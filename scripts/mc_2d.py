#!/usr/bin/env python
"""mc_2d.py
二维蜂窝磁体 经典 Heisenberg Monte Carlo, 计算序参量 vs 温度 -> 临界温度 (T_C / T_N)。

自包含: 只依赖 numpy + numba + matplotlib, 不依赖 spin_mc 共享引擎。
  - 自写 nopython 快速 kernel (xorshift64 RNG)。
  - 支持**总磁化** (铁磁, 算 T_C)、**交错磁化** (反铁磁, 算 T_N) 与
    **逐子晶格 Néel 矢量** (`--order sublattice:<mask>`, 交变磁体 P5) 三种序参量。
  - 序参量路径额外算 **Binder 累积量 U4** (第 4 列): 有序 -> 0.4, 高温无序 -> 0;
    不同 --size 的 U4(T) 交点是有限尺寸标度定 T_C 的稳健判据。
  - --seeds N: 跑 N 个独立种子 (不同初态/MC 状态) 并平均, 报告峰位标准差 (误差棒)。

用法:
  python mc_2d.py model.json --order auto --size 30 30 --tmin 5 --tmax 150 --nt 30 \
      --therm 2500 --meas 6000 --seeds 5 --out out_prefix
  --order auto: J1<0 -> 交错 (反铁磁), J1>0 -> 总磁化 (铁磁)
  --order sublattice:0     交变磁体: A = 子晶格 0, B = 其余 (与老的交错序参量同构)
  --order sublattice:0,1   交变磁体: A = 子晶格 {0,1}, B = 其余
"""

import argparse
import json
import math
import os

import numpy as np
from numba import njit


# --------------------------------------------------------------------------- lattice
def build_arrays(model, size, spin_seed=0):
    """从 poscar_to_mc3.py 的模型 JSON 直接构建 numba 所需数组。spin_seed 只决定初始自旋。"""
    prim = np.array(model["primitive_vectors"], dtype=float)
    basis = [np.array(f, dtype=float) for f in model["basis_frac"]]
    num_basis = len(basis)
    Lx, Ly = int(size[0]), int(size[1])
    N = Lx * Ly * num_basis

    spins = np.zeros((N, 3))
    onsite = np.zeros((N, 3, 3))
    for b, A in enumerate(model["onsite"]):
        onsite[b::num_basis] = np.array(A, dtype=float)

    # 邻接表 (先用 list 收集)
    nb = [[] for _ in range(N)]
    jm = [[] for _ in range(N)]
    for it in model["interactions"]:
        b1 = int(it["b1"]); b2 = int(it["b2"])
        off = np.array(it["offset"], dtype=int)
        J = np.array(it["J"], dtype=float)
        for cx in range(Lx):
            for cy in range(Ly):
                i = (cx * Ly + cy) * num_basis + b1
                tc = (np.array([cx, cy]) + off) % np.array([Lx, Ly])
                j = (tc[0] * Ly + tc[1]) * num_basis + b2
                nb[i].append(j); jm[i].append(J)
                nb[j].append(i); jm[j].append(J.T)

    maxn = max((len(x) for x in nb), default=0)
    neighbor_array = np.full((N, maxn), -1, dtype=np.int32)
    J_array = np.zeros((N, maxn, 3, 3), dtype=np.float64)
    for i in range(N):
        for k in range(len(nb[i])):
            neighbor_array[i, k] = nb[i][k]
            J_array[i, k] = jm[i][k]

    # 随机初始自旋
    rng = np.random.default_rng(spin_seed)
    phi = 2.0 * np.pi * rng.random(N)
    z = 2.0 * rng.random(N) - 1.0
    r = np.sqrt(1.0 - z * z)
    spins[:, 0] = r * np.cos(phi)
    spins[:, 1] = r * np.sin(phi)
    spins[:, 2] = z

    return spins, neighbor_array, J_array, onsite, num_basis


def stag_signs(neighbors, J_array, num_basis):
    """交错序参量的格点子晶格符号 (±1), 或 None (非二分/无法着色)。

    历史实现把「交错」写死成 `i % num_basis == 0` 为子晶格 A、其余为 B, 隐含
    **num_basis >= 2**。方形/矩形单离子晶格 num_basis == 1 ⇒ 所有格点都落在 A 上
    ⇒ cb == 0 ⇒ _run_T_fast 里 mbx/cb 抛 ZeroDivisionError (La2CuO4 是第一个
    J1 < 0 且 num_basis == 1 的体系, 所以这个坑一直没被踩到)。

    这里的规则:
      num_basis >= 2 -> 保持原样 (i % num_basis == 0 为 +1), 老体系数值**逐位不变**;
      num_basis == 1 -> 改用**图二着色**, 只用最强的一层键 (|J| >= 0.95*max|J|):
                        方形/矩形网的棋盘解就是真正的 Néel 子晶格。
                        奇数环 (三角/受挫) 着色冲突 -> None, 由调用方报错并建议 --order cv。
    """
    N = neighbors.shape[0]
    if num_basis >= 2:
        return np.where(np.arange(N) % num_basis == 0, 1.0, -1.0)
    M = neighbors.shape[1]
    jmax = 0.0
    for i in range(N):
        for k in range(M):
            if neighbors[i, k] < 0:
                break
            v = abs(J_array[i, k, 0, 0])
            if v > jmax:
                jmax = v
    if jmax <= 0.0:
        return None
    sign = np.zeros(N, dtype=np.float64)
    seen = np.zeros(N, dtype=bool)
    sign[0] = 1.0
    seen[0] = True
    stack = [0]
    while stack:
        i = stack.pop()
        for k in range(M):
            j = int(neighbors[i, k])
            if j < 0:
                break
            if abs(J_array[i, k, 0, 0]) < 0.95 * jmax:
                continue                      # 非主导壳, 不参与着色
            if not seen[j]:
                seen[j] = True
                sign[j] = -sign[i]
                stack.append(j)
            elif sign[j] == sign[i]:
                return None                   # 奇数环 -> 非二分
    if not seen.all():
        return None                           # 不连通 (不应发生)
    return sign


def sublattice_signs(num_basis, size, mask):
    """**逐子晶格 Néel 矢量**序参量的格点符号 (±1): mask 里的子晶格为 A (+1), 其余为 B (-1)。

    P5 (交变磁体) / NEW_TYPES_PLAN.md 第 24 项。与 `stag_signs` 的 `num_basis >= 2` 分支同构:
    **`mask = {0}` 就是老的交错序参量** ⇒ 双亚晶格体系(蜂窝/FGT/三角)逐位不变, 可作回归锚点。
    与 `stag_signs` 的区别是分组由调用方给(而不是写死 basis 0 vs 其余), 所以 basis >= 3
    (三角/kagome)也能取任意子晶格组合 —— 交变磁体要的就是「两套子晶格各自一个 Néel 矢量」。

    B 组为空(掩码覆盖全部子晶格)会让 `_run_T_fast` 里 `mbx/cb` 除零 ⇒ 这里直接报错(不静默)。
    """
    if num_basis < 2:
        raise ValueError(f"num_basis={num_basis}: 单亚晶格模型没有可分的子晶格 (需要 >= 2)")
    mask = sorted({int(b) for b in mask})
    if not mask:
        raise ValueError("子晶格掩码不能为空")
    bad = [b for b in mask if b < 0 or b >= num_basis]
    if bad:
        raise ValueError(f"子晶格编号越界: {bad} (num_basis={num_basis}, 允许 0..{num_basis - 1})")
    if len(mask) >= num_basis:
        raise ValueError(f"掩码 {mask} 覆盖了全部 {num_basis} 个子晶格 ⇒ B 组为空 (会让 mbx/cb 除零)")
    N = int(size[0]) * int(size[1]) * num_basis
    sel = np.zeros(num_basis, dtype=bool)
    sel[mask] = True
    return np.where(sel[np.arange(N) % num_basis], 1.0, -1.0)


ORDER_HELP = ("auto | staggered | total | cv | sublattice:<mask>  "
              "(mask = 逗号分隔的子晶格号, 例: sublattice:0 或 sublattice:0,1)")


def parse_order_spec(spec):
    """把 --order 的值解析成 (kind, mask)。kind ∈ auto/staggered/total/cv/sublattice。

    掩码解析是**严格**的: 空项 (如 `sublattice:0,` 或 `sublattice:0,,1`) 直接报错, 不静默丢掉 ——
    掩码是多打一个逗号就会把分组换掉的参数, 宁可让人重打一遍。
    """
    s = str(spec).strip()
    if s.split(":", 1)[0] == "sublattice":
        rest = s.split(":", 1)[1] if ":" in s else ""
        if not rest.strip():
            raise ValueError("--order sublattice:<mask> 需要掩码, 例如 sublattice:0 或 sublattice:0,1")
        toks = rest.replace(" ", "").split(",")
        if any(t == "" for t in toks):
            raise ValueError(f"--order {spec!r} 的掩码有空项 (写成 sublattice:0 或 sublattice:0,1)")
        try:
            mask = tuple(int(x) for x in toks)
        except ValueError:
            raise ValueError(f"--order {spec!r} 的掩码不是整数列表")
        return "sublattice", mask
    return s, None


@njit
def _xorshift64(state):
    x = state[0]
    x ^= x >> np.uint64(12)
    x ^= x << np.uint64(25)
    x ^= x >> np.uint64(27)
    state[0] = x
    return float(x * np.uint64(2685821657736338717)) * 5.421010862427522e-20


@njit
def _run_T_fast(spins, neighbors, J, onsite, T, n_therm, n_meas, state, order, stag):
    N = spins.shape[0]
    M = neighbors.shape[1]
    if T < 1e-6:
        T = 1e-6
    beta = 1.0 / (0.08617333262 * T)
    s1 = 0.0
    s2 = 0.0
    s4 = 0.0
    for it in range(n_therm + n_meas):
        for _ in range(N):
            i = int(_xorshift64(state) * N)
            if i >= N:
                i = N - 1
            oldx = spins[i, 0]
            oldy = spins[i, 1]
            oldz = spins[i, 2]
            phi = 6.283185307179586 * _xorshift64(state)
            cz = 2.0 * _xorshift64(state) - 1.0
            r = math.sqrt(1.0 - cz * cz)
            nx = r * math.cos(phi)
            ny = r * math.sin(phi)
            nz = cz
            dsx = nx - oldx
            dsy = ny - oldy
            dsz = nz - oldz
            dE = 0.0
            for k in range(M):
                j = neighbors[i, k]
                if j == -1:
                    break
                Jm = J[i, k]
                sx = spins[j, 0]
                sy = spins[j, 1]
                sz = spins[j, 2]
                dE += dsx * (Jm[0, 0] * sx + Jm[0, 1] * sy + Jm[0, 2] * sz) \
                    + dsy * (Jm[1, 0] * sx + Jm[1, 1] * sy + Jm[1, 2] * sz) \
                    + dsz * (Jm[2, 0] * sx + Jm[2, 1] * sy + Jm[2, 2] * sz)
            A = onsite[i]
            dE += (nx * A[0, 0] * nx + nx * A[0, 1] * ny + nx * A[0, 2] * nz
                   + ny * A[1, 0] * nx + ny * A[1, 1] * ny + ny * A[1, 2] * nz
                   + nz * A[2, 0] * nx + nz * A[2, 1] * ny + nz * A[2, 2] * nz) \
                - (oldx * A[0, 0] * oldx + oldx * A[0, 1] * oldy + oldx * A[0, 2] * oldz
                   + oldy * A[1, 0] * oldx + oldy * A[1, 1] * oldy + oldy * A[1, 2] * oldz
                   + oldz * A[2, 0] * oldx + oldz * A[2, 1] * oldy + oldz * A[2, 2] * oldz)
            if dE <= 0.0 or _xorshift64(state) < math.exp(-beta * dE):
                spins[i, 0] = nx
                spins[i, 1] = ny
                spins[i, 2] = nz
        if it >= n_therm:
            if order == 1:  # staggered: 子晶格符号由 stag 给出 (= +/-1)
                max_ = 0.0; may = 0.0; maz = 0.0
                mbx = 0.0; mby = 0.0; mbz = 0.0
                ca = 0; cb = 0
                for i in range(N):
                    if stag[i] > 0.0:
                        max_ += spins[i, 0]; may += spins[i, 1]; maz += spins[i, 2]; ca += 1
                    else:
                        mbx += spins[i, 0]; mby += spins[i, 1]; mbz += spins[i, 2]; cb += 1
                if ca == 0 or cb == 0:
                    # 单子晶格 (num_basis == 1 且未着色成功) -> 上游已拦, 这里兜底不崩
                    m = 0.0
                else:
                    stx = 0.5 * (max_ / ca - mbx / cb)
                    sty = 0.5 * (may / ca - mby / cb)
                    stz = 0.5 * (maz / ca - mbz / cb)
                    m = math.sqrt(stx * stx + sty * sty + stz * stz)
            else:  # total
                tx = 0.0; ty = 0.0; tz = 0.0
                for i in range(N):
                    tx += spins[i, 0]; ty += spins[i, 1]; tz += spins[i, 2]
                m = math.sqrt(tx * tx + ty * ty + tz * tz) / N
            s1 += m
            s2 += m * m
            s4 += m * m * m * m
    mean_m = s1 / n_meas
    m2_mean = s2 / n_meas
    var = m2_mean - mean_m * mean_m
    chi = var * N / (0.08617 * T)
    # Binder 累积量: U4 = 1 - <M^4> / ((1+2/n) <M^2>^2), n=3 (Heisenberg 序参量)。
    # M = 序参量模长 |ΣS|/N。有序(低温, M->1) -> 1-3/5 = 0.4; 高温无序(O(3) 高斯) -> 0。
    # 不同 L 的 U4(T) 交点是有限尺寸标度定 T_C 的稳健判据 (见 --seeds 误差棒与 README)。
    if m2_mean > 1e-16:
        u4 = 1.0 - (s4 / n_meas) / ((5.0 / 3.0) * m2_mean * m2_mean)
    else:
        u4 = 0.4  # 完全有序极限
    return mean_m, chi, u4


@njit
def _energy(spins, neighbors, J, onsite):
    """总能量/格点 (交换 + 单离子)。有向键和双计, 除 2。"""
    N = spins.shape[0]
    M = neighbors.shape[1]
    E = 0.0
    for i in range(N):
        for k in range(M):
            j = neighbors[i, k]
            if j == -1:
                break
            Jm = J[i, k]
            sx = spins[j, 0]; sy = spins[j, 1]; sz = spins[j, 2]
            E += spins[i, 0] * (Jm[0, 0] * sx + Jm[0, 1] * sy + Jm[0, 2] * sz) \
               + spins[i, 1] * (Jm[1, 0] * sx + Jm[1, 1] * sy + Jm[1, 2] * sz) \
               + spins[i, 2] * (Jm[2, 0] * sx + Jm[2, 1] * sy + Jm[2, 2] * sz)
    return 0.5 * E / N


@njit
def _run_T_cv(spins, neighbors, J, onsite, T, n_therm, n_meas, state):
    """单温度 Metropolis + 比热 Cv 测量 (能量涨落)。用于受挫/非二分晶格 (三角) 定 T_N。"""
    N = spins.shape[0]
    M = neighbors.shape[1]
    if T < 1e-6:
        T = 1e-6
    beta = 1.0 / (0.08617333262 * T)
    e1 = 0.0
    e2 = 0.0
    for it in range(n_therm + n_meas):
        for _ in range(N):
            i = int(_xorshift64(state) * N)
            if i >= N:
                i = N - 1
            oldx = spins[i, 0]; oldy = spins[i, 1]; oldz = spins[i, 2]
            phi = 6.283185307179586 * _xorshift64(state)
            cz = 2.0 * _xorshift64(state) - 1.0
            r = math.sqrt(1.0 - cz * cz)
            nx = r * math.cos(phi); ny = r * math.sin(phi); nz = cz
            dsx = nx - oldx; dsy = ny - oldy; dsz = nz - oldz
            dE = 0.0
            for k in range(M):
                j = neighbors[i, k]
                if j == -1:
                    break
                Jm = J[i, k]
                sx = spins[j, 0]; sy = spins[j, 1]; sz = spins[j, 2]
                dE += dsx * (Jm[0, 0] * sx + Jm[0, 1] * sy + Jm[0, 2] * sz) \
                    + dsy * (Jm[1, 0] * sx + Jm[1, 1] * sy + Jm[1, 2] * sz) \
                    + dsz * (Jm[2, 0] * sx + Jm[2, 1] * sy + Jm[2, 2] * sz)
            A = onsite[i]
            dE += (nx * A[0, 0] * nx + nx * A[0, 1] * ny + nx * A[0, 2] * nz
                   + ny * A[1, 0] * nx + ny * A[1, 1] * ny + ny * A[1, 2] * nz
                   + nz * A[2, 0] * nx + nz * A[2, 1] * ny + nz * A[2, 2] * nz) \
                - (oldx * A[0, 0] * oldx + oldx * A[0, 1] * oldy + oldx * A[0, 2] * oldz
                   + oldy * A[1, 0] * oldx + oldy * A[1, 1] * oldy + oldy * A[1, 2] * oldz
                   + oldz * A[2, 0] * oldx + oldz * A[2, 1] * oldy + oldz * A[2, 2] * oldz)
            if dE <= 0.0 or _xorshift64(state) < math.exp(-beta * dE):
                spins[i, 0] = nx; spins[i, 1] = ny; spins[i, 2] = nz
        if it >= n_therm:
            e = _energy(spins, neighbors, J, onsite)
            e1 += e
            e2 += e * e
    mean_e = e1 / n_meas
    cv = (e2 / n_meas - mean_e * mean_e) * N / (0.08617 * T * T)
    return mean_e, cv


def run_sweep(spins, neighbors, J, onsite, T_list, therm, meas, stag, order, seed=0):
    state = np.array([np.uint64(seed if seed != 0 else 88172645463325252)], dtype=np.uint64)
    results = []
    for T in T_list:
        mean_m, chi, u4 = _run_T_fast(spins, neighbors, J, onsite, float(T), therm, meas,
                                      state, order, stag)
        results.append([T, mean_m, chi, u4])
        print(f"  T={T:6.1f} K | m={mean_m:.4f} | chi={chi:.4f} | U4={u4:.4f}", flush=True)
    return np.array(results)


def run_sweep_cv(spins, neighbors, J, onsite, T_list, therm, meas, seed=0):
    state = np.array([np.uint64(seed if seed != 0 else 88172645463325252)], dtype=np.uint64)
    results = []
    for T in T_list:
        mean_e, cv = _run_T_cv(spins, neighbors, J, onsite, float(T), therm, meas, state)
        results.append([T, mean_e, cv])
        print(f"  T={T:6.1f} K | E/N={mean_e:.4f} | Cv={cv:.4f}", flush=True)
    return np.array(results)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("--order", default="auto", help=ORDER_HELP)
    ap.add_argument("--size", type=int, nargs=2, default=[30, 30])
    ap.add_argument("--tmin", type=float, default=5.0)
    ap.add_argument("--tmax", type=float, default=150.0)
    ap.add_argument("--nt", type=int, default=30)
    ap.add_argument("--therm", type=int, default=2500)
    ap.add_argument("--meas", type=int, default=6000)
    ap.add_argument("--cool", action="store_true")
    ap.add_argument("--seed", type=int, default=0,
                    help="MC 随机数种子; 0 = 默认 (可复现)。改种子可估计统计误差")
    ap.add_argument("--seeds", type=int, default=1,
                    help="独立种子数 N: 跑 N 次不同初态并平均, 报告峰位标准差 (误差棒)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    # 必须显式 utf-8: Windows 默认用区域编码 (GBK), 模型 JSON 里若含非 ASCII
    # (例如中文路径的 source_poscar, 或手改过的注释) 会被误解码而报
    # "Invalid \escape"。
    with open(args.model, encoding="utf-8") as f:
        model = json.load(f)
    meta = model["meta"]

    try:
        o_kind, o_mask = parse_order_spec(args.order)
    except ValueError as e:
        raise SystemExit(f"[error] {e}")
    if o_kind not in ("auto", "staggered", "total", "cv", "sublattice"):
        raise SystemExit(f"[error] 未知 --order {args.order!r}: 可用 {ORDER_HELP}")

    if o_kind == "cv":
        use_cv = True
        order = 0
        kind = "比热 Cv 峰 (受挫/三角晶格, T_N)"
    else:
        use_cv = False
        if o_kind == "auto":
            # 单壳层模型给 J1_meV; 多壳层模型 (FGT 双亚晶格) 没有这个键 ->
            # 用各壳层 J 之和的符号判断铁磁/反铁磁, 不能硬取 meta["J1_meV"]
            if "J1_meV" in meta:
                jsum = float(meta["J1_meV"])
            else:
                jsum = sum(-float(it["J"][0][0]) for it in model.get("interactions", []))
            order = 1 if jsum < 0 else 0
            kind = "交错 (反铁磁, T_N)" if order == 1 else "总磁化 (铁磁, T_C)"
        elif o_kind == "sublattice":
            order = 1
            kind = (f"逐子晶格 Néel 矢量 |L| (A = 子晶格 {sorted(set(o_mask))}, "
                    f"B = 其余; 交变磁体序参量, P5)")
        else:
            order = 1 if o_kind == "staggered" else 0
            kind = "交错 (反铁磁, T_N)" if order == 1 else "总磁化 (铁磁, T_C)"
    print(f"序参量: {kind}")

    prefix = args.out or os.path.splitext(os.path.basename(args.model))[0]
    # neighbors/J/onsite/num_basis 与种子无关; 先建一次拿 N/num_basis, 各种子只重建 spins。
    _spins, neighbors, J, onsite, num_basis = build_arrays(model, args.size, spin_seed=0)
    N = _spins.shape[0]
    # 单壳层模型 (poscar_to_mc3.py) 有 J1_meV/J2_meV; 多壳层模型 (如 FGT 双亚晶格)
    # 只给 shells_summary 或什么都没有 -> 不能硬取, 否则 KeyError。
    if "J1_meV" in meta:
        jinfo = (f"J1={meta['J1_meV']} J2={meta['J2_meV']} J3={meta.get('J3_meV', 0.0)}")
    elif meta.get("shells_summary"):
        jinfo = str(meta["shells_summary"])
    else:
        jinfo = f"{len(model.get('interactions', []))} 根键/胞 (多壳层)"
    print(f"格点: {N} (basis {num_basis} x {tuple(args.size)}) | "
          f"{jinfo} | MAE={meta.get('MAE_meV', 0.0)} meV | S={meta.get('S', 1.0)}")

    # 交错序参量需要子晶格 ±1 符号; num_basis == 1 (方形/矩形单离子) 必须图二着色,
    # 不能再用「i % num_basis」——那会让 cb == 0 而除零 (见 stag_signs docstring)。
    stag = np.zeros(N, dtype=np.float64)   # order == 0 (总磁化) 用不到; 占位避免 numba 的 None 类型
    if not use_cv and order == 1:
        if o_kind == "sublattice":
            # 交变磁体: 分组由掩码给 (不再写死 basis 0 vs 其余)
            try:
                stag = sublattice_signs(num_basis, args.size, o_mask)
            except ValueError as e:
                raise SystemExit(f"[error] --order {args.order}: {e}")
            n_a = int((stag > 0.0).sum())
            print(f"子晶格 (Néel 矢量序参量): A={n_a} B={N - n_a} | 掩码 "
                  f"{sorted(set(int(b) for b in o_mask))} / basis {num_basis}")
        else:
            stag = stag_signs(neighbors, J, num_basis)
            if stag is None:
                raise SystemExit(
                    f"序参量=交错 需要二分子晶格, 但该模型 (num_basis={num_basis}) 无法二分着色"
                    f" (受挫/非二分晶格) —— 请改用 `--order cv` 用比热峰定 T_N, "
                    f"或用 `--order total` 看总磁化")
            n_a = int((stag > 0.0).sum())
            print(f"子晶格 (交错序参量): A={n_a} B={N - n_a} | 图二着色, 主导壳判据 |J| >= 0.95*max|J|")
            if num_basis >= 3:
                print(f"[warn] num_basis={num_basis}: 交错序参量只区分 basis 0 与其余 "
                      f"(basis 1..{num_basis - 1} 被混成一组 B) —— 需要别的分组请用 "
                      f"`--order sublattice:<mask>`")

    T_list = np.linspace(args.tmin, args.tmax, args.nt)
    if args.cool:
        T_list = T_list[::-1]

    # 多种子: 每个种子独立初态 + 独立 MC 状态, 事后平均曲线 + 收集峰位 -> 误差棒。
    n_seeds = max(1, int(args.seeds))
    all_res = []
    peaks = {"chi": [], "dm": [], "cv": []}
    for k in range(n_seeds):
        spins, _, _, _, _ = build_arrays(model, args.size, spin_seed=args.seed + k)
        sd = args.seed + k
        if use_cv:
            res = run_sweep_cv(spins, neighbors, J, onsite, T_list, args.therm, args.meas,
                               seed=sd)
        else:
            res = run_sweep(spins, neighbors, J, onsite, T_list, args.therm, args.meas,
                            stag, order, seed=sd)
        res = res[np.argsort(res[:, 0])]
        all_res.append(res)
        if use_cv:
            peaks["cv"].append(res[int(np.argmax(res[:, 2])), 0])
        else:
            peaks["chi"].append(res[int(np.argmax(res[:, 2])), 0])
            peaks["dm"].append(res[int(np.argmax(-np.gradient(res[:, 1], res[:, 0]))), 0])
    res = np.mean(np.stack(all_res), axis=0)
    if n_seeds > 1:
        print(f">>> 已平均 {n_seeds} 个独立种子 (初态/MC 状态各不同)")

    if use_cv:
        fn = f"{prefix}_mt.csv"
        np.savetxt(fn, res, delimiter=",", header="T_Kelvin,E_per_site,Cv", comments="")
        print(f"\n>>> 数据存至 {fn}")
        i_cv = int(np.argmax(res[:, 2]))
        cv_peak = res[i_cv, 0]
        cv_edge = i_cv in (0, len(res) - 1)
        if cv_edge:
            print(f">>> [warn] 比热 Cv 峰落在温度窗口边界 T={cv_peak:.2f} K → T_N **未分辨**; "
                  f"这是窗口截断/网格跨过窄峰的假象, 请收窄 --tmin/--tmax、加密 --nt "
                  f"(低于 T_N 启动时再加 --cool)")
        else:
            msg = f">>> 比热 Cv 峰处 T_N ≈ {cv_peak:.2f} K"
            if n_seeds > 1:
                msg += f" ± {np.std(peaks['cv']):.2f} K ({n_seeds} 种子峰位标准差)"
            print(msg)
    else:
        fn = f"{prefix}_mt.csv"
        np.savetxt(fn, res, delimiter=",", header="T_Kelvin,order_param,Susceptibility,Binder_U4",
                   comments="")
        print(f"\n>>> 数据存至 {fn}")
        i_chi = int(np.argmax(res[:, 2]))
        dm = -np.gradient(res[:, 1], res[:, 0])
        i_dm = int(np.argmax(dm))
        chi_peak, tn_dm = res[i_chi, 0], res[i_dm, 0]
        chi_edge, dm_edge = i_chi in (0, len(res) - 1), i_dm in (0, len(res) - 1)
        se = (lambda name: f" ± {np.std(peaks[name]):.2f} K ({n_seeds} 种子峰位标准差)"
              if n_seeds > 1 else "")
        print(f">>> 磁化率峰处 T ≈ {chi_peak:.2f} K{se('chi')}"
              + (" [warn] 落在窗口边界 → 未分辨" if chi_edge else ""))
        print(f">>> |dm/dT| 最大处 T ≈ {tn_dm:.2f} K{se('dm')}"
              + (" [warn] 落在窗口边界 → 未分辨" if dm_edge else ""))
        print(f">>> Binder U4 见第 4 列: 有序 -> 0.4, 高温无序 -> 0; "
              f"不同 --size 的 U4(T) 交点是有限尺寸标度定 T_C 的稳健判据")
        if chi_edge or dm_edge:
            print(">>> [warn] 峰判据落在温度窗口边界: 请收窄 --tmin/--tmax、加密 --nt "
                  "(低温侧加 --cool); 该值不能当作 T_C/T_N")

    png = f"{prefix}_mt.png"
    tc_label = "T_N" if order == 1 else "T_C"
    # 两个版式各自只定义一半的峰量 (Cv 版没有 chi_*, 默认版没有 cv_*); 纯 Pillow 兜底
    # 函数按 use_cv 只取其中一半, 这里先补齐缺省值, 免得 UnboundLocalError 把出图吞掉。
    cv_peak = cv_edge = chi_peak = chi_edge = None
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        if use_cv:
            fig, ax1 = plt.subplots(figsize=(6.5, 4.5))
            ax1.plot(res[:, 0], res[:, 1], "o-", color="tab:blue", label="E/N (meV)")
            ax1.set_xlabel("T (K)")
            ax1.set_ylabel("energy per site", color="tab:blue")
            ax1.tick_params(axis="y", labelcolor="tab:blue")
            ax2 = ax1.twinx()
            ax2.plot(res[:, 0], res[:, 2], "s--", color="tab:red", alpha=0.7, label="Cv")
            ax2.set_ylabel("specific heat $C_v$", color="tab:red")
            ax2.tick_params(axis="y", labelcolor="tab:red")
            # 标题保持 ASCII: matplotlib 默认字体 (DejaVu Sans) 没有 CJK 字形, 中文会变豆腐块
            t_lbl = "T_N unresolved (peak at window edge)" if cv_edge \
                else f"$T_N\\approx${cv_peak:.0f} K"
            ax1.set_title(f"{meta['element']} MC (Cv)  ({t_lbl})")
        else:
            fig, (ax1, ax3) = plt.subplots(2, 1, figsize=(6.5, 6.5), sharex=True)
            ylab = "staggered $m_{stag}$" if order == 1 else "magnetization $M$"
            tc_label = "T_N" if order == 1 else "T_C"
            ax1.plot(res[:, 0], res[:, 1], "o-", color="tab:blue", label=ylab)
            ax1.set_ylabel(ylab, color="tab:blue")
            ax1.tick_params(axis="y", labelcolor="tab:blue")
            ax2 = ax1.twinx()
            ax2.plot(res[:, 0], res[:, 2], "s--", color="tab:red", alpha=0.7,
                     label="susceptibility")
            ax2.set_ylabel("susceptibility $\\chi$", color="tab:red")
            ax2.tick_params(axis="y", labelcolor="tab:red")
            t_lbl = (f"{tc_label} unresolved (peak at window edge)" if chi_edge
                     else f"${{{tc_label}}}\\approx${chi_peak:.0f} K")
            ax1.set_title(f"{meta['element']} MC  ({t_lbl})")
            # Binder 累积量子图 (有限尺寸标度判据): 有序 -> 0.4, 无序 -> 0
            ax3.plot(res[:, 0], res[:, 3], "^-", color="tab:green", label="Binder $U_4$")
            ax3.axhline(0.4, color="gray", ls=":", lw=0.8)
            ax3.set_xlabel("T (K)")
            ax3.set_ylabel("Binder $U_4$")
            ax3.grid(alpha=0.3)
        ax1.grid(alpha=0.3)
        fig.savefig(png, dpi=150, bbox_inches="tight")
        print(f">>> 图存至 {png}")
    except Exception as e:
        print(f"[warn] matplotlib 画图不可用 ({e})")
        try:
            _plot_mt_png(res, meta, png, order, use_cv, cv_peak, cv_edge,
                         chi_peak, chi_edge, tc_label)
            print(f">>> 图存至 {png} (纯 Pillow 兜底)")
        except Exception as e2:
            print(f"[warn] 画图跳过: {e2}")


def _pil_font(size_px):
    """mc_2d 自包含的字体加载: matplotlib 不在时也要能写中文 / 希腊字母。"""
    import os
    from PIL import ImageFont
    for p in (r"C:\Windows\Fonts\msyh.ttc", r"C:\Windows\Fonts\simhei.ttf",
              r"C:\Windows\Fonts\Deng.ttf", r"C:\Windows\Fonts\DejaVuSans.ttf"):
        if os.path.exists(p):
            try:
                return ImageFont.truetype(p, int(size_px))
            except OSError:
                continue
    return ImageFont.load_default()


def _dashed_poly(dr, pts, color, width, dash, gap):
    """手工虚/点折线 (Pillow 没有线型概念)。"""
    for i in range(len(pts) - 1):
        (x0, y0), (x1, y1) = pts[i], pts[i + 1]
        tot = ((x1 - x0) ** 2 + (y1 - y0) ** 2) ** 0.5
        if tot <= 0:
            continue
        ux, uy = (x1 - x0) / tot, (y1 - y0) / tot
        s = 0.0
        while s < tot:
            e = min(s + dash, tot)
            dr.line([(x0 + ux * s, y0 + uy * s), (x0 + ux * e, y0 + uy * e)],
                    fill=color, width=width)
            s = e + gap


def _plot_mt_png(res, meta, png, order, use_cv, cv_peak, cv_edge,
                 chi_peak, chi_edge, tc_label):
    """纯 Pillow 版 mc_mt.png (matplotlib 不可用时的兜底; 见本文件头 Pillow 说明)。

    版式与 matplotlib 版一致:
      - Cv 模式: 单面板, 左轴 E/N (蓝实线+圆点) / 右轴 Cv (红虚线+方点);
      - 默认  : 上=磁化(或交错磁化, 蓝) + 磁化率(红, 右轴), 下=Binder U4 (绿) + U4=0.4 参考线。
    列定义 (与 matplotlib 版相同): 0=T, 1=序参量/能量, 2=磁化率/Cv, 3=Binder U4。
    """
    import os
    from PIL import Image, ImageDraw
    SS, DPI = 2, 150
    CB, CR, CG, CGRID = (31, 119, 180), (214, 39, 40), (44, 160, 44), (226, 226, 226)
    fig_w, fig_h = (6.5, 4.5) if use_cv else (6.5, 6.5)
    W, H = int(fig_w * DPI), int(fig_h * DPI)
    Wd, Hd = W * SS, H * SS
    img = Image.new("RGB", (Wd, Hd), "white")
    dr = ImageDraw.Draw(img)
    L, R = 0.145 * Wd, 0.845 * Wd

    f_t = _pil_font(0.030 * Hd)
    f_l = _pil_font(0.0215 * Hd)
    f_k = _pil_font(0.0185 * Hd)

    def fmt(v):
        return f"{v:.3g}"

    xs = np.asarray(res[:, 0], dtype=float)
    x1, x2 = float(xs.min()), float(xs.max())
    if x2 - x1 <= 0:
        x2 = x1 + 1.0
    pad = 0.03 * (x2 - x1)
    xlo, xhi = x1 - pad, x2 + pad

    def rng(vals):
        v = np.asarray(vals, dtype=float)
        lo, hi = float(v.min()), float(v.max())
        if hi - lo < 1e-12:
            hi = lo + 1.0
        d = 0.08 * (hi - lo)
        return lo - d, hi + d

    def tk(lo, hi, n=5):
        return [lo + (hi - lo) * i / float(n - 1) for i in range(n)]

    def xpix(x):
        return L + (x - xlo) / (xhi - xlo) * (R - L)

    def make_ypix(box, lo, hi):
        y0, y1 = box
        return lambda v: y1 - (v - lo) / (hi - lo) * (y1 - y0)

    def frame(box, xlabel=None, show_x=True):
        y0, y1 = box
        for t in tk(xlo, xhi):
            xx = xpix(t)
            dr.line([(xx, y0), (xx, y1)], fill=CGRID, width=max(1, int(round(SS * 0.8))))
            if show_x:
                dr.text((xx, y1 + 0.016 * Hd), fmt(t), font=f_k, fill=(70, 70, 70),
                        anchor="ma")
        if xlabel:
            dr.text(((L + R) / 2.0, y1 + 0.050 * Hd), xlabel, font=f_l,
                    fill=(35, 35, 35), anchor="ma")

    def border(box):
        dr.rectangle([L, box[0], R, box[1]], outline=(130, 130, 130),
                     width=max(1, int(round(SS))))

    def yaxis(box, ypix, lo, hi, color, label, side="left"):
        y0, _y1 = box
        for t in tk(lo, hi):
            yy = ypix(t)
            dr.line([(L, yy), (R, yy)], fill=CGRID, width=max(1, int(round(SS * 0.8))))
            if side == "left":
                dr.text((L - 0.010 * Wd, yy), fmt(t), font=f_k, fill=(70, 70, 70),
                        anchor="rm")
            else:
                dr.text((R + 0.010 * Wd, yy), fmt(t), font=f_k, fill=color, anchor="lm")
        if side == "left":
            # 轴标题贴页边起排 (anchor=la/ra), 避免长标签被画布裁掉
            dr.text((0.008 * Wd, y0 - 0.034 * Hd), label, font=f_l, fill=color, anchor="la")
        else:
            dr.text((0.992 * Wd, y0 - 0.034 * Hd), label, font=f_l, fill=color, anchor="ra")

    def series(ypix, ys, color, marker, dashed=False):
        pts = [(xpix(float(a)), ypix(float(b))) for a, b in zip(xs, ys)]
        if dashed:
            _dashed_poly(dr, pts, color, max(1, int(round(SS * 1.4))), SS * 4.5, SS * 2.6)
        else:
            for i in range(len(pts) - 1):
                dr.line([pts[i], pts[i + 1]], fill=color,
                        width=max(1, int(round(SS * 1.4))))
        r = SS * 5.0
        for (px, py) in pts:
            if marker == "s":
                dr.rectangle([px - r, py - r, px + r, py + r], fill=color)
            elif marker == "^":
                dr.polygon([(px, py - r * 1.2), (px - r * 1.2, py + r),
                            (px + r * 1.2, py + r)], fill=color)
            else:
                dr.ellipse([px - r, py - r, px + r, py + r], fill=color)

    # 调用方只算自己版式的峰量: Cv 版给 cv_*, 默认版给 chi_*, 另一半是 None。
    # 这里按版式自己补算, 免得把 None 当成"峰落在窗口边界"或直接 format(None) 崩掉。
    if use_cv:
        if cv_peak is None:
            _i = int(np.argmax(res[:, 2]))
            cv_peak, cv_edge = float(res[_i, 0]), _i in (0, len(res) - 1)
    elif chi_peak is None:
        _i = int(np.argmax(res[:, 2]))
        chi_peak, chi_edge = float(res[_i, 0]), _i in (0, len(res) - 1)

    if use_cv:
        box = (0.12 * Hd, 0.845 * Hd)
        lo1, hi1 = rng(res[:, 1])
        lo2, hi2 = rng(res[:, 2])
        frame(box, xlabel="T (K)")
        yaxis(box, make_ypix(box, lo1, hi1), lo1, hi1, CB, "energy per site", "left")
        yaxis(box, make_ypix(box, lo2, hi2), lo2, hi2, CR, "Cv", "right")
        border(box)
        series(make_ypix(box, lo1, hi1), res[:, 1], CB, "o")
        series(make_ypix(box, lo2, hi2), res[:, 2], CR, "s", dashed=True)
        t_lbl = ("T_N unresolved (peak at window edge)" if cv_edge
                 else f"T_N≈{cv_peak:.0f} K")
        title = f"{meta['element']} MC (Cv)  ({t_lbl})"
    else:
        ylab = "staggered m_stag" if order == 1 else "magnetization M"
        box1 = (0.095 * Hd, 0.565 * Hd)
        box2 = (0.690 * Hd, 0.845 * Hd)
        lo1, hi1 = rng(res[:, 1])
        lo2, hi2 = rng(res[:, 2])
        lo3, hi3 = rng(list(res[:, 3]) + [0.4])
        frame(box1, show_x=False)
        yaxis(box1, make_ypix(box1, lo1, hi1), lo1, hi1, CB, ylab, "left")
        yaxis(box1, make_ypix(box1, lo2, hi2), lo2, hi2, CR, "susceptibility χ", "right")
        border(box1)
        series(make_ypix(box1, lo1, hi1), res[:, 1], CB, "o")
        series(make_ypix(box1, lo2, hi2), res[:, 2], CR, "s", dashed=True)
        frame(box2, xlabel="T (K)")
        yaxis(box2, make_ypix(box2, lo3, hi3), lo3, hi3, CG, "Binder U4", "left")
        border(box2)
        yp3 = make_ypix(box2, lo3, hi3)
        y04 = yp3(0.4)
        _dashed_poly(dr, [(L, y04), (R, y04)], (140, 140, 140),
                     max(1, int(round(SS * 0.8))), SS * 2.0, SS * 2.6)
        series(yp3, res[:, 3], CG, "^")
        t_lbl = (f"{tc_label} unresolved (peak at window edge)" if chi_edge
                 else f"{tc_label}≈{chi_peak:.0f} K")
        title = f"{meta['element']} MC  ({t_lbl})"

    dr.text((Wd / 2.0, 0.028 * Hd), title, font=f_t, fill="black", anchor="ma")
    foot = ("纯 Pillow 兜底渲染 (本机 matplotlib 缺依赖: pyparsing/cycler/contourpy/"
            "fontTools/kiwisolver)\n数据: 05_mc/mc_mt.csv")
    dr.multiline_text((Wd / 2.0, 0.982 * Hd), foot, font=_pil_font(0.016 * Hd),
                      fill=(120, 120, 120), anchor="md", align="center")
    img = img.resize((W, H), Image.LANCZOS)
    img.save(png, "PNG", optimize=True)


if __name__ == "__main__":
    main()
