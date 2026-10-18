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
import sys

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
                if j == i:
                    # 周期性回绕把这条相互作用映射回了它自己 (offset 是该轴尺寸的整数倍)。
                    # 不拦的话该格点的"自作用"会被当成近邻累加进能量/序参量, 静默错。
                    raise SystemExit(
                        f"[error] 相互作用 ({b1}->{b2}, offset="
                        f"({int(off[0])},{int(off[1])})) 在 --size {Lx} {Ly} 下回绕回同一"
                        f"格点 (cell {cx},{cy}); --size 太小。每个轴的尺寸必须严格大于该轴"
                        f"上 offset 的绝对值 (例如 offset (2,0) 需要 Lx>=3)。")
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
              "(mask = 逗号分隔的子晶格号, 例: sublattice:0 或 sublattice:0,1; "
              "auto = 先比 T=0 候选基态能量 (fm/neel/stripe) 再决定用哪个序参量)")


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
    acc = 0
    m_series = np.empty(n_meas, dtype=np.float64)
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
                acc += 1
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
            m_series[it - n_therm] = m
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
    return mean_m, chi, u4, acc / float((n_therm + n_meas) * N), m_series


@njit
def _energy(spins, neighbors, J, onsite):
    """总能量/格点 (交换 + 单离子各向异性), 单位 meV。

    两项的双计约定**不同**, 不能合并后统一乘 0.5:
      * 交换项: 有向键 (i->j 与 j->i 各存一条) 都数了一遍 -> 乘 0.5;
      * 单离子项: 每个格点只数一次 -> **不**乘 0.5。
    必须与 _run_T_cv 里 Metropolis 的 dE 用同一个哈密顿量 (那边 onsite 也是每格点
    一次), 否则抽样含各向异性、测量不含 -> 比热 Cv 的能量涨落少一块 -> T_N 系统
    偏移, 而且退出码仍然是 0 (静默错)。
    """
    N = spins.shape[0]
    M = neighbors.shape[1]
    E_ex = 0.0
    for i in range(N):
        for k in range(M):
            j = neighbors[i, k]
            if j == -1:
                break
            Jm = J[i, k]
            sx = spins[j, 0]; sy = spins[j, 1]; sz = spins[j, 2]
            E_ex += spins[i, 0] * (Jm[0, 0] * sx + Jm[0, 1] * sy + Jm[0, 2] * sz) \
                  + spins[i, 1] * (Jm[1, 0] * sx + Jm[1, 1] * sy + Jm[1, 2] * sz) \
                  + spins[i, 2] * (Jm[2, 0] * sx + Jm[2, 1] * sy + Jm[2, 2] * sz)
    E_on = 0.0
    for i in range(N):
        A = onsite[i]
        sx = spins[i, 0]; sy = spins[i, 1]; sz = spins[i, 2]
        E_on += sx * (A[0, 0] * sx + A[0, 1] * sy + A[0, 2] * sz) \
              + sy * (A[1, 0] * sx + A[1, 1] * sy + A[1, 2] * sz) \
              + sz * (A[2, 0] * sx + A[2, 1] * sy + A[2, 2] * sz)
    return 0.5 * E_ex / N + E_on / N


def _brute_energy(spins, neighbors, J, onsite):
    """纯 Python 暴力求值, 只当作 _energy 的对照 (同一组数组、同一套双计约定)。"""
    N = spins.shape[0]
    M = neighbors.shape[1]
    e_ex = 0.0
    for i in range(N):
        for k in range(M):
            j = neighbors[i, k]
            if j == -1:
                break
            e_ex += float(spins[i] @ J[i, k] @ spins[j])
    e_on = 0.0
    for i in range(N):
        e_on += float(spins[i] @ onsite[i] @ spins[i])
    return 0.5 * e_ex / N + e_on / N


def _candidate_orders(size, num_basis, neighbors, J, onsite):
    """`--order auto` 的判据: 在受支持的候选基态里比 T=0 能量 (经典共线 Heisenberg)。

    为什么不能只拿 J1 的符号猜 (旧实现在这里会错): J1>0 (AFM) + J2<0 (FM) 的经典基态是
    **条纹**反铁磁 —— 它既不是 FM 也不是 Néel; 用 ±1 交错序参量测出来的峰不是 T_N,
    而模型仍会正常打印一个温度 (静默错)。这里返回每种候选的每格点能量, 由调用方比较。
    只列该格点几何下真有定义的候选:
      fm          全部平行 (对应 T_C)
      neel        二分着色 ±1 (要求 stag_signs 能着色, 否则不出现)
      stripe_x/y  (-1)^cx / (-1)^cy, 只对 num_basis == 1 的简单方格/矩形有意义
    """
    Lx, Ly = int(size[0]), int(size[1])
    N = Lx * Ly * num_basis
    cell = np.arange(N) // num_basis
    cx = cell // Ly
    cy = cell % Ly
    cands = {"fm": np.ones(N)}
    stag = stag_signs(neighbors, J, num_basis)
    if stag is not None:
        cands["neel"] = np.asarray(stag, dtype=np.float64)
    if num_basis == 1:
        cands["stripe_x"] = np.where(cx % 2 == 0, 1.0, -1.0)
        cands["stripe_y"] = np.where(cy % 2 == 0, 1.0, -1.0)
    out = {}
    for name, sign in cands.items():
        s = np.zeros((N, 3))
        s[:, 2] = sign
        out[name] = float(_energy(s, neighbors, J, onsite))
    return out


def _resolve_auto(neighbors, J, onsite, size, num_basis):
    """auto 判定的**唯一**实现: 返回 (order, kind, use_cv, cand, best)。

    mc_2d.main 与 mc_fss 都走这里, 免得两个脚本的判据各自漂移 (旧版两边都只看 J1 符号)。
    """
    cand = _candidate_orders(size, num_basis, neighbors, J, onsite)
    ranked = sorted(cand.items(), key=lambda kv: kv[1])
    best = ranked[0][0]
    if best == "fm":
        return 0, "总磁化 (铁磁, T_C)", False, cand, best
    if best == "neel":
        return 1, "交错 (反铁磁, T_N)", False, cand, best
    tag = "x" if best == "stripe_x" else "y"
    kind = (f"比热 Cv 峰 (auto: T=0 基态是条纹反铁磁 (沿 {tag} 方向), "
            f"本引擎的 ±1 交错序参量不适用 -> 用比热峰定 T_N)")
    return 0, kind, True, cand, best


def resolve_auto_order(model, size):
    """给别的脚本 (mc_fss) 用的入口: 先按 size 建数组, 再判定 auto 序参量。"""
    _spins, neighbors, J, onsite, num_basis = build_arrays(model, size, spin_seed=0)
    return _resolve_auto(neighbors, J, onsite, size, num_basis)


def _tau_int(series):
    """积分自相关时间 τ_int (单位: 一个测量步), 用初始正序列 (Sokal) 截断。

    定义 τ_int = 1/2 + Σ_{t=1..W} ρ(t): 纯白噪声 ρ(t)=0 (t>=1) ⇒ τ=0.5, 即每次采样
    都算一个独立样本; 强关联时 τ 大, 有效独立样本数 n_eff = n / (2*τ) 远小于 n。
    截断规则: 在 ρ(t) 首次降到 0 以下**且 t >= 5** 时停 (Sokal 的自动窗口; 只取首个
    负值而不留窗口会系统性低估 τ)。

    返回 (tau, window)。样本 < 16 或方差为 0 时退化为 (0.5, 0) —— 此时 n_eff = n。
    """
    x = np.asarray(series, dtype=float).ravel()
    n = x.size
    if n < 16:
        return 0.5, 0
    x = x - x.mean()
    var = float(np.dot(x, x) / n)
    if var <= 1e-300:
        return 0.5, 0
    tmax = min(n // 2, 20000)
    tau = 0.5
    window = 0
    for t in range(1, tmax + 1):
        rho = float(np.dot(x[:n - t], x[t:]) / (n - t) / var)
        if rho <= 0.0 and t >= 5:
            break
        tau += rho
        window = t
    return max(0.5, float(tau)), window


def _binning_sem(series, n_bins=None):
    """把序列切成 n_bins 块、用块均值的标准误估 SEM (自动吸收相关性)。

    比 naive SEM (std/sqrt(n)) 稳健: 有相关性时 naive 会**低估**误差 (序列越相关
    低估越多), 而块均值在块长 >> τ_int 时近似独立。返回 (sem, n_bins, block_len);
    样本太少 (< 16) 时返回 (nan, 0, 0)。
    """
    x = np.asarray(series, dtype=float).ravel()
    n = x.size
    if n < 16:
        return float("nan"), 0, 0
    nb = int(n_bins) if n_bins else max(2, min(32, n // 8))
    blk = n // nb
    if blk < 2 or nb < 2:
        return float("nan"), 0, 0
    means = x[:blk * nb].reshape(nb, blk).mean(axis=1)
    sem = float(means.std(ddof=1) / np.sqrt(nb))
    return sem, nb, blk


def _naive_sem(series):
    """不扣相关性的 SEM (std/sqrt(n)) —— 只用来给「低估了多少倍」提供对照。"""
    x = np.asarray(series, dtype=float).ravel()
    if x.size < 2:
        return float("nan")
    return float(x.std(ddof=1) / np.sqrt(x.size))


def _diag_row(series, accept, therm, meas):
    """单温度的采样诊断: 接受率 / τ_int / 有效样本 / 两种 SEM。

    `series` 是测量期的序参量 (或能量) 时间序列, 步长 = 一次测量 (= N 次单自旋尝试)。
    """
    tau, win = _tau_int(series)
    n = len(series)
    neff = n / (2.0 * tau) if tau > 0 else float(n)
    sem_b, nb, blk = _binning_sem(series)
    return {"accept": float(accept), "tau_int": float(tau), "window": int(win),
            "n_eff": float(neff), "sem_naive": _naive_sem(series),
            "sem_binned": float(sem_b), "n_bins": int(nb), "block": int(blk),
            "n_meas": int(n), "n_therm": int(therm), "n_total": int(meas)}


def _selftest():
    """比热能量的自检: _energy 必须等于暴力求和, 且**含单离子各向异性**。

    回归的缺陷 (P0-1): 旧 _energy 只算交换项、把 onsite 整个丢掉; 而 Metropolis 的
    dE 含各向异性 -> 抽样含 A、测量不含 A -> Cv 的能量涨落少一块 -> T_N 系统偏移,
    且退出码仍是 0 (静默错)。这里用非对角/非对称的随机 J 与随机易轴 A 把两种算法
    摆在一起比, 三项都过才算通过。
    """
    rng = np.random.default_rng(20250516)
    L = 3
    N = L * L
    M = 4
    spins = rng.normal(size=(N, 3))
    spins /= np.linalg.norm(spins, axis=1)[:, None]
    neighbors = np.full((N, M), -1, dtype=np.int32)
    Jarr = np.zeros((N, M, 3, 3))
    on = np.zeros((N, 3, 3))
    for cx in range(L):
        for cy in range(L):
            i = cx * L + cy
            for k, (dx, dy) in enumerate([(1, 0), (0, 1), (2, 0), (0, 2)]):
                neighbors[i, k] = ((cx + dx) % L) * L + ((cy + dy) % L)
                Jarr[i, k] = rng.normal(size=(3, 3)) * 1.5
            on[i] = np.diag(rng.normal(size=3) * 1.0)   # 随机易轴, 各格点不同

    fails = []
    e_fast = float(_energy(spins, neighbors, Jarr, on))
    e_ref = _brute_energy(spins, neighbors, Jarr, on)
    if abs(e_fast - e_ref) > 1e-9:
        fails.append(f"含 A 时不一致: _energy={e_fast:.12f} != 暴力 {e_ref:.12f}")
    else:
        print(f"  [PASS] 含 A: _energy 与暴力求和一致 (E/N={e_fast:.6f} meV/格点)")

    zero = np.zeros_like(on)
    e0_fast = float(_energy(spins, neighbors, Jarr, zero))
    e0_ref = _brute_energy(spins, neighbors, Jarr, zero)
    if abs(e0_fast - e0_ref) > 1e-9:
        fails.append(f"A=0 时不一致: _energy={e0_fast:.12f} != 暴力 {e0_ref:.12f}")
    else:
        print(f"  [PASS] A=0: _energy 与暴力求和一致 (E/N={e0_fast:.6f} meV/格点)")

    d_e = e_fast - e0_fast
    if abs(d_e) < 1e-6:
        fails.append("A 与 A=0 的能量完全相同 -> onsite 项仍没被算进去")
    else:
        print(f"  [PASS] onsite 项确实进入能量: A=0 -> {e0_fast:.6f}, "
              f"含 A -> {e_fast:.6f} (差 {d_e:+.6f} meV/格点)")

    if fails:
        for msg in fails:
            print(f"  [FAIL] {msg}")
        return 1

    # ---- P1-4: 采样诊断的估计量本身也要自检 (τ_int / binning SEM) ----
    rng2 = np.random.default_rng(7)
    phi = 0.8
    n = 20000
    noise = rng2.normal(size=n)
    ar = np.empty(n)
    ar[0] = noise[0]
    for i in range(1, n):
        ar[i] = phi * ar[i - 1] + noise[i]
    tau_ar, _ = _tau_int(ar)
    tau_expect = 0.5 + phi / (1.0 - phi)   # AR(1): ρ(t)=φ^t ⇒ τ = 1/2 + φ/(1-φ) = 4.5
    ok = abs(tau_ar - tau_expect) < 1.5
    print(("  [PASS] " if ok else "  [FAIL] ")
          + f"AR(1) φ={phi} 的 τ_int {tau_ar:.2f} ≈ 理论 {tau_expect:.2f}")
    if not ok:
        fails.append(f"τ_int 估计偏差过大: {tau_ar:.2f} vs 理论 {tau_expect:.2f}")

    tau_wn, _ = _tau_int(rng2.normal(size=n))
    ok = 0.5 <= tau_wn <= 1.2
    print(("  [PASS] " if ok else "  [FAIL] ")
          + f"白噪声的 τ_int {tau_wn:.2f} ≈ 0.5 (每次采样近似独立)")
    if not ok:
        fails.append(f"白噪声 τ_int 应 ≈0.5, 实得 {tau_wn:.2f}")

    sem_ar, _, _ = _binning_sem(ar)
    ratio_ar = sem_ar / max(_naive_sem(ar), 1e-300)
    ok = ratio_ar > 2.0
    print(("  [PASS] " if ok else "  [FAIL] ")
          + f"强相关序列: binning SEM / naive SEM = {ratio_ar:.2f} (>2 ⇒ naive 确实低估)")
    if not ok:
        fails.append(f"相关序列上 binning SEM 应远大于 naive, 实得 {ratio_ar:.2f}x")

    wn = rng2.normal(size=n)
    ratio_wn = _binning_sem(wn)[0] / max(_naive_sem(wn), 1e-300)
    ok = 0.5 <= ratio_wn <= 2.0
    print(("  [PASS] " if ok else "  [FAIL] ")
          + f"白噪声: binning SEM / naive SEM = {ratio_wn:.2f} (≈1 ⇒ 两者一致)")
    if not ok:
        fails.append(f"白噪声上两种 SEM 应一致, 实得 {ratio_wn:.2f}x")

    if fails:
        for msg in fails:
            print(f"  [FAIL] {msg}")
        return 1
    print("mc_2d 自检: 比热能量 (交换 + 单离子各向异性) + 采样诊断 (τ_int/binning SEM) OK")
    return 0


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
    acc = 0
    e_series = np.empty(n_meas, dtype=np.float64)
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
                acc += 1
        if it >= n_therm:
            e = _energy(spins, neighbors, J, onsite)
            e1 += e
            e2 += e * e
            e_series[it - n_therm] = e
    mean_e = e1 / n_meas
    cv = (e2 / n_meas - mean_e * mean_e) * N / (0.08617 * T * T)
    return mean_e, cv, acc / float((n_therm + n_meas) * N), e_series


def run_sweep(spins, neighbors, J, onsite, T_list, therm, meas, stag, order, seed=0,
              diag=None):
    """逐温度扫一遍, 返回 (nt, 4) 数组 [T, m, chi, U4]。

    `diag` 给一个 dict 时会**顺手填**采样诊断 (P1-4: 接受率 / τ_int / 有效样本 / SEM),
    键为 "T"/"accept"/"tau_int"/"n_eff"/"sem_naive"/"sem_binned"/"window"/"n_bins",
    每项都是按 T_list 顺序的 list。不传就完全不算 (对 mc_fss 这类老调用方零开销)。
    """
    state = np.array([np.uint64(seed if seed != 0 else 88172645463325252)], dtype=np.uint64)
    results = []
    for T in T_list:
        mean_m, chi, u4, acc, mser = _run_T_fast(spins, neighbors, J, onsite, float(T), therm,
                                                  meas, state, order, stag)
        results.append([T, mean_m, chi, u4])
        if diag is not None:
            row = _diag_row(mser, acc, therm, meas)
            row["T"] = float(T)
            for key, val in row.items():
                diag.setdefault(key, []).append(val)
        print(f"  T={T:6.1f} K | m={mean_m:.4f} | chi={chi:.4f} | U4={u4:.4f}", flush=True)
    return np.array(results)


def run_sweep_cv(spins, neighbors, J, onsite, T_list, therm, meas, seed=0, diag=None):
    """比热路径的逐温度扫描, 返回 (nt, 3) 数组 [T, E/N, Cv]; `diag` 语义同 run_sweep。"""
    state = np.array([np.uint64(seed if seed != 0 else 88172645463325252)], dtype=np.uint64)
    results = []
    for T in T_list:
        mean_e, cv, acc, eser = _run_T_cv(spins, neighbors, J, onsite, float(T), therm, meas, state)
        results.append([T, mean_e, cv])
        if diag is not None:
            row = _diag_row(eser, acc, therm, meas)
            row["T"] = float(T)
            for key, val in row.items():
                diag.setdefault(key, []).append(val)
        print(f"  T={T:6.1f} K | E/N={mean_e:.4f} | Cv={cv:.4f}", flush=True)
    return np.array(results)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("--selftest", action="store_true",
                    help="只跑比热能量自检 (交换 + 单离子各向异性), 不需要 model")
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
    # --selftest 不读模型: 预扫描一次, 免得 argparse 因为缺 model 位置参数而报错
    if "--selftest" in sys.argv[1:]:
        return _selftest()
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
            # 真正的判据要等 build_arrays 之后: 在候选基态里比 T=0 能量
            # (见 _candidate_orders)。旧实现只拿 J1 的符号猜, 于是 J1>0 + J2<0 的
            # 条纹反铁磁会被当成 Néel —— 测出来的峰不是 T_N, 但模型照常打印一个温度。
            order = 1
            kind = "auto (待定: 依 T=0 候选基态能量)"
        elif o_kind == "sublattice":
            order = 1
            kind = (f"逐子晶格 Néel 矢量 |L| (A = 子晶格 {sorted(set(o_mask))}, "
                    f"B = 其余; 交变磁体序参量, P5)")
        else:
            order = 1 if o_kind == "staggered" else 0
            kind = "交错 (反铁磁, T_N)" if order == 1 else "总磁化 (铁磁, T_C)"
    if o_kind != "auto":
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

    if not use_cv and o_kind == "auto":
        order, kind, use_cv, cand, best = _resolve_auto(
            neighbors, J, onsite, args.size, num_basis)
        ranked = sorted(cand.items(), key=lambda kv: kv[1])
        e_best = ranked[0][1]
        print("auto 序参量判定: T=0 候选基态能量 (meV/格点) = "
              + " | ".join(f"{k}={v:+.4f}" for k, v in ranked))
        if len(ranked) > 1 and abs(ranked[1][1] - e_best) < 1e-6:
            print(f">>> [warn] 候选基态近简并 ({best} 与 {ranked[1][0]} 只差 "
                  f"{abs(ranked[1][1] - e_best):.2e} meV/格点) —— auto 序参量不可靠, "
                  f"请显式给 --order (staggered/total/cv/sublattice:<mask>)",
                  file=sys.stderr)
        if use_cv:
            # 条纹 (或别的非 fm/neel) 基态: ±1 交错序参量测的不是这个序, 改用比热峰。
            print(f">>> [warn] auto 判定的基态是 {best} (E={e_best:+.4f} meV/格点, "
                  f"低于 fm 的 {cand['fm']:+.4f}) —— 条纹反铁磁不能用 ±1 交错序参量测峰, "
                  f"已自动改用 --order cv (比热峰); 需要逐子晶格序参量请显式 "
                  f"--order sublattice:<mask>", file=sys.stderr)
        print(f"序参量: {kind}")

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
    diags = []
    peaks = {"chi": [], "dm": [], "cv": []}
    for k in range(n_seeds):
        spins, _, _, _, _ = build_arrays(model, args.size, spin_seed=args.seed + k)
        sd = args.seed + k
        dk = {}
        if use_cv:
            res = run_sweep_cv(spins, neighbors, J, onsite, T_list, args.therm, args.meas,
                               seed=sd, diag=dk)
        else:
            res = run_sweep(spins, neighbors, J, onsite, T_list, args.therm, args.meas,
                            stag, order, seed=sd, diag=dk)
        diags.append(dk)
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

    # ---- 采样诊断 (P1-4): 接受率 / 积分自相关时间 / 有效独立样本 / 两种 SEM ----
    # 没有这些量时, 「峰位差 1 K」「种子间 ±0.5 K」都无法判断是物理还是统计噪声。
    warns = []
    diag = {}
    if diags and diags[0]:
        for _key in diags[0]:
            diag[_key] = np.mean(np.array([d[_key] for d in diags], dtype=float), axis=0)
        Td = diag["T"]; dacc = diag["accept"]; dtau = diag["tau_int"]; dneff = diag["n_eff"]
        sn, sb = diag["sem_naive"], diag["sem_binned"]
        i_pk = int(np.argmax(res[:, 2]))
        # 峰落在窗口边界时, 上面已报「未分辨」; 再拿该点的 τ_int/接受率说事只是噪声
        # (端点温度往往太低, 自旋本来就冻住) —— 那一行诊断会被用户当成「再调参就能救」。
        pk_edge = i_pk in (0, len(res) - 1)
        ratio = (float(np.nanmedian(sb)) / max(float(np.nanmedian(sn)), 1e-300)
                 if np.isfinite(np.nanmedian(sn)) and np.isfinite(np.nanmedian(sb)) else float("nan"))
        print(f"\n>>> 采样诊断 (每温度 {args.meas} 次测量 x {N} 次单自旋尝试/步):")
        print(f"    接受率 {float(np.min(dacc)):.2f}..{float(np.max(dacc)):.2f} | "
              f"τ_int 最大 {float(np.max(dtau)):.2f} 步 (T={Td[int(np.argmax(dtau))]:.1f} K) | "
              f"有效独立样本 n_eff 最小 {float(np.min(dneff)):.0f}/{args.meas}")
        print(f"    SEM: naive 中位 {float(np.nanmedian(sn)):.3g} vs binning 中位 "
              f"{float(np.nanmedian(sb)):.3g} (naive 低估 {ratio:.2f}x)" if np.isfinite(ratio)
              else "    SEM: 样本不足, 未估计")
        if pk_edge:
            print(f"    (峰在窗口边界 T={Td[i_pk]:.2f} K, 该点的 τ_int/接受率不单独告警 —— "
                  f"先解决「峰未分辨」再谈采样量)")
        if dneff[i_pk] < 100 and not pk_edge:
            warns.append(f"峰处 T={Td[i_pk]:.2f} K 的有效独立样本只有 {dneff[i_pk]:.0f} < 100 "
                         f"(τ_int={dtau[i_pk]:.1f} 步): 峰位/峰高的统计误差大 —— 把 --meas 提到 "
                         f">= {int(200 * (1.0 + 2.0 * dtau[i_pk]))}, 或增加 --seeds 做多种子平均")
        if dacc[i_pk] < 0.05 and not pk_edge:
            warns.append(f"峰处 T={Td[i_pk]:.2f} K 的 Metropolis 接受率只有 {dacc[i_pk]:.2f}: "
                         f"自旋基本被冻在原地, 该点的序参量/比热不可信 (检查 --size/--tmin)")
        if np.isfinite(ratio) and ratio > 3.0:
            warns.append(f"naive SEM 比 binning SEM 小 {ratio:.1f}x: 序列强相关, "
                         f"误差棒必须用 binning 那个 (n_eff={float(np.min(dneff)):.0f})")
        for _w in warns:
            print(f">>> [warn] {_w}", file=sys.stderr)

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

    # ---- 诊断落盘: 机器可读 (供 driver/report 引用, 也让复核时不必重跑) ----
    if diag:
        dj = {"schema": "mc_diag/1", "script": "mc_2d.py", "order": str(args.order),
              "use_cv": bool(use_cv), "size": [int(x) for x in args.size],
              "n_mag": int(N), "n_therm": int(args.therm), "n_meas": int(args.meas),
              "seeds": int(n_seeds),
              "note": "tau_int/n_eff 按测量步计 (1 步 = N 次单自旋尝试); "
                      "误差棒请用 sem_binned (naive 在相关序列上会低估)",
              "warnings": warns, "columns": ["T", "accept", "tau_int", "n_eff",
                                             "sem_naive", "sem_binned", "window", "n_bins"]}
        for _k in dj["columns"]:
            dj[_k] = [float(x) for x in np.atleast_1d(diag[_k])]
        fn_diag = f"{prefix}_diag.json"
        with open(fn_diag, "w", encoding="utf-8") as f:
            json.dump(dj, f, ensure_ascii=False, indent=1)
        print(f">>> 采样诊断存至 {fn_diag}")


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
    # sys.exit 而不是裸 main(): 让 --selftest 的失败能变成非零退出码
    sys.exit(main())
