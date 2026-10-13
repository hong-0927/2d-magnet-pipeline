#!/usr/bin/env python
"""mc_2d.py
二维蜂窝磁体 经典 Heisenberg Monte Carlo, 计算序参量 vs 温度 -> 临界温度 (T_C / T_N)。

自包含: 只依赖 numpy + numba + matplotlib, 不依赖 spin_mc 共享引擎。
  - 自写 nopython 快速 kernel (xorshift64 RNG)。
  - 支持**总磁化** (铁磁, 算 T_C) 与**交错磁化** (反铁磁, 算 T_N) 两种序参量。
  - 序参量路径额外算 **Binder 累积量 U4** (第 4 列): 有序 -> 0.4, 高温无序 -> 0;
    不同 --size 的 U4(T) 交点是有限尺寸标度定 T_C 的稳健判据。
  - --seeds N: 跑 N 个独立种子 (不同初态/MC 状态) 并平均, 报告峰位标准差 (误差棒)。

用法:
  python mc_2d.py model.json --order auto --size 30 30 --tmin 5 --tmax 150 --nt 30 \
      --therm 2500 --meas 6000 --seeds 5 --out out_prefix
  --order auto: J1<0 -> 交错 (反铁磁), J1>0 -> 总磁化 (铁磁)
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


@njit
def _xorshift64(state):
    x = state[0]
    x ^= x >> np.uint64(12)
    x ^= x << np.uint64(25)
    x ^= x >> np.uint64(27)
    state[0] = x
    return float(x * np.uint64(2685821657736338717)) * 5.421010862427522e-20


@njit
def _run_T_fast(spins, neighbors, J, onsite, T, n_therm, n_meas, num_basis, state, order):
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
            if order == 1:  # staggered
                max_ = 0.0; may = 0.0; maz = 0.0
                mbx = 0.0; mby = 0.0; mbz = 0.0
                ca = 0; cb = 0
                for i in range(N):
                    if i % num_basis == 0:
                        max_ += spins[i, 0]; may += spins[i, 1]; maz += spins[i, 2]; ca += 1
                    else:
                        mbx += spins[i, 0]; mby += spins[i, 1]; mbz += spins[i, 2]; cb += 1
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


def run_sweep(spins, neighbors, J, onsite, T_list, therm, meas, num_basis, order, seed=0):
    state = np.array([np.uint64(seed if seed != 0 else 88172645463325252)], dtype=np.uint64)
    results = []
    for T in T_list:
        mean_m, chi, u4 = _run_T_fast(spins, neighbors, J, onsite, float(T), therm, meas,
                                      num_basis, state, order)
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
    ap.add_argument("--order", choices=["auto", "staggered", "total", "cv"], default="auto")
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

    if args.order == "cv":
        use_cv = True
        order = 0
        kind = "比热 Cv 峰 (受挫/三角晶格, T_N)"
    else:
        use_cv = False
        if args.order == "auto":
            # 单壳层模型给 J1_meV; 多壳层模型 (FGT 双亚晶格) 没有这个键 ->
            # 用各壳层 J 之和的符号判断铁磁/反铁磁, 不能硬取 meta["J1_meV"]
            if "J1_meV" in meta:
                jsum = float(meta["J1_meV"])
            else:
                jsum = sum(-float(it["J"][0][0]) for it in model.get("interactions", []))
            order = 1 if jsum < 0 else 0
            kind = "交错 (反铁磁, T_N)" if order == 1 else "总磁化 (铁磁, T_C)"
        else:
            order = 1 if args.order == "staggered" else 0
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
                            num_basis, order, seed=sd)
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
            ax1.set_title(f"{meta['element']} honeycomb MC  ({t_lbl})")
            # Binder 累积量子图 (有限尺寸标度判据): 有序 -> 0.4, 无序 -> 0
            ax3.plot(res[:, 0], res[:, 3], "^-", color="tab:green", label="Binder $U_4$")
            ax3.axhline(0.4, color="gray", ls=":", lw=0.8)
            ax3.set_xlabel("T (K)")
            ax3.set_ylabel("Binder $U_4$")
            ax3.grid(alpha=0.3)
        ax1.grid(alpha=0.3)
        png = f"{prefix}_mt.png"
        fig.savefig(png, dpi=150, bbox_inches="tight")
        print(f">>> 图存至 {png}")
    except Exception as e:
        print(f"[warn] 画图跳过: {e}")


if __name__ == "__main__":
    main()
