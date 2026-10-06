#!/usr/bin/env python
"""mc_2d.py
二维蜂窝磁体 经典 Heisenberg Monte Carlo, 计算序参量 vs 温度 -> 临界温度 (T_C / T_N)。

自包含: 只依赖 numpy + numba + matplotlib, 不依赖 spin_mc 共享引擎。
  - 自写 nopython 快速 kernel (xorshift64 RNG)。
  - 支持**总磁化** (铁磁, 算 T_C) 与**交错磁化** (反铁磁, 算 T_N) 两种序参量。

用法:
  python mc_2d.py model.json --order auto --size 30 30 --tmin 5 --tmax 150 --nt 30 \
      --therm 2500 --meas 6000 --out out_prefix
  --order auto: J1<0 -> 交错 (反铁磁), J1>0 -> 总磁化 (铁磁)
"""

import argparse
import json
import math
import os

import numpy as np
from numba import njit


# --------------------------------------------------------------------------- lattice
def build_arrays(model, size):
    """从 poscar_to_mc3.py 的模型 JSON 直接构建 numba 所需数组。"""
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
    rng = np.random.default_rng(0)
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
    mean_m = s1 / n_meas
    var = s2 / n_meas - mean_m * mean_m
    chi = var * N / (0.08617 * T)
    return mean_m, chi


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
        mean_m, chi = _run_T_fast(spins, neighbors, J, onsite, float(T), therm, meas,
                                  num_basis, state, order)
        results.append([T, mean_m, chi])
        print(f"  T={T:6.1f} K | m={mean_m:.4f} | chi={chi:.4f}", flush=True)
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
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    with open(args.model) as f:
        model = json.load(f)
    meta = model["meta"]

    if args.order == "cv":
        use_cv = True
        order = 0
        kind = "比热 Cv 峰 (受挫/三角晶格, T_N)"
    else:
        use_cv = False
        if args.order == "auto":
            order = 1 if meta["J1_meV"] < 0 else 0
            kind = "交错 (反铁磁, T_N)" if order == 1 else "总磁化 (铁磁, T_C)"
        else:
            order = 1 if args.order == "staggered" else 0
            kind = "交错 (反铁磁, T_N)" if order == 1 else "总磁化 (铁磁, T_C)"
    print(f"序参量: {kind}")

    prefix = args.out or os.path.splitext(os.path.basename(args.model))[0]
    spins, neighbors, J, onsite, num_basis = build_arrays(model, args.size)
    N = spins.shape[0]
    print(f"格点: {N} (basis {num_basis} x {tuple(args.size)}) | "
          f"J1={meta['J1_meV']} J2={meta['J2_meV']} J3={meta.get('J3_meV', 0.0)} "
          f"MAE={meta['MAE_meV']} meV | S={meta['S']}")

    T_list = np.linspace(args.tmin, args.tmax, args.nt)
    if args.cool:
        T_list = T_list[::-1]
    if use_cv:
        res = run_sweep_cv(spins, neighbors, J, onsite, T_list, args.therm, args.meas)
    else:
        res = run_sweep(spins, neighbors, J, onsite, T_list, args.therm, args.meas,
                        num_basis, order)
    res = res[np.argsort(res[:, 0])]

    if use_cv:
        fn = f"{prefix}_mt.csv"
        np.savetxt(fn, res, delimiter=",", header="T_Kelvin,E_per_site,Cv", comments="")
        print(f"\n>>> 数据存至 {fn}")
        cv_peak = res[np.argmax(res[:, 2]), 0]
        print(f">>> 比热 Cv 峰处 T_N ≈ {cv_peak:.1f} K")
    else:
        fn = f"{prefix}_mt.csv"
        np.savetxt(fn, res, delimiter=",", header="T_Kelvin,order_param,Susceptibility",
                   comments="")
        print(f"\n>>> 数据存至 {fn}")
        chi_peak = res[np.argmax(res[:, 2]), 0]
        dm = -np.gradient(res[:, 1], res[:, 0])
        tn_dm = res[np.argmax(dm), 0]
        print(f">>> 磁化率峰处 T ≈ {chi_peak:.1f} K")
        print(f">>> |dm/dT| 最大处 T ≈ {tn_dm:.1f} K")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax1 = plt.subplots(figsize=(6.5, 4.5))
        if use_cv:
            ax1.plot(res[:, 0], res[:, 1], "o-", color="tab:blue", label="E/N (meV)")
            ax1.set_xlabel("T (K)")
            ax1.set_ylabel("energy per site", color="tab:blue")
            ax1.tick_params(axis="y", labelcolor="tab:blue")
            ax2 = ax1.twinx()
            ax2.plot(res[:, 0], res[:, 2], "s--", color="tab:red", alpha=0.7, label="Cv")
            ax2.set_ylabel("specific heat $C_v$", color="tab:red")
            ax2.tick_params(axis="y", labelcolor="tab:red")
            ax1.set_title(f"{meta['element']} MC (Cv)  ($T_N\\approx${res[np.argmax(res[:, 2]), 0]:.0f} K)")
        else:
            ylab = "staggered $m_{stag}$" if order == 1 else "magnetization $M$"
            tc_label = "T_N" if order == 1 else "T_C"
            ax1.plot(res[:, 0], res[:, 1], "o-", color="tab:blue", label=ylab)
            ax1.set_xlabel("T (K)")
            ax1.set_ylabel(ylab, color="tab:blue")
            ax1.tick_params(axis="y", labelcolor="tab:blue")
            ax2 = ax1.twinx()
            ax2.plot(res[:, 0], res[:, 2], "s--", color="tab:red", alpha=0.7,
                     label="susceptibility")
            ax2.set_ylabel("susceptibility $\\chi$", color="tab:red")
            ax2.tick_params(axis="y", labelcolor="tab:red")
            chi_peak = res[np.argmax(res[:, 2]), 0]
            ax1.set_title(f"{meta['element']} honeycomb MC  (${{{tc_label}}}\\approx${chi_peak:.0f} K)")
        ax1.grid(alpha=0.3)
        png = f"{prefix}_mt.png"
        fig.savefig(png, dpi=150, bbox_inches="tight")
        print(f">>> 图存至 {png}")
    except Exception as e:
        print(f"[warn] 画图跳过: {e}")


if __name__ == "__main__":
    main()
