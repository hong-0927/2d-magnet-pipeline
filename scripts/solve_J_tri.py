#!/usr/bin/env python
"""solve_J_tri.py
从三角晶格各**共线**磁序总能量反解 J1/J2/J3 (三角 J1-J2-J3 模型)。

读 02_supercell/cases/magnetic_orders.json 里的 coefficients_per_cell + inversion_matrix_per_cell,
乘以各构型能量向量即可。约定 H = -J Σ S·S (J>0 = 铁磁), meV。

默认四态是 {fm, stripy_a, stripy_ab, uud} (全共线, 5x5x1 实测 rank=4, cond≈3.7e3)。
旧的 {fm, 120deg, stripy_a, stripy_ab} (cond≈5.7e3) 需要 LNONCOLLINEAR; 换成共线四态
后设计矩阵仍然满秩、条件数更好, 所以 120° 不再是反解的前提 (见 build_magnetic_cases_tri.py)。
反解本身只认 json 里的 selected_orders —— 换序后 builder 会连 cond 一起重算。

用法:
  python solve_J_tri.py <cases_dir>          # 读 cases/<order>/OUTCAR
  python solve_J_tri.py <cases_dir> --demo   # 合成 J 自检
"""
import argparse
import json
import os

import numpy as np


def read_toten(path):
    """取 OUTCAR 的总能量 (eV)。

    P3.17: 口径统一到 `dft_energy.py` —— 一律 `free energy TOTEN`、取**最后一次**出现、
    缺了就报错（不静默改用 energy(sigma->0)）。这里只转发, 不再各写一份。
    """
    from dft_energy import read_toten as _rt
    return _rt(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases")
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args()

    meta = json.load(open(os.path.join(args.cases, "magnetic_orders.json"), encoding="utf-8"))
    sel = meta["selected_orders"]
    S2 = meta["S2"]
    c = meta["coefficients_per_cell"]
    Ainv = np.array(meta["inversion_matrix_per_cell"], dtype=float)
    A = np.array([[1.0] + c[k] for k in sel])
    if len(sel) != 4:
        raise SystemExit(f"selected_orders 必须是 4 个构型 (e0+3 条 J), 实际 {sel}")
    spin_mode = meta.get("spin_mode") or ("noncollinear" if meta.get("noncollinear") else "collinear")
    rk = int(np.linalg.matrix_rank(A))
    cond = float(np.linalg.cond(A))
    if rk != 4 or cond > 1e4:
        raise SystemExit(f"设计矩阵不可解: rank={rk}/4, cond={cond:.4g} (>1e4); "
                         f"selected={sel} —— 请用 builder 重新生成 (它会重算 cond)。")

    if args.demo:
        J_true = [-0.3, 0.05, -0.02]
        Jt = [j * S2 for j in J_true]
        E = {k: sum(c[k][m] * Jt[m] for m in range(3)) for k in sel}
        b = np.array([E[k] for k in sel]) / 1000.0  # meV -> eV
        print("[demo] 合成 J (meV):", J_true, "-> 磁态能量 (meV):", [f"{E[k]:.4f}" for k in sel])
    else:
        b = []
        missing = []
        for k in sel:
            p = os.path.join(args.cases, k, "OUTCAR")
            if os.path.exists(p):
                b.append(read_toten(p))
            else:
                b.append(float("nan"))
                missing.append(k)
        if missing:
            raise SystemExit(f"缺少 {missing} 的 OUTCAR; 可加 --demo 自检")
        b = np.array(b)

    sol = Ainv @ b  # [e0, J1~, J2~, J3~] in eV
    J = [sol[1 + m] / S2 * 1000.0 for m in range(3)]
    print(f"\n构型顺序: {sel}   自旋模式: {spin_mode}")
    print(f"系数矩阵 (行=构型):\n{A}")
    print(f"可解性: rank={rk}/4, cond(A)={cond:.4g}")
    print(f"能量 (eV/胞): {[f'{x:+.6f}' for x in b]}")
    print(f"e0 = {sol[0]:+.6f} eV")
    print(f"J1 = {J[0]:+.4f} meV   J2 = {J[1]:+.4f} meV   J3 = {J[2]:+.4f} meV   (S^2={S2})")

    if args.demo:
        ok = all(abs(g - t) < 0.02 for g, t in zip(J, J_true))
        print(f"期望 (meV): {J_true} -> " + ("PASS" if ok else "FAIL"))


if __name__ == "__main__":
    main()
