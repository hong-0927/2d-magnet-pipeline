#!/usr/bin/env python
"""solve_J_rect.py
从四个共线磁态 (fm/afm_a/afm_b/afm_diag) 总能量反解矩形晶格 J1/J2/J3。

矩形晶格 (正交, a != b), 距离壳按升序: J1 = a 轴, J2 = b 轴, J3 = 对角 (|a+b|)。
四个共线磁态线性无关, 恰好解 {e0, J1, J2, J3}。

用法:
  python solve_J_rect.py <cases_dir>              # 读 cases/{fm,afm_a,afm_b,afm_diag}/OUTCAR
  python solve_J_rect.py <cases_dir> --demo        # 用 CrSBr 文献 J 生成的假能量自检

约定: H = -J1 Σ S·S - J2 Σ S·S - J3 Σ S·S, J>0 = 铁磁, meV。
参考: CrSBr J1(a)=2.79, J2(b)=2.12, J3(对角)=3.88 meV (Ke Yang & Hua Wu, PRB 104, 144416 (2021),
      按其论文命名 J1=对角 3.88 / J2=a 2.79 / J3=b 2.12, 此处按距离壳升序重排)。
"""
import argparse
import json
import os

import numpy as np

LIT = dict(J1=2.79, J2=2.12, J3=3.88)  # CrSBr 按距离排, 全铁磁


def read_toten(path):
    tot = None
    with open(path) as f:
        for line in f:
            if "TOTEN" in line and "free" in line:
                tot = float(line.split()[-2])
    if tot is None:
        raise ValueError(f"no 'free energy TOTEN' in {path}")
    return tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases")
    ap.add_argument("--demo", action="store_true")
    args = ap.parse_args()

    meta = json.load(open(os.path.join(args.cases, "magnetic_orders.json")))
    n = meta["n_mag"]
    S2 = meta["S2"]
    c = meta["coefficients_per_cell"]

    orders = ["fm", "afm_a", "afm_b", "afm_diag"]
    A = np.array([[1.0] + c[k] for k in orders])

    if args.demo:
        Jt = np.array([LIT["J1"], LIT["J2"], LIT["J3"]]) * S2  # J~ = J*S^2
        E = {k: float(c[k][0] * Jt[0] + c[k][1] * Jt[1] + c[k][2] * Jt[2]) for k in orders}
        print("[demo] 用 CrSBr 文献 J 生成的磁态能量 (eV):",
              {k: f"{v:.6f}" for k, v in E.items()})
        E = {k: v / 1000.0 for k, v in E.items()}  # meV -> eV
    else:
        E = {}
        for nm in orders:
            p = os.path.join(args.cases, nm, "OUTCAR")
            if os.path.exists(p):
                E[nm] = read_toten(p)
        if set(orders) != set(E):
            raise SystemExit("需要 fm/afm_a/afm_b/afm_diag 四态 OUTCAR; 可加 --demo 自检")

    b = np.array([E[k] for k in orders])
    e0, J1t, J2t, J3t = np.linalg.solve(A, b)

    J1 = J1t / S2 * 1000.0
    J2 = J2t / S2 * 1000.0
    J3 = J3t / S2 * 1000.0

    print(f"\n每胞系数: " + "  ".join(f"{k}{c[k]}" for k in orders))
    print(f"能量 (eV/胞): " + "  ".join(f"{k}={E[k]:+.6f}" for k in orders))
    print(f"J1~ = {J1t:+.6f} eV   J2~ = {J2t:+.6f} eV   J3~ = {J3t:+.6f} eV   (S^2 = {S2})")
    print(f"\nJ1 = {J1:+.4f} meV    J2 = {J2:+.4f} meV    J3 = {J3:+.4f} meV")
    print(f"   (J1=a 轴近邻, J2=b 轴近邻, J3=对角)")

    if args.demo:
        ok = (abs(J1 - LIT["J1"]) < 0.02 and abs(J2 - LIT["J2"]) < 0.02
              and abs(J3 - LIT["J3"]) < 0.02)
        print(f"\n期望 (CrSBr): J1={LIT['J1']}, J2={LIT['J2']}, J3={LIT['J3']} meV  ->  "
              + ("PASS" if ok else "FAIL"))


if __name__ == "__main__":
    main()
