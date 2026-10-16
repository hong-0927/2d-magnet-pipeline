#!/usr/bin/env python
"""solve_J3.py
从四个磁态 (FM/Néel/zigzag/stripy) 总能量反解 J1/J2/J3 (蜂窝 J1-J2-J3 模型)。

相对 CrI3/CrBr3 流水线的 solve_J.py (J1/J2) 扩展了第三近邻 J3。MnPS3 的 J3 不可忽略
(Sivadas PRB 91,235425: J3 > J2), 且 zigzag/stripy 两 M 点磁序在 J3 下不再简并, 所以
四个磁态正好解 {e0, J1, J2, J3}。

用法:
  python solve_J3.py <cases_dir>                # 读 cases/{fm,neel,zigzag,stripy}/OUTCAR
  python solve_J3.py <cases_dir> --demo          # 用文献 J 值生成的假能量自检

约定: H = -J1 sum_<ij> S.S - J2 sum_<<ij>> S.S - J3 sum_<<<ij>>> S.S, J>0 = 铁磁, meV。
"""
import argparse
import json
import os

import numpy as np

# 文献 J (Sivadas PRB 91,235425, 表 I, 单层 MnPS3), 转成 H=-J S.S 约定 (J>0 铁磁):
LIT = dict(J1=-0.79, J2=-0.04, J3=-0.23)


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

    meta = json.load(open(os.path.join(args.cases, "magnetic_orders.json")))
    n = meta["n_mag"]
    S2 = meta["S2"]
    c = meta["coefficients_per_cell"]  # 每胞系数 (8 Mn 超胞)

    A = np.array([[1.0] + c["fm"],
                  [1.0] + c["neel"],
                  [1.0] + c["zigzag"],
                  [1.0] + c["stripy"]])

    if args.demo:
        # 用文献 J 生成四个磁态能量 (每胞, meV), 自检反解
        Jt = np.array([LIT["J1"], LIT["J2"], LIT["J3"]]) * S2  # J~ = J*S^2
        E = dict(fm=0.0)
        for nm in ["fm", "neel", "zigzag", "stripy"]:
            E[nm] = float(c[nm][0] * Jt[0] + c[nm][1] * Jt[1] + c[nm][2] * Jt[2])
        print("[demo] 用文献 J 生成的磁态能量 (eV):",
              {k: f"{v:.6f}" for k, v in E.items()})
        E = {k: v / 1000.0 for k, v in E.items()}  # meV -> eV (与 OUTCAR 单位一致)
    else:
        E = {}
        for nm in ["fm", "neel", "zigzag", "stripy"]:
            p = os.path.join(args.cases, nm, "OUTCAR")
            if os.path.exists(p):
                E[nm] = read_toten(p)
        if set(["fm", "neel", "zigzag", "stripy"]) != set(E):
            raise SystemExit("需要 fm/neel/zigzag/stripy 四态 OUTCAR; 可加 --demo 自检")

    # 每胞能量 -> 反解 [e0, J1~, J2~, J3~] (eV)
    b = np.array([E["fm"], E["neel"], E["zigzag"], E["stripy"]])
    e0, J1t, J2t, J3t = np.linalg.solve(A, b)

    J1 = J1t / S2 * 1000.0
    J2 = J2t / S2 * 1000.0
    J3 = J3t / S2 * 1000.0

    print(f"\n每胞系数 (8 Mn): fm{c['fm']}  neel{c['neel']}  "
          f"zigzag{c['zigzag']}  stripy{c['stripy']}")
    print(f"能量 (eV/胞): fm={E['fm']:+.6f}  neel={E['neel']:+.6f}  "
          f"zigzag={E['zigzag']:+.6f}  stripy={E['stripy']:+.6f}")
    print(f"J1~ = {J1t:+.6f} eV   J2~ = {J2t:+.6f} eV   J3~ = {J3t:+.6f} eV   "
          f"(S^2 = {S2})")
    print(f"\nJ1 = {J1:+.4f} meV    J2 = {J2:+.4f} meV    J3 = {J3:+.4f} meV")

    if args.demo:
        ok = (abs(J1 - LIT["J1"]) < 0.02 and abs(J2 - LIT["J2"]) < 0.02
              and abs(J3 - LIT["J3"]) < 0.02)
        print(f"\n期望 (文献): J1={LIT['J1']}, J2={LIT['J2']}, J3={LIT['J3']} meV  ->  "
              + ("PASS" if ok else "FAIL"))


if __name__ == "__main__":
    main()
