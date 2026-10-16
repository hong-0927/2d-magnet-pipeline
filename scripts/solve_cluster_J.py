#!/usr/bin/env python
"""solve_cluster_J.py
从簇磁体各磁序总能量反解 {E0, J_intra, J1, J2, J3}（簇磁体能量映射模型）。

模型 (与 build_cluster_cases.py 一致, J>0=铁磁):
    E = E0 + J_intra*C_intra + J1*C1 + J2*C2 + J3*C3
其中 C_intra 由簇内 Nb-Nb 键和给出, C_k 由第 k 壳簇间净矩对给出。

读 02_supercell/cases/magnetic_orders.json 里 build_cluster_cases.py 存好的
selected_configs + pseudo_inverse_E0_Jintra_Jk, 乘以能量向量即可。

用法:
  python solve_cluster_J.py <cases_dir>          # 读 cases/<name>/OUTCAR
  python solve_cluster_J.py <cases_dir> --demo   # 用合成 J 生成能量自检

单位: 反解出的 J 为 eV/超胞, 打印时换算 meV。
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
    sel = meta["selected_configs"]
    coeffs = meta["coefficients"]
    pinv = np.array(meta["pseudo_inverse_E0_Jintra_Jk"], dtype=float)
    nbr = meta.get("nbr_shells", len(meta.get("inter_shells_angstrom", [])))

    # 设计矩阵行: [1, C_intra, C_inter...], 顺序与 sel 一致
    A = np.array([[1.0, coeffs[c]["C_intra"]] + list(coeffs[c]["C_inter"]) for c in sel])

    if args.demo:
        J_true = [-30.0] + [-2.0 - 0.3 * k for k in range(nbr)]  # J_intra + J1..Jk (meV)
        E0_true = 0.0
        Evec = np.array([E0_true + A[i, 1] * J_true[0]
                         + sum(A[i, 2 + k] * J_true[1 + k] for k in range(nbr))
                         for i in range(len(sel))])
        b = Evec / 1000.0  # meV -> eV
        print("[demo] 合成 J_intra/J:", [f"{x:g}" for x in J_true],
              "-> 能量 (meV/胞):", [f"{v:g}" for v in Evec])
    else:
        b = []
        missing = []
        for c in sel:
            p = os.path.join(args.cases, c, "OUTCAR")
            if os.path.exists(p):
                b.append(read_toten(p))
            else:
                b.append(float("nan"))
                missing.append(c)
        if missing:
            raise SystemExit(f"缺少 {missing} 的 OUTCAR; 可加 --demo 自检")
        b = np.array(b)

    sol = pinv @ b  # [E0, J_intra, J1..Jk] in eV
    names = ["E0", "J_intra"] + [f"J{k + 1}" for k in range(nbr)]
    print(f"\n构型顺序: {sel}")
    print(f"系数矩阵 (行=构型):\n{A}")
    print("\n反解结果:")
    for nm, v in zip(names, sol):
        print(f"  {nm:8s} = {v:+.6f} eV   = {v * 1000.0:+.4f} meV")

    if args.demo:
        exp = [E0_true] + J_true  # meV
        got = (sol * 1000.0).tolist()
        ok = all(abs(g - e) < 1e-6 for g, e in zip(got, exp))
        print(f"\n期望 (meV): {[f'{x:g}' for x in exp]}")
        print("-> " + ("PASS" if ok else "FAIL"))


if __name__ == "__main__":
    main()
