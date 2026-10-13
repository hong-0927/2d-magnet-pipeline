#!/usr/bin/env python
"""solve_fgt_J.py — Fe3GeTe2 双亚晶格反解: 从各磁态的静态总能量解出 J1/J2/J3/J4 (+E0)

原理 (与 build_fgt_cases.py 共用同一套系数矩阵):
    E(supercell) = N_cell · E0 + Σ_k C_k · J~_k ,   J~_k = J_k · S²
    M = [1 | C]  (行 = 磁态),  解  M · [E0, J~] = E_vec   ->  J_k = J~_k / S²
系数矩阵来自 cases/magnetic_orders.json (由 build_fgt_cases.py 写出并校验过条件数),
所以本脚本不自己数键 —— 保证"系数"和"反解"用的是同一份数据。

用法:
  # 真跑 DFT 之后 (runs/<order>/OUTCAR, 或退化成 runs/<order>/E.txt 里一个数字):
  python solve_fgt_J.py <cases_dir> --runs <runs_dir> --spin 1.0

  # 自检 (不需要 DFT: 用已知 J 正演能量再反解, 应精确还原):
  python solve_fgt_J.py <cases_dir> --demo
"""
import argparse
import json
import os
import re
import sys

import numpy as np

EV2MEV = 1000.0


def read_energy(run_dir):
    """优先 OUTCAR 的最后一个 'free  energy   TOTEN  ='; 否则读 E.txt 里的单个数字。"""
    outcar = os.path.join(run_dir, "OUTCAR")
    if os.path.exists(outcar):
        pat = re.compile(r"free\s+energy\s+TOTEN\s*=\s*(-?\d+\.\d+)")
        vals = pat.findall(open(outcar, encoding="utf-8", errors="ignore").read())
        if vals:
            return float(vals[-1]), "OUTCAR"
    e = os.path.join(run_dir, "E.txt")
    if os.path.exists(e):
        txt = open(e, encoding="utf-8").read().strip().split()
        if txt:
            return float(txt[0]), "E.txt"
    return None, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", help="cases 目录 (含 magnetic_orders.json)")
    ap.add_argument("--runs", default=None, help="runs/<order>/ 目录 (默认 cases 的同级 runs)")
    ap.add_argument("--spin", type=float, default=1.0, help="自旋 S (J = J~/S^2)")
    ap.add_argument("--demo", action="store_true", help="自检: 用已知 J 正演再反解")
    args = ap.parse_args()

    meta_path = os.path.join(args.cases, "magnetic_orders.json")
    meta = json.load(open(meta_path, encoding="utf-8"))
    rows = meta["coefficient_matrix"]["rows"]
    cols = meta["coefficient_matrix"]["cols"]          # ["E0", J1..]
    M = np.array(meta["coefficient_matrix"]["M"], float)
    n_cell = meta["n_cell"]
    labels = cols[1:]

    print(f"体系: {meta['system']} | 超胞 {meta['supercell']} | {n_cell} 胞 | "
          f"磁态 {len(rows)} 个 | 壳 {labels}")
    print(f"系数矩阵条件数 cond = {meta['coefficient_matrix']['cond']:.2f}")

    if args.demo:
        # 正演: 取一组"已知" J (meV) -> E (eV) -> 反解应精确还原
        truth_mev = {labels[i]: 10.0 * (i + 1) for i in range(len(labels))}
        truth = np.array([truth_mev[l] * args.spin ** 2 for l in labels])  # J~ (meV)
        rng = np.random.default_rng(7)
        E0_ev = -50.0
        noise = 0.0
        E_vec = M[:, 0] * E0_ev + M[:, 1:] @ (truth / EV2MEV) + noise
        src = "demo(解析)"
    else:
        if not args.runs:
            args.runs = os.path.join(os.path.dirname(os.path.abspath(args.cases)), "runs")
        if not os.path.isdir(args.runs):
            raise SystemExit(f"找不到 runs 目录: {args.runs}")
        E0_ev, E_vec, src = None, [], []
        for r in rows:
            d = os.path.join(args.runs, r)
            e, how = read_energy(d)
            if e is None:
                raise SystemExit(
                    f"runs/{r}/ 里没有 OUTCAR 也没有 E.txt。\n"
                    f"  跑完静态自洽后要么保留 OUTCAR, 要么把总能量 (eV, 一个数) 写进 E.txt。")
            E_vec.append(e)
            src.append(f"{r}:{how}")
            print(f"  {r:12s} E = {e:.6f} eV   [{how}]")
        E_vec = np.array(E_vec)
        if len(E_vec) < len(cols):
            raise SystemExit(f"磁态数 {len(E_vec)} < 未知数 {len(cols)} (E0 + {len(labels)} 个 J), 解不出")
        E0_ev = None

    # 最小二乘 (磁态数 > 未知数时更稳)
    sol, res, rank, sv = np.linalg.lstsq(M, E_vec, rcond=None)
    E0 = sol[0]
    J_sol_mev = sol[1:] * EV2MEV / args.spin ** 2

    print(f"\n反解结果 (S = {args.spin}, 约定 H = -Σ_k J_k Σ S·S, J>0 = FM):")
    print(f"  E0 (每胞) = {E0 * EV2MEV:.6f} meV")
    for l, v in zip(labels, J_sol_mev):
        kind = "铁磁 FM" if v > 0 else "反铁磁 AFM"
        print(f"  {l:3s} = {v:+9.4f} meV   ({kind})")
    resid = float(np.linalg.norm(M @ sol - E_vec))
    print(f"  残差 ||M·x - E|| = {resid:.3e} eV" +
          ("" if resid < 1e-6 else "  [注意: 残差偏大, 检查磁态是否收敛/是否同一结构]"))

    if args.demo:
        ok = all(abs(J_sol_mev[i] - truth_mev[labels[i]]) < 1e-6 for i in range(len(labels)))
        print(f"\n  已知真值: " + ", ".join(f"{l}={truth_mev[l]:.1f}" for l in labels))
        print("  DEMO " + ("PASS (精确还原)" if ok else "FAIL"))
        sys.exit(0 if ok else 1)

    out = {"E0_meV_per_cell": float(E0 * EV2MEV), "spin_S": args.spin,
           "J_meV": {l: float(v) for l, v in zip(labels, J_sol_mev)},
           "residual_eV": resid, "sources": src,
           "convention": "H = -sum_k J_k sum_<ij>_k S.S (J>0 = FM)"}
    with open(os.path.join(os.path.dirname(os.path.abspath(args.cases)), "J_solved.json"),
              "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"\n写出 J_solved.json")
    print("> 局限: 共线双亚晶格近似; 只到纳入的壳; 结构未弛豫; MAE/各向异性未纳入; "
          "U 值与 k 网格影响 J 的绝对值。")


if __name__ == "__main__":
    main()
