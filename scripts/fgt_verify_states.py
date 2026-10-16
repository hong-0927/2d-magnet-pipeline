#!/usr/bin/env python
"""fgt_verify_states.py — 反解前的最后一道闸: 核对每个磁态是否真收敛到指定的 sigma

为什么必须查:
    能量映射  E(supercell) = N_cell*E0 + sum_k C_k*J~_k  的前提是
    "每个磁态真的就是我们指定的那个自旋模式"。
    VASP 只照初始 MAGMOM 起跑, 收敛到哪个态由它自己定 —— 尤其 afm_b_vert 的 B-B 竖直对
    (2.55 A, 比 A-B 的 2.63 A 还短) 若层间耦合不足, 会把反向对翻回铁磁。
    一旦实际态 != 目标态, 代进系数矩阵的能量就是错的, J 会系统性偏掉, 而且**不会报错**。

做法:
    从 runs/<order>/OUTCAR 取**最后一个** magnetization (x) 块的每离子 tot 列,
    与 cases/<order>/magmom.txt 的目标 sigma 逐原子比符号 (两者都按 POSCAR 顺序)。

用法:
    python fgt_verify_states.py <cases_dir> [--runs <runs_dir>] [--tol 0.30] [--allow-missing]

退出码: 有磁态符号不一致 -> 1 (缺 OUTCAR 见 --allow-missing)
判定标记: VERDICT: PASS / VERDICT: FAIL (ASCII, 便于脚本抓)
"""
import argparse
import json
import os
import re
import sys

# P3.17: OUTCAR 的 TOTEN 正则与取法统一在 `dft_energy.py`（原地的 TOTEN 正则已删除）
MAGHEAD = "magnetization (x)"
SIGMA = re.compile(r"sigma=([+-]?\d+)")


def read_energy(run_dir):
    """OUTCAR 的最后一个 free energy TOTEN（口径统一在 `dft_energy.py`, P3.17）;
    否则 E.txt 第一个数（早期手工提交留下的兜底格式）。"""
    p = os.path.join(run_dir, "OUTCAR")
    if os.path.exists(p):
        from dft_energy import read_toten
        try:
            return read_toten(p)
        except ValueError:
            pass
    p = os.path.join(run_dir, "E.txt")
    if os.path.exists(p):
        t = open(p, encoding="utf-8", errors="ignore").read().split()
        if t:
            try:
                return float(t[0])
            except ValueError:
                return None
    return None


def read_final_moments(outcar):
    """取最后一个 magnetization (x) 块的每离子 tot 列 (POSCAR 顺序)。"""
    lines = open(outcar, encoding="utf-8", errors="ignore").read().splitlines()
    last = None
    for i, l in enumerate(lines):
        if MAGHEAD in l:
            last = i
    if last is None:
        return None
    vals = []
    for l in lines[last + 1:]:
        p = l.split()
        if not p:
            continue
        if p[0] == "tot":
            break
        if p[0].isdigit() and len(p) >= 5:
            try:
                vals.append(float(p[4]))
            except ValueError:
                pass
    return vals or None


def read_targets(magmom_path):
    """-> [(element, sigma, moment)] 按 POSCAR 顺序 (含非磁性原子, moment=0)"""
    out = []
    for l in open(magmom_path, encoding="utf-8", errors="ignore"):
        l = l.split("#")[0].strip()
        if not l:
            continue
        p = l.split()
        if len(p) < 4 or not p[1].isdigit():
            continue
        m = SIGMA.search(l)
        try:
            mom = float(p[-1])
        except ValueError:
            mom = 0.0
        out.append((p[0], int(m.group(1)) if m else 0, mom))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", help="cases 目录 (含 magnetic_orders.json 与 <order>/magmom.txt)")
    ap.add_argument("--runs", default=None, help="runs/<order>/ 目录 (默认 cases 同级 runs)")
    ap.add_argument("--tol", type=float, default=0.30,
                    help="|m| 相对目标磁矩允许的偏差 (muB), 超了也标出来")
    ap.add_argument("--allow-missing", action="store_true",
                    help="缺 OUTCAR 的磁态只警告, 不算失败")
    args = ap.parse_args()

    cases = os.path.abspath(args.cases)
    runs = args.runs or os.path.join(os.path.dirname(cases), "runs")
    meta = json.load(open(os.path.join(cases, "magnetic_orders.json"), encoding="utf-8"))
    rows = meta["coefficient_matrix"]["rows"]

    print(f"体系: {meta['system']} | cases {cases}")
    print(f"runs: {runs}")
    print(f"{'磁态':<14}{'状态':<12}{'符号一致':<12}{'|m| 均值':<14}{'SumM(Fe)':<12}{'E (eV)'}")
    print("-" * 78)

    n_bad, n_missing, n_ok = 0, 0, 0
    details = []
    for r in rows:
        tg = read_targets(os.path.join(cases, r, "magmom.txt"))
        idx = [i for i, (el, sg, mm) in enumerate(tg) if abs(mm) > 1e-9]
        outcar = os.path.join(runs, r, "OUTCAR")
        if not os.path.exists(outcar):
            n_missing += 1
            e = read_energy(os.path.join(runs, r))
            print(f"{r:<14}{'[MISSING]':<12}{'-':<12}{'-':<14}{'-':<12}"
                  f"{(f'{e:.6f}' if e is not None else '(无能量)')}")
            continue
        act = read_final_moments(outcar)
        if not act:
            n_bad += 1
            print(f"{r:<14}{'[NO-MAG]':<12}{'无法读取 magnetization 块 (非自洽/非共线? )'}")
            details.append(f"{r}: OUTCAR 里找不到 magnetization (x) 块")
            continue
        if len(act) < len(tg):
            n_bad += 1
            print(f"{r:<14}{'[BAD-ION]':<12}{f'离子数 {len(act)} < POSCAR {len(tg)}'}")
            details.append(f"{r}: magnetization 块离子数 {len(act)} != {len(tg)}")
            continue

        bad = []
        for i in idx:
            tgt = tg[i][2]
            a = act[i]
            if (a > 0) != (tgt > 0):
                bad.append((i + 1, tg[i][0], tgt, a))
            elif abs(abs(a) - abs(tgt)) > args.tol:
                bad.append((i + 1, tg[i][0], tgt, a))
        agree = len(idx) - len(bad)
        mags = [abs(act[i]) for i in idx]
        mean_m = sum(mags) / len(mags) if mags else 0.0
        sdev = (sum((m - mean_m) ** 2 for m in mags) / len(mags)) ** 0.5 if mags else 0.0
        sm = sum(act[i] for i in idx)
        e = read_energy(os.path.join(runs, r))
        tag = "[OK]" if not bad else "[MISMATCH]"
        if bad:
            n_bad += 1
        else:
            n_ok += 1
        print(f"{r:<14}{tag:<12}{f'{agree}/{len(idx)}':<12}"
              f"{f'{mean_m:.2f}+-{sdev:.2f}':<14}{f'{sm:+.2f}':<12}"
              f"{(f'{e:.6f}' if e is not None else '-')}")
        if bad:
            show = "  ".join(f"Fe#{i}(目标{t:+.0f} 实际{a:+.2f})" for i, el, t, a in bad[:8])
            more = f" ... 共 {len(bad)} 个" if len(bad) > 8 else ""
            details.append(f"{r}: {show}{more}")
            print(f"{'':<14}不一致原子: {show}{more}")

    print("-" * 78)
    print(f"可核对 {n_ok + n_bad} 个磁态: OK {n_ok}, 不一致 {n_bad}, 缺 OUTCAR {n_missing}")
    for d in details:
        print(f"  ! {d}")

    if n_bad or (n_missing and not args.allow_missing):
        print("\nVERDICT: FAIL")
        if n_bad:
            print("  -> 有磁态没有收敛到指定 sigma, 直接把它的能量代进系数矩阵会得到错的 J。")
            print("     可试: 先算好 FM 的 WAVECAR/CHGCAR 再逐个磁态起跑 (ISTART=1, ICHARG=1),")
            print("     或调大 NELM / 减小 AMIX, 或给反向对更强的初始 MAGMOM 对比。")
        if n_missing and not args.allow_missing:
            print("  -> 有磁态根本没跑 (缺 OUTCAR)。要么补齐, 要么显式加 --allow-missing。")
        sys.exit(1)
    print("\nVERDICT: PASS")


if __name__ == "__main__":
    main()
