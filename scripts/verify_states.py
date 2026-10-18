#!/usr/bin/env python
"""verify_states.py — 反解 J 之前, 逐位点核对每个磁态**真的**收敛到了目标 sigma

为什么必须有这一道 (通用版, 覆盖全部共线晶格):
    能量映射 E(supercell) = N_cell*E0 + sum_k C_k*J~_k 的前提是"每个磁态真的是我们
    指定的那个自旋模式"。VASP 只照初始 MAGMOM 起跑, 收敛到哪个态由它自己定 ——
    层间/亚晶格耦合不足时, 反向对会翻回铁磁, 于是:
      * 该态的能量其实是另一个态的能量 -> 系数矩阵对不上 -> J 系统性偏掉;
      * 而且**不会报错**: 退出码 0、残差也可能很小。更糟的是若两个磁态**互相坍缩**
        (都收敛到同一个态), 它们的能量相等 -> 残差 0 -> `j_result.confidence_of` 会
        因为 residual 小而给出 **high**(它从来看不到磁矩)。这是本项目最贵的一类静默错。

    以前的逐位点核对只存在于 FGT 那条链(fgt_verify_states.py, 由人手/agent 记得去调),
    其余 7 类晶格没有任何逐位点核对。本脚本是它的通用版本, 并被 `run_2d_magnet.py solve`
    当**硬闸门**调用(fail-close), 不再依赖"记得跑"。

做法:
    对 cases/<state>/magmom.txt 的目标磁矩 (POSCAR 顺序) 与 runs/<state>/OUTCAR 里
    **最后一个** magnetization (x) 块的每离子 tot 列逐原子比:
      (1) 符号必须与目标一致   (2) |m| 与目标之差 <= --tol (muB, 只提示不改判)
    再做**坍缩检测**: 两个磁态的实际磁矩签名(符号 + |m| 取整)若相同, 说明它们收敛到了
    同一个态 —— 这是能量映射彻底失效的情形, 单独 `[COLLAPSE]` 点名并判 FAIL。

用法:
    python verify_states.py <cases_dir> [--runs <runs_dir>] [--tol 0.30]
                            [--allow-missing] [--report-only] [--json <path>] [--quiet]
    python verify_states.py --selftest          # 合成夹具, 不需要 model/DFT

退出码: 0 PASS(或 --report-only / 非共线不适用), 1 FAIL(磁态不符/坍缩/缺 OUTCAR)
判定标记: VERDICT: PASS / VERDICT: FAIL / VERDICT: SKIP (ASCII, 便于脚本抓)

覆盖边界(诚实声明):
    * 只覆盖**共线**磁态 (magnetization (x) 的 tot 列)。非共线/SOC 体系 (P3) 的矢量和
      不能这样逐位点比, 本脚本打印 [SKIP] 并以 exit 0 结束, JSON 里 applicable=false
      —— 那条路径仍只有 ingest_dft.py 的净磁化软闸门, 不假装已经核对过。
"""
import argparse
import json
import os
import re
import shutil
import sys
import tempfile

MAGHEAD = "magnetization (x)"
SIGMA = re.compile(r"sigma=([+-]?\d+)")


# --------------------------------------------------------------------------- 读取
def read_energy(run_dir):
    """OUTCAR 的最后一个 free energy TOTEN (口径统一在 dft_energy.py, P3.17);
    否则 E.txt 第一个数 (早期手工提交留下的兜底格式)。"""
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
    """取最后一个 magnetization (x) 块的每离子 tot 列 (POSCAR 顺序)。

    真实 VASP 每行是 `index s p d tot` 五列; ingest_dft.synth_outcar 造的合成文件只写到
    `index s p d` 四列 (列名仍带 tot) —— 两种都认, 否则合成夹具会被误判成"读不到磁矩"。
    """
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
        if p[0].isdigit() and len(p) >= 4:
            try:
                vals.append(float(p[4] if len(p) >= 5 else p[3]))
            except ValueError:
                pass
    return vals or None


def outcar_is_noncollinear(outcar):
    """OUTCAR 自己声明了非共线/SOC 就跳过逐位点核对 (不靠用户命令行说了什么)。"""
    try:
        txt = open(outcar, encoding="utf-8", errors="ignore").read()
    except OSError:
        return False
    return bool(re.search(r"LNONCOLLINEAR\s*=\s*T", txt) or re.search(r"LSORBIT\s*=\s*T", txt))


def read_targets(magmom_path):
    """读目标磁矩 -> [(label, moment)] 每原子一项 (非磁性原子 moment=0)。

    两种格式都认 (仓库里历史原因并存):
      * 驱动 build 写的单行 `MAGMOM = m1 m2 ...`  (共线: 1 分量/原子; 非共线: 3 分量/原子)
      * FGT 链的表格行 `Fe 1 1 1.0 1.0 1.0 sigma=+1 2.0` (末列是目标磁矩)
    """
    out = []
    for l in open(magmom_path, encoding="utf-8", errors="ignore"):
        s = l.split("#")[0].strip()
        if not s:
            continue
        if s.upper().startswith("MAGMOM"):
            vals = []
            for tok in s.split("=", 1)[-1].split():
                try:
                    vals.append(float(tok))
                except ValueError:
                    pass
            for i, v in enumerate(vals):
                out.append((f"#{i + 1}", v))
            return out
        p = s.split()
        if len(p) < 4 or not p[1].isdigit():
            continue
        try:
            mom = float(p[-1])
        except ValueError:
            mom = 0.0
        out.append((p[0], mom))
    return out


def state_names(meta):
    """磁态名单: 共享规格层写 `selected_orders`; FGT 链写 `coefficient_matrix.rows`。"""
    if meta.get("selected_orders"):
        return list(meta["selected_orders"])
    rows = (meta.get("coefficient_matrix") or {}).get("rows")
    if rows:
        return list(rows)
    raise SystemExit("[error] magnetic_orders.json 里既没有 selected_orders 也没有 "
                     "coefficient_matrix.rows, 无法确定要核对哪些磁态")


# --------------------------------------------------------------------------- 核对
def run_check(cases, runs=None, tol=0.30, allow_missing=False, quiet=False):
    """核心核对, 返回 (report dict, exit code)。不打印时不写 stdout (quiet)。"""
    cases = os.path.abspath(cases)
    # OUTCAR 的规范位置是 cases/<state>/OUTCAR (ingest 收进来的那一份); runs 只是回收来源。
    # 默认去 <outdir>/runs 找, 其中 outdir = <outdir>/02_supercell/cases 往上两层。
    runs = runs or os.path.join(os.path.dirname(os.path.dirname(cases)), "runs")
    meta_path = os.path.join(cases, "magnetic_orders.json")
    if not os.path.exists(meta_path):
        raise SystemExit(f"[error] 缺 {meta_path} (先跑 build)")
    meta = json.load(open(meta_path, encoding="utf-8"))
    names = state_names(meta)
    element = meta.get("element") or "M"
    noncollinear = bool(meta.get("noncollinear"))

    def outcar_of(r):
        p = os.path.join(cases, r, "OUTCAR")
        return p if os.path.exists(p) else os.path.join(runs, r, "OUTCAR")

    _n_atoms = meta.get("n_atoms") or meta.get("n_ions")
    if not noncollinear and _n_atoms:
        try:
            _tg0 = read_targets(os.path.join(cases, names[0], "magmom.txt"))
        except OSError:
            _tg0 = []
        if _tg0 and len(_tg0) == 3 * int(_n_atoms):
            noncollinear = True   # MAGMOM 是 3 分量/原子 ⇒ 非共线
    if not noncollinear:
        # 以 OUTCAR 自己的声明为准 (非共线运行里 in-plane 的 120 度态不是"符号不一致")
        for r in names:
            _oc = outcar_of(r)
            if os.path.exists(_oc) and outcar_is_noncollinear(_oc):
                noncollinear = True
                break

    report = {
        "system": meta.get("system") or meta.get("input") or "",
        "cases": cases,
        "runs": os.path.abspath(runs),
        "element": element,
        "tol": tol,
        "noncollinear": noncollinear,
        "applicable": not noncollinear,
        "states": [],
        "collapse_pairs": [],
        "counts": {"states": len(names), "ok": 0, "mismatch": 0, "missing": 0, "collapse": 0},
        "verdict": "PASS",
    }

    def say(*a):
        if not quiet:
            print(*a)

    say(f"体系: {report['system'] or '(未命名)'} | 元素 {element} | 磁态 {len(names)}"
        f" | cases {cases}")
    say(f"runs: {report['runs']} | tol {tol:.2f} muB")
    if noncollinear:
        say("[SKIP] 非共线/SOC 体系: magnetization (x) 的 tot 列不能逐位点当矢量比, "
            "本闸门不适用 (仍靠 ingest_dft.py 的净磁化软闸门)")
        report["verdict"] = "SKIP"
        return report, 0

    say(f"{'磁态':<16}{'状态':<12}{'符号一致':<12}{'|m| 均值':<16}"
        f"{'SumM(' + element + ')':<14}{'E (eV)'}")
    say("-" * 86)

    sigs = {}
    for r in names:
        row = {"name": r, "status": "", "agree": "", "mean_m": None, "sd_m": None,
               "sum_m": None, "E_eV": None, "mismatch": []}
        tg = read_targets(os.path.join(cases, r, "magmom.txt"))
        idx = [i for i, (_lab, mm) in enumerate(tg) if abs(mm) > 1e-9]
        outcar = outcar_of(r)
        if not os.path.exists(outcar):
            report["counts"]["missing"] += 1
            e = read_energy(os.path.join(runs, r))
            row.update(status="[MISSING]", E_eV=e)
            say(f"{r:<16}{'[MISSING]':<12}{'-':<12}{'-':<16}{'-':<14}"
                f"{(f'{e:.6f}' if e is not None else '(无能量)')}")
            report["states"].append(row)
            continue
        act = read_final_moments(outcar)
        if not act:
            report["counts"]["mismatch"] += 1
            row.update(status="[NO-MAG]")
            say(f"{r:<16}{'[NO-MAG]':<12}无法读取 magnetization 块 (非自洽/非共线?)")
            report["states"].append(row)
            continue
        if len(act) < len(tg):
            report["counts"]["mismatch"] += 1
            row.update(status="[BAD-ION]")
            say(f"{r:<16}{'[BAD-ION]':<12}离子数 {len(act)} < POSCAR {len(tg)}")
            report["states"].append(row)
            continue

        bad = []
        for i in idx:
            tgt = tg[i][1]
            a = act[i]
            if (a > 0) != (tgt > 0):
                bad.append({"ion": i + 1, "element": tg[i][0], "target": tgt, "actual": a})
            elif abs(abs(a) - abs(tgt)) > tol:
                bad.append({"ion": i + 1, "element": tg[i][0], "target": tgt, "actual": a,
                            "note": "|m| 超差"})
        mags = [abs(act[i]) for i in idx]
        mean_m = sum(mags) / len(mags) if mags else 0.0
        sdev = (sum((m - mean_m) ** 2 for m in mags) / len(mags)) ** 0.5 if mags else 0.0
        sm = sum(act[i] for i in idx)
        e = read_energy(os.path.dirname(outcar))
        tag = "[OK]" if not bad else "[MISMATCH]"
        if bad:
            report["counts"]["mismatch"] += 1
        else:
            report["counts"]["ok"] += 1
        row.update(status=tag, agree=f"{len(idx) - len(bad)}/{len(idx)}",
                   mean_m=mean_m, sd_m=sdev, sum_m=sm, E_eV=e,
                   mismatch=[{k: v for k, v in b.items()} for b in bad])
        say(f"{r:<16}{tag:<12}{f'{len(idx) - len(bad)}/{len(idx)}':<12}"
            f"{f'{mean_m:.2f}+-{sdev:.2f}':<16}{f'{sm:+.2f}':<14}"
            f"{(f'{e:.6f}' if e is not None else '-')}")
        if row["agree"] != "":
            sig = tuple((1 if act[i] > 0 else -1, round(abs(act[i]), 1)) for i in idx)
            sigs.setdefault(sig, []).append(r)
        report["states"].append(row)

    # 坍缩: 两个磁态的实际磁矩签名完全相同 (四舍五入到 0.1 muB)
    for sig, group in sigs.items():
        if len(group) > 1:
            report["collapse_pairs"].append(group)
            report["counts"]["collapse"] += len(group)
    if report["collapse_pairs"]:
        say("")
        for g in report["collapse_pairs"]:
            es = [next((s["E_eV"] for s in report["states"] if s["name"] == n), None) for n in g]
            extra = ""
            if all(e is not None for e in es) and max(es) - min(es) < 1e-5:
                extra = f" (能量也相同: Δ={(max(es) - min(es)):.2e} eV)"
            say(f"  [COLLAPSE] {' / '.join(g)} 收敛到了同一个磁态 —— 它们的行在系数矩阵里"
                f"本该不同, 能量却同源{extra}")

    say("-" * 86)
    say(f"可核对 {report['counts']['ok'] + report['counts']['mismatch']} 个磁态: "
        f"OK {report['counts']['ok']}, 不一致 {report['counts']['mismatch']}, "
        f"缺 OUTCAR {report['counts']['missing']}")
    bad_total = report["counts"]["mismatch"] + len(report["collapse_pairs"])
    if bad_total or (report["counts"]["missing"] and not allow_missing):
        report["verdict"] = "FAIL"
        say("")
        say("VERDICT: FAIL")
        if report["counts"]["mismatch"]:
            say("  -> 有磁态没收敛到指定 sigma: 直接把它算出的能量代进系数矩阵会得到错的 J。")
            say("     可试: 先算好 FM 的 WAVECAR/CHGCAR 再逐态起跑 (ISTART=1, ICHARG=1)、"
                "调大 NELM / 减小 AMIX、给反向对更强的初始 MAGMOM 对比。")
        if report["collapse_pairs"]:
            say("  -> 有磁态互相坍缩(同一实际磁态): 这些行线性相关, 反解会严重病态甚至给出"
                "看似高置信度的结果。必须让它们真的收敛到不同磁序, 或从反解里剔除。")
        if report["counts"]["missing"] and not allow_missing:
            say("  -> 有磁态根本没跑(缺 OUTCAR)。补齐, 或显式加 --allow-missing。")
        return report, 1
    report["verdict"] = "PASS"
    say("")
    say("VERDICT: PASS")
    return report, 0


# --------------------------------------------------------------------------- 自检
def _selftest():
    """合成 4 个磁态: 正常 / 符号翻错 / 两态坍缩; 断言核对结果。"""
    td = tempfile.mkdtemp(prefix="_vs_selftest_")
    fails = []

    def outcar(path, moments):
        with open(path, "w", encoding="utf-8") as f:
            f.write("  free  energy   TOTEN  =      -10.00000000 eV\n")
            f.write(" magnetization (x)\n")
            f.write("\n")
            f.write("    # of ion       s       p       d       tot\n")
            f.write("    ------------------------------------------\n")
            for i, m in enumerate(moments):
                f.write(f"    {i + 1:<3}  0.000   0.000   0.000   {m:8.3f}\n")
            f.write("    ------------------------------------------\n")
            f.write("    tot          0.000   0.000   0.000   \n")

    try:
        cases = os.path.join(td, "02_supercell", "cases")
        runs = os.path.join(td, "02_supercell", "runs")
        meta = {"system": "selfcheck", "element": "Fe", "noncollinear": False,
                "selected_orders": ["fm", "afm_a", "afm_b", "afm_c"], "n_mag": 4}
        os.makedirs(cases)
        json.dump(meta, open(os.path.join(cases, "magnetic_orders.json"), "w",
                             encoding="utf-8"))
        targets = {"fm": [2, 2, 2, 2], "afm_a": [2, 2, -2, -2],
                   "afm_b": [2, -2, 2, -2], "afm_c": [2, -2, -2, 2]}
        actual = {"fm": [2.1, 2.1, 2.1, 2.1], "afm_a": [2.1, 2.1, -2.1, -2.1],
                  "afm_b": [2.1, 2.1, -2.1, -2.1], "afm_c": [2.1, -2.1, 2.1, -2.1]}
        for st, tgt in targets.items():
            os.makedirs(os.path.join(cases, st))
            with open(os.path.join(cases, st, "magmom.txt"), "w", encoding="utf-8") as f:
                for i, m in enumerate(tgt):
                    f.write(f"Fe {i + 1} {i + 1} 1.0 1.0 1.0  sigma={1 if m > 0 else -1} {m:.1f}\n")
            os.makedirs(os.path.join(runs, st))
            outcar(os.path.join(runs, st, "OUTCAR"), actual[st])

        # afm_b 的实际磁矩与自己的目标也不符 (它坍缩到了 afm_a 的态) —— 这正是要抓的东西
        rep, code = run_check(cases, runs, tol=0.30, quiet=True)
        got = [s["name"] for s in rep["states"] if s["status"] == "[MISMATCH]"]
        if code != 1:
            fails.append(f"应有不一致但 exit={code}")
        if got != ["afm_b", "afm_c"]:
            fails.append(f"不一致磁态应为 afm_b/afm_c, 实际 {got}")
        if rep["collapse_pairs"] != [["afm_a", "afm_b"]]:
            fails.append(f"坍缩对应为 [['afm_a','afm_b']], 实际 {rep['collapse_pairs']}")
        print(f"  [PASS] 符号不一致被点名 ({', '.join(got)}) / 坍缩对被点名 "
              f"(afm_a/afm_b) / VERDICT={rep['verdict']}")

        # 去掉两个坏磁态的 OUTCAR -> 只剩缺文件, 不该再判 FAIL
        for st in ("afm_b", "afm_c"):
            os.remove(os.path.join(runs, st, "OUTCAR"))
        rep2, code2 = run_check(cases, runs, tol=0.30, quiet=True)
        if code2 != 1 or rep2["counts"]["missing"] != 2:
            fails.append(f"缺 OUTCAR 应 FAIL 且 missing=2, 实际 exit={code2} "
                         f"missing={rep2['counts']['missing']}")
        rep3, code3 = run_check(cases, runs, tol=0.30, allow_missing=True, quiet=True)
        if code3 != 0 or rep3["verdict"] != "PASS":
            fails.append(f"--allow-missing 应放行 (exit 0/PASS), 实际 exit={code3} "
                         f"verdict={rep3['verdict']}")
        print(f"  [PASS] 缺 OUTCAR 默认 FAIL; --allow-missing 放行 "
              f"(exit {code2} -> {code3}, verdict {rep2['verdict']} / {rep3['verdict']})")

        meta_nc = dict(meta, noncollinear=True)
        json.dump(meta_nc, open(os.path.join(cases, "magnetic_orders.json"), "w",
                                encoding="utf-8"))
        rep4, code4 = run_check(cases, runs, quiet=True)
        if code4 != 0 or rep4["verdict"] != "SKIP" or rep4["applicable"]:
            fails.append(f"非共线应 SKIP/exit 0, 实际 exit={code4} verdict={rep4['verdict']}")
        print("  [PASS] 非共线体系判 SKIP 而不是假装通过")
    finally:
        shutil.rmtree(td, ignore_errors=True)

    if fails:
        for f in fails:
            print(f"  [FAIL] {f}")
        return 1
    print("verify_states 自检 OK")
    return 0


# --------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description="反解前的逐位点磁态核对 (通用版)")
    ap.add_argument("cases", nargs="?", help="cases 目录 (含 magnetic_orders.json)")
    ap.add_argument("--runs", default=None, help="runs/<state>/OUTCAR 目录 (默认 cases 同级 runs)")
    ap.add_argument("--tol", type=float, default=0.30,
                    help="|m| 相对目标磁矩的允许偏差 (muB), 超差只在表里标出")
    ap.add_argument("--allow-missing", action="store_true",
                    help="缺 OUTCAR 的磁态只警告, 不算失败")
    ap.add_argument("--report-only", action="store_true",
                    help="只报告不改退出码 (供 --allow-unverified-states 用; 会留 JSON 证据)")
    ap.add_argument("--json", default=None, help="把核对摘要写到这个路径")
    ap.add_argument("--quiet", action="store_true", help="不打印表 (只写 --json)")
    ap.add_argument("--selftest", action="store_true", help="跑合成夹具自检, 不需要 model/DFT")
    args = ap.parse_args()

    if args.selftest:
        return _selftest()
    if not args.cases:
        ap.error("缺 cases 目录 (或改用 --selftest)")

    rep, code = run_check(args.cases, args.runs, args.tol, args.allow_missing, args.quiet)
    if args.json:
        os.makedirs(os.path.dirname(os.path.abspath(args.json)), exist_ok=True)
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(rep, f, indent=2, ensure_ascii=False)
        if not args.quiet:
            print(f"[OK] 核对摘要 -> {args.json}")
    if args.report_only and code != 0:
        print("\n[warn] --report-only: 磁态核对没通过, 但按调用方要求不改退出码 "
              "(证据已写进 JSON; 用它算出的 J 视为未验证)", file=sys.stderr)
        return 0
    return code


if __name__ == "__main__":
    sys.exit(main())
