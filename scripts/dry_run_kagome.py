#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""无 VASP 的 P4.2 dry-run: 用**合成 OUTCAR** 把「装配 → 回收闸门 → 反解 → 凭证」整条链走通。

用途: 上真实机时之前, 确认失败只可能出在**物理**, 不会出在**管道**。
等价于 `knowledge/KAGOME_RUNBOOK.md` §5.2 的手工步骤, 但一条命令跑完。

**这个脚本产出的任何 J 都不是物理结论**:
  - 能量是 `E = Σ c_k · (J_k · S²) / 1000` 由 `--synthetic` 给的 J **合成**的;
  - 所以反解必然逐位还原那个 J —— 这只证明管道通, 不证明 kagome 的 J 是多少;
  - `J_result.json` 的 `source` 会被强制写成 `demo`, **不许**当 vasp 结果引用;
  - 产物文件名默认带 `.dryrun`, 免得混进真实结果目录。

    python scripts/dry_run_kagome.py <outdir> [--synthetic J1 J2 J3] [--cases <dir>]
                                              [--keep] [--mc] [--quiet]

`<outdir>` 必须先有 `02_supercell/cases/magnetic_orders.json`（先跑一次
`run_2d_magnet.py --model fe3sn --skip-dft --skip-mc --outdir <outdir>`）。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

SCRIPTS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPTS)

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors="replace")
    except Exception:
        pass

import ingest_dft as IG                       # noqa: E402
import j_result as JR                         # noqa: E402

BANNER = ("!" * 78 + "\n"
          "!! DRY-RUN: 下面所有能量都是**合成**的, 反解出的 J 只是你自己给的那个数。\n"
          "!! 这证明管道通, **不证明** kagome 的 J 或 T_C 是多少。\n"
          "!! 产物 source=demo, 文件名带 .dryrun —— 不要当物理结果引用。\n"
          + "!" * 78)


def say(msg, quiet=False):
    if not quiet:
        print(msg)


def cases_of(outdir):
    return os.path.join(outdir, "02_supercell", "cases")


def poscar_lattice(cases, order):
    """从 cases/<order>/POSCAR 取 (scale 后的) 晶格矢量, 供合成 OUTCAR 对齐。"""
    p = os.path.join(cases, order, "POSCAR")
    ln = [x for x in open(p, encoding="utf-8").read().splitlines() if x.strip()]
    sc = float(ln[1].split()[0])
    return [[float(v) * sc for v in ln[2 + i].split()[:3]] for i in range(3)]


def synth_runs(cases, orders, J, runs_dir, quiet=False):
    """按 coefficients_per_cell 合成四个 OUTCAR (只含闸门会查的标志)。"""
    with open(os.path.join(cases, "magnetic_orders.json"), encoding="utf-8-sig") as f:
        meta = json.load(f)
    coeff = meta["coefficients_per_cell"]
    S2 = float(meta["S2"])
    n_mag = int(meta["n_mag"])
    n_atoms = int(meta["n_atoms"])
    # 元素块: 从 POSCAR 的 species/counts 行取, 保证与 cases 的结构指纹一致
    ln = [x for x in open(os.path.join(cases, orders[0], "POSCAR"),
                          encoding="utf-8").read().splitlines() if x.strip()]
    counts = " ".join(ln[6].split())
    lat = poscar_lattice(cases, orders[0])
    missing = [o for o in orders if o not in coeff]
    if missing:
        raise SystemExit(f"magnetic_orders.json 的 coefficients_per_cell 缺这些序: {missing}\n"
                         f"  (是不是用了 --all-orders 之外的序集? 本脚本只吃 selected 四态)")
    E = {}
    for nm in orders:
        e = sum(coeff[nm][m] * (J[m] * S2) for m in range(3)) / 1000.0    # meV -> eV
        E[nm] = e
        os.makedirs(os.path.join(runs_dir, nm), exist_ok=True)
        IG.synth_outcar(os.path.join(runs_dir, nm, "OUTCAR"), e, steps=1,
                        n_ions=n_atoms, counts=counts, lattice=lat)
        say(f"  [synth] runs/{nm}/OUTCAR   E = {e:+.8f} eV   "
            f"(由 J=({J[0]}, {J[1]}, {J[2]}) meV 合成)", quiet)
    return E, meta, n_mag


def main(argv=None):
    ap = argparse.ArgumentParser(description="P4.2 无 VASP dry-run (合成 OUTCAR -> 闸门 -> 反解)")
    ap.add_argument("outdir", help="工程目录 (需先有 02_supercell/cases/magnetic_orders.json)")
    ap.add_argument("--cases", default=None, help="cases 目录 (默认 <outdir>/02_supercell/cases)")
    ap.add_argument("--synthetic", nargs=3, type=float, default=[5.0, -0.5, 0.2],
                    metavar=("J1", "J2", "J3"),
                    help="合成用的 J (meV), 只为验证管道; 默认 5.0 -0.5 0.2")
    ap.add_argument("--mc", action="store_true",
                    help="顺带用这份**合成** J 跑一次 mc_2d (只为验管道; T_C 无数值意义)")
    ap.add_argument("--keep", action="store_true", help="保留临时 runs/ 目录 (默认删)")
    ap.add_argument("--force", action="store_true",
                    help="放行 `not_shadowed` 闸门 (第二次跑 dry-run 时必需 —— 合成文件是"
                         "**确定性**的, 第二次必然与 cases/ 里的逐位相同, 会被防陈旧闸门拦下。"
                         "**只对合成数据用这个开关**; 真数据上被它拦下就说明忘了重算)")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args(argv)
    quiet = a.quiet

    print(BANNER)
    cases = a.cases or cases_of(a.outdir)
    mpath = os.path.join(cases, "magnetic_orders.json")
    if not os.path.isfile(mpath):
        raise SystemExit(f"找不到 {mpath}\n"
                         "  先跑: python run_2d_magnet.py --model fe3sn --skip-dft --skip-mc "
                         f"--outdir {a.outdir}")
    with open(mpath, encoding="utf-8-sig") as f:
        meta0 = json.load(f)
    orders = list(meta0.get("selected_orders") or [])
    if not orders:
        raise SystemExit(f"{mpath} 没有 selected_orders, 无法确定要合成哪几个序")
    J = list(a.synthetic)

    say(f"[1/5] 合成 OUTCAR: {len(orders)} 个序 -> runs/   J = {J} meV (合成, 非物理)", quiet)
    tmp_runs = tempfile.mkdtemp(prefix="_dryrun_runs_")
    E, meta, n_mag = synth_runs(cases, orders, J, tmp_runs, quiet)

    say("\n[2/5] 闸门预演 (--dry-run: 只校验不落盘)", quiet)
    ok, prob, warn, stats, _mp = IG.ingest(tmp_runs, cases, orders=orders, min_kb=0,
                                           dry_run=True, force=a.force, quiet=quiet)
    if not ok:
        print("  闸门没过 (这才是 dry-run 的价值: 在这里暴露):")
        for o, c, m in prob:
            print(f"    [{o or '-'}] {c}: {m}")
        raise SystemExit(1)
    say("  11 条闸门全过 (合成文件只含闸门会查的标志)", quiet)
    for o, c, m in warn:
        say(f"  [warn] [{o or '-'}] {c}: {m}", quiet)

    say("\n[3/5] 真收回 (写入 cases/ + ingest_manifest.json)", quiet)
    ok, prob, warn, stats, mpath2 = IG.ingest(tmp_runs, cases, orders=orders, min_kb=0,
                                              force=True, quiet=quiet)
    if not ok:
        for o, c, m in prob:
            print(f"    [{o or '-'}] {c}: {m}")
        raise SystemExit(1)
    say(f"  manifest: {mpath2}", quiet)

    say("\n[4/5] 反解 (J 应逐位还原上面的合成值)", quiet)
    jout = os.path.join(cases, "J_result.dryrun.json")
    cmd = [sys.executable, os.path.join(SCRIPTS, "solve_J_kagome.py"), cases,
           "--synthetic", str(J[0]), str(J[1]), str(J[2]), "--jfile-out", jout]
    r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace",
                       env=dict(os.environ, PYTHONIOENCODING="utf-8"))
    if not quiet:
        print(r.stdout[-2500:])
    if r.returncode != 0:
        print(r.stderr[-2000:])
        raise SystemExit(r.returncode)
    if not os.path.isfile(jout):
        raise SystemExit(f"反解没有写出 {jout}")
    obj = json.load(open(jout, encoding="utf-8"))
    err = max(abs(obj["J_meV"][f"J{i+1}"] - J[i]) for i in range(3))

    # dry-run 产物必须自曝身份: 不许伪装成 vasp 结果。**用 canonical builder 重建**而不是
    # 手改字段 —— 手改 source 会让 confidence 与 source 自相矛盾 (j_result.validate 会抓到,
    # 这正是它的职责: confidence 是**算出来的**)。重建成 source=demo ⇒ confidence=low。
    A = [[1.0] + list(meta["coefficients_per_cell"][o]) for o in orders]
    shells = [(nm, meta["shells_angstrom"][nm]) for nm in ("J1", "J2", "J3")]
    rebuilt = JR.build(lattice="kagome", shells=shells,
                       J=dict(J1=J[0], J2=J[1], J3=J[2]),
                       E=dict((o, E[o]) for o in orders), E0=obj["E0_eV"],
                       cond=obj["cond"], residual=obj["residual_eV"],
                       spin=meta["spin"], S2=meta["S2"], n_mag=meta["n_mag"],
                       orders=orders, A=A, supercell=meta.get("supercell"),
                       material=meta.get("material") or obj.get("material"),
                       source="demo", z_from="fm",
                       notes=["⚠ dry_run_kagome.py 产物: 能量由合成 J 生成, **不是 DFT**, "
                              "不构成任何物理结论", "往返误差 = %.3g meV" % err])
    rebuilt["dft_provenance"] = obj.get("dft_provenance")
    JR.write(jout, rebuilt)
    obj = rebuilt
    print(f"  往返误差 (反解 J vs 合成 J) = {err:.3g} meV")
    if err > 1e-6:
        raise SystemExit(f"dry-run 往返误差 {err:.3g} meV > 1e-6 —— 管道有问题")
    w = JR.validate(obj, strict=False)
    # confidence=low 是**预期**的 (source=demo), 不是缺陷; 其余警告才算问题。
    problems = [m for m in (w or []) if "confidence" not in m and "confidence" not in str(m)]
    if obj.get("confidence") != "low":
        raise SystemExit(f"dry-run 产物的 confidence 应为 low (source=demo), 实际 {obj.get('confidence')}")
    for m in problems:
        raise SystemExit(f"dry-run 产物没通过 j_result.validate: {m}")
    print(f"  j_result.validate: source=demo -> confidence=low (预期; 这才是真实 vasp 路径的对照)")
    print(f"  警告条数 (除 confidence 外) = {len(problems)}")

    say("\n[5/5] 结论", quiet)
    print(f"  管道: PASS (装配/闸门/收回/反解/凭证 全通)")
    print(f"  产物: {jout}  (source=demo, **不是物理结果**)")
    if a.mc:
        model = os.path.join(cases, "mc_model.dryrun.json")
        cmd = [sys.executable, os.path.join(SCRIPTS, "poscar_to_mc3.py"),
               os.path.join(cases, orders[0], "POSCAR"), "--element", meta["element"],
               "--S", str(meta["spin"]), "--J1", str(J[0]), "--J2", str(J[1]),
               "--J3", str(J[2]), "-o", model]
        r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                           errors="replace", env=dict(os.environ, PYTHONIOENCODING="utf-8"))
        if r.returncode != 0:
            print("  [mc] poscar_to_mc3 失败:"); print(r.stderr[-1500:])
        else:
            out = os.path.join(cases, "mc_dryrun")
            cmd = [sys.executable, os.path.join(SCRIPTS, "mc_2d.py"), model,
                   "--order", "auto", "--size", "12", "12", "--tmin", "5", "--tmax", "60",
                   "--nt", "6", "--out", out]
            r = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                               errors="replace",
                               env=dict(os.environ, PYTHONIOENCODING="utf-8"))
            print("  [mc] " + ("跑了 (T_C 是合成 J 的伪临界点, 无数值意义)"
                               if r.returncode == 0 else "失败: " + r.stderr[-800:]))
    if not a.keep:
        shutil.rmtree(tmp_runs, ignore_errors=True)
    else:
        print(f"  保留合成 OUTCAR: {tmp_runs}")
    print("\n  ⇒ 到此为止能证明的是**管道**, 不是物理。下一步需要真实机时 (见 runbook §5.0–5.1)。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
