#!/usr/bin/env python3
"""k 点收敛抽查 (P4.2 环节 2): 生成 KPOINTS 变体 -> 回收后判定"网格多密才够"。

它只做两件事:
  1) plan   : 把 <root>/runs/<order>/{INCAR,POSCAR} 复制成
              <root>/kconv/runs/<order>@k<a>_<b>_<c>/, **只换 KPOINTS**, 并写 kconv_plan.json;
  2) ingest : 读回每个 run 的 OUTCAR (唯一解析入口 scripts/ingest_dft.py), 算
              |E(该网格) - E(最密网格)| (meV/原子), 写 k_convergence.json,
              给出 recommended_kpoints = 偏差已 <= --tol 的**最粗**网格。

**它不能做** (说出来是为了不被误读):
  * 不能替你提交 VASP (本仓库从不代提交);
  * 不能证明 J 收敛 —— J 是能量**差**, 总能量收敛只是 J 可信的**必要**条件;
  * 不能判断 T_N / 基态 / 是否需要 SOC;
  * 没有真实 OUTCAR 时它只会写"待回收", 不会编数字。

退出码: 0 = 通过或还在等 OUTCAR; 1 = 已能判定且判为不收敛 / 数据质量不合格 / 晶胞不一致。
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import ingest_dft as ID      # noqa: E402  复用 OUTCAR 解析 (唯一实现)
import dft_energy as DE      # noqa: E402  复用能量口径 (free energy TOTEN)

VERSION = "1.0"
DEFAULT_GRIDS = "3,4,6,8"
KP_HEADER = "Gamma-centered k-mesh {a} {b} {c} (vacuum along c)"
NOTE = ("本文件只回答『总能量对 k 网格收敛到哪个网格够了』; **不**回答 J 是否收敛 "
        "(J 是能量差, 总能量收敛是必要条件而非充分条件), 也不回答 T_N / 基态 / SOC 是否需要。")
QUALITY_NOTE = ("每个 run 都要过 ingest_dft 的同款基本闸门: 体积 >= --min-kb、有 timing、"
                "SCF 收敛 (EDIFF reached); 不满足则该序**不参与比较**。")

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(errors="replace")
    except Exception:                                     # noqa: BLE001
        pass


# ---------------------------------------------------------------------------
# 网格解析 / KPOINTS 文本
# ---------------------------------------------------------------------------
def parse_grid(tok):
    """'6' -> (6,6,1); '6,6,1' 或 '6 6 1' -> (6,6,1)。非法输入硬报错。"""
    parts = [p for p in re.split(r"[,\s]+", str(tok).strip()) if p]
    if not parts:
        raise SystemExit("k 网格为空; 例: --grids 3,4,6,8 或 --grids '3 3 1;4 4 1'")
    try:
        vals = [int(p) for p in parts]
    except ValueError:
        raise SystemExit(f"k 网格 {tok!r} 不是整数: 只接受整数, 例 6 或 '6 6 1'")
    if len(vals) == 1:
        vals = [vals[0], vals[0], 1]
    if len(vals) != 3:
        raise SystemExit(f"k 网格 {tok!r} 需要 1 个或 3 个整数 (面内 k 或 kx ky kz)")
    if any(v <= 0 for v in vals):
        raise SystemExit(f"k 网格 {tok!r} 必须为正整数")
    return tuple(vals)


def parse_grids(spec):
    """'3,4,6,8' (面内 k, kz=1) 或 '3 3 1;4 4 1' (显式三元组)。去重 + 按面内 k 升序。"""
    spec = str(spec)
    if ";" in spec:
        grids = [parse_grid(x) for x in spec.split(";") if x.strip()]
    elif "," in spec:
        grids = [parse_grid(x) for x in spec.split(",") if x.strip()]
    else:
        grids = [parse_grid(spec)]
    uniq = sorted(set(grids), key=lambda g: (g[0] * g[1] * g[2], g))
    if not uniq:
        raise SystemExit(f"--grids {spec!r} 解析后为空")
    return uniq


def grid_label(g):
    return "k%d_%d_%d" % tuple(g)


def grid_text(g):
    return "%d %d %d" % tuple(g)


def kpoints_text(g):
    """与 scripts/assemble_dft.py 写出的 KPOINTS 同格式 (3 行: 注释 / 网格 / 0 0 0)。"""
    return KP_HEADER.format(a=g[0], b=g[1], c=g[2]) + "\n" + grid_text(g) + "\n0 0 0\n"


def label_to_grid(lbl):
    """'k6_6_1' -> (6,6,1)。"""
    return parse_grid(lbl[1:].replace("_", ","))


# ---------------------------------------------------------------------------
# plan: 复制 runs/ 并只换 KPOINTS
# ---------------------------------------------------------------------------
def order_names(root, orders=None):
    """默认取 <root>/runs/ 下有 INCAR 的序目录 (dft 阶段的产物)。"""
    runs = os.path.join(root, "runs")
    if orders:
        names = list(orders)
        missing = [n for n in names if not os.path.isfile(os.path.join(runs, n, "INCAR"))]
        if missing:
            raise SystemExit(f"--orders 里的 {missing} 在 {runs} 下没有 INCAR; "
                             f"可用的序: {sorted(d for d in os.listdir(runs) if os.path.isdir(os.path.join(runs, d))) if os.path.isdir(runs) else []}")
        return names
    if not os.path.isdir(runs):
        raise SystemExit(f"找不到 {runs}; 先跑 dft 阶段 (run_2d_magnet.py dft ...)")
    names = sorted(d for d in os.listdir(runs)
                   if os.path.isdir(os.path.join(runs, d)) and os.path.isfile(os.path.join(runs, d, "INCAR")))
    if not names:
        raise SystemExit(f"{runs} 下没有含 INCAR 的序目录; 先跑 dft 阶段")
    return names


def plan(root, grids, orders=None, quiet=False):
    """写 <root>/kconv/runs/<order>@<grid>/{INCAR,KPOINTS,POSCAR} + kconv_plan.json。"""
    names = order_names(root, orders)
    runs = os.path.join(root, "runs")
    base = os.path.join(root, "kconv")
    made = []
    for nm in names:
        src = os.path.join(runs, nm)
        with open(os.path.join(src, "INCAR"), encoding="utf-8", errors="replace") as f:
            incar = f.read()
        for g in grids:
            dst = os.path.join(base, "runs", f"{nm}@{grid_label(g)}")
            os.makedirs(dst, exist_ok=True)
            with open(os.path.join(dst, "INCAR"), "w", encoding="utf-8") as f:
                f.write(incar.rstrip("\n")
                        + f"\n! kconv: KPOINTS = {grid_text(g)} (k 点收敛抽查, P4.2)\n")
            with open(os.path.join(dst, "KPOINTS"), "w", encoding="utf-8") as f:
                f.write(kpoints_text(g))
            for extra in ("POSCAR",):
                p = os.path.join(src, extra)
                if os.path.isfile(p):
                    shutil.copy(p, os.path.join(dst, extra))
            made.append({"order": nm, "grids": list(g), "dir": dst})
    os.makedirs(base, exist_ok=True)
    man = {
        "version": VERSION,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "root": os.path.abspath(root),
        "grids": [list(g) for g in grids],
        "orders": names,
        "runs": made,
        "incar_note": "每个 run 的 INCAR 与 runs/<order>/INCAR 逐字节相同, 只多一行 '!' 注释; 只换 KPOINTS。",
        "next": "把每个目录**各提交一次静态 SCF** (同一个 INCAR, 只差 KPOINTS), "
                "把 OUTCAR 放回同目录, 再跑同一条命令 (或加 --ingest-only) 判定。",
        "note": NOTE,
        "quality": QUALITY_NOTE,
    }
    with open(os.path.join(base, "kconv_plan.json"), "w", encoding="utf-8") as f:
        json.dump(man, f, ensure_ascii=False, indent=2)
    if not quiet:
        print(f"[kconv] plan: {len(names)} 个序 x {len(grids)} 个网格 = {len(made)} 个 run")
        print(f"        目录: {os.path.join(base, 'runs')}")
        print(f"        清单: {os.path.join(base, 'kconv_plan.json')}")
        print(f"        网格: {', '.join(grid_text(g) for g in grids)}; 序: {', '.join(names)}")
    return man


# ---------------------------------------------------------------------------
# ingest: 读回 OUTCAR 并判定
# ---------------------------------------------------------------------------
def _scan_runs(root):
    """返回 {(order, grid): dir}。目录名形如 <order>@k<a>_<b>_<c>。"""
    d = os.path.join(root, "kconv", "runs")
    out = {}
    if not os.path.isdir(d):
        return out
    for nm in sorted(os.listdir(d)):
        p = os.path.join(d, nm)
        if not os.path.isdir(p) or "@" not in nm:
            continue
        order, lbl = nm.rsplit("@", 1)
        try:
            g = label_to_grid(lbl)
        except SystemExit:
            continue
        out[(order, g)] = p
    return out


def ingest(root, grids=None, orders=None, tol=1.0, min_kb=200, quiet=False):
    """逐序比较 |E(网格) - E(最密网格)| (meV/原子), 写 k_convergence.json 并返回结果 dict。"""
    found = _scan_runs(root)
    if not found:
        raise SystemExit(f"找不到 {os.path.join(root, 'kconv', 'runs')}; 先跑 plan (不加 --ingest-only)")
    by_order = {}
    for (order, g) in found:
        by_order.setdefault(order, {})[g] = found[(order, g)]
    if orders:
        by_order = {k: v for k, v in by_order.items() if k in set(orders)}

    res = {
        "version": VERSION,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        "root": os.path.abspath(root),
        "tol_meV_per_atom": float(tol),
        "min_kb": int(min_kb),
        "orders": {},
        "note": NOTE,
        "quality": QUALITY_NOTE,
    }
    for order in sorted(by_order):
        cells = by_order[order]
        g_sorted = sorted(cells, key=lambda g: (g[0] * g[1] * g[2], g), reverse=True)  # 密 -> 粗
        rows, problems = [], []
        for g in g_sorted:
            dirp = cells[g]
            op = ID._find_outcar(dirp)
            row = {"grids": list(g), "dir": dirp, "outcar": op}
            if not op:
                row["status"] = "missing_outcar"
                rows.append(row)
                continue
            info = ID.parse_outcar(op, min_kb=min_kb)
            row.update({"energy_eV": info.get("energy_eV"), "n_ions": info.get("n_ions"),
                        "size_kb": info.get("size_kb"), "size_ok": info.get("size_ok"),
                        "has_timing": info.get("has_timing"),
                        "scf_converged": info.get("scf_converged"),
                        "n_kpts": info.get("n_kpts"), "encut": info.get("encut"),
                        "energy_error": info.get("energy_error")})
            bad = []
            if not row["size_ok"]:
                bad.append("volume_too_small")
            if not row["has_timing"]:
                bad.append("no_timing")
            if not row["scf_converged"]:
                bad.append("scf_not_converged")
            if row["energy_eV"] is None:
                bad.append("no_energy")
            by_kpts = row["n_kpts"]
            n_irr_note = None
            if by_kpts is not None:
                n_irr_note = by_kpts          # 只记录, 不当判据
            row["n_kpts_outcar"] = n_irr_note
            row["problems"] = bad
            row["status"] = "ok" if not bad else "quality_fail"
            rows.append(row)
        good = [r for r in rows if r["status"] == "ok"]
        for r in rows:
            if r["status"] == "quality_fail":
                problems.append(f"{order}@{grid_label(r['grids'])}: {','.join(r['problems'])}")
        verdict = {
            "rows": rows,
            "n_usable": len(good),
            "problems": problems,
        }
        if len(good) < 2:
            verdict["status"] = "insufficient"
            verdict["n_note"] = ("至少需要 2 个网格的可用能量才能比较 "
                                 f"(当前可用 {len(good)} / 共 {len(rows)} 个网格)")
        else:
            n_ions_set = {r["n_ions"] for r in good}
            if len(n_ions_set) > 1 or None in n_ions_set:
                verdict["status"] = "cell_mismatch"
                verdict["n_note"] = f"各网格的原子数不一致 {sorted(x for x in n_ions_set if x is not None)}; 能量不可比"
            else:
                n_ions = good[0]["n_ions"] or 0
                ref = good[0]                        # 最密网格
                for r in good:
                    d_eV = r["energy_eV"] - ref["energy_eV"]
                    r["delta_eV"] = d_eV
                    r["delta_meV_per_atom"] = (d_eV * 1000.0 / n_ions) if n_ions else None
                ok_rows = [r for r in good
                           if r["delta_meV_per_atom"] is not None and abs(r["delta_meV_per_atom"]) <= tol]
                conv = all(r["delta_meV_per_atom"] is not None and abs(r["delta_meV_per_atom"]) <= tol
                           for r in good[1:]) if len(good) >= 2 else False
                verdict["status"] = "converged" if conv else "not_converged"
                verdict["reference_grids"] = list(ref["grids"])
                verdict["reference_energy_eV"] = ref["energy_eV"]
                if ok_rows:
                    verdict["recommended_grids"] = list(min(ok_rows, key=lambda r: r["grids"][0] * r["grids"][1])["grids"])
                verdict["max_abs_delta_meV_per_atom"] = max(
                    abs(r["delta_meV_per_atom"]) for r in good[1:])
        res["orders"][order] = verdict

    statuses = {v["status"] for v in res["orders"].values()}
    if not res["orders"]:
        res["status"] = "insufficient"
    elif statuses == {"converged"}:
        res["status"] = "ok"
    elif "not_converged" in statuses:
        res["status"] = "not_converged"
    elif "cell_mismatch" in statuses:
        res["status"] = "cell_mismatch"
    elif "quality_fail" in statuses or "insufficient" in statuses:
        res["status"] = "insufficient" if statuses == {"insufficient"} else "partial"
    else:
        res["status"] = "insufficient"
    res["recommended_kpoints"] = {o: v.get("recommended_grids")
                                  for o, v in res["orders"].items() if v.get("recommended_grids")}
    os.makedirs(os.path.join(root, "kconv"), exist_ok=True)
    with open(os.path.join(root, "kconv", "k_convergence.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    if not quiet:
        render(res)
    return res


def render(res):
    print("")
    print(f"[kconv] 判定: {res['status']}  (判据 |ΔE| <= {res['tol_meV_per_atom']} meV/原子, "
          f"参照 = 每个序最密的可用网格)")
    for order, v in res["orders"].items():
        print(f"  - {order}: {v['status']}"
              + (f" | 建议网格 {grid_text(v['recommended_grids'])}" if v.get("recommended_grids") else "")
              + (f" | {v['n_note']}" if v.get("n_note") else ""))
        for r in v["rows"]:
            e = r.get("energy_eV")
            d = r.get("delta_meV_per_atom")
            print("      %-12s E=%s  Δ=%s  %s"
                  % (grid_text(r["grids"]),
                     f"{e:+.6f} eV" if isinstance(e, float) else "—",
                     f"{d:+.3f} meV/原子" if isinstance(d, float) else "—",
                     r["status"] if r["status"] == "ok" else r["status"] + "(" + ",".join(r.get("problems") or []) + ")"))
    if res["status"] == "insufficient":
        _any_out = any(r.get("outcar") for v in res["orders"].values() for r in v["rows"])
        if _any_out:
            print("      有 OUTCAR 但**没过闸门** (看上面的 quality_fail 行): 体积/timing/SCF 收敛标志"
                  "缺一个就不参与比较; 截断或半成品 OUTCAR 不能拿来判收敛。")
        else:
            print("      还没法判定: 把每个 kconv/runs/<序>@<网格>/ 各提交一次静态 SCF, "
                  "OUTCAR 放回同目录后重跑本命令。")
    print(f"      明细: {os.path.join(res['root'], 'kconv', 'k_convergence.json')}")
    print(f"      注意: {NOTE}")


# ---------------------------------------------------------------------------
# 自检
# ---------------------------------------------------------------------------
def _selftest():
    import tempfile
    ok = fail = 0

    def ck(name, cond, extra=""):
        nonlocal ok, fail
        if cond:
            ok += 1
            print(f"  [ok ] {name}" + (f"  {extra}" if extra else ""))
        else:
            fail += 1
            print(f"  [FAIL] {name}" + (f"  {extra}" if extra else ""))

    # A. 解析
    ck("parse_grid '6' -> (6,6,1)", parse_grid("6") == (6, 6, 1))
    ck("parse_grid '6,6,1' / '6 6 1'", parse_grid("6,6,1") == (6, 6, 1) and parse_grid("6 6 1") == (6, 6, 1))
    for bad in ("0", "6,6", "abc", ""):
        try:
            parse_grid(bad)
            ck(f"parse_grid {bad!r} 应硬报错", False)
        except SystemExit:
            ck(f"parse_grid {bad!r} 硬报错", True)
    ck("parse_grids '3,4,6,8' -> 4 个 (kz=1)",
       parse_grids("3,4,6,8") == [(3, 3, 1), (4, 4, 1), (6, 6, 1), (8, 8, 1)])
    ck("parse_grids '3 3 1;4 4 1' 显式三元组", parse_grids("3 3 1;4 4 1") == [(3, 3, 1), (4, 4, 1)])
    ck("parse_grids 去重 '4,4' -> 1 个", parse_grids("4,4") == [(4, 4, 1)])
    ck("kpoints_text 三行且与 assemble_dft 同格式",
       kpoints_text((6, 6, 1)) == "Gamma-centered k-mesh 6 6 1 (vacuum along c)\n6 6 1\n0 0 0\n")
    ck("label_to_grid 'k3_3_1' -> (3,3,1)", label_to_grid("k3_3_1") == (3, 3, 1))

    tmp = tempfile.mkdtemp(prefix="kconv_selftest_")
    try:
        # B. plan
        root = os.path.join(tmp, "proj")
        for nm in ("fm", "afm_c"):
            d = os.path.join(root, "runs", nm)
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, "INCAR"), "w", encoding="utf-8") as f:
                f.write("ENCUT    = 520\nNPAR     = 4\nMAGMOM = 1 2 3\n")
            with open(os.path.join(d, "KPOINTS"), "w", encoding="utf-8") as f:
                f.write(kpoints_text((3, 3, 1)))
            with open(os.path.join(d, "POSCAR"), "w", encoding="utf-8") as f:
                f.write("POSCAR\n1.0\n")
        grids = parse_grids("3,4,6")
        man = plan(root, grids, quiet=True)
        ck("plan: 2 序 x 3 网格 = 6 个 run", len(man["runs"]) == 6)
        d1 = os.path.join(root, "kconv", "runs", "fm@k6_6_1")
        ck("plan: 目录/INCAR/KPOINTS/POSCAR 都在",
           os.path.isfile(os.path.join(d1, "INCAR")) and os.path.isfile(os.path.join(d1, "KPOINTS"))
           and os.path.isfile(os.path.join(d1, "POSCAR")))
        ck("plan: KPOINTS 内容正确", open(os.path.join(d1, "KPOINTS"), encoding="utf-8").read()
           == kpoints_text((6, 6, 1)))
        body = open(os.path.join(d1, "INCAR"), encoding="utf-8").read()
        ck("plan: INCAR 只多一行 '!' 注释, 参数逐字不变",
           body.splitlines()[:3] == ["ENCUT    = 520", "NPAR     = 4", "MAGMOM = 1 2 3"]
           and body.rstrip().endswith("(k 点收敛抽查, P4.2)"))
        ck("plan: 清单 JSON 键齐全",
           all(k in json.load(open(os.path.join(root, "kconv", "kconv_plan.json"), encoding="utf-8"))
               for k in ("grids", "orders", "runs", "next", "note", "quality")))
        try:
            order_names(root, ["bogus"])
            ck("plan: --orders 里的未知序名应硬报错", False)
        except SystemExit:
            ck("plan: --orders 里的未知序名硬报错", True)

        # C. ingest: 收敛序列
        def write_all(energy_of, counts="6 6", scf=True, root=root, gl=grids):
            for nm in ("fm", "afm_c"):
                for g in gl:
                    dd = os.path.join(root, "kconv", "runs", f"{nm}@{grid_label(g)}")
                    os.makedirs(dd, exist_ok=True)
                    cts = counts(g) if callable(counts) else counts
                    ID.synth_outcar(os.path.join(dd, "OUTCAR"), energy_of(g, nm), counts=cts,
                                    scf=scf, n_kpts=g[0] * g[1])

        write_all(lambda g, nm: -100.0 + 0.01 / (g[0] ** 2))
        res = ingest(root, grids, tol=1.0, min_kb=0, quiet=True)
        ck("ingest: 收敛序列 -> status=ok", res["status"] == "ok", res["status"])
        ck("ingest: 建议网格 = 最粗的合格网格 (3 3 1)",
           res["recommended_kpoints"].get("fm") == [3, 3, 1], str(res["recommended_kpoints"].get("fm")))
        _rows_fm = res["orders"]["fm"]["rows"]
        ck("ingest: 参照行 (最密网格) Δ=0 且行按 密->粗 排列",
           abs(_rows_fm[0]["delta_meV_per_atom"]) < 1e-12 and _rows_fm[0]["grids"] == [6, 6, 1]
           and _rows_fm[-1]["grids"] == [3, 3, 1])
        ck("ingest: Δ 是 meV/原子口径 (k=3 vs k=6 手算 0.0694 meV/原子)",
           abs(_rows_fm[-1]["delta_meV_per_atom"] - 0.0694444) < 1e-4,
           f"{_rows_fm[-1]['delta_meV_per_atom']:.6f}")

        write_all(lambda g, nm: -100.0 + 5.0 / (g[0] ** 2))
        res = ingest(root, grids, tol=1.0, min_kb=0, quiet=True)
        ck("ingest: 未收敛序列 -> status=not_converged", res["status"] == "not_converged", res["status"])
        ck("ingest: 未收敛时仍给出已合格的粗网格", res["recommended_kpoints"].get("fm") == [6, 6, 1],
           str(res["recommended_kpoints"].get("fm")))

        # D. 闸门
        write_all(lambda g, nm: -100.0)
        res = ingest(root, grids, tol=1.0, min_kb=200, quiet=True)
        ck("ingest: 体积闸门 (--min-kb 200 对合成小文件) -> 不参与比较",
           res["status"] in ("insufficient", "partial"), res["status"])
        write_all(lambda g, nm: -100.0, scf=False)
        res = ingest(root, grids, tol=1.0, min_kb=0, quiet=True)
        ck("ingest: SCF 未收敛 -> 不参与比较", res["status"] == "insufficient", res["status"])
        write_all(lambda g, nm: -100.0)
        os.remove(os.path.join(root, "kconv", "runs", "fm@k4_4_1", "OUTCAR"))
        res = ingest(root, grids, tol=1.0, min_kb=0, quiet=True)
        ck("ingest: 缺一个 OUTCAR 仍能判定 (missing 行记 missing_outcar)",
           res["status"] == "ok"
           and any(r["status"] == "missing_outcar" for r in res["orders"]["fm"]["rows"]), res["status"])
        write_all(lambda g, nm: -100.0)
        write_all(lambda g, nm: -100.0, counts=lambda g: "6 6" if g != (4, 4, 1) else "6 8")
        res = ingest(root, grids, tol=1.0, min_kb=0, quiet=True)
        ck("ingest: 原子数不一致 -> cell_mismatch",
           res["orders"]["fm"]["status"] == "cell_mismatch", res["orders"]["fm"]["status"])
        ck("ingest: 明细 JSON 落盘", os.path.isfile(os.path.join(root, "kconv", "k_convergence.json")))
        ck("ingest: 结论里带免责说明 (J/T_N 不算)",
           "J 是能量差" in res["note"] and "不" in res["note"])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print(f"\n--selftest: {'ALL PASS' if not fail else 'FAIL'} ({ok}/{ok + fail} 项)")
    return 0 if not fail else 1


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="k 点收敛抽查 (P4.2): 生成 KPOINTS 变体 + 回收判定")
    ap.add_argument("--root", default=None, help="项目目录 (含 runs/<order>/)")
    ap.add_argument("--grids", default=DEFAULT_GRIDS,
                    help=f"k 网格: '3,4,6,8' (面内 k, kz=1) 或 '3 3 1;4 4 1' (默认 {DEFAULT_GRIDS})")
    ap.add_argument("--tol", type=float, default=1.0, help="判据: |ΔE| <= 该值 (meV/原子, 默认 1.0)")
    ap.add_argument("--orders", nargs="+", default=None, help="只抽查这些序 (默认 runs/ 下全部)")
    ap.add_argument("--min-kb", type=int, default=200, help="OUTCAR 最小体积 KB (默认 200)")
    ap.add_argument("--plan-only", action="store_true", help="只生成网格变体, 不判定")
    ap.add_argument("--ingest-only", action="store_true", help="只判定 (不重新生成)")
    ap.add_argument("--json", action="store_true", help="把结果 JSON 打到 stdout")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)

    if args.selftest:
        return _selftest()
    if not args.root:
        raise SystemExit("需要 --root <项目目录> (含 runs/<order>/); 例: "
                         "python scripts/k_convergence.py --root ./Fe3Sn_pipeline --grids 3,4,6,8")

    root = os.path.abspath(args.root)
    grids = parse_grids(args.grids)
    if args.tol <= 0:
        raise SystemExit(f"--tol {args.tol} 必须为正 (meV/原子)")
    if not args.ingest_only:
        plan(root, grids, orders=args.orders, quiet=args.quiet)
    if args.plan_only:
        return 0
    res = ingest(root, grids, orders=args.orders, tol=args.tol, min_kb=args.min_kb, quiet=args.quiet)
    if args.json:
        print(json.dumps(res, ensure_ascii=False, indent=2))
    return 1 if res["status"] in ("not_converged", "cell_mismatch", "partial") else 0


if __name__ == "__main__":
    sys.exit(main())
