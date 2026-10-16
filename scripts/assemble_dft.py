#!/usr/bin/env python
"""assemble_dft.py
从超胞磁序 cases 生成 VASP 输入: INCAR.static / INCAR.relax / KPOINTS / POTCAR.spec,
并把每个磁态的 MAGMOM 写入对应 INCAR, 组装成可直接提交的 runs/<order>/ 目录。

两种自旋模式（P3.16 起统一入口）:
  * 默认（共线）: `ISPIN = 2` + 1 分量 MAGMOM。
  * `--noncollinear`: `LNONCOLLINEAR = .TRUE.` + **3 分量 MAGMOM**（**绝不写 `ISPIN`**）;
    加 `--lsorbit` 再开 `LSORBIT = .TRUE.` + `GGA_COMPAT = .FALSE.`（SOC，弱铁磁倾斜/
    DMI/kagome 所需）; `--saxis x y z` 可指定自旋量子化轴。

不静默降级（NEW_TYPES_PLAN §115）——下面三条都是**硬报错**，不是警告:
  1. 磁序是 3 分量（非共线）却没给 `--noncollinear` ⇒ 拒绝按共线近似继续;
  2. 磁序需要 SOC（`cant`/`tilt`/`spiral` 系，见 `mag_vectors.needs_soc`）却没给
     `--lsorbit` ⇒ 拒绝算（没有 SOC 就没有弱铁磁倾斜这回事）;
  3. `--lsorbit` 而没给 `--noncollinear` ⇒ 报错（VASP 的 LSORBIT 依赖 LNONCOLLINEAR）。

元素顺序: `LDAUL/LDAUU/LDAUJ`、`POTCAR.spec`、`MAGMOM` 一律按 **POSCAR 的元素块顺序**
渲染（P2 抓到的真 bug: 磁性元素并不总在第一个块，CrOCl 就是 Cl Cr O）。`--species`
只用于一致性检查，不改变顺序。

角度规则全部来自 `mag_vectors.py`（P3.15 的单一来源）; 能量口径见 `dft_energy.py`。

用法:
  python assemble_dft.py --root <proj> --element Mn --nonmag P S --U 5.0 --moment 5.0 --encut 520
  python assemble_dft.py --root <proj> --element Cr --nonmag O Cl --noncollinear --lsorbit
  python assemble_dft.py --selftest
"""
import argparse
import json
import os
import shutil
import sys

import numpy as np

import mag_vectors as MV
from element_data import ldaul_for

# 元素 -> 建议 PAW（POTCAR.spec 用）
PAW_HINT = {"Mn": "Mn_pv (含 3s3p 半芯)", "Cr": "Cr_pv", "V": "V_pv", "Fe": "Fe_pv",
            "Co": "Co", "Ni": "Ni_pv", "P": "P", "S": "S", "Se": "Se", "Te": "Te",
            "Cl": "Cl", "Br": "Br", "I": "I"}

# 每个 INCAR 变体的"体"（不含自旋块与 MAGMOM），数据驱动，避免复制两份模板
INCAR_BODY = {
    "static": """PREC     = Normal
ENCUT    = {encut}
EDIFF    = 1E-6
ALGO     = Normal
NELM     = 200
ISMEAR   = 0
SIGMA    = 0.05
""",
    "relax": """PREC     = Normal
ENCUT    = {encut}
EDIFF    = 1E-5
EDIFFG   = -0.02
IBRION   = 2
NSW      = 300
NELM     = 200
ISIF     = 2
ISYM     = 0
ISMEAR   = 0
SIGMA    = 0.05
""",
}

# k 网格默认值。P4.2 起网格可以显式给（`--kpoints`），且 header 必须写**实际**网格：
# 旧 header 写死了 "for 2x2x1 monolayer supercell"，既是对超胞尺寸的错误陈述，也掩盖了
# 「这个网格从没做过收敛抽查」这件事。
KPOINTS_MESH = (3, 3, 1)


def kpoints_text(mesh):
    """Gamma 心 k 网格文本（monolayer + 真空沿 c）。"""
    kx, ky, kz = (int(x) for x in mesh)
    return (f"Gamma-centered k-mesh {kx} {ky} {kz} (vacuum along c)\n"
            f"{kx} {ky} {kz}\n"
            "0 0 0\n")


# ---------------------------------------------------------------------------
# POSCAR / 元素
# ---------------------------------------------------------------------------
def read_species_counts(poscar_path):
    """-> (species, counts)。同时给出 read_poscar_head 的 (species, counts, n_atoms)。"""
    lines = [ln for ln in open(poscar_path, encoding="utf-8").read().splitlines() if ln.strip()]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    return species, counts


def read_poscar_head(poscar_path):
    species, counts = read_species_counts(poscar_path)
    return species, counts, sum(counts)


def magnetic_block(species, counts, element):
    """磁性元素在 POSCAR 原子序列里的**下标区间**（按元素块顺序，不是"排在最前"）。

    这是 P2 抓到的一类真 bug 的根: MAGMOM 必须按 POSCAR 的元素块顺序排，
    而 CrOCl 这类体系磁元素并不在第一个块 (Cl Cr O)。
    """
    off = 0
    for s, c in zip(species, counts):
        if s == element:
            return off, off + c
        off += c
    raise SystemExit(f"元素 {element!r} 不在 POSCAR 的 {species}")


def mag_idx_of(species, counts, element):
    lo, hi = magnetic_block(species, counts, element)
    return list(range(lo, hi))


def build_potcar_spec(species):
    lines = ["# POTCAR 组装说明", "", "元素顺序 **必须与 POSCAR/magmom.txt 一致**:  " + " ".join(species), ""]
    lines.append("| 顺序 | 元素 | 建议 PAW | 说明 |")
    lines.append("|---|---|---|---|")
    for i, s in enumerate(species):
        lines.append(f"| {i + 1} | {s} | {PAW_HINT.get(s, s)} | |")
    lines.append("")
    lines.append("实际使用的 PAW 版本/日期务必回填。")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# INCAR 渲染
# ---------------------------------------------------------------------------
def spin_block(noncollinear, lsorbit, saxis, magmom_line):
    """自旋相关的那几行。共线与非共线**互斥**（非共线块里绝不出现 ISPIN）。"""
    if noncollinear:
        out = ["LNONCOLLINEAR = .TRUE."]
        if lsorbit:
            out.append("LSORBIT  = .TRUE.")
            out.append("GGA_COMPAT = .FALSE.      ! SOC 必须关掉 GGA 的兼容模式")
        if saxis is not None:
            out.append("SAXIS    = {:.6f} {:.6f} {:.6f}".format(*[float(x) for x in saxis]))
        out.append(magmom_line)
        return "\n".join(out) + "\n"
    return "ISPIN    = 2\n" + magmom_line + "\n"


def _with_isym(body, isym):
    """显式指定 ISYM 时覆盖/插入（缺省 isym=None 一律不动模板）。

    动机：relax 模板自带 `ISYM = 0`，而 static 模板**没有** ISYM（VASP 默认 2）。
    对磁序比较敏感的计算，用户要能显式关掉对称性而不必手改 INCAR。
    """
    if isym is None:
        return body
    want = f"ISYM     = {int(isym)}"
    lines = body.splitlines()
    if any(ln.strip().upper().startswith("ISYM") for ln in lines):
        return "\n".join(want if ln.strip().upper().startswith("ISYM") else ln
                         for ln in lines) + "\n"
    return body.rstrip("\n") + "\n" + want + "\n"


def build_incar(kind, species, counts, element, U, moment, encut, magmom_line,
                noncollinear=False, lsorbit=False, saxis=None, npar=4, isym=None):
    """渲染一个 INCAR（kind ∈ {'static','relax'}）。"""
    if kind not in INCAR_BODY:
        raise ValueError(f"未知 INCAR 变体 {kind!r}")
    ldaul = " ".join(ldaul_for(s) if s == element else "-1" for s in species)
    ldauu = " ".join(f"{U:.1f}" if s == element else "0.0" for s in species)
    ldauj = " ".join("0.0" for _ in species)
    tag = "noncollinear-SOC" if (noncollinear and lsorbit) else \
          ("noncollinear" if noncollinear else "collinear")
    system = " ".join(species)
    return (f"SYSTEM   = {system} monolayer magnetic order ({tag})\n"
            + _with_isym(INCAR_BODY[kind].format(encut=encut), isym)
            + "\n"
            + spin_block(noncollinear, lsorbit, saxis, magmom_line)
            + "\nLDAU     = .TRUE.\n"
            + "LDAUTYPE = 2\n"
            + f"LDAUL    = {ldaul}\n"
            + f"LDAUU    = {ldauu}\n"
            + f"LDAUJ    = {ldauj}\n"
            + "LMAXMIX  = 4\n\n"
            + f"LREAL    = {'Auto' if kind == 'relax' else '.FALSE.'}\n"
            + "LWAVE    = .FALSE.\n"
            + "LCHARG   = .FALSE.\n"
            + f"NPAR     = {npar}\n")


def default_magmom_line(species, counts, element, moment, noncollinear):
    """模板里的默认 MAGMOM 行（**按 POSCAR 元素块顺序**；非共线写全 3N 个值）。"""
    n_atoms = sum(counts)
    mag_idx = mag_idx_of(species, counts, element)
    if noncollinear:
        spins = np.zeros((n_atoms, 3))
        spins[mag_idx, 2] = 1.0
        return MV.render_magmom(spins * float(moment))
    vals = np.zeros(n_atoms)
    vals[mag_idx] = float(moment)
    return "MAGMOM   = " + " ".join(f"{v:.1f}" for v in vals)


# ---------------------------------------------------------------------------
# 每个磁序的 MAGMOM
# ---------------------------------------------------------------------------
def magmom_line_for_case(case_dir, n_atoms, noncollinear, order, warn=None):
    """读 cases/<order>/magmom.txt -> 可直接写进 INCAR 的 MAGMOM 行。

    共线模式遇到 3 分量 -> 报错（不静默按共线近似继续）;
    非共线模式遇到 1 分量 -> 按"自旋沿 z"展开并**显式提示**（共线磁态在非共线设置下
    计算的原意，不是降级）。
    """
    path = os.path.join(case_dir, "magmom.txt")
    if not os.path.exists(path):
        raise SystemExit(f"磁序 {order!r} 缺 magmom.txt ({path})")
    text = open(path, encoding="utf-8").read().strip()
    vals, ncomp = MV.parse_magmom(text, n_atoms=n_atoms)
    if ncomp == 0 or vals.size == 0:
        raise SystemExit(f"磁序 {order!r} 的 magmom.txt 解析不出数值 ({path})")
    if noncollinear:
        if ncomp == 3:
            line = text if text.lstrip().upper().startswith("MAGMOM") else MV.render_magmom(vals.reshape(-1, 3))
        else:
            m = np.asarray(vals, dtype=float).reshape(n_atoms, 1)
            line = MV.render_magmom(np.hstack([np.zeros((n_atoms, 2)), m]))
            if warn:
                warn(f"[注意] {order}: magmom.txt 是共线 (1 分量) -> 非共线设置下按自旋沿 z 展开")
        ok, _nc, probs = MV.check_magmom(line, n_atoms, want_ncomp=3)
        if not ok:
            raise SystemExit(f"磁序 {order!r} 的 MAGMOM 不合法: {'; '.join(probs)}")
        return line
    if ncomp == 3:
        raise SystemExit(
            f"磁序 {order!r} 的 magmom.txt 是 3 分量 (非共线), 但当前是共线模式 (ISPIN=2) "
            f"=> 拒绝按共线近似继续; 请加 --noncollinear"
            + (" --lsorbit" if MV.needs_soc(order) else ""))
    return text if text.lstrip().upper().startswith("MAGMOM") else MV.render_magmom(np.asarray(vals, dtype=float).reshape(-1, 1))


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="生成 VASP 输入 (共线 / 非共线+SOC)")
    ap.add_argument("--root", required=True)
    ap.add_argument("--element", required=True)
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--species", nargs="+", default=None,
                    help="元素顺序; 只做一致性检查, 渲染一律按 POSCAR 顺序 (推荐不给)")
    ap.add_argument("--U", type=float, default=5.0)
    ap.add_argument("--moment", type=float, default=5.0)
    ap.add_argument("--encut", type=float, default=520)
    ap.add_argument("--kpoints", type=int, nargs=3, default=list(KPOINTS_MESH),
                    help=f"Gamma 心 k 网格 (默认 {' '.join(str(x) for x in KPOINTS_MESH)}, 真空沿 c)")
    ap.add_argument("--npar", type=int, default=4, help="NPAR (默认 4)")
    ap.add_argument("--isym", type=int, default=None,
                    help="显式写 ISYM (0 = 关对称性); 缺省: static 不写, relax 用模板的 0")
    ap.add_argument("--cases", default="02_supercell/cases")
    ap.add_argument("--noncollinear", action="store_true", help="LNONCOLLINEAR + 3 分量 MAGMOM")
    ap.add_argument("--lsorbit", action="store_true", help="SOC (需配合 --noncollinear)")
    ap.add_argument("--saxis", nargs=3, type=float, default=None, help="自旋量子化轴, 如 0 0 1")
    ap.add_argument("--all-orders", action="store_true",
                    help="组装 cases/ 下**所有**磁序目录 (默认只用 magnetic_orders.json 的 selected_orders)")
    ap.add_argument("--selftest", action="store_true")
    _argv = list(sys.argv[1:] if argv is None else argv)
    if "--selftest" in _argv:          # 自检不需要 --root/--element, 先短路
        return _selftest()
    args = ap.parse_args(_argv)

    if args.lsorbit and not args.noncollinear:
        raise SystemExit("--lsorbit 需要 --noncollinear (VASP 的 LSORBIT 依赖 LNONCOLLINEAR)")
    if any(int(k) < 1 for k in args.kpoints):
        raise SystemExit(f"--kpoints 三个分量都必须 >= 1, 实际 {list(args.kpoints)}")
    if int(args.npar) < 1:
        raise SystemExit(f"--npar 必须 >= 1, 实际 {args.npar}")

    explicit_kmesh = tuple(int(x) for x in args.kpoints) != tuple(KPOINTS_MESH)

    root = os.path.abspath(args.root)
    cases_dir = os.path.join(root, args.cases)
    fm_poscar = os.path.join(cases_dir, "fm", "POSCAR")
    if not os.path.exists(fm_poscar):
        orders_any = sorted(d for d in os.listdir(cases_dir)
                            if os.path.isfile(os.path.join(cases_dir, d, "POSCAR"))) if os.path.isdir(cases_dir) else []
        if not orders_any:
            raise SystemExit(f"{cases_dir} 下没有可用磁序目录")
        fm_poscar = os.path.join(cases_dir, orders_any[0], "POSCAR")
    sp_poscar, counts, n_atoms = read_poscar_head(fm_poscar)
    if args.species and list(args.species) != sp_poscar:
        if sorted(args.species) != sorted(sp_poscar):
            raise SystemExit(f"--species {args.species} 与 POSCAR 的 {sp_poscar} 不是同一组元素")
        print(f"[注意] --species {args.species} 与 POSCAR 顺序 {sp_poscar} 不同 "
              f"-> LDAU/POTCAR/MAGMOM 一律按 POSCAR 顺序渲染")
    species = list(sp_poscar)

    # 磁序列表: 默认用 magnetic_orders.json 的 selected_orders; --all-orders 用全部磁序目录
    if args.all_orders:
        orders = sorted(d for d in os.listdir(cases_dir)
                        if os.path.isfile(os.path.join(cases_dir, d, "POSCAR"))) \
            if os.path.isdir(cases_dir) else []
    else:
        try:
            meta = json.load(open(os.path.join(cases_dir, "magnetic_orders.json"), encoding="utf-8"))
            orders = list(meta.get("selected_orders") or meta.get("coefficients_per_cell", {}))
        except (OSError, ValueError):
            orders = ["fm", "neel", "zigzag", "stripy"]

    # §115-2: 需要 SOC 的序必须显式开 LSORBIT
    if not args.lsorbit:
        soc_orders = [o for o in orders if MV.needs_soc(o)]
        if soc_orders:
            raise SystemExit(
                f"磁序 {soc_orders} 需要 SOC (弱铁磁倾斜/螺旋), 但没给 --lsorbit => 拒绝按共线近似继续 "
                f"(NEW_TYPES_PLAN §115)")

    default_line = default_magmom_line(species, counts, args.element, args.moment, args.noncollinear)

    dft_dir = os.path.join(root, "03_dft")
    os.makedirs(dft_dir, exist_ok=True)
    for kind, fname in (("static", "INCAR.static"), ("relax", "INCAR.relax")):
        text = build_incar(kind, species, counts, args.element, args.U, args.moment,
                           args.encut, default_line, args.noncollinear, args.lsorbit, args.saxis,
                           npar=args.npar, isym=args.isym)
        with open(os.path.join(dft_dir, fname), "w", encoding="utf-8") as f:
            f.write(text)
    with open(os.path.join(dft_dir, "KPOINTS"), "w", encoding="utf-8") as f:
        f.write(kpoints_text(args.kpoints))
    with open(os.path.join(dft_dir, "POTCAR.spec"), "w", encoding="utf-8") as f:
        f.write(build_potcar_spec(species))

    runs_dir = os.path.join(root, "runs")
    os.makedirs(runs_dir, exist_ok=True)
    n_written = 0
    for order in orders:
        d = os.path.join(cases_dir, order)
        if not os.path.isdir(d):
            continue
        line = magmom_line_for_case(d, n_atoms, args.noncollinear, order,
                                    warn=lambda m: print(m))
        incar = build_incar("static", species, counts, args.element, args.U, args.moment,
                            args.encut, line, args.noncollinear, args.lsorbit, args.saxis,
                            npar=args.npar, isym=args.isym)
        rd = os.path.join(runs_dir, order)
        os.makedirs(rd, exist_ok=True)
        with open(os.path.join(rd, "INCAR"), "w", encoding="utf-8") as f:
            f.write(incar)
        shutil.copy(os.path.join(dft_dir, "KPOINTS"), os.path.join(rd, "KPOINTS"))
        shutil.copy(os.path.join(d, "POSCAR"), os.path.join(rd, "POSCAR"))
        n_written += 1

    mode = ("非共线" + ("+SOC (LSORBIT)" if args.lsorbit else "")) if args.noncollinear else "共线 (ISPIN=2)"
    print(f"写出 03_dft/ (INCAR.static/.relax, KPOINTS, POTCAR.spec) 与 runs/<order>/ 组装目录")
    print(f"k 网格: {' '.join(str(x) for x in args.kpoints)}"
          + ("" if explicit_kmesh else " (默认值, **未做收敛抽查**)")
          + f"; NPAR = {args.npar}"
          + (f"; ISYM = {args.isym}" if args.isym is not None else
             " (static 不写 ISYM, relax = 0)"))
    if not explicit_kmesh and n_atoms >= 60:
        print(f"[提示] 体系 {n_atoms} 个原子却用默认 k 网格 "
              f"{' '.join(str(x) for x in args.kpoints)}: 上真实机时前先做一次 k 收敛抽查 "
              f"(P4_2_KAGOME_RUNBOOK.md §1), 或用 --kpoints 指定")
    print(f"自旋模式: {mode}; 磁序 {n_written} 个: {orders}")
    print(f"MAGMOM 分量数: {3 if args.noncollinear else 1} / 原子; 原子总数 {n_atoms}"
          + (f"; SAXIS {args.saxis}" if args.saxis else ""))
    print("POTCAR 顺序:", " ".join(species), "(与 POSCAR/magmom.txt 一致)")
    if args.noncollinear and not args.lsorbit:
        print("[提示] 未开 SOC: 非共线几何阻挫 (120 度等) 可算; 弱铁磁倾斜/DMI 需 --lsorbit")
    return 0


# ---------------------------------------------------------------------------
# selftest（无 DFT: 只渲染 INCAR 文本与 MAGMOM，逐条查）
# ---------------------------------------------------------------------------
# 故意让磁性元素**不在第一个元素块**（CrOCl 型: Cl Cr O）—— 这是 P2 抓到的
# "MAGMOM 与 POSCAR 元素分组错位" bug 的回归位。
_SPECIES = ["Cl", "Cr", "O"]
_COUNTS = [2, 2, 2]
_POSCAR_TXT = ("probe\n1.0\n3.2 0 0\n0 3.9 0\n0 0 20\nCl Cr O\n2 2 2\nDirect\n"
               + "".join(f"{0.1 * i:.6f} {0.1 * i:.6f} 0.5\n" for i in range(6)))


def _selftest():
    checks = []

    def ck(name, cond, extra=""):
        checks.append(bool(cond))
        print(f"  [{'ok ' if cond else 'FAIL'}] {name}" + (f"  {extra}" if extra else ""))

    base = dict(species=_SPECIES, counts=_COUNTS, element="Cr", U=4.0, moment=3.0, encut=520)

    d_col = default_magmom_line(_SPECIES, _COUNTS, "Cr", 3.0, False)
    vals_c, nc_c = MV.parse_magmom(d_col, n_atoms=6)
    ck("共线默认 MAGMOM: 1 分量/原子, 共 6 个值", nc_c == 1 and vals_c.size == 6)
    ck("共线默认 MAGMOM: 按 POSCAR 元素块顺序 (Cl 块=0, Cr 块=3, O 块=0)",
       np.allclose(vals_c[:2], 0.0) and np.allclose(vals_c[2:4], 3.0) and np.allclose(vals_c[4:], 0.0),
       f"vals={np.round(vals_c, 3).tolist()}")

    d_nc = default_magmom_line(_SPECIES, _COUNTS, "Cr", 3.0, True)
    v_nc, nc_nc = MV.parse_magmom(d_nc, n_atoms=6)
    m_nc = v_nc.reshape(6, 3)
    ck("非共线默认 MAGMOM: 3 分量/原子, 共 18 个值", nc_nc == 3 and v_nc.size == 18)
    ck("非共线默认 MAGMOM: 磁性块沿 z = moment, 其余为 0",
       np.allclose(m_nc[2:4], [0, 0, 3.0]) and np.allclose(np.delete(m_nc, [2, 3], axis=0), 0.0))

    inc_col = build_incar(kind="static", **base, magmom_line=d_col)
    ck("共线 INCAR: 有 ISPIN=2, 无 LNONCOLLINEAR/LSORBIT",
       "ISPIN    = 2" in inc_col and "LNONCOLLINEAR" not in inc_col and "LSORBIT" not in inc_col)
    ck("共线 INCAR: LDAUL 按 POSCAR 顺序 (Cl Cr O -> -1 2 -1)",
       "LDAUL    = -1 2 -1" in inc_col)

    inc_nc = build_incar(kind="static", **base, magmom_line=d_nc, noncollinear=True)
    ck("非共线 INCAR: 有 LNONCOLLINEAR = .TRUE. 且**不含 ISPIN**",
       "LNONCOLLINEAR = .TRUE." in inc_nc and "ISPIN" not in inc_nc)
    ck("非共线 INCAR: 未开 --lsorbit 时无 LSORBIT/GGA_COMPAT",
       "LSORBIT" not in inc_nc and "GGA_COMPAT" not in inc_nc)

    inc_soc = build_incar(kind="static", **base, magmom_line=d_nc, noncollinear=True,
                          lsorbit=True, saxis=(0.0, 0.0, 1.0))
    ck("SOC INCAR: LSORBIT + GGA_COMPAT=.FALSE. + SAXIS, 仍无 ISPIN",
       "LSORBIT  = .TRUE." in inc_soc and "GGA_COMPAT = .FALSE." in inc_soc
       and "SAXIS    = 0.000000 0.000000 1.000000" in inc_soc and "ISPIN" not in inc_soc)
    ck("SOC INCAR: SYSTEM 标注 noncollinear-SOC", "(noncollinear-SOC)" in inc_soc.splitlines()[0])

    inc_relax = build_incar(kind="relax", **base, magmom_line=d_nc, noncollinear=True)
    ck("relax 变体: 有 IBRION/NSW/EDIFFG 且非共线块完整",
       "IBRION   = 2" in inc_relax and "EDIFFG" in inc_relax and "LNONCOLLINEAR" in inc_relax)

    # --- k 网格 / NPAR / ISYM（P4.2 参数化）---
    _kp_def = kpoints_text(KPOINTS_MESH)
    ck("k 网格: 默认 header 写实际网格 (不再写死 2x2x1)",
       _kp_def.splitlines()[0] == "Gamma-centered k-mesh 3 3 1 (vacuum along c)"
       and _kp_def.splitlines()[1] == "3 3 1")
    _kp6 = kpoints_text((6, 6, 1))
    ck("k 网格: --kpoints 6 6 1 -> header 与网格行同步",
       _kp6.splitlines()[0].endswith("6 6 1 (vacuum along c)")
       and _kp6.splitlines()[1] == "6 6 1" and _kp6.splitlines()[2] == "0 0 0")
    _static_def = build_incar(kind="static", **base, magmom_line=d_col)
    _relax_def = build_incar(kind="relax", **base, magmom_line=d_col)
    ck("ISYM 缺省: static 不写 ISYM, relax 保持模板的 ISYM=0",
       "ISYM" not in _static_def and "ISYM     = 0" in _relax_def)
    ck("--isym 0: static 里插入 ISYM=0",
       "ISYM     = 0" in build_incar(kind="static", **base, magmom_line=d_col, isym=0))
    _relax_i2 = build_incar(kind="relax", **base, magmom_line=d_col, isym=2)
    ck("--isym 2: relax 的模板值被覆盖且只有一条 ISYM",
       _relax_i2.count("ISYM") == 1 and "ISYM     = 2" in _relax_i2)
    ck("--npar 8: INCAR 落 NPAR = 8",
       "NPAR     = 8" in build_incar(kind="static", **base, magmom_line=d_col, npar=8))

    import tempfile
    td = tempfile.mkdtemp(prefix="asmdft_")
    try:
        cases = os.path.join(td, "02_supercell", "cases")
        spins = {"fm": MV.magmom_values(MV.spin_120deg(np.array([0, 1])), 3.0, [2, 3], 6),
                 "fm_col": MV.magmom_values(MV.spin_collinear([1, -1]), 3.0, [2, 3], 6),
                 "cant_a": MV.magmom_values(MV.spin_cant([1, -1], 8.0), 3.0, [2, 3], 6)}
        for name, sp in spins.items():
            d = os.path.join(cases, name)
            os.makedirs(d, exist_ok=True)
            with open(os.path.join(d, "magmom.txt"), "w", encoding="utf-8") as f:
                f.write(MV.render_magmom(sp))
            with open(os.path.join(d, "POSCAR"), "w", encoding="utf-8") as f:
                f.write(_POSCAR_TXT)

        def set_orders(orders):
            with open(os.path.join(cases, "magnetic_orders.json"), "w", encoding="utf-8") as f:
                json.dump({"selected_orders": orders, "n_mag": 2, "n_atoms": 6, "moment": 3.0}, f)

        def try_main(extra):
            try:
                main(["--root", td, "--element", "Cr", "--nonmag", "O", "Cl",
                      "--U", "4.0", "--moment", "3.0"] + extra)
                return None
            except SystemExit as e:
                return str(e)

        set_orders(["fm"])
        msg = try_main([])
        ck("硬报错 1: 3 分量磁序 + 共线模式 -> 拒绝 (提示 --noncollinear)",
           msg is not None and "--noncollinear" in msg, (msg or "")[:64])

        set_orders(["fm_col", "cant_a"])
        msg = try_main(["--noncollinear"])
        ck("硬报错 2: cant 序 + 无 --lsorbit -> 拒绝",
           msg is not None and "SOC" in msg and "--lsorbit" in msg, (msg or "")[:64])

        msg = try_main(["--lsorbit"])
        ck("硬报错 3: --lsorbit 而没 --noncollinear -> 拒绝",
           msg is not None and "--noncollinear" in msg, (msg or "")[:64])

        rc = main(["--root", td, "--element", "Cr", "--nonmag", "O", "Cl", "--U", "4.0",
                   "--moment", "3.0", "--noncollinear", "--lsorbit", "--saxis", "0", "0", "1"])
        inc = open(os.path.join(td, "runs", "cant_a", "INCAR"), encoding="utf-8").read()
        ck("全链路: --noncollinear --lsorbit 组装成功, runs/<order>/INCAR 带 SOC 且无 ISPIN",
           rc == 0 and "LSORBIT  = .TRUE." in inc and "ISPIN" not in inc and "SAXIS" in inc)
        ck("共线 magmom.txt 在非共线模式下按自旋沿 z 展开 (3 分量)",
           MV.parse_magmom(inc, n_atoms=6)[1] == 3)
        ck("非共线模式下 runs/ 覆盖全部选定序",
           sorted(os.listdir(os.path.join(td, "runs"))) == ["cant_a", "fm_col"])
        ck("03_dft/ 四个文件齐全",
           all(os.path.exists(os.path.join(td, "03_dft", n)) for n in
               ("INCAR.static", "INCAR.relax", "KPOINTS", "POTCAR.spec")))

        # --- --kpoints / --npar / --isym 的全链路落盘 ---
        rc2 = main(["--root", td, "--element", "Cr", "--nonmag", "O", "Cl", "--U", "4.0",
                    "--moment", "3.0", "--noncollinear", "--lsorbit", "--saxis", "0", "0", "1",
                    "--kpoints", "5", "5", "1", "--npar", "8", "--isym", "0"])
        kp = open(os.path.join(td, "03_dft", "KPOINTS"), encoding="utf-8").read()
        inc8 = open(os.path.join(td, "runs", "cant_a", "INCAR"), encoding="utf-8").read()
        ck("全链路: --kpoints 5 5 1 落进 03_dft/KPOINTS (header 同步)",
           rc2 == 0 and kp.splitlines()[1].strip() == "5 5 1" and "5 5 1" in kp.splitlines()[0])
        ck("全链路: --npar 8 / --isym 0 落进 runs/<order>/INCAR",
           "NPAR     = 8" in inc8 and "ISYM     = 0" in inc8)
        ck("全链路: runs/<order>/KPOINTS 与 03_dft/KPOINTS 逐字一致",
           open(os.path.join(td, "runs", "cant_a", "KPOINTS"), encoding="utf-8").read() == kp)
        ck("参数域: --kpoints 0 3 1 被硬拒",
           try_main(["--noncollinear", "--lsorbit", "--kpoints", "0", "3", "1"]) is not None)
        ck("参数域: --npar 0 被硬拒",
           try_main(["--noncollinear", "--lsorbit", "--npar", "0"]) is not None)
    finally:
        shutil.rmtree(td, ignore_errors=True)

    n_ok = sum(checks)
    print(f"\n--selftest: {'ALL PASS' if n_ok == len(checks) else 'FAIL'} ({n_ok}/{len(checks)} 项)")
    return 0 if n_ok == len(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main())
