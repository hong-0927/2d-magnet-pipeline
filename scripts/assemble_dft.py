#!/usr/bin/env python
"""assemble_dft.py
从超胞磁序 cases 生成 VASP 输入: INCAR.static / INCAR.relax / KPOINTS / POTCAR.spec,
并把每个磁态的 MAGMOM 写入对应 INCAR, 组装成可直接提交的 runs/<order>/ 目录。

用法:
  python assemble_dft.py --root <proj> --element Mn --nonmag P S --U 5.0 --moment 5.0 \
      --species "Mn P S" --encut 520
"""
import argparse
import json
import os
import shutil

import numpy as np

from element_data import ldaul_for


def read_species_counts(poscar_path):
    lines = [ln for ln in open(poscar_path, encoding="utf-8").read().splitlines() if ln.strip()]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    return species, counts


def build_incar_static(species, counts, element, U, moment, encut, npar=4):
    total_mag = counts[species.index(element)] if element in species else 0
    n_nonmag = sum(counts) - total_mag
    magmom_default = f"{total_mag}*{moment:.1f} {n_nonmag}*0.0"
    ldaul = " ".join(ldaul_for(s) if s == element else "-1" for s in species)
    ldauu = " ".join(f"{U:.1f}" if s == element else "0.0" for s in species)
    ldauj = " ".join("0.0" for _ in species)
    return f"""SYSTEM   = {" ".join(species)} monolayer magnetic order
PREC     = Normal
ENCUT    = {encut}
EDIFF    = 1E-6
ALGO     = Normal
NELM     = 200
ISMEAR   = 0
SIGMA    = 0.05

ISPIN    = 2
MAGMOM   = {magmom_default}     ! <-- 各磁态用 runs/<order>/INCAR 里的实际值

LDAU     = .TRUE.
LDAUTYPE = 2
LDAUL    = {ldaul}
LDAUU    = {ldauu}
LDAUJ    = {ldauj}
LMAXMIX  = 4

LREAL    = .FALSE.
LWAVE    = .FALSE.
LCHARG   = .FALSE.
NPAR     = {npar}
"""


def build_incar_relax(species, counts, element, U, moment, encut, npar=4):
    total_mag = counts[species.index(element)] if element in species else 0
    n_nonmag = sum(counts) - total_mag
    magmom_default = f"{total_mag}*{moment:.1f} {n_nonmag}*0.0"
    ldaul = " ".join(ldaul_for(s) if s == element else "-1" for s in species)
    ldauu = " ".join(f"{U:.1f}" if s == element else "0.0" for s in species)
    ldauj = " ".join("0.0" for _ in species)
    return f"""SYSTEM   = {" ".join(species)} monolayer relax (FM)
PREC     = Normal
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

ISPIN    = 2
MAGMOM   = {magmom_default}

LDAU     = .TRUE.
LDAUTYPE = 2
LDAUL    = {ldaul}
LDAUU    = {ldauu}
LDAUJ    = {ldauj}
LMAXMIX  = 4

LREAL    = Auto
LWAVE    = .FALSE.
LCHARG   = .FALSE.
NPAR     = {npar}
"""


KPOINTS_TMPL = """Gamma-centered k-mesh for 2x2x1 monolayer supercell (vacuum along c)
3 3 1
0 0 0
"""


def build_potcar_spec(species):
    lines = ["# POTCAR 组装说明", "", "元素顺序 **必须与 POSCAR/magmom.txt 一致**:  " + " ".join(species), ""]
    lines.append("| 顺序 | 元素 | 建议 PAW | 说明 |")
    lines.append("|---|---|---|---|")
    pa = {"Mn": "Mn_pv (含 3s3p 半芯)", "Cr": "Cr_pv", "V": "V_pv", "Fe": "Fe_pv",
          "Co": "Co", "Ni": "Ni_pv", "P": "P", "S": "S", "Se": "Se", "Te": "Te",
          "Cl": "Cl", "Br": "Br", "I": "I"}
    for i, s in enumerate(species):
        lines.append(f"| {i + 1} | {s} | {pa.get(s, s)} | |")
    lines.append("")
    lines.append("实际使用的 PAW 版本/日期务必回填。")
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--element", required=True)
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--species", nargs="+", required=True, help="元素顺序, 如 Mn P S")
    ap.add_argument("--U", type=float, default=5.0)
    ap.add_argument("--moment", type=float, default=5.0)
    ap.add_argument("--encut", type=float, default=520)
    ap.add_argument("--cases", default="02_supercell/cases")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    species = args.species
    cases_dir = os.path.join(root, args.cases)
    # 从 fm 的 POSCAR 拿 species/counts (若未显式给 counts)
    fm_poscar = os.path.join(cases_dir, "fm", "POSCAR")
    sp, counts = read_species_counts(fm_poscar)

    dft_dir = os.path.join(root, "03_dft")
    os.makedirs(dft_dir, exist_ok=True)
    with open(os.path.join(dft_dir, "INCAR.static"), "w") as f:
        f.write(build_incar_static(sp, counts, args.element, args.U, args.moment, args.encut))
    with open(os.path.join(dft_dir, "INCAR.relax"), "w") as f:
        f.write(build_incar_relax(sp, counts, args.element, args.U, args.moment, args.encut))
    with open(os.path.join(dft_dir, "KPOINTS"), "w") as f:
        f.write(KPOINTS_TMPL)
    with open(os.path.join(dft_dir, "POTCAR.spec"), "w") as f:
        f.write(build_potcar_spec(sp))

    # 组装 runs/<order>/ (INCAR 带实际 MAGMOM + KPOINTS + POSCAR)
    runs_dir = os.path.join(root, "runs")
    os.makedirs(runs_dir, exist_ok=True)
    static_tmpl = build_incar_static(sp, counts, args.element, args.U, args.moment, args.encut)
    try:
        with open(os.path.join(cases_dir, "magnetic_orders.json"), encoding="utf-8") as _f:
            _meta = json.load(_f)
        orders = _meta.get("selected_orders") or list(_meta.get("coefficients_per_cell", {}))
    except (OSError, ValueError):
        orders = ["fm", "neel", "zigzag", "stripy"]
    for order in orders:
        d = os.path.join(cases_dir, order)
        if not os.path.isdir(d):
            continue
        magmom_line = open(os.path.join(d, "magmom.txt"), encoding="utf-8").read().strip()
        incar = static_tmpl.replace(f"MAGMOM   = {counts[sp.index(args.element)]}*{args.moment:.1f} "
                                    f"{sum(counts) - counts[sp.index(args.element)]}*0.0     ! <-- 各磁态用 runs/<order>/INCAR 里的实际值",
                                    magmom_line)
        rd = os.path.join(runs_dir, order)
        os.makedirs(rd, exist_ok=True)
        with open(os.path.join(rd, "INCAR"), "w") as f:
            f.write(incar)
        shutil.copy(os.path.join(dft_dir, "KPOINTS"), os.path.join(rd, "KPOINTS"))
        shutil.copy(os.path.join(d, "POSCAR"), os.path.join(rd, "POSCAR"))

    print("写出 03_dft/ (INCAR.static/.relax, KPOINTS, POTCAR.spec) 与 runs/<order>/ 组装目录")
    print("POTCAR 顺序:", " ".join(sp))


if __name__ == "__main__":
    main()
