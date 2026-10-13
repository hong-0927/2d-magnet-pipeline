#!/usr/bin/env python
"""assemble_dft_tri.py
从三角磁序 cases 生成 VASP 非共线输入: INCAR / KPOINTS / POTCAR.spec, 并组装 runs/<order>/。

非共线 (含 120°): LNONCOLLINEAR = .TRUE., MAGMOM 为每原子 3 分量 (自旋沿 x/y/z)。
共线态 (fm/stripy/uud) 也用同一非共线设置跑 (自旋沿 z), 保证全组可比。

用法:
  python assemble_dft_tri.py --root <proj> --element Mn --nonmag I --U 4.0 --moment 5.0
"""
import argparse
import os
import shutil


def read_species_counts(poscar_path):
    lines = [ln for ln in open(poscar_path, encoding="utf-8").read().splitlines() if ln.strip()]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    return species, counts


def build_incar(species, counts, element, U, moment, encut, npar=4):
    total_mag = counts[species.index(element)] if element in species else 0
    n_nonmag = sum(counts) - total_mag
    ldaul = " ".join("2" if s == element else "-1" for s in species)
    ldauu = " ".join(f"{U:.1f}" if s == element else "0.0" for s in species)
    ldauj = " ".join("0.0" for _ in species)
    return f"""SYSTEM   = {" ".join(species)} monolayer triangular magnetic order
PREC     = Normal
ENCUT    = {encut}
EDIFF    = 1E-6
ALGO     = Normal
NELM     = 200
ISMEAR   = 0
SIGMA    = 0.05

LNONCOLLINEAR = .TRUE.
MAGMOM   = {total_mag}*{moment:.1f} 0.0 0.0 {n_nonmag}*0.0 0.0 0.0     ! <-- 实际用 runs/<order>/INCAR

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


KPOINTS_TMPL = """Gamma-centered k-mesh for triangular monolayer supercell (vacuum along c)
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
    ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--moment", type=float, default=5.0)
    ap.add_argument("--encut", type=float, default=520)
    ap.add_argument("--cases", default="02_supercell/cases")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    cases_dir = os.path.join(root, args.cases)
    orders = sorted(d for d in os.listdir(cases_dir)
                    if os.path.isdir(os.path.join(cases_dir, d))
                    and os.path.exists(os.path.join(cases_dir, d, "POSCAR")))
    fm_poscar = os.path.join(cases_dir, orders[0], "POSCAR")
    sp, counts = read_species_counts(fm_poscar)

    dft_dir = os.path.join(root, "03_dft")
    os.makedirs(dft_dir, exist_ok=True)
    static_tmpl = build_incar(sp, counts, args.element, args.U, args.moment, args.encut)
    with open(os.path.join(dft_dir, "INCAR.static"), "w") as f:
        f.write(static_tmpl)
    with open(os.path.join(dft_dir, "KPOINTS"), "w") as f:
        f.write(KPOINTS_TMPL)
    with open(os.path.join(dft_dir, "POTCAR.spec"), "w") as f:
        f.write(build_potcar_spec(sp))

    nm = counts[sp.index(args.element)]
    n_nonmag = sum(counts) - nm
    default_line = f"MAGMOM   = {nm}*{args.moment:.1f} 0.0 0.0 {n_nonmag}*0.0 0.0 0.0     ! <-- 实际用 runs/<order>/INCAR"

    runs_dir = os.path.join(root, "runs")
    os.makedirs(runs_dir, exist_ok=True)
    for order in orders:
        d = os.path.join(cases_dir, order)
        magmom_line = open(os.path.join(d, "magmom.txt"), encoding="utf-8").read().strip()
        incar = static_tmpl.replace(default_line, magmom_line)
        rd = os.path.join(runs_dir, order)
        os.makedirs(rd, exist_ok=True)
        with open(os.path.join(rd, "INCAR"), "w") as f:
            f.write(incar)
        shutil.copy(os.path.join(dft_dir, "KPOINTS"), os.path.join(rd, "KPOINTS"))
        shutil.copy(os.path.join(d, "POSCAR"), os.path.join(rd, "POSCAR"))

    print("写出 03_dft/ (INCAR.static, KPOINTS, POTCAR.spec) 与 runs/<order>/")
    print("非共线 (LNONCOLLINEAR=.TRUE.), 磁序:", orders)
    print("POTCAR 顺序:", " ".join(sp))


if __name__ == "__main__":
    main()
