#!/usr/bin/env python
"""build_rectangular_cases.py
矩形晶格单离子磁体: 超胞 + 共线磁序枚举 + J1/J2/J3 能量映射系数。

矩形晶格 (正交, a != b) 与三角/蜂窝的根本区别:
  - 三角/蜂窝是 60° (a=b) 晶格, 近邻配位 6/3, 有 120° 非共线 (受挫) 磁序;
  - 矩形晶格 a 垂直 b 且 |a| != |b|, 磁原子配位 = 2 (a 轴) + 2 (b 轴) + 4 (对角),
    是**非受挫**晶格, 低能态全是共线, 无需 120° 非共线态。
  典型材料: CrSBr (a=3.54 A, b=4.76 A; 三项各向异性铁磁 J, 见 Ke Yang & Hua Wu, PRB 104, 144416 (2021))。

模型: H = -J1 Σ_<ij> S·S - J2 Σ_<<ij>> S·S - J3 Σ_<<<ij>>> S·S  (J>0 = 铁磁)
  距离壳按升序: J1 = a 轴近邻, J2 = b 轴近邻, J3 = 对角 (|a+b|)。
能量映射 (每超胞): E = e0 + c1·J1~ + c2·J2~ + c3·J3~,  J~ = J·S²
  c_k = -0.5 * Σ_{第 k 壳有向键} σ_i·σ_j  (σ = ±1 共线)。

四个共线磁序 (恰好线性无关, 解 {e0, J1, J2, J3}):
  fm      : 全同向
  afm_a   : σ = (-1)^n0   (沿 a 反平行条纹)
  afm_b   : σ = (-1)^n1   (沿 b 反平行条纹)
  afm_diag: σ = (-1)^(n0+n1) (棋盘 / 对角 AFM)

超胞: 默认 4x4x1 (**偶数**; 条纹态 afm_a/afm_b/afm_diag 才能跨周期边界干净闭合,
     奇数超胞会在边界多出一条「同号键」使条纹态带缺陷, 系数偏离理想值)。

用法:
  python build_rectangular_cases.py POSCAR --supercell 5 5 1 --outdir cases \
      --element Cr --nonmag S Br --moment 3.0 --spin 1.5 [--selftest]
"""
import argparse
import json
import os
from itertools import product

import numpy as np


def parse_poscar(path):
    with open(path) as f:
        lines = [ln for ln in f.read().splitlines() if ln.strip()]
    scale = float(lines[1].strip())
    lat = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)]) * scale
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    idx = 7
    if lines[idx].strip().lower()[:1] == "s":
        idx += 1
    mode = lines[idx].strip().lower()
    idx += 1
    coords, atoms = [], []
    for sp, n in zip(species, counts):
        for _ in range(n):
            coords.append([float(x) for x in lines[idx].split()[:3]])
            atoms.append(sp)
            idx += 1
    coords = np.array(coords)
    if mode[:1] in "ck":
        coords = coords @ np.linalg.inv(lat)
    coords -= np.floor(coords)
    return lat, atoms, coords, mode


def write_poscar(path, lat, atoms, coords, comment="rectangular magnetic-order supercell"):
    order = []
    for a in atoms:
        if a not in order:
            order.append(a)
    counts = [atoms.count(sp) for sp in order]
    lines = [comment, "1.00000000000000"]
    for row in lat:
        lines.append(" ".join(f"{x:20.12f}" for x in row))
    lines.append("  ".join(order))
    lines.append("  ".join(str(c) for c in counts))
    lines.append("Direct")
    for sp in order:
        for a, c in zip(atoms, coords):
            if a == sp:
                lines.append(" ".join(f"{x:20.12f}" for x in c))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def make_supercell(lat, atoms, coords, sc):
    sc = np.asarray(sc, dtype=int)
    new_atoms, new_coords, cell_idx = [], [], []
    for n0 in range(sc[0]):
        for n1 in range(sc[1]):
            for n2 in range(sc[2]):
                for a, c in zip(atoms, coords):
                    fc = (c + np.array([n0, n1, n2], dtype=float)) / sc
                    new_atoms.append(a)
                    new_coords.append(fc)
                    cell_idx.append((n0, n1, n2))
    new_lat = lat.copy()
    new_lat[0] *= sc[0]
    new_lat[1] *= sc[1]
    new_lat[2] *= sc[2]
    return new_lat, new_atoms, np.array(new_coords), np.array(cell_idx, dtype=int)


def magnetic_shells(idxs, coords, lat, max_shells=3, nshift=2, tol=0.2):
    """显式平移搜索识别磁性亚晶格前 max_shells 个距离壳 (不做最小镜像截断)。"""
    shifts = list(product(range(-nshift, nshift + 1), repeat=3))
    ds = []
    n = len(idxs)
    for a in range(n):
        ci = coords[idxs[a]]
        for b in range(n):
            base = coords[idxs[b]] - ci
            for sh in shifts:
                cart = (base + np.array(sh, dtype=float)) @ lat
                d = float(np.linalg.norm(cart))
                if d < 1e-6:
                    continue
                ds.append(d)
    ds = sorted(ds)
    shells = []
    for d in ds:
        if not shells or d - shells[-1] > tol:
            shells.append(d)
    return shells[:max_shells]


def classify_bonds(idxs, coords, lat, sigma, shells, nshift=1):
    """有向键和 Σ σ_i σ_j 按壳归类, 返回与 shells 等长的列表。"""
    shifts = list(product(range(-nshift, nshift + 1), repeat=3))
    sums = [0.0] * len(shells)
    n = len(idxs)
    for i in range(n):
        ci = coords[idxs[i]]
        for j in range(n):
            base = coords[idxs[j]] - ci
            for sh in shifts:
                cart = (base + np.array(sh, dtype=float)) @ lat
                d = float(np.linalg.norm(cart))
                if d < 1e-6:
                    continue
                k = int(np.argmin([abs(d - s) for s in shells]))
                if abs(d - shells[k]) < 0.3:
                    sums[k] += sigma[i] * sigma[j]
    return sums


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poscar")
    ap.add_argument("--supercell", type=int, nargs=3, default=[4, 4, 1])
    ap.add_argument("--outdir", default="cases")
    ap.add_argument("--element", default="Cr")
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--moment", type=float, default=3.0)
    ap.add_argument("--spin", type=float, default=1.5)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    lat, atoms, coords, mode = parse_poscar(args.poscar)  # parse_poscar 已把 Cartesian 归一为 Direct

    sc = args.supercell
    lat_sc, atoms_sc, coords_sc, cell_idx = make_supercell(lat, atoms, coords, sc)
    mag = [i for i, a in enumerate(atoms_sc) if a == args.element]
    n = len(mag)
    shells = magnetic_shells(mag, coords_sc, lat_sc, max_shells=3)
    if len(shells) < 3:
        raise SystemExit(f"只找到 {len(shells)} 个磁距壳, 需要 3 (J1/J2/J3); 请加大 --supercell")
    d_nn, d_nnn, d_3rd = shells

    n0 = cell_idx[mag][:, 0]
    n1 = cell_idx[mag][:, 1]

    orders = ["fm", "afm_a", "afm_b", "afm_diag"]
    sigmas = {
        "fm": np.ones(n),
        "afm_a": np.where(n0 % 2 == 0, 1.0, -1.0),
        "afm_b": np.where(n1 % 2 == 0, 1.0, -1.0),
        "afm_diag": np.where((n0 + n1) % 2 == 0, 1.0, -1.0),
    }

    sums = {k: classify_bonds(mag, coords_sc, lat_sc, s, shells) for k, s in sigmas.items()}
    # E_mag = e0 + c1*J1~ + c2*J2~ + c3*J3~,  c_k = -0.5 * 有向键和 (有向双计)
    coeff = {k: [-0.5 * s[0], -0.5 * s[1], -0.5 * s[2]] for k, s in sums.items()}

    A = np.array([[1.0] + coeff[k] for k in orders])
    Ainv = np.linalg.inv(A)
    cond = np.linalg.cond(A)

    os.makedirs(args.outdir, exist_ok=True)
    nonmag_counts = {el: atoms_sc.count(el) for el in args.nonmag}
    for name in orders:
        d = os.path.join(args.outdir, name)
        os.makedirs(d, exist_ok=True)
        write_poscar(os.path.join(d, "POSCAR"), lat_sc, atoms_sc, coords_sc,
                     comment=f"{args.element} rectangular {name}")
        mag_parts = [f"{args.moment * s:.1f}" for s in sigmas[name]]
        for el, cnt in nonmag_counts.items():
            mag_parts += ["0.0"] * cnt
        with open(os.path.join(d, "magmom.txt"), "w") as f:
            f.write(f"MAGMOM = {' '.join(mag_parts)}\n")

    report = {
        "input": args.poscar, "supercell": sc, "n_atoms": len(atoms_sc),
        "n_mag": n, "element": args.element, "nonmag": args.nonmag,
        "moment": args.moment, "spin": args.spin, "S2": args.spin ** 2,
        "shells_angstrom": {"J1": d_nn, "J2": d_nnn, "J3": d_3rd},
        "selected_orders": orders,
        "noncollinear": False,
        "coefficients_per_cell": {k: list(v) for k, v in coeff.items()},
        "coefficients_per_mag": {k: [c / n for c in v] for k, v in coeff.items()},
        "inversion_matrix_per_cell": Ainv.tolist(),
    }
    with open(os.path.join(args.outdir, "magnetic_orders.json"), "w") as f:
        json.dump(report, f, indent=2)

    print(f"模型: E = e0 + c1*J1~ + c2*J2~ + c3*J3~  (J~=J*S^2, S={args.spin})")
    print(f"超胞 {tuple(sc)} -> {len(atoms_sc)} 原子 ({n} {args.element}), "
          f"nonmag {nonmag_counts}")
    print(f"距离壳: J1={d_nn:.3f} A  J2={d_nnn:.3f} A  J3={d_3rd:.3f} A")
    print(f"四共线磁序 {orders}, 条件数 cond(A)={cond:.3g}")
    for k in orders:
        print(f"  {k:10s} c1={coeff[k][0]:8.3f} c2={coeff[k][1]:8.3f} c3={coeff[k][2]:8.3f}"
              f"   (per mag: {coeff[k][0]/n:+.3f}, {coeff[k][1]/n:+.3f}, {coeff[k][2]/n:+.3f})")
    print("\n期望 (矩形晶格, 每 mag): fm(-1,-1,-2) afm_a(+1,-1,+2) "
          "afm_b(-1,+1,+2) afm_diag(+1,+1,-2)  [x|S|^2 系数]")

    if args.selftest:
        # 1) 系数对照理想值 (验证超胞为偶数 + 够大, 条纹态纯净, 无周期边界缺陷)
        ideal = {"fm": (-1, -1, -2), "afm_a": (+1, -1, +2),
                 "afm_b": (-1, +1, +2), "afm_diag": (+1, +1, -2)}
        max_err = max(abs(coeff[k][m] / n - ideal[k][m]) for k in orders for m in range(3))
        coef_ok = max_err < 0.15
        # 2) 合成能量 -> 反解 J (验证 A 可逆 + 系数自洽)
        J_true = [2.79, 2.12, 3.88]  # CrSBr 按距离排 (a 轴/b 轴/对角), 全铁磁
        Jt = [j * args.spin ** 2 for j in J_true]
        E = {k: sum(coeff[k][m] * Jt[m] for m in range(3)) for k in orders}
        Evec = np.array([E[k] for k in orders]) / 1000.0  # meV -> eV
        sol = Ainv @ Evec
        got = [sol[1 + m] / args.spin ** 2 * 1000.0 for m in range(3)]
        inv_ok = all(abs(g - t) < 0.02 for g, t in zip(got, J_true))
        print(f"\n--selftest: 系数最大偏差={max_err:.3f} {'PASS' if coef_ok else 'FAIL'}")
        print(f"--selftest: 合成 J={J_true} -> 反解 {[f'{x:.4f}' for x in got]}  "
              f"{'PASS' if inv_ok else 'FAIL'}")


if __name__ == "__main__":
    main()
