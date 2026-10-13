#!/usr/bin/env python
"""build_magnetic_cases_tri.py
三角晶格单离子磁体: 超胞 + 磁序枚举 (含 120° 非共线) + J1/J2/J3 能量映射系数。

模型: H = -J1 Σ_<ij> S·S - J2 Σ_<<ij>> S·S - J3 Σ_<<<ij>>> S·S  (J>0 = 铁磁)
能量映射 (每超胞): E = e0 + c1·J1~ + c2·J2~ + c3·J3~,  J~ = J·S²
  c_k = -Σ_{第 k 壳无向键} (s_i · s_j),  s_i 为单位自旋矢量 (共线 ±z 或非共线面内 120°)。

磁序池 (三角 3 子格色标 color = (n0 + n1) mod 3, 使 NN 矢量 a1/a2/a1+a2 都变色):
  fm       : 全同向 (共线)
  120deg   : 色标 0/1/2 各转 120° (非共线, 面内) —— 反解 J1 必需
  stripy_a : σ = (-1)^n0    (共线)
  stripy_b : σ = (-1)^n1    (共线)
  stripy_ab: σ = (-1)^(n0+n1) (共线)
  uud      : 两子格 ↑ 一子格 ↓ (共线 3 子格)
从池中必含 fm + 120deg, 贪心选 4 个良态构型 (SVD) 解 {e0, J1, J2, J3}。

超胞: 默认 5x5x1 (J3=2a 需 smax>=2 -> 2*2+1=5, 避免 J1/J3 周期像混叠; 见 cluster-magnets.md 规则)。

用法:
  python build_magnetic_cases_tri.py POSCAR --supercell 5 5 1 --outdir cases \
      --element Mn --nonmag I --moment 5.0 --spin 2.5 [--selftest]
"""
import argparse
import json
import os

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


def write_poscar(path, lat, atoms, coords, comment="triangular magnetic-order supercell"):
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


def min_image_frac(c1, c2):
    d = c2 - c1
    d -= np.floor(d + 0.5)
    return d


def magnetic_shells(idxs, coords, lat, max_shells=3, tol=0.2):
    ds = []
    n = len(idxs)
    for i in range(n):
        for j in range(i + 1, n):
            d = min_image_frac(coords[idxs[i]], coords[idxs[j]])
            dist = float(np.linalg.norm(d @ lat))
            if dist > 1e-6:
                ds.append(dist)
    ds = sorted(ds)
    shells = []
    for d in ds:
        if not shells or d - shells[-1] > tol:
            shells.append(d)
    return shells[:max_shells]


def bond_coeffs(idxs, coords, lat, spins, shells, tol=0.3):
    """返回 c_k = -Σ_{第 k 壳无向键} (s_i · s_j)。"""
    c = [0.0] * len(shells)
    n = len(idxs)
    for i in range(n):
        ci = coords[idxs[i]]
        for j in range(i + 1, n):
            d = min_image_frac(ci, coords[idxs[j]])
            dist = float(np.linalg.norm(d @ lat))
            k = int(np.argmin([abs(dist - s) for s in shells]))
            if abs(dist - shells[k]) < tol:
                c[k] -= float(spins[i] @ spins[j])
    return c


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poscar")
    ap.add_argument("--supercell", type=int, nargs=3, default=[5, 5, 1])
    ap.add_argument("--outdir", default="cases")
    ap.add_argument("--element", default="Mn")
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--moment", type=float, default=5.0)
    ap.add_argument("--spin", type=float, default=2.5)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    lat, atoms, coords, mode = parse_poscar(args.poscar)
    assert mode[:1] == "d", f"POSCAR 需 Direct, 得到 {mode!r}"

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
    color = (n0 + n1) % 3

    def spins_120():
        theta = color * (2.0 * np.pi / 3.0)
        s = np.zeros((n, 3))
        s[:, 0] = np.cos(theta)
        s[:, 1] = np.sin(theta)
        return s

    def spins_collinear(sigma):
        s = np.zeros((n, 3))
        s[:, 2] = sigma
        return s

    pool = {
        "fm": spins_collinear(np.ones(n)),
        "120deg": spins_120(),
        "stripy_a": spins_collinear(np.where(n0 % 2 == 0, 1.0, -1.0)),
        "stripy_b": spins_collinear(np.where(n1 % 2 == 0, 1.0, -1.0)),
        "stripy_ab": spins_collinear(np.where((n0 + n1) % 2 == 0, 1.0, -1.0)),
        "uud": spins_collinear(np.where(color == 0, -1.0, 1.0)),
    }

    coeff = {k: bond_coeffs(mag, coords_sc, lat_sc, s, shells) for k, s in pool.items()}

    required = ["fm", "120deg"]
    selected = list(required)
    candidates = [k for k in pool if k not in selected]

    def svd_score(keys):
        A = np.array([[1.0] + coeff[k] for k in keys])
        s = np.linalg.svd(A, compute_uv=False)
        if s.size == 0:
            return 0.0
        tol = s.max() * max(A.shape) * np.finfo(float).eps * 100
        nz = s[s > tol]
        rank = int(nz.size)
        return rank * 1e12 + float(nz.min()) if rank else 0.0

    while len(selected) < 4 and candidates:
        best, bs = None, -1.0
        for c in candidates:
            sc_ = svd_score(selected + [c])
            if sc_ > bs:
                bs, best = sc_, c
        if bs <= 1e-10:
            break
        selected.append(best)
        candidates.remove(best)
    if len(selected) < 4:
        raise SystemExit(f"构型不足/退化: 只能选 {len(selected)}/4 个独立构型")

    A = np.array([[1.0] + coeff[k] for k in selected])
    Ainv = np.linalg.inv(A)
    cond = np.linalg.cond(A)
    print(f"模型: E = e0 + c1*J1~ + c2*J2~ + c3*J3~  (J~=J*S^2, S={args.spin})")
    print(f"超胞 {tuple(sc)} -> {len(atoms_sc)} 原子 ({n} {args.element})")
    print(f"距离壳: J1={d_nn:.3f} A  J2={d_nnn:.3f} A  J3={d_3rd:.3f} A")
    print(f"选取 {len(selected)} 构型 {selected}, 条件数 cond(A)={cond:.3g}")
    for k in selected:
        print(f"  {k:10s} c1={coeff[k][0]:8.3f} c2={coeff[k][1]:8.3f} c3={coeff[k][2]:8.3f}")

    os.makedirs(args.outdir, exist_ok=True)
    nonmag_counts = {el: atoms_sc.count(el) for el in args.nonmag}
    for name in selected:
        d = os.path.join(args.outdir, name)
        os.makedirs(d, exist_ok=True)
        write_poscar(os.path.join(d, "POSCAR"), lat_sc, atoms_sc, coords_sc,
                     comment=f"{args.element} triangular {name}")
        mag_parts = []
        for s in pool[name]:
            mag_parts.append(f"{args.moment * s[0]:.6f} {args.moment * s[1]:.6f} {args.moment * s[2]:.6f}")
        for _ in range(sum(nonmag_counts.values())):
            mag_parts.append("0.0 0.0 0.0")
        with open(os.path.join(d, "magmom.txt"), "w") as f:
            f.write(f"MAGMOM = {' '.join(mag_parts)}\n")

    report = {
        "input": args.poscar, "supercell": sc, "n_atoms": len(atoms_sc),
        "n_mag": n, "element": args.element, "nonmag": args.nonmag,
        "moment": args.moment, "spin": args.spin, "S2": args.spin ** 2,
        "shells_angstrom": {"J1": d_nn, "J2": d_nnn, "J3": d_3rd},
        "selected_orders": selected,
        "noncollinear": "120deg" in selected,
        "coefficients_per_cell": {k: list(v) for k, v in coeff.items() if k in selected},
        "coefficients_per_mag": {k: [c / n for c in v] for k, v in coeff.items() if k in selected},
        "inversion_matrix_per_cell": Ainv.tolist(),
    }
    with open(os.path.join(args.outdir, "magnetic_orders.json"), "w") as f:
        json.dump(report, f, indent=2)

    if args.selftest:
        J_true = [-0.3, 0.05, -0.02]
        Jt = [j * args.spin ** 2 for j in J_true]
        E = {k: sum(coeff[k][m] * Jt[m] for m in range(3)) for k in selected}
        Evec = np.array([E[k] for k in selected]) / 1000.0  # meV -> eV
        sol = Ainv @ Evec
        got = [sol[1 + m] / args.spin ** 2 * 1000.0 for m in range(3)]
        ok = all(abs(g - t) < 0.02 for g, t in zip(got, J_true))
        print(f"\n--selftest: 合成 J={J_true} -> 反解 {[f'{x:.4f}' for x in got]}  "
              f"{'PASS' if ok else 'FAIL'}")


if __name__ == "__main__":
    main()
