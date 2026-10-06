#!/usr/bin/env python
"""build_magnetic_cases3.py
蜂窝磁序超胞生成 + J1/J2/J3 能量映射系数 (相对 CrI3 流水线的 J1/J2 版, 扩展到第三近邻 J3)。

MnPS3 是蜂窝 Néel 反铁磁 (S=5/2), 其 J3 不可忽略 (Sivadas 2015: J3>J2)。所以这里把
"cr_shells" 扩展到识别 3 个距离壳, "directed_bond_sums" 扩展到统计 3 类键和, 并
把 zigzag/stripy 两个 M 点磁序区分开 (它们在 J1-J2 模型简并, 但 J3 会破简并, 见
Sivadas PRB 91,235425 式(2): E_zigzag/stripy = E0 + (±J1 -2J2 ∓3J3)|S|^2)。

用法:
  python build_magnetic_cases3.py POSCAR --supercell 2 2 1 --outdir cases \
      --element Mn --nonmag "P S" --moment 5.0 --spin 2.5
"""

import argparse
import json
import os
from collections import deque
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


def write_poscar(path, lat, atoms, coords, comment="magnetic-order supercell"):
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
    """用显式平移搜索 (不做最小镜像截断) 识别磁性亚晶格前 max_shells 个距离壳。"""
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
    """统计各类键的 directed sum sigma_i*sigma_j, 返回与 shells 等长的列表。"""
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
                # 归类到最近壳
                k = int(np.argmin([abs(d - s) for s in shells]))
                if abs(d - shells[k]) < 0.3:
                    sums[k] += sigma[i] * sigma[j]
    return sums


def min_image_distance(c1, c2, lat):
    d = c2 - c1
    d -= np.floor(d + 0.5)
    return float(np.linalg.norm(d @ lat))


def sublattice_colors(idxs, coords, lat, d_nn):
    n = len(idxs)
    color = np.full(n, -1, dtype=int)
    color[0] = 0
    q = deque([0])
    while q:
        i = q.popleft()
        for j in range(n):
            if i == j:
                continue
            if min_image_distance(coords[idxs[i]], coords[idxs[j]], lat) < d_nn + 0.5:
                if color[j] == -1:
                    color[j] = 1 - color[i]
                    q.append(j)
    return color


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poscar")
    ap.add_argument("--supercell", type=int, nargs=3, default=[2, 2, 1])
    ap.add_argument("--outdir", default="cases")
    ap.add_argument("--element", default="Mn")
    ap.add_argument("--nonmag", nargs="+", default=["P", "S"])
    ap.add_argument("--moment", type=float, default=5.0)
    ap.add_argument("--spin", type=float, default=2.5)
    args = ap.parse_args()

    lat, atoms, coords, mode = parse_poscar(args.poscar)
    assert mode[:1] == "d", f"POSCAR 需 Direct, 得到 {mode!r}"

    sc = args.supercell
    lat_sc, atoms_sc, coords_sc, cell_idx = make_supercell(lat, atoms, coords, sc)
    cr = [i for i, a in enumerate(atoms_sc) if a == args.element]
    n = len(cr)
    shells = magnetic_shells(cr, coords_sc, lat_sc, max_shells=3)
    if len(shells) < 3:
        raise SystemExit(f"只找到 {len(shells)} 个磁距壳, 需要 3 (J1/J2/J3); "
                         f"请加大 --supercell")
    d_nn, d_nnn, d_3rd = shells
    color = sublattice_colors(cr, coords_sc, lat_sc, d_nn)

    # 磁序 sigma
    def m_pattern(which):
        n0 = cell_idx[cr][:, 0]
        n1 = cell_idx[cr][:, 1]
        phase = [n0, n1, n0 + n1][which]
        return np.where(phase % 2 == 0, 1.0, -1.0)

    sigmas = {
        "fm": np.ones(n),
        "neel": np.where(color == 0, 1.0, -1.0),
        "M0": m_pattern(0),
        "M1": m_pattern(1),
        "M01": m_pattern(2),
    }

    sums = {}
    for k, s in sigmas.items():
        sums[k] = classify_bonds(cr, coords_sc, lat_sc, s, shells)

    # zigzag: M 点里 NN 键和最大 (+8, 即 J1 系数最小); stripy: 另一个不同的 M 点成员
    m_keys = ["M0", "M1", "M01"]
    zig_key = min(m_keys, key=lambda k: sums[k][0])
    rest = [k for k in m_keys if k != zig_key]
    # stripy 取与 zigzag NN 键和相同但第三近邻键和不同的那个
    strip_key = None
    for k in rest:
        if abs(sums[k][0] - sums[zig_key][0]) < 1e-9 and abs(sums[k][2] - sums[zig_key][2]) > 1e-9:
            strip_key = k
            break
    if strip_key is None:
        strip_key = rest[0]
    zig_sig = sigmas[zig_key]
    strip_sig = sigmas[strip_key]

    orders = {"fm": sigmas["fm"], "neel": sigmas["neel"],
              "zigzag": zig_sig, "stripy": strip_sig}
    sums_std = {
        "fm": sums["fm"], "neel": sums["neel"],
        "zigzag": sums[zig_key], "stripy": sums[strip_key],
    }
    # E_mag = e0 + c1*J1~ + c2*J2~ + c3*J3~,  c_k = -0.5 * bond_sum_k (directed 双计)
    coeff = {k: [-0.5 * s[0], -0.5 * s[1], -0.5 * s[2]] for k, s in sums_std.items()}

    # 4x4 反解 [e0, J1~, J2~, J3~] = Ainv @ [E_fm, E_neel, E_zigzag, E_stripy]
    A = np.array([[1.0] + coeff["fm"],
                  [1.0] + coeff["neel"],
                  [1.0] + coeff["zigzag"],
                  [1.0] + coeff["stripy"]])
    Ainv = np.linalg.inv(A)

    os.makedirs(args.outdir, exist_ok=True)
    nonmag_counts = {el: atoms_sc.count(el) for el in args.nonmag}
    for name, sig in orders.items():
        d = os.path.join(args.outdir, name)
        os.makedirs(d, exist_ok=True)
        write_poscar(os.path.join(d, "POSCAR"), lat_sc, atoms_sc, coords_sc,
                     comment=f"MnPS3 monolayer {name}")
        mag = " ".join(f"{args.moment * s:.1f}" for s in sig)
        for el, cnt in nonmag_counts.items():
            mag += " " + " ".join(["0.0"] * cnt)
        with open(os.path.join(d, "magmom.txt"), "w") as f:
            f.write(f"MAGMOM = {mag}\n")

    report = {
        "input": args.poscar, "supercell": sc, "n_atoms": len(atoms_sc),
        "n_mag": n, "element": args.element, "nonmag": args.nonmag,
        "moment": args.moment, "spin": args.spin, "S2": args.spin ** 2,
        "shells_angstrom": {"J1": d_nn, "J2": d_nnn, "J3": d_3rd},
        "zigzag_key": zig_key, "stripy_key": strip_key,
        "coefficients_per_cell": {k: list(v) for k, v in coeff.items()},
        "coefficients_per_mag": {k: [c / n for c in v] for k, v in coeff.items()},
        "inversion_matrix_per_cell": Ainv.tolist(),
    }
    with open(os.path.join(args.outdir, "magnetic_orders.json"), "w") as f:
        json.dump(report, f, indent=2)

    print(f"supercell {sc} -> {len(atoms_sc)} atoms ({n} {args.element}), "
          f"nonmag {nonmag_counts}")
    print(f"距离壳: J1={d_nn:.3f} A  J2={d_nnn:.3f} A  J3={d_3rd:.3f} A")
    print(f"zigzag 取 {zig_key}, stripy 取 {strip_key}")
    print("\nE_mag [per cell] = c1*J1~ + c2*J2~ + c3*J3~ , J~ = J*S^2, S=", args.spin)
    for name in ["fm", "neel", "zigzag", "stripy"]:
        c1, c2, c3 = coeff[name]
        print(f"  {name:8s}: c1={c1:8.3f}  c2={c2:8.3f}  c3={c3:8.3f}   "
              f"(per Mn: {c1/n:+.3f}, {c2/n:+.3f}, {c3/n:+.3f})")
    print("\n期望 (Sivadas 式(2), 每 Mn): fm(+3,+6,+3) neel(-3,+6,-3) "
          "zigzag(+1,-2,-3) stripy(-1,-2,+3)  [x|S|^2 系数]")
    print("  注: 本表 c 是 -0.5*bond_sum, 对每 Mn 除以 n 后应与上式一致。")


if __name__ == "__main__":
    main()
