#!/usr/bin/env python
"""cluster_sites.py
把磁性原子按成键距离聚成"簇"（如 Nb3I8 的 Nb3 三聚体），输出簇质心作为有效磁位。

簇磁体（cluster Mott insulator）的模型：每个磁性簇 = 一个有效自旋（Nb3 → S=1/2），
簇中心构成有效晶格（Nb3I8 是三角晶格）。本脚本做簇识别 + 质心，供下游 MC / DFT 用。

用法:
  python cluster_sites.py POSCAR --element Nb --cutoff 3.3 -o reduced.vasp
输出:
  - reduced.vasp : 磁性原子替换为簇质心（元素仍用原名），非磁原子保留
  - 打印簇数 / 每簇原子数 / 簇质心
"""

import argparse
import os
import sys

import numpy as np


def parse_poscar(path):
    lines = [ln for ln in open(path, encoding="utf-8").read().splitlines() if ln.strip()]
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
    return lat, atoms, coords


def min_image_frac(c1, c2):
    d = c2 - c1
    d -= np.floor(d + 0.5)
    return d


def find_clusters(atoms, coords, lat, element, cutoff):
    """把 element 原子按 min-image 距离 < cutoff 聚成簇。返回簇列表（每簇 = 原子索引列表）。"""
    idxs = [i for i, a in enumerate(atoms) if a == element]
    n = len(idxs)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    for i in range(n):
        for j in range(i + 1, n):
            d = min_image_frac(coords[idxs[i]], coords[idxs[j]]) @ lat
            if np.linalg.norm(d) < cutoff:
                union(i, j)

    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(idxs[i])
    return list(groups.values())


def cluster_centroids(coords, lat, groups):
    """每簇质心（笛卡尔平均，正确处理 PBC 解缠），返回 frac 坐标列表。"""
    cents = []
    for g in groups:
        ref = coords[g[0]]
        cart = []
        for gi in g:
            d = min_image_frac(ref, coords[gi])
            cart.append((ref + d) @ lat)
        c = np.mean(cart, axis=0)
        frac = c @ np.linalg.inv(lat)
        frac -= np.floor(frac)
        cents.append(frac)
    return np.array(cents)


def write_poscar(path, lat, atoms, coords, comment):
    order = []
    for a in atoms:
        if a not in order:
            order.append(a)
    counts = [atoms.count(sp) for sp in order]
    out = [comment, "1.00000000000000"]
    for row in lat:
        out.append(" ".join(f"{x:20.12f}" for x in row))
    out.append("  ".join(order))
    out.append("  ".join(str(c) for c in counts))
    out.append("Direct")
    for sp in order:
        for a, c in zip(atoms, coords):
            if a == sp:
                out.append(" ".join(f"{x:20.12f}" for x in c))
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


def main():
    ap = argparse.ArgumentParser(description="磁性簇识别 -> 簇质心有效磁位")
    ap.add_argument("poscar")
    ap.add_argument("--element", required=True, help="磁性元素, 如 Nb")
    ap.add_argument("--cutoff", type=float, default=3.3, help="簇内成键阈值 (A)")
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    lat, atoms, coords = parse_poscar(args.poscar)
    groups = find_clusters(atoms, coords, lat, args.element, args.cutoff)
    cents = cluster_centroids(coords, lat, groups)

    print(f"{args.element} 原子 {sum(1 for a in atoms if a == args.element)} 个 -> "
          f"{len(groups)} 个簇")
    for i, (g, c) in enumerate(zip(groups, cents)):
        sizes = [f"{np.linalg.norm(min_image_frac(coords[g[0]], coords[gi]) @ lat):.2f}A"
                 for gi in g[1:]]
        print(f"  簇{i}: {len(g)} 个 {args.element} (质心 frac {c[0]:.4f}, {c[1]:.4f}, {c[2]:.4f})"
              f" 内距 {sizes}")

    # 写 reduced POSCAR: 磁性原子替换为簇质心, 非磁原子保留
    keep_atoms = [a for a in atoms if a != args.element]
    keep_coords = [c for a, c in zip(atoms, coords) if a != args.element]
    new_atoms = keep_atoms + [args.element] * len(cents)
    new_coords = np.vstack([np.array(keep_coords), cents]) if keep_coords else cents
    if not len(keep_atoms):
        new_coords = cents

    out = args.out or (args.poscar + ".clusters.vasp")
    write_poscar(out, lat, new_atoms, new_coords,
                 f"{args.element} clusters -> {len(cents)} effective sites")
    print(f"写出 reduced POSCAR: {out}")


if __name__ == "__main__":
    main()
