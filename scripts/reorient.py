#!/usr/bin/env python
"""reorient.py
把数据库给出的二维单层结构整理成"层法线 = +z、面内基矢在 xy 平面、沿 z 加真空"的标准取向。

动机: Materials Project 存的 MnPS3 单层是单斜设置 (C2/m), 面内基矢 a,b 带 z 分量,
且 c 与层法线不平行、真空不足 (~3.8 A)。本脚本:
  1. 求层平面法线 n = a x b (归一化);
  2. Rodrigues 旋转使 n -> +z (层铺平到 xy);
  3. 以原子在 z 的极差为层厚, 重设 c = (层厚 + vacuum) * z_hat, 层中心放到 z=0.5;
  4. 重写分数坐标, 输出标准 POSCAR。

用法:
  python reorient.py POSCAR --element Mn --vacuum 15.0 -o POSCAR.reorient
"""

import argparse
import numpy as np


def parse_poscar(path):
    lines = [ln for ln in open(path).read().splitlines() if ln.strip()]
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


def rodrigues(n, target):
    """返回把单位向量 n 旋转到 target 的 3x3 旋转矩阵。"""
    n = n / np.linalg.norm(n)
    target = np.asarray(target, dtype=float)
    target = target / np.linalg.norm(target)
    c = float(np.dot(n, target))
    if c > 1 - 1e-12:
        return np.eye(3)
    if c < -1 + 1e-12:
        return -np.eye(3)
    k = np.cross(n, target)
    k = k / np.linalg.norm(k)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(np.arccos(c)) * K + (1 - c) * (K @ K)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poscar")
    ap.add_argument("--vacuum", type=float, default=15.0)
    ap.add_argument("-o", "--out", default=None)
    # POSCAR 第 1 行只是注释 (VASP 忽略), 但它会进产物被人读 —— 默认值不再写死材料名,
    # 由 driver 按 --model 的化学式传入 (否则 CrOCl 的 POSCAR 会自称 MnPS3)。
    ap.add_argument("--comment", default="monolayer (reoriented, vacuum)")
    args = ap.parse_args()

    lat, atoms, coords = parse_poscar(args.poscar)

    # 面内基矢: 取模长最短的两条 (层在它们张成的平面内); 对蜂窝单层就是 a,b。
    # 层法线 = a x b (归一化)。
    n = np.cross(lat[0], lat[1])
    n = n / np.linalg.norm(n)
    print(f"层法线 n = ({n[0]:+.4f}, {n[1]:+.4f}, {n[2]:+.4f})")

    R = rodrigues(n, [0, 0, 1])
    lat_r = lat @ R.T
    cart = coords @ lat          # 原始笛卡尔
    cart_r = cart @ R.T          # 旋转后笛卡尔

    # 沿层法线的周期 = 旋转后 c 的 z 分量 (即 n·c_orig)。分数坐标沿 c 会 wrap,
    # 需先按该周期把原子 z 解缠到对称区间, 才能正确量层厚 (否则"下方"原子被当"上方")。
    period = float(lat_r[2, 2])
    z = cart_r[:, 2]
    zref = float(np.median(z))
    z = ((z - zref) + period / 2.0) % period - period / 2.0
    thick = float(z.max() - z.min())
    print(f"层厚 = {thick:.3f} A (周期 {period:.3f} A), 真空 = {args.vacuum:.1f} A "
          f"-> 新 c = {thick + args.vacuum:.3f} A")

    # 新晶格: a,b 用旋转后 (应已纯面内), c = (厚+真空) z_hat
    a2 = lat_r[0].copy(); a2[2] = 0.0
    b2 = lat_r[1].copy(); b2[2] = 0.0
    c2 = np.array([0.0, 0.0, thick + args.vacuum])
    new_lat = np.array([a2, b2, c2])

    # 层中心移到 z=0.5
    zc = (z.max() + z.min()) / 2.0
    z -= zc
    z += 0.5 * (thick + args.vacuum)
    cart_r[:, 2] = z

    frac = cart_r @ np.linalg.inv(new_lat)
    frac -= np.floor(frac)

    # 检查: 面内基矢应无 z 分量, 层内原子 z 应集中在中部
    print(f"新 a = ({a2[0]:.4f}, {a2[1]:.4f}, {a2[2]:.4f})  |a|={np.linalg.norm(a2):.4f}")
    print(f"新 b = ({b2[0]:.4f}, {b2[1]:.4f}, {b2[2]:.4f})  |b|={np.linalg.norm(b2):.4f}")
    print(f"新 c = ({c2[0]:.4f}, {c2[1]:.4f}, {c2[2]:.4f})  |c|={np.linalg.norm(c2):.4f}")

    out = args.out or (args.poscar + ".reorient")
    write_poscar(out, new_lat, atoms, frac, args.comment)
    print(f"写出: {out}")


if __name__ == "__main__":
    main()
