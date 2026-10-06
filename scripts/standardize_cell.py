#!/usr/bin/env python
"""standardize_cell.py
把二维六方单层的原胞标准化成"紧凑蜂窝原胞"——使胞内两个磁性原子互为最近邻。

为什么需要: 数据库(如 Materials Project)给出的原胞常常是"非标准设置"
(胞内两磁性原子相距 2a/sqrt3 而不是 a/sqrt3)。在这种原胞里, 教科书里的
zigzag/stripy 磁序无法写成纯晶胞索引模式 (+1)^n0, 会算错能量映射系数,
进而算错 J2。标准化成紧凑原胞后, 标准定义与 PPT/文献一致。

做法: 在输入晶格矢 a1,a2 的整数组合中搜索 (det = +/-1), 取使胞内
磁性原子间距最小的那一组作为新原胞。晶格本身不变, 只是换了一组基矢。

用法:
  python standardize_cell.py POSCAR --element Cr -o POSCAR.std
"""

import argparse
import numpy as np


def _is_float(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def parse_poscar(path):
    lines = [ln.rstrip("\n") for ln in open(path) if ln.strip()]
    scale = float(lines[1].split()[0])
    L = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]) * scale
    idx = 5
    tok = lines[idx].split()
    if all(_is_float(t) for t in tok):
        symbols, counts, idx = None, [int(t) for t in tok], idx + 1
    else:
        symbols, counts, idx = tok, [int(t) for t in lines[idx + 1].split()], idx + 2
    if lines[idx].strip()[:1] in ("s", "S"):
        idx += 1
    mode = lines[idx].strip()[:1]
    idx += 1
    nat = sum(counts)
    coords = np.array([[float(x) for x in lines[idx + i].split()[:3]] for i in range(nat)])
    if mode in ("c", "C", "k", "K"):
        coords = coords @ np.linalg.inv(L)
    coords -= np.floor(coords)
    if symbols is None:
        symbols = ["X"] * nat
    else:
        symbols = [s for s, c in zip(symbols, counts) for _ in range(c)]
    return L, symbols, coords


def write_poscar(path, L, symbols, coords, comment):
    order = []
    for s in symbols:
        if s not in order:
            order.append(s)
    counts = [symbols.count(s) for s in order]
    out = [comment, "1.00000000000000"]
    for v in L:
        out.append(" ".join(f"{x:20.12f}" for x in v))
    out.append("  ".join(order))
    out.append("  ".join(str(c) for c in counts))
    out.append("Direct")
    for s in order:
        for a, c in zip(symbols, coords):
            if a == s:
                out.append(" ".join(f"{x:20.12f}" for x in c))
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


def main():
    ap = argparse.ArgumentParser(description="标准化二维蜂窝原胞")
    ap.add_argument("poscar")
    ap.add_argument("--element", required=True, help="磁性元素, 如 Cr")
    ap.add_argument("--search", type=int, default=2, help="整数组合搜索范围")
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    L, symbols, coords = parse_poscar(args.poscar)
    mag = [i for i, a in enumerate(symbols) if a == args.element]
    if len(mag) != 2:
        raise SystemExit(f"期望原胞里有 2 个 {args.element}, 实际 {len(mag)} 个")

    r = args.search
    best = None
    for m1 in range(-r, r + 1):
        for n1 in range(-r, r + 1):
            for m2 in range(-r, r + 1):
                for n2 in range(-r, r + 1):
                    if abs(m1 * n2 - n1 * m2) != 1:
                        continue
                    na1 = m1 * L[0] + n1 * L[1]
                    na2 = m2 * L[0] + n2 * L[1]
                    newL = np.array([na1, na2, L[2]])
                    if np.linalg.det(newL) < 0:
                        newL = newL[[1, 0, 2]]
                    nf = coords @ L @ np.linalg.inv(newL)
                    nf -= np.floor(nf)
                    # 胞内两磁性原子距离 (用实际分数坐标, 非最小镜像)
                    d = float(np.linalg.norm((nf[mag[1]] - nf[mag[0]]) @ newL))
                    # 字典序目标: 先要最短基矢(|a1|+|a2| 最小), 再要胞内最近邻
                    span = float(np.linalg.norm(newL[0]) + np.linalg.norm(newL[1]))
                    key = (round(span, 1), round(d, 4))
                    if best is None or key < best[0]:
                        best = (key, d, newL, nf, (m1, n1, m2, n2))

    _, d, newL, nf, combo = best
    d_orig = float(np.linalg.norm((coords[mag[1]] - coords[mag[0]]) @ L))
    print(f"原胞内 {args.element}-{args.element}: {d_orig:.3f} A -> 标准化后 {d:.3f} A")
    print(f"新基矢组合 (a1,a2) = ({combo[0]},{combo[1]}) , ({combo[2]},{combo[3]})")
    print(f"|a1|={np.linalg.norm(newL[0]):.4f}  |a2|={np.linalg.norm(newL[1]):.4f}  "
          f"|c|={np.linalg.norm(newL[2]):.4f}")

    out = args.out or (args.poscar + ".std")
    write_poscar(out, newL, symbols, nf, f"{args.element} standard compact cell")
    print(f"写出: {out}")


if __name__ == "__main__":
    main()
