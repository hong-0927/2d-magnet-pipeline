#!/usr/bin/env python
"""render_supercell.py -- 画 2x2x1 超胞四磁序的俯视图 (matplotlib, 无需 OVITO)。

读 02_supercell/cases/<order>/{POSCAR,magmom.txt}, 投影到 xy 俯视:
  Cr 按自旋上/下染色 (红/蓝), 非磁元素浅灰, 画 Cr-Cr 最近邻键 + 超胞边界。

用法:
  python render_supercell.py --root <CrBr3_pipeline>
"""
import argparse
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

C_UP = "#d62728"
C_DN = "#1f77b4"
C_NON = "#c8c8c8"
C_BOND = "#9a9a9a"


def _isfloat(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def parse_poscar(path):
    lines = [l.rstrip("\n") for l in open(path, encoding="utf-8") if l.strip()]
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]) * scale
    idx = 5
    tok = lines[idx].split()
    if all(_isfloat(t) for t in tok):
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
        coords = coords @ np.linalg.inv(lat)
    coords -= np.floor(coords)
    if symbols is None:
        symbols = ["X"] * nat
    else:
        symbols = [s for s, c in zip(symbols, counts) for _ in range(c)]
    return lat, symbols, coords


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--orders", nargs="+", default=["fm", "neel", "zigzag", "stripy"])
    ap.add_argument("--element", default="Mn")
    ap.add_argument("--nonmag", nargs="+", default=["P", "S"])
    ap.add_argument("--label", default="MnPS3")
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    titles = {"fm": "FM (铁磁)", "neel": "Néel (反铁磁)", "zigzag": "zigzag", "stripy": "stripy"}

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 9.4), dpi=200)
    axes = axes.ravel()

    for ax, order in zip(axes, args.orders):
        d = os.path.join(root, "02_supercell", "cases", order)
        lat, symbols, coords = parse_poscar(os.path.join(d, "POSCAR"))
        lat2 = lat[:2, :2]
        xy = coords[:, :2] @ lat2

        mom = np.array([float(t) for t in
                        open(os.path.join(d, "magmom.txt"), encoding="utf-8").read().split("=", 1)[1].split()])
        mag_idx = [i for i, s in enumerate(symbols) if s == args.element]
        non_idx = [i for i, s in enumerate(symbols) if s in args.nonmag]
        cr_mom = mom[:len(mag_idx)]

        crf = coords[mag_idx, :2]
        ds = []
        for i in range(len(mag_idx)):
            for j in range(i + 1, len(mag_idx)):
                dd = crf[j] - crf[i]
                dd -= np.round(dd)
                ds.append(np.linalg.norm(dd @ lat2))
        d_nn = min(ds)

        a1, a2 = lat2[0], lat2[1]
        crxy = xy[mag_idx]
        up = [k for k, m in enumerate(cr_mom) if m > 0]
        dn = [k for k, m in enumerate(cr_mom) if m < 0]

        # 背景周期像 (淡, 只画 Cr): 让蜂窝在超胞边界处看起来完整
        for oi in (-1, 0, 1):
            for oj in (-1, 0, 1):
                if oi == 0 and oj == 0:
                    continue
                sh = oi * a1 + oj * a2
                ax.scatter(crxy[up, 0] + sh[0], crxy[up, 1] + sh[1],
                           s=150, c=C_UP, alpha=0.16, edgecolors="none", zorder=1)
                ax.scatter(crxy[dn, 0] + sh[0], crxy[dn, 1] + sh[1],
                           s=150, c=C_DN, alpha=0.16, edgecolors="none", zorder=1)

        # Cr-Cr 最近邻键 (中心胞 -> 周期像, 不裁剪)
        for i in range(len(mag_idx)):
            for j in range(i + 1, len(mag_idx)):
                dd = crf[j] - crf[i]
                dd -= np.round(dd)
                if abs(np.linalg.norm(dd @ lat2) - d_nn) < 0.25:
                    p0 = xy[mag_idx[i]]
                    p1 = p0 + dd @ lat2
                    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=C_BOND, lw=1.5, zorder=2)

        # 中心胞: Br + Cr
        ax.scatter(xy[non_idx, 0], xy[non_idx, 1], s=72, c=C_NON, edgecolors="none", zorder=3)
        ax.scatter(crxy[up, 0], crxy[up, 1], s=245, c=C_UP, edgecolors="white", linewidths=1.6, zorder=4)
        ax.scatter(crxy[dn, 0], crxy[dn, 1], s=245, c=C_DN, edgecolors="white", linewidths=1.6, zorder=4)

        # 超胞边界
        ax.plot([0, a1[0], a1[0] + a2[0], a2[0], 0],
                [0, a1[1], a1[1] + a2[1], a2[1], 0], color="black", lw=1.1, ls="--", zorder=5)

        corners = np.array([[0, 0], a1, a1 + a2, a2])
        m = 3.2
        ax.set_xlim(corners[:, 0].min() - m, corners[:, 0].max() + m)
        ax.set_ylim(corners[:, 1].min() - m, corners[:, 1].max() + m)

        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(titles.get(order, order), fontsize=13)
        for sp in ax.spines.values():
            sp.set_visible(False)

    handles = [
        Line2D([], [], marker="o", ls="", mfc=C_UP, mec="white", ms=12, label=f"{args.element} 自旋 ↑"),
        Line2D([], [], marker="o", ls="", mfc=C_DN, mec="white", ms=12, label=f"{args.element} 自旋 ↓"),
        Line2D([], [], marker="o", ls="", mfc=C_NON, mec="none", ms=9, label=f"{'/'.join(args.nonmag)} (非磁)"),
        Line2D([], [], color=C_BOND, lw=1.6, label=f"{args.element}–{args.element} 最近邻键"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=10, bbox_to_anchor=(0.5, 0.045))
    fig.suptitle(f"{args.label} 单层 2×2×1 超胞：四种共线磁序俯视图", fontsize=15, y=0.965)
    fig.text(0.5, 0.018,
             f"俯视沿 c (z)。红/蓝 = {args.element} 自旋上/下；浅灰 = {'/'.join(args.nonmag)}。"
             f"数据: 02_supercell/cases/*/POSCAR + magmom.txt",
             ha="center", fontsize=8.5, color="0.45")
    fig.subplots_adjust(left=0.03, right=0.97, top=0.93, bottom=0.085, wspace=0.05, hspace=0.06)

    out = args.out or os.path.join(root, "figures", "supercell_orders.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out)
    print("写出:", os.path.abspath(out))


if __name__ == "__main__":
    main()
