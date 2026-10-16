#!/usr/bin/env python
"""render_cluster_orders.py -- 簇磁体磁序俯视图 (matplotlib, 无需 OVITO)。

簇磁体（如 Nb3X8）把每个磁性簇（Nb3）当作 1 个有效自旋。本脚本:
  - 自动扫描 02_supercell/cases/<order>/{POSCAR,magmom.txt}（order 数可变: fm + 波矢态）;
  - 按 --cutoff 把磁性原子聚成簇, 取簇净矩符号给质心着色 (红=↑, 蓝=↓);
  - 画簇心之间 J1/J2/J3 壳键 (距离读 magnetic_orders.json 的 inter_shells_angstrom);
  - 非磁原子浅灰背景。

用法:
  python render_cluster_orders.py --root <proj> --element Nb --nonmag I --cutoff 3.3 \
      --label Nb3I8 -o figures/supercell_orders.png
"""
import argparse
import json
import os

import numpy as np


def _plt():
    """延迟 import matplotlib: 让本模块的几何辅助函数 (parse_poscar / find_clusters /
    cluster_centroids / SHELL_STYLES ...) 能在**没有 matplotlib** 的环境里被
    render_cluster_orders_pil.py 复用 (本机 matplotlib 缺依赖, 见该文件头)。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False
    return plt, Line2D

C_UP = "#d62728"
C_DN = "#1f77b4"
C_NON = "#c8c8c8"

SHELL_STYLES = [
    ("J1", "#1f4e79", "-", 1.8, "NN"),
    ("J2", "#2ca02c", "--", 1.5, "NNN"),
    ("J3", "#ff7f0e", ":", 1.6, "NNNN"),
]
SHELL_TOL = 0.4  # 簇心键距归壳容差 (A)


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


def _isfloat(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def min_image_frac(c1, c2):
    d = c2 - c1
    d -= np.floor(d + 0.5)
    return d


def min_image_2d(c1, c2):
    """面内 (x,y) 分数坐标最小镜像差。"""
    d = c2[:2] - c1[:2]
    d -= np.round(d)
    return d


def find_clusters(atoms, coords, lat, element, cutoff):
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
    cents = []
    for g in groups:
        ref = coords[g[0]]
        cart = []
        for gi in g:
            cart.append((ref + min_image_frac(ref, coords[gi])) @ lat)
        c = np.mean(cart, axis=0)
        frac = c @ np.linalg.inv(lat)
        frac -= np.floor(frac)
        cents.append(frac)
    return np.array(cents)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--element", required=True)
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--cutoff", type=float, default=3.3)
    ap.add_argument("--label", default=None)
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    plt, Line2D = _plt()
    root = os.path.abspath(args.root)
    cases_dir = os.path.join(root, "02_supercell", "cases")

    shells = None
    mp = os.path.join(cases_dir, "magnetic_orders.json")
    if os.path.exists(mp):
        with open(mp, encoding="utf-8") as f:
            meta = json.load(f)
        shells = meta.get("inter_shells_angstrom")

    orders = []
    for d in sorted(os.listdir(cases_dir)):
        p = os.path.join(cases_dir, d)
        if os.path.isdir(p) and os.path.exists(os.path.join(p, "POSCAR")) \
                and os.path.exists(os.path.join(p, "magmom.txt")):
            orders.append(d)
    # fm 置顶
    orders = sorted(orders, key=lambda o: (o != "fm", o))
    n = len(orders)
    cols = int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4.2, rows * 4.0), dpi=200)
    axes = np.atleast_1d(axes).ravel()

    titles = {"fm": "FM (铁磁)", "fm_trim": "fm_trim (S=1/2)"}

    resolved = None
    for ax, order in zip(axes, orders):
        lat, symbols, coords = parse_poscar(os.path.join(cases_dir, order, "POSCAR"))
        lat2 = lat[:2, :2]
        magmom = np.array([float(t) for t in
                           open(os.path.join(cases_dir, order, "magmom.txt"), encoding="utf-8")
                          .read().split("=", 1)[1].split()])
        mag_idx = [i for i, s in enumerate(symbols) if s == args.element]
        non_idx = [i for i, s in enumerate(symbols) if s in args.nonmag]
        mom_mag = magmom[:len(mag_idx)]  # magmom 与 POSCAR 顺序一致(按物种分组)

        groups = find_clusters(symbols, coords, lat, args.element, args.cutoff)
        cents = cluster_centroids(coords, lat, groups)
        net = np.array([sum(mom_mag[[mag_idx.index(a) for a in g]]) for g in groups])

        # 统一壳距: 优先 json, 否则按簇心距自动聚壳
        if resolved is None:
            if shells:
                resolved = [float(s) for s in shells]
            else:
                ds = []
                for i in range(len(cents)):
                    for j in range(i + 1, len(cents)):
                        dd = min_image_2d(cents[i], cents[j])
                        ds.append(float(np.linalg.norm(dd @ lat2)))
                resolved = []
                for d in sorted(ds):
                    if not resolved or d - resolved[-1] > SHELL_TOL:
                        resolved.append(d)
                resolved = resolved[:3]

        a1, a2 = lat2[0], lat2[1]
        cen_xy = cents[:, :2] @ lat2

        # 簇心间壳键
        pairs_by_shell = [[] for _ in range(len(resolved))]
        for i in range(len(cents)):
            for j in range(i + 1, len(cents)):
                dd = min_image_2d(cents[i], cents[j])
                dist = float(np.linalg.norm(dd @ lat2))
                k = int(np.argmin([abs(dist - s) for s in resolved]))
                if abs(dist - resolved[k]) <= SHELL_TOL:
                    pairs_by_shell[k].append((i, j, dd))
        for k in range(len(resolved) - 1, -1, -1):
            _, color, ls, lw, _alias = SHELL_STYLES[k % len(SHELL_STYLES)]
            for i, j, dd in pairs_by_shell[k]:
                p0 = cen_xy[i]
                p1 = p0 + dd @ lat2
                ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, ls=ls, lw=lw, zorder=2)

        # 非磁原子浅灰
        nm_xy = coords[non_idx, :2] @ lat2
        if len(non_idx):
            ax.scatter(nm_xy[:, 0], nm_xy[:, 1], s=45, c=C_NON, edgecolors="none", zorder=1)

        # 簇心按净矩符号着色
        up = [i for i, v in enumerate(net) if v > 0]
        dn = [i for i, v in enumerate(net) if v < 0]
        ax.scatter(cen_xy[up, 0], cen_xy[up, 1], s=420, c=C_UP, edgecolors="white",
                   linewidths=1.6, zorder=4)
        ax.scatter(cen_xy[dn, 0], cen_xy[dn, 1], s=420, c=C_DN, edgecolors="white",
                   linewidths=1.6, zorder=4)

        ax.plot([0, a1[0], a1[0] + a2[0], a2[0], 0],
                [0, a1[1], a1[1] + a2[1], a2[1], 0], color="black", lw=1.0, ls="--", zorder=5)
        corners = np.array([[0, 0], a1, a1 + a2, a2])
        m = 3.0
        ax.set_xlim(corners[:, 0].min() - m, corners[:, 0].max() + m)
        ax.set_ylim(corners[:, 1].min() - m, corners[:, 1].max() + m)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(titles.get(order, order), fontsize=11)
        for sp in ax.spines.values():
            sp.set_visible(False)

    for ax in axes[n:]:
        ax.axis("off")

    handles = [
        Line2D([], [], marker="o", ls="", mfc=C_UP, mec="white", ms=12, label=f"{args.element} 簇 自旋 ↑"),
        Line2D([], [], marker="o", ls="", mfc=C_DN, mec="white", ms=12, label=f"{args.element} 簇 自旋 ↓"),
        Line2D([], [], marker="o", ls="", mfc=C_NON, mec="none", ms=9,
               label=f"{'/'.join(args.nonmag)} (非磁)" if args.nonmag else ""),
    ]
    handles = [h for h in handles if h.get_label()]
    for k in range(len(resolved)):
        _, color, ls, lw, alias = SHELL_STYLES[k % len(SHELL_STYLES)]
        name = f"J{k + 1}"
        handles.append(Line2D([], [], color=color, ls=ls, lw=2.0,
                              label=f"{name} ({alias}) {resolved[k]:.2f} Å"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False,
               fontsize=8.5, bbox_to_anchor=(0.5, 0.02))
    label = args.label or args.element
    fig.suptitle(f"{label} 簇磁体磁序俯视图（每簇 = 1 个有效自旋）", fontsize=14, y=0.99)
    fig.text(0.5, 0.012,
             f"俯视沿 c (z)。红/蓝 = 簇净矩 ↑/↓；键按 J1/J2/J3 簇间壳着色。"
             f"数据: 02_supercell/cases/*/POSCAR + magmom.txt + magnetic_orders.json",
             ha="center", fontsize=8, color="0.45")
    fig.subplots_adjust(left=0.03, right=0.97, top=0.93, bottom=0.10, wspace=0.06, hspace=0.12)

    out = args.out or os.path.join(root, "figures", "supercell_orders.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out)
    print("写出:", os.path.abspath(out))


if __name__ == "__main__":
    main()
