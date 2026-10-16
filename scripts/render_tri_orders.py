#!/usr/bin/env python
"""render_tri_orders.py -- 三角磁序俯视图 (matplotlib, 无需 OVITO)。

读 02_supercell/cases/<order>/{POSCAR,magmom.txt}:
  - 共线态 (fm/stripy/uud): 磁性原子按自旋上/下染色 (红/蓝);
  - 非共线 120°: 磁性原子画面内自旋矢量箭头;
  - 画 J1/J2/J3 壳键 (距离读 magnetic_orders.json 的 shells_angstrom);
  - 非磁原子浅灰。

用法:
  python render_tri_orders.py --root <proj> --element Mn --nonmag I --label MnI2 \
      -o figures/supercell_orders.png
"""
import argparse
import json
import os

import numpy as np


def _plt():
    """延迟 import matplotlib: 让本模块的几何辅助函数 (parse_poscar / min_image_2d /
    SHELL_STYLES / TITLES ...) 能在**没有 matplotlib** 的环境里被
    render_tri_orders_pil.py 复用 (本机 matplotlib 缺依赖, 见该文件头)。"""
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
SHELL_TOL = 0.4

TITLES = {"fm": "FM (铁磁)", "120deg": "120° (非共线)", "stripy_a": "stripy_a",
          "stripy_b": "stripy_b", "stripy_ab": "stripy_ab", "uud": "uud"}


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


def min_image_2d(c1, c2):
    d = c2[:2] - c1[:2]
    d -= np.round(d)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--element", required=True)
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--label", default=None)
    ap.add_argument("--title", default=None,
                    help="整条标题覆盖（默认 '<label> 三角晶格磁序俯视图（含 120° 非共线）'）。"
                         "kagome 等复用本渲染器的晶格由驱动传自己的标题，避免标题张冠李戴。")
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
        sa = meta.get("shells_angstrom")
        if sa:
            shells = [sa.get(k) for k in ("J1", "J2", "J3")]
            shells = [s for s in shells if s is not None]

    orders = sorted(d for d in os.listdir(cases_dir)
                    if os.path.isdir(os.path.join(cases_dir, d))
                    and os.path.exists(os.path.join(cases_dir, d, "POSCAR")))
    n = len(orders)
    cols = int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4.2, rows * 4.0), dpi=200)
    axes = np.atleast_1d(axes).ravel()

    resolved = None
    for ax, order in zip(axes, orders):
        lat, symbols, coords = parse_poscar(os.path.join(cases_dir, order, "POSCAR"))
        lat2 = lat[:2, :2]
        magmom = np.array([float(t) for t in
                           open(os.path.join(cases_dir, order, "magmom.txt"), encoding="utf-8")
                          .read().split("=", 1)[1].split()])
        mag_idx = [i for i, s in enumerate(symbols) if s == args.element]
        non_idx = [i for i, s in enumerate(symbols) if s in args.nonmag]
        n_mag = len(mag_idx)
        spins = magmom[:3 * n_mag].reshape(n_mag, 3)  # 每 Mn 3 分量

        crf = coords[mag_idx, :2]
        if resolved is None:
            if shells:
                resolved = [float(s) for s in shells]
            else:
                ds = []
                for i in range(n_mag):
                    for j in range(i + 1, n_mag):
                        dd = min_image_2d(crf[i], crf[j])
                        ds.append(float(np.linalg.norm(dd @ lat2)))
                resolved = []
                for d in sorted(ds):
                    if not resolved or d - resolved[-1] > SHELL_TOL:
                        resolved.append(d)
                resolved = resolved[:3]

        a1, a2 = lat2[0], lat2[1]
        mn_xy = crf @ lat2

        # 壳键
        for k in range(len(resolved) - 1, -1, -1):
            _, color, ls, lw, _alias = SHELL_STYLES[k % len(SHELL_STYLES)]
            for i in range(n_mag):
                for j in range(i + 1, n_mag):
                    dd = min_image_2d(crf[i], crf[j])
                    dist = float(np.linalg.norm(dd @ lat2))
                    if abs(dist - resolved[k]) <= SHELL_TOL:
                        p0 = mn_xy[i]
                        p1 = p0 + dd @ lat2
                        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, ls=ls, lw=lw, zorder=2)

        # 非磁原子
        if non_idx:
            nm_xy = coords[non_idx, :2] @ lat2
            ax.scatter(nm_xy[:, 0], nm_xy[:, 1], s=45, c=C_NON, edgecolors="none", zorder=1)

        # 磁性原子: 共线 -> 上/下点; 非共线 -> 面内箭头
        inplane = np.hypot(spins[:, 0], spins[:, 1])
        collinear = inplane < 1e-6
        up = np.where(collinear & (spins[:, 2] > 0))[0]
        dn = np.where(collinear & (spins[:, 2] < 0))[0]
        if len(up):
            ax.scatter(mn_xy[up, 0], mn_xy[up, 1], s=280, c=C_UP, edgecolors="white",
                       linewidths=1.4, zorder=4)
        if len(dn):
            ax.scatter(mn_xy[dn, 0], mn_xy[dn, 1], s=280, c=C_DN, edgecolors="white",
                       linewidths=1.4, zorder=4)
        ncz = np.where(~collinear)[0]
        if len(ncz):
            # 面内自旋只画**方向** (统一长度 0.30*|a|), 不按 |S| 缩放:
            # 原作 quiver(scale_units="xy", scale=0.9) 在 |S|=5 时会画出 5/0.9 ≈ 5.6 A
            # 的长箭头, 在 a≈4.2 A 的三角晶格上互相纠缠。PIL 版同几何 (见 render_tri_orders_pil.py)。
            l_arr = 0.30 * float(np.linalg.norm(a1))
            for t in ncz:
                u = np.array([spins[t, 0], spins[t, 1]]) / max(float(inplane[t]), 1e-12)
                ax.annotate("", xy=mn_xy[t] + u * l_arr, xytext=mn_xy[t],
                            arrowprops=dict(arrowstyle="-|>", color="#1a1a1a", lw=1.4,
                                            mutation_scale=12), zorder=5)

        ax.plot([0, a1[0], a1[0] + a2[0], a2[0], 0],
                [0, a1[1], a1[1] + a2[1], a2[1], 0], color="black", lw=1.0, ls="--", zorder=5)
        corners = np.array([[0, 0], a1, a1 + a2, a2])
        m = 3.0
        ax.set_xlim(corners[:, 0].min() - m, corners[:, 0].max() + m)
        ax.set_ylim(corners[:, 1].min() - m, corners[:, 1].max() + m)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(TITLES.get(order, order), fontsize=11)
        for sp in ax.spines.values():
            sp.set_visible(False)

    for ax in axes[n:]:
        ax.axis("off")

    handles = [
        Line2D([], [], marker="o", ls="", mfc=C_UP, mec="white", ms=12, label=f"{args.element} 自旋 ↑"),
        Line2D([], [], marker="o", ls="", mfc=C_DN, mec="white", ms=12, label=f"{args.element} 自旋 ↓"),
    ]
    if args.nonmag:
        handles.append(Line2D([], [], marker="o", ls="", mfc=C_NON, mec="none", ms=9,
                              label=f"{'/'.join(args.nonmag)} (非磁)"))
    for k in range(len(resolved)):
        _, color, ls, lw, alias = SHELL_STYLES[k % len(SHELL_STYLES)]
        handles.append(Line2D([], [], color=color, ls=ls, lw=2.0,
                              label=f"J{k + 1} ({alias}) {resolved[k]:.2f} Å"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False,
               fontsize=8.5, bbox_to_anchor=(0.5, 0.02))
    label = args.label or args.element
    fig.suptitle(args.title or f"{label} 三角晶格磁序俯视图（含 120° 非共线）", fontsize=14, y=0.99)
    fig.text(0.5, 0.012,
             f"俯视沿 c (z)。红/蓝 = 自旋 ↑/↓；箭头 = 面内自旋；键按 J1/J2/J3 壳着色。"
             f"数据: 02_supercell/cases/*/POSCAR + magmom.txt + magnetic_orders.json",
             ha="center", fontsize=8, color="0.45")
    fig.subplots_adjust(left=0.03, right=0.97, top=0.93, bottom=0.10, wspace=0.06, hspace=0.12)

    out = args.out or os.path.join(root, "figures", "supercell_orders.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out)
    print("写出:", os.path.abspath(out))


if __name__ == "__main__":
    main()
