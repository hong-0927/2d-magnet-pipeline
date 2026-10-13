#!/usr/bin/env python
"""render_supercell.py -- 画 2x2x1 超胞四磁序的俯视图 (matplotlib, 无需 OVITO)。

读 02_supercell/cases/<order>/{POSCAR,magmom.txt}, 投影到 xy 俯视:
  磁性元素按自旋上/下染色 (红/蓝), 非磁元素浅灰, 画磁性元素间最近邻键 + 超胞边界。

增强 (方案 B):
  - 从 magnetic_orders.json 读 shells_angstrom (J1/J2/J3 距离) 与 coefficients_per_mag,
    把磁性元素间的键按壳 J1/J2/J3 用不同颜色+线型画出 (J1 实线 / J2 虚线 / J3 点线)。
  - 图例列出每壳距离; 每个面板副标题标注该磁序的 (J1,J2,J3) 系数。
  - 无 magnetic_orders.json 时回退为自动聚壳 (最多 --shells 个)。

用法:
  python render_supercell.py --root <proj> --element Mn --nonmag P S --label MnPS3 \
      -o figures/supercell_orders.png [--shells 3]
"""
import argparse
import json
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

# 壳 -> (颜色, 线型, 线宽, 图例名)
SHELL_STYLES = [
    ("J1", "#1f4e79", "-", 1.8, "NN"),
    ("J2", "#2ca02c", "--", 1.5, "NNN"),
    ("J3", "#ff7f0e", ":", 1.6, "NNNN"),
]

SHELL_TOL = 0.3  # 键距归壳容差 (A)


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


def pairwise_distances(crf, lat2):
    ds = []
    for i in range(len(crf)):
        for j in range(i + 1, len(crf)):
            dd = crf[j] - crf[i]
            dd -= np.round(dd)
            ds.append(float(np.linalg.norm(dd @ lat2)))
    return np.array(sorted(ds)) if ds else np.array([])


def cluster_shells(ds, tol=SHELL_TOL):
    shells = []
    for d in ds:
        if not shells or d - shells[-1][-1] > tol:
            shells.append([d])
        else:
            shells[-1].append(d)
    return [float(np.mean(s)) for s in shells]


def load_order_meta(cases_dir):
    """返回 (shells dict, coefficients dict), 缺失时为 (None, None)。

    兼容两代 magnetic_orders.json:
      - 新格式: shells_angstrom {J1,J2,J3} + coefficients_per_mag
      - 旧格式: d_NN_angstrom / d_NNN_angstrom + coefficients_per_cr
    """
    p = os.path.join(cases_dir, "magnetic_orders.json")
    if not os.path.exists(p):
        return None, None
    with open(p, encoding="utf-8") as f:
        d = json.load(f)
    shells = d.get("shells_angstrom")
    if not shells:
        nn = d.get("d_NN_angstrom")
        nnn = d.get("d_NNN_angstrom")
        if nn is not None:
            shells = {"J1": nn, "J2": nnn}
    coeffs = d.get("coefficients_per_mag") or d.get("coefficients_per_cr")
    return shells, coeffs


def coef_label(coeff):
    if not coeff:
        return ""
    names = ["J1", "J2", "J3"][:len(coeff)]
    return "(" + ",".join(names) + ")=(" + ",".join(f"{c:g}" for c in coeff) + ")"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--orders", nargs="+", default=["fm", "neel", "zigzag", "stripy"])
    ap.add_argument("--element", default="Mn")
    ap.add_argument("--nonmag", nargs="+", default=["P", "S"])
    ap.add_argument("--label", default="MnPS3")
    ap.add_argument("--shells", type=int, default=3, help="自动聚壳时最多取的壳数 (默认 3)")
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    cases_dir = os.path.join(root, "02_supercell", "cases")
    titles = {"fm": "FM (铁磁)", "neel": "Néel (反铁磁)", "zigzag": "zigzag", "stripy": "stripy"}

    # 壳距离: 优先 magnetic_orders.json, 否则自动聚壳
    shells_angstrom, coeffs = load_order_meta(cases_dir)
    if shells_angstrom:
        shell_dists = [shells_angstrom.get(k) for k in ("J1", "J2", "J3")]
        shell_dists = [d for d in shell_dists if d is not None]
    else:
        shell_dists = None

    fig, axes = plt.subplots(2, 2, figsize=(10.6, 9.8), dpi=200)
    axes = axes.ravel()

    # 跨面板统一壳距: 若未从 json 得到, 从第一个 order 的结构自动聚壳
    resolved = shell_dists
    for ax, order in zip(axes, args.orders):
        d = os.path.join(cases_dir, order)
        lat, symbols, coords = parse_poscar(os.path.join(d, "POSCAR"))
        lat2 = lat[:2, :2]
        xy = coords[:, :2] @ lat2

        mom = np.array([float(t) for t in
                        open(os.path.join(d, "magmom.txt"), encoding="utf-8").read().split("=", 1)[1].split()])
        mag_idx = [i for i, s in enumerate(symbols) if s == args.element]
        non_idx = [i for i, s in enumerate(symbols) if s in args.nonmag]
        cr_mom = mom[:len(mag_idx)]

        crf = coords[mag_idx, :2]
        ds = pairwise_distances(crf, lat2)
        d_nn = min(ds) if len(ds) else 0.0

        if resolved is None:
            resolved = cluster_shells(ds)[:args.shells]

        a1, a2 = lat2[0], lat2[1]
        crxy = xy[mag_idx]
        up = [k for k, m in enumerate(cr_mom) if m > 0]
        dn = [k for k, m in enumerate(cr_mom) if m < 0]

        # 背景周期像 (淡, 只画磁性元素)
        for oi in (-1, 0, 1):
            for oj in (-1, 0, 1):
                if oi == 0 and oj == 0:
                    continue
                sh = oi * a1 + oj * a2
                ax.scatter(crxy[up, 0] + sh[0], crxy[up, 1] + sh[1],
                           s=150, c=C_UP, alpha=0.16, edgecolors="none", zorder=1)
                ax.scatter(crxy[dn, 0] + sh[0], crxy[dn, 1] + sh[1],
                           s=150, c=C_DN, alpha=0.16, edgecolors="none", zorder=1)

        # 按壳画键: 远壳先画, 最近邻最后画 (NN 压上)
        pairs_by_shell = [[] for _ in range(len(resolved))]
        for i in range(len(mag_idx)):
            for j in range(i + 1, len(mag_idx)):
                dd = crf[j] - crf[i]
                dd -= np.round(dd)
                dist = float(np.linalg.norm(dd @ lat2))
                k = int(np.argmin([abs(dist - s) for s in resolved]))
                if abs(dist - resolved[k]) <= SHELL_TOL:
                    pairs_by_shell[k].append((i, j, dd))
        for k in range(len(resolved) - 1, -1, -1):
            _, color, ls, lw, _alias = SHELL_STYLES[k % len(SHELL_STYLES)]
            for i, j, dd in pairs_by_shell[k]:
                p0 = xy[mag_idx[i]]
                p1 = p0 + dd @ lat2
                ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=color, ls=ls, lw=lw, zorder=2)

        # 中心胞: 非磁 + 磁性元素
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
        ttl = titles.get(order, order)
        if coeffs and order in coeffs:
            ttl += "\n" + coef_label(coeffs[order])
        ax.set_title(ttl, fontsize=12)
        for sp in ax.spines.values():
            sp.set_visible(False)

    handles = [
        Line2D([], [], marker="o", ls="", mfc=C_UP, mec="white", ms=12, label=f"{args.element} 自旋 ↑"),
        Line2D([], [], marker="o", ls="", mfc=C_DN, mec="white", ms=12, label=f"{args.element} 自旋 ↓"),
        Line2D([], [], marker="o", ls="", mfc=C_NON, mec="none", ms=9, label=f"{'/'.join(args.nonmag)} (非磁)"),
    ]
    for k in range(len(resolved)):
        _, color, ls, lw, alias = SHELL_STYLES[k % len(SHELL_STYLES)]
        name = f"J{k + 1}" if k < 3 else f"J{k + 1}"
        handles.append(Line2D([], [], color=color, ls=ls, lw=2.0,
                              label=f"{name} ({alias}) {resolved[k]:.2f} Å"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 0.055))
    fig.suptitle(f"{args.label} 单层 2×2×1 超胞：四种共线磁序俯视图", fontsize=15, y=0.975)
    fig.text(0.5, 0.018,
             f"俯视沿 c (z)。红/蓝 = {args.element} 自旋上/下；浅灰 = {'/'.join(args.nonmag)}；"
             f"键按 J1/J2/J3 近邻壳着色。副标题为该磁序 (J1,J2,J3) 系数。\n"
             f"数据: 02_supercell/cases/*/POSCAR + magmom.txt + magnetic_orders.json",
             ha="center", fontsize=8.5, color="0.45")
    fig.subplots_adjust(left=0.03, right=0.97, top=0.94, bottom=0.115, wspace=0.05, hspace=0.10)

    out = args.out or os.path.join(root, "figures", "supercell_orders.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out)
    print("写出:", os.path.abspath(out))


if __name__ == "__main__":
    main()
