#!/usr/bin/env python
"""render_cluster_orders_pil.py -- 纯 Pillow 版簇磁体磁序俯视图
(render_cluster_orders.py 的离线替代)。

为什么存在这条路径: 本机 matplotlib 3.11.1 缺 pyparsing / cycler / contourpy /
fontTools / kiwisolver, 且无网络、缓存里也没有这些 wheel => `import matplotlib`
必然失败。Pillow 可用, 因此每条 matplotlib 渲染路径都配一条 Pillow 等价路径
(同 CLI、同输入、同输出、同配色/线型)。见 pil_fig.py 头部说明。

几何/文本规则复用 render_cluster_orders.py (该模块已改成不在顶层 import matplotlib):
  parse_poscar / min_image_2d / find_clusters / cluster_centroids / _isfloat /
  C_UP / C_DN / C_NON / SHELL_STYLES / SHELL_TOL。

用法 (与 render_cluster_orders.py 完全一致):
  python render_cluster_orders_pil.py --root <proj> --element Nb --nonmag I \
      --cutoff 3.3 --label Nb3I8 -o figures/supercell_orders.png
"""
import argparse
import os
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pil_fig as P          # noqa: E402
import render_cluster_orders as R  # noqa: E402  (纯辅助函数, 顶层不 import matplotlib)

TITLES = {"fm": "FM (铁磁)", "fm_trim": "fm_trim (S=1/2)"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--element", required=True)
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--cutoff", type=float, default=3.3)
    ap.add_argument("--label", default=None)
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    cases_dir = os.path.join(root, "02_supercell", "cases")

    shells = None
    mp = os.path.join(cases_dir, "magnetic_orders.json")
    if os.path.exists(mp):
        import json
        with open(mp, encoding="utf-8") as f:
            meta = json.load(f)
        shells = meta.get("inter_shells_angstrom")

    orders = []
    for d in sorted(os.listdir(cases_dir)):
        p = os.path.join(cases_dir, d)
        if os.path.isdir(p) and os.path.exists(os.path.join(p, "POSCAR")) \
                and os.path.exists(os.path.join(p, "magmom.txt")):
            orders.append(d)
    orders = sorted(orders, key=lambda o: (o != "fm", o))  # fm 置顶
    n = len(orders)
    cols = int(np.ceil(np.sqrt(n)))
    rows = int(np.ceil(n / cols))

    W = int(round(cols * 4.2 * P.DPI))
    H = int(round(rows * 4.0 * P.DPI))
    Wd, Hd = W * P.SS, H * P.SS
    img = Image.new("RGB", (Wd, Hd), "white")
    dr = ImageDraw.Draw(img)

    boxes = P.panel_grid(rows, cols, Wd, Hd)
    LEG_Y0, LEG_Y1 = 0.855 * Hd, 0.895 * Hd

    resolved = None
    for k, order in enumerate(orders):
        x0, y0, pw, ph = boxes[k]
        lat, symbols, coords = R.parse_poscar(os.path.join(cases_dir, order, "POSCAR"))
        lat2 = lat[:2, :2]
        with open(os.path.join(cases_dir, order, "magmom.txt"), encoding="utf-8") as f:
            magmom = np.array([float(t) for t in f.read().split("=", 1)[1].split()])
        mag_idx = [i for i, s in enumerate(symbols) if s == args.element]
        non_idx = [i for i, s in enumerate(symbols) if s in args.nonmag]
        mom_mag = magmom[:len(mag_idx)]  # magmom 与 POSCAR 顺序一致(按物种分组)

        groups = R.find_clusters(symbols, coords, lat, args.element, args.cutoff)
        cents = R.cluster_centroids(coords, lat, groups)
        net = np.array([sum(mom_mag[[mag_idx.index(a) for a in g]]) for g in groups])

        # 统一壳距: 优先 json 的 inter_shells_angstrom, 否则按簇心距自动聚壳
        if resolved is None:
            if shells:
                resolved = [float(s) for s in shells]
            else:
                ds = []
                for i in range(len(cents)):
                    for j in range(i + 1, len(cents)):
                        dd = R.min_image_2d(cents[i], cents[j])
                        ds.append(float(np.linalg.norm(dd @ lat2)))
                resolved = []
                for d in sorted(ds):
                    if not resolved or d - resolved[-1] > R.SHELL_TOL:
                        resolved.append(d)
                resolved = resolved[:3]

        a1, a2 = lat2[0], lat2[1]
        cen_xy = np.asarray(cents)[:, :2] @ lat2
        corners = np.array([[0, 0], a1, a1 + a2, a2])
        T = P.panel_transform(x0, y0, pw, ph, corners, 3.0)

        # 非磁原子浅灰 (matplotlib 版 zorder=1, 画在键之下)
        if len(non_idx):
            nm_xy = coords[non_idx, :2] @ lat2
            r_non = P.marker_radius(45)
            for t in range(len(non_idx)):
                P.circle(dr, T(nm_xy[t]), r_non, fill=R.C_NON)

        # 簇心间壳键 (matplotlib 版 zorder=2, 压住非磁原子)
        pairs_by_shell = [[] for _ in range(len(resolved))]
        for i in range(len(cents)):
            for j in range(i + 1, len(cents)):
                dd = R.min_image_2d(cents[i], cents[j])
                dist = float(np.linalg.norm(dd @ lat2))
                kk = int(np.argmin([abs(dist - s) for s in resolved]))
                if abs(dist - resolved[kk]) <= R.SHELL_TOL:
                    pairs_by_shell[kk].append((i, j, dd))
        for kk in range(len(resolved) - 1, -1, -1):
            _n, color, ls, lw, _alias = R.SHELL_STYLES[kk % len(R.SHELL_STYLES)]
            for i, j, dd in pairs_by_shell[kk]:
                p0 = cen_xy[i]
                p1 = p0 + dd @ lat2
                P.draw_seg(dr, T(p0), T(p1), color, ls, lw)

        # 簇心按净矩符号着色
        r_cl = P.marker_radius(420)
        for i in range(len(net)):
            if net[i] > 0:
                P.circle(dr, T(cen_xy[i]), r_cl, fill=R.C_UP, outline="white",
                         width=P.pt2px(1.6))
            elif net[i] < 0:
                P.circle(dr, T(cen_xy[i]), r_cl, fill=R.C_DN, outline="white",
                         width=P.pt2px(1.6))

        # 超胞边界 (虚线黑框)
        box = [T(corners[0]), T(corners[1]), T(corners[2]), T(corners[3])]
        for i in range(4):
            P.draw_seg(dr, box[i], box[(i + 1) % 4], "black", "--", 1.1)

        # 面板标题
        dr.multiline_text((x0 + pw / 2.0, y0 - P.pt2px(7)), TITLES.get(order, order),
                          font=P.font(P.pt2px(12)), fill="black", anchor="md",
                          align="center", spacing=P.pt2px(3))

    # ---- 图例 ----
    items = P.legend_items(f"{args.element} 簇", args.nonmag, resolved, R.SHELL_STYLES)
    P.draw_legend(dr, Wd, 0.03 * Wd, 0.97 * Wd, items, y0=LEG_Y0, y1=LEG_Y1)

    label = args.label or args.element
    foot = ("俯视沿 c (z)。红/蓝 = 簇净矩 ↑/↓；键按 J1/J2/J3 簇间壳着色。"
            "每簇 = 1 个有效自旋。\n"
            "数据: 02_supercell/cases/*/POSCAR + magmom.txt + magnetic_orders.json\n"
            "(纯 Pillow 渲染: 本机 matplotlib 缺依赖, 由 render_cluster_orders_pil.py 出图)")
    P.draw_titles(dr, Wd, Hd, f"{label} 簇磁体磁序俯视图（每簇 = 1 个有效自旋）", foot)

    img = img.resize((W, H), Image.LANCZOS)
    out = args.out or os.path.join(root, "figures", "supercell_orders.png")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    img.save(out, "PNG", optimize=True)
    print("写出:", os.path.abspath(out), "(Pillow)")


if __name__ == "__main__":
    main()
