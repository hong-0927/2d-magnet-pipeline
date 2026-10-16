#!/usr/bin/env python
"""render_tri_orders_pil.py -- 纯 Pillow 版三角磁序俯视图 (render_tri_orders.py 的离线替代)。

为什么存在这条路径: 本机 matplotlib 3.11.1 缺 pyparsing / cycler / contourpy /
fontTools / kiwisolver, 且无网络、缓存里也没有这些 wheel => `import matplotlib`
必然失败。Pillow 可用, 因此每条 matplotlib 渲染路径都配一条 Pillow 等价路径
(同 CLI、同输入、同输出、同配色/线型)。见 pil_fig.py 头部说明。

共线态 -> 红/蓝点 (自旋上/下); 非共线 120° -> 面内定长方向箭头
(长度 = 0.30*|a|, 与原 matplotlib 版的 quiver 同几何; 只表方向不按 |S| 缩放)。

几何/文本规则复用 render_tri_orders.py (该模块已改成不在顶层 import matplotlib):
  parse_poscar / min_image_2d / _isfloat / TITLES / C_UP / C_DN / C_NON /
  SHELL_STYLES / SHELL_TOL。

用法 (与 render_tri_orders.py 完全一致):
  python render_tri_orders_pil.py --root <proj> --element Mn --nonmag I --label MnI2 \
      -o figures/supercell_orders.png
"""
import argparse
import json
import os
import sys

import numpy as np
from PIL import Image, ImageDraw

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pil_fig as P        # noqa: E402
import render_tri_orders as R   # noqa: E402  (纯辅助函数, 顶层不 import matplotlib)

ARROW_FRAC = 0.30   # 面内自旋箭头长度 = ARROW_FRAC * |a| (只表方向, 不按 |S| 缩放)
C_ARROW = "#1a1a1a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--element", required=True)
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--label", default=None)
    ap.add_argument("--title", default=None,
                    help="整条标题覆盖（与 render_tri_orders.py 同义）。kagome 等复用本渲染器的"
                         "晶格由驱动传自己的标题，避免标题张冠李戴。")
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

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

    W = int(round(cols * 4.2 * P.DPI))
    H = int(round(rows * 4.0 * P.DPI))
    Wd, Hd = W * P.SS, H * P.SS
    img = Image.new("RGB", (Wd, Hd), "white")
    dr = ImageDraw.Draw(img)
    boxes = P.panel_grid(rows, cols, Wd, Hd)
    LEG_Y0, LEG_Y1 = 0.855 * Hd, 0.895 * Hd

    resolved = None
    n_arrow = 0
    for k, order in enumerate(orders):
        x0, y0, pw, ph = boxes[k]
        lat, symbols, coords = R.parse_poscar(os.path.join(cases_dir, order, "POSCAR"))
        lat2 = lat[:2, :2]
        with open(os.path.join(cases_dir, order, "magmom.txt"), encoding="utf-8") as f:
            magmom = np.array([float(t) for t in f.read().split("=", 1)[1].split()])
        mag_idx = [i for i, s in enumerate(symbols) if s == args.element]
        non_idx = [i for i, s in enumerate(symbols) if s in args.nonmag]
        n_mag = len(mag_idx)
        spins = magmom[:3 * n_mag].reshape(n_mag, 3)   # 每磁性原子 3 分量

        crf = coords[mag_idx, :2]
        if resolved is None:
            if shells:
                resolved = [float(s) for s in shells]
            else:
                ds = []
                for i in range(n_mag):
                    for j in range(i + 1, n_mag):
                        dd = R.min_image_2d(crf[i], crf[j])
                        ds.append(float(np.linalg.norm(dd @ lat2)))
                resolved = []
                for d in sorted(ds):
                    if not resolved or d - resolved[-1] > R.SHELL_TOL:
                        resolved.append(d)
                resolved = resolved[:3]

        a1, a2 = lat2[0], lat2[1]
        mn_xy = crf @ lat2
        corners = np.array([[0, 0], a1, a1 + a2, a2])
        T = P.panel_transform(x0, y0, pw, ph, corners, 3.0)

        # 非磁原子浅灰 (matplotlib 版 zorder=1, 画在键之下)
        if non_idx:
            nm_xy = coords[non_idx, :2] @ lat2
            r_non = P.marker_radius(45)
            for t in range(len(non_idx)):
                P.circle(dr, T(nm_xy[t]), r_non, fill=R.C_NON)

        # 壳键 (远壳先画; matplotlib 版 zorder=2, 压住非磁原子)
        for kk in range(len(resolved) - 1, -1, -1):
            _n, color, ls, lw, _alias = R.SHELL_STYLES[kk % len(R.SHELL_STYLES)]
            for i in range(n_mag):
                for j in range(i + 1, n_mag):
                    dd = R.min_image_2d(crf[i], crf[j])
                    dist = float(np.linalg.norm(dd @ lat2))
                    if abs(dist - resolved[kk]) <= R.SHELL_TOL:
                        p0 = mn_xy[i]
                        p1 = p0 + dd @ lat2
                        P.draw_seg(dr, T(p0), T(p1), color, ls, lw)

        # 磁性原子: 共线 -> 上/下点; 非共线 -> 面内箭头 (定长方向标记, 见 ARROW_FRAC)
        inplane = np.hypot(spins[:, 0], spins[:, 1])
        collinear = inplane < 1e-6
        l_arr = ARROW_FRAC * float(np.linalg.norm(a1))
        r_mag = P.marker_radius(280)
        for t in range(n_mag):
            if collinear[t]:
                if spins[t, 2] > 0:
                    P.circle(dr, T(mn_xy[t]), r_mag, fill=R.C_UP, outline="white",
                             width=P.pt2px(1.4))
                elif spins[t, 2] < 0:
                    P.circle(dr, T(mn_xy[t]), r_mag, fill=R.C_DN, outline="white",
                             width=P.pt2px(1.4))
            else:
                u = np.array([spins[t, 0], spins[t, 1]]) / max(float(inplane[t]), 1e-12)
                P.arrow(dr, T(mn_xy[t]), T(mn_xy[t] + u * l_arr), C_ARROW,
                        lw_pt=1.6, head_pt=7.0)
                n_arrow += 1

        # 超胞边界 (虚线黑框)
        box = [T(corners[0]), T(corners[1]), T(corners[2]), T(corners[3])]
        for i in range(4):
            P.draw_seg(dr, box[i], box[(i + 1) % 4], "black", "--", 1.1)

        dr.multiline_text((x0 + pw / 2.0, y0 - P.pt2px(7)), R.TITLES.get(order, order),
                          font=P.font(P.pt2px(12)), fill="black", anchor="md",
                          align="center", spacing=P.pt2px(3))

    # ---- 图例 ----
    items = P.legend_items(args.element, args.nonmag, resolved, R.SHELL_STYLES)
    if n_arrow:
        items.append(("bond", C_ARROW, "面内自旋 (120°)", ("-", 2.0)))
    P.draw_legend(dr, Wd, 0.03 * Wd, 0.97 * Wd, items, y0=LEG_Y0, y1=LEG_Y1)

    label = args.label or args.element
    foot = ("俯视沿 c (z)。红/蓝 = 自旋 ↑/↓；箭头 = 面内自旋；键按 J1/J2/J3 壳着色。\n"
            "数据: 02_supercell/cases/*/POSCAR + magmom.txt + magnetic_orders.json\n"
            "(纯 Pillow 渲染: 本机 matplotlib 缺依赖, 由 render_tri_orders_pil.py 出图)")
    P.draw_titles(dr, Wd, Hd,
                  args.title or f"{label} 三角晶格磁序俯视图（含 120° 非共线）", foot)

    img = img.resize((W, H), Image.LANCZOS)
    out = args.out or os.path.join(root, "figures", "supercell_orders.png")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    img.save(out, "PNG", optimize=True)
    print("写出:", os.path.abspath(out), "(Pillow)")


if __name__ == "__main__":
    main()
