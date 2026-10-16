#!/usr/bin/env python
"""render_supercell_pil.py -- 纯 Pillow 版 2x2x1 超胞四磁序俯视图 (render_supercell.py 的离线替代)。

为什么存在这条路径:
  本机 matplotlib 3.11.1 已装但其依赖 (pyparsing / cycler / contourpy / fontTools /
  kiwisolver) 全部缺失, 且无网络、pip/uv 缓存里也没有这些 wheel =>
  `import matplotlib` 必然失败 => render_supercell.py 永远出不了图。
  Pillow (12.3.0) 可用, 因此这里提供一条不依赖 matplotlib 的等价渲染路径:
  同样的 CLI、同样的输入文件、同样的输出文件、同样的配色与壳线型。

几何/文本规则全部复用 render_supercell.py (该模块已改成不在顶层 import matplotlib):
  parse_poscar / pairwise_distances / cluster_shells / load_order_meta / coef_label /
  SHELL_STYLES / SHELL_TOL / C_UP / C_DN / C_NON。

用法 (与 render_supercell.py 完全一致):
  python render_supercell_pil.py --root <proj> --element Mn --nonmag P S --label MnPS3 \
      -o figures/supercell_orders.png [--orders fm neel zigzag stripy] [--shells 3]
"""
import argparse
import os
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import render_supercell as R  # noqa: E402  (纯辅助函数, 顶层不 import matplotlib)

DPI = 200          # 与 matplotlib 版 figsize=(10.6, 9.8), dpi=200 对齐
SS = 2             # 超采样: 先画 SS 倍大图再 LANCZOS 缩小 => 抗锯齿
FIG_W_IN, FIG_H_IN = 10.6, 9.8
W = int(round(FIG_W_IN * DPI))
H = int(round(FIG_H_IN * DPI))

FONT_CANDIDATES = [
    r"C:\Windows\Fonts\msyh.ttc",      # 微软雅黑 (中文)
    r"C:\Windows\Fonts\msyhbd.ttc",
    r"C:\Windows\Fonts\simhei.ttf",
    r"C:\Windows\Fonts\Deng.ttf",
    r"C:\Windows\Fonts\simsun.ttc",
]
_font_cache = {}


def font(size_px):
    """按像素字号取中文字体; 找不到任何中文字体时回退 Pillow 内置位图字体。"""
    if size_px in _font_cache:
        return _font_cache[size_px]
    f = None
    for p in FONT_CANDIDATES:
        if os.path.exists(p):
            try:
                f = ImageFont.truetype(p, size_px)
                break
            except OSError:
                continue
    if f is None:
        f = ImageFont.load_default()
    _font_cache[size_px] = f
    return f


def pt2px(pt):
    """磅 -> 超采样画布上的像素 (matplotlib 的 pt 口径)。"""
    return pt / 72.0 * DPI * SS


def dash_gap(ls, lw):
    """matplotlib 默认虚线/点线模式 (lines.dashed_pattern / dotted_pattern) 按线宽缩放。"""
    if ls == "--":
        return 3.7 * lw, 1.6 * lw
    return 1.0 * lw, 1.65 * lw


def dashed_line(draw, p0, p1, color, width, dash, gap):
    x0, y0 = p0
    x1, y1 = p1
    tot = ((x1 - x0) ** 2 + (y1 - y0) ** 2) ** 0.5
    if tot <= 0:
        return
    ux, uy = (x1 - x0) / tot, (y1 - y0) / tot
    s = 0.0
    while s < tot:
        e = min(s + dash, tot)
        draw.line([(x0 + ux * s, y0 + uy * s), (x0 + ux * e, y0 + uy * e)],
                  fill=color, width=max(1, int(round(width))))
        s = e + gap


def mix(hex_color, alpha, bg=(255, 255, 255)):
    """hex 颜色按 alpha 混到背景上 (Pillow 无全局 alpha, 用等价实色近似)。"""
    c = hex_color.lstrip("#")
    rgb = tuple(int(c[i:i + 2], 16) for i in (0, 2, 4))
    return tuple(int(round(bg[k] + (rgb[k] - bg[k]) * alpha)) for k in range(3))


def circle(draw, center, radius, fill=None, outline=None, width=1):
    """实心圆 (+可选描边)。描边用"大白圆 + 内层色圆"实现: Pillow 的 ellipse(width>1)
    在同时给 fill+outline 时会画出十字状伪影, 不能用。"""
    x, y = center
    if outline is not None and width and width > 0:
        draw.ellipse([x - radius, y - radius, x + radius, y + radius], fill=outline)
        r2 = max(1.0, radius - float(width))
        draw.ellipse([x - r2, y - r2, x + r2, y + r2], fill=fill)
    else:
        draw.ellipse([x - radius, y - radius, x + radius, y + radius], fill=fill)


def marker_radius(s_pt2):
    """matplotlib scatter 的 s (磅^2, 面积) -> 半径像素。"""
    return pt2px(0.5 * float(s_pt2) ** 0.5)


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
    titles = {"fm": "FM (铁磁)", "neel": "Néel (反铁磁)", "zigzag": "zigzag", "stripy": "stripy",
              "afm_a": "AFM-a (沿 a 反平行)", "afm_b": "AFM-b (沿 b 反平行)",
              "afm_diag": "AFM-对角 (棋盘)", "afm_c": "AFM-c (C 心网 Néel, 胞内 basis 反平行)"}

    shells_angstrom, coeffs = R.load_order_meta(cases_dir)
    shell_dists = None
    if shells_angstrom:
        shell_dists = [shells_angstrom.get(k) for k in ("J1", "J2", "J3")]
        shell_dists = [d for d in shell_dists if d is not None]

    resolved = shell_dists
    spec_alt = [None] * (len(shell_dists) if shell_dists else 0)
    aligned = shell_dists is None

    Wd, Hd = W * SS, H * SS
    img = Image.new("RGB", (Wd, Hd), "white")
    dr = ImageDraw.Draw(img)

    # 版面: 4 条互不重叠的横带 (总标题 / 面板 / 图例两行 / 脚注)。
    # matplotlib 版靠 subplots_adjust + fig.legend 自适应, 这里显式给死带高, 避免中文标题压字。
    P_TOP, P_BOT = 0.10, 0.835          # 面板区 (从页顶 / 从页底算的比例)
    LEG_Y0, LEG_Y1 = 0.855, 0.895       # 图例第 1 / 第 2 行
    FOOT_Y = 0.985                      # 脚注文本块底部
    TITLE_Y = 0.018                     # 总标题 (从页顶算)

    L, Rr = 0.03 * Wd, 0.97 * Wd
    TOPY, BOTY = P_TOP * Hd, P_BOT * Hd
    wspace, hspace = 0.05, 0.10
    gapx = wspace * ((Rr - L) / 2.0)
    gapy = hspace * ((BOTY - TOPY) / 2.0)
    pw = ((Rr - L) - gapx) / 2.0
    ph = ((BOTY - TOPY) - gapy) / 2.0
    for k, order in enumerate(args.orders[:4]):
        row, col = divmod(k, 2)
        x0 = L + col * (pw + gapx)
        y0 = TOPY + row * (ph + gapy)

        d = os.path.join(cases_dir, order)
        lat, symbols, coords = R.parse_poscar(os.path.join(d, "POSCAR"))
        lat2 = lat[:2, :2]
        xy = coords[:, :2] @ lat2

        with open(os.path.join(d, "magmom.txt"), encoding="utf-8") as f:
            mom = np.array([float(t) for t in f.read().split("=", 1)[1].split()])
        mag_idx = [i for i, s in enumerate(symbols) if s == args.element]
        non_idx = [i for i, s in enumerate(symbols) if s in args.nonmag]
        # magmom.txt 逐原子 (与 POSCAR 同序) 时按索引取; 只写磁性原子时退化为前 N 个 (老行为)
        cr_mom = mom[mag_idx] if len(mom) == len(symbols) else mom[:len(mag_idx)]

        crf = coords[mag_idx, :2]
        ds = R.pairwise_distances(crf, lat2)

        if resolved is None:
            resolved = R.cluster_shells(ds)[:args.shells]
        elif not aligned:
            # 规格层是 3D 键长, 这里是俯视投影: 按次序对齐 (见 match_projected_shells)
            resolved, spec_alt = R.match_projected_shells(resolved, R.cluster_shells(ds, tol=R.SHELL_TOL))
            aligned = True

        a1, a2 = lat2[0], lat2[1]
        crxy = xy[mag_idx]
        up = [t for t, m in enumerate(cr_mom) if m > 0]
        dn = [t for t, m in enumerate(cr_mom) if m < 0]

        # 数据坐标 -> 画布像素 (等比缩放 + 居中, y 轴翻转)
        corners = np.array([[0, 0], a1, a1 + a2, a2])
        m = 3.2
        xmin, xmax = corners[:, 0].min() - m, corners[:, 0].max() + m
        ymin, ymax = corners[:, 1].min() - m, corners[:, 1].max() + m
        scale = min(pw / (xmax - xmin), ph / (ymax - ymin))
        ox = x0 + (pw - (xmax - xmin) * scale) / 2.0
        oy = y0 + ph - (ph - (ymax - ymin) * scale) / 2.0

        def T(p, _ox=ox, _oy=oy, _s=scale, _xmin=xmin, _ymin=ymin):
            return (_ox + (p[0] - _xmin) * _s, _oy - (p[1] - _ymin) * _s)

        # 背景周期像 (淡, 只画磁性元素)
        r_bg = marker_radius(150)
        for oi in (-1, 0, 1):
            for oj in (-1, 0, 1):
                if oi == 0 and oj == 0:
                    continue
                sh = oi * a1 + oj * a2
                for t in up:
                    circle(dr, T(crxy[t] + sh), r_bg, fill=mix(R.C_UP, 0.16))
                for t in dn:
                    circle(dr, T(crxy[t] + sh), r_bg, fill=mix(R.C_DN, 0.16))

        # 按壳画键: 远壳先画, 最近邻最后画 (NN 压上)
        pairs_by_shell = [[] for _ in range(len(resolved))]
        for i in range(len(mag_idx)):
            for j in range(i + 1, len(mag_idx)):
                dd = crf[j] - crf[i]
                dd -= np.round(dd)
                dist = float(np.linalg.norm(dd @ lat2))
                kk = int(np.argmin([abs(dist - s) for s in resolved]))
                if abs(dist - resolved[kk]) <= R.SHELL_TOL:
                    pairs_by_shell[kk].append((i, j, dd))
        for kk in range(len(resolved) - 1, -1, -1):
            _, color, ls, lw, _alias = R.SHELL_STYLES[kk % len(R.SHELL_STYLES)]
            width = pt2px(lw)
            dash, gap = dash_gap(ls, lw)
            for i, j, dd in pairs_by_shell[kk]:
                p0 = xy[mag_idx[i]]
                p1 = p0 + dd @ lat2
                if ls == "-":
                    dr.line([T(p0), T(p1)], fill=color, width=max(1, int(round(width))))
                else:
                    dashed_line(dr, T(p0), T(p1), color, width, pt2px(dash), pt2px(gap))

        # 中心胞: 非磁 + 磁性元素
        r_non = marker_radius(72)
        r_mag = marker_radius(245)
        for t in non_idx:
            circle(dr, T(xy[t]), r_non, fill=R.C_NON)
        for t in up:
            circle(dr, T(crxy[t]), r_mag, fill=R.C_UP, outline="white", width=pt2px(1.6))
        for t in dn:
            circle(dr, T(crxy[t]), r_mag, fill=R.C_DN, outline="white", width=pt2px(1.6))

        # 超胞边界 (虚线黑框)
        box = [T(np.array([0, 0])), T(a1), T(a1 + a2), T(a2)]
        for i in range(4):
            p0, p1 = box[i], box[(i + 1) % 4]
            dash, gap = dash_gap("--", 1.1)
            dashed_line(dr, p0, p1, "black", pt2px(1.1), pt2px(dash), pt2px(gap))

        # 面板标题 (序名 + 该序的 (J1,J2,J3) 系数)
        ttl = titles.get(order, order)
        if coeffs and order in coeffs:
            ttl += "\n" + R.coef_label(coeffs[order])
        dr.multiline_text((x0 + pw / 2.0, y0 - pt2px(7)), ttl, font=font(pt2px(12)),
                          fill="black", anchor="md", align="center", spacing=pt2px(3))

    # ---- 图例 (等价的 fig.legend(loc="lower center", ncol=4, frameon=False)) ----
    items = [("mag_up", R.C_UP, f"{args.element} 自旋 ↑", pt2px(6)),
             ("mag_dn", R.C_DN, f"{args.element} 自旋 ↓", pt2px(6)),
             ("non", R.C_NON, f"{'/'.join(args.nonmag)} (非磁)", pt2px(4.5))]
    for kk in range(len(resolved)):
        _, color, ls, lw, alias = R.SHELL_STYLES[kk % len(R.SHELL_STYLES)]
        lab = f"J{kk + 1} ({alias}) {resolved[kk]:.2f} Å"
        if kk < len(spec_alt) and spec_alt[kk]:
            lab += f" [规格层 {spec_alt[kk]:.2f} Å]"
        items.append(("bond", color, lab, (ls, lw)))

    ncol = 4
    rows = [items[i:i + ncol] for i in range(0, len(items), ncol)]
    f_leg = font(pt2px(9))
    row_ys = [LEG_Y0 * Hd, LEG_Y1 * Hd]
    for ri, row_items in enumerate(rows):
        y = row_ys[min(ri, len(row_ys) - 1)]
        colw = (Rr - L) / float(ncol)
        for ci, (kind, color, label, extra) in enumerate(row_items):
            cx = L + ci * colw + colw * 0.06
            if kind == "mag_up" or kind == "mag_dn":
                circle(dr, (cx + extra, y), extra, fill=color, outline="white", width=pt2px(1.2))
                tx = cx + 2 * extra + pt2px(4)
            elif kind == "non":
                circle(dr, (cx + extra, y), extra, fill=color)
                tx = cx + 2 * extra + pt2px(4)
            else:
                ls, lw = extra
                dash, gap = dash_gap(ls, lw)
                p0 = (cx, y)
                p1 = (cx + pt2px(26), y)
                if ls == "-":
                    dr.line([p0, p1], fill=color, width=max(1, int(round(pt2px(2.0)))))
                else:
                    dashed_line(dr, p0, p1, color, pt2px(2.0), pt2px(dash), pt2px(gap))
                tx = cx + pt2px(30)
            dr.text((tx, y), label, font=f_leg, fill=(51, 51, 51), anchor="lm")

    # ---- 总标题 + 脚注 ----
    dr.text((Wd / 2.0, TITLE_Y * Hd), f"{args.label} 单层 2×2×1 超胞：四种共线磁序俯视图",
            font=font(pt2px(15)), fill="black", anchor="ma")
    foot = (f"俯视沿 c (z)。红/蓝 = {args.element} 自旋上/下；浅灰 = {'/'.join(args.nonmag)}；"
            f"键按 J1/J2/J3 近邻壳着色。副标题为该磁序 (J1,J2,J3) 系数。\n"
            f"数据: 02_supercell/cases/*/POSCAR + magmom.txt + magnetic_orders.json")
    if any(spec_alt):
        foot += "\nC 心矩形网: 规格层壳距是 3D 键长 (J1 含面外起伏), 图中线段是俯视投影, 图例同时给出两值。"
    foot += "\n(纯 Pillow 渲染: 本机 matplotlib 缺依赖, 由 render_supercell_pil.py 出图)"
    dr.multiline_text((Wd / 2.0, FOOT_Y * Hd), foot, font=font(pt2px(8.5)),
                      fill=(115, 115, 115), anchor="md", align="center", spacing=pt2px(2))

    img = img.resize((W, H), Image.LANCZOS)
    out = args.out or os.path.join(root, "figures", "supercell_orders.png")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    img.save(out, "PNG", optimize=True)
    print("写出:", os.path.abspath(out), "(Pillow)")


if __name__ == "__main__":
    main()
