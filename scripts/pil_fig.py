#!/usr/bin/env python
"""pil_fig.py -- Pillow 兜底渲染器共用的纯 Pillow 绘图原语。

为什么有这条路径: 本机 matplotlib 3.11.1 已装但其依赖 (pyparsing / cycler /
contourpy / fontTools / kiwisolver) 全部缺失, 且无网络、pip/uv 缓存里也没有这些
wheel => `import matplotlib` 必然失败。Pillow (12.3.0) 可用, 于是本仓库为每个
matplotlib 渲染器配一条 Pillow 等价路径:

    render_supercell.py       -> render_supercell_pil.py
    render_cluster_orders.py  -> render_cluster_orders_pil.py
    render_tri_orders.py      -> render_tri_orders_pil.py
    mc_2d.py 的 mc_mt.png     -> mc_2d.py 内部 _plot_mt_png()
    mc_fss.py 的 fss.png      -> 本模块 chart_* 原语 + mc_fss.py 内部 _plot_fss_png()
    build_fgt_cases.py 的 FGT 图 -> 本模块 chart_* 原语 + build_fgt_cases.py 内部 _plot_fgt_png()

render_supercell_pil.py 早于本模块, 自带一份等价拷贝; 之后新增的 Pillow 渲染器
一律复用本模块 (render_cluster_orders_pil.py / render_tri_orders_pil.py 即是)。
mc_2d.py 按设计要自包含(可单独拷走运行), 因此它内部另有一份极简字体加载与折线绘制。

约定 (三个渲染器完全一致):
  - 画布 = figsize(inch) * DPI * SS, 画完再 LANCZOS 缩到 figsize*DPI => 抗锯齿;
  - 字号/线宽用 pt2px(pt) 换算, 口径与 matplotlib 的 pt 相同;
  - 颜色/壳线型仍从各自 matplotlib 模块 import (那里是单一事实来源);
  - Pillow 的 ellipse(width>1) 在同时给 fill+outline 时会画十字伪影, 一律用
    circle() 的 "大白圆 + 内层色圆" 两段式描边。
"""
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont  # noqa: F401  (Image 供调用方使用)

DPI = 200          # 与 matplotlib 版 figsize*, dpi=200 对齐
SS = 2             # 超采样倍率

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
                f = ImageFont.truetype(p, int(round(size_px)))
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
    """matplotlib 默认虚线/点线模式 (lines.dashed_pattern / dotted_pattern) 按线宽缩放。

    返回 (dash, gap), 单位磅。
    """
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


def draw_seg(draw, p0, p1, color, ls, lw):
    """按 matplotlib 线型字符串画一段 (实线/虚线/点线)。"""
    width = pt2px(lw)
    if ls == "-":
        draw.line([p0, p1], fill=color, width=max(1, int(round(width))))
    else:
        dash, gap = dash_gap(ls, lw)
        dashed_line(draw, p0, p1, color, width, pt2px(dash), pt2px(gap))


def arrow(draw, p0, p1, color, lw_pt=1.4, head_pt=5.0):
    """带箭头的面内自旋矢量 (等价 matplotlib quiver)。"""
    x0, y0 = p0
    x1, y1 = p1
    draw.line([p0, p1], fill=color, width=max(1, int(round(pt2px(lw_pt)))))
    dx, dy = x1 - x0, y1 - y0
    tot = (dx * dx + dy * dy) ** 0.5
    if tot <= 0:
        return
    ux, uy = dx / tot, dy / tot
    hl = pt2px(head_pt)
    hw = hl * 0.5
    bx, by = x1 - ux * hl, y1 - uy * hl
    px, py = -uy, ux
    draw.polygon([(x1, y1), (bx + px * hw, by + py * hw), (bx - px * hw, by - py * hw)],
                 fill=color)


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


def panel_grid(rows, cols, Wd, Hd, left=0.03, right=0.97, top=0.10, bottom=0.835,
               wspace=0.05, hspace=0.10):
    """把画布切 rows*cols 个面板, 返回按行优先排列的 (x0, y0, pw, ph) 列表。

    top/bottom/left/right 是相对画布的比例; 面板区之外的横带留给总标题/图例/脚注。
    """
    L, R = left * Wd, right * Wd
    TOPY, BOTY = top * Hd, bottom * Hd
    gapx = wspace * ((R - L) / float(cols))
    gapy = hspace * ((BOTY - TOPY) / float(rows))
    pw = ((R - L) - (cols - 1) * gapx) / float(cols)
    ph = ((BOTY - TOPY) - (rows - 1) * gapy) / float(rows)
    boxes = []
    for r in range(rows):
        for c in range(cols):
            boxes.append((L + c * (pw + gapx), TOPY + r * (ph + gapy), pw, ph))
    return boxes


def panel_transform(x0, y0, pw, ph, corners, margin):
    """构造 数据坐标(面内 xy, 单位 A) -> 画布像素 的变换 (等比缩放+居中, y 翻转)。"""
    corners = np.asarray(corners, dtype=float)
    xmin, xmax = corners[:, 0].min() - margin, corners[:, 0].max() + margin
    ymin, ymax = corners[:, 1].min() - margin, corners[:, 1].max() + margin
    scale = min(pw / (xmax - xmin), ph / (ymax - ymin))
    ox = x0 + (pw - (xmax - xmin) * scale) / 2.0
    oy = y0 + ph - (ph - (ymax - ymin) * scale) / 2.0

    def T(p):
        return (ox + (p[0] - xmin) * scale, oy - (p[1] - ymin) * scale)

    return T


def draw_legend(draw, Wd, L, R, items, y0=0.855, y1=0.895, ncol=4, sep=0.06):
    """等价 fig.legend(loc="lower center", ncol=n, frameon=False)。

    y0 / y1 是**画布像素**的 y 坐标 (第 1 / 第 2 行); 超过 2 行时后面的行沿用 y1。
    items: (kind, color, label, extra) 序列, kind ∈ {"mag_up","mag_dn","non","bond"};
      - mag_*: extra = 圆半径像素 (白描边);
      - non  : extra = 圆半径像素;
      - bond : extra = (ls, lw_pt), 画一小段线型样例。
    """
    rows = [items[i:i + ncol] for i in range(0, len(items), ncol)]
    f_leg = font(pt2px(9))
    for ri, row_items in enumerate(rows):
        y = (y0 if ri == 0 else y1)
        colw = (R - L) / float(ncol)
        for ci, (kind, color, label, extra) in enumerate(row_items):
            cx = L + ci * colw + colw * sep
            if kind in ("mag_up", "mag_dn"):
                circle(draw, (cx + extra, y), extra, fill=color, outline="white",
                       width=pt2px(1.2))
                tx = cx + 2 * extra + pt2px(4)
            elif kind == "non":
                circle(draw, (cx + extra, y), extra, fill=color)
                tx = cx + 2 * extra + pt2px(4)
            else:
                ls, lw = extra
                draw_seg(draw, (cx, y), (cx + pt2px(26), y), color, ls, 2.0)
                tx = cx + pt2px(30)
            draw.text((tx, y), label, font=f_leg, fill=(51, 51, 51), anchor="lm")


def draw_titles(draw, Wd, Hd, title, foot, title_pt=15, foot_pt=8.5,
                title_y=0.018, foot_y=0.985):
    """总标题 (顶部) + 脚注 (底部, 可多行, 居中)。"""
    draw.text((Wd / 2.0, title_y * Hd), title, font=font(pt2px(title_pt)),
              fill="black", anchor="ma")
    draw.multiline_text((Wd / 2.0, foot_y * Hd), foot, font=font(pt2px(foot_pt)),
                        fill=(115, 115, 115), anchor="md", align="center",
                        spacing=pt2px(2))


def legend_items(element, nonmag, resolved, shell_styles, aliases_updn=True,
                 spec_alt=None):
    """生成三个渲染器通用的图例条目 (自旋上/下 + 非磁 + 各壳)。"""
    up_lab = f"{element} 自旋 ↑" if aliases_updn else f"{element} ↑"
    dn_lab = f"{element} 自旋 ↓" if aliases_updn else f"{element} ↓"
    items = [("mag_up", "#d62728", up_lab, pt2px(6)),
             ("mag_dn", "#1f77b4", dn_lab, pt2px(6))]
    if nonmag:
        items.append(("non", "#c8c8c8", f"{'/'.join(nonmag)} (非磁)", pt2px(4.5)))
    for k in range(len(resolved)):
        _name, color, ls, lw, alias = shell_styles[k % len(shell_styles)]
        lab = f"J{k + 1} ({alias}) {resolved[k]:.2f} Å"
        if spec_alt and k < len(spec_alt) and spec_alt[k]:
            lab += f" [规格层 {spec_alt[k]:.2f} Å]"
        items.append(("bond", color, lab, (ls, lw)))
    return items


# ---------------------------------------------------------------------------
# 通用直角坐标系原语 (折线图 / 散点图)
#
# panel_grid + panel_transform 那套是给"面内 xy、等比缩放、不显示刻度"的俯视图用的,
# 拿它画 T 依赖曲线会得到错误的坐标轴 (T 与 U4 量纲不同, 必须独立标度 + 刻度)。
# mc_fss.py 的 Binder U4 标度图与 build_fgt_cases.py 的 FGT 磁序散点图都需要普通
# matplotlib 风格坐标系, 于是这里补一层极简 chart 原语, 两个渲染器共用。
# ---------------------------------------------------------------------------

def nice_ticks(lo, hi, n=5):
    """给出 [lo, hi] 内约 n 个"好看"的刻度值 (1/2/2.5/5 阶梯, 等价 MaxNLocator 粗版)。"""
    lo, hi = float(lo), float(hi)
    if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
        return [lo]
    raw = (hi - lo) / max(1, int(n))
    mag = 10.0 ** np.floor(np.log10(raw))
    step = 10.0 * mag
    for m in (1.0, 2.0, 2.5, 5.0, 10.0):
        if m * mag >= raw:
            step = m * mag
            break
    start = np.ceil(lo / step - 1e-9) * step
    out = []
    v = start
    while v <= hi + 1e-9 * max(1.0, abs(hi)):
        out.append(float(round(v, 10)) or 0.0)   # 归一化 -0.0 -> 0.0
        v += step
    return out or [lo]


def rotated_text(img, xy, text, f, fill=(0, 0, 0), angle=90):
    """把文字旋转 angle 度后以 xy 为中心贴到 img 上 (等价 mpl rotation=90 的 y 轴标签)。"""
    bbox = f.getbbox(text)
    pad = 4
    tmp = Image.new("RGBA", (bbox[2] - bbox[0] + 2 * pad,
                             bbox[3] - bbox[1] + 2 * pad), (0, 0, 0, 0))
    ImageDraw.Draw(tmp).text((pad - bbox[0], pad - bbox[1]), text, font=f, fill=fill)
    rot = tmp.rotate(angle, expand=True)
    img.paste(rot, (int(round(xy[0] - rot.size[0] / 2.0)),
                    int(round(xy[1] - rot.size[1] / 2.0))), rot)


def marker_shape(draw, center, r, shape, fill=None, edge=None, width=1):
    """按 mpl 记号风格画点: "o" 圆 / "s" 方 / "^" 上三角 (edge 用两段式描边)。"""
    x, y = center
    if shape in ("^", "v"):
        s = 1.0 if shape == "^" else -1.0
        pts = [(x, y - s * r * 1.15), (x + r, y + s * r * 0.75), (x - r, y + s * r * 0.75)]
        if edge is not None and width and width > 0:
            draw.polygon(pts, fill=edge)
            rr = max(1.0, r - float(width))
            pts2 = [(x, y - s * rr * 1.15), (x + rr, y + s * rr * 0.75),
                    (x - rr, y + s * rr * 0.75)]
            draw.polygon(pts2, fill=fill)
        else:
            draw.polygon(pts, fill=fill)
        return
    if shape == "s":
        if edge is not None and width and width > 0:
            draw.rectangle([x - r, y - r, x + r, y + r], fill=edge)
            r2 = max(1.0, r - float(width))
            draw.rectangle([x - r2, y - r2, x + r2, y + r2], fill=fill)
        else:
            draw.rectangle([x - r, y - r, x + r, y + r], fill=fill)
        return
    circle(draw, center, r, fill=fill, outline=edge, width=width)


def plot_series(draw, T, xs, ys, color, ls="-", lw=1.5, marker=None, ms_pt=4.0,
                marker_face=None, marker_edge=None):
    """折线 (+可选数据点标记); T 是 chart_axes 返回的 数据->像素 变换。"""
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    good = np.isfinite(xs) & np.isfinite(ys)
    xs, ys = xs[good], ys[good]
    for i in range(len(xs) - 1):
        draw_seg(draw, T((xs[i], ys[i])), T((xs[i + 1], ys[i + 1])), color, ls, lw)
    if marker:
        r = pt2px(ms_pt / 2.0)
        for i in range(len(xs)):
            marker_shape(draw, T((xs[i], ys[i])), r, marker,
                         fill=marker_face or color, edge=marker_edge,
                         width=pt2px(0.9) if marker_edge else 1)


def chart_axes(img, box, xlim, ylim, xlabel="", ylabel="", title=None,
               nticks=(6, 5), xfmt="{:.6g}", yfmt="{:.3f}", grid=True,
               grid_color=(229, 229, 229), spine=(0, 0, 0), spine_pt=0.8,
               tick_pt=9, lab_pt=10, title_pt=12, tick_len_pt=3.5, tick_pad_pt=3.0,
               grid_pt=0.8, lab_pad_pt=4.0):
    """画一个 matplotlib 风格的直角坐标系, 返回 数据坐标->画布像素 的变换 T(p)。

    box = (x0, y0, pw, ph) 为画布像素, y0 是**上**边缘 (与 panel_grid 的输出同约定)。
    只画四条轴脊 + 向外刻度 (mpl 默认风格); 范围按 nice 刻度外扩, 数据落在刻度内。
    """
    draw = ImageDraw.Draw(img)
    x0, y0, pw, ph = box
    xt = nice_ticks(xlim[0], xlim[1], nticks[0])
    yt = nice_ticks(ylim[0], ylim[1], nticks[1])
    xmin, xmax = min(xt[0], float(xlim[0])), max(xt[-1], float(xlim[1]))
    ymin, ymax = min(yt[0], float(ylim[0])), max(yt[-1], float(ylim[1]))
    if xmax <= xmin:
        xmax = xmin + 1.0
    if ymax <= ymin:
        ymax = ymin + 1.0
    fx0, fy0, fx1, fy1 = x0, y0 + ph, x0 + pw, y0

    def T(p):
        return (x0 + (p[0] - xmin) / (xmax - xmin) * pw,
                fy0 - (p[1] - ymin) / (ymax - ymin) * ph)

    if grid:
        for v in xt:
            px = T((v, ymin))[0]
            draw.line([(px, fy0), (px, fy1)], fill=grid_color,
                      width=max(1, int(round(pt2px(grid_pt)))))
        for v in yt:
            py = T((xmin, v))[1]
            draw.line([(fx0, py), (fx1, py)], fill=grid_color,
                      width=max(1, int(round(pt2px(grid_pt)))))
    sw = max(1, int(round(pt2px(spine_pt))))
    draw.rectangle([fx0, fy1, fx1, fy0], outline=spine, width=sw)
    tl = pt2px(tick_len_pt)
    tp = pt2px(tick_pad_pt)
    f_tick, f_lab = font(pt2px(tick_pt)), font(pt2px(lab_pt))
    for v in xt:
        px = T((v, ymin))[0]
        draw.line([(px, fy0), (px, fy0 + tl)], fill=spine, width=sw)
        draw.text((px, fy0 + tl + tp), xfmt.format(v), font=f_tick,
                  fill=(38, 38, 38), anchor="mt")
    maxw = 0
    for v in yt:
        py = T((xmin, v))[1]
        draw.line([(fx0 - tl, py), (fx0, py)], fill=spine, width=sw)
        s = yfmt.format(v)
        draw.text((fx0 - tl - tp, py), s, font=f_tick, fill=(38, 38, 38), anchor="rm")
        maxw = max(maxw, f_tick.getbbox(s)[2])
    if xlabel:
        draw.text((fx0 + pw / 2.0, fy0 + tl + tp + pt2px(lab_pad_pt + 11)),
                  xlabel, font=f_lab, fill=(0, 0, 0), anchor="mt")
    if ylabel:
        rotated_text(img, (fx0 - tl - tp - maxw - pt2px(lab_pad_pt + 10), fy0 - ph / 2.0),
                     ylabel, f_lab, fill=(0, 0, 0), angle=90)
    if title:
        draw.text((fx0 + pw / 2.0, fy1 - pt2px(lab_pad_pt + 3)), title, font=font(pt2px(title_pt)),
                  fill=(0, 0, 0), anchor="mb")
    return T


def legend_size(entries, f_pt=9, pad_pt=5, seg_pt=22.0):
    """图例框的 (宽, 高) 像素 —— 调用方据此把图例摆到想要的位置。"""
    f = font(pt2px(f_pt))
    pad = pt2px(pad_pt)
    seg = pt2px(seg_pt)
    w = max(f.getbbox(lab)[2] for _c, _ls, lab in entries)
    lh = pt2px(f_pt * 1.5)
    return pad * 2 + seg + pt2px(4) + w, pad * 2 + lh * len(entries)


def legend_box(draw, p0, entries, f_pt=9, pad_pt=5, edge=(204, 204, 204), lw_pt=1.6,
               seg_pt=22.0):
    """在 p0 (左上角像素) 处画一个 mpl 风格的内嵌图例 (白底细边框)。

    entries: (color, ls, label) 序列; 返回图例框的 (宽, 高) 像素。
    """
    f = font(pt2px(f_pt))
    pad = pt2px(pad_pt)
    seg = pt2px(seg_pt)
    lh = pt2px(f_pt * 1.5)
    bw, bh = legend_size(entries, f_pt, pad_pt, seg_pt)
    draw.rectangle([p0[0], p0[1], p0[0] + bw, p0[1] + bh], fill="white", outline=edge,
                   width=max(1, int(round(pt2px(0.8)))))
    for i, (color, ls, lab) in enumerate(entries):
        y = p0[1] + pad + lh * (i + 0.5)
        draw_seg(draw, (p0[0] + pad, y), (p0[0] + pad + seg, y), color, ls, lw_pt)
        draw.text((p0[0] + pad + seg + pt2px(4), y), lab, font=f, fill=(0, 0, 0),
                  anchor="lm")
    return bw, bh
