#!/usr/bin/env python
"""mc_fss.py
多尺寸有限尺寸标度 (Finite-Size Scaling): 用 Binder 累积量 U4(T) 的尺寸间交点自动定 T_C/T_N。

原理: Binder 累积量 U4(L,T) 在临界温度处**与系统尺寸无关** (标度不变), 所以不同 L 的
U4(T) 曲线在 T_C 附近交于一点; 交点温度随 L 增大收敛到真实 T_C。这比单一尺寸的
χ 峰 / |dm/dT| 判据更稳健 —— 那些判据的峰位随 L 有系统漂移 (有限尺寸伪影)。

做法: 对每个 --sizes L 复用 mc_2d.py 的 kernel 跑同一套 MC, 取 U4(T) 曲线; 两两插值
求交点 (限制在 U4 ∈ [0.05, 0.35] 的过渡带内, 避开有序 0.4 / 无序 0 的平坦区噪声),
交点温度均值 ± 标准差 = T_C 估计。

要求 --order ∈ {auto, staggered, total} (U4 是序参量判据, 比热 Cv 路径没有 U4)。

用法:
  python mc_fss.py model.json --sizes 12 18 24 30 --order auto \
      --tmin 5 --tmax 60 --nt 30 --therm 2500 --meas 6000 --seeds 2 --cool --out fss

输出: <out>_fss.csv (T 与各尺寸 U4 宽表) + <out>_fss.png + 控制台 T_C 估计。
      matplotlib 不可用时 <out>_fss.png 自动改由纯 Pillow 兜底渲染 (见 pil_fig.py)。
"""
import argparse
import json
import os

import numpy as np

import mc_2d


COLOR_CYCLE = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
               "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")   # matplotlib tab10
C_REF = "#808080"   # 参考线 (gray)
C_TC = "#d62728"    # T_C 竖线 (tab:red)


def _plot_fss_png(curves, sizes, meta, tc_mark, png):
    """纯 Pillow 兜底渲染 Binder U4(T) 标度图 (本机 matplotlib 缺依赖时走这里)。

    与 matplotlib 版同图同信息: 每尺寸一条圆点折线 + U4=0.4/0.0 灰点线参考 +
    T_C/T_N 红虚线 + 图内图例。画布按 pil_fig 约定 (figsize*200dpi, 2x 超采样后
    LANCZOS 缩小)——mpl 版此处用 dpi=150, 故像素尺寸更大但版面一致。
    唯一文本差异: mpl 的 LaTeX 下标 "Binder $U_4$" 写成 "Binder U4" (Pillow 无排版)。
    """
    import sys as _sys
    _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import pil_fig as P
    from PIL import Image, ImageDraw

    W, H = int(round(7 * P.DPI)), int(round(5 * P.DPI))
    Wd, Hd = W * P.SS, H * P.SS
    img = Image.new("RGB", (Wd, Hd), "white")
    dr = ImageDraw.Draw(img)

    T_grid = np.asarray(curves[sizes[0]][0], dtype=float)
    ulo = min([0.0] + [float(np.min(np.asarray(curves[L][1], dtype=float))) for L in sizes])
    uhi = max([0.4] + [float(np.max(np.asarray(curves[L][1], dtype=float))) for L in sizes])
    xr, xR = float(T_grid.min()), float(T_grid.max())
    if xR <= xr:
        xr, xR = xr - 1.0, xr + 1.0
    px, py = 0.03 * (xR - xr), 0.05 * max(uhi - ulo, 1e-6)
    box = (0.125 * Wd, 0.12 * Hd, 0.775 * Wd, 0.77 * Hd)   # mpl 默认 axes rect
    T = P.chart_axes(img, box, (xr - px, xR + px), (ulo - py, uhi + py),
                     xlabel="T (K)", ylabel="Binder U4",
                     title=f"{meta['element']} FSS (Binder U4 crossing)",
                     nticks=(7, 5), xfmt="{:.0f}", yfmt="{:.1f}")
    for yv in (0.4, 0.0):
        P.draw_seg(dr, T((xr - px, yv)), T((xR + px, yv)), C_REF, ":", 0.8)
    for i, L in enumerate(sizes):
        P.plot_series(dr, T, curves[L][0], curves[L][1], COLOR_CYCLE[i % len(COLOR_CYCLE)],
                      ls="-", lw=1.5, marker="o", ms_pt=4.0)
    if tc_mark is not None:
        P.draw_seg(dr, T((tc_mark, ulo - py)), T((tc_mark, uhi + py)), C_TC, "--", 1.0)
    entries = [(COLOR_CYCLE[i % len(COLOR_CYCLE)], "-", f"L={L}") for i, L in enumerate(sizes)]
    if tc_mark is not None:
        entries.append((C_TC, "--", f"T_C={tc_mark:.1f} K"))
    # mpl 的 loc="best": U4 曲线从左上降到右下, 右上角必空 => 固定摆轴内右上角
    bw, _bh = P.legend_size(entries)
    P.legend_box(dr, (box[0] + box[2] - bw - P.pt2px(8), box[1] + P.pt2px(8)), entries)
    img = img.resize((W, H), Image.LANCZOS)
    os.makedirs(os.path.dirname(os.path.abspath(png)) or ".", exist_ok=True)
    img.save(png, "PNG", optimize=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("model")
    ap.add_argument("--sizes", type=int, nargs="+", required=True,
                    help="至少两个线性尺寸 L (方晶格 LxL), 如 12 18 24 30")
    ap.add_argument("--order", default="auto", help=mc_2d.ORDER_HELP)
    ap.add_argument("--tmin", type=float, default=5.0)
    ap.add_argument("--tmax", type=float, default=150.0)
    ap.add_argument("--nt", type=int, default=30)
    ap.add_argument("--therm", type=int, default=2500)
    ap.add_argument("--meas", type=int, default=6000)
    ap.add_argument("--cool", action="store_true",
                    help="从高温往低温退火扫描 (低温启动时避免被亚稳态冻住)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--seeds", type=int, default=1,
                    help="每个尺寸独立种子数 (平均降噪, 交点更稳)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if len(args.sizes) < 2:
        raise SystemExit("--sizes 至少给 2 个线性尺寸 (如 --sizes 12 18)")

    with open(args.model, encoding="utf-8") as f:
        model = json.load(f)
    meta = model["meta"]

    # 序参量选择: 同 mc_2d.py 的 auto 逻辑 (J 求和符号定铁磁/反铁磁)
    try:
        o_kind, o_mask = mc_2d.parse_order_spec(args.order)
    except ValueError as e:
        raise SystemExit(f"[error] {e}")
    if o_kind not in ("auto", "staggered", "total", "cv", "sublattice"):
        raise SystemExit(f"[error] 未知 --order {args.order!r}: 可用 {mc_2d.ORDER_HELP}")
    if o_kind == "cv":
        raise SystemExit("[error] FSS 用 U4 交点定临界温度, 没有 U4 的比热路径不适用 "
                         "(与 mc_2d.py 的 --order cv 不同)")
    if o_kind == "auto":
        if "J1_meV" in meta:
            jsum = float(meta["J1_meV"])
        else:
            jsum = sum(-float(it["J"][0][0]) for it in model.get("interactions", []))
        order = 1 if jsum < 0 else 0
        kind = "交错 (反铁磁, T_N)" if order == 1 else "总磁化 (铁磁, T_C)"
    elif o_kind == "sublattice":
        order = 1
        kind = f"逐子晶格 Néel 矢量 |L| (A = 子晶格 {sorted(set(o_mask))}, B = 其余)"
    else:
        order = 1 if o_kind == "staggered" else 0
        kind = "交错 (反铁磁, T_N)" if order == 1 else "总磁化 (铁磁, T_C)"
    print(f"序参量: {kind}")

    prefix = args.out or (os.path.splitext(os.path.basename(args.model))[0] + "_fss")
    T_list = np.linspace(args.tmin, args.tmax, args.nt)
    if args.cool:
        T_list = T_list[::-1]
    n_seeds = max(1, int(args.seeds))

    curves = {}  # L -> (T, U4)
    for L in args.sizes:
        # neighbors/J/onsite/num_basis 与种子无关, 先建一次; 各种子只重建初始自旋。
        # 交错序参量要子晶格 ±1 符号: num_basis==1 的方形/矩形必须图二着色 ——
        # 老代码按 `i % num_basis` 分会退化成单侧为空 (cb==0) 而除零。
        _s, neighbors, J, onsite, num_basis = mc_2d.build_arrays(model, (L, L), spin_seed=0)
        stag = np.zeros(_s.shape[0], dtype=np.float64)
        if order == 1:
            if o_kind == "sublattice":
                try:
                    stag = mc_2d.sublattice_signs(num_basis, (L, L), o_mask)
                except ValueError as e:
                    raise SystemExit(f"[error] --order {args.order}: {e}")
            else:
                stag = mc_2d.stag_signs(neighbors, J, num_basis)
                if stag is None:
                    raise SystemExit(
                        f"序参量=交错 需要二分子晶格, 但该模型 (num_basis={num_basis}) 无法二分着色"
                        f" (受挫/非二分晶格) —— 请改用 `--order total` 用 U4 交点定 T_C")
        all_res = []
        for k in range(n_seeds):
            spins, _, _, _, _ = mc_2d.build_arrays(model, (L, L), spin_seed=args.seed + k)
            res = mc_2d.run_sweep(spins, neighbors, J, onsite, T_list, args.therm,
                                  args.meas, stag, order, seed=args.seed + k)
            all_res.append(res[np.argsort(res[:, 0])])
        res = np.mean(np.stack(all_res), axis=0)
        curves[L] = (res[:, 0], res[:, 3])  # T, U4
        print(f"L={L:3d}: U4 [{res[:, 3].min():.3f}, {res[:, 3].max():.3f}] "
              f"(平均 {n_seeds} 种子)", flush=True)

    sizes = sorted(curves)
    T_grid = curves[sizes[0]][0]

    # 平滑窗口: 压掉单种子的 Binder 采样噪声 (Binder 累积量在高 T 无序相会因有限样本
    # 在 0 附近抖动、甚至掉到负值; 窗口越大越平滑, 但别吃掉真峰)。≥7 个温度点才平滑。
    sw = 3 if len(T_grid) >= 7 else 1

    def smooth(y):
        if sw < 2 or len(y) < sw:
            return y
        k = np.ones(sw) / sw
        return np.convolve(y, k, mode="same")

    # U4 曲线矩阵: 所有尺寸共享同一温度网格 (T_grid), 可直接按列做统计。
    Umat = np.stack([np.interp(T_grid, curves[L][0], curves[L][1]) for L in sizes], axis=0)
    mean_U = Umat.mean(axis=0)

    # 主判据 A (collapse): 在标度不变点 T_C 处 U4 与尺寸无关 -> 尺寸间方差最小。
    # 过渡带 (平均 U4 ∈ [0.05, 0.35]) 内找方差最小; 平坦区 (全 0.4 或全 0) 方差也小,
    # 但不在过渡带内, 被排除。这比两两交点更稳: 一次性用上所有尺寸, 不受单对噪声影响。
    band = (mean_U >= 0.05) & (mean_U <= 0.35)
    T_collapse = None
    if band.any():
        vs = smooth(Umat.var(axis=0))
        T_collapse = float(T_grid[band][int(np.argmin(vs[band]))])

    # 主判据 B (两两交点): 平滑后, 在 U4 过渡带 [0.02, 0.38] 内取**最陡**的那个符号
    # 变化 —— 真相变点在 U4 急降处 (斜率大), 平坦区 (有序 0.4 / 无序 0) 的噪声小抖动
    # 被过渡带过滤掉, 残余抖动也因斜率小而不被选中。每对尺寸只取 1 个交点。
    crossings = []
    for i in range(len(sizes)):
        for j in range(i + 1, len(sizes)):
            T1, U1 = curves[sizes[i]]
            T2, U2 = curves[sizes[j]]
            U1s = smooth(U1)
            U2i = smooth(np.interp(T1, T2, U2))
            diff = U1s - U2i
            cand = []
            for k in range(len(T1) - 1):
                d0, d1 = diff[k], diff[k + 1]
                if d0 * d1 < 0.0:
                    tc = T1[k] - d0 * (T1[k + 1] - T1[k]) / (d1 - d0)
                    uc = np.interp(tc, T1, U1s)
                    if 0.02 <= uc <= 0.38:
                        slope = abs(d1 - d0) / (T1[k + 1] - T1[k])
                        cand.append((tc, uc, slope))
            if cand:
                cand.sort(key=lambda x: -x[2])  # 最陡优先
                crossings.append((sizes[i], sizes[j], cand[0][0]))

    header = "T_Kelvin," + ",".join(f"U4_L{L}" for L in sizes)
    cols = [T_grid] + [Umat[k] for k in range(len(sizes))]
    out_csv = prefix + "_fss.csv"
    np.savetxt(out_csv, np.column_stack(cols), delimiter=",", header=header, comments="")
    print(f"\n>>> 数据存至 {out_csv}")

    # 综合两种判据: 优先两两交点均值 (用户要的 U4 交点), collapse 作为交叉验证。
    tc = sd = None
    if crossings:
        tc_vals = np.array([c[2] for c in crossings])
        tc = float(tc_vals.mean())
        sd = float(tc_vals.std(ddof=1)) if len(tc_vals) > 1 else 0.0
        print(f">>> U4 交点 ({len(crossings)} 个尺寸对):")
        for a, b, t in crossings:
            print(f"    L={a} x L={b} -> {t:.2f} K")
        if T_collapse is not None:
            print(f">>> U4 方差最小 (collapse): {T_collapse:.2f} K")
        print(f">>> 有限尺寸标度 T_C/T_N ≈ {tc:.2f} K ± {sd:.2f} K "
              f"({len(crossings)} 个尺寸对交点标准差)")
    elif T_collapse is not None:
        tc, sd = T_collapse, None
        print(f">>> 两两交点未找到干净的, 但 U4 方差最小 (collapse) = {T_collapse:.2f} K; "
              f"以此作为 T_C/T_N 估计 (建议加 --seeds / 加密温度点以得到两两交点)")
    else:
        print(">>> [warn] 过渡带内既无 U4 交点也无 collapse 点: "
              "窗口可能没覆盖相变, 或尺寸差太小/噪声太大; "
              "请调整 --tmin/--tmax/--nt, 增大 --sizes 跨度, 或加 --seeds 降噪")

    # 机器可读结果 (driver --fss 读它; 独立跑时也方便下游解析)。始终写出, tc 可能为 None。
    with open(prefix + "_fss.json", "w", encoding="utf-8") as f:
        json.dump({"method": "binder_u4_fss", "sizes": sizes, "tc": tc, "sd": sd,
                   "tc_collapse": T_collapse, "n_crossings": len(crossings),
                   "crossings": [[int(a), int(b), float(t)] for a, b, t in crossings]},
                  f, ensure_ascii=False, indent=2)

    # 图上标出的转变温度口径与上面打印的一致 (两两交点优先, 否则 collapse)。
    # 必须在 try 之前算出来: 否则 matplotlib 一失败, 兜底渲染器就拿不到它
    # (mc_2d.py 的 cv_peak 踩过同一个坑)。
    tc_mark = tc if crossings else T_collapse
    png = prefix + "_fss.png"
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 5))
        for L in sizes:
            ax.plot(curves[L][0], curves[L][1], "o-", ms=4, label=f"L={L}")
        ax.axhline(0.4, color="gray", ls=":", lw=0.8)
        ax.axhline(0.0, color="gray", ls=":", lw=0.8)
        if tc_mark is not None:
            ax.axvline(tc_mark, color="tab:red", ls="--", lw=1,
                       label=f"$T_C$={tc_mark:.1f} K")
        ax.set_xlabel("T (K)")
        ax.set_ylabel("Binder $U_4$")
        ax.set_title(f"{meta['element']} FSS (Binder U4 crossing)")
        ax.legend()
        ax.grid(alpha=0.3)
        fig.savefig(png, dpi=150, bbox_inches="tight")
        print(f">>> 图存至 {png}")
    except Exception as e:
        print(f"[warn] matplotlib 画图不可用 ({e}); 改用纯 Pillow 兜底")
        try:
            _plot_fss_png(curves, sizes, meta, tc_mark, png)
            print(f">>> 图存至 {png} (纯 Pillow 兜底)")
        except Exception as e2:
            print(f"[warn] 画图跳过: {e2}")


if __name__ == "__main__":
    main()
