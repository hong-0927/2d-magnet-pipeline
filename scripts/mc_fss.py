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
import sys

import numpy as np

import mc_2d


COLOR_CYCLE = ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
               "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")   # matplotlib tab10
C_REF = "#808080"   # 参考线 (gray)
C_TC = "#d62728"    # T_C 竖线 (tab:red)


def _pick_crossing(cand, ref):
    """从一对尺寸的候选交点里选一个: 取离 `ref` (全局 collapse 估计) 最近的。

    旧实现取「斜率最陡」的候选 (P1-1): 平滑后残余的噪声抖动只要有局部陡坡就可能被选中,
    这是选择偏差, 而且会让 T_C 系统性偏向斜率大的位置。collapse 估计是尺寸无关的全局判据
    (一次用上所有尺寸、不受单对噪声影响), 用它当参照更稳。
    """
    return sorted(cand, key=lambda x: abs(x[0] - ref))[0]


def _extrapolate_1overL(crossings, tc_vals):
    """标度区假设下按 T_c(L) = T_c(∞) + a/L 做最小二乘, 返回 (T_c(∞), a)。

    x 取每个尺寸对的平均 L 的倒数; 只有 >=3 个交点才谈得上拟合。斜率 a 的大小本身是
    「有没有进标度区」的体检指标 —— |a| 大说明交点随 L 还在明显漂移, 外推不可信。
    """
    xs = np.array([1.0 / (0.5 * (a + b)) for a, b, _ in crossings], dtype=float)
    A = np.column_stack([np.ones_like(xs), xs])
    coef = np.linalg.lstsq(A, np.asarray(tc_vals, dtype=float), rcond=None)[0]
    return float(coef[0]), float(coef[1])


def _selftest():
    """本脚本自己的不变量: 交点选择不取「最陡」+ 1/L 外推能还原真值。

    与 smoke §5f 同源 (那里跑同一组断言, 这里给手动 `--selftest` 用)。
    """
    fails = []

    def ck(name, cond):
        print(("  [PASS] " if cond else "  [FAIL] ") + name)
        if not cond:
            fails.append(name)

    # 1. 交点选择: 候选 (8.0, 斜率 1) 与 (9.0, 斜率 100), 参照 8.1 -> 必须选 8.0
    pick = _pick_crossing([(8.0, 0.30, 1.0), (9.0, 0.30, 100.0)], 8.1)[0]
    ck(f"交点选择取近参照 (期望 8.0, 实得 {pick})", abs(pick - 8.0) < 1e-9)
    # 2. 1/L 外推: 合成「严格在 T = 10 + 40/L 上」的交点 -> 截距 10, 斜率 40
    cr = [(4, 4, 20.0), (4, 8, 16.6666666667), (8, 8, 15.0)]
    t_inf, drift = _extrapolate_1overL(cr, [c[2] for c in cr])
    ck(f"1/L 外推截距 (期望 10.0, 实得 {t_inf:.4f})", abs(t_inf - 10.0) < 1e-6)
    ck(f"1/L 外推斜率 (期望 40.0, 实得 {drift:.4f})", abs(drift - 40.0) < 1e-6)
    # 3. 报告措辞: ± 必须是跨尺寸离散度, 不能写成统计标准误
    #    被禁的短语在本文件里拼出来, 否则自检会读到自己的源码而永远失败
    src = open(os.path.abspath(__file__), encoding="utf-8").read()
    bad = "交点" + "标准差"
    ck("± 的措辞是「跨尺寸对离散」而不是「统计标准误」",
       "不是统计标准误" in src and bad not in src)
    if fails:
        print(f"  [FAIL] {len(fails)} 项失败: {fails}")
        return 1
    print("  ALL PASS")
    return 0


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
    if "--selftest" in sys.argv[1:]:
        return _selftest()
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
        # 判据与 mc_2d.py 同源 (P1-2): 比 T=0 候选基态能量, 不看 J1 符号。
        _L = max(3, min(int(x) for x in args.sizes))
        order, kind, use_cv, cand, best = mc_2d.resolve_auto_order(model, (_L, _L))
        ranked = sorted(cand.items(), key=lambda kv: kv[1])
        print(f"auto 序参量判定 (按 L={_L} 的 T=0 候选基态能量): "
              + " | ".join(f"{k}={v:+.4f}" for k, v in ranked))
        if use_cv:
            raise SystemExit(
                f"[error] auto 判定该模型 T=0 基态是 {best} (条纹反铁磁): ±1 交错序参量测的"
                f"不是这个序, 而 FSS 用 U4 交点定临界温度、没有比热路径 —— 请显式给 "
                f"--order sublattice:<mask> (按条纹方向分子晶格), 或先用 "
                f"`mc_2d.py --order cv` 看比热峰")
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
                # 一对尺寸在过渡带里可能出现多个交点 (平滑后的残余抖动)。旧实现取「最陡」,
                # 那是**选择偏差**: 噪声小抖动只要有局部陡坡就可能被选中, 而且系统性地偏向
                # 斜率大的位置 (P1-1)。改用尺寸无关的全局判据当参照: collapse 估计 (它一次
                # 用上所有尺寸、不依赖单对噪声) 离得最近的那个交点。
                ref = T_collapse if T_collapse is not None else float(np.median([c[0] for c in cand]))
                crossings.append((sizes[i], sizes[j], _pick_crossing(cand, ref)[0]))

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
              f"({len(crossings)} 个尺寸对交点的离散度, **不是统计标准误**: 这些交点是"
              f"相关的估计, 不是独立样本; 随 L 单调漂移时应当读作系统误差/标度区未达)")
        if len(crossings) >= 3:
            # 1/L 外推: 标度区里 T_c(L) = T_c(∞) + a/L, 用各尺寸对的平均 L 作 x。
            t_inf, drift = _extrapolate_1overL(crossings, tc_vals)
            print(f">>> 1/L 线性外推 (参考值, 前提是已进标度区): "
                  f"T_C/T_N(L->∞) ≈ {t_inf:.2f} K (斜率 {drift:+.1f} K·L; "
                  f"|斜率| 大 = 尺寸还没进标度区, 此时该外推不可信)")
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
    sys.exit(main())
