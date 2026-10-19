#!/usr/bin/env python3
"""交变磁体 (altermagnet) 识别层检查清单 —— P5 第 25 项 (knowledge/NEW_TYPES_PLAN.md)。

**只出报告, 不下结论。** 本脚本判定的是**几何清单**:

 1. 补偿: 磁性原子分成**两套亚晶格**(按分数坐标的等价类), 原子数相等, 自旋反向 ⇒ 净磁化严格为零;
 2. 连接操作必须是**晶体的对称操作** —— 即把**全部**原子按元素映射回原子位置(容差 --tol);
    这一条是判据的关键: 「两套亚晶格之间有个平移」这种说法只有在那个平移**同时**是晶体对称操作时才算数,
    否则它只是一个坐标差(例: 金红石结构里 (1/2,1/2,1/2) 把两个磁性亚晶格互相映射, 但它不是晶体对称操作,
    因为非磁原子的位置对不上 ⇒ 真正的连接操作是 4_2 螺旋, 旋转部分非平凡 ⇒ 交变磁体)。
 3. **旋转部分**决定分类:
      M = I 且平移非格矢  → 平移连接 → 常规共线反铁磁 (自旋劈裂被对称性禁止)
      M = -I 且非共面胞   → 反演连接 → 宇称禁止 (也不是交变磁体)
                             (共面 2D 胞例外: 那里的「反演」与面内 C2 是同一个操作,
                              真二维交变磁体天然落在「旋转连接」这一支)
      M ≠ ±I               → 旋转连接 → **交变磁体几何候选** (真旋转 / 镜像 / 旋转反演都算)
    反演只要**在场**就判宇称禁止 —— 同一 t 下同时存在反演与旋转时, 不许因为「还找到了旋转」改判候选。

为什么这条分类不是拍脑袋 (「连接操作 g + 自旋翻转」是磁性结构的对称操作 gT, 自旋 1/2 的
双群表示决定 (gT)² 的符号, 也就是有没有 Kramers 型的强制简并):
      g = 平移   → (gT)² = -1 ⇒ 所有 k 强制二重简并 ⇒ 自旋劈裂被禁止 (常规反铁磁)
      g = 反演   → (gT)² = -1 ⇒ 同样强制简并   ⇒ 禁止 (宇称禁止)
      g = 旋转/镜像/旋转反演 → 代数的 ² 给出 +1 ⇒ 没有强制简并 ⇒ 劈裂允许 (交变磁体的几何前提)
    所以「镜像连接也算候选」是有代价的: 一个**常规**蜂窝 Néel 反铁磁只要两套亚晶格恰好由镜面
    相连, 也会落进「候选」这一支 (实测 MnPS3 的 fixture 就是这样, 残差 0.012 mÅ)。「候选」
    只表示**没有被几何判据排除**, 不等于「是交变磁体」—— 请照 `说明` 与 `边界` 两行读。

**容差按 Å 判, 不按「四舍六入查元组」**: 判据的唯一容差来源是映射后原子与同元素最近原子的
笛卡尔距离 (`--tol`, 默认 1e-3 Å)。真实结构文件里同一条对称关系两侧坐标常在第 6~7 位小数上
不同 (实测 MnPS3/CrOCl fixture 的残差只有 0.01~0.5 mÅ, 按 6 位四舍五入查表会全部落空 —— 那是
把「精确对称」误报成「没有对称」的静默错误)。无连接时报告仍给出**最接近的候选操作与残差**,
让人自己判断是「真的没有」还是「弛豫/写文件精度差一点」。

**不可判定**: k 分辨自旋劈裂(能带层面) —— 那需要自旋分辨能带计算, 属研究问题, 不在本流水线契约内
(见 knowledge/LIMITATIONS.md 10 / knowledge/WORKFLOW.md 局限 1g)。本脚本永远不会输出「这是交变磁体」这种结论。

用法:
    python scripts/altermagnet_check.py --poscar tests/fixtures/Fe3Sn_POSCAR --element Fe
    python scripts/altermagnet_check.py --poscar X_POSCAR --element Mn --spins +1,-1 --json out.json
    python scripts/altermagnet_check.py --selftest
"""
import argparse
import itertools
import json
import os
import sys

import numpy as np

# 报告文本含 '⇒'/'²' 等非 GBK 字符: cp936 控制台直跑会 UnicodeEncodeError。
# 与 run_2d_magnet.py 顶部同款处理: 只放宽 errors, 不改编码。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors="replace")
    except Exception:
        pass

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import magnetic_orders as MO  # noqa: E402  (parse_poscar / _sublattice_keys, 与规格层同一套约定)

VERDICT_CONVENTIONAL = "常规共线反铁磁（两套亚晶格由平移连接 ⇒ 自旋劈裂被对称性禁止）"
VERDICT_PARITY = "反演连接（宇称禁止自旋劈裂 ⇒ 也不是交变磁体）"
VERDICT_CANDIDATE = "交变磁体几何候选（两套亚晶格由旋转连接 ⇒ 满足几何检查清单）"
VERDICT_NONE = "不成立：没有任何晶体对称操作连接两套亚晶格"
VERDICT_NA = "不适用：磁性原子不是「两套原子数相等的亚晶格」"
DISCLAIMER = ("不可判定：k 分辨自旋劈裂（能带层面）—— 本清单只做几何分类, "
              "不构成「这是交变磁体」的判定结论。")


def metric(lat):
    """分数坐标的度规 G = A A^T (行向量约定: r = f @ A, |r|^2 = f G f^T)。"""
    return np.asarray(lat, dtype=float) @ np.asarray(lat, dtype=float).T


def point_ops(lat, tol=1e-6):
    """晶格的**点群**在分数坐标下的整数矩阵 M (行约定 f' = f @ M)。

    枚举 {-1,0,1} 的 3x3 整数矩阵并要求 M G M^T = G (保长) 且 det = ±1。
    常规三维晶格(含 120° 六角胞)的点群都能被覆盖。
    """
    G = metric(lat)
    out = []
    for vals in itertools.product((-1, 0, 1), repeat=9):
        M = np.array(vals, dtype=float).reshape(3, 3)
        if abs(np.linalg.det(M)) < 0.5:
            continue
        if np.max(np.abs(M @ G @ M.T - G)) <= tol * max(1.0, np.max(np.abs(G))):
            out.append(M.astype(int))
    return out


def rot_kind(M):
    """按旋转部分给操作分类。"""
    I = np.eye(3, dtype=int)
    if np.array_equal(M, I):
        return "translation"
    if np.array_equal(M, -I):
        return "inversion"
    if int(round(np.linalg.det(M))) == 1:
        return "proper"
    return "improper"


def cart_rotation(lat, M):
    """分数坐标算子 -> 笛卡尔行向量算子 Rc (r' = r @ Rc)。"""
    A = np.asarray(lat, dtype=float)
    return np.linalg.inv(A) @ np.asarray(M, dtype=float) @ A


def rot_axis_angle(lat, M, tol=1e-8):
    """真旋转的轴(笛卡尔, 单位化)与角度(度); 非真旋转返回 (None, None)。"""
    if rot_kind(M) != "proper":
        return None, None
    Rc = cart_rotation(lat, M)
    ang = float(np.degrees(np.arccos(np.clip((np.trace(Rc) - 1.0) / 2.0, -1.0, 1.0))))
    w, v = np.linalg.eig(Rc.T)
    axis = None
    for k in range(len(w)):
        if abs(w[k] - 1.0) < 1e-6:
            cand = np.real(v[:, k])
            n = np.linalg.norm(cand)
            if n > tol:
                axis = cand / n
                break
    return axis, ang


def _frac_to_cart(lat, d):
    """分数坐标差 -> 笛卡尔位移 (行向量约定: r = f @ A)。"""
    return np.asarray(d, dtype=float) @ np.asarray(lat, dtype=float)


def nearest_same_element(lat, frac, atoms, a, g):
    """把映射后的分数坐标 g 落到「同元素最近原子」的笛卡尔距离 (最小镜像, 单位 Å)。

    这是判据的**唯一**容差来源: 不用「四舍六入后查元组」——真实结构文件里同一条对称
    关系两侧的坐标常在**第 6~7 位小数**上不同 (MP/DFT 弛豫/写文件精度), 查元组会把
    精确的对称操作判成不存在 (实测: MnPS3/CrOCl fixture 的残差只有 0.01~0.5 mÅ,
    但按 6 位四舍五入查表全部落空)。
    """
    best = None
    for a2, f2 in zip(atoms, frac):
        if a2 != a:
            continue
        d = g - f2
        d = d - np.round(d)
        r = float(np.linalg.norm(_frac_to_cart(lat, d)))
        if best is None or r < best:
            best = r
    return float("inf") if best is None else best


def sublattices(atoms, frac, element, dec=6):
    """磁性原子的等价类(与规格层 _sublattice_keys 同约定: 四舍六入后字典序)。

    返回 (mag_idx, sub_of_atom, n_sub, sub_sizes)。
    """
    idxs, key_map, _keys = MO._sublattice_keys(atoms, frac, element)
    sub_of = {i: key_map[tuple(np.round(frac[i], 6))] for i in idxs}
    sizes = [sum(1 for i in idxs if sub_of[i] == k) for k in range(len(key_map))]
    return idxs, sub_of, len(key_map), sizes


def check(lat, atoms, frac, element, spins=None, tol=1e-3, max_pairs=4000):
    """跑一遍检查清单, 返回结果 dict (不抛异常, 便于测试与 JSON 输出)。

    tol: 判定「(M, t) 是晶体对称操作」的容差, **单位 Å** (映射后原子与同元素最近原子
         的笛卡尔距离)。理想结构用默认 1e-3 Å; DFT 弛豫/实验结构可放宽到 0.01~0.05 Å。
         报告里同时给出最接近的候选操作及其残差, 便于判断「差一点」是不是精度问题。
    """
    idxs, sub_of, n_sub, sizes = sublattices(atoms, frac, element)
    res = {
        "element": element, "n_sub": int(n_sub), "sub_sizes": [int(s) for s in sizes],
        "n_mag": int(len(idxs)), "compensated": False, "ops": [],
        "verdict": VERDICT_NA, "note": "", "disclaimer": DISCLAIMER,
    }
    if n_sub != 2:
        res["note"] = f"磁性原子 {element} 有 {n_sub} 套分数坐标等价类（需要恰好 2 套）"
        return res
    if sizes[0] != sizes[1]:
        res["note"] = f"两套亚晶格原子数不等: {sizes}"
        return res
    if spins is None:
        spins = [1, -1]
    spins = [int(s) for s in spins]
    if len(spins) != 2:
        res["note"] = "spins 必须给 2 个值"
        return res
    res["spins"] = spins
    res["compensated"] = bool(sum(spins) == 0 and sizes[0] == sizes[1])

    A_idx = [i for i in idxs if sub_of[i] == 0]
    B_idx = [i for i in idxs if sub_of[i] == 1]
    # 共面性: 全部原子 z 相同 ⇒ 这是二维单层胞, 此时「反演」与「面内 C2」是同一个操作,
    # 本清单永远落不到「反演连接」分支(二维交变磁体因此天然是「旋转连接」)。
    res["coplanar"] = bool(np.max(np.abs(frac[:, 2] - frac[0, 2])) < 1e-3)
    res["tol_ang"] = float(tol)

    found = []
    best = None            # (残差 Å, M, kind, t) —— 即使不成立也记下最接近的一条
    pairs = 0

    def _note_candidate(M, kind, t):
        nonlocal best
        ok, worst = _is_symmetry(lat, M, t, atoms, frac, tol)
        if best is None or worst < best[0]:
            best = (worst, M.tolist(), kind, [float(x) for x in np.round(t, 6)])
        return ok

    for M in point_ops(lat):
        kind = rot_kind(M)
        if kind == "translation":
            # 平移: 只需要看坐标差(A 的某点 + t = B 的某点), 与 A 里哪个点无关
            for i in A_idx[:1]:
                for j in B_idx:
                    t = np.mod(frac[j] - frac[i] @ M, 1.0)
                    tw = t - np.round(t)
                    if float(np.linalg.norm(_frac_to_cart(lat, tw))) <= tol:
                        continue  # 格矢平移(恒等): 会把 A 映到 A, 不是连接操作
                    if _note_candidate(M, kind, t):
                        found.append({"M": M.tolist(), "kind": "translation",
                                      "t": [float(x) for x in np.round(t, 6)]})
                    break
            continue
        for i in A_idx:
            for j in B_idx:
                pairs += 1
                if pairs > max_pairs:
                    break
                t = np.mod(frac[j] - frac[i] @ M, 1.0)
                if _note_candidate(M, kind, t):
                    axis, ang = rot_axis_angle(lat, M)
                    found.append({"M": M.tolist(), "kind": kind,
                                  "t": [float(x) for x in np.round(t, 6)],
                                  "angle_deg": (None if ang is None else round(ang, 4)),
                                  "axis": (None if axis is None else [round(float(x), 6) for x in axis])})

    # 去重(同一个 (M,t) 可能被多个 (i,j) 命中)
    uniq = []
    seen = set()
    for o in found:
        key = (tuple(o["M"][0]), tuple(o["M"][1]), tuple(o["M"][2]), tuple(o["t"]))
        if key not in seen:
            seen.add(key)
            uniq.append(o)
    res["ops"] = uniq
    if best is not None:
        res["best_mismatch_ang"] = round(float(best[0]), 8)
        res["best_candidate"] = {"M": best[1], "kind": best[2], "t": best[3]}
    if uniq:
        # 已接受的操作里最差的那条残差 (给「这个对称性有多干净」一个数)
        _op_worst = 0.0
        for o in uniq:
            _ok, _w = _is_symmetry(lat, np.array(o["M"], dtype=float), np.array(o["t"], dtype=float),
                                   atoms, frac, tol)
            _op_worst = max(_op_worst, _w)
        res["worst_mismatch_ang"] = round(float(_op_worst), 8)

    kinds = {o["kind"] for o in uniq}
    if not uniq:
        res["verdict"] = VERDICT_NONE
        _bm = res.get("best_mismatch_ang")
        _bc = res.get("best_candidate")
        if _bm is not None and _bc is not None:
            res["note"] = (f"没有晶体对称操作在容差 {tol:g} Å 内连接两套亚晶格; "
                           f"最接近的一条是 {_bc['kind']} t={_bc['t']}, 残差 {_bm * 1000:.2f} mÅ "
                           f"(放宽 --tol 可看到它是否只是弛豫/写文件精度问题)")
        else:
            res["note"] = "两套亚晶格之间没有任何候选连接操作"
    elif "translation" in kinds:
        res["verdict"] = VERDICT_CONVENTIONAL
        res["note"] = "平移连接优先于旋转: 只要存在一条**晶体对称的**非格矢平移连接两套亚晶格, 自旋劈裂就被禁止"
    elif "inversion" in kinds and not res["coplanar"]:
        # 反演连接 ⇒ 宇称禁止自旋劈裂。共面 2D 胞是例外: 那里的「反演」就是面内 C2,
        # 属旋转连接 (真二维交变磁体天然落在这一支), 故只在非共面时判 PARITY。
        res["verdict"] = VERDICT_PARITY
        _others = sorted(kinds - {"inversion"})
        res["note"] = ("连接操作含反演 ⇒ 宇称禁止自旋劈裂"
                       + (f"; 同一 t 下另有 {_others} —— 反演在场即已禁止" if _others else ""))
    else:
        res["verdict"] = VERDICT_CANDIDATE
        res["note"] = "连接操作含非平凡旋转部分; 判据到此为止 —— 是否真有 k 分辨自旋劈裂要能带计算"
    return res


def _is_symmetry(lat, M, t, atoms, frac, tol_ang):
    """(M, t) 是否把**全部**原子按元素映射回原子位置 (容差 tol_ang, 单位 Å)。

    返回 (是否成立, 最差原子的残差 Å)。残差即使不成立也要返回 —— 报告里要写清
    「没有任何对称操作」与「最接近的操作差 0.02 Å」是两件完全不同的事。
    """
    worst = 0.0
    for a, f in zip(atoms, frac):
        g = np.mod(f @ M + t, 1.0)
        r = nearest_same_element(lat, frac, atoms, a, g)
        if r > worst:
            worst = r
        if r > tol_ang:
            return False, worst
    return True, worst


def fmt_report(res, poscar=None):
    lines = []
    if poscar:
        lines.append(f"POSCAR      : {poscar}")
    lines.append(f"磁性元素    : {res['element']}   (磁原子 {res['n_mag']} 个 / 胞)")
    lines.append(f"亚晶格      : {res['n_sub']} 套, 原子数 {res['sub_sizes']}"
                 + (f", 自旋 {res['spins']}" if "spins" in res else ""))
    lines.append(f"补偿(净磁化=0): {'是' if res['compensated'] else '否'}")
    if res.get("coplanar"):
        lines.append("共面 2D 胞  : 是（此时「反演」与「面内 C2」是同一操作 ⇒ 不会出现「反演连接」分支）")
    if res["ops"]:
        lines.append(f"连接操作    : {len(res['ops'])} 条 (下列为去重后的晶体对称操作, "
                     f"容差 {res.get('tol_ang', 0):g} Å"
                     + (f", 最差残差 {res['worst_mismatch_ang'] * 1000:.3f} mÅ"
                        if "worst_mismatch_ang" in res else "") + ")")
        for o in res["ops"][:8]:
            extra = ""
            if o.get("angle_deg") is not None:
                extra = f"  轴={o['axis']} 角={o['angle_deg']}°"
            elif o["kind"] == "inversion":
                extra = "  (M = -I)"
            lines.append(f"  - {o['kind']:<12} t={o['t']}{extra}")
        if len(res["ops"]) > 8:
            lines.append(f"  … 其余 {len(res['ops']) - 8} 条见 --json")
    else:
        lines.append(f"连接操作    : 无 (容差 {res.get('tol_ang', 0):g} Å)")
        if "best_mismatch_ang" in res and "best_candidate" in res:
            bc = res["best_candidate"]
            lines.append(f"  最接近者  : {bc['kind']} t={bc['t']}, 残差 "
                         f"{res['best_mismatch_ang'] * 1000:.2f} mÅ"
                         f" (放宽 --tol 再看是否只是弛豫/写文件精度)")
    lines.append("")
    lines.append(f"结论(几何)  : {res['verdict']}")
    if res["note"]:
        lines.append(f"说明        : {res['note']}")
    lines.append(f"边界        : {res['disclaimer']}")
    return "\n".join(lines)


# ---------------------------------------------------------------- selftest

def _write_poscar(path, lat, species, counts, coords, comment="selftest"):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"{comment}\n1.0\n")
        for v in lat:
            fh.write("  ".join(f"{x:.10f}" for x in v) + "\n")
        fh.write(" ".join(species) + "\n")
        fh.write(" ".join(str(c) for c in counts) + "\n")
        fh.write("Direct\n")
        for c in coords:
            fh.write("  ".join(f"{x:.10f}" for x in c) + "\n")


def _selftest():
    import tempfile
    ok = True
    td = tempfile.mkdtemp(prefix="altermag_selftest_")
    sq = [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 12.0]]

    def run(name, lat, species, counts, coords, element, expect_verdict, spins=None):
        nonlocal ok
        p = os.path.join(td, name + ".vasp")
        _write_poscar(p, lat, species, counts, coords)
        lat2, atoms, frac, _ = MO.parse_poscar(p)
        res = check(lat2, atoms, frac, element, spins=spins)
        good = res["verdict"] == expect_verdict
        print(f"  [{'PASS' if good else 'FAIL'}] {name:<22} -> {res['verdict'][:34]}"
              f"{'' if good else '   期望: ' + expect_verdict[:34]}")
        ok = ok and good
        return res

    # (a) 平移连接: 非磁原子让 (0,1/2,0) 成为晶体对称操作 -> 常规共线反铁磁
    run("conv_translation", sq, ["Fe", "Sn"], [2, 2],
        [[0, 0, 0.5], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.5]], "Fe",
        VERDICT_CONVENTIONAL)

    # (b) 旋转连接(金红石式 4_2 螺旋): (1/2,1/2,1/2) 不是晶体对称操作, 螺旋才是
    tet = [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 6.0]]
    run("screw_rotation", tet, ["Fe", "O"], [2, 4],
        [[0, 0, 0], [0.5, 0.5, 0.5],
         [0.3, 0.3, 0.0], [0.7, 0.7, 0.0], [0.8, 0.2, 0.5], [0.2, 0.8, 0.5]], "Fe",
        VERDICT_CANDIDATE)

    # (c) 反演连接: 三斜胞(点群只有 ±I)里唯一连接操作是反演 -> 宇称禁止
    #     注意: 共面 2D 胞里「反演」与「面内 C2」是同一个操作 ⇒ 二维单层几乎总落到「旋转连接」,
    #     PARITY 分支只在 3D(非共面)输入上出现(见文件头说明)。
    tri = [[3.0, 0, 0], [0.5, 3.0, 0], [0.3, 0.4, 10.0]]
    run("inversion_only_3D", tri, ["Fe", "Sn"], [2, 2],
        [[0.1, 0.2, 0.3], [0.9, 0.8, 0.7], [0.4, 0.1, 0.6], [0.6, 0.9, 0.4]], "Fe",
        VERDICT_PARITY)

    # (c') 反演与面内 C2 **同时**连接(全部原子 z=0.5 的共面胞): 旋转连接优先 -> 几何候选
    run("inversion_plus_C2", sq, ["Fe", "Sn"], [2, 2],
        [[0.2, 0.3, 0.5], [0.8, 0.7, 0.5], [0.4, 0.1, 0.5], [0.6, 0.9, 0.5]], "Fe",
        VERDICT_CANDIDATE)

    # (d) 无连接: 两个磁位点不靠任何晶体对称操作相连
    run("no_connection", sq, ["Fe", "Sn"], [2, 2],
        [[0.2, 0.3, 0.5], [0.55, 0.7, 0.5], [0.4, 0.1, 0.5], [0.9, 0.9, 0.5]], "Fe",
        VERDICT_NONE)

    # (e) 不适用: 三套等价类(kagome 的真实 fixture)
    fix = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "fixtures",
                       "Fe3Sn_POSCAR")
    if os.path.exists(fix):
        lat3, atoms3, frac3, _ = MO.parse_poscar(fix)
        r3 = check(lat3, atoms3, frac3, "Fe")
        good = r3["verdict"] == VERDICT_NA and r3["n_sub"] == 3
        print(f"  [{'PASS' if good else 'FAIL'}] {'Fe3Sn(kagome)':<22} -> n_sub={r3['n_sub']}, "
              f"{r3['verdict'][:26]}")
        ok = ok and good
    else:
        print("  [SKIP] Fe3Sn_POSCAR 不存在")
        ok = False

    # (f) 自旋给反了(不补偿)也要能看出来
    r = run("non_compensated", sq, ["Fe", "Sn"], [2, 2],
            [[0, 0, 0.5], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.5]], "Fe",
            VERDICT_CONVENTIONAL, spins=[1, 1])
    good = (r["compensated"] is False)
    print(f"  [{'PASS' if good else 'FAIL'}] {'compensation flag':<22} -> compensated={r['compensated']}")
    ok = ok and good

    # (g) 点群自检: 正方晶格要有 16 个保长整数矩阵(含 ±I)
    ops = point_ops(sq)
    good = len(ops) == 16 and any(np.array_equal(M, -np.eye(3, dtype=int)) for M in ops)
    print(f"  [{'PASS' if good else 'FAIL'}] {'point group (square)':<22} -> {len(ops)} 个点操作")
    ok = ok and good

    # (h) **第 6~7 位小数的坐标抖动不能毁掉判据**: 同一个平移连接, 把 B 组坐标
    #     抖 3e-6 (真实结构文件/DFT 弛豫的常见量级) 后仍必须认出平移连接。
    #     这是把「查四舍六入元组」换成「按 Å 距离判定」的回归锚点。
    jit = 3e-6
    r = run("jit_translation", sq, ["Fe", "Sn"], [2, 2],
            [[0, 0, 0.5], [0, 0.5, 0.5 + jit], [0.5, 0, 0.5], [0.5, 0.5, 0.5]], "Fe",
            VERDICT_CONVENTIONAL)
    good = r.get("worst_mismatch_ang", 1.0) < 1e-2
    print(f"  [{'PASS' if good else 'FAIL'}] {'jitter residue':<22} -> "
          f"worst={r.get('worst_mismatch_ang')} Å")
    ok = ok and good

    # (i) 3D 胞里反演与真旋转**同时**连接 ⇒ 反演在场即宇称禁止 (不能因为「还找到了旋转」就改判候选)
    r = run("inv_and_rot_3D", tri, ["Fe", "Sn"], [2, 2],
            [[0.1, 0.2, 0.3], [0.9, 0.8, 0.7], [0.4, 0.1, 0.6], [0.6, 0.9, 0.4]], "Fe",
            VERDICT_PARITY)
    good = ("inversion" in {o["kind"] for o in r["ops"]})
    print(f"  [{'PASS' if good else 'FAIL'}] {'inversion in ops':<22} -> "
          f"{sorted({o['kind'] for o in r['ops']})}")
    ok = ok and good

    # (j) 无连接时也要给出「最接近者 + 残差」, 供人判断是不是精度问题
    r = run("no_connection", sq, ["Fe", "Sn"], [2, 2],
            [[0.2, 0.3, 0.5], [0.55, 0.7, 0.5], [0.4, 0.1, 0.5], [0.9, 0.9, 0.5]], "Fe",
            VERDICT_NONE)
    good = ("best_mismatch_ang" in r and r["best_mismatch_ang"] > 0.1
            and r.get("best_candidate", {}).get("kind"))
    print(f"  [{'PASS' if good else 'FAIL'}] {'best candidate report':<22} -> "
          f"{r.get('best_candidate', {}).get('kind')} "
          f"{r.get('best_mismatch_ang')} Å")
    ok = ok and good

    import shutil
    shutil.rmtree(td, ignore_errors=True)
    print("ALL PASS" if ok else "FAIL")
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--poscar", help="POSCAR 路径")
    ap.add_argument("--element", help="磁性元素符号 (默认取 POSCAR 里第一个元素)")
    ap.add_argument("--spins", help="两个亚晶格的自旋符号, 例如 +1,-1 (默认 +1,-1)")
    ap.add_argument("--tol", type=float, default=1e-3,
                    help="对称性判据容差, **单位 Å** (映射后原子到同元素最近原子的笛卡尔距离); "
                         "默认 1e-3, DFT 弛豫/实验结构可放宽到 0.01~0.05")
    ap.add_argument("--json", help="把结果写成 JSON")
    ap.add_argument("--selftest", action="store_true", help="跑合成结构的自检")
    args = ap.parse_args()

    if args.selftest:
        print("[altermagnet_check] selftest")
        return _selftest()
    if not args.poscar:
        ap.error("需要 --poscar (或用 --selftest)")

    lat, atoms, frac, _ = MO.parse_poscar(args.poscar)
    element = args.element or atoms[0]
    if element not in atoms:
        print(f"[error] POSCAR 里没有元素 {element} (有: {sorted(set(atoms))})", file=sys.stderr)
        return 2
    spins = None
    if args.spins:
        try:
            spins = [int(x) for x in args.spins.replace(" ", "").split(",") if x]
        except ValueError:
            print("[error] --spins 形如 +1,-1", file=sys.stderr)
            return 2
    res = check(lat, atoms, frac, element, spins=spins, tol=args.tol)
    print(fmt_report(res, poscar=args.poscar))
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(res, fh, ensure_ascii=False, indent=2)
        print(f"\n>>> JSON 已写 {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
