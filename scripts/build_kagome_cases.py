#!/usr/bin/env python
"""build_kagome_cases.py
六角 kagome 晶格 (三个磁位/胞, 最近邻 z=4) 单离子磁体: 超胞 + 磁序枚举 + J1/J2/J3 能量映射系数。

与其它晶格的区别:
  - 三角/蜂窝是 60° 晶格, 1~2 个磁位/胞; 方形/矩形是 90°, 1~2 个磁位/胞;
  - kagome 是 **3 个磁位/胞的受挫晶格**: 每个磁位 4 个最近邻, 由**上/下两种三角形**拼成,
    低能态可以是共面 120° 非共线序 (q=0 或 √3×√3)。
  典型材料: Fe3Sn (S≈1, T_C≈660 K 的金属 kagome) / 以及各类 M3X 家族。

模型: H = -J1 Σ_<ij> S·S - J2 Σ_<<ij>> S·S - J3 Σ_<<<ij>>> S·S   (J>0 = 铁磁)
  距离壳按升序 (以 |a1| = a 为单位):
    J1 = a/2     每格点 4 根 (三角形内最近邻, 上/下三角同距离)
    J2 = a√3/2   每格点 4 根 (跨六角)
    J3 = a       每格点 6 根 (同亚晶格/对角)
  理想 kagome 里 J1(上三角) == J1'(下三角) —— 对称性使然, 故只有 3 个独立距离壳.
  ⚠ 要分辨 J1 与 J1' 必须用 `pair=` 亚晶格过滤的显式壳 (见 magnetic_orders.Shell.pair),
     本脚本不做 (那是另一个 J 自由度, 需要畸变/各向异性结构), 见文件末尾"已知局限"。

能量映射 (每超胞): E = e0 + c1·J1~ + c2·J2~ + c3·J3~,  J~ = J·S²
  c_k = -Σ_{第 k 壳无向键} S_i·S_j = -0.5 · Σ_{第 k 壳有向键} S_i·S_j

⚠ 本晶格的核心陷阱 —— q=0 一族分不开 J1 与 J2:
  q=0 的四个序 (fm / q0_120 / q0_ud / q0_duu) 在设计矩阵里**秩只有 2**:
  它们的 (c1, c2) 恒定成比例 (c1 : c2 = 1 : 1), J1 与 J2 两列线性相关,
  {e0, J1, J2} 不可分离 —— 与方形的"周期 2 序分不开 J3"同构。
  必须引入 **√3×√3** 调制 (sq3_120 / sq3_updntri / sq3_alt) 才能把 J1/J2 分开。
  实测 (3x3x1): 只取 q=0 一族 rank=2, cond≈2.3e34 (数值奇异);
  SELECTED = (fm, q0_120, sq3_updntri, sq3_alt) -> rank=4, cond≈98.5 (< 1e4 闸门)。
  §122 纪律: 先算系数矩阵条件数, 再决定提交哪些构型。

⚠ 第二个陷阱 —— 小超胞里的最小镜像折叠:
  用**最小镜像**枚举键时, 超胞宽度越小, 远壳的像越容易被折叠到近处而被壳容差滤掉。
  本脚本一律用**显式平移搜索** (nshift=1, 与矩形/蜂窝/方形同理), 规格层也必须
  policy="shift_search"; 3x3 超胞里两种策略恰好一致 (已实测), 但不要依赖这一点。

七个磁序 (SELECTED 的 4 个线性无关, 恰好解 {e0, J1, J2, J3}):
  fm          : S = (0,0,+1)              (共线铁磁)
  q0_120      : θ = sub*120°, 面内          (q=0 共面 120° 非共线, 上/下三角同型)
  q0_ud       : σ = (+1,+1,-1) 按亚晶格      (q=0 共线 2up1down)
  q0_duu      : σ = (+1,-1,-1) 按亚晶格      (q=0 共线 1up2down)
  sq3_120     : θ = (sub + ℓ)*120°, ℓ=(n0+n1)%3 (√3×√3 共面 120°)
  sq3_updntri : σ = (-1)^ℓ                  (√3×√3 共线, 胞一相翻转)
  sq3_alt     : σ = (-1)^(sub + ℓ)          (亚晶格 + √3 胞交替)

超胞: 默认 3x3x1。**n0 与 n1 都必须是 3 的倍数**, 否则 √3×√3 的相位 ℓ=(n0+n1)%3
      跨周期边界闭合不干净 (q=0 的序不受影响, 但 sq3_* 会偏)。

用法:
  python build_kagome_cases.py POSCAR --supercell 3 3 1 --outdir cases \
      --element Fe --nonmag Sn --moment 2.0 --spin 1.0 [--shell-tol 0.2] [--selftest]

已知局限 (NEW_TYPES_PLAN P4.1):
  1. 只有 3 个距离壳 ⇒ 拿不到 J1/J1' 的拆分 (需要 pair= 显式壳 + 畸变结构);
  2. S 未在仓库定值: --spin 只影响 J~ 的口径, 反解出的 J 绝对值随 S² 缩放;
  3. 从不共线的角度: q0_120/sq3_120 是**共面** 120° 序, 不需要 SOC (soc=False);
     面外倾斜/螺旋 (cant*/spiral) 才需要 LSORBIT, 那是后续阶段的事。
"""
import argparse
import json
import os
import sys
from itertools import product

import numpy as np

ORDERS = ["fm", "q0_120", "q0_ud", "q0_duu", "sq3_120", "sq3_updntri", "sq3_alt"]
# 4 个线性无关的序 (3x3x1 实测 rank=4, cond≈98.5, 远低于 1e4 闸门); 其余作对照/校验
SELECTED = ["fm", "q0_120", "sq3_updntri", "sq3_alt"]

RULES = {
    "fm": "S = (0,0,+1) (共线铁磁)",
    "q0_120": "kagome q=0: theta = sub*120deg; S=(cos,sin,0)",
    "q0_ud": "kagome q=0: sigma = (+1,+1,-1) 按亚晶格 (2up1down)",
    "q0_duu": "kagome q=0: sigma = (+1,-1,-1) 按亚晶格 (1up2down)",
    "sq3_120": "kagome sqrt3xsqrt3: theta = (sub + (n0+n1)%3)*120deg; S=(cos,sin,0)",
    "sq3_updntri": "kagome sqrt3xsqrt3: sigma = (-1)^((n0+n1)%3)",
    "sq3_alt": "kagome sqrt3xsqrt3: sigma = (-1)^(sub + (n0+n1)%3)",
}

# 每格点键数 (配位数): kagome 的几何指纹 —— 用来抓"壳被静默合并/别名折叠"
EXPECT_Z = (4.0, 4.0, 6.0)
# 3x3x1 下系数的解析理想值 (**per mag**, 与下方 coeff[k]/n 同口径; 每胞值 = 这里的 3 倍)。
# 已用独立键枚举逐位复核 (2026-xx-xx): 每胞键数 (6,6,9), z=(4,4,6)。
IDEAL_33 = {
    "fm": (-2.0, -2.0, -3.0),
    "q0_120": (+1.0, +1.0, -3.0),
    "q0_ud": (+2.0 / 3.0, +2.0 / 3.0, -3.0),
    "q0_duu": (+2.0 / 3.0, +2.0 / 3.0, -3.0),
    "sq3_120": (0.0, 0.0, +1.5),
    "sq3_updntri": (-2.0 / 9.0, -2.0 / 9.0, +1.0),
    "sq3_alt": (-2.0 / 9.0, +2.0 / 3.0, +1.0),
}
Q0_FAMILY = ["fm", "q0_120", "q0_ud", "q0_duu"]


def sign_pow(base):
    """(-1)^base, 对整数数组。"""
    return np.where(np.asarray(base) % 2 == 0, 1.0, -1.0)


def spins_of(n0, n1, sub):
    """七个磁序的自旋矢量 (n,3)。共线序沿 z, 120° 序在面内。"""
    n0 = np.asarray(n0, dtype=int)
    n1 = np.asarray(n1, dtype=int)
    sub = np.asarray(sub, dtype=int)
    m = len(n0)
    zero = np.zeros(m)
    one = np.ones(m)
    lab = (n0 + n1) % 3

    def planar(theta):
        return np.stack([np.cos(theta), np.sin(theta), zero], axis=1)

    def along_z(s):
        return np.stack([zero, zero, s], axis=1)

    return {
        "fm": along_z(one),
        "q0_120": planar(2.0 * np.pi * sub / 3.0),
        "q0_ud": along_z(np.array([1.0, 1.0, -1.0])[sub]),
        "q0_duu": along_z(np.array([1.0, -1.0, -1.0])[sub]),
        "sq3_120": planar(2.0 * np.pi * (sub + lab) / 3.0),
        "sq3_updntri": along_z(sign_pow(lab)),
        "sq3_alt": along_z(sign_pow(sub + lab)),
    }


def parse_poscar(path):
    with open(path, encoding="utf-8") as f:
        lines = [ln for ln in f.read().splitlines() if ln.strip()]
    scale = float(lines[1].strip())
    lat = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)]) * scale
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    idx = 7
    if lines[idx].strip().lower()[:1] == "s":
        idx += 1
    mode = lines[idx].strip().lower()
    idx += 1
    coords, atoms = [], []
    for sp, n in zip(species, counts):
        for _ in range(n):
            coords.append([float(x) for x in lines[idx].split()[:3]])
            atoms.append(sp)
            idx += 1
    coords = np.array(coords)
    if mode[:1] in "ck":
        coords = coords @ np.linalg.inv(lat)
    coords -= np.floor(coords)
    return lat, atoms, coords, mode


def write_poscar(path, lat, atoms, coords, comment="kagome magnetic-order supercell"):
    order = []
    for a in atoms:
        if a not in order:
            order.append(a)
    counts = [atoms.count(sp) for sp in order]
    lines = [comment, "1.00000000000000"]
    for row in lat:
        lines.append(" ".join(f"{x:20.12f}" for x in row))
    lines.append("  ".join(order))
    lines.append("  ".join(str(c) for c in counts))
    lines.append("Direct")
    for sp in order:
        for a, c in zip(atoms, coords):
            if a == sp:
                lines.append(" ".join(f"{x:20.12f}" for x in c))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def make_supercell(lat, atoms, coords, sc):
    sc = np.asarray(sc, dtype=int)
    new_atoms, new_coords, cell_idx = [], [], []
    for n0 in range(sc[0]):
        for n1 in range(sc[1]):
            for n2 in range(sc[2]):
                for a, c in zip(atoms, coords):
                    fc = (c + np.array([n0, n1, n2], dtype=float)) / sc
                    new_atoms.append(a)
                    new_coords.append(fc)
                    cell_idx.append((n0, n1, n2))
    new_lat = lat.copy()
    new_lat[0] *= sc[0]
    new_lat[1] *= sc[1]
    new_lat[2] *= sc[2]
    return new_lat, new_atoms, np.array(new_coords), np.array(cell_idx, dtype=int)


def sublattice_labels(atoms, coords, element):
    """原胞内该元素的亚晶格标签 (按分数坐标字典序, 与 magnetic_orders._sublattice_keys 同序)。

    返回 (原始胞内索引列表, {round 坐标: 标签}, 标签数)。
    """
    idxs = [i for i, a in enumerate(atoms) if a == element]
    keys = sorted({tuple(np.round(coords[i], 6)) for i in idxs})
    return idxs, {k: j for j, k in enumerate(keys)}, len(keys)


def magnetic_shells(idxs, coords, lat, max_shells=3, nshift=2, tol=0.2):
    """显式平移搜索识别磁性亚晶格前 max_shells 个距离壳 (不做最小镜像截断)。"""
    shifts = list(product(range(-nshift, nshift + 1), repeat=3))
    ds = []
    n = len(idxs)
    for a in range(n):
        ci = coords[idxs[a]]
        for b in range(n):
            base = coords[idxs[b]] - ci
            for sh in shifts:
                cart = (base + np.array(sh, dtype=float)) @ lat
                d = float(np.linalg.norm(cart))
                if d < 1e-6:
                    continue
                ds.append(d)
    ds = sorted(ds)
    shells = []
    for d in ds:
        if not shells or d - shells[-1] > tol:
            shells.append(d)
    return shells[:max_shells]


def classify_bonds(idxs, coords, lat, spins, shells, nshift=1, tol=0.3):
    """有向键和 Σ S_i·S_j 按壳归类 (显式平移搜索, 保留像的重数)。"""
    shifts = list(product(range(-nshift, nshift + 1), repeat=3))
    sums = [0.0] * len(shells)
    n = len(idxs)
    for i in range(n):
        ci = coords[idxs[i]]
        for j in range(n):
            base = coords[idxs[j]] - ci
            for sh in shifts:
                cart = (base + np.array(sh, dtype=float)) @ lat
                d = float(np.linalg.norm(cart))
                if d < 1e-6:
                    continue
                k = int(np.argmin([abs(d - s) for s in shells]))
                if abs(d - shells[k]) < tol:
                    sums[k] += float(spins[i] @ spins[j])
    return sums


def _cross_check(spec_kwargs, coeff):
    """与 P0.1 规格层 (magnetic_orders.py) 交叉校验: 两套独立实现必须逐位一致。"""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import magnetic_orders as MO

    spec = MO.Spec(
        key="kagome_selftest", label="kagome (selftest)",
        structure=MO.StructureSpec(
            poscar=spec_kwargs["poscar"], element=spec_kwargs["element"],
            nonmag=tuple(spec_kwargs["nonmag"]), supercell=tuple(spec_kwargs["supercell"]),
            sublattice_names=("A", "B", "C")),
        shells=tuple(MO.Shell(f"J{i + 1}", d) for i, d in enumerate(spec_kwargs["shells"])),
        orders=tuple(MO.Order(nm, (lambda c, nm=nm: spins_of(c.n0, c.n1, c.sub_label)[nm]))
                     for nm in ORDERS),
        selected=tuple(SELECTED), spin=spec_kwargs["spin"],
        # 显式给 shift_search: 与 build_* 的显式平移搜索同口径 (默认 min_image 在小超胞
        # 会把远壳的像折叠到近处)。kagome 的 √3×√3 序依赖胞相位, 更要口径一致。
        policy="shift_search", nshift=1,
    )
    res, ctx = MO.compute(spec)
    worst, where = 0.0, ""
    for nm in ORDERS:
        a = np.asarray(coeff[nm], dtype=float)
        b = np.asarray(res[nm].coeffs, dtype=float)
        if a.shape != b.shape:
            return False, f"{nm}: 形状不符 {a.shape} vs {b.shape}"
        err = float(np.max(np.abs(a - b)))
        if err > worst:
            worst, where = err, nm
    return worst <= 1e-9, f"最大逐位偏差 {worst:.3g} ({where})"


def _order_specs(coeff):
    return [dict(name=nm, kind="planar120" if nm.endswith("120") else "collinear",
                 rule=RULES[nm], coeffs=[float(x) for x in coeff[nm]]) for nm in ORDERS]


def _selftest(args, sc, n, n_atoms, coeff, shells, A, Ainv, cond, sigmas, outdir):
    """无 DFT 的自检: 几何指纹 + q=0 简并 + 设计矩阵 + 合成反解 + 规格层交叉校验。"""
    bad = []
    per = lambda k, m: coeff[k][m] / n  # noqa: E731
    z = tuple(-2.0 * coeff["fm"][m] / n for m in range(3))

    # 1) 几何指纹: 每格点键数
    dz = max(abs(a - b) for a, b in zip(z, EXPECT_Z))
    print(f"\n--selftest: 配位数 z={tuple(round(x, 6) for x in z)} == {EXPECT_Z} "
          f"{'PASS' if dz < 1e-9 else 'FAIL'}")
    if dz >= 1e-9:
        bad.append(f"配位数 {z} 与 kagome 几何不符 (期望 {EXPECT_Z})")

    # 2) 壳距离比 = 1 : √3 : 2 (以 J1 = a/2 为单位)
    ratio = (shells[1] / shells[0], shells[2] / shells[0])
    want = (np.sqrt(3.0), 2.0)
    dr = max(abs(x - y) for x, y in zip(ratio, want))
    print(f"--selftest: 壳距离比 {ratio[0]:.6f}, {ratio[1]:.6f} == (√3, 2) "
          f"{'PASS' if dr < 1e-6 else 'FAIL'}")
    if dr >= 1e-6:
        bad.append(f"壳距离比 {ratio} != (√3, 2): 壳被合并或结构不是 kagome")

    # 3) 3x3x1 解析理想值 (per cell)
    if tuple(sc) == (3, 3, 1):
        err = max(abs(per(k, m) - IDEAL_33[k][m]) for k in IDEAL_33 for m in range(3))
        print(f"--selftest: 3x3x1 解析理想值最大偏差 = {err:.3g} "
              f"{'PASS' if err < 1e-9 else 'FAIL'}")
        if err >= 1e-9:
            bad.append(f"3x3x1 系数偏离解析理想值 (最大 {err:.3g})")

    # 4) q=0 一族的 J1/J2 简并: (c1, c2) 恒定成比例 —— 这才是"分不开 J1/J2"的机制
    r = [coeff[k][0] / coeff[k][1] for k in Q0_FAMILY if abs(coeff[k][1]) > 1e-12]
    deg = max(r) - min(r) if r else 1.0
    print(f"--selftest: q=0 一族 c1/c2 比值离散度 = {deg:.3g} "
          f"{'PASS' if deg < 1e-9 else 'FAIL'}")
    if deg >= 1e-9:
        bad.append("q=0 一族的 c1/c2 不再恒定 —— 简并结构变了, 请复核壳/序定义")

    # 5) 只取 q=0 一族时设计矩阵必然秩亏 (2 而不是 4)
    A_q0 = np.array([[1.0] + coeff[k] for k in Q0_FAMILY])
    rk_q0 = int(np.linalg.matrix_rank(A_q0))
    print(f"--selftest: q=0 一族只给 rank={rk_q0} (期望 2: J1/J2 不可分离) "
          f"{'PASS' if rk_q0 == 2 else 'FAIL'}")
    if rk_q0 != 2:
        bad.append(f"q=0 一族 rank={rk_q0} != 2 —— 简并预期被破坏, 请复核")

    # 6) SELECTED 必须满秩
    rk = int(np.linalg.matrix_rank(A))
    print(f"--selftest: SELECTED rank={rk}/{len(SELECTED)} "
          f"{'PASS' if rk == len(SELECTED) else 'FAIL'}")
    if rk != len(SELECTED):
        bad.append(f"SELECTED 设计矩阵秩亏: rank={rk}")

    # 7) cond / det 闸门
    det = float(np.linalg.det(A))
    print(f"--selftest: cond(A)={cond:.4g} < 1e4 {'PASS' if cond < 1e4 else 'FAIL'}; "
          f"|det(A)|={abs(det):.4g} > 1e-6 {'PASS' if abs(det) > 1e-6 else 'FAIL'}")
    if not (cond < 1e4) or not (abs(det) > 1e-6):
        bad.append(f"设计矩阵闸门不过: cond={cond:.4g}, det={det:.4g}")

    # 8) 七个磁序必须两两不同
    keys = [tuple(np.round(sigmas[k].ravel(), 9).tolist()) for k in ORDERS]
    uniq = len(set(keys))
    print(f"--selftest: 磁序互异 {uniq}/{len(ORDERS)} "
          f"{'PASS' if uniq == len(ORDERS) else 'FAIL'}")
    if uniq != len(ORDERS):
        bad.append(f"磁序数组有重复: 只有 {uniq} 个不同 (期望 {len(ORDERS)})")

    # 9) √3×√3 提供 q=0 一族**没有**的两条独立信息 (这才是打破简并的机制):
    #    (a) sq3_alt 的 c1 != c2  -> 把 J1 与 J2 两列分开;
    #    (b) sq3_updntri 的 c3 与 q=0 一族的 c3 不同 -> 把 J3 与 e0 分开。
    d12 = abs(coeff["sq3_alt"][0] - coeff["sq3_alt"][1])
    d3 = abs(coeff["sq3_updntri"][2] - coeff["fm"][2])
    print(f"--selftest: sq3_alt 的 |c1-c2| = {d12:.6g} > 0 {'PASS' if d12 > 1e-9 else 'FAIL'}; "
          f"sq3_updntri 的 |c3 - c3(q=0)| = {d3:.6g} > 0 {'PASS' if d3 > 1e-9 else 'FAIL'}")
    if d12 <= 1e-9:
        bad.append("sq3_alt 的 c1 == c2: √3×√3 序没有把 J1/J2 两列分开")
    if d3 <= 1e-9:
        bad.append("sq3_updntri 的 c3 与 q=0 一族相同: √3×√3 序没有把 J3/e0 分开")

    # 10) 合成 J -> 反解闭合 (验证 A 可逆 + 系数自洽)
    J_true = [3.0, -0.4, 0.15]  # 合成值 (AFM J1 为主 + 小 J2/J3)
    Jt = [j * args.spin ** 2 for j in J_true]
    E = {k: sum(coeff[k][m] * Jt[m] for m in range(3)) for k in SELECTED}
    Evec = np.array([E[k] for k in SELECTED]) / 1000.0  # meV -> eV
    sol = Ainv @ Evec
    got = [sol[1 + m] / args.spin ** 2 * 1000.0 for m in range(3)]
    inv_ok = all(abs(g - t) < 1e-6 for g, t in zip(got, J_true))
    print(f"--selftest: 合成 J={J_true} -> 反解 {[f'{x:.6f}' for x in got]}  "
          f"{'PASS' if inv_ok else 'FAIL'}")
    if not inv_ok:
        bad.append("合成 J 反解不闭合")

    # 11) 与 P0.1 规格层交叉校验 (两套独立实现必须逐位一致)
    try:
        ok, msg = _cross_check(
            dict(poscar=os.path.abspath(args.poscar), element=args.element,
                 nonmag=args.nonmag, supercell=sc, shells=list(shells),
                 spin=args.spin), coeff)
        print(f"--selftest: 与 P0.1 规格层交叉校验 {'PASS' if ok else 'FAIL'} ({msg})")
        if not ok:
            bad.append(f"与规格层不一致: {msg}")
    except ImportError as exc:
        print(f"--selftest: 跳过与规格层交叉校验 (无法 import magnetic_orders: {exc})")

    # 12) 落盘产物自洽: magmom.txt 单行 + 全原子 3 分量 + 恰好 n 个磁位非零 + POSCAR 注释行
    for nm in ORDERS:
        p = os.path.join(outdir, nm, "magmom.txt")
        with open(p, encoding="utf-8") as f:
            txt = f.read()
        vals = txt.strip().split("=", 1)[1].split()
        nz = sum(1 for k in range(0, len(vals), 3)
                 if max(abs(float(vals[k + c])) for c in range(3)) > 1e-12)
        ok_txt = (txt.count("MAGMOM") == 1 and len(txt.strip().splitlines()) == 1
                  and len(vals) == 3 * n_atoms and nz == n)
        if not ok_txt:
            bad.append(f"{nm}: magmom.txt 格式不符 (单行/值数 {len(vals)} != 3*{n_atoms}"
                       f"/非零磁位 {nz} != {n})")
        c1 = open(os.path.join(outdir, nm, "POSCAR"), encoding="utf-8").readline().strip()
        if c1 != f"{args.element} kagome {nm}":
            bad.append(f"{nm}: POSCAR 注释行 {c1!r} 不符")
    print(f"--selftest: 落盘产物 (magmom.txt 单行 {3 * n_atoms} 值 / {n} 个非零磁位 / "
          f"POSCAR 注释行) "
          f"{'PASS' if not [b for b in bad if 'magmom.txt' in b or 'POSCAR' in b] else 'FAIL'}")

    # 13) 落盘 JSON 的 order_spec 系数与现场重算逐位一致
    with open(os.path.join(outdir, "magnetic_orders.json"), encoding="utf-8") as f:
        rep = json.load(f)
    dj = max(abs(float(rep["order_spec"][i]["coeffs"][m]) - coeff[ORDERS[i]][m])
             for i in range(len(ORDERS)) for m in range(3))
    print(f"--selftest: 落盘 order_spec 系数偏差 = {dj:.3g} "
          f"{'PASS' if dj < 1e-12 else 'FAIL'}")
    if dj >= 1e-12:
        bad.append(f"落盘 order_spec 与现场重算不一致 ({dj:.3g})")

    # 14) per_site_bonds 与超胞无关 (几何不变量)
    pb = rep["per_site_bonds"]
    dpb = max(abs(pb[nm] - EXPECT_Z[m]) for m, nm in enumerate(("J1", "J2", "J3")))
    print(f"--selftest: 落盘 per_site_bonds={ {k: round(v, 6) for k, v in pb.items()} } "
          f"== {EXPECT_Z} {'PASS' if dpb < 1e-9 else 'FAIL'}")
    if dpb >= 1e-9:
        bad.append(f"per_site_bonds 与 kagome 几何不符 (最大偏差 {dpb:.3g})")

    print(f"\n--selftest: {'ALL PASS' if not bad else 'FAIL'}")
    for b in bad:
        print(f"  !! {b}")
    if bad:
        raise SystemExit(1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poscar")
    ap.add_argument("--supercell", type=int, nargs=3, default=[3, 3, 1])
    ap.add_argument("--outdir", default="cases")
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--moment", type=float, default=2.0)
    ap.add_argument("--spin", type=float, default=1.0)
    ap.add_argument("--shell-tol", type=float, default=0.2)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    sc = [int(x) for x in args.supercell]
    if sc[0] % 3 or sc[1] % 3:
        raise SystemExit(
            f"--supercell {sc} 的前两维必须是 3 的倍数: √3×√3 序的相位 ℓ=(n0+n1)%3 "
            "只在 3 的倍数超胞里闭合, 否则 sq3_* 跨周期边界不干净 (J1/J2 分不开)。")

    lat, atoms, coords, mode = parse_poscar(args.poscar)
    _idxs, sub_map, n_sub = sublattice_labels(atoms, coords, args.element)
    if n_sub != 3:
        raise SystemExit(
            f"原胞里有 {n_sub} 个 {args.element} 磁位, kagome 必须恰好 3 个/胞 "
            f"(三角形内三个亚晶格)。这不是 kagome —— 请改用别的 builder "
            f"(build_square_cases.py / build_magnetic_cases_tri.py / build_rectangular_cases.py)。")

    lat_sc, atoms_sc, coords_sc, cell_idx = make_supercell(lat, atoms, coords, sc)
    mag = [i for i, a in enumerate(atoms_sc) if a == args.element]
    n = len(mag)
    if n == 0:
        raise SystemExit(f"超胞里找不到元素 {args.element!r}; 胞内元素 = {sorted(set(atoms_sc))}")
    shells = magnetic_shells(mag, coords_sc, lat_sc, max_shells=3, tol=args.shell_tol)
    if len(shells) < 3:
        raise SystemExit(f"只找到 {len(shells)} 个磁距壳, 需要 3 (J1/J2/J3); 请加大 --supercell")
    gaps = [shells[i + 1] - shells[i] for i in range(len(shells) - 1)]
    if min(gaps) < 0.3:
        print(f"[warn] 相邻壳间距最小 {min(gaps):.4f} A < 0.3 A: --shell-tol "
              f"{args.shell_tol} 可能把两个壳静默合并, 请复核")
    d_nn, d_2nd, d_3rd = shells

    # 亚晶格标签: 按原胞内分数坐标字典序 (与规格层同序)。
    # 超胞分数坐标必须先乘回超胞尺寸才是原胞分数坐标 —— 否则每个位点都映射到同一个标签。
    sc_f = np.asarray(sc, dtype=float)
    sub = np.array([sub_map[tuple(np.round((coords_sc[i] * sc_f) % 1.0, 6))] for i in mag],
                   dtype=int)
    n0 = cell_idx[mag][:, 0]
    n1 = cell_idx[mag][:, 1]
    spins = spins_of(n0, n1, sub)
    coeff = {k: [-0.5 * x for x in classify_bonds(mag, coords_sc, lat_sc, s, shells)]
             for k, s in spins.items()}

    # 几何指纹闸门: kagome 的配位数必须是 (4,4,6)
    z = tuple(-2.0 * coeff["fm"][m] / n for m in range(3))
    if max(abs(a - b) for a, b in zip(z, EXPECT_Z)) >= 1e-6:
        raise SystemExit(
            f"配位数 z={tuple(round(x, 6) for x in z)} 与 kagome 几何不符 (期望 {EXPECT_Z}): "
            "壳很可能被 --shell-tol 静默合并, 或被超胞别名折叠。请显式指定壳距重建。")

    A = np.array([[1.0] + coeff[k] for k in SELECTED])
    Ainv = np.linalg.inv(A)
    rk = int(np.linalg.matrix_rank(A))
    cond = float(np.linalg.cond(A))
    det = float(np.linalg.det(A))
    if rk != len(SELECTED):
        raise SystemExit(
            f"选定的 {SELECTED} 设计矩阵秩亏 (rank={rk} < {len(SELECTED)}): J1/J2 分不开。\n"
            "  kagome 的 q=0 一族天然分不开 J1/J2, 必须让 √3×√3 序 (sq3_*) 进入 selected; "
            "见本文件 docstring 的简并说明与 NEW_TYPES_PLAN §122。")
    if abs(det) < 1e-6:
        raise SystemExit(f"设计矩阵接近奇异: |det(A)|={abs(det):.4g} < 1e-6; 拒绝反解。")
    if cond > 1e4:
        raise SystemExit(
            f"选定的 4 个磁序设计矩阵病态: cond(A)={cond:.4g} > 1e4; 拒绝反解。\n"
            "  提示: 先算系数矩阵条件数再决定提交哪些构型 (NEW_TYPES_PLAN §122)。")

    os.makedirs(args.outdir, exist_ok=True)
    nonmag_counts = {el: atoms_sc.count(el) for el in args.nonmag}
    sp_order, blocks = [], []
    for a in atoms_sc:
        if a not in sp_order:
            sp_order.append(a)
    for sp in sp_order:
        blocks.append((sp, [i for i, a in enumerate(atoms_sc) if a == sp]))
    pos_in_mag = {i: k for k, i in enumerate(mag)}
    for name in ORDERS:
        d = os.path.join(args.outdir, name)
        os.makedirs(d, exist_ok=True)
        write_poscar(os.path.join(d, "POSCAR"), lat_sc, atoms_sc, coords_sc,
                     comment=f"{args.element} kagome {name}")
        vals = []
        for sp, idxs in blocks:
            if sp == args.element:
                # 全部写成 3 分量 (共线序沿 z) —— 非共线链路口径统一
                vals += [f"{args.moment * spins[name][pos_in_mag[i]][c]:.6f}"
                         for i in idxs for c in range(3)]
            else:
                vals += ["0.0"] * (3 * len(idxs))
        with open(os.path.join(d, "magmom.txt"), "w") as f:
            f.write(f"MAGMOM = {' '.join(vals)}\n")

    # P0.3 自描述块 (让「系数行 ↔ 磁序 ↔ 壳定义」可被下游机器复查)
    order_spec = _order_specs(coeff)
    model_spec = dict(
        schema="model_spec/1", lattice="kagome", label="六角 kagome 单离子",
        element=args.element, nonmag=list(args.nonmag), spin=args.spin, S2=args.spin ** 2,
        supercell=sc, policy="shift_search", nshift=1,
        noncollinear=True, soc=False, n_mag=n, n_basis_mag=3,
        shells=[dict(name="J1", d=d_nn, tol=args.shell_tol, pair=None),
                dict(name="J2", d=d_2nd, tol=args.shell_tol, pair=None),
                dict(name="J3", d=d_3rd, tol=args.shell_tol, pair=None)],
    )

    report = {
        "input": args.poscar, "supercell": sc, "n_atoms": len(atoms_sc),
        "n_mag": n, "element": args.element, "nonmag": args.nonmag,
        "moment": args.moment, "spin": args.spin, "S2": args.spin ** 2,
        "shells_angstrom": {"J1": d_nn, "J2": d_2nd, "J3": d_3rd},
        "selected_orders": SELECTED,
        "all_orders": ORDERS,
        "noncollinear": True,
        "soc": False,
        "coefficients_per_cell": {k: list(v) for k, v in coeff.items()},
        "coefficients_per_mag": {k: [c / n for c in v] for k, v in coeff.items()},
        # 跨口径不变量 (每格点键数): poscar_to_mc3.py 在 2D 投影里按距离把 J 绑到壳上,
        # 而本文件写的是 3D 键长 —— 距离对不上, 但这个三元组必须一致。
        "per_site_bonds": {"J1": z[0], "J2": z[1], "J3": z[2]},
        "inversion_matrix_per_cell": Ainv.tolist(),
        "design_cond": cond,
        "design_rank": rk,
        "design_det": det,
        "shell_tol": args.shell_tol,
        "n_cells": len(atoms_sc) // len(atoms),
        "q0_degeneracy_note": ("q=0 一族 (fm/q0_120/q0_ud/q0_duu) 的设计矩阵秩只有 2: "
                               "J1 与 J2 在该族内比值恒定, 必须靠 √3×√3 序分离 "
                               "(NEW_TYPES_PLAN P4 风险条)"),
        "order_spec": order_spec,
        "model_spec": model_spec,
    }
    with open(os.path.join(args.outdir, "magnetic_orders.json"), "w") as f:
        json.dump(report, f, indent=2)

    print(f"模型: E = e0 + c1*J1~ + c2*J2~ + c3*J3~  (J~=J*S^2, S={args.spin})")
    print(f"超胞 {tuple(sc)} -> {len(atoms_sc)} 原子 ({n} {args.element}), nonmag {nonmag_counts}")
    print(f"距离壳: J1={d_nn:.4f} A (a/2)  J2={d_2nd:.4f} A (a√3/2)  J3={d_3rd:.4f} A (a)")
    print(f"配位数 z = {tuple(round(x, 6) for x in z)}")
    print(f"选定 {SELECTED}, rank={rk}, cond(A)={cond:.4g}, |det|={abs(det):.4g}")
    for k in ORDERS:
        tag = " *" if k in SELECTED else "  "
        print(f" {tag}{k:12s} c1={coeff[k][0]:8.3f} c2={coeff[k][1]:8.3f} c3={coeff[k][2]:8.3f}"
              f"   (per mag: {coeff[k][0] / n:+.3f}, {coeff[k][1] / n:+.3f}, {coeff[k][2] / n:+.3f})")
    print("\n注意: q=0 一族 (fm/q0_120/q0_ud/q0_duu) 的 J1 与 J2 列线性相关 (rank=2) —— "
          "J1/J2 只能靠 √3×√3 序 (sq3_*) 分离.")

    if args.selftest:
        _selftest(args, sc, n, len(atoms_sc), coeff, shells, A, Ainv, cond, spins, args.outdir)


if __name__ == "__main__":
    main()
