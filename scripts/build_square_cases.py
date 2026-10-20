#!/usr/bin/env python
"""build_square_cases.py
方形晶格 (a ⊥ b 且 |a| = |b|) 单离子磁体: 超胞 + 共线磁序枚举 + J1/J2/J3 能量映射系数。

与三角/蜂窝/矩形的区别:
  - 三角/蜂窝是 60° 晶格, 近邻配位 6/3, 受挫, 有 120° 非共线磁序;
  - 矩形是 90° 但 |a| != |b| (三个距离壳互不相同);
  - 方形是 90° 且 |a| = |b|, 也是**非受挫**晶格 (低能态全是共线), 但它有一个
    矩形没有的陷阱: 一切**周期 2** 的序给出的 J3 系数完全相同 (见下)。
  典型材料: K2NiF4 (S=1, 2D 反铁磁, T_N=97.1 K) / K2CuF4 (S=1/2, 2D 铁磁, T_C=6.25 K)
            / La2CuO4 (S=1/2, J≈130 meV 的 2D 反铁磁, 高温超导母体)。

模型: H = -J1 Σ_<ij> S·S - J2 Σ_<<ij>> S·S - J3 Σ_<<<ij>>> S·S  (J>0 = 铁磁)
  距离壳按升序: J1 = a 近邻,  J2 = 对角 |a+b| = a√2,  J3 = 轴向 2a。
  ⚠ 文献里 K2NiF4 / 铜氧化物的 "J2" 指的是**对角**壳 (= 本管线的 J2), 不是轴向 2a;
     而 "J3" 在文献里往往指别的路径, 照抄论文数值前先核对壳的定义。
能量映射 (每超胞): E = e0 + c1·J1~ + c2·J2~ + c3·J3~,  J~ = J·S²
  c_k = -0.5 * Σ_{第 k 壳有向键} σ_i·σ_j   (σ = ±1 共线)

⚠ 本晶格的核心陷阱 —— 为什么必须引入周期 4 的序:
  任何**周期 2** 的序 (fm / g_type / stripe_a / stripe_b) 在 2a 壳上给出的系数
  **完全相同** (因为对周期 2 的 σ, C(2,0) = N² 恒成立)。于是 J3 与 e0 的列
  线性相关, {e0, J3} 不可分离 —— 与 kagome 的 q=0 简并同构, 只是成因不同。
  必须引入沿 a / b 方向**周期 4** 的 p4_a / mix_ab 才能解开。
  (这也是 §122 那条纪律「先算系数矩阵条件数, 再决定提交哪些构型」的由来。)

⚠ 第二个陷阱 —— 2a 壳恰好是 4 宽超胞的半周期:
  4x4x1 里 J3 = 2a 的分数位移正好是 (0.5, 0, 0) 与 (0, 0.5, 0)。用**最小镜像**
  枚举键时 +2a 与 -2a 两个像折叠成**同一个位点**, 每胞只数到 1 根 2a 键, 而无限
  晶格每胞有 2 根 —— 于是 c3 少算一半 (fm: -16 而非 -32), J3 被系统性低估 2 倍。
  所以本脚本的 `classify_bonds` 用**显式平移搜索** (nshift=1) 保留像的重数, 与
  矩形/蜂窝脚本同理; `magnetic_orders.py` 里方形规格也必须 policy="shift_search"。
  J1 = a 与 J2 = a√2 的分数位移是 0.25 / 0.3535, 不到半周期, 两种策略一致。

六个共线磁序 (selected 的 4 个线性无关, 恰好解 {e0, J1, J2, J3}):
  fm      : σ = +1                 (全同向)
  g_type  : σ = (-1)^(n0+n1)       (G 型 / Néel, 棋盘)
  stripe_a: σ = (-1)^n0            (C 型条纹, 沿 a 交替)
  stripe_b: σ = (-1)^n1            (a=b 时与 stripe_a 完全简并, 只作对照/校验)
  p4_a    : σ = (-1)^floor(n0/2)   (周期 4, 打破 J3/e0 简并)
  mix_ab  : σ = (-1)^n0 · (-1)^floor(n1/2)

超胞: 默认 4x4x1。**n0 与 n1 都必须是 4 的倍数**, 否则 p4_a / mix_ab 跨周期边界
      闭合不干净, 系数会偏离理想值。

用法:
  python build_square_cases.py POSCAR --supercell 4 4 1 --outdir cases \
      --element Ni --nonmag K F --moment 2.0 --spin 1.0 [--selftest]
"""
import argparse
import json
import os
import sys
from itertools import product

import numpy as np

ORDERS = ["fm", "g_type", "stripe_a", "stripe_b", "p4_a", "mix_ab"]
# 4 个线性无关的序 (cond ≈ 38 @ 4x4x1, 远低于 1e4 闸门); 其余作对照/校验
SELECTED = ["fm", "stripe_a", "p4_a", "mix_ab"]

RULES = {
    "fm": "sigma_i = +1 (全同向)",
    "g_type": "sigma = (-1)^(n0+n1) (方形 G 型/Néel, 棋盘)",
    "stripe_a": "sigma = (-1)^n0 (方形 C 型条纹, 沿 a 交替)",
    "stripe_b": "sigma = (-1)^n1 (方形 C 型条纹, 沿 b 交替)",
    "p4_a": "sigma = (-1)^floor(n0/2) (方形周期 4; 打破 J3/e0 简并的必要条件)",
    "mix_ab": "sigma = (-1)^n0 * (-1)^floor(n1/2) (周期 2 x 周期 4 混合)",
}


def sign_pow(base):
    """(-1)^base, 对整数数组。"""
    return np.where(np.asarray(base) % 2 == 0, 1.0, -1.0)


def sigmas_of(n0, n1):
    """六个磁序的 σ 数组 (±1)。"""
    n0 = np.asarray(n0, dtype=int)
    n1 = np.asarray(n1, dtype=int)
    return {
        "fm": np.ones(len(n0)),
        "g_type": sign_pow(n0 + n1),
        "stripe_a": sign_pow(n0),
        "stripe_b": sign_pow(n1),
        "p4_a": sign_pow(n0 // 2),
        "mix_ab": sign_pow(n0) * sign_pow(n1 // 2),
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


def write_poscar(path, lat, atoms, coords, comment="square magnetic-order supercell"):
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


def classify_bonds(idxs, coords, lat, sigma, shells, nshift=1):
    """有向键和 Σ σ_i σ_j 按壳归类, 返回与 shells 等长的列表。"""
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
                if abs(d - shells[k]) < 0.3:
                    sums[k] += sigma[i] * sigma[j]
    return sums


def _cross_check(spec_kwargs, coeff, n):
    """与 P0.1 规格层 (magnetic_orders.py) 交叉校验: 两套独立实现必须逐位一致。

    build_* 脚本自带几何实现 (是被回归的**独立基准**), 规格层是共享语义层;
    两者一旦漂移必须立刻炸出来, 而不是等某个产物悄悄错掉。
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import magnetic_orders as MO

    spec = MO.Spec(
        key="sq_selftest", label="方形 (selftest)",
        structure=MO.StructureSpec(
            poscar=spec_kwargs["poscar"], element=spec_kwargs["element"],
            nonmag=tuple(spec_kwargs["nonmag"]), supercell=tuple(spec_kwargs["supercell"])),
        shells=tuple(MO.Shell(f"J{i + 1}", d) for i, d in enumerate(spec_kwargs["shells"])),
        orders=tuple(MO.Order(nm, (lambda c, nm=nm: sigmas_of(c.n0, c.n1)[nm]))
                     for nm in ORDERS),
        selected=tuple(SELECTED), spin=spec_kwargs["spin"],
        # 必须显式给 shift_search: 规格层默认 min_image, 而 4x4 超胞里 2a 恰是半周期,
        # min_image 会把 ±2a 两个像折叠成同一位点, 每胞少算一半 2a 键。
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poscar")
    ap.add_argument("--supercell", type=int, nargs=3, default=[4, 4, 1])
    ap.add_argument("--outdir", default="cases")
    ap.add_argument("--element", default="Ni")
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--moment", type=float, default=2.0)
    ap.add_argument("--spin", type=float, default=1.0)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    sc = [int(x) for x in args.supercell]
    if sc[0] % 4 or sc[1] % 4:
        raise SystemExit(
            f"--supercell {sc} 的前两维必须是 4 的倍数: p4_a/mix_ab 靠周期 4 打破 "
            "J3/e0 简并, 超胞宽度不是 4 的倍数时它们跨周期边界闭合不干净。")
    if sc[0] != sc[1]:
        raise SystemExit(f"--supercell {sc} 的 n0 与 n1 必须相等 (方形晶格 |a|=|b|)。")

    lat, atoms, coords, mode = parse_poscar(args.poscar)  # 已把 Cartesian 归一为 Direct
    lat_sc, atoms_sc, coords_sc, cell_idx = make_supercell(lat, atoms, coords, sc)
    mag = [i for i, a in enumerate(atoms_sc) if a == args.element]
    n = len(mag)
    if n == 0:
        raise SystemExit(f"超胞里找不到元素 {args.element!r}; 胞内元素 = {sorted(set(atoms_sc))}")
    shells = magnetic_shells(mag, coords_sc, lat_sc, max_shells=3)
    if len(shells) < 3:
        raise SystemExit(f"只找到 {len(shells)} 个磁距壳, 需要 3 (J1/J2/J3); 请加大 --supercell")
    d_nn, d_nnn, d_3rd = shells

    n0 = cell_idx[mag][:, 0]
    n1 = cell_idx[mag][:, 1]
    sigmas = sigmas_of(n0, n1)

    sums = {k: classify_bonds(mag, coords_sc, lat_sc, s, shells) for k, s in sigmas.items()}
    # E_mag = e0 + c1*J1~ + c2*J2~ + c3*J3~,  c_k = -0.5 * 有向键和 (有向双计)
    coeff = {k: [-0.5 * s[0], -0.5 * s[1], -0.5 * s[2]] for k, s in sums.items()}

    A = np.array([[1.0] + coeff[k] for k in SELECTED])
    Ainv = np.linalg.inv(A)
    cond = float(np.linalg.cond(A))
    if cond > 1e4:
        raise SystemExit(
            f"选定的 4 个磁序设计矩阵病态: cond(A)={cond:.4g} > 1e4; 拒绝反解。\n"
            "  提示: 先算系数矩阵条件数再决定提交哪些构型 (LATTICE_ROADMAP §122)。")

    os.makedirs(args.outdir, exist_ok=True)
    nonmag_counts = {el: atoms_sc.count(el) for el in args.nonmag}
    # MAGMOM 必须与 POSCAR 的**元素分组顺序**一致 (write_poscar 按首次出现顺序分组)。
    # 不能假设磁性元素排在第一位 —— K2NiF4/La2CuO4 恰好是的, 所以旧写法看不出问题,
    # 但换个 POSCAR (例如 JARVIS 原胞) 就会让 MAGMOM 与原子错位而 DFT 照样跑完。
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
                     comment=f"{args.element} square {name}")
        vals = []
        for sp, idxs in blocks:
            if sp == args.element:
                vals += [f"{args.moment * sigmas[name][pos_in_mag[i]]:.1f}" for i in idxs]
            else:
                vals += ["0.0"] * len(idxs)
        with open(os.path.join(d, "magmom.txt"), "w") as f:
            f.write(f"MAGMOM = {' '.join(vals)}\n")

    # P0.3 自描述块: 让「系数行 ↔ 磁序 ↔ 壳定义」可被下游机器复查 (poscar_to_mc3 --jfile 用它)
    order_spec = [dict(name=nm, kind="collinear", rule=RULES[nm],
                       coeffs=[float(x) for x in coeff[nm]]) for nm in ORDERS]
    model_spec = dict(
        schema="model_spec/1", lattice="square", label="方形单离子",
        element=args.element, nonmag=list(args.nonmag), spin=args.spin, S2=args.spin ** 2,
        supercell=sc, policy="shift_search", nshift=1,
        noncollinear=False, n_mag=n,
        shells=[dict(name="J1", d=d_nn, tol=0.3, pair=None),
                dict(name="J2", d=d_nnn, tol=0.3, pair=None),
                dict(name="J3", d=d_3rd, tol=0.3, pair=None)],
    )

    report = {
        "input": args.poscar, "supercell": sc, "n_atoms": len(atoms_sc),
        "n_mag": n, "element": args.element, "nonmag": args.nonmag,
        "moment": args.moment, "spin": args.spin, "S2": args.spin ** 2,
        "shells_angstrom": {"J1": d_nn, "J2": d_nnn, "J3": d_3rd},
        "selected_orders": SELECTED,
        "all_orders": ORDERS,
        "noncollinear": False,
        "coefficients_per_cell": {k: list(v) for k, v in coeff.items()},
        "coefficients_per_mag": {k: [c / n for c in v] for k, v in coeff.items()},
        # 跨口径不变量 (每格点键数): poscar_to_mc3.py 在 2D 投影里搜索并按距离排序把
        # J1/J2/J3 绑到壳上, 而本文件写的是 3D 键长 —— 距离对不上, 但这个三元组必须
        # 一致。驱动层用 run_2d_magnet.py:verify_mc_shells 据此校验绑定没被换序。
        "per_site_bonds": {"J1": -2.0 * coeff["fm"][0] / n, "J2": -2.0 * coeff["fm"][1] / n,
                           "J3": -2.0 * coeff["fm"][2] / n},
        "inversion_matrix_per_cell": Ainv.tolist(),
        "design_cond": cond,
        "order_spec": order_spec,
        "model_spec": model_spec,
    }
    with open(os.path.join(args.outdir, "magnetic_orders.json"), "w") as f:
        json.dump(report, f, indent=2)

    print(f"模型: E = e0 + c1*J1~ + c2*J2~ + c3*J3~  (J~=J*S^2, S={args.spin})")
    print(f"超胞 {tuple(sc)} -> {len(atoms_sc)} 原子 ({n} {args.element}), nonmag {nonmag_counts}")
    print(f"距离壳: J1={d_nn:.4f} A (a)  J2={d_nnn:.4f} A (a*sqrt2)  J3={d_3rd:.4f} A (2a)")
    print(f"选定 {SELECTED}, cond(A)={cond:.4g}")
    for k in ORDERS:
        tag = " *" if k in SELECTED else "  "
        print(f" {tag}{k:9s} c1={coeff[k][0]:8.3f} c2={coeff[k][1]:8.3f} c3={coeff[k][2]:8.3f}"
              f"   (per mag: {coeff[k][0] / n:+.3f}, {coeff[k][1] / n:+.3f}, {coeff[k][2] / n:+.3f})")
    print("\n注意: 周期 2 的序 (fm/g_type/stripe_*) 三项 c3 完全相同 -> 必须靠 p4_a/mix_ab 分离 J3.")

    if args.selftest:
        bad = []

        # 1) 周期 2 简并不变式: 这是本晶格的核心物理事实, 也是必须用周期 4 序的理由
        c3 = [coeff[k][2] for k in ("fm", "g_type", "stripe_a", "stripe_b")]
        deg = max(c3) - min(c3)
        print(f"\n--selftest: 周期 2 序的 c3 简并度 = {deg:.3g} "
              f"{'PASS' if abs(deg) < 1e-9 else 'FAIL'}")
        if abs(deg) >= 1e-9:
            bad.append("周期 2 序 c3 不再简并 (壳定义或磁序定义被改动?)")

        # 2) a=b => stripe_a 与 stripe_b 完全简并
        ds = max(abs(a - b) for a, b in zip(coeff["stripe_a"], coeff["stripe_b"]))
        print(f"--selftest: stripe_a vs stripe_b 简并度 = {ds:.3g} "
              f"{'PASS' if ds < 1e-9 else 'FAIL'}")
        if ds >= 1e-9:
            bad.append("stripe_a != stripe_b, 但方形晶格 |a|=|b| 下两者应当简并")

        # 3) 六个磁序必须两两不同 (防命名/切片写重复)
        keys = [tuple(sigmas[k].tolist()) for k in ORDERS]
        uniq = len(set(keys))
        print(f"--selftest: 磁序互异 {uniq}/{len(ORDERS)} "
              f"{'PASS' if uniq == len(ORDERS) else 'FAIL'}")
        if uniq != len(ORDERS):
            bad.append(f"磁序数组有重复: 只有 {uniq} 个不同 (期望 {len(ORDERS)})")

        # 4) 合成 J -> 反解 (验证 A 可逆 + 系数自洽)
        J_true = [5.0, -0.5, 0.2]  # 制造一组合成值 (AFM J1 为主 + 小 J2/J3)
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

        # 5) 4x4x1 下的解析理想值 (per mag): 每胞 2 个 NN / 2 个对角 / 2 个 2a 键。
        #    ⚠ 2a 键在 4 宽超胞里恰好是半周期, ±2a 两个像折叠到同一位点;
        #      必须用 shift_search 保留像的重数, 否则每胞只数到 1 根 (c3 少算一半)。
        if tuple(sc) == (4, 4, 1):
            ideal = {"fm": (-2.0, -2.0, -2.0), "g_type": (+2.0, -2.0, -2.0),
                     "stripe_a": (0.0, +2.0, -2.0), "p4_a": (-1.0, 0.0, 0.0),
                     "mix_ab": (+1.0, 0.0, 0.0)}
            err = max(abs(coeff[k][m] / n - ideal[k][m])
                      for k in ideal for m in range(3))
            print(f"--selftest: 4x4x1 解析理想值最大偏差 = {err:.3g} "
                  f"{'PASS' if err < 1e-9 else 'FAIL'}")
            if err >= 1e-9:
                bad.append(f"4x4x1 系数偏离解析理想值 (最大 {err:.3g})")

        # 6) 与 P0.1 规格层交叉校验 (两套独立实现必须逐位一致)
        try:
            ok, msg = _cross_check(
                dict(poscar=os.path.abspath(args.poscar), element=args.element,
                     nonmag=args.nonmag, supercell=sc, shells=list(shells),
                     spin=args.spin), coeff, n)
            print(f"--selftest: 与 P0.1 规格层交叉校验 {'PASS' if ok else 'FAIL'} ({msg})")
            if not ok:
                bad.append(f"与规格层不一致: {msg}")
        except ImportError as exc:
            print(f"--selftest: 跳过与规格层交叉校验 (无法 import magnetic_orders: {exc})")

        # 7) 条件数闸门
        print(f"--selftest: cond(A)={cond:.4g} < 1e4 {'PASS' if cond < 1e4 else 'FAIL'}")

        print(f"\n--selftest: {'ALL PASS' if not bad else 'FAIL'}")
        for b in bad:
            print(f"  !! {b}")
        if bad:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
