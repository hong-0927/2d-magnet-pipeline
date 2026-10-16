#!/usr/bin/env python
"""build_fgt_cases.py — Fe3GeTe2 双亚晶格反解: 超胞 + 磁序枚举 + 系数矩阵 + DFT 输入 + 图 + 哈密顿卡片

为什么需要这个脚本: FGT 单层胞有 3 个 Fe、分属两个不等价位 (Fe I@2c 与 Fe II@4e),
流水线自带的 蜂窝(2 位/胞)/三角(1 位/胞)/三聚体簇 三种模型都不适用 (见 SKILL.md
「多亚晶格铁磁网络」)。这里按**双亚晶格**建模并用总能量反解 J。

模型 (约定与流水线一致: H = -Σ_k J_k Σ_{<ij>_k} S_i·S_j , J>0 = FM):
  A = Fe I  2c : 1 个/胞, 面内 (1/3,2/3)
  B = Fe II 4e : 2 个/胞, 面内 (0,0), 胞内两者相距 2.55 A (**竖直对**)
壳层 (从结构实测):
  J1  A-B  2.63 A  (每个 A 位 6 个)
  J2  B-B  2.55 A  (每个 B 位 1 个, 即胞内竖直对)
  J3  A-A  3.99 A  (每个 A 位 6 个)
  J4  B-B  3.99 A  (每个 B 位 6 个)
  J5  A-B  4.78 A  (每个 A 位 6 个, 默认不纳入, --shells 控制)

能量映射:  E(supercell) = N_cell·E0 + Σ_k C_k·J~_k ,  J~ = J·S²
           C_k = -Σ_{bonds∈shell k} σ_iσ_j            (σ = ±1, 共线)
反解:      solve_fgt_J.py 吃 N_config 个静态总能量解这个线性方程组 (C 矩阵可逆性脚本内校验)。
配套:      fgt_relax_backfill.py  弛豫后的 CONTCAR -> 5 个磁态共用同一结构
           fgt_verify_states.py   反解前核对每态 OUTCAR 最终磁矩 vs 目标 sigma
           (端到端 mock 测试: tests/test_fgt_e2e.py, 不需要 VASP)
"""
import argparse
import json
import os
from itertools import product

import numpy as np

from element_data import ldaul_for

# 壳层定义: (标签, 配对, 距离 A, 描述)  —— 顺序即 J1..J5
SHELL_DEF = [
    ("J1", "AB", 2.63, "A-B 最近邻 (主体铁磁耦合)"),
    ("J2", "BB", 2.55, "B-B 胞内竖直对"),
    ("J3", "AA", 3.99, "A-A 面内三角"),
    ("J4", "BB", 3.99, "B-B 面内三角"),
    ("J5", "AB", 4.78, "A-B 次近邻"),
]
TOL = 0.25  # 归类到壳的距离容差 (A)


def parse_poscar(path):
    lines = [l for l in open(path, encoding="utf-8").read().splitlines() if l.strip()]
    scale = float(lines[1])
    lat = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)]) * scale
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    idx = 7
    if lines[idx].strip().lower()[:1] == "s":
        idx += 1
    mode = lines[idx].strip()
    idx += 1
    atoms = []
    for sp, n in zip(species, counts):
        atoms += [sp] * n
    frac = np.array([[float(x) for x in lines[idx + i].split()[:3]]
                     for i in range(len(atoms))])
    if mode[:1] in "cCkK":
        frac = frac @ np.linalg.inv(lat)
    return lat, atoms, frac, species, counts


def write_poscar(path, lat, atoms, frac, comment):
    order, counts = [], []
    for a in atoms:
        if a not in order:
            order.append(a)
            counts.append(0)
        counts[order.index(a)] += 1
    L = [comment, "1.0"]
    for v in lat:
        L.append("  ".join(f"{x:.10f}" for x in v))
    L.append(" ".join(order))
    L.append(" ".join(str(c) for c in counts))
    L.append("Direct")
    for f in frac:
        L.append("  ".join(f"{x:.10f}" for x in f))
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")


def make_supercell(lat, atoms, frac, sc):
    """显式超胞。返回 (new_lat, atoms, frac, cell_index, is_A)。"""
    new_lat = lat * np.array(sc)[:, None]
    out_a, out_f, cell, isA = [], [], [], []
    for ia, a in enumerate(atoms):
        for n in product(*[range(s) for s in sc]):
            f = (frac[ia] + np.array(n)) / np.array(sc)
            out_a.append(a)
            out_f.append(f)
            cell.append(n)
            isA.append(None)
    return new_lat, out_a, np.array(out_f), np.array(cell), isA


def min_image(d, lat):
    d = d - np.floor(d + 0.5)
    return float(np.linalg.norm(d @ lat))


def _fgt_visibility_offsets(mag_sc, sub_sc, frac, fe_prim, n_cell):
    """B 位"竖直对"的面内可见性偏移: 俯视图里同一面内位置的两个 B 完全重叠, 按上/下
    成员横向错开, 让各自的 ↑/↓ 都看得见 (仅为可见性, 不代表真实面内位移)。"""
    dx = np.zeros(len(mag_sc))
    for k in range(len(mag_sc)):
        if sub_sc[k] == "B":
            base = int(mag_sc[k]) // n_cell
            dx[k] = 0.45 if frac[fe_prim[base]][2] > 0.5 else -0.45
    return dx


def _plot_fgt_png(png, configs, sigma, xy, sub_sc, dx, cell=None):
    """纯 Pillow 兜底渲染 FGT 磁序俯视图 (本机 matplotlib 缺依赖时走这里)。

    与 matplotlib 版同图同信息: 每个磁序一格; 圆 = A(2c), 方 = B(4e); 格内 ↑/↓ = 自旋;
    B 竖直对按 dx 做面内错位以便同时可见。画布走 pil_fig 约定 (figsize*200dpi, 2x 超采样
    后 LANCZOS 缩小)。与 mpl 版的两点差异: 总标题/脚注改成 draw_titles 的换行排版
    (mpl 版是一行超长 suptitle), 点内箭头字号按点径取 7pt (mpl 的 11pt 明显溢出色点)。
    """
    import sys as _sys
    _sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import pil_fig as P
    from PIL import Image, ImageDraw

    n = max(1, len(configs))
    W, H = int(round(3.1 * n * P.DPI)), int(round(3.4 * P.DPI))
    Wd, Hd = W * P.SS, H * P.SS
    img = Image.new("RGB", (Wd, Hd), "white")
    dr = ImageDraw.Draw(img)

    xy = np.asarray(xy, dtype=float)
    pts = xy + np.column_stack([np.asarray(dx, dtype=float), np.zeros(len(xy))])
    corners = np.asarray(cell, dtype=float) if cell is not None else pts
    boxes = P.panel_grid(1, n, Wd, Hd, left=0.02, right=0.98, top=0.17, bottom=0.80,
                         wspace=0.05)
    r_mag = P.marker_radius(34)          # mpl scatter s=34
    f_glyph = P.font(P.pt2px(7.0))
    for b, kind in enumerate(configs):
        x0, y0, pw, ph = boxes[b]
        T = P.panel_transform(x0, y0, pw, ph, corners, 0.3)
        s = sigma(kind)
        dr.rectangle([x0, y0, x0 + pw, y0 + ph], outline=(204, 204, 204),
                     width=max(1, int(round(P.pt2px(0.8)))))
        for k in range(len(pts)):
            up = s[k] > 0
            P.marker_shape(dr, T(tuple(pts[k])), r_mag,
                           ("o" if sub_sc[k] == "A" else "s"),
                           fill=("#1f77b4" if sub_sc[k] == "A" else "#d62728"),
                           edge="black", width=P.pt2px(0.4))
            dr.text(T(tuple(pts[k])), "↑" if up else "↓", font=f_glyph,
                    fill=("white" if up else "black"), anchor="mm")
        dr.text((x0 + pw / 2.0, y0 - P.pt2px(6)), kind, font=P.font(P.pt2px(11)),
                fill="black", anchor="md")
    foot = ("圆 = A (2c) Fe I, 方 = B (4e) Fe II; 点内 ↑/↓ = 自旋; "
            "B 的竖直对按面内偏移画开以便同时可见 (不代表真实面内位移)。\n"
            "数据: 02_supercell/cases/*/POSCAR + magmom.txt + magnetic_orders.json\n"
            "(纯 Pillow 渲染: 本机 matplotlib 缺依赖, 由 build_fgt_cases.py 出图)")
    P.draw_titles(dr, Wd, Hd, "Fe3GeTe2 two-sublattice magnetic orders", foot)
    img = img.resize((W, H), Image.LANCZOS)
    os.makedirs(os.path.dirname(os.path.abspath(png)) or ".", exist_ok=True)
    img.save(png, "PNG", optimize=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poscar",
                    help="单层原胞 POSCAR (必须含两个 Fe 亚晶格: 1 个 2c + 2 个 4e)")
    ap.add_argument("--supercell", type=int, nargs=3, default=[3, 3, 1],
                    help="默认 3x3: 2x2 (边长 8 A) 会让 ~4 A 的同亚晶格键**别名坍缩** -> J3/J4 少算")
    ap.add_argument("--outdir", default=os.path.abspath("Fe3GeTe2_pipeline"),
                    help="项目目录 (默认 ./Fe3GeTe2_pipeline)")
    ap.add_argument("--root", default=None,
                    help="产物根目录 (默认 == --outdir)")
    ap.add_argument("--element", default="Fe")
    ap.add_argument("--nonmag", nargs="+", default=["Ge", "Te"])
    ap.add_argument("--moment", type=float, default=2.0, help="每个 Fe 的初始磁矩 (muB)")
    ap.add_argument("--spin", type=float, default=1.0, help="Fe 的自旋 S (统一 S=1)")
    ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--encut", type=float, default=500.0)
    ap.add_argument("--shells", nargs="+", default=["J1", "J2", "J3", "J4"],
                    help="纳入反解的壳 (默认 J1..J4; 加 J5 需要更多磁态)")
    args = ap.parse_args()
    if args.root is None:
        args.root = args.outdir

    lat, atoms, frac, species, counts = parse_poscar(args.poscar)
    sc = tuple(args.supercell)
    n_cell = sc[0] * sc[1] * sc[2]

    # --- 原胞里区分 A(2c) / B(4e): 面内 (1/3,2/3) 为 A, (0,0) 为 B ---
    fe_prim = [i for i, a in enumerate(atoms) if a == args.element]
    sub_prim = []
    for i in fe_prim:
        f = frac[i]
        a_like = abs(f[0] - 1 / 3) < 1e-3 and abs(f[1] - 2 / 3) < 1e-3
        sub_prim.append("A" if a_like else "B")
    print(f"原胞磁性原子: {len(fe_prim)} 个 -> 亚晶格 {sub_prim}")
    assert sorted(sub_prim) == ["A", "B", "B"], "本脚本只处理 Fe I@2c ×1 + Fe II@4e ×2 的 FGT 型胞"
    # 唯一的"位点键": 两个 B 同属亚晶格 B (共用 J), 但映射到超胞时必须区分开,
    # 否则 (亚晶格, 胞) 作键会互相覆盖 -> 键张冠李戴。
    site_prim = []
    for s in sub_prim:
        if s == "A":
            site_prim.append("A")
        else:
            site_prim.append(f"B{sum(1 for x in site_prim if x.startswith('B'))}")

    # --- 壳层: 从原胞按显式平移统计, 按 (配对, 距离) 归类 ---
    shift = list(product(range(-2, 3), repeat=3))
    shell_bonds = {name: [] for name, _, _, _ in SHELL_DEF}  # (i, j, shift) 代表键
    for ii, i in enumerate(fe_prim):
        for jj, j in enumerate(fe_prim):
            pair = "".join(sorted([sub_prim[ii], sub_prim[jj]]))
            if pair == "AA":
                pair = "AA"
            for sh in shift:
                cart = (frac[j] - frac[i] + np.array(sh, float)) @ lat
                d = float(np.linalg.norm(cart))
                if d < 1e-6:
                    continue
                for name, pr, dist, _ in SHELL_DEF:
                    if pr == pair and abs(d - dist) < TOL:
                        shell_bonds[name].append((ii, jj, sh, d))
    # 去重 (i,j,sh) 与其反向 (-j,-i,-sh) 是同一根键
    for name in shell_bonds:
        seen, uniq = set(), []
        for (ii, jj, sh, d) in shell_bonds[name]:
            key = (min(ii, jj), max(ii, jj), tuple(sh) if ii <= jj else tuple(-np.array(sh)))
            if key in seen:
                continue
            seen.add(key)
            uniq.append((ii, jj, sh, d))
        shell_bonds[name] = uniq

    use = [s for s in args.shells]
    print("\n壳层 (原胞内唯一键数, 已按配对与距离分离):")
    for name, pr, dist, desc in SHELL_DEF:
        mark = "纳入" if name in use else "跳过"
        print(f"  {name} {pr} {dist:.2f} A : {len(shell_bonds[name]):2d} 根/胞   [{mark}] {desc}")

    # --- 超胞 ---
    new_lat, sc_atoms, sc_frac, sc_cell, _ = make_supercell(lat, atoms, frac, sc)
    mag_sc = [i for i, a in enumerate(sc_atoms) if a == args.element]
    # 每个超胞磁性原子 -> (亚晶格, 原胞索引)
    sub_sc, cell_of = [], []
    for i in mag_sc:
        base = i // n_cell  # make_supercell 的循环顺序: 外层原子, 内层平移
        sub_sc.append(sub_prim[base])
        cell_of.append(tuple(sc_cell[i]))
    sub_sc = np.array(sub_sc)
    print(f"\n超胞 {sc} -> {len(sc_atoms)} 原子, 其中 {len(mag_sc)} 个 {args.element} "
          f"(A {int((sub_sc == 'A').sum())} / B {int((sub_sc == 'B').sum())}), {n_cell} 个胞")

    # --- 磁序构型 (σ = ±1) ---
    def sigma(kind):
        s = np.ones(len(mag_sc))
        for k, i in enumerate(mag_sc):
            c = cell_of[k]
            if kind == "fm":
                pass
            elif kind == "afm_ab":                      # A↑ B↓  (亚铁磁)
                s[k] = +1 if sub_sc[k] == "A" else -1
            elif kind == "afm_b_vert":                  # 胞内竖直 B-B 对反向
                base = int(mag_sc[k]) // n_cell
                z = frac[fe_prim[base]][2]
                s[k] = -1 if (sub_sc[k] == "B" and z > 0.5) else +1
            elif kind == "afm_a":                       # A 面内棋盘反向
                s[k] = (-1) ** (c[0] + c[1]) if sub_sc[k] == "A" else +1
            elif kind == "afm_b_plane":                 # B 面内棋盘反向
                s[k] = (-1) ** (c[0] + c[1]) if sub_sc[k] == "B" else +1
            else:
                raise SystemExit(f"未知磁序 {kind}")
        return s

    CONFIGS = ["fm", "afm_ab", "afm_b_vert", "afm_a", "afm_b_plane"]
    DESC = {"fm": "全部 Fe 同向 (参考态, 定 E0)",
            "afm_ab": "A↑B↓ 亚铁磁 (主要定 J1)",
            "afm_b_vert": "A↑, B 的胞内竖直对反向 (定 J2)",
            "afm_a": "A 面内棋盘反向, B↑ (定 J3)",
            "afm_b_plane": "A↑, B 面内棋盘反向 (定 J4)"}

    # --- 系数矩阵: C_k = -Σ_{唯一键} σ_iσ_j (超胞内, 最小镜像) ---
    def bonds_in_supercell(name):
        """把原胞壳键映射到超胞, 返回**唯一无向键**列表 (每根键只算一次)。

        必须去重: 原胞壳表里 (ii,jj,sh) 与 (jj,ii,-sh) 都代表同一类键, 逐胞展开时
        会把同一根键生成两次 -> 系数翻倍 -> 反解出的 J 少一半。
        """
        idx_of = {}
        for k, i in enumerate(mag_sc):
            base = int(i) // n_cell
            idx_of[(site_prim[base],) + tuple(cell_of[k])] = k
        bonds = set()
        for c in product(*[range(s) for s in sc]):
            for (ii, jj, sh, d) in shell_bonds[name]:
                cj = tuple(int(x) for x in
                           (np.array(c, float) + np.array(sh, float)) % np.array(sc, float))
                ki = idx_of.get((site_prim[ii],) + tuple(c))
                kj = idx_of.get((site_prim[jj],) + cj)
                if ki is None or kj is None or ki == kj:
                    continue
                bonds.add((min(ki, kj), max(ki, kj)))
        return sorted(bonds)

    shell_bonds_sc = {n: bonds_in_supercell(n) for n in use}

    # 别名自检 (几何判据): 任何纳入的键长都必须 < 超胞面内边长的一半, 否则周期像会把
    # **不同的键折到同一对原子**上, 系数被系统性少算 -> 反解出的 J 偏大。
    # 实测: 2x2 (边长 8.0 A) 对 3.99 A 的 J3/J4 已经临界坍缩, 对 4.78 A 的 J5 必然坍缩。
    Lmin = min(sc[0] * float(np.linalg.norm(new_lat[0])),
               sc[1] * float(np.linalg.norm(new_lat[1])))
    rmax = max(d for (n, pr, d, de) in SHELL_DEF if n in use)
    if rmax >= Lmin / 2.0:
        raise SystemExit(
            f"[build_fgt_cases] 超胞 {sc} 太小: 最远纳入壳 r_max = {rmax:.2f} A "
            f">= 面内边长一半 {Lmin/2:.2f} A -> 键会别名坍缩, 反解出的 J 不可信。\n"
            f"  请加大 --supercell (边长需 > 2*r_max = {2*rmax:.1f} A, 即 3x3 及以上), "
            f"或减少 --shells。")
    bonds_per_cell = {}
    for n in use:
        tot = len(shell_bonds_sc[n])
        if tot % n_cell:
            raise SystemExit(f"[build_fgt_cases] 壳 {n} 的唯一键数 {tot} 不是胞数 {n_cell} 的整数倍")
        bonds_per_cell[n] = tot // n_cell
    print(f"\n纳入的壳 (超胞 {sc} 内唯一键 / 胞):")
    for n in use:
        pr = [p for (nm, p, d, de) in SHELL_DEF if nm == n][0]
        dist = [d for (nm, p, d, de) in SHELL_DEF if nm == n][0]
        print(f"  {n} {pr} {dist:.2f} A : {bonds_per_cell[n]:2d} 根/胞, 超胞内 {len(shell_bonds_sc[n])} 根")
    C = np.zeros((len(CONFIGS), len(use)))
    for r, kind in enumerate(CONFIGS):
        s = sigma(kind)
        for q, n in enumerate(use):
            C[r, q] = -sum(s[i] * s[j] for (i, j) in shell_bonds_sc[n])

    # 可逆性/条件数校验
    M = np.column_stack([np.ones(len(CONFIGS)), C])   # [E0, J1..Jk]
    cond = float(np.linalg.cond(M))
    print(f"\n系数矩阵 M (行=磁态, 列=[E0, {', '.join(use)}]):")
    for r, kind in enumerate(CONFIGS):
        print(f"  {kind:12s} " + "  ".join(f"{v:8.2f}" for v in M[r]))
    print(f"  条件数 cond(M) = {cond:.2f} " + ("[可逆 OK]" if cond < 1e6 else "[WARN 病态]"))
    if abs(np.linalg.det(M)) < 1e-9:
        raise SystemExit("系数矩阵奇异: 磁态不能区分这些壳, 请增删 --shells 或补磁态")

    # --- 写 cases/ ---
    cases = os.path.join(args.outdir, "02_supercell", "cases")
    os.makedirs(cases, exist_ok=True)
    mag_list = []
    for r, kind in enumerate(CONFIGS):
        s = sigma(kind)
        d = os.path.join(cases, kind)
        os.makedirs(d, exist_ok=True)
        write_poscar(os.path.join(d, "POSCAR"), new_lat, sc_atoms, sc_frac,
                     f"Fe3GeTe2 {sc} supercell, order={kind}")
        mom = np.zeros(len(sc_atoms))
        j = 0
        for i in range(len(sc_atoms)):
            if sc_atoms[i] == args.element:
                mom[i] = s[j] * args.moment
                j += 1
        with open(os.path.join(d, "magmom.txt"), "w", encoding="utf-8") as fh:
            fh.write(f"# {kind}: {DESC[kind]}\n")
            fh.write("# 顺序与 POSCAR 一致; 每行 元素 序号 亚晶格 sigma 磁矩(muB)\n")
            j = 0
            for i, a in enumerate(sc_atoms):
                if a == args.element:
                    fh.write(f"{a} {i+1} {sub_sc[j]} sigma={int(s[j]):+d} {mom[i]:+.2f}\n")
                    j += 1
                else:
                    fh.write(f"{a} {i+1} - - {mom[i]:+.2f}\n")
        mag_list.append({"name": kind, "desc": DESC[kind], "coeffs": list(C[r]),
                         "m_total_muB": float(mom.sum())})

    report = {
        "system": "Fe3GeTe2 monolayer (two-sublattice)",
        "model": "H = -sum_k J_k sum_<ij>_k S_i.S_j  (J>0 = FM); E(supercell) = N_cell*E0 + sum_k C_k*J~_k, J~=J*S^2",
        "sublattices": {"A": "Fe I 2c (1/cell)", "B": "Fe II 4e (2/cell)"},
        "supercell": list(sc),
        "n_cell": n_cell,
        "n_mag": len(mag_sc),
        "spin_S": args.spin,
        "moment_muB": args.moment,
        "shells": [{"label": n, "pair": pr, "d_A": dist,
                    "bonds_per_cell": bonds_per_cell.get(n), "desc": desc, "used": n in use}
                   for (n, pr, dist, desc) in SHELL_DEF],
        "configurations": mag_list,
        "coefficient_matrix": {"rows": CONFIGS, "cols": ["E0"] + use,
                               "M": M.tolist(), "cond": cond},
        "poscar": os.path.abspath(args.poscar),
    }
    with open(os.path.join(cases, "magnetic_orders.json"), "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
    print(f"\n写出 {cases}/<order>/POSCAR + magmom.txt, magnetic_orders.json")

    # --- DFT 输入 ---
    dft = os.path.join(args.outdir, "03_dft")
    os.makedirs(dft, exist_ok=True)
    ldaul = " ".join(ldaul_for(s) if s == args.element else "-1" for s in species)
    ldauu = " ".join(f"{args.U:.1f}" if s == args.element else "0.0" for s in species)
    ldauj = " ".join("0.0" for _ in species)

    def incar(magmom_line, system, ediff):
        return f"""SYSTEM   = {system}
PREC     = Normal
ENCUT    = {args.encut:.0f}
EDIFF    = {ediff}
ALGO     = Normal
NELM     = 200
ISMEAR   = 0
SIGMA    = 0.05

ISPIN    = 2
MAGMOM   = {magmom_line}

LDAU     = .TRUE.
LDAUTYPE = 2
LDAUL    = {ldaul}
LDAUU    = {ldauu}
LDAUJ    = {ldauj}
LMAXMIX  = 4

LREAL    = .FALSE.
LWAVE    = .FALSE.
LCHARG   = .FALSE.
LORBIT   = 11
NPAR     = 4
"""

    nmag = len(mag_sc)
    fm_line = " ".join([f"{args.moment:.1f}"] * nmag +
                       ["0.0"] * (len(sc_atoms) - nmag))
    with open(os.path.join(dft, "INCAR.static"), "w", encoding="utf-8") as fh:
        fh.write(incar(fm_line, "Fe3GeTe2 monolayer static (per-order MAGMOM in runs/)", "1E-6"))
    relax_txt = (incar(fm_line, "Fe3GeTe2 monolayer FM relax (先用 CONTCAR 再跑各磁态)", "1E-5")
                 + "IBRION   = 2\nNSW      = 60\nEDIFFG   = -0.02\nISIF     = 2\n")
    with open(os.path.join(dft, "INCAR.relax"), "w", encoding="utf-8") as fh:
        fh.write(relax_txt)
    with open(os.path.join(dft, "KPOINTS"), "w", encoding="utf-8") as fh:
        fh.write(f"Gamma-centered {sc[0]}x{sc[1]}x1 supercell, vacuum along c\n0\nGamma\n"
                 f"{max(3, sc[0])} {max(3, sc[1])} 1\n0 0 0\n")
    pa = {"Fe": "Fe_pv", "Ge": "Ge_d", "Te": "Te"}
    with open(os.path.join(dft, "POTCAR.spec"), "w", encoding="utf-8") as fh:
        fh.write("# POTCAR 顺序必须与 POSCAR/magmom.txt 一致: " + " ".join(species) + "\n\n")
        fh.write("| 顺序 | 元素 | 建议 PAW |\n|---|---|---|\n")
        for i, s in enumerate(species):
            fh.write(f"| {i+1} | {s} | {pa.get(s, s)} |\n")

    runs = os.path.join(args.outdir, "runs")
    for r, kind in enumerate(CONFIGS):
        rd = os.path.join(runs, kind)
        os.makedirs(rd, exist_ok=True)
        s = sigma(kind)
        j = 0
        line = []
        for i, a in enumerate(sc_atoms):
            if a == args.element:
                line.append(f"{s[j] * args.moment:+.1f}")
                j += 1
            else:
                line.append("0.0")
        with open(os.path.join(rd, "INCAR"), "w", encoding="utf-8") as fh:
            fh.write(incar(" ".join(line), f"Fe3GeTe2 {kind}", "1E-6"))
        for f in ("POSCAR", "magmom.txt"):
            src = os.path.join(cases, kind, f)
            if os.path.exists(src):
                with open(src, encoding="utf-8") as a, open(os.path.join(rd, f), "w",
                                                            encoding="utf-8") as b:
                    b.write(a.read())
        with open(os.path.join(dft, "KPOINTS"), encoding="utf-8") as a, \
                open(os.path.join(rd, "KPOINTS"), "w", encoding="utf-8") as b:
            b.write(a.read())
    # runs/fm/ 里同时放一份弛豫 INCAR: 该目录的 INCAR 是**静态**的 (IBRION 默认 -1),
    # 直接跑它不会弛豫 —— 这是最容易踩的坑, 所以把弛豫输入摆在旁边并标明。
    fm_relax = os.path.join(runs, "fm", "INCAR.relax")
    if os.path.isdir(os.path.dirname(fm_relax)):
        with open(fm_relax, "w", encoding="utf-8") as fh:
            fh.write(relax_txt)
    print(f"写出 {dft}/{{INCAR.static,INCAR.relax,KPOINTS,POTCAR.spec}} 与 runs/<order>/ "
          f"({len(CONFIGS)} 个磁态, INCAR 里是**逐原子** MAGMOM; "
          f"runs/fm/ 另存 INCAR.relax, 因为 runs/fm/INCAR 是静态的)")

    # --- 图 ---
    figd = os.path.join(args.outdir, "figures")
    os.makedirs(figd, exist_ok=True)
    png = os.path.join(figd, "fgt_supercell_orders.png")
    cart = sc_frac @ new_lat
    xy_mag = cart[mag_sc][:, :2]                       # 磁性原子的面内坐标 (A)
    dx_mag = _fgt_visibility_offsets(mag_sc, sub_sc, frac, fe_prim, n_cell)
    cell2 = np.array([[0.0, 0.0], new_lat[0, :2], new_lat[0, :2] + new_lat[1, :2],
                      new_lat[1, :2]])
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, len(CONFIGS), figsize=(3.1 * len(CONFIGS), 3.4))
        for ax, kind in zip(np.atleast_1d(axes), CONFIGS):
            s = sigma(kind)
            ax.set_title(kind, fontsize=11)
            for k, i in enumerate(mag_sc):
                up = s[k] > 0
                # B 位是"竖直对"(同一面内位置), 俯视图会完全重叠 -> 按上/下成员做面内偏移
                # 让两个成员与各自的 ↑/↓ 都看得见 (仅为可见性, 不代表真实面内位移)。
                x, y = cart[i, 0] + dx_mag[k], cart[i, 1]
                ax.scatter(x, y, s=34,
                           c=("tab:blue" if sub_sc[k] == "A" else "tab:red"),
                           marker=("o" if sub_sc[k] == "A" else "s"),
                           edgecolor="k", linewidth=0.4, zorder=3)
                ax.annotate("↑" if up else "↓", (x, y),
                            ha="center", va="center", fontsize=11, zorder=4,
                            color="white" if up else "black")
            ax.set_aspect("equal")
            ax.set_xticks([]); ax.set_yticks([])
            ax.grid(alpha=0.2)
        fig.suptitle("Fe3GeTe2 two-sublattice magnetic orders "
                     "(circle=A(2c) square=B(4e); up/down arrow = spin; "
                     "the two B of each vertical pair are offset in-plane for visibility)",
                     fontsize=11)
        fig.tight_layout()
        fig.savefig(png, dpi=150, bbox_inches="tight")
        print(f"图 -> {png}")
    except Exception as e:
        print(f"[warn] matplotlib 画图不可用 ({e}); 改用纯 Pillow 兜底")
        try:
            _plot_fgt_png(png, CONFIGS, sigma, xy_mag, sub_sc, dx_mag, cell2)
            print(f"图 -> {png} (纯 Pillow 兜底)")
        except Exception as e2:  # pragma: no cover
            print("[warn] 画图跳过:", e2)

    # --- 哈密顿卡片 ---
    h = [f"# Fe3GeTe2 单层 双亚晶格 Heisenberg 哈密顿量", "",
         "> 本卡片由 `build_fgt_cases.py` 从结构 + 磁序系数矩阵生成，数值均来自 "
         "`02_supercell/cases/magnetic_orders.json`。", "",
         "## 符号式哈密顿量", "", "```text",
         "H = - J1 Σ_{A-B, 2.63Å} S·S - J2 Σ_{B-B, 2.55Å} S·S",
         "    - J3 Σ_{A-A, 3.99Å} S·S - J4 Σ_{B-B, 3.99Å} S·S", "```", "",
         f"亚晶格: **A = Fe I (2c), 1 个/胞**; **B = Fe II (4e), 2 个/胞** "
         f"(胞内两者相距 2.55 Å, 构成竖直对)。超胞 {sc}, {n_cell} 个胞, "
         f"{len(mag_sc)} 个 Fe。", "",
         "## 壳层", "", "| 壳 | 配对 | d (Å) | 键数/胞 | 是否纳入反解 | 说明 |",
         "|---|---|---|---|---|---|"]
    for (n, pr, dist, desc) in SHELL_DEF:
        bpc = bonds_per_cell.get(n)
        h.append(f"| {n} | {pr} | {dist:.2f} | {bpc if bpc else '—'} | "
                 f"{'✅' if n in use else '—'} | {desc} |")
    h += ["", "## 能量映射与反解", "",
          "对每个磁态 (σ=±1 共线) 有", "",
          "```text",
          "E(supercell) = N_cell · E0 + Σ_k C_k · J~_k ,   J~_k = J_k · S²",
          "C_k = - Σ_{bonds ∈ shell k} σ_i σ_j",
          "```", "",
          f"系数矩阵 M = [1 | C] ({len(CONFIGS)} 行 × {len(use) + 1} 列), "
          f"条件数 cond = {cond:.2f}; 反解 = M⁻¹ · E:", "",
          "| 磁态 | " + " | ".join(["E0"] + use) + " | 说明 |",
          "|" + "---|" * (len(use) + 3)]
    for r, kind in enumerate(CONFIGS):
        h.append(f"| {kind} | " + " | ".join(f"{v:.2f}" for v in M[r]) +
                 f" | {DESC[kind]} |")
    h += ["", "## 用法 (5 步, 只有第 2 步需要 VASP)", "",
          "```powershell",
          "# 1) 先弛豫 FM: 用 03_dft/INCAR.relax (== runs/fm/INCAR.relax, 含 IBRION=2/NSW=60)",
          "#    !!! runs/fm/INCAR 是**静态**的 (IBRION 默认 -1), 直接跑它不会弛豫 !!!",
          "#    拿到 CONTCAR 后回填成 5 个磁态共用的弛豫后结构 (原坐标存 POSCAR.pristine):",
          f'python scripts/fgt_relax_backfill.py "{os.path.join(args.outdir, "02_supercell", "cases")}"',
          "# 2) 每个磁态静态自洽: runs/<order>/ (INCAR 里已写逐原子 MAGMOM)",
          "#    —— 只有这一步要跑 VASP (INCAR.static + KPOINTS + POTCAR.spec)",
          "# 3) 反解前核对每个磁态真的收敛到指定 sigma (符号/大小逐原子比):",
          f'python scripts/fgt_verify_states.py "{os.path.join(args.outdir, "02_supercell", "cases")}"',
          "# 4) 反解 (注意第一个位置参数是 cases 目录, 不是 runs):",
          f'python scripts/solve_fgt_J.py "{os.path.join(args.outdir, "02_supercell", "cases")}" '
          f'--runs "{os.path.join(args.outdir, "runs")}" --spin {args.spin}',
          "# 自检 (不需要 DFT):",
          f'python scripts/solve_fgt_J.py "{os.path.join(args.outdir, "02_supercell", "cases")}" --demo',
          "```", "",
          "> 局限: 共线双亚晶格近似; 只到 J4; S=1 统一的等效自旋; MAE/各向异性未纳入;",
          "> 单层结构未弛豫 (实验体相 a=3.991 Å); k 网格 3x3x1 偏粗 (金属建议加密复核)。",
          "> **第 1 步别跑错 INCAR**: `runs/<order>/INCAR` 全是**静态**输入 (IBRION 默认 -1), 弛豫要用",
          "> `03_dft/INCAR.relax`; 跑错的话 VASP 照样输出 CONTCAR (坐标与 POSCAR 逐位相同), 回填脚本",
          "> 会报 `[warn] ... 最大坐标位移 0.0000 A` —— 见到这句就说明结构根本没弛豫、反解前提是废的。",
          "> **反解前必查 (第 3 步)**: 每个磁态 OUTCAR 里**最终**每个 Fe 的磁矩是否与 magmom.txt 的 ± 目标一致",
          "> (VASP 可能从一个初始 MAGMOM 弛豫到别的态, 尤其 afm_b_vert 是 2.55 A 强键合二聚体的反向态);",
          "> 不一致就换初始磁矩/换磁态, 否则代入系数矩阵的能量是错的。",
          "> 没有 OUTCAR 时 (只留了总能量) 可把总能量 (eV, 一个数) 写进 runs/<order>/E.txt。", ""]
    with open(os.path.join(args.root, "hamiltonian.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(h))
    print(f"哈密顿卡片 -> {os.path.join(args.root, 'hamiltonian.md')}")


if __name__ == "__main__":
    main()
