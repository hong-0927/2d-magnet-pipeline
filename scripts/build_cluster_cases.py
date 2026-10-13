#!/usr/bin/env python
"""build_cluster_cases.py
簇磁体（如 Nb3X8 的 Nb3 三聚体）的超胞 + 磁序枚举 + 能量映射系数。

模型（共线 Ising 映射，J>0=铁磁）:
    E = E0 + J_intra*C_intra + J1*C1 + J2*C2 + J3*C3 ...
    C_intra = -Σ_{簇内 Nb-Nb 键} σ_i σ_j                (微观簇内键, 分离 J_intra)
    C_k     = -Σ_{第 k 个簇间壳的簇对} m_a m_b,  m_a = 簇 a 净磁矩 = Σ_{i∈a} σ_i

簇间壳按**簇心晶格**距离分类（不是 Nb-Nb 键长）:
    1NN=7.7 A (J1), 2NN=13.33 A (J2), 3NN=15.4 A (J3), ...

关键点（见 references/cluster-magnets.md）:
  - 次近邻 (J2) / 次次近邻 (J3) **必须考虑**，不能只取最近邻 J1。本脚本用
    --nbr-shells 控制取到第几个簇间壳，默认 3（含 J1/J2/J3）。
  - 超胞尺寸要求: 解到第 k 个壳需 sc >= 2*smax_k + 1 (smax_k = 该壳最大平移分量)。
    1NN (smax=1) -> >=3x3; 2NN/3NN (smax=2) -> >=5x5。脚本会打印并警告过小的超胞。

用法:
  python build_cluster_cases.py POSCAR --element Nb --nonmag I \
      --cutoff 3.5 --supercell 5 5 1 --nbr-shells 3 --outdir cases \
      --U 3.0 --moment 1.0 --encut 520 --root <proj> [--selftest]
输出:
  cases/<name>/{POSCAR, magmom.txt}   磁序超胞
  cases/magnetic_orders.json         键壳 + 键和 + 反解矩阵 + 超胞尺寸自检
  runs/<name>/                        组装好的 VASP 提交目录 (INCAR/KPOINTS/POSCAR)
  --selftest: 前向模型合成能量 -> 反解, 验证映射良态
"""
import argparse
import json
import os
import shutil
from itertools import product

import numpy as np


# ---------- POSCAR / supercell ----------

def parse_poscar(path):
    lines = [ln for ln in open(path, encoding="utf-8").read().splitlines() if ln.strip()]
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


def write_poscar(path, lat, atoms, coords, comment="magnetic-order supercell"):
    order = []
    for a in atoms:
        if a not in order:
            order.append(a)
    counts = [atoms.count(sp) for sp in order]
    out = [comment, "1.00000000000000"]
    for row in lat:
        out.append(" ".join(f"{x:20.12f}" for x in row))
    out.append("  ".join(order))
    out.append("  ".join(str(c) for c in counts))
    out.append("Direct")
    for sp in order:
        for a, c in zip(atoms, coords):
            if a == sp:
                out.append(" ".join(f"{x:20.12f}" for x in c))
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


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


def min_image_frac(c1, c2):
    d = c2 - c1
    d -= np.floor(d + 0.5)
    return d


def min_image_distance(c1, c2, lat):
    return float(np.linalg.norm(min_image_frac(c1, c2) @ lat))


# ---------- cluster identification ----------

def find_clusters(atoms, coords, lat, element, cutoff):
    idxs = [i for i, a in enumerate(atoms) if a == element]
    n = len(idxs)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if min_image_distance(coords[idxs[i]], coords[idxs[j]], lat) < cutoff:
                ra, rb = find(i), find(j)
                if ra != rb:
                    parent[rb] = ra
    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(idxs[i])
    return list(groups.values())


# ---------- trimer-lattice shells (簇心晶格, 非 Nb-Nb 键长) ----------

def trimer_shells(lat_p, R=3, tol=0.5):
    """在**原胞晶格**上枚举整数平移 s, 按 |s@lat_p| 聚成簇间壳。

    返回 (centers, shell_shifts):
      centers[k]      = 第 k 壳的簇心距离 (A)
      shell_shifts[k] = 该壳的所有平移矢量 (整数分量)
    """
    by_dist = {}
    for s in product(range(-R, R + 1), range(-R, R + 1), [0]):
        if s[0] == 0 and s[1] == 0:
            continue
        d = float(np.linalg.norm(np.array(s, dtype=float) @ lat_p))
        by_dist.setdefault(round(d, 4), []).append(s)
    dists = sorted(by_dist)
    centers = []
    for d in dists:
        if not centers or d - centers[-1] > tol:
            centers.append(d)
    shell_shifts = [[] for _ in centers]
    for d in dists:
        k = int(np.argmin([abs(d - c) for c in centers]))
        shell_shifts[k].extend(by_dist[d])
    for k in range(len(shell_shifts)):
        shell_shifts[k] = sorted(set(shell_shifts[k]))
    return centers, shell_shifts


def shell_smax(shifts):
    """壳内平移的最大 |分量| -> 决定超胞最小尺寸。"""
    return max(max(abs(s[0]), abs(s[1])) for s in shifts)


# ---------- 磁序构型 ----------

def sigma_from_tau(tau, trimers, atoms_sc, coords_sc):
    """由每簇净自旋 tau[g] (S=1/2: +/-1) 构造微观 σ (Nb=+/-1, 非磁=0)。

    簇内固定 2up1down 构型: 排序后前 2 个 +1, 第 3 个 -1, 再乘 tau[g]。
    """
    sig = np.zeros(len(atoms_sc))
    for g, t in enumerate(trimers):
        ordered = sorted(t, key=lambda i: tuple(coords_sc[i]))
        for k, i in enumerate(ordered):
            sig[i] = tau[g] * (-1.0 if k == 2 else 1.0)
    return sig


def q_pattern(q, cell_coords):
    """按波矢 q (倒空间分数坐标) 生成每簇 tau = sign(cos(2π q·n))。"""
    tau = np.empty(len(cell_coords))
    for g, (n0, n1, _n2) in enumerate(cell_coords):
        ph = 2 * np.pi * (q[0] * n0 + q[1] * n1)
        c = np.cos(ph)
        tau[g] = 1.0 if c >= -1e-9 else -1.0
    return tau


# ---------- INCAR ----------

def build_incar_static(species, counts, element, U, moment, encut, npar=4):
    nm = counts[species.index(element)]
    n_nonmag = sum(counts) - nm
    ldaul = " ".join("2" if s == element else "-1" for s in species)
    ldauu = " ".join(f"{U:.1f}" if s == element else "0.0" for s in species)
    ldauj = " ".join("0.0" for _ in species)
    return f"""SYSTEM   = {" ".join(species)} monolayer cluster magnet order
PREC     = Normal
ENCUT    = {encut}
EDIFF    = 1E-6
ALGO     = Normal
NELM     = 200
ISMEAR   = 0
SIGMA    = 0.05

ISPIN    = 2
MAGMOM   = {nm}*{moment:.1f} {n_nonmag}*0.0     ! <-- 实际用 runs/<order>/INCAR

LDAU     = .TRUE.
LDAUTYPE = 2
LDAUL    = {ldaul}
LDAUU    = {ldauu}
LDAUJ    = {ldauj}
LMAXMIX  = 4

LREAL    = .FALSE.
LWAVE    = .FALSE.
LCHARG   = .FALSE.
NPAR     = {npar}
"""


KPOINTS_TMPL = """Gamma-centered k-mesh for monolayer supercell (vacuum along c)
3 3 1
0 0 0
"""


def build_potcar_spec(species):
    pa = {"Nb": "Nb_pv (含 4s4p 半芯)", "Cl": "Cl", "Br": "Br", "I": "I",
          "V": "V_pv", "Ti": "Ti_pv", "Ta": "Ta_pv"}
    lines = ["# POTCAR 组装说明", "",
             "元素顺序 **必须与 POSCAR/magmom.txt 一致**:  " + " ".join(species), ""]
    lines.append("| 顺序 | 元素 | 建议 PAW |")
    lines.append("|---|---|---|")
    for i, s in enumerate(species):
        lines.append(f"| {i + 1} | {s} | {pa.get(s, s)} |")
    lines.append("")
    lines.append("实际 PAW 版本/日期务必回填。")
    return "\n".join(lines) + "\n"


# ---------- main ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("poscar")
    ap.add_argument("--element", required=True, help="磁性元素, 如 Nb")
    ap.add_argument("--nonmag", nargs="+", default=[], help="非磁元素, 如 I")
    ap.add_argument("--cutoff", type=float, default=3.5, help="簇内成键阈值 (A)")
    ap.add_argument("--supercell", type=int, nargs=3, default=None,
                    help="超胞 (默认按 --nbr-shells 自动定最小尺寸)")
    ap.add_argument("--nbr-shells", type=int, default=3,
                    help="簇间壳数 (1=只 J1; 3=J1/J2/J3 含次近邻/次次近邻)")
    ap.add_argument("--outdir", default="cases")
    ap.add_argument("--root", default=None, help="项目根, 用于组装 runs/")
    ap.add_argument("--U", type=float, default=3.0)
    ap.add_argument("--moment", type=float, default=1.0)
    ap.add_argument("--encut", type=float, default=520)
    ap.add_argument("--selftest", action="store_true",
                    help="前向模型合成能量 -> 反解, 验证映射良态")
    args = ap.parse_args()

    lat, atoms, coords, mode = parse_poscar(args.poscar)
    assert mode[:1] == "d", f"POSCAR 需 Direct, 得到 {mode!r}"

    # 簇间壳 (原胞晶格)
    centers, shell_shifts = trimer_shells(lat)
    n_shells_avail = len(centers)
    if args.nbr_shells > n_shells_avail:
        raise SystemExit(f"--nbr-shells {args.nbr_shells} 超过可分辨壳数 {n_shells_avail}")
    inter_centers = centers[:args.nbr_shells]
    inter_shifts = shell_shifts[:args.nbr_shells]

    # 超胞尺寸: 默认按最大平移分量自动定, 过小则警告
    smax = max(shell_smax(sh) for sh in inter_shifts)
    req = 2 * smax + 1
    if args.supercell is None:
        sc = [req, req, 1]
    else:
        sc = list(args.supercell)
        if sc[0] < req or sc[1] < req:
            print(f"[WARN] 超胞 {tuple(sc)} 过小: 解到这些壳需 >= {req}x{req}x1 "
                  f"(smax={smax}), 平移会重叠导致键和**双重计数**。建议 {req} {req} 1。")

    lat_sc, atoms_sc, coords_sc, cell_idx = make_supercell(lat, atoms, coords, sc)
    mag_idx = [i for i, a in enumerate(atoms_sc) if a == args.element]
    n = len(mag_idx)

    groups = find_clusters(atoms_sc, coords_sc, lat_sc, args.element, args.cutoff)
    n_clusters = len(groups)

    # 分子簇模型的两个前提, 不满足就必须**报错退出**(原先只 warn, 一路走到 compute_coeffs
    # 才以 KeyError 崩掉 —— FGT 就是这样: 整个 5x5 超胞的 75 个 Fe 连成 1 个簇)。
    #   1) 超胞内 >= 2 个簇, 否则"簇间耦合"无从定义 (c2g 查不到周期像 -> KeyError);
    #   2) 每簇恰为 3 原子, 因为簇内系数 loop 写死 range(3) (Nb3I8 三聚体)。
    if n_clusters < 2:
        raise SystemExit(
            f"[cluster] 超胞内只有 {n_clusters} 个簇 (共 {n} 个 {args.element} 原子) -> "
            f"簇间耦合无从定义, 建不出「分子簇磁体」模型。\n"
            f"  多半是该体系**不是分子簇磁体, 而是连通的磁性网络**: 磁性近邻 < --cutoff "
            f"{args.cutoff} A 把整个超胞连成了 1 个簇 (典型: Fe3GeTe2 / Fe3GaTe2 这类多亚晶格"
            f"铁磁金属)。\n"
            f"  本流水线不支持多亚晶格模型; 需另建双亚晶格/有效自旋模型 "
            f"(见 SKILL.md「多亚晶格铁磁网络」)。")
    sizes = sorted({len(g) for g in groups})
    if sizes != [3]:
        raise SystemExit(
            f"[cluster] 簇大小 = {sizes}, 不是 [3]: 本脚本的簇内系数只对 **3 原子簇** 成立 "
            f"(Nb3I8 三聚体), 簇更大时多出的原子会被静默忽略、更小时直接越界。\n"
            f"  说明 --cutoff {args.cutoff} A 没圈出正确的分子单元, 或该体系本就不是三聚体簇磁体。")
    print(f"supercell {tuple(sc)} -> {len(atoms_sc)} atoms ({n} {args.element}), "
          f"{n_clusters} 簇, 每簇 {len(groups[0])} 原子")

    trimers = [sorted(g, key=lambda i: tuple(coords_sc[i])) for g in groups]
    # 每簇所属原胞 (n0,n1,n2), 用于簇心晶格 + 波矢构型
    cell_of = np.array([tuple(cell_idx[t[0]]) for t in trimers], dtype=int)
    c2g = {tuple(cell_of[g]): g for g in range(len(trimers))}

    # 打印簇间壳
    print("簇间壳 (簇心晶格):")
    for k, (c, sh) in enumerate(zip(inter_centers, inter_shifts)):
        print(f"  shell{k + 1} (J{k + 1}): {c:.3f} A, {len(sh)} 个平移, smax={shell_smax(sh)}")

    def compute_coeffs(sig):
        # 簇内
        intra_sum = 0.0
        for t in trimers:
            for a in range(3):
                for b in range(a + 1, 3):
                    intra_sum += sig[t[a]] * sig[t[b]]
        C_intra = -intra_sum
        m = np.array([sum(sig[i] for i in t) for t in trimers])
        Cs = []
        for sh in inter_shifts:
            direct = 0.0
            for g in range(len(trimers)):
                for s in sh:
                    tc = tuple((cell_of[g][d] + s[d]) % sc[d] for d in range(3))
                    h = c2g[tc]
                    direct += m[g] * m[h]
            Cs.append(-direct / 2.0)
        return C_intra, Cs

    # 构型: fm (全同, 分离 J_intra) + 一组波矢磁序 (S=1/2 流形, 分离 J1/J2/J3)
    configs = []

    sig_fm = np.zeros(len(atoms_sc))
    sig_fm[mag_idx] = 1.0
    ci, cs = compute_coeffs(sig_fm)
    configs.append({"name": "fm", "sig": sig_fm, "C_intra": ci, "Cs": cs,
                    "note": "全部 Nb 同向 (簇 S=3/2), 用于 J_intra"})

    q_cands = [
        (0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5),
        (1 / 3, 0.0), (0.0, 1 / 3), (1 / 3, 1 / 3), (1 / 3, 2 / 3),
        (2 / 3, 1 / 3), (2 / 3, 2 / 3),
        (1 / 4, 0.0), (0.0, 1 / 4), (1 / 4, 1 / 4), (1 / 4, 1 / 2), (1 / 2, 1 / 4),
    ]
    seen = set()
    for q in q_cands:
        tau = q_pattern(q, cell_of)
        key = tuple(int(round(x)) for x in tau)
        if key in seen:
            continue
        seen.add(key)
        sig = sigma_from_tau(tau, trimers, atoms_sc, coords_sc)
        ci, cs = compute_coeffs(sig)
        qname = "fm_trim" if q == (0.0, 0.0) else "q_" + "_".join(f"{v:g}" for v in q)
        configs.append({"name": qname, "sig": sig, "C_intra": ci, "Cs": cs,
                        "note": f"波矢 q=({q[0]:g},{q[1]:g}) S=1/2 磁序"})

    n_unknowns = 2 + args.nbr_shells  # E0 + J_intra + J1..J{nbr}
    n_select = min(n_unknowns + 2, len(configs))  # 多选 2 个 -> 超定最小二乘, 降噪

    def design(indices):
        return np.array([[1.0] + [configs[i]["C_intra"]] + configs[i]["Cs"]
                         for i in indices])

    def svd_score(indices):
        """优先满秩(rank), 再最大化最小非零奇异值(良态)。"""
        A = design(indices)
        s = np.linalg.svd(A, compute_uv=False)
        if s.size == 0:
            return 0.0
        tol = s.max() * max(A.shape) * np.finfo(float).eps * 100
        nz = s[s > tol]
        rank = int(nz.size)
        return rank * 1e12 + float(nz.min()) if rank else 0.0

    # 良态选取: 必含 fm(0) + fm_trim(1), 再贪心加满秩
    fm_trim_idx = next(i for i, c in enumerate(configs) if c["name"] == "fm_trim")
    selected = [0, fm_trim_idx]
    candidates = [i for i in range(len(configs)) if i not in selected]
    while len(selected) < n_select and candidates:
        best, bs = None, -1.0
        for c in candidates:
            sc_ = svd_score(selected + [c])
            if sc_ > bs:
                bs, best = sc_, c
        if bs <= 1e-10:
            break
        selected.append(best)
        candidates.remove(best)
    if len(selected) < n_unknowns:
        raise SystemExit(f"构型不足/退化: 只能选 {len(selected)}/{n_unknowns} 个独立构型, "
                         f"请增大超胞或增加波矢构型")

    A = design(selected)
    pinv = np.linalg.pinv(A)
    cond = np.linalg.cond(A)
    sv = np.linalg.svd(A, compute_uv=False)
    print(f"\n模型: E = E0 + J_intra*C_intra + " +
          " + ".join(f"J{k + 1}*C{k + 1}" for k in range(args.nbr_shells)))
    print(f"选取 {len(selected)} 个构型 (未知数 {n_unknowns}), 条件数 cond(A)={cond:.3g}, "
          f"奇异值 {[f'{x:.3g}' for x in sv]}")
    for i in selected:
        c = configs[i]
        print(f"  {c['name']:12s} C_intra={c['C_intra']:8.2f} "
              + "  ".join(f"C{k + 1}={c['Cs'][k]:8.2f}" for k in range(args.nbr_shells)))

    # 输出 cases
    os.makedirs(args.outdir, exist_ok=True)

    def grouped_order(atoms_list):
        order = []
        for a in atoms_list:
            if a not in order:
                order.append(a)
        return order, [[i for i, a in enumerate(atoms_list) if a == sp] for sp in order]

    order, groups_atom = grouped_order(atoms_sc)

    for i in selected:
        c = configs[i]
        d = os.path.join(args.outdir, c["name"])
        os.makedirs(d, exist_ok=True)
        write_poscar(os.path.join(d, "POSCAR"), lat_sc, atoms_sc, coords_sc,
                     comment=f"{args.element}3 cluster {c['name']}")
        mag_parts = []
        for sp, gidx in zip(order, groups_atom):
            if sp == args.element:
                mag_parts += [f"{args.moment * c['sig'][j]:.1f}" for j in gidx]
            else:
                mag_parts += ["0.0"] * len(gidx)
        with open(os.path.join(d, "magmom.txt"), "w") as f:
            f.write(f"MAGMOM = {' '.join(mag_parts)}\n")

    report = {
        "input": args.poscar, "supercell": sc, "n_atoms": len(atoms_sc),
        "n_mag": n, "n_clusters": n_clusters, "cluster_size": len(groups[0]),
        "element": args.element, "nonmag": args.nonmag,
        "moment": args.moment, "U": args.U, "nbr_shells": args.nbr_shells,
        "model": "E = E0 + J_intra*C_intra + sum_k Jk*Ck (J>0=铁磁)",
        "inter_shells_angstrom": inter_centers,
        "inter_shell_shifts": [[list(s) for s in sh] for sh in inter_shifts],
        "required_supercell": [req, req, 1],
        "supercell_ok": (sc[0] >= req and sc[1] >= req),
        "selected_configs": [configs[i]["name"] for i in selected],
        "coefficients": {configs[i]["name"]: {"C_intra": configs[i]["C_intra"],
                                              "C_inter": configs[i]["Cs"]}
                         for i in selected},
        "pseudo_inverse_E0_Jintra_Jk": pinv.tolist(),
    }
    with open(os.path.join(args.outdir, "magnetic_orders.json"), "w") as f:
        json.dump(report, f, indent=2)

    # 组装 runs/ + INCAR (若给了 --root)
    if args.root:
        root = os.path.abspath(args.root)
        sp, counts = parse_poscar_counts(os.path.join(args.outdir, configs[selected[0]]["name"], "POSCAR"))
        dft_dir = os.path.join(root, "03_dft")
        os.makedirs(dft_dir, exist_ok=True)
        static_tmpl = build_incar_static(sp, counts, args.element, args.U, args.moment, args.encut)
        with open(os.path.join(dft_dir, "INCAR.static"), "w") as f:
            f.write(static_tmpl)
        with open(os.path.join(dft_dir, "KPOINTS"), "w") as f:
            f.write(KPOINTS_TMPL)
        with open(os.path.join(dft_dir, "POTCAR.spec"), "w") as f:
            f.write(build_potcar_spec(sp))
        nm = counts[sp.index(args.element)]
        n_nonmag = sum(counts) - nm
        default_line = f"MAGMOM   = {nm}*{args.moment:.1f} {n_nonmag}*0.0     ! <-- 实际用 runs/<order>/INCAR"
        runs_dir = os.path.join(root, "runs")
        os.makedirs(runs_dir, exist_ok=True)
        for i in selected:
            c = configs[i]
            d = os.path.join(args.outdir, c["name"])
            if not os.path.isdir(d):
                continue
            magmom_line = open(os.path.join(d, "magmom.txt"), encoding="utf-8").read().strip()
            incar = static_tmpl.replace(default_line, magmom_line)
            rd = os.path.join(runs_dir, c["name"])
            os.makedirs(rd, exist_ok=True)
            with open(os.path.join(rd, "INCAR"), "w") as f:
                f.write(incar)
            shutil.copy(os.path.join(dft_dir, "KPOINTS"), os.path.join(rd, "KPOINTS"))
            shutil.copy(os.path.join(d, "POSCAR"), os.path.join(rd, "POSCAR"))
        print("\n写出 03_dft/ (INCAR.static, KPOINTS, POTCAR.spec) 与 runs/<order>/")
        print("POTCAR 顺序:", " ".join(sp))

    # 自检: 前向模型合成能量 -> 反解, 验证恢复 J
    if args.selftest:
        J_true = [-30.0] + [-2.0 - 0.3 * k for k in range(args.nbr_shells)]  # J_intra + J1..Jk
        E0_true = 0.0
        Evec = np.array([E0_true + configs[i]["C_intra"] * J_true[0]
                         + sum(configs[i]["Cs"][k] * J_true[1 + k]
                               for k in range(args.nbr_shells))
                         for i in selected])
        sol = pinv @ Evec
        ok = np.allclose(sol, [E0_true] + J_true, atol=1e-8)
        print(f"\n--selftest: 合成 J_intra={J_true[0]}, "
              f"J={J_true[1:]} -> 反解 {[f'{x:.3f}' for x in sol]}  "
              f"{'PASS' if ok else 'FAIL'}")
        if not ok:
            raise SystemExit("自检失败: 映射矩阵良态性有问题")


def parse_poscar_counts(path):
    lines = [ln for ln in open(path, encoding="utf-8").read().splitlines() if ln.strip()]
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    return species, counts


if __name__ == "__main__":
    main()
