#!/usr/bin/env python
"""magnetic_orders.py — 磁序规格层 (P0.1)

把 5 条既有链路 (三角 / 矩形 / 蜂窝 / 簇 / FGT) 里各自复制的
「超胞展开 + 磁序构造 + 键和系数」逻辑收敛到**一个**模块。
新晶格 (方形 / kagome / oxyhalide / 4f) 只需往 REGISTRY 里加**规格数据**,
不再新写 build_*_cases.py。

核心不变式 (与 5 条旧路径数值一致)
----------------------------------
位点-位点模型::

    E = e0 + Σ_k c_k · J~_k ,   J~ = J·S²
    c_k = -Σ_{第 k 壳的无向键 i<j} (s_i · s_j)

旧代码里出现的「有向计数 × (-0.5)」写法与之**恒等**: 每根无向键在有向和里
出现两次且 σ_iσ_j = σ_jσ_i, 故 -0.5·Σ_{有向} ≡ -Σ_{i<j}。
(build_magnetic_cases_tri.py:118 / build_rectangular_cases.py:14 /
 build_magnetic_cases3.py:219 / build_fgt_cases.py:19)

簇磁体 (Nb3I8) 额外一层「簇净矩」粗粒化::

    C_intra = -Σ_{簇内键} σ_i σ_j
    C_k     = -0.5 · Σ_{第 k 簇间壳的簇对 (a,b)} (m_a m_b) ,  m_a = Σ_{i∈a} σ_i

三个可复用部件
--------------
StructureSpec : POSCAR + 元素 + 超胞        -> 几何 (lat/atoms/frac/cell_idx)
Shell         : 壳中心距离 + 容差 + (可选)亚晶格对过滤
Order         : 磁序 = 一个 (cell_idx, 亚晶格, 簇) -> 自旋矢量的纯函数

用法
----
  python scripts/magnetic_orders.py --list
  python scripts/magnetic_orders.py --spec tri_mni2
  python scripts/magnetic_orders.py --regress          # 逐位比对 5 条旧链路的 JSON
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from itertools import product
from typing import Callable

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))  # scripts -> 2d_magnet_pipeline -> 项目根


def resolve(path: str) -> str:
    return path if os.path.isabs(path) else os.path.join(ROOT, path)


# =====================================================================
# 1. 几何: POSCAR -> 超胞
# =====================================================================

def parse_poscar(path: str):
    """读 POSCAR, 返回 (lat, atoms, frac, mode)。frac 一律归一化到 [0,1) 分数坐标。"""
    with open(path, encoding="utf-8") as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()]
    scale = float(lines[1].split()[0])
    lat = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)]) * scale
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    idx = 7
    if lines[idx].strip().lower()[:1] == "s":  # Selective dynamics
        idx += 1
    mode = lines[idx].strip()
    idx += 1
    frac, atoms = [], []
    for sp, n in zip(species, counts):
        for _ in range(n):
            frac.append([float(x) for x in lines[idx].split()[:3]])
            atoms.append(sp)
            idx += 1
    frac = np.asarray(frac, dtype=float)
    if mode[:1].lower() in ("c", "k"):  # Cartesian / Kpoints-line
        frac = frac @ np.linalg.inv(lat)
    frac -= np.floor(frac)
    return lat, atoms, frac, mode


def make_supercell(lat, atoms, frac, sc):
    """显式超胞 (与旧链路同约定: 原胞索引外循环, 原子内循环)。

    返回 (lat_sc, atoms_sc, frac_sc, cell_idx, prim_idx)。
    """
    sc = np.asarray(sc, dtype=int)
    out_a, out_f, cell, prim = [], [], [], []
    for n0 in range(int(sc[0])):
        for n1 in range(int(sc[1])):
            for n2 in range(int(sc[2])):
                shift = np.array([n0, n1, n2], dtype=float)
                for ia, (a, f) in enumerate(zip(atoms, frac)):
                    out_a.append(a)
                    out_f.append((f + shift) / sc)
                    cell.append((n0, n1, n2))
                    prim.append(ia)
    lat_sc = np.asarray(lat, dtype=float).copy()
    for d in range(3):
        lat_sc[d] = lat_sc[d] * sc[d]
    return (lat_sc, out_a, np.asarray(out_f, dtype=float),
            np.asarray(cell, dtype=int), np.asarray(prim, dtype=int))


def min_image_frac(c1, c2):
    d = np.asarray(c2, dtype=float) - np.asarray(c1, dtype=float)
    return d - np.floor(d + 0.5)


def min_image_distance(c1, c2, lat):
    return float(np.linalg.norm(min_image_frac(c1, c2) @ lat))


def magnetic_shells(idxs, frac, lat, max_shells=3, tol=0.2):
    """最小镜像距离聚类, 自动发现前 max_shells 个磁距壳 (供新晶格用)。

    注意: 该启发式与 build_rectangular_cases.py 的显式平移搜索**不等价**,
    正式规格一律**显式声明**壳距离 (见 REGISTRY), 本函数只用于快速探路。
    """
    ds = []
    n = len(idxs)
    for i in range(n):
        for j in range(i + 1, n):
            dist = min_image_distance(frac[idxs[i]], frac[idxs[j]], lat)
            if dist > 1e-6:
                ds.append(dist)
    ds.sort()
    shells: list[float] = []
    for d in ds:
        if not shells or d - shells[-1] > tol:
            shells.append(d)
    return shells[:max_shells]


def cluster_center_shells(lat_prim, R=3, tol=0.5):
    """在**原胞晶格**上枚举整数平移, 按 |s·lat| 聚成簇间壳。

    返回 (centers, shift_sets): 壳心距离与对应的整数平移矢量集。
    与 build_cluster_cases.py:141 trimer_shells 同算法。
    """
    by_dist: dict[float, list] = {}
    for s in product(range(-R, R + 1), range(-R, R + 1), [0]):
        if s[0] == 0 and s[1] == 0:
            continue
        d = float(np.linalg.norm(np.asarray(s, dtype=float) @ lat_prim))
        by_dist.setdefault(round(d, 4), []).append(s)
    dists = sorted(by_dist)
    centers: list[float] = []
    for d in dists:
        if not centers or d - centers[-1] > tol:
            centers.append(d)
    shift_sets: list[list] = [[] for _ in centers]
    for d in dists:
        k = int(np.argmin([abs(d - c) for c in centers]))
        shift_sets[k].extend(by_dist[d])
    return centers, [sorted(set(s)) for s in shift_sets]


# =====================================================================
# 2. 规格数据
# =====================================================================

@dataclass(frozen=True)
class Shell:
    """一个距离壳。pair 为 None 表示不按亚晶格过滤。"""
    name: str
    d: float
    tol: float = 0.3
    pair: str | None = None       # "AB" / "AA" / "BB" / "A0A1" ...


@dataclass(frozen=True)
class StructureSpec:
    poscar: str
    element: str
    supercell: tuple = (1, 1, 1)
    nonmag: tuple = ()
    # 按原胞分数坐标**排序**后的亚晶格命名; 默认 A, B, C...
    sublattice_names: tuple | None = None
    # 每个亚晶格所属「族」(用于 Shell.pair 匹配); 默认 == names
    sublattice_family: tuple | None = None


@dataclass(frozen=True)
class ClusterModel:
    """簇净矩模式 (Nb3I8 型三聚体): 簇间壳按**簇心晶格**定义。"""
    cutoff: float                  # 簇内连通判据 (Nb-Nb 距离, A)
    intra: Shell
    inter: tuple = ()
    R: int = 3
    tol: float = 0.5


@dataclass(frozen=True)
class Order:
    """一个磁序构型: build(ctx) -> (n,) 共线标量 或 (n,3) 自旋矢量。"""
    name: str
    build: Callable


@dataclass(frozen=True)
class Spec:
    key: str
    label: str
    structure: StructureSpec
    orders: tuple
    shells: tuple = ()
    spin: float = 1.0
    policy: str = "min_image"             # 键枚举策略: "min_image" | "shift_search"
    nshift: int = 1
    cluster: ClusterModel | None = None
    post: Callable | None = None          # 派生构型 (如蜂窝 zigzag/stripy)
    fixture: str | None = None
    fixture_kind: str = "coefficients_per_cell"
    selected: tuple | None = None         # 旧链路实际选用的构型 (信息用)


# =====================================================================
# 3. 磁序构造小工具 (供 Order.build 使用)
# =====================================================================

def collinear(sigma):
    sigma = np.asarray(sigma, dtype=float)
    s = np.zeros((len(sigma), 3))
    s[:, 2] = sigma
    return s


def sign_pow(base):
    """(-1)^base, base 为整数数组 (返回 ±1)。"""
    return np.where(np.asarray(base) % 2 == 0, 1.0, -1.0)


def planar120(cell, period=3):
    """三角/蜂窝晶格 120° 共面非共线序: θ = (n0+n1) mod period · 120°。

    注意: 相位只依赖**晶胞索引**, 对同一胞内的多个磁性亚晶格给同一自旋。
    因此它**不能**用于 kagome (3 磁位/胞): 那样会退化成铁磁。
    kagome 用 `planar120_sub` (q=0) 与 `planar120_sqrt3` (√3×√3)。
    """
    color = (np.asarray(cell)[:, 0] + np.asarray(cell)[:, 1]) % period
    theta = color * (2.0 * np.pi / period)
    s = np.zeros((len(cell), 3))
    s[:, 0] = np.cos(theta)
    s[:, 1] = np.sin(theta)
    return s


def planar120_sub(sub_label, period=3):
    """kagome q=0 120° 共面非共线序: θ = sub_label · 120° (相位按**亚晶格**)。

    同一胞内 A/B/C 三个磁位分别取 0°/120°/240°, 所有胞等价 (q=0)。
    """
    sub_label = np.asarray(sub_label)
    theta = sub_label * (2.0 * np.pi / period)
    s = np.zeros((len(sub_label), 3))
    s[:, 0] = np.cos(theta)
    s[:, 1] = np.sin(theta)
    return s


def sqrt3_label(n0, n1, period=3):
    """kagome √3×√3 超胞的 Z3 标签 ℓ = (n0+n1) mod 3。

    在 √3×√3 的平移矢量 a1-a2 与 a1+2a2 下不变 (二者都把 n0+n1 改变 3 的倍数),
    故 ℓ 在商群 Z²/L(√3) ≅ Z3 上良定义 —— 这正是打破 q=0 简并所需的调制。
    """
    return (np.asarray(n0) + np.asarray(n1)) % period


def planar120_sqrt3(sub_label, n0, n1, period=3):
    """kagome √3×√3 120° 共面非共线序: θ = (sub_label + ℓ) · 120°。"""
    label = np.asarray(sub_label) + sqrt3_label(n0, n1, period)
    return planar120_sub(label, period)


def sublattice_signs(sub_label, signs):
    """按亚晶格标签取 ±1 (kagome 的 2up1down / 1up2down 共线序)。"""
    return np.asarray(signs, dtype=float)[np.asarray(sub_label, dtype=int)]


def sublattice_colors(mag_idx, frac, lat, d_nn):
    """按最近邻二分着色 (BFS), 返回 0/1 标签数组 (蜂窝 neel 用)。"""
    from collections import deque
    n = len(mag_idx)
    color = np.full(n, -1, dtype=int)
    if n == 0:
        return color
    color[0] = 0
    q = deque([0])
    while q:
        i = q.popleft()
        for j in range(n):
            if i == j:
                continue
            if min_image_distance(frac[mag_idx[i]], frac[mag_idx[j]], lat) < d_nn + 0.5:
                if color[j] == -1:
                    color[j] = 1 - color[i]
                    q.append(j)
    return color


def q_pattern(q, cell):
    """波矢型 tau[g] = sign(cos 2π(q·n)) (簇净矩序)。"""
    cell = np.asarray(cell)
    ph = 2.0 * np.pi * (q[0] * cell[:, 0] + q[1] * cell[:, 1])
    return np.where(np.cos(ph) >= -1e-9, 1.0, -1.0)


def cluster_2up1down(ctx, tau):
    """簇内固定 2up1down (按坐标排序后前 2 个 +1, 第 3 个 -1), 再乘每簇 tau。"""
    sig = np.zeros(ctx.n)
    pos = {int(a): p for p, a in enumerate(ctx.mag_idx)}
    for g, grp in enumerate(ctx.groups):
        ordered = sorted(grp, key=lambda i: tuple(ctx.frac_sc[i]))
        for k, i in enumerate(ordered):
            sig[pos[int(i)]] = tau[g] * (-1.0 if k == 2 else 1.0)
    return sig


# =====================================================================
# 4. 上下文
# =====================================================================

@dataclass
class OrderContext:
    lat_prim: np.ndarray
    lat_sc: np.ndarray
    atoms_sc: list
    frac_sc: np.ndarray
    cell_idx: np.ndarray
    prim_idx: np.ndarray
    mag_idx: np.ndarray
    n: int
    sub_label: np.ndarray        # (n,) 亚晶格序号
    sub_name: np.ndarray         # (n,) 亚晶格名 (str)
    sub_family: np.ndarray       # (n,) 亚晶格族 (str)
    supercell: tuple
    d_nn: float | None = None          # 最近邻壳距离 (蜂窝二分着色用)
    groups: list | None = None
    group_cell: np.ndarray | None = None

    @property
    def mcell(self):
        return self.cell_idx[self.mag_idx]

    @property
    def n0(self):
        return self.mcell[:, 0]

    @property
    def n1(self):
        return self.mcell[:, 1]

    @property
    def n_groups(self):
        return 0 if self.groups is None else len(self.groups)


def _sublattice_keys(atoms, frac, element):
    idxs = [i for i, a in enumerate(atoms) if a == element]
    coords = sorted({tuple(np.round(frac[i], 6)) for i in idxs})
    return idxs, {c: k for k, c in enumerate(coords)}, coords


def _find_clusters(atoms_sc, frac_sc, lat_sc, element, cutoff):
    """并查集: 距离 < cutoff 的同类原子连成一簇 (build_cluster_cases.py:116)。"""
    idxs = [i for i, a in enumerate(atoms_sc) if a == element]
    n = len(idxs)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n):
        for j in range(i + 1, n):
            if min_image_distance(frac_sc[idxs[i]], frac_sc[idxs[j]], lat_sc) < cutoff:
                ra, rb = find(i), find(j)
                if ra != rb:
                    parent[rb] = ra
    groups: dict[int, list] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(idxs[i])
    return [sorted(v) for _, v in sorted(groups.items())]


def build_context(spec: Spec) -> OrderContext:
    st = spec.structure
    lat, atoms, frac, _mode = parse_poscar(resolve(st.poscar))
    lat_sc, atoms_sc, frac_sc, cell_idx, prim_idx = make_supercell(lat, atoms, frac, st.supercell)
    mag_idx = np.array([i for i, a in enumerate(atoms_sc) if a == st.element], dtype=int)
    _idxs, mapping, coord_keys = _sublattice_keys(atoms, frac, st.element)

    names = st.sublattice_names
    if names is None:
        names = tuple(chr(ord("A") + k) for k in range(len(coord_keys)))
    assert len(names) == len(coord_keys), (
        f"{spec.key}: sublattice_names 给了 {len(names)} 个, 但原胞里有 {len(coord_keys)} 个磁性亚晶格")
    family = st.sublattice_family or names
    assert len(family) == len(names), f"{spec.key}: sublattice_family 长度不符"

    sub_label = np.array([mapping[tuple(np.round(frac[prim_idx[i]], 6))] for i in mag_idx], dtype=int)
    ctx = OrderContext(
        lat_prim=lat, lat_sc=lat_sc, atoms_sc=atoms_sc, frac_sc=frac_sc,
        cell_idx=cell_idx, prim_idx=prim_idx, mag_idx=mag_idx, n=len(mag_idx),
        sub_label=sub_label,
        sub_name=np.array([names[k] for k in sub_label]),
        sub_family=np.array([family[k] for k in sub_label]),
        supercell=tuple(int(x) for x in st.supercell),
        d_nn=(spec.shells[0].d if spec.shells else None),
    )
    if spec.cluster is not None:
        ctx.groups = _find_clusters(atoms_sc, frac_sc, lat_sc, st.element, spec.cluster.cutoff)
        ctx.group_cell = np.array(
            [cell_idx[g[0]] for g in ctx.groups], dtype=int)
    return ctx


# =====================================================================
# 5. 系数引擎
# =====================================================================

def _pair_ok(pair: str | None, fa: str, fb: str) -> bool:
    if pair is None:
        return True
    return sorted(pair) == sorted((fa, fb))


def _assign_shell(dist: float, shells, fa: str, fb: str):
    """先按亚晶格对过滤再取最近壳 (否则同距离的 AA/BB 壳会互相遮蔽)。"""
    best, bd = None, float("inf")
    for k, s in enumerate(shells):
        if not _pair_ok(s.pair, fa, fb):
            continue
        dd = abs(dist - s.d)
        if dd < bd:
            bd, best = dd, k
    if best is not None and bd < shells[best].tol:
        return best
    return None


def site_pair_coeffs(ctx: OrderContext, spins, shells, policy="min_image", nshift=1) -> np.ndarray:
    """c_k = -Σ_{第 k 壳无向键 i<j} (s_i·s_j)。

    policy="min_image"   : 无向对 i<j + 最小镜像距离 (三角 / FGT)
    policy="shift_search": 有向对 (i,j) + 显式平移集, 整体乘 -0.5 (矩形 / 蜂窝)
                           两种写法**数学恒等**, 但小超胞下最小镜像会把远壳别名坍缩到
                           近处而被容差滤掉, 显式平移搜索不会 —— 故旧链路各异, 这里保留为参数。
    """
    n = ctx.n
    if policy == "min_image":
        c = np.zeros(len(shells))
        for i in range(n):
            ci = ctx.frac_sc[ctx.mag_idx[i]]
            fi = ctx.sub_family[i]
            for j in range(i + 1, n):
                dist = min_image_distance(ci, ctx.frac_sc[ctx.mag_idx[j]], ctx.lat_sc)
                if dist < 1e-6:
                    continue
                k = _assign_shell(dist, shells, fi, ctx.sub_family[j])
                if k is None:
                    continue
                c[k] -= float(spins[i] @ spins[j])
        return c

    if policy == "shift_search":
        c = np.zeros(len(shells))
        shifts = list(product(range(-nshift, nshift + 1), repeat=3))
        for i in range(n):
            ci = ctx.frac_sc[ctx.mag_idx[i]]
            fi = ctx.sub_family[i]
            for j in range(n):
                base = ctx.frac_sc[ctx.mag_idx[j]] - ci
                fj = ctx.sub_family[j]
                for sh in shifts:
                    dist = float(np.linalg.norm((base + np.asarray(sh, dtype=float)) @ ctx.lat_sc))
                    if dist < 1e-6:
                        continue
                    k = _assign_shell(dist, shells, fi, fj)
                    if k is None:
                        continue
                    c[k] += float(spins[i] @ spins[j])
        return -0.5 * c

    raise ValueError(f"未知键枚举策略 policy={policy!r}")


def cluster_coeffs(ctx: OrderContext, cm: ClusterModel, sigma) -> np.ndarray:
    """返回 [C_intra, C_1, ..., C_k]。sigma 为 (n,) 共线标量 (±1)。"""
    pos = {int(a): p for p, a in enumerate(ctx.mag_idx)}
    intra = 0.0
    for grp in ctx.groups:
        for x in range(len(grp)):
            for y in range(x + 1, len(grp)):
                intra += sigma[pos[int(grp[x])]] * sigma[pos[int(grp[y])]]
    m = np.array([sum(sigma[pos[int(i)]] for i in grp) for grp in ctx.groups])

    centers, shift_sets = cluster_center_shells(ctx.lat_prim, cm.R, cm.tol)
    nsh = len(cm.inter)
    if len(centers) < nsh:
        raise ValueError(f"簇心壳只有 {len(centers)} 个, 规格要 {nsh} 个; 请加大 R")
    for k in range(nsh):
        if abs(centers[k] - cm.inter[k].d) > 0.05:
            raise ValueError(f"簇间壳 {cm.inter[k].name}: 实测 {centers[k]:.4f} A "
                             f"vs 规格 {cm.inter[k].d:.4f} A")
    sc = ctx.supercell
    c2g = {tuple(int(v) for v in ctx.group_cell[g]): g for g in range(len(ctx.groups))}
    out = [float(-intra)]
    for k in range(nsh):
        # 超胞周期性: 同一目标胞的平移只算一次 (s 与 s-sc 等价)
        seen, uniq = set(), []
        for s in shift_sets[k]:
            key = tuple(int(s[d]) % sc[d] for d in range(3))
            if key not in seen:
                seen.add(key)
                uniq.append(s)
        direct = 0.0
        for g in range(len(ctx.groups)):
            for s in uniq:
                tc = tuple((int(ctx.group_cell[g][d]) + int(s[d])) % sc[d] for d in range(3))
                h = c2g[tc]
                direct += m[g] * m[h]
        out.append(float(-direct / 2.0))
    return np.array(out)


# =====================================================================
# 6. 计算 / 设计矩阵 / 条件数闸门
# =====================================================================

@dataclass
class OrderResult:
    name: str
    spins: np.ndarray
    coeffs: np.ndarray


def compute(spec: Spec) -> tuple[dict, OrderContext]:
    ctx = build_context(spec)
    results: dict[str, OrderResult] = {}
    for od in spec.orders:
        spins = np.asarray(od.build(ctx), dtype=float)
        if spins.ndim == 1:
            spins = collinear(spins)
        if spec.cluster is not None:
            coef = cluster_coeffs(ctx, spec.cluster, spins[:, 2])
        else:
            coef = site_pair_coeffs(ctx, spins, spec.shells, spec.policy, spec.nshift)
        results[od.name] = OrderResult(od.name, spins, coef)
    if spec.post is not None:
        results.update(spec.post(results))
    return results, ctx


def design_matrix(names, coeffs) -> np.ndarray:
    """A = [[1, c_1, ..., c_k]] —— 每行一个磁序, 首列定 E0。"""
    return np.array([[1.0] + [float(x) for x in coeffs[nm]] for nm in names])


def conditioning(A) -> float:
    return float(np.linalg.cond(np.asarray(A, dtype=float)))


def check_design(A, cond_max: float = 1e4, label: str = "") -> float:
    """条件数闸门: 病态/退化就 raise, **绝不静默出 J**。"""
    cd = conditioning(A)
    if not np.isfinite(cd) or cd > cond_max:
        raise ValueError(f"{label} 设计矩阵病态: cond(A)={cd:.4g} > {cond_max:.4g}; 拒绝反解")
    return cd


# =====================================================================
# 7. 注册表 (5 条既有链路 + 新晶格从这里长出来)
# =====================================================================

_Q_CANDS = [
    (0.0, 0.0), (0.5, 0.0), (0.0, 0.5), (0.5, 0.5),
    (1 / 3, 0.0), (0.0, 1 / 3), (1 / 3, 1 / 3), (1 / 3, 2 / 3),
    (2 / 3, 1 / 3), (2 / 3, 2 / 3),
    (0.25, 0.0), (0.0, 0.25), (0.25, 0.25), (0.25, 0.5), (0.5, 0.25),
]


def _cluster_qname(q) -> str:
    if q == (0.0, 0.0):
        return "fm_trim"
    return "q_" + "_".join(f"{v:g}" for v in q)


def _cluster_orders() -> tuple:
    """簇净矩序: fm (全 Nb 同向, 定 J_intra) + 一组波矢序 (S=1/2 流形, 定 J1/J2/J3)。

    tau 按**簇心晶格**的波矢生成, 簇内固定 2up1down。
    """
    orders = [Order("fm", lambda c: np.ones(c.n))]
    for q in _Q_CANDS:
        orders.append(Order(
            _cluster_qname(q),
            (lambda qq: (lambda c: cluster_2up1down(c, q_pattern(qq, c.group_cell))))(q),
        ))
    return tuple(orders)


def _honeycomb_zigzag_stripy(results):
    """zigzag = 三个 M 点序里 J1 系数**最大**的; stripy = 与之 J1 同值但 J3 不同者。

    复刻 build_magnetic_cases3.py:198-211 的判据 —— 那里比较的是 classify_bonds 返回的
    **有向键和** sums, 而 coeff = -0.5·sums, 故 "min(sums[0])" 等价于 "max(coeff[0])"。
    """
    m = {k: np.asarray(results[k].coeffs, dtype=float) for k in ("M0", "M1", "M01")}
    zig = max(m, key=lambda k: m[k][0])
    rest = [k for k in m if k != zig]
    strip = next((k for k in rest
                  if abs(m[k][0] - m[zig][0]) < 1e-9 and abs(m[k][2] - m[zig][2]) > 1e-9),
                 rest[0])
    return {"zigzag": results[zig], "stripy": results[strip]}


def _basis_sign(c):
    """C 心网 (2 个磁性原子/胞): 同胞两个磁性原子取相反符号 (basis0 -> +1, basis1 -> -1)。

    这是 C 心网的 Néel 态, 在矩形/方形/蜂窝的规格里都不存在 —— 那些晶格每胞只有
    一个磁性位点, 磁序完全由 (n0, n1) 决定, 表达不出「亚晶格反平行」。

    ⚠ 关键结构事实: C 心矢量 (1/2,1/2) 连接的是**同一胞内**的两个 basis, 所以这个
      符号对一切只依赖 (n0,n1) 奇偶的序**完全不可见** ——
          basis*(-1)^n0 == (-1)^n0,  basis*(-1)^(n0+n1) == (-1)^(n0+n1)
      也就是说 basis 翻倍版本与其母序逐位相同, 二者**不能同时选入设计矩阵**
      (会奇异)。只有它自己 (afm_c) 才带来新信息。
    """
    pid = c.prim_idx[c.mag_idx]
    return np.where(pid == pid.min(), 1.0, -1.0)


REGISTRY: dict[str, Spec] = {s.key: s for s in (

    # ---------------- 三角 (MnI2): 1 磁位/胞, 60° 晶格; 默认全共线四态 ----------------
    # 120° 序仍保留在 orders 中 (需要非共线对照/文献复现时显式启用), 但**不再是反解必需**:
    # fm + stripy_a + stripy_ab + uud 已给出 rank=4 (cond≈3.7e3)。
    Spec(
        key="tri_mni2",
        label="三角单离子 (MnI2)",
        structure=StructureSpec(
            poscar="MnI2_pipeline/01_structure/POSCAR", element="Mn",
            nonmag=("I",), supercell=(5, 5, 1)),
        shells=(
            Shell("J1", 4.168580934901683),
            Shell("J2", 7.220193659999998),
            Shell("J3", 8.337161869803365),
        ),
        orders=(
            Order("fm", lambda c: np.ones(c.n)),
            Order("120deg", lambda c: planar120(c.mcell, 3)),
            Order("stripy_a", lambda c: sign_pow(c.n0)),
            Order("stripy_b", lambda c: sign_pow(c.n1)),
            Order("stripy_ab", lambda c: sign_pow(c.n0 + c.n1)),
            Order("uud", lambda c: np.where((c.n0 + c.n1) % 3 == 0, -1.0, 1.0)),
        ),
        selected=("fm", "stripy_a", "stripy_ab", "uud"),
        spin=2.5,
        fixture="MnI2_pipeline/02_supercell/cases/magnetic_orders.json",
        fixture_kind="coefficients_per_cell",
    ),

    # ---------------- 矩形 (CrSBr): a≠b, 非受挫, 全共线序 ----
    Spec(
        key="rect_crsbr",
        label="矩形单离子 (CrSBr)",
        structure=StructureSpec(
            poscar="CrSBr_pipeline_v2/01_structure/POSCAR", element="Cr",
            nonmag=("S", "Br"), supercell=(4, 4, 1)),
        shells=(
            Shell("J1", 3.5404561399999994),
            Shell("J2", 4.75547185),
            Shell("J3", 5.928688066967777),
        ),
        orders=(
            Order("fm", lambda c: np.ones(c.n)),
            Order("afm_a", lambda c: sign_pow(c.n0)),
            Order("afm_b", lambda c: sign_pow(c.n1)),
            Order("afm_diag", lambda c: sign_pow(c.n0 + c.n1)),
        ),
        selected=("fm", "afm_a", "afm_b", "afm_diag"),
        spin=1.5,
        policy="shift_search",
        fixture="CrSBr_pipeline_v2/02_supercell/cases/magnetic_orders.json",
        fixture_kind="coefficients_per_cell",
    ),

    # ---------------- 方形 (K2NiF4): a=b, 1 磁位/胞, 非受挫; 需要周期 4 的序 ----
    # 关键 1: 一切周期 2 的序 (fm/g_type/stripe) 给出的 J3 系数**完全相同**
    #         (C(2,0)=N^2 恒成立) ⇒ 不引入周期 4 的序, J3 与 E0 无法分离。
    #         这与 kagome 的 q=0 简并同构, 只是成因不同。见 _probe_square 记录。
    # 关键 2: J3 = 2a 恰好是 4x4 超胞的半周期, 最小镜像会把 ±2a 两个像折叠成
    #         **同一个位点**, 于是每胞只数到 1 根 2a 键而正确值是 2 根 (少算一半)
    #         ⇒ 必须 policy="shift_search" 保留像的重数 (与矩形/蜂窝同样的理由)。
    #         证据: 4x4x1 下 fm 的 c3: min_image=-16 (错), shift_search=-32 (对,
    #         与无限晶格每胞 2 根 2a 键一致); J1=a 与 J2=a√2 分数位移 0.25/0.3535,
    #         不到半周期, 两种策略一致。
    Spec(
        key="sq_k2nif4",
        label="方形单离子 (K2NiF4)",
        structure=StructureSpec(
            poscar="K2NiF4_pipeline/01_structure/POSCAR", element="Ni",
            nonmag=("K", "F"), supercell=(4, 4, 1)),
        shells=(
            Shell("J1", 4.0172),
            Shell("J2", 4.0172 * np.sqrt(2.0)),
            Shell("J3", 2.0 * 4.0172),
        ),
        orders=(
            Order("fm", lambda c: np.ones(c.n)),
            Order("g_type", lambda c: sign_pow(c.n0 + c.n1)),
            Order("stripe_a", lambda c: sign_pow(c.n0)),
            Order("stripe_b", lambda c: sign_pow(c.n1)),
            Order("p4_a", lambda c: sign_pow(np.floor(c.n0 / 2.0))),
            Order("mix_ab", lambda c: sign_pow(c.n0) * sign_pow(np.floor(c.n1 / 2.0))),
        ),
        selected=("fm", "stripe_a", "p4_a", "mix_ab"),
        spin=1.0,
        policy="shift_search", nshift=1,
        fixture="K2NiF4_pipeline/02_supercell/cases/magnetic_orders.json",
        fixture_kind="coefficients_per_cell",
    ),

    # ---------------- 方形 (La2CuO4): 同一方形格子的**定量参考**材料 ----
    # K2NiF4 给的是无 Jahn-Teller 畸变/无八面体倾转的理想方形胞（结构原型），
    # 但它自己的 magnon 拟合 J 未能从可获取来源核实 ⇒ literature_J.json 里 J 留
    # null + _note（只有 T_N = 97.1 K 作实验锚点，符合 WORKFLOW §63 惯例）。
    # La2CuO4 是同一 K2NiF4 型方形格子的**定量锚点**（G5 闸门）:
    #   J(近邻, 反铁磁) = 111.8 ± 4 meV (295 K) / 104.1 ± 4 meV (10 K)
    #   J'(对角, 铁磁)  = -11.4 ± 3 meV (295 K) / -18   ± 3 meV (10 K)
    #   T_N = 325 ± 5 K
    #   (Coldea et al., "Spin Waves and Electronic Interactions in La2CuO4",
    #    PRL 86, 5377 (2001); arXiv:cond-mat/0006384 —— 原文约定 H = +J Σ_<ij> S·S,
    #    J>0 = 反铁磁。换算到本仓库约定 H = -J1 Σ S·S (J>0 = 铁磁):
    #    J1 = -111.8 meV (295 K 值) 或 -104.1 meV (10 K 值), J2 = +11.4 或 +18 meV。)
    # 结构与 K2NiF4 完全同型 (shells 形状与 6 个磁序的定义逐字相同, 只换材料/自旋),
    # 所以两个 spec 的 coefficients_per_cell 必然逐位相等 —— 这不是冗余, 而是
    # 「规格层只吃几何、材料信息不进系数」这条设计的回归证据。
    Spec(
        key="sq_la2cuo4",
        label="方形单离子 (La2CuO4)",
        structure=StructureSpec(
            poscar="La2CuO4_pipeline/01_structure/POSCAR", element="Cu",
            nonmag=("La", "O"), supercell=(4, 4, 1)),
        shells=(
            Shell("J1", 3.7992),
            Shell("J2", 3.7992 * np.sqrt(2.0)),
            Shell("J3", 2.0 * 3.7992),
        ),
        orders=(
            Order("fm", lambda c: np.ones(c.n)),
            Order("g_type", lambda c: sign_pow(c.n0 + c.n1)),
            Order("stripe_a", lambda c: sign_pow(c.n0)),
            Order("stripe_b", lambda c: sign_pow(c.n1)),
            Order("p4_a", lambda c: sign_pow(np.floor(c.n0 / 2.0))),
            Order("mix_ab", lambda c: sign_pow(c.n0) * sign_pow(np.floor(c.n1 / 2.0))),
        ),
        selected=("fm", "stripe_a", "p4_a", "mix_ab"),
        spin=0.5,
        policy="shift_search", nshift=1,
        fixture="La2CuO4_pipeline/02_supercell/cases/magnetic_orders.json",
        fixture_kind="coefficients_per_cell",
    ),

    # ---------------- C 心矩形 (CrOCl): 2 磁位/胞, 90° 晶格, a≠b, 非受挫 ----
    # 与矩形 (CrSBr) 的区别是**两点物理差异**, 不是风格差异:
    #  (1) 磁亚晶格是 2 个 Cr/胞的 C 心网 (Cr 在 (0,0,z) 与 (1/2,1/2,z')), 所以磁序
    #      多一个「同胞两 Cr 反平行」的 basis 自由度 (afm_c), 矩形规格表达不出来;
    #  (2) 最近两个 Cr-Cr 壳只差 0.163 A:
    #         J1 = 3.034707132522032 A  x4  <- (1/2,1/2) C 心键, **连同胞两个 Cr**
    #         J2 = 3.19795066          A  x2  <- 沿 a (同一亚晶格)
    #         J3 = 3.86026201          A  x2  <- 沿 b (同一亚晶格)
    #      build_rectangular_cases.py 的 magnetic_shells(tol=0.2) 会把 J1/J2 **静默合并**
    #      成一个 J1 (两个物理上不同的交换路径算成一个, 反解照样收敛、不报错)。
    #      所以这里必须显式声明壳距离, 并把 tol 压到 0.05。
    # 核心结构事实 (也是回归不变式): C 心矢量连的是同一胞内的两个 basis ⇒ basis 符号
    # 对一切 (n0,n1)-奇偶序不可见, afm_c_a == afm_a, afm_c_diag == afm_diag;
    # 反过来 fm 与 afm_c 只差 c1 (c2/c3 完全相同)。selected 取物理最直白的四个:
    # 4x4x1 下 cond(A)=64.02, det=524288 (非奇异, 远低于 1e4 闸门)。
    # 与方形不同, 这里**不需要**周期 4 的序: 周期 2 的四个序已经给出互不相同的
    # (c1,c2,c3), J3 与 e0 不简并。
    Spec(
        key="crect_crocl",
        label="C 心矩形多亚晶格 (CrOCl)",
        structure=StructureSpec(
            poscar="CrOCl_pipeline/01_structure/POSCAR", element="Cr",
            nonmag=("O", "Cl"), supercell=(4, 4, 1)),
        shells=(
            Shell("J1", 3.034707132522032, tol=0.05),
            Shell("J2", 3.19795066, tol=0.05),
            Shell("J3", 3.86026201, tol=0.05),
        ),
        orders=(
            Order("fm", lambda c: np.ones(c.n)),
            Order("afm_c", _basis_sign),
            Order("afm_a", lambda c: sign_pow(c.n0)),
            Order("afm_b", lambda c: sign_pow(c.n1)),
            Order("afm_diag", lambda c: sign_pow(c.n0 + c.n1)),
            Order("afm_c_a", lambda c: _basis_sign(c) * sign_pow(c.n0)),
            Order("afm_c_diag", lambda c: _basis_sign(c) * sign_pow(c.n0 + c.n1)),
        ),
        selected=("fm", "afm_c", "afm_a", "afm_b"),
        spin=1.5,
        policy="shift_search", nshift=1,
        fixture="CrOCl_pipeline/02_supercell/cases/magnetic_orders.json",
        fixture_kind="coefficients_per_cell",
    ),

    # ---------------- 蜂窝 (MnPS3): 2 磁位/胞, 60° 晶格, neel + 两个 M 点序 ----
    Spec(
        key="honeycomb_mnps3",
        label="蜂窝单离子 (MnPS3)",
        structure=StructureSpec(
            poscar="MnPS3_pipeline/01_structure/POSCAR", element="Mn",
            nonmag=("P", "S"), supercell=(2, 2, 1)),
        shells=(
            Shell("J1", 3.506189543925007),
            Shell("J2", 6.072013676948339),
            Shell("J3", 7.007530851884672),
        ),
        orders=(
            Order("fm", lambda c: np.ones(c.n)),
            Order("neel", lambda c: np.where(
                sublattice_colors(c.mag_idx, c.frac_sc, c.lat_sc, c.d_nn) == 0, 1.0, -1.0)),
            Order("M0", lambda c: sign_pow(c.n0)),
            Order("M1", lambda c: sign_pow(c.n1)),
            Order("M01", lambda c: sign_pow(c.n0 + c.n1)),
        ),
        post=_honeycomb_zigzag_stripy,
        spin=2.5,
        policy="shift_search",
        fixture="MnPS3_pipeline/02_supercell/cases/magnetic_orders.json",
        fixture_kind="coefficients_per_cell",
    ),

    # ---------------- 簇 (Nb3I8): 三聚体净矩粗粒化, 簇心晶格定壳 ----
    Spec(
        key="cluster_nb3i8",
        label="簇磁体 (Nb3I8 三聚体)",
        structure=StructureSpec(
            poscar="Nb3I8_pipeline/01_structure/POSCAR", element="Nb",
            nonmag=("I",), supercell=(5, 5, 1)),
        cluster=ClusterModel(
            cutoff=3.5,
            intra=Shell("intra", 2.9911229298393174),
            inter=(Shell("J1", 7.6988), Shell("J2", 13.3348), Shell("J3", 15.3977)),
            R=3, tol=0.5),
        orders=_cluster_orders(),
        spin=1.0,
        fixture="Nb3I8_pipeline/02_supercell/cases_j123/magnetic_orders.json",
        fixture_kind="cluster",
    ),

    # ---------------- FGT: 双亚晶格, 同距离的 AA/BB 壳必须按 pair 区分 ----
    Spec(
        key="fgt_fe3gete2",
        label="双亚晶格铁磁网络 (Fe3GeTe2)",
        structure=StructureSpec(
            poscar="Fe3GeTe2_pipeline/01_structure/POSCAR", element="Fe",
            nonmag=("Ge", "Te"), supercell=(3, 3, 1),
            sublattice_names=("B0", "B1", "A"),
            sublattice_family=("B", "B", "A")),
        shells=(
            Shell("J1", 2.63, tol=0.25, pair="AB"),
            Shell("J2", 2.55, tol=0.25, pair="BB"),
            Shell("J3", 3.99, tol=0.25, pair="AA"),
            Shell("J4", 3.99, tol=0.25, pair="BB"),
        ),
        orders=(
            Order("fm", lambda c: np.ones(c.n)),
            Order("afm_ab", lambda c: np.where(c.sub_family == "A", 1.0, -1.0)),
            Order("afm_b_vert", lambda c: np.where(c.sub_name == "B0", 1.0, -1.0)),
            Order("afm_a", lambda c: np.where(c.sub_family == "A",
                                              sign_pow(c.n0 + c.n1), 1.0)),
            Order("afm_b_plane", lambda c: np.where(c.sub_family == "A", 1.0,
                                                    sign_pow(c.n0 + c.n1))),
        ),
        selected=("fm", "afm_ab", "afm_b_vert", "afm_a", "afm_b_plane"),
        spin=1.0,
        fixture="2d_magnet_pipeline/tests/fixtures/fgt_cases/magnetic_orders.json",
        fixture_kind="configurations",
    ),

    # ---------------- 六角 kagome (Fe3Sn): 3 磁位/胞, 最近邻 z=4 ----------------
    # 3 个 Fe 亚晶格 = kagome 的 A/B/C。壳按理想 kagome 的距离比 1 : √3 : 2 (单位 a/2):
    #   J1 d=a/2    (z=4, 三角形内最近邻)
    #   J2 d=a√3/2  (z=4, 跨六角)
    #   J3 d=a      (z=6, 同亚晶格/对角)
    # 理想 kagome 里 J1 == J1' (对称性), 故只有 3 个独立壳, 对应 P4 计划的三自由度。
    # **简并警告**: q=0 一族 (fm / q0_120 / q0_ud / q0_duu) 的设计矩阵秩只有 2
    # (J1 与 J2 两列在该族内比值恒定), 分不开 J1/J2; 必须靠 √3×√3 序 (sq3_*) 才能分开。
    # P6 起默认不生成 120° (q0_120/sq3_120), 由 q0_ud 充当 q=0 代表;
    # selected 实测秩 = 4、cond≈1.0e2 (见 tests/test_magnetic_orders.py 的 COND_MAX 闸门)。
    Spec(
        key="kagome_fe3sn",
        label="六角 kagome (Fe3Sn)",
        structure=StructureSpec(
            poscar="Fe3Sn_pipeline/01_structure/POSCAR", element="Fe",
            nonmag=("Sn",), supercell=(3, 3, 1),
            sublattice_names=("A", "B", "C")),
        shells=(
            Shell("J1", 2.70, tol=0.30),
            Shell("J2", 4.676537, tol=0.30),
            Shell("J3", 5.40, tol=0.30),
        ),
        orders=(
            Order("fm", lambda c: np.ones(c.n)),
            Order("q0_120", lambda c: planar120_sub(c.sub_label, 3)),
            Order("q0_ud", lambda c: sublattice_signs(c.sub_label, (1.0, 1.0, -1.0))),
            Order("q0_duu", lambda c: sublattice_signs(c.sub_label, (1.0, -1.0, -1.0))),
            Order("sq3_120", lambda c: planar120_sqrt3(c.sub_label, c.n0, c.n1, 3)),
            Order("sq3_updntri", lambda c: sign_pow(sqrt3_label(c.n0, c.n1, 3))),
            Order("sq3_alt", lambda c: sign_pow(np.asarray(c.sub_label)
                                                + sqrt3_label(c.n0, c.n1, 3))),
        ),
        selected=("fm", "q0_ud", "sq3_updntri", "sq3_alt"),
        # S=1 占位: Fe3Sn 的 S 未在仓库定值, 反解出的 J 绝对值随 S² 缩放 (见 knowledge/LIMITATIONS.md)。
        spin=1.0,
        policy="shift_search",
        nshift=1,
        fixture="Fe3Sn_pipeline/02_supercell/cases/magnetic_orders.json",
        fixture_kind="coefficients_per_cell",
    ),
)}


# =====================================================================
# 7b. P0.3: 自描述块 (order_spec[] + model_spec)
#     写进 cases/magnetic_orders.json, 让「系数行 ↔ 磁序 ↔ 壳定义」可被机器复查
# =====================================================================

# key -> 晶格家族名 (写进产物; 新晶格往这里加: P1 已加 square, P2 已加 rect_centered, P4 已加 kagome)
LATTICE_OF = {
    "tri_mni2": "triangular",
    "rect_crsbr": "rectangular",
    "sq_k2nif4": "square",
    "sq_la2cuo4": "square",
    "crect_crocl": "rect_centered",
    "honeycomb_mnps3": "honeycomb",
    "cluster_nb3i8": "cluster",
    "fgt_fe3gete2": "fgt_bilayer",
    "kagome_fe3sn": "kagome",
}

# 磁序的**文字规则** (纯文档, 便于人复查产物; 不在计算路径上)
ORDER_RULES = {
    "fm": "sigma_i = +1 (全同向)",
    "120deg": "color=(n0+n1)%3; theta=color*120deg; S=(cos,sin,0)",
    "stripy_a": "sigma = (-1)^n0",
    "stripy_b": "sigma = (-1)^n1",
    "stripy_ab": "sigma = (-1)^(n0+n1)",
    "uud": "color=(n0+n1)%3; color==0 -> -1, 其余 +1",
    "afm_a": "sigma = (-1)^n0 (沿 a 反平行条纹)",
    "afm_b": "sigma = (-1)^n1 (沿 b 反平行条纹)",
    "afm_diag": "sigma = (-1)^(n0+n1) (棋盘)",
    "afm_c": "sigma = basis (C 心网 Néel: 2 磁原子/胞, 同胞两个磁原子反平行)",
    "afm_c_a": "sigma = basis * (-1)^n0 (与 afm_a 逐位相同: C 心矢量连同一胞内两个 basis)",
    "afm_c_diag": "sigma = basis * (-1)^(n0+n1) (与 afm_diag 逐位相同; 不可与 afm_diag 同选)",
    "g_type": "sigma = (-1)^(n0+n1) (方形 G 型/Néel, 棋盘)",
    "stripe_a": "sigma = (-1)^n0 (方形 C 型条纹, 沿 a 交替)",
    "stripe_b": "sigma = (-1)^n1 (方形 C 型条纹, 沿 b 交替)",
    "p4_a": "sigma = (-1)^floor(n0/2) (方形周期 4; 打破 J3/E0 简并的必要条件)",
    "mix_ab": "sigma = (-1)^n0 * (-1)^floor(n1/2) (周期 2 x 周期 4 混合)",
    "neel": "最近邻 BFS 二着色 (sublattice_colors): A->+1, B->-1",
    "M0": "sigma = (-1)^n0",
    "M1": "sigma = (-1)^n1",
    "M01": "sigma = (-1)^(n0+n1)",
    "zigzag": "三个 M 点序中 J1 系数最大者 (等价旧代码 min(有向键和))",
    "stripy": "与 zigzag 的 J1 同值、J3 不同的 M 点序",
    "fm_trim": "簇内 2 up 1 down; 簇净矩 m_a 全同向",
    "afm_ab": "A(2c) 与 B(4e) 反平行",
    "afm_b_vert": "B 胞内竖直对反平行, A 不动",
    "afm_b_plane": "A 不动, B 面内棋盘 (-1)^(n0+n1)",
    # --- P4.1 六角 kagome (Fe3Sn): 相位按**亚晶格**取, 不是按晶胞 ---
    "q0_120": "kagome q=0: theta = sub_label*120deg; S=(cos,sin,0) (同胞 A/B/C 互差 120°)",
    "q0_ud": "kagome q=0: sigma = (+1,+1,-1) 按亚晶格 (2up1down, 每个三角形同型)",
    "q0_duu": "kagome q=0: sigma = (+1,-1,-1) 按亚晶格 (1up2down)",
    "sq3_120": "kagome sqrt3xsqrt3: theta = (sub_label + (n0+n1)%3)*120deg; S=(cos,sin,0)",
    "sq3_updntri": "kagome sqrt3xsqrt3: sigma = (-1)^((n0+n1)%3) (√3 胞一相翻转)",
    "sq3_alt": "kagome sqrt3xsqrt3: sigma = (-1)^(sub_label + (n0+n1)%3) (亚晶格 + √3 胞交替)",
}


def _is_noncollinear(spins) -> bool:
    s = np.asarray(spins, dtype=float)
    return s.ndim == 2 and s.shape[1] == 3 and bool(np.any(np.abs(s[:, :2]) > 1e-9))


def _rule_of(name: str) -> str:
    if name in ORDER_RULES:
        return ORDER_RULES[name]
    if name.startswith("q_"):
        return "簇净矩 m_a 沿 q 相位旋转 (平面螺旋)"
    if name == "s12":
        return "12 个 Nb 反号 (簇内 1 up 2 down)"
    if name in ("interA", "interB"):
        return "簇间构型 (单胞粗粒化, 见 build_cluster_cases.py)"
    return "(未登记规则)"


def order_spec_block(spec: Spec, results, ctx) -> dict:
    """产出 order_spec[] (逐磁序) + model_spec (晶格/壳定义)。纯读, 不改任何计算。"""
    if spec.cluster is not None:
        cm = spec.cluster
        shells = [dict(name=cm.intra.name, d=float(cm.intra.d), tol=float(cm.intra.tol),
                       pair=None)]
        for s in cm.inter:
            shells.append(dict(name=s.name, d=float(s.d), tol=float(s.tol), pair=s.pair))
    else:
        shells = [dict(name=s.name, d=float(s.d), tol=float(s.tol), pair=s.pair)
                  for s in spec.shells]

    order_spec = []
    for nm, r in results.items():
        order_spec.append(dict(
            name=nm,
            kind="noncollinear" if _is_noncollinear(r.spins) else "collinear",
            rule=_rule_of(nm),
            coeffs=[float(x) for x in r.coeffs],
        ))

    model_spec = dict(
        schema="model_spec/1",
        lattice=LATTICE_OF.get(spec.key, spec.key),
        label=spec.label,
        element=spec.structure.element,
        nonmag=list(spec.structure.nonmag),
        spin=float(spec.spin),
        S2=float(spec.spin) ** 2,
        supercell=list(spec.structure.supercell),
        policy=spec.policy,
        nshift=int(spec.nshift),
        noncollinear=any(o["kind"] == "noncollinear" for o in order_spec),
        n_mag=int(ctx.n),
        shells=shells,
    )
    if spec.cluster is not None:
        model_spec["cluster"] = dict(cutoff=float(spec.cluster.cutoff),
                                     R=int(spec.cluster.R),
                                     tol=float(spec.cluster.tol),
                                     n_clusters=int(ctx.n_groups),
                                     cluster_size=len(ctx.groups[0]) if ctx.groups else None)
    return dict(order_spec=order_spec, model_spec=model_spec)


def emit_spec_into(path: str, key: str, tol: float = 1e-9):
    """把 order_spec[]/model_spec 注入既有 magnetic_orders.json, 并**顺手复核**系数。

    返回 (block, messages)。系数与规格层算不出一致时, message 里带 "!!" —— 这条路径
    同时是 P0.5 「脚本对契约」的一环: 产物一旦与规格漂移, 立刻可见。
    """
    spec = REGISTRY.get(key)
    if spec is None:
        raise KeyError(f"未知规格 {key!r}; 可选 {', '.join(REGISTRY)}")
    results, ctx = compute(spec)
    block = order_spec_block(spec, results, ctx)

    with open(path, encoding="utf-8") as fh:
        data = json.load(fh)

    msgs = []
    have = data.get("coefficients_per_cell") or {}
    if not have:
        have = _expected_from_fixture(data, spec.fixture_kind)
    for o in block["order_spec"]:
        want = have.get(o["name"])
        if want is None:
            msgs.append(f"!! {o['name']}: 产物里没有这一行, 已补 (规格层算得 {o['coeffs']})")
            continue
        err = float(np.max(np.abs(np.asarray(o["coeffs"], float) - np.asarray(want, float))))
        msgs.append(f"{'ok' if err <= tol else '!!'} {o['name']}: 与产物偏差 {err:.3g}")

    data["order_spec"] = block["order_spec"]
    data["model_spec"] = block["model_spec"]
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    return block, msgs


# =====================================================================
# 8. 回归: 与 5 条旧链路的 JSON 逐位比对
# =====================================================================

def _expected_from_fixture(data: dict, kind: str) -> dict:
    if kind == "coefficients_per_cell":
        return {k: list(v) for k, v in data["coefficients_per_cell"].items()}
    if kind == "configurations":
        return {c["name"]: list(c["coeffs"]) for c in data["configurations"]}
    if kind == "cluster":
        out = {}
        for k, v in data["coefficients"].items():
            out[k] = [v["C_intra"]] + list(v["C_inter"]) if isinstance(v, dict) else list(v)
        return out
    raise ValueError(f"未知 fixture_kind={kind!r}")


def regress(spec: Spec, tol: float = 1e-9):
    results, ctx = compute(spec)
    if spec.fixture is None:
        return None, ctx, "无 fixture"
    fx = resolve(spec.fixture)
    if not os.path.exists(fx):
        # 产物未生成 (新晶格刚进规格层, build_*.py 还没跑过) -> 明确报 SKIP, 不静默
        return None, ctx, f"SKIP (产物未生成: {spec.fixture})"
    with open(fx, encoding="utf-8") as fh:
        data = json.load(fh)
    exp = _expected_from_fixture(data, spec.fixture_kind)
    rows, nbad = [], 0
    for name, want in exp.items():
        if name not in results:
            rows.append((name, None, want, "缺磁序定义"))
            nbad += 1
            continue
        got = np.asarray(results[name].coeffs, dtype=float)
        want = np.asarray(want, dtype=float)
        if got.shape != want.shape:
            rows.append((name, got, want, f"长度不符 {got.shape} vs {want.shape}"))
            nbad += 1
            continue
        err = float(np.max(np.abs(got - want))) if got.size else 0.0
        ok = err <= tol
        nbad += 0 if ok else 1
        rows.append((name, got, want, "ok" if ok else f"偏差 {err:.3g}"))
    return rows, ctx, ("PASS" if nbad == 0 else f"FAIL ({nbad}/{len(exp)})")


def _fmt(v) -> str:
    if v is None:
        return "-"
    return "[" + ", ".join(f"{float(x):+.6g}" for x in v) + "]"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="磁序规格层 (P0.1)")
    ap.add_argument("--list", action="store_true", help="列出注册表")
    ap.add_argument("--spec", help="只算某条规格, 打印系数")
    ap.add_argument("--regress", action="store_true", help="与旧链路 JSON 逐位回归")
    ap.add_argument("--emit-spec", metavar="KEY",
                    help="把 order_spec[]/model_spec 注入该链路的 magnetic_orders.json (P0.3)")
    ap.add_argument("--into", default=None,
                    help="配合 --emit-spec: 目标 magnetic_orders.json (默认用注册表里的 fixture)")
    ap.add_argument("--dry-run", action="store_true", help="配合 --emit-spec: 只复核不写回")
    ap.add_argument("--tol", type=float, default=1e-9)
    args = ap.parse_args(argv)

    if args.list or not (args.spec or args.regress or args.emit_spec):
        print(f"{'key':22s} {'label':30s} {'shells':>7s} {'orders':>7s}  fixture")
        for s in REGISTRY.values():
            nsh = len(s.cluster.inter) + 1 if s.cluster else len(s.shells)
            print(f"{s.key:22s} {s.label:30s} {nsh:7d} {len(s.orders):7d}  "
                  f"{os.path.basename(s.fixture) if s.fixture else '-'}")
        if not (args.spec or args.regress or args.emit_spec):
            return 0

    if args.emit_spec:
        spec = REGISTRY.get(args.emit_spec)
        if spec is None:
            print(f"未知规格 {args.emit_spec!r}; 可选: {', '.join(REGISTRY)}", file=sys.stderr)
            return 2
        target = args.into or (resolve(spec.fixture) if spec.fixture else None)
        if not target:
            print(f"{spec.key} 没有 fixture, 必须用 --into 指定目标文件", file=sys.stderr)
            return 2
        if args.dry_run:
            results, ctx = compute(spec)
            block = order_spec_block(spec, results, ctx)
            msgs = [f"(dry-run) 将写入 {len(block['order_spec'])} 条 order_spec, "
                    f"lattice={block['model_spec']['lattice']}"]
        else:
            block, msgs = emit_spec_into(target, spec.key, args.tol)
        print(f"[{spec.key}] {spec.label}  ->  {target}")
        for m in msgs:
            print(f"  {m}")
        print(f"  model_spec: lattice={block['model_spec']['lattice']} "
              f"policy={block['model_spec']['policy']} "
              f"shells={[s['name'] for s in block['model_spec']['shells']]}")
        if not args.dry_run and not args.spec and not args.regress:
            return 0

    if args.spec:
        spec = REGISTRY.get(args.spec)
        if spec is None:
            print(f"未知规格 {args.spec!r}; 可选: {', '.join(REGISTRY)}", file=sys.stderr)
            return 2
        results, ctx = compute(spec)
        print(f"[{spec.key}] {spec.label}")
        print(f"  原胞磁性位点 {len(set(ctx.prim_idx[ctx.mag_idx].tolist()))} 个 -> "
              f"亚晶格 {sorted(set(ctx.sub_name.tolist()))}")
        print(f"  超胞 {ctx.supercell} -> {len(ctx.atoms_sc)} 原子, 磁位点 {ctx.n}"
              + (f", 簇 {ctx.n_groups}" if ctx.groups else ""))
        for nm, r in results.items():
            print(f"  {nm:12s} {_fmt(r.coeffs)}")

    if args.regress:
        print("=" * 96)
        print(f"回归: 磁序规格层 vs 各链路的 magnetic_orders.json ({len(REGISTRY)} 条规格)")
        print("=" * 96)
        bad = 0
        skipped = []
        for spec in REGISTRY.values():
            try:
                rows, ctx, verdict = regress(spec, args.tol)
            except FileNotFoundError as exc:
                # 规格指向工作区里的项目目录时 (如 MnI2 链路), 独立分发的 skill 包里没有这些目录。
                # 结构缺失 = 本条无法回归, 明确 SKIP, 不把整条命令判失败。
                skipped.append(f"{spec.key}: SKIP (结构不在本包内: {exc.filename})")
                continue
            if rows is None:
                skipped.append(f"{spec.key}: {verdict}")
                continue
            print(f"\n[{spec.key}] {spec.label}  ->  {verdict}")
            print(f"  几何: 超胞 {ctx.supercell}, 磁位点 {ctx.n}"
                  + (f", 簇 {ctx.n_groups}×{len(ctx.groups[0])}" if ctx.groups else ""))
            for name, got, want, note in rows:
                mark = "  " if note == "ok" else "!!"
                print(f" {mark} {name:12s} got={_fmt(got):34s} want={_fmt(want):34s} {note}")
            if verdict != "PASS":
                bad += 1
        for s in skipped:
            print(f"\n[skip] {s}")
        print("\n" + "=" * 96)
        print("全部 PASS" if bad == 0 else f"{bad} 条规格未通过")
        print("=" * 96)
        return 0 if bad == 0 else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
