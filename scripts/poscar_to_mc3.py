#!/usr/bin/env python
"""
poscar_to_mc.py
把 POSCAR + J1/J2/MAE 转成 MnSTe_repro 引擎 (spin_mc.py) 能直接读的 Monte Carlo 模型输入文件 (JSON)。

物理约定 (务必对齐):
  笔记/DFT 侧的哈密顿:  H = -J1 * sum_<ij> S_i.S_j  -  J2 * sum_<<ij>> S_i.S_j     (J>0 = 铁磁)
  spin_mc 引擎侧能量:    E = sum_<ij> S_i . J_mat . S_j  +  sum_i S_i . A . S_i  -  g muB B . sum S_i
  引擎用单位矢量 |S|=1 做经典 MC, 所以交换要乘 S^2 才能对上量子能量尺度:
      J_mat(NN)  = -(J1 * S^2) * I3
      J_mat(NNN) = -(J2 * S^2) * I3
  单离子各向异性 (易轴 z, 每离子 MAE 单位 meV):
      A = diag(0, 0, -MAE)

用法:
  python poscar_to_mc.py POSCAR --element Cr --S 1.5 --J1 3.05 --J2 0.73 --MAE 0.6 -o cri3_mc_model.json
  python poscar_to_mc.py POSCAR --element Cr --S 1.5 --J1 3.05 -o cri3_mc_model.json   # 只取 J1

  # P0.4: 显式壳定义 (绕开「距离排序硬映射 J1/J2/J3」), 或直接吃 P0.2 的 J_result.json
  python poscar_to_mc3.py POSCAR --element K Cu --S 1.0 0.5 \
      --shell-def shells.json -o k2cuf4_mc_model.json
  python poscar_to_mc3.py POSCAR --element Cr --S 1.5 \
      --jfile 02_supercell/cases/J_result.json -o crsbr_mc_model.json

参数:
  --element      磁性元素符号, **可多个** (顺序决定 basis 分组顺序)
  --S            自旋大小; 给一个值通吃, 或逐元素给 (与 --element 等长)
  --J1/--J2/--J3 各向同性交换 (meV, H=-J sum S.S, J>0 铁磁); 用了 --jfile 可省
  --jfile        P0.2 统一产物 J_result.json: 同时提供 **壳定义 (d/tol/pair) + J**
  --shell-def    显式壳定义 JSON (文件或内联); 每项 {name, d, tol?, pair?, J_meV?, J_mat?}
                 - pair: 亚晶格对 ("AB"/"AA"/...) 或显式 basis 索引对 [[b1,b2], ...]
                 - J_mat: 3x3 矩阵, 给各向异性交换 / DMI 留门 (省略则用 J_meV 各向同性)
  (!) --jfile/--shell-def 一旦给出, 壳就**不再按距离自动排序** —— MC 模型与 J 拟合必须同壳,
      否则同一个 d 落到不同壳里, 能量对不上。这是 P0.4 存在的唯一理由。
  --MAE          单离子各向异性能 (meV/离子, 易轴 z), 可省略 (默认 0)
  --no-s2        不做 S^2 归一 (若你给的 J1/J2 已是经典单位自旋耦合)
  --nrange       找近邻时平移搜索范围 (默认 2)
  -o/--out       输出 JSON 路径
"""

import argparse
import json
import os
import sys

import numpy as np

I3 = np.eye(3)


# ----------------------------------------------------------------------------- POSCAR
def _is_float(tok):
    try:
        float(tok)
        return True
    except ValueError:
        return False


def parse_poscar(path):
    # encoding 必须显式 utf-8: 中文 Windows 的默认编码是 GBK, 而 POSCAR 注释行
    # 常带非 ASCII (材料名/说明), 不指定就会 UnicodeDecodeError 崩在注释行上。
    # 与本仓库其它 parse_poscar (build_square_cases.py:92 等) 保持一致。
    with open(path, encoding="utf-8") as f:
        lines = [ln.rstrip("\n") for ln in f]

    comment = lines[0]
    scale = float(lines[1].split()[0])
    lattice = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)])
    if scale < 0:  # negative scale = target volume
        vol = abs(np.linalg.det(lattice))
        lattice *= (abs(scale) / vol) ** (1.0 / 3.0)
    else:
        lattice *= scale

    idx = 5
    tok = lines[idx].split()
    if all(_is_float(t) for t in tok):
        symbols = None
        counts = [int(t) for t in tok]
        idx += 1
    else:
        symbols = tok
        counts = [int(t) for t in lines[idx + 1].split()]
        idx += 2

    if lines[idx].strip()[:1] in ("s", "S"):
        idx += 1
    mode = lines[idx].strip()[:1]
    idx += 1

    nat = sum(counts)
    coords = np.array([[float(x) for x in lines[idx + i].split()[:3]] for i in range(nat)])
    if mode in ("c", "C", "k", "K"):
        coords = coords @ np.linalg.inv(lattice)
    coords -= np.floor(coords)

    if symbols is None:
        symbols = ["X"] * nat
    else:
        symbols = [s for s, c in zip(symbols, counts) for _ in range(c)]

    return dict(comment=comment, lattice=lattice, symbols=symbols, coords=coords)


# ----------------------------------------------------------------------------- 2D 磁格
def get_magnetic_sites(pos, element):
    frac = np.array([pos["coords"][i, :2] for i, s in enumerate(pos["symbols"]) if s == element])
    if len(frac) == 0:
        raise SystemExit(f"POSCAR 里找不到元素 {element}")
    return frac


def is_symmetry(frac, t, tol=1e-3):
    shifted = frac + np.asarray(t)
    diff = frac[None, :, :] - shifted[:, None, :]
    diff -= np.round(diff)
    return bool(np.all(np.any(np.all(np.abs(diff) < tol, axis=2), axis=1)))


def reduce_to_primitive(frac, cell, tol=1e-3, max_den=4):
    """把磁性亚晶格约化到最小重复单元。cell 为 2x2 行矢矩阵。
    平移按分数坐标搜索 (分母 up to max_den), 以覆盖超胞输入。"""
    cands = {}
    for d in range(1, max_den + 1):
        for i in range(-2 * d, 2 * d + 1):
            for j in range(-2 * d, 2 * d + 1):
                if i == 0 and j == 0:
                    continue
                t = np.array([i, j], dtype=float) / d
                key = (round(t[0], 6), round(t[1], 6))
                if key in cands:
                    continue
                if is_symmetry(frac, t, tol):
                    cands[key] = (np.linalg.norm(t @ cell), t)
    ordered = sorted(cands.values(), key=lambda x: x[0])

    if not ordered:
        return cell, frac, False

    t1 = ordered[0][1]
    t2 = None
    for _, t in ordered[1:]:
        # 2D 叉积的 z 分量 = t1[0]*t[1] - t1[1]*t[0], 判断两平移矢量是否非共线
        if abs(t1[0] * t[1] - t1[1] * t[0]) > tol:
            t2 = t
            break
    if t2 is None:
        return cell, frac, False

    new_cell = np.array([t1, t2]) @ cell
    if abs(np.linalg.det(new_cell)) < 1e-6:
        return cell, frac, False
    if abs(np.linalg.det(new_cell)) >= abs(np.linalg.det(cell)) - 1e-6:
        return cell, frac, False
    if np.linalg.det(new_cell) < 0:
        new_cell = new_cell[::-1]

    red = frac @ cell @ np.linalg.inv(new_cell)
    red -= np.floor(red + 1e-9)
    uniq = []
    for r in red:
        if not any(np.all(np.abs(((r - u + 0.5) % 1.0) - 0.5) < tol) for u in uniq):
            uniq.append(r)
    return new_cell, np.array(uniq), True


def find_bonds(basis, prim, nrange=2, tol=1e-3):
    """返回 {canonical_key: distance}, key=(k,l,i,j) 表示 basis[l] 平移 (i,j) 与 basis[k] 成键。"""
    n = len(basis)
    bonds = {}
    for k in range(n):
        rk = basis[k] @ prim
        for l in range(n):
            for i in range(-nrange, nrange + 1):
                for j in range(-nrange, nrange + 1):
                    if k == l and i == 0 and j == 0:
                        continue
                    rl = (basis[l] + np.array([i, j])) @ prim
                    d = float(np.linalg.norm(rk - rl))
                    if d < tol:
                        continue
                    key = (k, l, i, j)
                    rev = (l, k, -i, -j)
                    ck = min(key, rev)
                    if ck not in bonds:
                        bonds[ck] = d
    return bonds


def shell_of(distances, tol=1e-2):
    ds = sorted(set(round(d, 5) for d in distances))
    shells = []
    for d in ds:
        if not shells or d - shells[-1][-1] > tol:
            shells.append([d])
        else:
            shells[-1].append(d)
    return [float(np.mean(s)) for s in shells]


# ============================================================================= P0.4
# 显式壳定义 / 多磁性元素 / 逐壳 J —— 「壳」不再由距离排序猜出来
# =============================================================================
def magnetic_sites_multi(pos, elements):
    """按 elements 的顺序取磁性位点。多元素时逐元素分组 —— 分组顺序决定 basis 排序。"""
    if isinstance(elements, str):
        elements = [elements]
    idxs = []
    for el in elements:
        sub = [i for i, s in enumerate(pos["symbols"]) if s == el]
        if not sub:
            raise SystemExit(f"POSCAR 里找不到元素 {el}")
        idxs.extend(sub)
    return np.array([pos["coords"][i, :2] for i in idxs])


def load_shell_defs(spec):
    """读壳定义: 支持 P0.2 的 J_result.json / magnetic_orders.json / 内联 JSON / 纯列表。

    返回 (defs, origin)。defs 每项: {name, d, tol, pair, J_meV, J_mat}
    """
    text = spec.strip()
    inline = text.startswith("{") or text.startswith("[")
    if inline:
        data = json.loads(text)
        origin = "<inline>"
    else:
        if not os.path.exists(spec):
            raise SystemExit(f"壳定义文件不存在: {spec}")
        # BOM 分派: utf-8-sig 吃 Windows 编辑器/Set-Content 的 UTF-8 BOM;
        # PowerShell 的 `python ... > x.json` 会写成 **UTF-16LE+ BOM** (0xFF 0xFE),
        # 用 utf-8-sig 直接 UnicodeDecodeError -> 按前两字节显式选编码, 不猜。
        with open(spec, "rb") as f:
            head = f.read(2)
        enc = "utf-16" if head in (b"\xff\xfe", b"\xfe\xff") else "utf-8-sig"
        with open(spec, encoding=enc) as f:
            data = json.load(f)
        origin = os.path.abspath(spec) + ("" if enc == "utf-8-sig" else f" [{enc}]")

    defs = []
    if isinstance(data, dict) and isinstance(data.get("shells"), list) \
            and data["shells"] and isinstance(data["shells"][0], dict):
        jmap = data.get("J_meV") or {}
        for s in data["shells"]:
            jv = s.get("J_meV")
            defs.append(dict(name=str(s["name"]), d=float(s["d"]),
                             tol=float(s.get("tol", 0.3)), pair=s.get("pair"),
                             J_meV=(jv if jv is not None else jmap.get(str(s["name"]))),
                             J_mat=s.get("J_mat")))
    elif isinstance(data, dict) and isinstance(data.get("shells_angstrom"), (dict, list)):
        # 注意口径: shells_angstrom 是**规格层的 3D 键长** (CrOCl 的 C 心键 3.0347 A),
        # 而本脚本把 POSCAR 投到 2D 再搜索, 同一个键在这里是 2.506 A。所以拿
        # magnetic_orders.json 的 shells_angstrom 当 --jfile 会一条键都匹配不上
        # (下面的空壳检查会拦)。驱动层因此改用 --J1/--J2/--J3 按序传 J, 再用
        # magnetic_orders.json 的 per_site_bonds 与本模型的 per_site 交叉校验。
        sa = data["shells_angstrom"]
        items = list(sa.items()) if isinstance(sa, dict) else \
            [(f"J{i + 1}", v) for i, v in enumerate(sa)]
        for k, v in items:
            defs.append(dict(name=str(k), d=float(v), tol=0.3, pair=None,
                             J_meV=None, J_mat=None))
    elif isinstance(data, list):
        for i, s in enumerate(data):
            if isinstance(s, dict):
                defs.append(dict(name=str(s.get("name", f"J{i + 1}")), d=float(s["d"]),
                                 tol=float(s.get("tol", 0.3)), pair=s.get("pair"),
                                 J_meV=s.get("J_meV"), J_mat=s.get("J_mat")))
            else:
                defs.append(dict(name=f"J{i + 1}", d=float(s), tol=0.3, pair=None,
                                 J_meV=None, J_mat=None))
    else:
        raise SystemExit(f"认不出壳定义格式 (前 80 字): {spec[:80]}")
    if not defs:
        raise SystemExit("壳定义为空")
    return defs, origin


def parse_families(items, n_basis):
    """--family A=0 B=1,2 -> {"A": [0], "B": [1, 2]}。供 shell 的 pair="AB" 用。"""
    if not items:
        return None
    fam = {}
    for it in items:
        if "=" not in it:
            raise SystemExit(f"--family 需要 NAME=i,j 形式, 收到 {it!r}")
        nm, idxs = it.split("=", 1)
        vals = [int(x) for x in idxs.replace(",", " ").split() if x.strip()]
        for i in vals:
            if not 0 <= i < n_basis:
                raise SystemExit(f"--family {it!r}: basis 索引 {i} 超出 0..{n_basis - 1}")
        fam[nm.strip()] = vals
    return fam


def pair_basis_pairs(pair, fam, label=""):
    """把 shell 的 pair 归一成 {(b1,b2)} 允许集; None = 不过滤。"""
    if pair is None:
        return None
    if isinstance(pair, (list, tuple)) and pair and isinstance(pair[0], (list, tuple)):
        return {(min(int(a), int(b)), max(int(a), int(b))) for a, b in pair}
    if isinstance(pair, str) and len(pair) == 2 and fam and pair[0] in fam and pair[1] in fam:
        out = set()
        for a in fam[pair[0]]:
            for b in fam[pair[1]]:
                out.add((min(a, b), max(a, b)))
        return out
    raise SystemExit(
        f"壳 {label} 的 pair={pair!r} 无法解析。两种写法:\n"
        f"  (a) 显式 basis 索引对: \"pair\": [[0,1],[0,2]]\n"
        f"  (b) 亚晶格族字符串 \"AB\"/\"AA\"/\"BB\" + CLI --family A=0 B=1,2")


def assign_to_shells(bonds, defs, fam):
    """把 bonds 分到**显式**壳。返回 (per_shell_keys, unassigned)。

    规则: 每条键在「距离在该壳 tol 内 **且** 通过该壳 pair 过滤」的候选壳里挑**最近**的一个。
    - 挑最近 -> 壳之间不会互相抢键, 结果与 defs 的书写顺序无关
    - pair 过滤逐壳各判 -> FGT 那种 J1=AB@2.63 / J2=BB@2.55 的密集壳也能分开
    """
    allowed = [pair_basis_pairs(d.get("pair"), fam, d.get("name", i))
               for i, d in enumerate(defs)]
    out = [[] for _ in defs]
    unassigned = {}
    for key, d in sorted(bonds.items(), key=lambda kv: kv[1]):
        best, best_dev = None, None
        for si, sd in enumerate(defs):
            dev = abs(d - sd["d"])
            if dev > sd["tol"]:
                continue
            if allowed[si] is not None:
                k, l = key[0], key[1]
                if (min(k, l), max(k, l)) not in allowed[si]:
                    continue
            if best is None or dev < best_dev:
                best, best_dev = si, dev
        if best is None:
            unassigned[key] = d
        else:
            out[best].append(key)
    return out, unassigned


def jmat_of(sd, idx, args, s2):
    """逐壳 J 矩阵: J_mat (显式 3x3, 给各向异性/DMI 留门) 优先, 否则 J_meV 各向同性。"""
    if sd.get("J_mat") is not None:
        m = np.asarray(sd["J_mat"], dtype=float)
        if m.shape != (3, 3):
            raise SystemExit(f"壳 {sd['name']} 的 J_mat 形状 {m.shape} != (3,3)")
        return m
    j = sd.get("J_meV")
    if j is None and idx < 3:
        j = getattr(args, f"J{idx + 1}", None)
    if j is None:
        raise SystemExit(f"壳 {sd['name']} (d={sd['d']:.4f}) 没有 J: "
                         f"既无 J_meV/J_mat, 也没有 --J{idx + 1}")
    return (-float(j) * s2) * I3


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="POSCAR + J1/J2/MAE -> Monte Carlo 模型输入文件")
    ap.add_argument("poscar")
    ap.add_argument("--element", required=True, nargs="+",
                    help="磁性元素符号, 如 Cr; 可多个 (--element K Cu), 顺序 = basis 排序")
    ap.add_argument("--S", type=float, required=True, nargs="+", help="自旋大小, 如 Cr3+ = 1.5")
    ap.add_argument("--J1", type=float, default=None, help="最近邻 J1 (meV); 无 --jfile/--shell-def 时必给")
    ap.add_argument("--J2", type=float, default=0.0, help="次近邻 J2 (meV)")
    ap.add_argument("--J3", type=float, default=0.0, help="第三近邻 J3 (meV)")
    ap.add_argument("--shell-def", default=None,
                    help='显式壳定义 (JSON 文件路径或内联), 如 \'{"shells":[{"name":"J1",'
                         '"d":5.3,"tol":0.3,"pair":[[0,1]]}]}\'; 给了就不再按距离排序猜 J1/J2/J3')
    ap.add_argument("--jfile", default=None, help="j_result.py 写出的 J_result.json (P0.2 产物)")
    ap.add_argument("--family", nargs="+", default=None,
                    help="亚晶格族映射, 如 --family A=0 B=1,2; 供 shell 的 pair='AB' 用")
    ap.add_argument("--MAE", type=float, default=0.0, help="单离子各向异性 (meV/离子, 易轴 z)")
    ap.add_argument("--no-s2", action="store_true", help="J1/J2 已是经典耦合, 不再乘 S^2")
    ap.add_argument("--nrange", type=int, default=2)
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    elements = list(args.element)
    spins = list(args.S)
    if len(spins) == 1 and len(elements) > 1:
        spins = spins * len(elements)
    if len(spins) != len(elements):
        raise SystemExit(f"--S 给了 {len(spins)} 个值, --element 给了 {len(elements)} 个: "
                         f"要么等长, 要么 --S 只给一个 (自动套用)")
    if len(set(elements)) != len(elements):
        raise SystemExit(f"--element 有重复元素: {elements}")
    if len(set(spins)) > 1:
        raise SystemExit(
            "多个不同的 --S 需要逐壳的 S^2 归一 (壳两端自旋之积), 那是 P3 的工作。\n"
            f"当前引擎只有一个全局 S^2, 收到 {dict(zip(elements, spins))} 会静默算错, 拒绝继续。\n"
            "请用统一的 --S, 或等 P3 的逐壳 S_eff。")
    S = spins[0]

    pos = parse_poscar(args.poscar)
    cell3 = pos["lattice"]
    if max(abs(cell3[2, 0]), abs(cell3[2, 1]), abs(cell3[0, 2]), abs(cell3[1, 2])) > 0.1:
        print("[warn] 晶格 c 似乎不在 z 方向, 2D 投影可能不准", file=sys.stderr)
    cell2 = np.array([cell3[0, :2], cell3[1, :2]])

    frac = magnetic_sites_multi(pos, elements)
    n_raw = len(frac)

    prim, basis, reduced = reduce_to_primitive(frac, cell2)
    print(f"磁性元素 {'+'.join(elements)}: {n_raw} 个/输入胞 -> 约化后 {len(basis)} 个/原胞 "
          f"({'已约化' if reduced else '未约化'})")

    bonds = find_bonds(basis, prim, nrange=args.nrange)
    if not bonds:
        raise SystemExit("找不到任何近邻键")

    s2 = 1.0 if args.no_s2 else S ** 2

    src = args.jfile or args.shell_def
    if src:
        defs, origin = load_shell_defs(src)
        fam = parse_families(args.family, len(basis))
        per_shell, unassigned = assign_to_shells(bonds, defs, fam)
        names = [d["name"] for d in defs]
        shells_d = [float(d["d"]) for d in defs]
        # 吸收了一个离壳心较远的键 = 可能把不存在于该壳的键算进了 J。允许, 但必须说出来。
        for si, keys in enumerate(per_shell):
            devs = [abs(bonds[k] - shells_d[si]) for k in keys]
            if devs and max(devs) > 0.05:
                far = sorted({round(bonds[k], 4) for k in keys
                              if abs(bonds[k] - shells_d[si]) > 0.05})
                print(f"[warn] 壳 {names[si]} (d={shells_d[si]:.4f}) 吸收了偏离 {max(devs):.3f} A "
                      f"的键: {far[:6]} —— 该壳可能掺进了别的近邻壳, 请核对壳定义",
                      file=sys.stderr)
        reach = max(shells_d) + max(d["tol"] for d in defs)
        stray = {k: v for k, v in unassigned.items() if v <= reach}
        if stray:
            uniq = sorted(set(round(v, 4) for v in stray.values()))
            print(f"[warn] {len(stray)} 条键没落进任何显式壳 (虽在 {reach:.3f} A 内): {uniq[:8]}\n"
                  f"       壳定义漏了就会丢键 -> 能量对不上。补壳或放宽 tol。", file=sys.stderr)
        empty = [names[si] for si, keys in enumerate(per_shell) if not keys]
        if empty:
            raise SystemExit(
                f"壳 {empty} 一条键都没匹配到 -> 相应的 J 会落空 (MC 会少一整条耦合)。\n"
                f"  最可能的原因: 壳定义的距离口径与本脚本的几何口径不一致 —— 本脚本把\n"
                f"  POSCAR 投到 2D 找键, 规格层 (magnetic_orders.json) 写的是 3D 键长。")
        print(f"壳来源 ({'jfile' if args.jfile else 'shell-def'}): {origin}")
        print(f"{len(defs)} 个显式壳, 不做距离自动排序: "
              + ", ".join(f"{names[i]}@{shells_d[i]:.4f}A" for i in range(len(defs))))
    else:
        if args.J1 is None:
            raise SystemExit("没有 --jfile/--shell-def 时必须给 --J1")
        all_shells = shell_of(list(bonds.values()))
        shells = all_shells[:3]          # 旧路径只认最近 3 个壳 (J1/J2/J3)
        if len(all_shells) < 2:
            print("[warn] 只找到 1 个近邻壳, J2/J3 将不被使用", file=sys.stderr)
        if len(all_shells) < 3:
            print(f"[warn] 只找到 {len(all_shells)} 个近邻壳, J3 将不被使用 "
                  f"(若 J3 重要请加大 --nrange)", file=sys.stderr)
        if len(all_shells) > 3:
            print(f"[note] 距离上共有 {len(all_shells)} 个壳, 旧路径只用最近 3 个 (J1/J2/J3); "
                  f"要更多壳请用 --shell-def/--jfile", file=sys.stderr)
        names = ["J1", "J2", "J3"][:len(shells)]
        shells_d = [float(x) for x in shells]
        defs = [dict(name=names[i], d=shells_d[i], tol=1e-2, pair=None,
                     J_meV=(getattr(args, f"J{i + 1}", 0.0) or 0.0), J_mat=None)
                for i in range(len(shells_d))]
        per_shell = [[k for k, d in bonds.items() if abs(d - shells_d[i]) < 1e-2]
                     for i in range(len(shells_d))]
        origin = "距离自动排序 (旧路径)"

    jmats = [jmat_of(defs[i], i, args, s2) for i in range(len(defs))]

    interactions = []
    for si, keys in enumerate(per_shell):
        for k, l, i, j in keys:
            interactions.append(dict(b1=int(k), b2=int(l), offset=[int(i), int(j)],
                                     shell=names[si], dist=shells_d[si],
                                     J=jmats[si].tolist()))

    onsite = np.zeros((len(basis), 3, 3))
    onsite[:, 2, 2] = -args.MAE

    per_site = [len(keys) * 2 / len(basis) for keys in per_shell]
    for si in range(len(defs)):
        print(f"{names[si]} 壳: d = {shells_d[si]:.4f} A, {len(per_shell[si])} 条键/原胞 "
              f"(~{per_site[si]:.1f} 近邻/格点), J_mat diag = {np.diag(jmats[si])[0]:.4f} meV")
    print(f"各向异性: A_zz = {-args.MAE:.4f} meV | S = {S} | S^2 = {s2} | "
          f"共 {len(interactions)} 条相互作用")

    # 标量 J: 仅当 J_mat 各向同性时存在; 含 DMI/各向异性 -> NaN (本来就没有单一标量)。
    # 走 -m[0,0]/s2 反推会经 J*S^2 往返, 留下 3.8800000000000003 这类尾差 -> 圆整到 10 位。
    jres = {}
    for i, nm in enumerate(names):
        m = jmats[i]
        jres[nm] = (round(float(-m[0, 0] / s2), 10)
                    if (s2 and np.allclose(m, m[0, 0] * I3)) else float("nan"))

    # mc_2d.py 的 if "J1_meV" in meta 分支会一起读 J2_meV/J3_meV: 要么三个都在, 要么都不在
    jmeta = {}
    if "J1" in jres and np.isfinite(jres["J1"]):
        def _f(nm):
            v = jres.get(nm, float("nan"))
            return float(v) if np.isfinite(v) else 0.0
        jmeta = {"J1_meV": jres["J1"], "J2_meV": _f("J2"), "J3_meV": _f("J3")}

    model = dict(
        meta=dict(
            source_poscar=os.path.abspath(args.poscar),
            element="+".join(elements), elements=elements, spins=spins, S=S, S2=s2,
            MAE_meV=args.MAE, easy_axis="z", dim=2,
            shell_source=origin,
            jfile=(os.path.abspath(args.jfile) if args.jfile else None),
            shells_summary=" ".join(
                f"{names[i]}={shells_d[i]:.3f}A x{len(per_shell[i])}" for i in range(len(defs))),
            convention="H=-sum_k J_k*sum_<ij,k> S.S; engine J=-(J*S^2); A_zz=-MAE",
            **jmeta,
        ),
        primitive_vectors=prim.tolist(),
        basis_frac=basis.tolist(),
        onsite=onsite.tolist(),
        interactions=interactions,
        shells=dict(
            d_nn=shells_d[0], d_nnn=(shells_d[1] if len(shells_d) > 1 else None),
            d_3rd=(shells_d[2] if len(shells_d) > 2 else None),
            nn_per_site=per_site[0],
            nnn_per_site=(per_site[1] if len(per_site) > 1 else 0.0),
            nnnn_per_site=(per_site[2] if len(per_site) > 2 else 0.0),
            n_shells=len(defs),
            all=[dict(name=names[i], d=shells_d[i], tol=defs[i]["tol"],
                      n_bonds=len(per_shell[i]), per_site=per_site[i],
                      J_meV=(float(jres[names[i]]) if np.isfinite(jres[names[i]]) else None),
                      J_diag=float(np.diag(jmats[i])[0])) for i in range(len(defs))],
        ),
    )

    out = args.out or (os.path.splitext(os.path.basename(args.poscar))[0] + "_mc_model.json")
    with open(out, "w") as f:
        json.dump(model, f, indent=2)
    print(f"\n已写出 Monte Carlo 模型输入文件: {os.path.abspath(out)}")
    print(f"用 mc_run.py 运行:  python mc_run.py {os.path.basename(out)} --task mt")


if __name__ == "__main__":
    main()
