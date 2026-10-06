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

参数:
  --element      磁性元素符号 (如 Cr)
  --S            自旋大小 S (Cr3+ d3 -> 1.5)
  --J1           最近邻交换 J1 (meV, H=-J sum S.S, J>0 铁磁)
  --J2           次近邻交换 J2 (meV), 可省略
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
    with open(path) as f:
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


# ----------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description="POSCAR + J1/J2/MAE -> Monte Carlo 模型输入文件")
    ap.add_argument("poscar")
    ap.add_argument("--element", required=True, help="磁性元素符号, 如 Cr")
    ap.add_argument("--S", type=float, required=True, help="自旋大小, 如 Cr3+ = 1.5")
    ap.add_argument("--J1", type=float, required=True, help="最近邻 J1 (meV)")
    ap.add_argument("--J2", type=float, default=0.0, help="次近邻 J2 (meV)")
    ap.add_argument("--J3", type=float, default=0.0, help="第三近邻 J3 (meV)")
    ap.add_argument("--MAE", type=float, default=0.0, help="单离子各向异性 (meV/离子, 易轴 z)")
    ap.add_argument("--no-s2", action="store_true", help="J1/J2 已是经典耦合, 不再乘 S^2")
    ap.add_argument("--nrange", type=int, default=2)
    ap.add_argument("-o", "--out", default=None)
    args = ap.parse_args()

    pos = parse_poscar(args.poscar)
    cell3 = pos["lattice"]
    if max(abs(cell3[2, 0]), abs(cell3[2, 1]), abs(cell3[0, 2]), abs(cell3[1, 2])) > 0.1:
        print("[warn] 晶格 c 似乎不在 z 方向, 2D 投影可能不准", file=sys.stderr)
    cell2 = np.array([cell3[0, :2], cell3[1, :2]])

    frac = get_magnetic_sites(pos, args.element)
    n_raw = len(frac)

    prim, basis, reduced = reduce_to_primitive(frac, cell2)
    print(f"磁性元素 {args.element}: {n_raw} 个/输入胞 -> 约化后 {len(basis)} 个/原胞 "
          f"({'已约化' if reduced else '未约化'})")

    bonds = find_bonds(basis, prim, nrange=args.nrange)
    shells = shell_of(list(bonds.values()))
    if len(shells) < 1:
        raise SystemExit("找不到任何近邻键")
    if len(shells) < 2:
        print("[warn] 只找到 1 个近邻壳, J2/J3 将不被使用", file=sys.stderr)
    if len(shells) < 3:
        print("[warn] 只找到 %d 个近邻壳, J3 将不被使用 (若 J3 重要请加大 --nrange)" % len(shells),
              file=sys.stderr)

    d1 = shells[0]
    d2 = shells[1] if len(shells) > 1 else None
    d3 = shells[2] if len(shells) > 2 else None
    shell1 = [k for k, d in bonds.items() if abs(d - d1) < 1e-2]
    shell2 = [k for k, d in bonds.items() if d2 is not None and abs(d - d2) < 1e-2]
    shell3 = [k for k, d in bonds.items() if d3 is not None and abs(d - d3) < 1e-2]

    s2 = 1.0 if args.no_s2 else args.S ** 2
    jmat1 = (-args.J1 * s2) * I3
    jmat2 = (-args.J2 * s2) * I3
    jmat3 = (-args.J3 * s2) * I3

    interactions = []
    for k, l, i, j in shell1:
        interactions.append(dict(b1=int(k), b2=int(l), offset=[int(i), int(j)],
                                 shell="J1", dist=d1, J=jmat1.tolist()))
    for k, l, i, j in shell2:
        interactions.append(dict(b1=int(k), b2=int(l), offset=[int(i), int(j)],
                                 shell="J2", dist=d2, J=jmat2.tolist()))
    for k, l, i, j in shell3:
        interactions.append(dict(b1=int(k), b2=int(l), offset=[int(i), int(j)],
                                 shell="J3", dist=d3, J=jmat3.tolist()))

    onsite = np.zeros((len(basis), 3, 3))
    onsite[:, 2, 2] = -args.MAE

    nn_per_site = len(shell1) * 2 / len(basis)
    nnn_per_site = len(shell2) * 2 / len(basis) if shell2 else 0.0
    nnnn_per_site = len(shell3) * 2 / len(basis) if shell3 else 0.0
    print(f"J1 壳: d = {d1:.4f} A, {len(shell1)} 条键/原胞 (~{nn_per_site:.1f} 近邻/格点)")
    if d2 is not None:
        print(f"J2 壳: d = {d2:.4f} A, {len(shell2)} 条键/原胞 (~{nnn_per_site:.1f} 近邻/格点)")
    if d3 is not None:
        print(f"J3 壳: d = {d3:.4f} A, {len(shell3)} 条键/原胞 (~{nnnn_per_site:.1f} 近邻/格点)")
    print(f"引擎耦合: J1_mat diag = {np.diag(jmat1)[0]:.4f} meV, "
          f"J2_mat diag = {np.diag(jmat2)[0]:.4f} meV, "
          f"J3_mat diag = {np.diag(jmat3)[0]:.4f} meV, A_zz = {-args.MAE:.4f} meV")

    model = dict(
        meta=dict(
            source_poscar=os.path.abspath(args.poscar),
            element=args.element, S=args.S, S2=s2,
            J1_meV=args.J1, J2_meV=args.J2, J3_meV=args.J3, MAE_meV=args.MAE,
            easy_axis="z", dim=2,
            convention="H=-J1*sum_<ij>S.S - J2*sum_<<ij>>S.S - J3*sum_<<<ij>>>S.S; "
                       "engine J=-(J*S^2); A_zz=-MAE",
        ),
        primitive_vectors=prim.tolist(),
        basis_frac=basis.tolist(),
        onsite=onsite.tolist(),
        interactions=interactions,
        shells=dict(d_nn=d1, d_nnn=d2, d_3rd=d3, nn_per_site=nn_per_site,
                    nnn_per_site=nnn_per_site, nnnn_per_site=nnnn_per_site),
    )

    out = args.out or (os.path.splitext(os.path.basename(args.poscar))[0] + "_mc_model.json")
    with open(out, "w") as f:
        json.dump(model, f, indent=2)
    print(f"\n已写出 Monte Carlo 模型输入文件: {os.path.abspath(out)}")
    print(f"用 mc_run.py 运行:  python mc_run.py {os.path.basename(out)} --task mt")


if __name__ == "__main__":
    main()
