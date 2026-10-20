#!/usr/bin/env python3
"""DMI 的 3x3 交换矩阵等价形式: 把 D 矢量转成 `poscar_to_mc3.py --shell-def` 能吃的 `J_mat`。

## 为什么需要它

反称交换（Dzyaloshinskii-Moriya）在 MC 模型里没有独立字段，但它**可以精确**写进逐壳 3x3 `J_mat`：

    物理形式:  H_DM = +D · (S_i x S_j)                    （本脚本默认约定, --convention plus）
    等价形式:  H_DM = S_i^T J_D S_j,   J_D = -[D]_x       （[D]_x v = D x v）

    [D]_x = [[0, -Dz, +Dy], [+Dz, 0, -Dx], [-Dy, +Dx, 0]]
    J_D   = [[0, +Dz, -Dy], [-Dz, 0, +Dx], [+Dy, -Dx, 0]]

推导: `S_i^T (-[D]_x) S_j = -(D x S_j)·S_i = -(S_j x S_i)·D = +D·(S_i x S_j)` ✓

`--convention minus` 对应 `H_DM = -D · (S_i x S_j)` ⇒ `J_D = +[D]_x`。**符号约定必须与 J 的符号约定
一起写在结论里**（本仓: `H = -J Σ S·S`，`J > 0 = 铁磁`，见 `skill/2d-magnet/SKILL.md`「J 符号」）。

## 与 `J_mat` 的既有约定一致

`poscar_to_mc3.py` 的模型约定（文件头行 8-11）: `E = Σ_<ij> S_i·J_mat·S_j + Σ_i S_i·A·S_i - gμ_B B·Σ S_i`，
且各向同性时 `J_mat = -(J·S²)·I₃` —— **S² 被吸收进 `J_mat`**。DMI 项与交换同为双线性，所以同样吸收:

    J_mat(壳) = -(J_iso · S²) · I₃  +  J_D · S²        （`--J` 省略时只给 DMI 部分）

## 陷阱（实测踩过，务必知道）

`mc_2d._energy()` 返回的是**每格点**能量（源码: `0.5 * Σ_有向 / N`），而物理模型是
`E_cell = Σ_<ij> S_i·J_mat·S_j`（每键**一次**）⇒ `E_cell = N_sites × mc_2d._energy(...)`。
写对拍脚本时把 `_energy` 当成键能量会得到**差 4 倍**的数（本脚本 `--selftest` 的第 4 条就是钉这个）。

## 用法

    python scripts/dmi_matrix.py --D 0.31 -0.77 0.55 --J 1.9 --S 1.5      # 人读
    python scripts/dmi_matrix.py --D 0.31 -0.77 0.55 --J 1.9 --json       # 给 --shell-def 用
    python scripts/dmi_matrix.py --selftest                               # 无 DFT 自检 (含与 mc_2d 对拍)

产物只是**矩阵**：`D` 本身仍必须来自文献或你自己的 DFT（本仓不提取 D，见 `knowledge/DMI_MTI_MEMO.md` §3.2）。
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)


def skew(D):
    """[D]_x : 满足 [D]_x v = D x v 的反称矩阵。"""
    Dx, Dy, Dz = (float(x) for x in D)
    return np.array([[0.0, -Dz, Dy],
                     [Dz, 0.0, -Dx],
                     [-Dy, Dx, 0.0]])


def dmi_matrix(D, S=1.0, convention="plus"):
    """D 矢量 -> 等价反称交换矩阵 J_D（已含 S² 吸收, 与 poscar_to_mc3 的 J_mat 同尺度）。"""
    if convention not in ("plus", "minus"):
        raise SystemExit(f"--convention 只能是 plus 或 minus, 收到 {convention!r}")
    J = -skew(D) if convention == "plus" else skew(D)
    return J * (float(S) ** 2)


def jmat_shell(D, J_iso=None, S=1.0, convention="plus"):
    """整壳 J_mat = 各向同性部分 + DMI 部分（都按 poscar_to_mc3 的 S² 吸收约定）。"""
    M = np.zeros((3, 3))
    if J_iso is not None:
        M += -(float(J_iso) * float(S) ** 2) * np.eye(3)
    M += dmi_matrix(D, S=S, convention=convention)
    return M


# --------------------------------------------------------------------------- selftest
def _selftest():
    import mc_2d  # 只在自检里用 (对拍能量约定)

    fails = []
    rng = np.random.default_rng(7)

    def rand_unit():
        v = rng.normal(size=3)
        return v / np.linalg.norm(v)

    def model_energy(Jmat, s0, s1):
        """单键 1x1 胞: mc_2d._energy 是每格点能量, 乘回 N 得键能量。"""
        model = {
            "primitive_vectors": np.eye(3) * 10.0,
            "basis_frac": [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]],
            "onsite": [[[0.0] * 3 for _ in range(3)] for _ in range(2)],
            "interactions": [{"b1": 0, "b2": 1, "offset": [0, 0], "J": np.asarray(Jmat).tolist()}],
        }
        spins, nb, J, onsite, _ = mc_2d.build_arrays(model, (1, 1))
        spins[0] = s0
        spins[1] = s1
        return mc_2d._energy(spins, nb, J, onsite) * spins.shape[0]

    D = np.array([0.31, -0.77, 0.55])

    # 1) 反称矩阵性质: J_D^T = -J_D, 对角 0
    M = dmi_matrix(D, S=1.0)
    if not np.allclose(M.T, -M, atol=1e-15) or not np.allclose(np.diag(M), 0.0, atol=1e-15):
        fails.append("J_D 不是反称矩阵 / 对角非零")

    # 2) 等价性 (plus 约定): S0^T J_D S1 == +D.(S0 x S1),  mc_2d 实测
    for _ in range(5):
        a, b = rand_unit(), rand_unit()
        got, want = model_energy(M, a, b), float(D @ np.cross(a, b))
        if abs(got - want) > 1e-12:
            fails.append(f"plus 约定对拍: got={got:.12f} want={want:.12f}")

    # 3) minus 约定给出相反符号
    Mm = dmi_matrix(D, S=1.0, convention="minus")
    a, b = rand_unit(), rand_unit()
    if abs(model_energy(Mm, a, b) + float(D @ np.cross(a, b))) > 1e-12:
        fails.append("minus 约定对拍失败")

    # 4) 与各向同性部分线性可分 + S² 只做整体缩放
    Jiso, S = 1.9, 1.5
    Jmix = jmat_shell(D, J_iso=Jiso, S=S)
    for _ in range(5):
        a, b = rand_unit(), rand_unit()
        got = model_energy(Jmix, a, b)
        want = -(Jiso * S ** 2) * float(a @ b) + (S ** 2) * float(D @ np.cross(a, b))
        if abs(got - want) > 1e-12:
            fails.append(f"混合壳对拍: got={got:.12f} want={want:.12f}")
    if not np.allclose(dmi_matrix(D, S=2.0), 4.0 * dmi_matrix(D, S=1.0), atol=1e-15):
        fails.append("S² 缩放不是纯倍数")

    # 5) 反向键用 J.T: 正反两个方向一致（键能量与 i<j 的取向无关）
    Janiso = np.array([[1.2, 0.4, -0.3], [-0.9, 0.7, 0.2], [0.5, -0.1, 1.1]])
    for _ in range(3):
        a, b = rand_unit(), rand_unit()
        got = model_energy(Janiso, a, b)
        want = 0.5 * (float(a @ Janiso @ b) + float(b @ Janiso.T @ a))
        if abs(got - want) > 1e-12:
            fails.append(f"双向一致性: got={got:.12f} want={want:.12f}")

    print(f"[dmi] 反称矩阵性质            : {'PASS' if len(fails) == 0 else '见下'}")
    print(f"[dmi] plus 约定 vs mc_2d 对拍 : {'PASS' if not any('plus' in f for f in fails) else 'FAIL'}")
    print(f"[dmi] minus 约定符号          : {'PASS' if not any('minus' in f for f in fails) else 'FAIL'}")
    print(f"[dmi] J_iso + J_D 线性可分    : {'PASS' if not any('混合' in f or '缩放' in f for f in fails) else 'FAIL'}")
    print(f"[dmi] 反向键 J.T 双向一致     : {'PASS' if not any('双向' in f for f in fails) else 'FAIL'}")
    if fails:
        print("\n".join("  " + f for f in fails))
        print("FAIL")
        return 1
    print("ALL PASS")
    return 0


# --------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="D 矢量 -> 3x3 J_mat (给 poscar_to_mc3.py --shell-def 用); 只做矩阵换算, 不产生 D")
    ap.add_argument("--D", nargs=3, type=float, metavar=("Dx", "Dy", "Dz"),
                    help="DMI 矢量 (meV)。必须来自文献或你自己的 DFT —— 本仓不提取 D")
    ap.add_argument("--J", type=float, default=None, help="同一壳的各向同性 J (meV), 省略则只输出 DMI 部分")
    ap.add_argument("--S", type=float, default=1.0, help="自旋量子数 (吸收进 J_mat, 默认 1.0)")
    ap.add_argument("--convention", choices=("plus", "minus"), default="plus",
                    help="plus: H_DM = +D.(Si x Sj) -> J_D = -[D]_x (默认); minus: 反号")
    ap.add_argument("--name", default="J1", help="壳名 (写进 JSON 片段, 默认 J1)")
    ap.add_argument("--d", type=float, default=None, help="壳距离 (Å), 写进 JSON 片段")
    ap.add_argument("--tol", type=float, default=None, help="壳容差 (Å), 写进 JSON 片段")
    ap.add_argument("--json", action="store_true", help="只输出 JSON 片段 (可直接贴进 --shell-def)")
    ap.add_argument("--out", default=None,
                    help="把 JSON 片段写到文件 (UTF-8)。**优先用它而不是 shell 的 `>`** —— "
                         "PowerShell 的 `python ... > x.json` 会写成 UTF-16LE+ BOM")
    ap.add_argument("--selftest", action="store_true", help="无 DFT 自检 (含与 mc_2d._energy 对拍)")
    args = ap.parse_args(argv)

    if args.selftest:
        return _selftest()
    if args.D is None:
        ap.error("需要 --D Dx Dy Dz（或 --selftest）")

    D = np.array(args.D, dtype=float)
    M = jmat_shell(D, J_iso=args.J, S=args.S, convention=args.convention)
    shell = {"name": args.name}
    if args.d is not None:
        shell["d"] = float(args.d)
    if args.tol is not None:
        shell["tol"] = float(args.tol)
    if args.J is not None:
        shell["J_meV"] = float(args.J)
    shell["J_mat"] = [[float(x) for x in row] for row in M]

    payload = {"shells": [shell]}
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
            f.write("\n")
        print(f"> 壳定义已写 {os.path.abspath(args.out)} (UTF-8) —— 直接给 "
              f"`poscar_to_mc3.py --shell-def {args.out}`")
        print(f"> 同一个 --S ({args.S:g}) 要一起传给 poscar_to_mc3: J_mat 路径不会再乘 S^2")
        return 0

    if args.json:
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return 0

    sign_txt = "+D.(Si x Sj)" if args.convention == "plus" else "-D.(Si x Sj)"
    print(f"D = ({D[0]:+.4f}, {D[1]:+.4f}, {D[2]:+.4f}) meV | |D| = {np.linalg.norm(D):.4f} meV | "
          f"S = {args.S:g} | H_DM = {sign_txt}")
    print(f"J_iso = {'(未给)' if args.J is None else f'{args.J:+.4f} meV'}")
    print("J_mat (3x3, 给 poscar_to_mc3.py --shell-def):")
    for row in M:
        print("  [ " + "  ".join(f"{x:+.6f}" for x in row) + " ]")
    print("> 各向异性/DMI 存在时没有单一标量 J（poscar_to_mc3 会把标量 J 记成 NaN），"
          "下游只读标量 J 的地方要能容忍。")
    print("> 本脚本只换算矩阵；D 必须来自文献或你自己的 DFT（见 knowledge/DMI_MTI_MEMO.md §3.2）。")
    return 0


if __name__ == "__main__":
    sys.exit(main())
