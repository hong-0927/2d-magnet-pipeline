#!/usr/bin/env python
"""write_hamiltonian.py
从 poscar_to_mc3.py 生成的 Monte Carlo 模型 JSON 写出一份可读的 Heisenberg 哈密顿卡片 (Markdown)。

输出内容 (全部来自 model JSON, 无臆造):
  - 符号式 H (物理约定, 用 meta 的 J1/J2/J3 + MAE)
  - 模型参数表 (element / S / S^2 / 维度 / 易轴 / MAE / J)
  - 引擎 vs 物理约定说明 (经典 |S|=1, J_mat = -(J·S^2) I3, A_zz = -MAE)
  - 近邻壳摘要表 (shell / J / d / 配位数 / 每胞键数)
  - 逐键明细表 (b1 / b2 / offset / shell / dist / J 张量)
  - 晶格 (primitive_vectors + basis_frac)

用法:
  python write_hamiltonian.py mc_model.json [-o hamiltonian.md] [--label Cr]
"""

import argparse
import json
import os

import numpy as np


def fmt_num(x, nd=4):
    if x is None:
        return "-"
    if abs(x) < 1e-12:
        return "0"
    return f"{x:.{nd}f}".rstrip("0").rstrip(".")


def jmatrix_repr(J):
    """J (3x3) -> 简短可读表示。各向同性 -> 标量; 否则 -> 完整 3x3。"""
    a = np.asarray(J, dtype=float)
    diag = a[0, 0]
    off = a.copy()
    np.fill_diagonal(off, 0.0)
    isotropic = (np.allclose(off, 0.0, atol=1e-12)
                 and np.allclose(np.diag(a), diag, atol=1e-12))
    if isotropic:
        return f"{fmt_num(diag)}·I₃"
    rows = []
    for r in range(3):
        rows.append("[" + ", ".join(fmt_num(x, 3) for x in a[r]) + "]")
    return " ; ".join(rows)


def load_shells(model):
    meta = model.get("meta", {})
    shells = model.get("shells", {})
    return [
        ("J1", meta.get("J1_meV"), shells.get("d_nn"), shells.get("nn_per_site")),
        ("J2", meta.get("J2_meV", 0.0), shells.get("d_nnn"), shells.get("nnn_per_site")),
        ("J3", meta.get("J3_meV", 0.0), shells.get("d_3rd"), shells.get("nnnn_per_site")),
    ]


def symbolic_terms(shell_rows, mae):
    """按非零 J 生成符号式项列表; 返回 (terms, note)。"""
    terms = []
    for name, J, _, _ in shell_rows:
        if J is None or abs(J) == 0.0:
            continue
        terms.append((name, J))
    if mae is not None and abs(mae) > 0.0:
        terms.append(("MAE", mae))
    return terms


def build_markdown(model, label):
    meta = model.get("meta", {})
    element = meta.get("element", label or "?")
    S = meta.get("S")
    S2 = meta.get("S2")
    J1 = meta.get("J1_meV")
    J2 = meta.get("J2_meV", 0.0)
    J3 = meta.get("J3_meV", 0.0)
    mae = meta.get("MAE_meV", 0.0)
    easy_axis = meta.get("easy_axis", "z")
    dim = meta.get("dim", 2)
    source = meta.get("source_poscar", "-")

    prim = model.get("primitive_vectors", [])
    basis = model.get("basis_frac", [])
    interactions = model.get("interactions", [])
    shell_rows = load_shells(model)
    terms = symbolic_terms(shell_rows, mae)

    L = []
    L.append(f"# {label or element} 单层 Heisenberg 哈密顿量\n")
    L.append("> 本卡片由 `write_hamiltonian.py` 从 Monte Carlo 模型 JSON 生成，全部数值源自模型文件，无额外假设。\n")

    # ---- 符号式 ----
    L.append("## 符号式哈密顿量\n")
    if terms:
        parts = []
        for name, J in terms:
            if name == "MAE":
                parts.append(f"- MAE · Σ_i (S_i^{{z}})^2")
            elif name == "J1":
                parts.append(f"- J1 · Σ_{{⟨ij⟩}} S_i·S_j")
            elif name == "J2":
                parts.append(f"- J2 · Σ_{{⟨⟨ij⟩⟩}} S_i·S_j")
            else:
                parts.append(f"- J3 · Σ_{{⟨⟨⟨ij⟩⟩⟩}} S_i·S_j")
        L.append("```text")
        L.append("H = " + "\n      ".join(parts))
        L.append("```")
        L.append("")
    else:
        L.append("(无可用的非零交换项)\n")

    # ---- 参数表 ----
    L.append("## 模型参数\n")
    L.append("| 项 | 值 |")
    L.append("|---|---|")
    L.append(f"| 磁性元素 | {element} |")
    L.append(f"| 自旋 S / S² | {fmt_num(S)} / {fmt_num(S2)} |")
    L.append(f"| 维度 | {dim}D |")
    L.append(f"| 易轴 | {easy_axis} |")
    L.append(f"| MAE (meV/离子) | {fmt_num(mae)} |")
    L.append(f"| J1 (meV) | {fmt_num(J1)} |")
    L.append(f"| J2 (meV) | {fmt_num(J2)} |")
    L.append(f"| J3 (meV) | {fmt_num(J3)} |")
    L.append("")
    L.append("符号约定: `H = -J Σ S·S`，**J > 0 = 铁磁，J < 0 = 反铁磁**。")
    if ".clusters" in str(source):
        L.append("")
        L.append("> 注: 本卡片为**簇间有效三角模型**（每磁性簇 = 1 个有效自旋）。真实簇磁体为两层模型"
                 " `H = -J_intra Σ_簇内 S·S - Σ_k J_k Σ_第k壳簇间 S·S`；本卡片只含其中簇间部分，"
                 " J_intra 需另由 DFT/文献给出。")
    L.append("")

    # ---- 引擎约定说明 ----
    L.append("## 引擎 ↔ 物理约定\n")
    L.append("- 经典 Heisenberg，单位自旋 |S|=1；交换张量 `J_mat = -(J·S²)·I₃`（把量子自旋尺度吸收进耦合）。")
    L.append("- 单离子各向异性 `A_zz = -MAE`（易轴 z，MAE>0 时 (S^z)² 被压低、有利于自旋沿 z）。")
    L.append("")

    # ---- 近邻壳摘要 ----
    L.append("## 近邻壳摘要\n")
    L.append("| 壳 | J (meV) | d (Å) | 配位数/位点 |")
    L.append("|---|---|---|---|")
    for name, J, d, coord in shell_rows:
        if d is None:
            continue
        L.append(f"| {name} | {fmt_num(J)} | {fmt_num(d, 4)} | {fmt_num(coord)} |")
    L.append("")

    # ---- 逐键明细 ----
    L.append("## 逐键明细\n")
    L.append("| # | b1 | b2 | offset (i,j) | 壳 | d (Å) | J 张量 (meV, 引擎) |")
    L.append("|---|---|---|---|---|---|---|")
    for idx, it in enumerate(interactions, 1):
        b1 = it.get("b1")
        b2 = it.get("b2")
        off = it.get("offset")
        shell = it.get("shell")
        dist = it.get("dist")
        Jm = it.get("J")
        off_s = "(" + ", ".join(str(x) for x in off) + ")" if off is not None else "-"
        L.append(f"| {idx} | {b1} | {b2} | {off_s} | {shell} | {fmt_num(dist, 4)} | {jmatrix_repr(Jm)} |")
    L.append("")

    # ---- 晶格 ----
    L.append("## 晶格\n")
    L.append(f"- 原胞矢量 (Å): `{prim}`")
    L.append(f"- 基矢位点 (分数坐标): `{basis}`")
    L.append(f"- 源 POSCAR: `{source}`")
    L.append("")
    L.append("> 局限: 纯各向同性 Heisenberg + 单离子 MAE；二维 Mermin-Wagner 下临界温度为有限尺寸伪相变估计。")
    L.append("")

    return "\n".join(L)


def main():
    ap = argparse.ArgumentParser(description="模型 JSON -> Heisenberg 哈密顿 Markdown 卡片")
    ap.add_argument("model")
    ap.add_argument("-o", "--out", default=None)
    ap.add_argument("--label", default=None, help="标题用化学式/元素 (默认取 meta.element)")
    args = ap.parse_args()

    with open(args.model, encoding="utf-8") as f:
        model = json.load(f)

    label = args.label or model.get("meta", {}).get("element")
    md = build_markdown(model, label)

    out = args.out or os.path.join(os.path.dirname(os.path.abspath(args.model)), "hamiltonian.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write(md)
    print(f"已写出哈密顿卡片: {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
