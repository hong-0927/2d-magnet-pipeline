#!/usr/bin/env python
"""solve_J_square.py
从四个共线磁态 (fm / stripe_a / p4_a / mix_ab) 总能量反解方形晶格 J1/J2/J3。

方形晶格 (a = b, 原胞 1 个磁原子), 距离壳按升序: J1 = a (近邻), J2 = a√2 (对角),
J3 = 2a (同轴向次近邻)。四个共线磁态线性无关, 恰好解 {e0, J1, J2, J3}。

为什么这四个态 (P1 的关键物理论证, 见 build_square_cases.py 注释):
  一切「周期 2」的序 (fm / g_type / stripe_a / stripe_b) 给出的 **J3 系数完全相同**
  (因为 C(2,0) = N^2 对周期 2 的序恒成立) ⇒ 只用周期 2 的序, J3 与 E0 无法分离。
  必须引入「周期 4」的序 p4_a = (-1)^floor(n0/2) 才能把 J3 单独解出来。这与 kagome
  的 q=0 简并同构, 只是成因不同。P0.1 规格层里 square 的 selected 就是这四个。
  副产物: a=b 时 stripe_a 与 stripe_b 严格简并, 只留一个代表。

用法:
  python solve_J_square.py <cases_dir>                 # 读 cases/{fm,stripe_a,p4_a,mix_ab}/OUTCAR
  python solve_J_square.py <cases_dir> --demo          # 用文献 J 生成的假能量做无 DFT 自检
  python solve_J_square.py <cases_dir> --demo --ref la2cuo4   # 从 literature_J.json 取锚点

约定: H = -J1 Σ S·S - J2 Σ S·S - J3 Σ S·S, J>0 = 铁磁, meV。

参考锚点 (方形格子族的定量来源; 原文约定 H = +J Σ_<ij> S·S, J>0 = 反铁磁):
  R. Coldea et al., "Spin Waves and Electronic Interactions in La2CuO4",
  PRL 86, 5377 (2001); arXiv:cond-mat/0006384
    295 K 只用 J,J' 自由拟合 (Z_c=1.18): J = 111.8 ± 4 meV (反铁磁), J' = -11.4 ± 3 meV (对角铁磁)
     10 K 同上:                          J = 104.1 ± 4 meV,           J' = -18   ± 3 meV
    Hubbard/环交换拟合给出 J'' = J_c/20 = 2 ± 0.5 meV (295 K, 反铁磁)
    T_N = 325 ± 5 K
  换算到本脚本约定 (符号取反): J1 = -111.8, J2 = +11.4, J3 = -2.0 meV。
  K2NiF4 只有 T_N = 97.1 K 是核实的实验锚点, 其 magnon 拟合 J 未核实 ⇒
  literature_J.json 里 J 留 null; 对它加 --ref 会被明确拒绝 (不静默降级)。
"""
import argparse
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import j_result as JR  # noqa: E402

ORDERS = ["fm", "stripe_a", "p4_a", "mix_ab"]

# 方形格子族的定量锚点 = La2CuO4 (Coldea et al. 2001), 按本仓库约定 (J>0 = 铁磁)
LIT = dict(J1=-111.8, J2=11.4, J3=-2.0)
# 同一份文献的低温值, 留作对照 (只用近邻+对角两壳自由拟合)
LIT_10K = dict(J1=-104.1, J2=18.0, J3=-2.0)
LIT_SOURCE = "La2CuO4: Coldea et al., PRL 86, 5377 (2001)  [295 K: J=111.8, J'=-11.4 meV + J''=-2 meV]"

LIT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "literature_J.json")


def read_toten(path):
    """取 OUTCAR 的总能量 (eV)。

    P3.17: 口径统一到 `dft_energy.py` —— 一律 `free energy TOTEN`、取**最后一次**出现、
    缺了就报错（不静默改用 energy(sigma->0)）。这里只转发, 不再各写一份。
    """
    from dft_energy import read_toten as _rt
    return _rt(path)


def guess_material(path):
    """从 cases 目录向上找 `*_pipeline` / `*_pipeline_v<N>` 并剥后缀 -> "La2CuO4"。

    不直接用 cases 的 basename: 那是 "cases" / "cases_j123" / "02_supercell",
    当材料名没有任何信息量。找不到就返回 None (产物里留 null 比写错好)。
    """
    p = os.path.abspath(path)
    for _ in range(6):
        m = re.match(r"^(?P<name>.+?)_pipeline(?:_v\d+)?$", os.path.basename(p))
        if m:
            return m.group("name")
        parent = os.path.dirname(p)
        if parent == p:
            break
        p = parent
    return None


def _lit_key(name):
    return re.sub(r"[^a-z0-9]", "", (name or "").lower())


def load_ref(ref):
    """从 literature_J.json 取某个材料的 J1/J2/J3 (meV)。

    找不到、或 J 为 null 一律 SystemExit —— 「规格层没有可靠 J 就当没有」,
    绝不用 0 或别的材料的值顶上 (WORKFLOW §99 不静默降级)。
    """
    src = ref if os.path.exists(ref) else LIT_FILE
    if not os.path.exists(src):
        raise SystemExit(f"--ref 指向的文献库不存在: {src}")
    with open(src, encoding="utf-8-sig") as f:
        lit = json.load(f)
    want = _lit_key(ref if src == LIT_FILE else os.path.splitext(os.path.basename(ref))[0])
    hit = None
    for k, v in lit.items():
        if not isinstance(v, dict):        # 顶层 _comment / _convention 是字符串
            continue
        if _lit_key(k) == want or _lit_key(v.get("formula", "")) == want:
            hit = (k, v)
            break
    if hit is None:
        raise SystemExit(f"literature_J.json 里没有材料 {ref!r}; 现有条目: {sorted(lit)}")
    key, ent = hit
    missing = [nm for nm in ("J1", "J2", "J3") if ent.get(nm) is None]
    if missing:
        raise SystemExit(
            f"文献库条目 {key!r} 的 {'/'.join(missing)} 是 null（未核实）⇒ 不能当锚点。\n"
            f"  note: {ent.get('_note', '(无)')}\n"
            "  该材料仍可进主链（几何/规格层可用）, 但 --demo 需要一个已核实的 J;"
            " 换成已核实的同晶型材料, 或显式接受这一点后自己给 J。")
    return dict(J1=float(ent["J1"]), J2=float(ent["J2"]), J3=float(ent["J3"])), \
        f"{key}: {ent.get('source', '(literature_J.json)')}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases")
    ap.add_argument("--demo", action="store_true",
                    help="用文献锚点 J 合成磁态能量自检（不需要 DFT）")
    ap.add_argument("--ref", default=None,
                    help="从 literature_J.json 取锚点的材料名（或直接给该 json 的路径）")
    ap.add_argument("--material", default=None, help="材料名 (写进 J_result.json)")
    ap.add_argument("--jfile-out", default=None,
                    help="J_result.json 输出路径 (默认 <cases>/J_result.json)")
    ap.add_argument("--no-jresult", action="store_true", help="不写统一产物 (只打印)")
    args = ap.parse_args()

    meta = json.load(open(os.path.join(args.cases, "magnetic_orders.json"),
                          encoding="utf-8-sig"))
    n = meta["n_mag"]
    S2 = meta["S2"]
    c = meta["coefficients_per_cell"]
    missing = [k for k in ORDERS if k not in c]
    if missing:
        raise SystemExit(f"magnetic_orders.json 缺磁序 {missing}; 现有 {sorted(c)}")
    A = np.array([[1.0] + c[k] for k in ORDERS])
    cond = float(np.linalg.cond(A))

    lit, lit_src = LIT, LIT_SOURCE
    if args.ref:
        lit, lit_src = load_ref(args.ref)

    if args.demo:
        Jt = np.array([lit["J1"], lit["J2"], lit["J3"]]) * S2   # J~ = J*S^2
        E = {k: float(c[k][0] * Jt[0] + c[k][1] * Jt[1] + c[k][2] * Jt[2])
             for k in ORDERS}
        print("[demo] 参考 J 来源: " + lit_src)
        print("[demo] 用该 J 生成的磁态能量 (eV/胞):",
              {k: f"{v:.6f}" for k, v in E.items()})
        E = {k: v / 1000.0 for k, v in E.items()}               # meV -> eV
    else:
        E = {}
        for nm in ORDERS:
            p = os.path.join(args.cases, nm, "OUTCAR")
            if os.path.exists(p):
                E[nm] = read_toten(p)
        if set(ORDERS) != set(E):
            raise SystemExit(f"需要 {'/'.join(ORDERS)} 四态 OUTCAR; 可加 --demo 自检")

    b = np.array([E[k] for k in ORDERS])
    e0, J1t, J2t, J3t = np.linalg.solve(A, b)

    J1 = J1t / S2 * 1000.0
    J2 = J2t / S2 * 1000.0
    J3 = J3t / S2 * 1000.0

    # 反解自洽性: 把解回代, 看与输入能量差多少 —— 差值大说明设计矩阵/能量不自洽
    residual = float(np.max(np.abs(A @ np.array([e0, J1t, J2t, J3t]) - b)))

    print(f"\n每胞系数: " + "  ".join(f"{k}{c[k]}" for k in ORDERS))
    print(f"能量 (eV/胞): " + "  ".join(f"{k}={E[k]:+.6f}" for k in ORDERS))
    print(f"J1~ = {J1t:+.6f} eV   J2~ = {J2t:+.6f} eV   J3~ = {J3t:+.6f} eV   (S^2 = {S2})")
    print(f"\nJ1 = {J1:+.4f} meV    J2 = {J2:+.4f} meV    J3 = {J3:+.4f} meV")
    print(f"   (J1 = a 近邻, J2 = a√2 对角, J3 = 2a 同轴向次近邻)")
    print(f"拟合质量: cond(A) = {cond:.4g}   residual = {residual:.3g} eV")
    if cond > JR.COND_WARN:
        print(f"  [warn] cond(A) = {cond:.4g} > {JR.COND_WARN:.4g} —— 设计矩阵接近病态,"
              " 反解结果对能量噪声敏感")

    material = args.material or guess_material(args.cases)
    if args.demo and _lit_key(material) != _lit_key(lit_src.split(":")[0]):
        print(f"  [note] 几何来自 {material!r}, 但 demo 用的锚点 J 来自 "
              f"{lit_src.split(':')[0]!r} —— 系数只吃几何, J 只吃材料; 这是无 DFT 自检,"
              " 产物里的 material 以 --material 为准")

    if not args.no_jresult:
        d_ang = meta.get("shells_angstrom") or {}
        shells = [dict(name=nm, d=float(d_ang[nm]))
                  for nm in ("J1", "J2", "J3") if nm in d_ang]
        if len(shells) != 3:
            raise SystemExit(f"magnetic_orders.json 缺 shells_angstrom 的 J1/J2/J3: {d_ang}")
        notes = [f"coefficients_per_cell source = {os.path.basename(args.cases)}"]
        if args.demo:
            notes.append("demo 能量由锚点 J 合成: " + lit_src)
        obj = JR.build(
            lattice="square",
            material=material,
            shells=shells, J=dict(J1=J1, J2=J2, J3=J3), E=E, E0=e0,
            cond=cond, residual=residual, spin=meta["spin"], n_mag=n,
            orders=ORDERS, A=A, supercell=meta.get("supercell"),
            source="demo" if args.demo else "vasp", z_from="fm",
            notes=notes,
        )
        out = args.jfile_out or os.path.join(args.cases, "J_result.json")
        warns = JR.write(out, obj)
        print(f"\n已写出统一产物: {os.path.abspath(out)}  (confidence={obj['confidence']})")
        for m in warns:
            print(f"  [warn] {m}")

    if args.demo:
        ok = all(abs({"J1": J1, "J2": J2, "J3": J3}[k] - lit[k]) < 0.02
                 for k in ("J1", "J2", "J3"))
        print(f"\n期望 ({lit_src.split(':')[0]}): J1={lit['J1']}, J2={lit['J2']}, "
              f"J3={lit['J3']} meV  ->  " + ("PASS" if ok else "FAIL"))


if __name__ == "__main__":
    main()
