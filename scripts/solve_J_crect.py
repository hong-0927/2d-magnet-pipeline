#!/usr/bin/env python
"""solve_J_crect.py
从四个共线磁态 (fm / afm_c / afm_a / afm_b) 总能量反解 **C 心矩形** 晶格 J1/J2/J3。

## 这个晶格是什么 (CrOCl)
C 心矩形: 正交 (a != b, 90°), 但**每胞 2 个磁原子** —— Cr 在 (0,0,z) 与 (1/2,1/2,z')。
所以它既不是 `rectangular`(1 磁位/胞, CrSBr) 也不是 `square`(a=b)。CrOCl 单层实测:
a = 3.1995 A, b = 3.8603 A (JARVIS JVASP-27906, Pmmn No.59)。

距离壳按升序 (每磁位点的键数):
  J1 = 3.0347 A  x4   <- **C 心矢量 (1/2,1/2)**; 注意它连的是**同一胞内的两个 Cr**
  J2 = 3.1980 A  x2   <- 沿 a 轴 (同一亚晶格内)
  J3 = 3.8603 A  x2   <- 沿 b 轴 (同一亚晶格内)

⚠ 两个坑, 都会**静默**给出错答案:
  (1) J1 与 J2 只差 0.1632 A。`build_rectangular_cases.py` 的 `magnetic_shells(tol=0.2)`
      会把它们**合并成一个 J1** —— 两个物理上不同的交换路径算成一个, 反解照样收敛、
      残差照样是 0、不报任何错。所以本晶格必须用**显式壳距** (build_crect_cases.py 用
      tol=0.05), 且产物里的 tol 由 j_result 按「相邻壳间距的一半」自动推 (=0.0816)。
  (2) C 心矢量连胞内两个 basis ⇒ **basis 符号对一切 (n0,n1) 奇偶序不可见**:
      afm_c_a 与 afm_a 逐位相同, afm_c_diag 与 afm_diag 逐位相同。所以序集合里必须
      同时有 pure basis (afm_c) 和 pure 格点奇偶 (afm_a/afm_b), 但不能混着选。

## 为什么恰好是这四个序
设计矩阵 A 的第 k+1 列 = 该壳在各磁序下的**每胞系数和** c_k(σ)。实测 (4x4x1):
    fm      c1=-64 c2=-32 c3=-32      (per mag: -2, -1, -1)
    afm_c   c1=+64 c2=-32 c3=-32      (per mag: +2, -1, -1)
    afm_a   c1=  0 c2=+32 c3=-32
    afm_b   c1=  0 c2=-32 c3=+32
  * fm 与 afm_c 只在 c1 上反号 ⇒ 这一对唯一的作用是把 **J1 从 E0 里分出来**。
    没有 afm_c, J1 与 E0 就不可分 (这正是「必须支持多亚晶格」的物理必要性)。
  * afm_a / afm_b 把 c2 / c3 分别拉到 + 侧 ⇒ 分开 J2 与 J3。
  * 不可选 (fm, afm_c, p4_a, afm_c_p4): 该组合 **det(A)=0** (P2 侦察时实测), 因为
    一切周期 2 与周期 4 的格点奇偶序 c1 恒为 0, 两行成比例。
  * 与方形不同, 这里**不需要周期 4 的序**: 方形是 a=b 使 stripe_a/stripe_b 简并、
    且周期 2 的序 J3 与 E0 简并, 才被迫引入 p4_a; C 心矩形 a≠b, afm_a/afm_b 的
    (c2,c3) 已经独立, J3 与 E0 不再简并。
  教训 (换晶格必须重算 cond/det): 方形那组 (fm,stripe_a,p4_a,mix_ab) 的 cond=4.62,
  这里换成周期 2 的四个序反而 cond=64.02 —— 更大但仍远低于 1e4 闸门, 且 det=524288。

用法:
  python solve_J_crect.py <cases_dir>                          # 读 cases/{fm,afm_c,afm_a,afm_b}/OUTCAR
  python solve_J_crect.py <cases_dir> --demo --synthetic 5 -0.5 0.2   # 合成 J, 无 DFT 自检
  python solve_J_crect.py <cases_dir> --demo --ref CrOCl        # 从 literature_J.json 取锚点

⚠ CrOCl 的文献 J **未核实** (literature_J.json 里 J1/J2/J3 = null) ⇒ `--ref CrOCl` 会
  明确拒绝, 且 `--demo` 必须显式给 `--ref` 或 `--synthetic` 之一 —— 绝不用 0 或别的材料
  的值顶上 (WORKFLOW §63/§99 不静默降级)。

约定: H = -J1 Σ S·S - J2 Σ S·S - J3 Σ S·S, J>0 = 铁磁, meV。
"""
import argparse
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import j_result as JR  # noqa: E402

ORDERS = ["fm", "afm_c", "afm_a", "afm_b"]

# 本晶格家族的**定量**锚点 = 没有 (CrOCl 的文献 J 未核实) ⇒ LIT 故意留 None。
LIT = None
LIT_SOURCE = "(无: CrOCl 的文献 J 未核实, 见 literature_J.json 的 CrOCl/_note)"

LIT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "literature_J.json")

# fm 行给出的配位数 (每磁位点的无向键数), 这是纯几何量, 与超胞大小无关。
# 反解出的 J 必须与它对上; 对不上说明壳被合并/切错 (坑 1)。
EXPECT_Z_N = (4.0, 2.0, 2.0)


def read_toten(path):
    """取 OUTCAR 的总能量 (eV)。

    P3.17: 口径统一到 `dft_energy.py` —— 一律 `free energy TOTEN`、取**最后一次**出现、
    缺了就报错（不静默改用 energy(sigma->0)）。这里只转发, 不再各写一份。
    """
    from dft_energy import read_toten as _rt
    return _rt(path)


def guess_material(path):
    """从 cases 目录向上找 `*_pipeline` / `*_pipeline_v<N>` 并剥后缀 -> "CrOCl"。

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
            "  C 心矩形仍可进主链（几何/规格层可用）, 但 --demo 需要一个已核实的 J;"
            " 换成已核实的同晶型材料, 或用 --synthetic 显式给一组合成 J 做无 DFT 自检。")
    return dict(J1=float(ent["J1"]), J2=float(ent["J2"]), J3=float(ent["J3"])), \
        f"{key}: {ent.get('source', '(literature_J.json)')}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases")
    ap.add_argument("--demo", action="store_true",
                    help="用锚点 J 合成磁态能量自检（不需要 DFT）; 需配 --ref 或 --synthetic")
    ap.add_argument("--ref", default=None,
                    help="从 literature_J.json 取锚点的材料名（或直接给该 json 的路径）")
    ap.add_argument("--synthetic", nargs=3, type=float, default=None, metavar=("J1", "J2", "J3"),
                    help="显式给一组**合成** J (meV) 做无 DFT 自检; 不是文献值, 产物标 demo")
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
    if meta.get("lattice") not in (None, "rect_centered"):
        print(f"  [warn] magnetic_orders.json 的 lattice={meta['lattice']!r}, "
              f"本脚本期望 'rect_centered'（系数是按 C 心矩形算的吗?）")
    A = np.array([[1.0] + c[k] for k in ORDERS])
    cond = float(np.linalg.cond(A))

    # ---- 锚点解析: --ref 与 --synthetic 互斥, --demo 必须给其中之一 ------------
    if args.ref and args.synthetic:
        raise SystemExit("--ref 与 --synthetic 互斥: 一个是文献锚点, 一个是合成值, 不许混用")
    lit, lit_src = LIT, LIT_SOURCE
    if args.ref:
        lit, lit_src = load_ref(args.ref)
    elif args.synthetic:
        lit = dict(J1=args.synthetic[0], J2=args.synthetic[1], J3=args.synthetic[2])
        lit_src = "synthetic (显式给的合成 J, **不是文献值**)"
    if args.demo and lit is None:
        raise SystemExit(
            "本晶格的文献 J 未核实 ⇒ --demo 需要一个显式锚点, 二者之一:\n"
            f"  --ref <材料名>           从 literature_J.json 取已核实的文献 J\n"
            "  --synthetic J1 J2 J3     显式给合成 J 做无 DFT 自检 (例: --synthetic 5 -0.5 0.2)\n"
            "拒绝用 0 或别的材料的值顶上 (WORKFLOW §99 不静默降级)。")

    if args.demo:
        Jt = np.array([lit["J1"], lit["J2"], lit["J3"]]) * S2   # J~ = J*S^2
        E = {k: float(c[k][0] * Jt[0] + c[k][1] * Jt[1] + c[k][2] * Jt[2])
             for k in ORDERS}
        print("[demo] 锚点 J 来源: " + lit_src)
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
    print("   (J1 = 3.0347 A C 心矢量 x4 [连胞内两个 Cr], "
          "J2 = 3.1980 A 沿 a x2, J3 = 3.8603 A 沿 b x2)")

    # 配位数自检: 由 fm 行直接给出, 纯几何量 ⇒ 必须是 4/2/2。
    z = tuple(round(2.0 * (-c["fm"][i]) / n, 6) for i in range(3))
    ok_z = all(abs(z[i] - EXPECT_Z_N[i]) < 1e-9 for i in range(3))
    print(f"配位数 (由 fm 行导出, 每磁位点): J1={z[0]:g} J2={z[1]:g} J3={z[2]:g}"
          f"   期望 {EXPECT_Z_N}  ->  " + ("PASS" if ok_z else "FAIL"))
    if not ok_z:
        raise SystemExit(
            "配位数与 C 心矩形几何不符 (期望 4/2/2) —— 壳很可能被 magnetic_shells(tol=0.2) "
            "静默合并了 (J1 与 J2 只差 0.1632 A)。请用显式壳距重建 magnetic_orders.json。")

    print(f"拟合质量: cond(A) = {cond:.4g}   residual = {residual:.3g} eV")
    if cond > JR.COND_WARN:
        print(f"  [warn] cond(A) = {cond:.4g} > {JR.COND_WARN:.4g} —— 设计矩阵接近病态,"
              " 反解结果对能量噪声敏感")

    material = args.material or guess_material(args.cases)
    if args.demo and lit and _lit_key(material) != _lit_key(lit_src.split(":")[0]):
        print(f"  [note] 几何来自 {material!r}, 但 demo 用的锚点 J 来自 "
              f"{lit_src.split(':')[0]!r} —— 系数只吃几何, J 只吃材料; 这是无 DFT 自检,"
              " 产物里的 material 以 --material 为准")

    if not args.no_jresult:
        d_ang = meta.get("shells_angstrom") or {}
        shells = [dict(name=nm, d=float(d_ang[nm]))
                  for nm in ("J1", "J2", "J3") if nm in d_ang]
        if len(shells) != 3:
            raise SystemExit(f"magnetic_orders.json 缺 shells_angstrom 的 J1/J2/J3: {d_ang}")
        notes = [f"coefficients_per_cell source = {os.path.basename(args.cases)}",
                 "C 心矩形 (2 磁位/胞): J1 = C 心矢量, 连胞内两个 basis; J2 沿 a; J3 沿 b"]
        if args.synthetic:
            notes.append("demo 能量由**合成** J 生成 (非文献值): "
                         f"J1={lit['J1']}, J2={lit['J2']}, J3={lit['J3']} meV")
        elif args.demo:
            notes.append("demo 能量由锚点 J 生成: " + lit_src)
        obj = JR.build(
            lattice="rect_centered",
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
        # 往返误差: E 由锚点 J 经同一个 A 合成, 反解必须逐位还原。
        err = max(abs({"J1": J1, "J2": J2, "J3": J3}[k] - lit[k]) for k in ("J1", "J2", "J3"))
        ok = err < 1e-6 and ok_z and cond < 1e4
        print(f"\n往返误差 (J1,J2,J3 vs 锚点) = {err:.3g} meV  ->  " + ("PASS" if ok else "FAIL"))
        if not ok:
            raise SystemExit("--demo 自检 FAIL")


if __name__ == "__main__":
    main()
