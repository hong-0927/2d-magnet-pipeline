#!/usr/bin/env python
"""solve_J_kagome.py
从四个**共线**磁态 (fm / q0_ud / sq3_updntri / sq3_alt) 总能量反解 **kagome** 晶格 J1/J2/J3。

## 为什么是共线四态
J 反解只要求设计矩阵 [1,c1,c2,c3] 满秩; 与磁序是否共线无关。旧的组合
(fm/q0_120/sq3_updntri/sq3_alt) 里 q0_120 是面内 120° 非共线序, 需要 LNONCOLLINEAR。
换成共线的 q0_ud 后 rank 仍然是 4、cond 同量级 (实测 101 vs 98), 但整个 DFT 链路
只需要标量 MAGMOM —— 这就是默认走共线的原因。120° 序仍可由
`build_kagome_cases.py --include-noncollinear` 生成作对照, 但不再是反解的前提。

## 这个晶格是什么 (Fe3Sn)
kagome: 六角 (a=b, 120°), **每胞 3 个磁原子** (Fe 在 (0,1/2), (1/2,0), (1/2,1/2)),
非磁 Sn 在 (0,0)。每磁位点的键数 (纯几何, 与超胞无关):
  J1 = a/2      x4   <- 角共享三角形内, **同一三角形**上的两个近邻 (两条不等价键 J1 / J1')
  J2 = a√3/2    x4   <- 跨六角环的近邻 (与 J1 不等价, 经典 kagome 的 J1 vs J1')
  J3 = a        x6   <- 同亚晶格的对角 (再下一壳)
  第四壳 = a√7/2 ≈ 1.3229a (x8) 本脚本不采用。

## ⚠ 两个坑, 都会**静默**给出错答案
  (1) **q=0 一族对 J1/J2 是简并的**: fm / q0_ud / q0_duu (以及 q0_120) 的系数列满足
      c1 == c2 且 c3 相同 ⇒ 该子块 rank=2 (实测 cond=2.3e34, 即数值上完全奇异)。
      只用 q=0 序反解, 你会得到一个**看似正常**的 J1 与 J2 —— 它们只是
      (J1+J2)/2 的任意分配。必须引入 **√3×√3** 序 才能分离:
      sq3_alt 给出 c1 != c2 (分开 J1/J2 两列), sq3_updntri 的 c3 与 q=0 一族不同
      (把 J3 与 E0 分开)。本脚本对 rank<4 **直接报错**, 不做静默近似 (§115)。
  (2) **J1 与 J1' 同距离**: 角共享三角形的两条最近邻键距离都是 a/2, 但分属不同亚晶格对。
      本脚本 (与 builder) 的 3 壳模型把它们并成**一条** J1 ⇒ 反解出的 J1 是两者的
      **平均**。要把 J1 与 J1' 分开, 需要降对称的超胞 (P4.2 待做)。
      另外 `poscar_to_mc3 --jfile` 必须用 `pair=` 显式指定 (旧的 `shells[:3]` 自动路径
      会把同距离的两条并成一条, 见 poscar_to_mc3.py 的壳分配逻辑)。

## 为什么恰好是这四个序 (3x3x1 实测每胞系数)
    fm            c1=-54 c2=-54 c3=-81   (per mag: -2.00, -2.00, -3.00)
    q0_ud         c1=+18 c2=+18 c3=-81   (per mag: +0.67, +0.67, -3.00)
    sq3_updntri   c1= -6 c2= -6 c3=+27   (per mag: -0.22, -0.22, +1.00)
    sq3_alt       c1= -6 c2=+18 c3=+27   (per mag: -0.22, +0.67, +1.00)
  * rank=4, cond(A)≈101, |det|=1.87e5 —— 远低于 1e4 闸门。
  * q=0 一族单独 rank=2 (上面坑 1); √3×√3 两序提供两条独立信息 ⇒ 才可解。
  * 换序/换晶格必须**重算** cond/det (方形那组在 kagome 上 rank 会掉到 3)。

用法:
  python solve_J_kagome.py <cases_dir>                              # 读 cases/<序>/OUTCAR
  python solve_J_kagome.py <cases_dir> --demo --synthetic 3 -0.4 0.15   # 合成 J, 无 DFT 自检
  python solve_J_kagome.py <cases_dir> --demo --ref Fe3Sn              # 从 literature_J.json 取锚点
  python solve_J_kagome.py <cases_dir> --demo --synthetic 1 1 1 --jfile-out /tmp/J.json

⚠ Fe3Sn 的文献 J **未核实** (literature_J.json 里 Fe3Sn 条目已有, 但 J1/J2/J3 = null) ⇒ `--ref Fe3Sn` 会明确拒绝,
  且 `--demo` 必须显式给 `--ref` 或 `--synthetic` 之一 —— 绝不用 0 或别的材料的值顶上
  (WORKFLOW §63/§99 不静默降级)。

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

# 与 build_kagome_cases.py 的 SELECTED 必须一致 (builder 把 selected_orders 写进 json)
ORDERS = ["fm", "q0_ud", "sq3_updntri", "sq3_alt"]
Q0_FAMILY = ["fm", "q0_120", "q0_ud", "q0_duu"]

# fm 行给出的配位数 (每磁位点的无向键数): kagome 的几何指纹 4/4/6。
# 反解出的 J 必须与它对上; 对不上说明壳被合并/切错或 json 不是 kagome 的产物。
EXPECT_Z_N = (4.0, 4.0, 6.0)

# 本晶格家族的**定量**锚点 = 没有 (Fe3Sn 的文献 J 未核实) ⇒ LIT 故意留 None。
LIT = None
LIT_SOURCE = "(无: Fe3Sn 的文献 J 未核实, literature_J.json 里该条目的 J1/J2/J3 = null)"

LIT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "literature_J.json")


def read_toten(path):
    """取 OUTCAR 的总能量 (eV)。

    P3.17: 口径统一到 `dft_energy.py` —— 一律 `free energy TOTEN`、取**最后一次**出现、
    缺了就报错（不静默改用 energy(sigma->0)）。这里只转发, 不再各写一份。
    """
    from dft_energy import read_toten as _rt
    return _rt(path)


def guess_material(path):
    """从 cases 目录向上找 `*_pipeline` / `*_pipeline_v<N>` 并剥后缀 -> "Fe3Sn"。

    不直接用 cases 的 basename (那是 "cases"), 找不到就返回 None (产物里留 null 比写错好)。
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
            "  kagome 仍可进主链（几何/规格层可用）, 但 --demo 需要一个已核实的 J;"
            " 换成已核实的同晶型材料, 或用 --synthetic 显式给一组合成 J 做无 DFT 自检。")
    return dict(J1=float(ent["J1"]), J2=float(ent["J2"]), J3=float(ent["J3"])), \
        f"{key}: {ent.get('source', '(literature_J.json)')}"


def check_selected(meta):
    """校验 selected_orders 与脚本期望的四个序一致 (防 builder/solver 漂移)。"""
    sel = meta.get("selected_orders")
    if sel is None:
        print(f"  [note] magnetic_orders.json 没有 selected_orders, 按脚本内置的 {ORDERS} 反解")
        return
    if list(sel) != ORDERS:
        raise SystemExit(
            f"magnetic_orders.json 的 selected_orders={list(sel)} 与本脚本期望的 {ORDERS} 不一致。\n"
            "  kagome 的可解组合不是随便四个序: 必须含 √3×√3 序才能把 J1/J2 (q=0 简并) 分开,\n"
            "  且默认四态全部共线 (q0_ud 取代旧的 q0_120 非共线序)。\n"
            "  要么用 build_kagome_cases.py 重新生成, 要么显式改本脚本的 ORDERS 并**重算 cond/rank**。")


def read_ingest_manifest(cases, orders=None):
    """读 ingest_dft.py 写下的 cases/ingest_manifest.json (收回凭证)。

    返回 (prov, problem): prov 为 None 时 problem 说明为什么没有凭证。
    **凭证只证明「OUTCAR 过了格式/收敛/结构/口径闸门」, 不证明 J/T_N 的物理正确性。**
    """
    p = os.path.join(cases, "ingest_manifest.json")
    if not os.path.exists(p):
        return None, (f"{p} 不存在 —— 这批 OUTCAR 不是经 ingest_dft.py 收回的 "
                      "(手工 cp 的 OUTCAR 没有收敛/口径/结构/防旧文件闸门, "
                      "未收敛或混口径的能量会让 J 看起来正常却是错的)")
    try:
        with open(p, encoding="utf-8-sig") as f:
            man = json.load(f)
    except Exception as e:                                     # noqa: BLE001
        return None, f"{p} 读不动: {type(e).__name__}: {e}"
    if man.get("schema") != "dft_ingest/v1":
        return None, f"{p} 的 schema={man.get('schema')!r}, 本脚本认 'dft_ingest/v1'"
    o = man.get("orders") or {}
    want = list(orders or ORDERS)
    missing = [k for k in want if k not in o]
    if missing:
        return None, f"{p} 里没有这些序的记录: {missing}"
    files = {os.path.basename((f.get("outcar") or "")): f for f in (man.get("files") or [])}
    prov = {
        "schema": man.get("schema"),
        "generated_utc": man.get("generated_utc"),
        "source_runs": man.get("source_runs"),
        "energy_convention": man.get("energy_convention"),
        "lenient": [c["name"] for c in (man.get("checks") or []) if not c.get("hard")],
        "orders": {k: {kk: o[k].get(kk) for kk in
                       ("energy_eV", "n_ionic_steps", "scf_converged", "n_ions",
                        "n_kpts", "encut", "noncollinear", "lsorbit", "vasp_version")}
                   for k in want},
        "outcar_sha256": {k: (files.get("OUTCAR", {}).get("sha256") if len(want) == 1 else None)
                          for k in want},
        "manifest": os.path.abspath(p),
        "note": ("凭证只证明 OUTCAR 过了格式/收敛/结构/口径闸门, "
                 "**不证明 J 与 T_C 的物理正确性**"),
    }
    # 逐序 sha256 存在 manifest 的 files 列表里 (按 outcar 相对路径区分序)
    sha = {}
    for f in (man.get("files") or []):
        rel = (f.get("outcar") or "").replace("\\", "/")
        if rel:
            sha[rel.split("/")[0] if "/" in rel else rel] = f.get("sha256")
    if sha:
        prov["outcar_sha256"] = sha
    return prov, None


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
    ap.add_argument("--no-provenance-check", action="store_true",
                    help="放行没有 cases/ingest_manifest.json 的 OUTCAR (默认拒绝: 手工 cp 的 "
                         "OUTCAR 没有收敛/口径/结构/防旧文件闸门); 放行后 J_result.json 会**没有**收回凭证")
    args = ap.parse_args()

    meta = json.load(open(os.path.join(args.cases, "magnetic_orders.json"),
                          encoding="utf-8-sig"))
    n = meta["n_mag"]
    S2 = meta["S2"]
    c = meta["coefficients_per_cell"]
    missing = [k for k in ORDERS if k not in c]
    if missing:
        raise SystemExit(f"magnetic_orders.json 缺磁序 {missing}; 现有 {sorted(c)}")
    check_selected(meta)
    if meta.get("lattice") not in (None, "kagome"):
        print(f"  [warn] magnetic_orders.json 的 lattice={meta['lattice']!r}, "
              f"本脚本期望 'kagome'（系数是按 kagome 算的吗?）")

    A = np.array([[1.0] + c[k] for k in ORDERS])
    rk = int(np.linalg.matrix_rank(A))
    cond = float(np.linalg.cond(A))
    det = float(abs(np.linalg.det(A)))

    # ---- 可解性硬闸门 (§115 不静默降级): rank 与 cond 不达标就**不**给 J ----------
    if rk < len(ORDERS):
        q0_only = [k for k in ORDERS if k in Q0_FAMILY]
        raise SystemExit(
            f"设计矩阵 rank={rk} < {len(ORDERS)} ⇒ 这几个磁序**无法唯一确定** J1/J2/J3, 拒绝给 J。\n"
            f"  q=0 一族 {Q0_FAMILY} 满足 c1 == c2 且 c3 相同 ⇒ 其子块 rank=2 (cond≈2.3e34),\n"
            "  J1 与 J2 只能靠 √3×√3 序 (sq3_alt 的 c1 != c2) 分开。\n"
            f"  当前选中 {ORDERS} (其中 q=0 成员: {q0_only})。\n"
            "  修法: 用 build_kagome_cases.py 生成 (它内置 selected_orders), 或换含 √3×√3 的组合。\n"
            "  绝不允许对 rank 不足的矩阵求解后把结果当 J 报出去 (NEW_TYPES_PLAN §143)。")
    if cond > 1e4 or det < 1e-6:
        raise SystemExit(
            f"设计矩阵病态: cond={cond:.4g} (>1e4) 或 |det|={det:.3g} (<1e-6) ⇒ 拒绝给 J。\n"
            "  反解对能量噪声的放大倍数就是这个 cond; 请换壳/换序后**重算** cond。")

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
            "  --ref <材料名>           从 literature_J.json 取已核实的文献 J\n"
            "  --synthetic J1 J2 J3     显式给合成 J 做无 DFT 自检 (例: --synthetic 3 -0.4 0.15)\n"
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
    print("   (J1 = a/2 角共享三角形内 x4 [J1 与 J1' 的同距离平均], "
          "J2 = a√3/2 跨六角 x4, J3 = a 同亚晶格对角 x6)")
    print("   [note] J1 是三角形两条不等价最近邻键 (J1/J1') 的**平均**; 分开它们需要"
          " 降对称超胞 (P4.2)。")

    # 配位数自检: 由 fm 行直接给出, 纯几何量 ⇒ 必须是 4/4/6 (kagome 指纹)。
    z = tuple(round(2.0 * (-c["fm"][i]) / n, 6) for i in range(3))
    ok_z = all(abs(z[i] - EXPECT_Z_N[i]) < 1e-9 for i in range(3))
    print(f"配位数 (由 fm 行导出, 每磁位点): J1={z[0]:g} J2={z[1]:g} J3={z[2]:g}"
          f"   期望 {EXPECT_Z_N}  ->  " + ("PASS" if ok_z else "FAIL"))
    if not ok_z:
        raise SystemExit(
            "配位数与 kagome 几何不符 (期望 4/4/6) —— 壳很可能被合并/切错, 或这份 json "
            "不是 kagome 的产物。请用 build_kagome_cases.py 显式壳距重建 magnetic_orders.json。")

    print(f"可解性: rank(A) = {rk}/{len(ORDERS)}   cond(A) = {cond:.4g}   |det| = {det:.3g}")
    if cond > JR.COND_WARN:
        print(f"  [warn] cond(A) = {cond:.4g} > {JR.COND_WARN:.4g} —— 设计矩阵接近病态,"
              " 反解结果对能量噪声敏感")
    print(f"拟合质量: residual = {residual:.3g} eV")

    material = args.material or guess_material(args.cases)
    if args.demo and lit and material and _lit_key(material) != _lit_key(lit_src.split(":")[0]):
        print(f"  [note] 几何来自 {material!r}, 但 demo 用的锚点 J 来自 "
              f"{lit_src.split(':')[0]!r} —— 系数只吃几何, J 只吃材料; 这是无 DFT 自检,"
              " 产物里的 material 以 --material 为准")

    # ---- 收回凭证: 真值路径必须能说清「这个 J 是哪批 OUTCAR 算出来的」 ----
    prov, prov_problem = None, None
    if not args.demo:
        prov, prov_problem = read_ingest_manifest(args.cases)
        if prov is None:
            msg = (f"[provenance] {prov_problem}\n"
                   "  J 的下游 (mc_2d/mc_fss) 与任何物理结论都需要这份凭证: 没有它就分不清"
                   "「真实 DFT 能量」和「看起来正常的坏能量」。\n"
                   "  正确做法: python run_2d_magnet.py ingest --outdir <outdir> "
                   "(或 scripts/ingest_dft.py <runs> <cases>)。\n"
                   "  确要用未经闸门的 OUTCAR (只在万不得已、且已人工逐条核对时): "
                   "--no-provenance-check")
            if args.no_provenance_check:
                for line in msg.splitlines():
                    print("  [warn] " + line)
                print("  [warn] 已放行: J_result.json 将**没有**收回凭证")
            else:
                raise SystemExit(msg)
        else:
            print(f"收回凭证: {prov['manifest']}")
            print(f"  来源 runs: {prov['source_runs']}   口径: {prov['energy_convention']}")

    if not args.no_jresult:
        d_ang = meta.get("shells_angstrom") or {}
        shells = [dict(name=nm, d=float(d_ang[nm]))
                  for nm in ("J1", "J2", "J3") if nm in d_ang]
        if len(shells) != 3:
            raise SystemExit(f"magnetic_orders.json 缺 shells_angstrom 的 J1/J2/J3: {d_ang}")
        notes = [f"coefficients_per_cell source = {os.path.basename(args.cases)}",
                 "kagome (3 磁位/胞): J1 = a/2 三角形内 (J1/J1' 同距离, 本模型取其平均); "
                 "J2 = a√3/2 跨六角; J3 = a 同亚晶格对角",
                 f"可解性: rank(A)={rk}/{len(ORDERS)}, cond(A)={cond:.4g} "
                 "(q=0 一族对 J1/J2 简并, 必须含 √3×√3 序)"]
        if args.synthetic:
            notes.append("demo 能量由**合成** J 生成 (非文献值): "
                         f"J1={lit['J1']}, J2={lit['J2']}, J3={lit['J3']} meV")
        elif args.demo:
            notes.append("demo 能量由锚点 J 生成: " + lit_src)
        if prov:
            notes.append(f"收回凭证: {prov['manifest']} "
                         f"(source_runs={prov['source_runs']}, 生成于 {prov['generated_utc']})")
            if prov.get("lenient"):
                notes.append(f"⚠ 以下闸门被降级为警告: {prov['lenient']} —— "
                             "J 的可信度随之下降, 引用前先看 ingest_report.md")
        elif not args.demo:
            notes.append("⚠ **没有收回凭证** (--no-provenance-check): 这批 OUTCAR 未过 "
                         "ingest_dft.py 的收敛/口径/结构/防旧文件闸门")
        obj = JR.build(
            lattice="kagome",
            material=material,
            shells=shells, J=dict(J1=J1, J2=J2, J3=J3), E=E, E0=e0,
            cond=cond, residual=residual, spin=meta["spin"], n_mag=n,
            orders=ORDERS, A=A, supercell=meta.get("supercell"),
            source="demo" if args.demo else "vasp", z_from="fm",
            notes=notes,
            dft_provenance=prov,
        )
        out = args.jfile_out or os.path.join(args.cases, "J_result.json")
        warns = JR.write(out, obj)
        print(f"\n已写出统一产物: {os.path.abspath(out)}  (confidence={obj['confidence']})")
        for m in warns:
            print(f"  [warn] {m}")

    if args.demo:
        # 往返误差: E 由锚点 J 经同一个 A 合成, 反解必须逐位还原。
        err = max(abs({"J1": J1, "J2": J2, "J3": J3}[k] - lit[k]) for k in ("J1", "J2", "J3"))
        ok = err < 1e-6 and ok_z and cond < 1e4 and rk == len(ORDERS)
        print(f"\n往返误差 (J1,J2,J3 vs 锚点) = {err:.3g} meV  ->  " + ("PASS" if ok else "FAIL"))
        if not ok:
            raise SystemExit("--demo 自检 FAIL")


if __name__ == "__main__":
    main()
