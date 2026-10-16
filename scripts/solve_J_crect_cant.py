#!/usr/bin/env python
"""solve_J_crect_cant.py — C 心矩形（CrOCl 类）**弱铁磁倾斜（±θ）**构型的 J 反解器（P2′，系数层）。

一句话定位: 本脚本**不反解 DMI**, 也**不能**反解 DMI。它做的是三件事:

  ① **θ 感知的设计矩阵**: 系数不写死共线公式, 而是直接吃 `magnetic_orders.json` 里
     builder 按倾斜构型算出的 `coefficients_per_cell`（含倾斜序）⇒ 倾斜构型的能量也能
     进同一个 `E = e0 + Σ_k c_k·J_k~`（`J~ = J·S²`）线性系统。
  ② **独立复核 θ 标度**: 倾斜系数必须满足闭式关系
         c_k(θ) = sin²θ·c_k(fm) + cos²θ·c_k(B)        (B = 母序, 本脚本**不信任名字**, 用该式反查)
     这一条是**独立实现**: 它只用 `coefficients_per_cell` 的 fm 行 + 母序行 + 声明的 θ,
     不调用 builder 的任何函数。对不上就硬报错（系数与声明的 θ 不一致 ⇒ 构型文件被人手改过,
     或者 builder 的符号约定变了）。
  ③ **把「为什么反解不出 D」写成输出**, 而不是让残留悄悄过去（见下面的边界）。

物理与边界（读这段再改代码）
-------------------------------------------------------------------------------
* 倾斜对的键上: `S_i·S_j = sin²θ + σ_iσ_j·cos²θ` ⇒ **面内/面外同向**的子晶格键（J2/J3）
  与 θ 无关, **两套亚晶格之间**的 C 心键（J1）按 `−cos2θ` 变。
* 因此 **`cant_p8` 与 `cant_m8` 的系数逐位相同**（cos2θ 是偶函数）⇒ **θ 的符号不可由总能量辨识**,
  脚本会显式点名这类退化对; 面内方位角 φ 更是**几何/总能都定不出来**（只有 SAXIS + 磁各向异性
  能定, 那要真实 SOC 计算）。
* **倾斜构型不给设计矩阵增加秩**: 由 ② 的闭式可见, θ 已知时倾斜行与**它自己的母序行等价**
  （`c_B = (c(θ) − sin²θ·c_fm) / cos²θ`）⇒ 它能提供的秩 = rank{fm 行} ∪ {各母序行},
  **不引入任何新方向**（同一个母序的倾斜家族自身只有 rank 2）。所以倾斜构型**可以用来**
  反解同一组 J1/J2/J3（本脚本的 `--demo` 就是证明这件事）, 但**不能**分离 DMI,
  也不能分离 J1 与 J1′ ——
  真要有 D 的信息, 需要 (a) 一个真正独立的通道（换 φ / 换 DM 方向）, (b) DM 矢量方向的对称性
  信息, 二者都要**真实 SOC DFT**。本脚本对这种「残留不为零」只报告、**绝不猜 D**。
* 本脚本**没有** DFT: `--demo` 用锚点 J 合成能量自检（`source=demo`）, 真实路径读
  `cases/<order>/OUTCAR` 的 `free energy TOTEN`（口径统一走 `solve_J_crect.read_toten`）。

用法:
  python scripts/solve_J_crect_cant.py <cases>                       # 默认用 selected_orders（4 共线序）
  python scripts/solve_J_crect_cant.py <cases> --orders cant          # 4 个共线序 + 倾斜家族（额外自洽校验）
  python scripts/solve_J_crect_cant.py <cases> --orders all           # 最小二乘, 用上所有序
  python scripts/solve_J_crect_cant.py <cases> --audit-only           # 只看系数审计表, 不要 OUTCAR
  python scripts/solve_J_crect_cant.py <cases> --demo --synthetic 5 -0.5 0.2

⚠ CrOCl 的文献 J **未核实**（`literature_J.json` 里 J1/J2/J3 = null）⇒ `--ref CrOCl` 会明确拒绝,
  且 `--demo` 必须显式给 `--ref` 或 `--synthetic` 之一 —— 绝不用 0 或别的材料的值顶上。

约定: `H = -J1 Σ S·S - J2 Σ S·S - J3 Σ S·S`, J>0 = 铁磁, meV; 系数约定与
`build_crect_cases.py` 一致（`c_k = -Σ_bonds S_i·S_j`, 每胞）。
"""
import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import j_result as JR  # noqa: E402
import mag_vectors as MV  # noqa: E402
import solve_J_crect as C  # noqa: E402

# 本晶格家族的**定量**锚点 = 没有 (CrOCl 的文献 J 未核实) ⇒ 故意留 None (同 solve_J_crect.py)。
LIT = None
LIT_SOURCE = C.LIT_SOURCE

# 系数标度身份的容差。系数是 O(1)~O(10) 的量, 闭式关系是**精确**的代数恒等式,
# 所以容差可以给得很紧; 放宽它只会掩盖「系数与声明的 θ 不一致」这类真错。
AUDIT_TOL = 1e-9

RANK_GATE = 4                 # 未知量 = {e0, J1~, J2~, J3~}
COND_GATE = 1e4               # 与 kagome 的闸门同量级 (JR.COND_WARN)

SHELL_NAMES = ("J1", "J2", "J3")


# --------------------------------------------------------------------------- 元数据
def load_meta(cases):
    """读 `cases/magnetic_orders.json`（builder 产物）。缺文件/缺系数一律硬报错。"""
    p = os.path.join(cases, "magnetic_orders.json")
    if not os.path.exists(p):
        raise SystemExit(
            f"找不到 {p}\n"
            "  先跑 builder:  python scripts/build_crect_cases.py <POSCAR> --element Cr "
            "--nonmag Cl O --outdir <cases> [--cant-theta \"8 16\"]")
    with open(p, encoding="utf-8-sig") as f:
        meta = json.load(f)
    if not meta.get("coefficients_per_cell"):
        raise SystemExit(f"{p} 里没有 coefficients_per_cell —— 不是 build_crect_cases.py 的产物?")
    return p, meta


def order_names(meta):
    """全部磁序名, 顺序 = all_orders（保 builder 的稳定顺序）, 再补 cant_orders。"""
    names = [k for k in (meta.get("all_orders") or []) if k in meta["coefficients_per_cell"]]
    for k in meta["coefficients_per_cell"]:
        if k not in names:
            names.append(k)
    return names


def cant_orders(meta):
    """倾斜序名。优先用显式名单; 没有就按前缀兜底（并保持稳定顺序）。"""
    listed = [k for k in (meta.get("cant_orders") or []) if k in meta["coefficients_per_cell"]]
    if listed:
        return listed
    return [k for k in order_names(meta) if k.startswith("cant")]


def collinear_orders(meta):
    cant = set(cant_orders(meta))
    return [k for k in order_names(meta) if k not in cant]


def theta_of(meta, name):
    """倾斜构型的声明 θ（度, 带符号）。两级来源: 顶层 cant_theta_deg → cant_report。"""
    t = (meta.get("cant_theta_deg") or {}).get(name)
    if t is None:
        ent = (meta.get("cant_report") or {}).get(name) or {}
        t = ent.get("theta_deg")
    if t is None:
        raise SystemExit(
            f"倾斜序 {name!r} 在 magnetic_orders.json 里没有声明的 θ "
            f"(既不在 cant_theta_deg 也不在 cant_report) ⇒ 无法复核 θ 标度, 拒绝继续")
    return float(t)


def declared_base(meta, name):
    """builder 声明的母序（只用于报错信息; 判定一律走闭式反查）。"""
    ent = (meta.get("cant_report") or {}).get(name) or {}
    if ent.get("base"):
        return ent["base"]
    for spec in (meta.get("order_spec") or []):
        if spec.get("name") == name and spec.get("base"):
            return spec["base"]
    return None


# --------------------------------------------------------------------------- 审计
def canted_coeff(coeff_fm, coeff_base, theta_deg):
    """闭式 θ 标度: c_k(θ) = sin²θ·c_k(fm) + cos²θ·c_k(B)。

    builder 的系数约定: c_k = -Σ_bonds S_i·S_j（每胞）。倾斜 ±θ 后键上
    S_i·S_j = sin²θ + σ_iσ_j·cos²θ, 而 Σ_bonds 1 = -c_k(fm) 且 Σ_bonds σ_iσ_j = -c_k(B)
    ⇒ c_k(θ) = sin²θ·c_k(fm) + cos²θ·c_k(B)。θ 与 -θ 给出同一个式子（cos2θ 偶）。
    """
    s2 = float(np.sin(np.radians(theta_deg)) ** 2)
    c2 = float(np.cos(np.radians(theta_deg)) ** 2)
    return [s2 * float(a) + c2 * float(b) for a, b in zip(coeff_fm, coeff_base)]


def identify_base(meta, name, tol=AUDIT_TOL):
    """用闭式 θ 标度**反查**母序（不信任名字/元数据里的 base 声明）。

    返回 (base, max_abs_dev, predicted)。若没有任何共线序能在 tol 内解释这条倾斜行, 就
    SystemExit —— 这通常意味着构型目录被人手改过, 或系数约定变了（绝不当成「差不多」放过去）。
    """
    th = theta_of(meta, name)
    coeff = meta["coefficients_per_cell"]
    if "fm" not in coeff:
        raise SystemExit("magnetic_orders.json 缺 fm 行 ⇒ 无法复核 θ 标度（fm 行是闭式的参照）")
    c_name = [float(x) for x in coeff[name]]
    hits, best = [], (None, float("inf"), None)
    for b in collinear_orders(meta):
        if b == name:
            continue
        pred = canted_coeff(coeff["fm"], coeff[b], th)
        dev = max(abs(p - q) for p, q in zip(pred, c_name))
        if dev < best[1]:
            best = (b, dev, pred)
        if dev <= tol:
            hits.append((b, dev, pred))
    if not hits:
        dec = declared_base(meta, name)
        raise SystemExit(
            f"倾斜序 {name!r} (θ={th:g}°) 的系数不满足闭式标度 c_k(θ) = sin²θ·c_k(fm) + cos²θ·c_k(B):\n"
            f"  最接近的共线序 {best[0]!r} 的最大偏差 = {best[1]:.3g} (容差 {tol:g})\n"
            f"  builder 声明的母序 = {dec!r}\n"
            "  可能原因: 构型目录被人手改过 / 系数不是按倾斜自旋算的 / 系数约定变了。"
            "不要放宽 --tol 硬过 —— 系数错了 J 就是错的。")
    dec = declared_base(meta, name)
    if dec:
        for b, dev, pred in hits:
            if b == dec:
                return b, dev, pred
    hits.sort(key=lambda h: h[1])
    return hits[0]


def coefficient_audit(meta, tol=AUDIT_TOL):
    """逐条复核倾斜系数: 母序、闭式偏差、退化对（θ 与 −θ 不可分辨）。"""
    audit = {"tol": tol, "orders": {}, "degenerate_pairs": [], "noncollinear": bool(meta.get("noncollinear"))}
    for nm in cant_orders(meta):
        base, dev, pred = identify_base(meta, nm, tol=tol)
        th = theta_of(meta, nm)
        audit["orders"][nm] = {
            "theta_deg": th,
            "base": base,
            "max_dev": float(dev),
            "predicted": [float(x) for x in pred],
            "coeff": [float(x) for x in meta["coefficients_per_cell"][nm]],
            # 与母序的「差异」本身就是信息: 差异为零 ⇒ 这个 θ 白跑了（构型无信息量）
            "delta_vs_base": [float(float(a) - float(b)) for a, b in
                              zip(meta["coefficients_per_cell"][nm], meta["coefficients_per_cell"][base])],
        }

    # 系数逐位相同的序 ⇒ 总能量不可区分。对倾斜家族这就是 (θ, −θ) 对（cos2θ 偶函数）。
    groups = {}
    for nm in order_names(meta):
        key = tuple(round(float(x), 9) for x in meta["coefficients_per_cell"][nm])
        groups.setdefault(key, []).append(nm)
    for key, names in groups.items():
        if len(names) > 1 and any(n in audit["orders"] for n in names):
            audit["degenerate_pairs"].append(sorted(names))
    audit["degenerate_pairs"].sort()
    return audit


def print_audit(meta, audit):
    print(f"\n[审计] 倾斜系数 vs 闭式标度 c_k(θ) = sin²θ·c_k(fm) + cos²θ·c_k(母序)  (容差 {audit['tol']:g})")
    if not audit["orders"]:
        print("  没有倾斜序（magnetic_orders.json 不含 cant_orders）—— 只有共线序, 走老路径即可")
    for nm, ent in audit["orders"].items():
        d = ent["delta_vs_base"]
        print(f"  {nm:<12} θ={ent['theta_deg']:+8.4g}°  母序={ent['base']:<10} "
              f"最大偏差={ent['max_dev']:.3g}  与母序系数差={['%+.4g' % x for x in d]}")
    for pair in audit["degenerate_pairs"]:
        print(f"  [退化] {' / '.join(pair)} 的系数**逐位相同** ⇒ 总能量完全不可区分"
              f"（cos2θ 是偶函数 ⇒ θ 的符号不可辨识）")
    if audit["orders"]:
        print("  [边界] 每一条倾斜行都是「fm 行」与「母序行」的仿射组合 ⇒ θ 已知时它**等价于母序行**"
              "（不引入新方向; 同一个母序的倾斜家族自身 rank = 2）⇒ **不给设计矩阵增加秩**: "
              "能反解同一组 J, 但**不能**分离 DMI, 也不能分离 J1/J1′")


# --------------------------------------------------------------------------- 设计矩阵
def choose_orders(meta, spec):
    """把 --orders 解析成 (orders, note)。规则写死在这里, 便于测试与报错引用。"""
    sel = [k for k in (meta.get("selected_orders") or []) if k in meta["coefficients_per_cell"]]
    cant = cant_orders(meta)
    if spec in (None, "", "default", "collinear"):
        if not sel:
            raise SystemExit("magnetic_orders.json 没有 selected_orders ⇒ 用 --orders 显式点名")
        return sel, "selected_orders（4 个共线序, 与 solve_J_crect.py 同一条路径）"
    if spec == "all":
        return order_names(meta), "全部磁序（最小二乘）"
    if spec == "cant":
        # 倾斜家族**必须和共线序一起用**: 倾斜行只是「fm 行 + 母序行」的仿射组合,
        # 单独一堆倾斜行 rank 不足（见审计的边界一行）。所以这里默认把 selected_orders
        # 一起带上, 让倾斜构型变成**同一组 J 的额外自洽校验**, 而不是假装它们带来新信息。
        picks = list(sel)
        for nm in cant:
            if nm not in picks:
                picks.append(nm)
        if not picks:
            raise SystemExit(
                "magnetic_orders.json 既没有 selected_orders 也没有倾斜序 ⇒ 用 --orders 显式点名")
        return picks, "selected_orders + 倾斜家族（倾斜行不给秩, 是同一组 J 的额外自洽校验）"
    picks = [s.strip() for s in str(spec).split(",") if s.strip()]
    if not picks:
        raise SystemExit("--orders 给了空列表")
    unknown = [p for p in picks if p not in meta["coefficients_per_cell"]]
    if unknown:
        raise SystemExit(f"--orders 里有未知磁序 {unknown}; 现有 {order_names(meta)}")
    return picks, "显式点名"


def design_matrix(meta, orders):
    """A = [[1, c1, c2, c3], ...]; 未知量 = [e0, J1~, J2~, J3~]。"""
    rows = []
    for nm in orders:
        c = [float(x) for x in meta["coefficients_per_cell"][nm]]
        if len(c) != len(SHELL_NAMES):
            raise SystemExit(f"磁序 {nm!r} 的系数长度 {len(c)} != {len(SHELL_NAMES)}（壳数对不上）")
        rows.append([1.0] + c)
    return np.array(rows)


def fit(A, b, orders):
    """解 / 最小二乘 + rank 与条件数闸门。返回 (x, cond, residual, method)。"""
    rank = int(np.linalg.matrix_rank(A))
    cond = float(np.linalg.cond(A))
    if rank < RANK_GATE:
        raise SystemExit(
            f"设计矩阵 rank={rank} < {RANK_GATE} ⇒ 这组磁序**无法唯一确定** e0/J1/J2/J3: {orders}\n"
            "  倾斜构型只是「fm 行 + 母序行」的仿射组合 ⇒ 光靠倾斜序凑不满秩; "
            "至少要凑够 4 条线性无关的行（例: --orders all, 或 4 个共线序 + 倾斜序一起用）。")
    if len(orders) == RANK_GATE:
        x = np.linalg.solve(A, b)
        method = "solve(4x4)"
    else:
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        method = f"lstsq({A.shape[0]}x{A.shape[1]})"
    residual = float(np.max(np.abs(A @ x - b)))
    return x, cond, residual, method


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", help="build_crect_cases.py 的 --outdir（内含 magnetic_orders.json）")
    ap.add_argument("--orders", default=None,
                    help="default|collinear|all|cant 或逗号分隔的磁序名 "
                         "(默认 default = selected_orders; cant = 倾斜家族 + 张成所需的 fm/母序行)")
    ap.add_argument("--audit-only", action="store_true",
                    help="只打印系数审计表（不需要 OUTCAR, 不下结论）")
    ap.add_argument("--demo", action="store_true",
                    help="用锚点 J 合成磁态能量自检（不需要 DFT）; 需配 --ref 或 --synthetic")
    ap.add_argument("--ref", default=None, help="从 literature_J.json 取锚点的材料名（或该 json 的路径）")
    ap.add_argument("--synthetic", nargs=3, type=float, default=None, metavar=("J1", "J2", "J3"),
                    help="显式给一组**合成** J (meV) 做无 DFT 自检; 不是文献值, 产物标 demo")
    ap.add_argument("--tol", type=float, default=AUDIT_TOL,
                    help=f"闭式 θ 标度的判定容差（默认 {AUDIT_TOL:g}; 放宽容差只会掩盖真错）")
    ap.add_argument("--material", default=None, help="材料名 (写进 J_result.json)")
    ap.add_argument("--jfile-out", default=None, help="J_result.json 输出路径 (默认 <cases>/J_result.json)")
    ap.add_argument("--no-jresult", action="store_true", help="不写统一产物 (只打印)")
    args = ap.parse_args()

    _p, meta = load_meta(args.cases)
    n = meta["n_mag"]
    S2 = meta["S2"]
    coeff = meta["coefficients_per_cell"]

    audit = coefficient_audit(meta, tol=args.tol)
    print_audit(meta, audit)

    if args.audit_only:
        # 审计模式不拟合, --orders 不参与计算; 但**仍然校验**: 打错一个序名不该被静默忽略。
        if args.orders:
            choose_orders(meta, args.orders)
        print("\n[audit-only] 只做系数审计, 没读能量、没下任何 J 的结论")
        return 0

    orders, why = choose_orders(meta, args.orders)
    print(f"\n用 {len(orders)} 个磁序反解 ({why}): {orders}")
    _need_soc = [o for o in orders if MV.needs_soc(o)]
    if _need_soc:
        print(f"  [硬要求] {_need_soc} 是 SOC 构型 ⇒ 上游 assemble_dft 必须 "
              f"`--noncollinear` **且** `--lsorbit`（缺了会硬报错, 不会静默降级）")

    # ---- 锚点解析: --ref 与 --synthetic 互斥, --demo 必须给其中之一 ------------
    if args.ref and args.synthetic:
        raise SystemExit("--ref 与 --synthetic 互斥: 一个是文献锚点, 一个是合成值, 不许混用")
    lit, lit_src = LIT, LIT_SOURCE
    if args.ref:
        lit, lit_src = C.load_ref(args.ref)
    elif args.synthetic:
        lit = dict(J1=args.synthetic[0], J2=args.synthetic[1], J3=args.synthetic[2])
        lit_src = "synthetic (显式给的合成 J, **不是文献值**)"
    if args.demo and lit is None:
        raise SystemExit(
            "本晶格的文献 J 未核实 ⇒ --demo 需要一个显式锚点, 二者之一:\n"
            "  --ref <材料名>           从 literature_J.json 取已核实的文献 J\n"
            "  --synthetic J1 J2 J3     显式给合成 J 做无 DFT 自检 (例: --synthetic 5 -0.5 0.2)\n"
            "拒绝用 0 或别的材料的值顶上 (WORKFLOW §99 不静默降级)。")

    A = design_matrix(meta, orders)

    if args.demo:
        Jt = np.array([lit["J1"], lit["J2"], lit["J3"]]) * S2      # J~ = J*S^2 (meV)
        E = {k: (float(np.dot(coeff[k], Jt)) / 1000.0) for k in orders}
        print("[demo] 锚点 J 来源: " + lit_src)
        print("[demo] 用该 J 生成的磁态能量 (eV/胞):",
              {k: f"{v:.6f}" for k, v in E.items()})
        # 倾斜构型的能量与母序不同吗? 相同就说明这批构型**没有信息量**（例如拿 fm 当母序）。
        bases = {audit["orders"][k]["base"] for k in orders if k in audit["orders"]}
        for b in sorted(bases):
            if b in orders:
                same = [k for k in orders if k in audit["orders"]
                        and audit["orders"][k]["base"] == b and abs(E[k] - E[b]) < 1e-12]
                if same:
                    print(f"  [warn] {same} 的能量与母序 {b} **逐位相同** ⇒ 这些构型对该 J 无信息量"
                          f"（典型原因: 母序取了 fm —— 整体转动不改能量）")
    else:
        E = {}
        for nm in orders:
            p = os.path.join(args.cases, nm, "OUTCAR")
            if os.path.exists(p):
                E[nm] = C.read_toten(p)
        if len(E) != len(orders):
            miss = [nm for nm in orders if nm not in E]
            raise SystemExit(
                f"缺 {miss} 的 OUTCAR（需要 {orders}）; 可加 --demo 做无 DFT 自检, "
                f"或 --audit-only 只看系数审计")

    b = np.array([E[k] for k in orders])
    x, cond, residual, method = fit(A, b, orders)
    e0, J1t, J2t, J3t = x

    J1, J2, J3 = (J1t / S2 * 1000.0, J2t / S2 * 1000.0, J3t / S2 * 1000.0)

    print(f"\n每胞系数: " + "  ".join(f"{k}{['%+.4g' % v for v in coeff[k]]}" for k in orders))
    print(f"能量 (eV/胞): " + "  ".join(f"{k}={E[k]:+.6f}" for k in orders))
    print(f"解: {method}   cond(A) = {cond:.4g}")
    print(f"J1~ = {J1t:+.6f} eV   J2~ = {J2t:+.6f} eV   J3~ = {J3t:+.6f} eV   (S^2 = {S2})")
    print(f"\nJ1 = {J1:+.4f} meV    J2 = {J2:+.4f} meV    J3 = {J3:+.4f} meV")
    print("   (J1 = 3.0347 A C 心矢量 x4 [连胞内两个 Cr], "
          "J2 = 3.1980 A 沿 a x2, J3 = 3.8603 A 沿 b x2)")

    # 配位数自检: 由 fm 行直接给出, 纯几何量 ⇒ 必须是 4/2/2（与 solve_J_crect.py 同一条闸门）。
    if "fm" in coeff:
        z = tuple(round(2.0 * (-float(coeff["fm"][i])) / n, 6) for i in range(len(SHELL_NAMES)))
        ok_z = all(abs(z[i] - C.EXPECT_Z_N[i]) < 1e-9 for i in range(len(SHELL_NAMES)))
        print(f"配位数 (由 fm 行导出, 每磁位点): J1={z[0]:g} J2={z[1]:g} J3={z[2]:g}"
              f"   期望 {C.EXPECT_Z_N}  ->  " + ("PASS" if ok_z else "FAIL"))
        if not ok_z:
            raise SystemExit(
                "配位数与 C 心矩形几何不符 (期望 4/2/2) —— 壳很可能被 magnetic_shells(tol=0.2) "
                "静默合并了 (J1 与 J2 只差 0.1632 A)。请用显式壳距重建 magnetic_orders.json。")
    else:
        z, ok_z = None, None

    print(f"拟合质量: residual = {residual:.3g} eV   cond = {cond:.4g}")
    if cond > COND_GATE:
        print(f"  [warn] cond(A) = {cond:.4g} > {COND_GATE:.4g} —— 设计矩阵接近病态,"
              " 反解结果对能量噪声敏感")
    if residual > 1e-6:
        print(f"  [warn] residual = {residual:.3g} eV **不为零** ⇒ 这个 Heisenberg 三壳模型"
              "解释不了这批能量: 可能是 DMI / 磁各向异性 / 壳合并 / J1 与 J1′ 没分开。"
              "本脚本**不反解 D**（倾斜构型不给设计矩阵增加秩, 且 DM 方向要对称性信息 + 真实 SOC DFT）"
              "—— 残留只当「模型缺项」的指纹报告, 不许当成 D 的值。")
    if args.demo and lit is not None:
        dev = max(abs(J1 - lit["J1"]), abs(J2 - lit["J2"]), abs(J3 - lit["J3"]))
        print(f"[demo] 往返误差 max|ΔJ| = {dev:.3g} meV  ->  " + ("PASS" if dev < 1e-6 else "FAIL"))
        if dev >= 1e-6:
            return 1

    material = args.material or C.guess_material(args.cases)
    if args.demo and lit is not None and material and \
            C._lit_key(material) != C._lit_key(lit_src.split(":")[0]):
        print(f"  [note] 几何来自 {material!r}, 但 demo 用的锚点 J 来自 "
              f"{lit_src.split(':')[0]!r} —— 系数只吃几何, J 只吃材料; 这是无 DFT 自检,"
              " 产物里的 material 以 --material 为准")

    if args.no_jresult:
        return 0

    shells = []
    for nm in SHELL_NAMES:
        d = (meta.get("shells_angstrom") or {}).get(nm)
        shells.append(dict(name=nm, d=(float(d) if d is not None else None),
                           tol=meta.get("shell_tol")))

    notes = [
        "P2′ 系数层: 倾斜构型的系数由 builder 按声明的 θ 算出; 本脚本用闭式 "
        "c_k(θ)=sin²θ·c_k(fm)+cos²θ·c_k(B) **独立复核**（不调用 builder 的任何函数）",
        "倾斜构型是 fm 行与母序行的仿射组合 ⇒ **不给设计矩阵增加秩**: 能反解同一组 J, "
        "但**不能**分离 DMI, 也不能分离 J1/J1′",
        "θ 的符号不可由总能量辨识（cos2θ 偶函数 ⇒ cant_p/cant_m 系数逐位相同）; "
        "面内方位角 φ 由几何与总能都定不出来（需 SAXIS + 磁各向异性 ⇒ 真实 SOC DFT）",
        "与共线路径共享同一条能量口径（free energy TOTEN, dft_energy.read_toten）",
    ]
    if not meta.get("noncollinear"):
        notes.append("magnetic_orders.json 标 noncollinear=False ⇒ 这批构型是共线的（没有走 SOC 通道）")

    obj = JR.build(
        lattice="rect_centered", shells=shells,
        J=dict(J1=J1, J2=J2, J3=J3), E=E, E0=float(e0),
        cond=cond, residual=residual, spin=meta["spin"], S2=S2, n_mag=n,
        orders=orders, A=A, supercell=meta.get("supercell"), material=material,
        source=("demo" if args.demo else "vasp"), z_from="fm", notes=notes,
        cant_audit=audit, orders_spec=(args.orders or "default"),
    )
    out = args.jfile_out or os.path.join(args.cases, "J_result.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=2)
    print(f"\n产物: {out}  (source={obj['source']}, confidence={obj['confidence']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
