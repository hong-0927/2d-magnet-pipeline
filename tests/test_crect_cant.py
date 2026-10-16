#!/usr/bin/env python
"""test_crect_cant.py — P2′（C 心矩形 弱铁磁倾斜 ±θ）反解器的契约测试

**不依赖 builder、不依赖 DFT**: 本文件自己合成一份 `magnetic_orders.json`
（4 条共线序 + 5 条倾斜序), 系数用**测试内独立实现**的闭式生成, 并用一个**硬编码锚点**钉住公式。

覆盖:
  1. **闭式 θ 标度**: `c_k(θ) = sin²θ·c_k(fm) + cos²θ·c_k(B)` 与测试内独立实现逐位一致;
     θ→0 回到母序、θ→90° 回到 fm; θ 与 −θ 系数逐位相同（cos2θ 偶函数）;
     另有一个手算硬编码锚点 (8°/afm_c 的 J1 = 3.84504880…), 防止「公式抄错还自洽」。
  2. **名字无关的母序识别**: 只吃系数 + 声明的 θ 反查母序; 系数被人手改 1e-3 必须硬报错
     （绝不放宽容差蒙过去）。
  3. **退化与边界**: (θ, −θ) 对与「母序=fm」的行被点名; 倾斜行只是 fm 与母序行的仿射组合
     ⇒ 单独一堆倾斜行 rank 不足必须被闸门拦住（这正是「反解不出 DMI」的算术原因）。
  4. **设计矩阵与反解**: 4 条共线行的 rank=4、J 恢复到 1e-9; 2 条倾斜行 rank=2 ⇒ 报错。
  5. **CLI 契约**: --audit-only 不需要 OUTCAR; --demo 必须有锚点（--ref 或 --synthetic）;
     --ref CrOCl（文献 J 为 null）必须拒绝; --ref 与 --synthetic 互斥; 缺 OUTCAR 报错;
     配位数闸门（fm 行导出 4/2/2）; 产物字段（source=demo / confidence=low / cant_audit / shells.tol）。

运行:
  python tests/test_crect_cant.py
  pytest tests/test_crect_cant.py -q
"""
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(HERE, "scripts")
PY = sys.executable
SOLVER = os.path.join(SCRIPTS, "solve_J_crect_cant.py")
CHILD_ENV = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONDONTWRITEBYTECODE="1")

sys.path.insert(0, SCRIPTS)

# 报告/报错文本含 ° / ² / ⇒ 等非 GBK 字符: 只放宽 errors, 不改编码。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors="replace")
    except Exception:
        pass

import solve_J_crect_cant as X  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + (("   " + extra) if extra else ""))
    if not cond:
        fails.append(name)


def run_cli(*cmd, cwd=None):
    return subprocess.run([PY] + list(cmd), capture_output=True, text=True,
                          encoding="utf-8", errors="replace", env=CHILD_ENV,
                          cwd=cwd or HERE)


# ------------------------------------------------------------------ 合成元数据
# 合成模型的bond 计数（每胞）: J1 是 A-B 键 x4, J2 是 A-A 键 x2, J3 是 B-B 键 x2
# （与真 CrOCl 的 z_n = 4/2/2 同形, 这样配位数闸门也走真路径）。
COLLINEAR = {
    "fm":       (-4.0, -2.0, -2.0),
    "afm_c":    (+4.0, -2.0, -2.0),
    "afm_a":    (+0.0, +2.0, -2.0),
    "afm_diag": (-4.0, +2.0, +2.0),
}
CANT = [("cant_p8", "afm_c", +8.0), ("cant_m8", "afm_c", -8.0),
        ("cant_p16", "afm_a", +16.0), ("cant_m16", "afm_a", -16.0),
        ("cant_fm16", "fm", +16.0)]
SPIN, N_MAG = 1.5, 2
TRUE_J = {"J1": 5.0, "J2": -0.5, "J3": 0.2}


def cant_row(c_fm, c_base, theta_deg):
    """测试内**独立实现**的闭式（刻意不复用 solve_J_crect_cant.canted_coeff）:
    每根键上 S_i·S_j = sin²θ + σ_iσ_j·cos²θ, 而 c_k = -Σ_bonds S_i·S_j。"""
    s2 = math.sin(math.radians(theta_deg)) ** 2
    c2 = math.cos(math.radians(theta_deg)) ** 2
    return tuple(s2 * a + c2 * b for a, b in zip(c_fm, c_base))


def make_meta(path, collinear=None):
    coll = dict(collinear or COLLINEAR)
    coeff = {k: list(v) for k, v in coll.items()}
    theta = {}
    report = {}
    for name, base, th in CANT:
        coeff[name] = list(cant_row(coll["fm"], coll[base], th))
        theta[name] = th
        report[name] = {"base": base, "theta_deg": th}
    meta = {
        "n_mag": N_MAG, "spin": SPIN, "S2": SPIN ** 2, "element": "Cr",
        "noncollinear": True, "supercell": [4, 4, 1],
        "shell_tol": 0.05,
        "shells_angstrom": {"J1": 3.0347, "J2": 3.1980, "J3": 3.8603},
        "selected_orders": ["fm", "afm_c", "afm_a", "afm_diag"],
        "all_orders": list(coll) + [c[0] for c in CANT],
        "coefficients_per_cell": coeff,
        "cant_orders": [c[0] for c in CANT],
        "cant_theta_deg": theta,
        "cant_report": report,
    }
    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, "magnetic_orders.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    return meta


TD = tempfile.mkdtemp(prefix="crect_cant_")
CASES = os.path.join(TD, "cases")
META = make_meta(CASES)

# ================================================================== 1. 闭式 θ 标度
print("== 1. 闭式 θ 标度（与测试内独立实现比）")
worst = 0.0
for name, base, th in CANT + [("cant_p8", "afm_c", 90.0), ("cant_p8", "afm_c", 0.0)]:
    mine = X.canted_coeff(COLLINEAR["fm"], COLLINEAR[base], th)
    ref = cant_row(COLLINEAR["fm"], COLLINEAR[base], th)
    worst = max(worst, max(abs(a - b) for a, b in zip(mine, ref)))
check("canted_coeff 与独立实现逐位一致 (<=1e-15)", worst <= 1e-15, "max|Δ|=%.2g" % worst)

r8 = X.canted_coeff(COLLINEAR["fm"], COLLINEAR["afm_c"], 8.0)
# 手算锚点: c_J1(8°) = 4 - 8·sin²8° = 3.8450467837532756（sin8° = 0.13917310096006544,
# **平方**才是权重: 若有人把 sin²θ 写成 sinθ, 这里会得到 2.8866…, 差了 1 eV 一眼就看出来）
check("硬编码锚点: 8°/afm_c 的 J1 = 3.8450467837532756",
      abs(r8[0] - 3.8450467837532756) < 1e-12, "得到 %.13f" % r8[0])
check("硬编码锚点: 8°/afm_c 的 J2=J3=-2（母序与 fm 在这两壳上同值 ⇒ 与 θ 无关）",
      abs(r8[1] + 2.0) < 1e-15 and abs(r8[2] + 2.0) < 1e-15)
check("θ→0 回到母序行", X.canted_coeff(COLLINEAR["fm"], COLLINEAR["afm_c"], 0.0)
      == list(COLLINEAR["afm_c"]))
check("θ→90° 回到 fm 行", X.canted_coeff(COLLINEAR["fm"], COLLINEAR["afm_c"], 90.0)
      == list(COLLINEAR["fm"]))
p8 = X.canted_coeff(COLLINEAR["fm"], COLLINEAR["afm_c"], 8.0)
m8 = X.canted_coeff(COLLINEAR["fm"], COLLINEAR["afm_c"], -8.0)
check("θ 与 −θ 系数逐位相同（cos2θ 偶 ⇒ 符号不可辨识）", p8 == m8)

# ================================================================== 2. 母序识别
print("\n== 2. 名字无关的母序识别 + 改系数必须报错")
for name, base, _th in CANT:
    got, dev, _pred = X.identify_base(META, name)
    check(f"identify_base({name}) = {base}", got == base and dev < 1e-12,
          f"得到 {got} (dev={dev:.2g})")

bad = json.loads(json.dumps(META))
bad["coefficients_per_cell"]["cant_p8"][2] += 1e-3
try:
    X.identify_base(bad, "cant_p8")
    check("系数被改动 1e-3 -> 必须硬报错", False, "居然通过了")
except SystemExit as e:
    check("系数被改动 1e-3 -> 必须硬报错", "闭式标度" in str(e))

# ================================================================== 3. 退化与边界
print("\n== 3. 退化对 + 秩边界")
audit = X.coefficient_audit(META)
pairs = [list(p) for p in audit["degenerate_pairs"]]
check("点名 (cant_p8, cant_m8) 退化对", ["cant_m8", "cant_p8"] in pairs, str(pairs))
check("点名 (cant_p16, cant_m16) 退化对", ["cant_m16", "cant_p16"] in pairs)
check("点名 (cant_fm16, fm) 退化对（母序=fm ⇒ 该行无信息量）",
      ["cant_fm16", "fm"] in pairs)
check("审计记录母序与偏差", all(audit["orders"][n]["base"] == b and
                          audit["orders"][n]["max_dev"] < 1e-12 for n, b, _ in CANT))

picks, _why = X.choose_orders(META, "cant")
check("--orders cant = selected_orders + 倾斜家族（9 行, 共线序在前）",
      picks[:4] == META["selected_orders"] and len(picks) == 9
      and set(picks[4:]) == {c[0] for c in CANT}, str(picks))
A = X.design_matrix(META, META["selected_orders"])
check("4 条共线行的设计矩阵 rank = 4", int(np.linalg.matrix_rank(A)) == 4)

A2 = X.design_matrix(META, ["cant_p8", "cant_m8"])
try:
    X.fit(A2, np.zeros(2), ["cant_p8", "cant_m8"])
    check("只有 2 条倾斜行 -> rank 闸门必须拦住", False, "居然通过了")
except SystemExit as e:
    check("只有 2 条倾斜行 -> rank 闸门必须拦住", "无法唯一确定" in str(e))

# 「倾斜行不给秩」的**算术证明**（θ 已知 ⇒ 倾斜行 ⇔ 它自己的母序行）:
#   a) 同一个母序的倾斜家族 = rank 2（fm + 母序两根方向, 第 3/4 个未知量无从谈起）;
#   b) 倾斜行只会把它**自己的母序方向**带进来, 不引入任何新方向:
#      rank{fm, afm_c, afm_a, 全部倾斜行} == rank{fm, afm_c, afm_a} == 3。
r_same_base = ["fm", "afm_c", "cant_p8", "cant_m8"]
A_same = X.design_matrix(META, r_same_base)
check("同一母序的倾斜家族 rank = 2（fm + 母序, 没有第三个方向）",
      int(np.linalg.matrix_rank(A_same)) == 2, "rank=%d" % int(np.linalg.matrix_rank(A_same)))
try:
    X.fit(A_same, np.zeros(len(r_same_base)), r_same_base)
    check("这套 rank=2 的构型必须被闸门拦住", False, "居然通过了")
except SystemExit as e:
    check("这套 rank=2 的构型必须被闸门拦住",
          "无法唯一确定" in str(e) and "仿射组合" in str(e))

r_bases = ["fm", "afm_c", "afm_a"]
A_b = X.design_matrix(META, r_bases)
A_all = X.design_matrix(META, r_bases + [c[0] for c in CANT])
check("倾斜行只带来自己的母序方向（加不加倾斜行, rank 一样）",
      int(np.linalg.matrix_rank(A_b)) == int(np.linalg.matrix_rank(A_all)) == 3,
      "rank(共线)=%d rank(+倾斜)=%d" % (int(np.linalg.matrix_rank(A_b)),
                                    int(np.linalg.matrix_rank(A_all))))

# ================================================================== 4. 反解恢复
print("\n== 4. 反解（合成能量 -> J）")
Jt = np.array([TRUE_J["J1"], TRUE_J["J2"], TRUE_J["J3"]]) * SPIN ** 2
orders = META["selected_orders"]
b = np.array([float(np.dot(META["coefficients_per_cell"][k], Jt)) / 1000.0 for k in orders])
x, cond, residual, method = X.fit(A, b, orders)
J = [x[1] / SPIN ** 2 * 1000.0, x[2] / SPIN ** 2 * 1000.0, x[3] / SPIN ** 2 * 1000.0]
devj = max(abs(J[0] - TRUE_J["J1"]), abs(J[1] - TRUE_J["J2"]), abs(J[2] - TRUE_J["J3"]))
check("4x4 反解恢复合成 J (<=1e-9 meV)", devj <= 1e-9,
      "max|ΔJ|=%.2g method=%s cond=%.3g" % (devj, method, cond))
check("残差 ~ 0", residual < 1e-15, "residual=%.2g" % residual)
check("e0 = 0（合成能量没有常数项 ⇒ 常数项不许吸进 J）", abs(x[0]) < 1e-12,
      "e0=%.2g" % x[0])

# ================================================================== 5. CLI 契约
print("\n== 5. CLI 契约")
p = run_cli(SOLVER, CASES, "--audit-only")
check("--audit-only 不需要 OUTCAR 也不需锚点, exit 0", p.returncode == 0)
check("--audit-only 打印退化与秩边界说明",
      "退化" in p.stdout and "不给设计矩阵增加秩" in p.stdout)

out_json = os.path.join(TD, "J_all.json")
p = run_cli(SOLVER, CASES, "--orders", "all", "--demo",
            "--synthetic", "5", "-0.5", "0.2", "--material", "CrOCl",
            "--jfile-out", out_json)
check("--orders all --demo 全流程 exit 0", p.returncode == 0, p.stdout[-300:] if p.returncode else "")
check("demo 打印配位数 PASS (4/2/2 闸门走真路径)", "配位数" in p.stdout and "PASS" in p.stdout)
check("demo 往返误差 PASS", "往返误差" in p.stdout and "PASS" in p.stdout)
check("demo 提示无信息量构型 (母序=fm)", "[warn]" in p.stdout and "无信息量" in p.stdout)
if os.path.exists(out_json):
    obj = json.load(open(out_json, encoding="utf-8"))
    check("产物 source=demo / confidence=low",
          obj["source"] == "demo" and obj["confidence"] == "low")
    check("产物 lattice=rect_centered / material 透传",
          obj["lattice"] == "rect_centered" and obj["material"] == "CrOCl")
    check("产物 J_meV 恢复合成值 (<=1e-6)",
          max(abs(obj["J_meV"][k] - TRUE_J[k]) for k in TRUE_J) <= 1e-6,
          str(obj["J_meV"]))
    check("产物 shells 带 builder 的 tol=0.05", all(s["tol"] == 0.05 for s in obj["shells"]))
    check("产物带 cant_audit（母序识别结论进了档案）",
          obj.get("cant_audit", {}).get("orders", {}).get("cant_p8", {}).get("base") == "afm_c")
    check("产物记录 orders_spec = all", obj.get("orders_spec") == "all")
else:
    check("产物 J_result.json 存在", False, out_json)

p = run_cli(SOLVER, CASES, "--orders", "cant", "--demo", "--synthetic", "5", "-0.5", "0.2",
            "--no-jresult")
check("--orders cant = 4 共线 + 5 倾斜 = 9 行最小二乘",
      p.returncode == 0 and "lstsq(9x4)" in p.stdout,
      (p.stdout + p.stderr)[-400:] if p.returncode else "")

p = run_cli(SOLVER, CASES, "--demo", "--orders", "all")
check("--demo 缺锚点 -> 拒绝（提到文献 J 未核实）",
      p.returncode != 0 and "未核实" in (p.stdout + p.stderr))
p = run_cli(SOLVER, CASES, "--demo", "--ref", "CrOCl", "--orders", "all")
check("--ref CrOCl（文献 J 为 null）-> 拒绝",
      p.returncode != 0 and "null" in (p.stdout + p.stderr))
p = run_cli(SOLVER, CASES, "--demo", "--ref", "CrOCl", "--synthetic", "5", "-0.5", "0.2")
check("--ref 与 --synthetic 互斥", p.returncode != 0 and "互斥" in (p.stdout + p.stderr))
p = run_cli(SOLVER, CASES, "--orders", "all")
check("非 demo 且缺 OUTCAR -> 报错并说清缺哪些",
      p.returncode != 0 and "OUTCAR" in (p.stdout + p.stderr))
p = run_cli(SOLVER, CASES, "--orders", "bogus")
check("未知磁序名 -> 报错", p.returncode != 0 and "未知磁序" in (p.stdout + p.stderr))

# 配位数闸门: 把 fm 行改成 z=3.7 的坏几何（cant 行随之重算, 保证审计先通过）
bad_cases = os.path.join(TD, "bad_cases")
make_meta(bad_cases, collinear=dict(COLLINEAR, fm=(-3.7, -2.0, -2.0)))
p = run_cli(SOLVER, bad_cases, "--orders", "all", "--demo", "--synthetic", "5", "-0.5", "0.2")
check("配位数闸门: fm 行导出 3.7 != 4 -> 报错",
      p.returncode != 0 and "配位数" in (p.stdout + p.stderr),
      (p.stdout + p.stderr)[-200:] if p.returncode == 0 else "")

# ------------------------------------------------------------------ 收尾
shutil.rmtree(TD, ignore_errors=True)
print()
if fails:
    print("FAILED %d 项:" % len(fails))
    for f in fails:
        print("  - " + f)
    sys.exit(1)
print("全部通过")
sys.exit(0)
