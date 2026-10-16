#!/usr/bin/env python
"""test_altermagnet.py — 交变磁体 (P5) 识别层 + 序参量的契约测试

覆盖:
  1. **纯函数**: 点群枚举 (point_ops) / 操作分类 (rot_kind) / 笛卡尔旋转 (cart_rotation) /
     轴角 (rot_axis_angle) / 最近同元素距离 (nearest_same_element);
  2. **判定四分支 + 不适用**: 平移连接 / 旋转连接 / 反演连接 / 无连接 / 非两套亚晶格,
     全部用合成 POSCAR (写进临时目录), 不下发任何 DFT;
  3. **容差语义**: 判据容差是**Å 距离**而不是「四舍六入查元组」—— 坐标抖 3e-6 后仍必须
     认出对称操作 (真 fixture 曾被这条静默误报「没有对称操作」), 而把 --tol 收到 1e-6
     又必须真的判不出来 (否则容差是摆设);
  4. **反演在场即宇称禁止**: 同一 t 下同时存在反演与真旋转时不许改判候选; 共面 2D 胞例外
     (那里「反演」与面内 C2 是同一操作);
  5. **真 fixture 的落盘结论与数字**: MnPS3 = 几何候选 (镜面连接, 残差 ~0.01 mÅ),
     CrOCl = 反演禁止, Fe3Sn = 不适用 (3 套等价类);
  6. **CLI 契约**: --selftest / --json / --element / --spins / 错误码;
  7. **与 MC 侧的接线**: `mc_2d.py --order sublattice:<mask>` 在二亚晶格模型上与老
     `--order staggered` **逐位一致**; 空掩码 / 越界 / 覆盖全部亚晶格必须报错;
     驱动有 `--mc-order` 透传; `mc_fss.py` 复用同一套解析。

运行:
  python tests/test_altermagnet.py
  pytest tests/test_altermagnet.py -q
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(HERE, "scripts")
FIXTURES = os.environ.get("MAGNET_FIXTURES", os.path.join(HERE, "tests", "fixtures"))
PY = sys.executable
CHILD_ENV = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONDONTWRITEBYTECODE="1")

POSCAR_MNPS3 = os.path.join(FIXTURES, "MnPS3_POSCAR")
POSCAR_CROCL = os.path.join(FIXTURES, "CrOCl_POSCAR")
POSCAR_FE3SN = os.path.join(FIXTURES, "Fe3Sn_POSCAR")

sys.path.insert(0, SCRIPTS)

# 报告文本含 ⇒ / ² 等非 GBK 字符 (报错的 extra 也会带): 只放宽 errors, 不改编码。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors="replace")
    except Exception:
        pass

import altermagnet_check as A  # noqa: E402
import mc_2d as M  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + (("   " + extra) if extra else ""))
    if not cond:
        fails.append(name)


def run_cli(*cmd, cwd=None):
    return subprocess.run([PY] + list(cmd), capture_output=True, text=True,
                          encoding="utf-8", errors="replace", env=CHILD_ENV,
                          cwd=cwd or HERE)


# ------------------------------------------------------------------ 合成结构

SQ = [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 12.0]]
TET = [[3.0, 0, 0], [0, 3.0, 0], [0, 0, 6.0]]
TRI = [[3.0, 0, 0], [0.5, 3.0, 0], [0.3, 0.4, 10.0]]
HEX = [[5.4, 0, 0], [-2.7, 4.6765371, 0], [0, 0, 20.0]]

TD = tempfile.mkdtemp(prefix="altermag_test_")
_counter = [0]


def synth(lat, species, counts, coords, element, spins=None, tol=1e-3):
    """写一个合成 POSCAR 再跑判定 (返回结果 dict)。"""
    _counter[0] += 1
    p = os.path.join(TD, "synth%d.vasp" % _counter[0])
    A._write_poscar(p, lat, species, counts, coords)
    lat2, atoms, frac, _ = A.MO.parse_poscar(p)
    return A.check(lat2, atoms, frac, element, spins=spins, tol=tol)


# ============================================================ (A) 纯函数
print("--- (A) 纯函数 ---")

G = A.metric(SQ)
check("metric 对称", np.allclose(G, G.T))
check("metric 正定", np.all(np.linalg.eigvalsh(G) > 0))

ops_sq = A.point_ops(SQ)
check("point_ops(正方) = 16", len(ops_sq) == 16, "得到 %d" % len(ops_sq))
check("point_ops 含 ±I",
      any(np.array_equal(m, np.eye(3, dtype=int)) for m in ops_sq)
      and any(np.array_equal(m, -np.eye(3, dtype=int)) for m in ops_sq))

ops_hex = A.point_ops(HEX)
check("point_ops(120° 六角) = 24", len(ops_hex) == 24, "得到 %d" % len(ops_hex))
check("point_ops(三斜) = 2", len(A.point_ops(TRI)) == 2)
check("point_ops 全部保长",
      all(np.allclose(m @ A.metric(SQ) @ m.T, A.metric(SQ)) for m in ops_sq))
check("point_ops 行列式都是 ±1",
      all(abs(round(np.linalg.det(m))) == 1 for m in ops_sq))

C2Z = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=int)
MIRR = np.array([[1, 0, 0], [0, -1, 0], [0, 0, 1]], dtype=int)
check("rot_kind(I) = translation", A.rot_kind(np.eye(3, dtype=int)) == "translation")
check("rot_kind(-I) = inversion", A.rot_kind(-np.eye(3, dtype=int)) == "inversion")
check("rot_kind(C2z) = proper", A.rot_kind(C2Z) == "proper")
check("rot_kind(镜面) = improper", A.rot_kind(MIRR) == "improper")

Rc = A.cart_rotation(SQ, C2Z)
check("cart_rotation 是正交矩阵", np.allclose(Rc @ Rc.T, np.eye(3), atol=1e-12))
check("cart_rotation 保 det=+1", abs(np.linalg.det(Rc) - 1.0) < 1e-12)
ax, ang = A.rot_axis_angle(SQ, C2Z)
check("rot_axis_angle(C2z) = 90° 绕 z",
      abs(ang - 90.0) < 1e-6 and abs(abs(ax[2]) - 1.0) < 1e-9,
      "轴=%s 角=%.4f" % (np.round(ax, 6).tolist(), ang))
check("rot_axis_angle(镜面) 无轴角", A.rot_axis_angle(SQ, MIRR) == (None, None))

check("nearest_same_element 精确命中 ~ 0",
      A.nearest_same_element(SQ, [[0, 0, 0.5]], ["Fe"], "Fe", np.array([0.0, 0.0, 0.5])) < 1e-12)
check("nearest_same_element 位移 0.1 胞 -> 0.3 Å",
      abs(A.nearest_same_element(SQ, [[0, 0, 0.5]], ["Fe"], "Fe",
                                 np.array([0.1, 0.0, 0.5])) - 0.3) < 1e-12)
check("nearest_same_element 只认同元素",
      A.nearest_same_element(SQ, [[0, 0, 0.5]], ["Sn"], "Fe", np.array([0.0, 0.0, 0.5])) == float("inf"))

# ============================================================ (B) 判定分支
print("--- (B) 判定分支 (合成结构) ---")

r = synth(SQ, ["Fe", "Sn"], [2, 2],
          [[0, 0, 0.5], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.5]], "Fe")
check("平移连接 -> 常规共线反铁磁", r["verdict"] == A.VERDICT_CONVENTIONAL)
check("平移连接: 命中 translation (优先于反演/旋转)",
      "translation" in {o["kind"] for o in r["ops"]},
      str(sorted({o["kind"] for o in r["ops"]})))

r_screw = synth(TET, ["Fe", "O"], [2, 4],
                [[0, 0, 0], [0.5, 0.5, 0.5],
                 [0.3, 0.3, 0.0], [0.7, 0.7, 0.0], [0.8, 0.2, 0.5], [0.2, 0.8, 0.5]], "Fe")
check("旋转连接(4_2 螺旋) -> 几何候选", r_screw["verdict"] == A.VERDICT_CANDIDATE)
check("旋转连接: 给出轴与角",
      any(o.get("angle_deg") for o in r_screw["ops"]),
      str([o.get("angle_deg") for o in r_screw["ops"] if o.get("angle_deg")]))
check("候选的说明写明要能带计算", "能带计算" in r_screw["note"], r_screw["note"])

r = synth(TRI, ["Fe", "Sn"], [2, 2],
          [[0.1, 0.2, 0.3], [0.9, 0.8, 0.7], [0.4, 0.1, 0.6], [0.6, 0.9, 0.4]], "Fe")
check("3D 反演连接 -> 宇称禁止", r["verdict"] == A.VERDICT_PARITY)
check("3D 反演连接: kinds = {inversion}", {o["kind"] for o in r["ops"]} == {"inversion"})

r = synth(SQ, ["Fe", "Sn"], [2, 2],
          [[0.2, 0.3, 0.5], [0.8, 0.7, 0.5], [0.4, 0.1, 0.5], [0.6, 0.9, 0.5]], "Fe")
check("共面胞: 反演与面内 C2 是同一操作 -> 候选 (不是宇称禁止)",
      r["verdict"] == A.VERDICT_CANDIDATE and r["coplanar"] is True)

r = synth(SQ, ["Fe", "Sn"], [2, 2],
          [[0.2, 0.3, 0.5], [0.55, 0.7, 0.5], [0.4, 0.1, 0.5], [0.9, 0.9, 0.5]], "Fe")
check("无连接 -> 不成立", r["verdict"] == A.VERDICT_NONE)
check("无连接也要报「最接近者 + 残差」",
      "best_mismatch_ang" in r and r["best_mismatch_ang"] > 0.1 and r["best_candidate"]["kind"],
      "%s %.4f Å" % (r["best_candidate"]["kind"], r["best_mismatch_ang"]))
check("无连接的 note 指明容差与放宽办法",
      ("Å" in r["note"]) and ("--tol" in r["note"]))

r = synth(SQ, ["Fe", "Sn"], [3, 2],
          [[0, 0, 0.5], [0.5, 0.5, 0.5], [0.25, 0.25, 0.5], [0.3, 0.7, 0.5], [0.8, 0.2, 0.5]], "Fe")
check("3 套等价类 -> 不适用", r["verdict"] == A.VERDICT_NA and r["n_sub"] == 3, r["note"])

r = synth(SQ, ["Fe", "Sn"], [3, 2],
          # 两个磁位点给**同一**组分数坐标 (合成退化胞, 只为触发「原子数不等」这条守卫):
          # 等价类 {0:2 个, 1:1 个} ⇒ 两套亚晶格但原子数不等。
          [[0, 0, 0.5], [0, 0, 0.5], [0.5, 0.5, 0.5], [0.3, 0.7, 0.5], [0.8, 0.2, 0.5]], "Fe")
check("两套亚晶格但原子数不等 -> 不适用",
      r["verdict"] == A.VERDICT_NA and r["sub_sizes"] == [2, 1], str(r["sub_sizes"]))

r = synth(SQ, ["Fe", "Sn"], [2, 2],
          [[0, 0, 0.5], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.5]], "Fe", spins=[1, 1])
check("自旋同号 -> compensated = False", r["compensated"] is False)
check("自旋同号仍走几何分支", r["verdict"] == A.VERDICT_CONVENTIONAL)

r = synth(SQ, ["Fe", "Sn"], [2, 2],
          [[0, 0, 0.5], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.5]], "Fe", spins=[1, -1, 1])
check("spins 给 3 个 -> 不适用 + 说明", r["verdict"] == A.VERDICT_NA and "spins" in r["note"])

r = synth(SQ, ["Fe", "Sn"], [2, 2],
          [[0, 0, 0.5], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0.5]], "Fe")
check("结果永远带边界声明", r["disclaimer"] == A.DISCLAIMER and "不可判定" in r["disclaimer"])
check("候选措辞是「几何候选」而不是结论",
      "几何候选" in A.VERDICT_CANDIDATE and "这是交变磁体" not in A.VERDICT_CANDIDATE)

# ---- 容差语义: 判据按 Å 距离, 不按四舍六入查元组
print("--- (B') 容差语义 ---")
JIT = 3e-6
r = synth(SQ, ["Fe", "Sn"], [2, 2],
          [[0, 0, 0.5], [0, 0.5, 0.5 + JIT], [0.5, 0, 0.5], [0.5, 0.5, 0.5]], "Fe")
check("坐标抖 3e-6 仍认出平移连接 (旧版查元组会落空)",
      r["verdict"] == A.VERDICT_CONVENTIONAL)
check("抖动的残差被量化报告 (< 0.01 Å)",
      r.get("worst_mismatch_ang", 1.0) < 0.01, "%s Å" % r.get("worst_mismatch_ang"))

r_tight = synth(SQ, ["Fe", "Sn"], [2, 2],
                [[0, 0, 0.5], [0, 0.5, 0.5 + JIT], [0.5, 0, 0.5], [0.5, 0.5, 0.5]], "Fe",
                tol=1e-9)
check("--tol 收紧到 1e-9 后同一结构判不出来 (容差真的在用)",
      r_tight["verdict"] == A.VERDICT_NONE
      and r_tight["best_mismatch_ang"] > 1e-9)

# ============================================================ (C) 真 fixture
print("--- (C) 真 fixture ---")
if os.path.exists(POSCAR_MNPS3):
    lat, atoms, frac, _ = A.MO.parse_poscar(POSCAR_MNPS3)
    r_mn = A.check(lat, atoms, frac, "Mn")
    check("MnPS3: 两套亚晶格 (1,1)", r_mn["n_sub"] == 2 and r_mn["sub_sizes"] == [1, 1])
    check("MnPS3: 镜面连接 -> 几何候选 (不是「没有对称操作」)",
          r_mn["verdict"] == A.VERDICT_CANDIDATE,
          "%s ops=%s" % (r_mn["verdict"][:12], sorted({o["kind"] for o in r_mn["ops"]})))
    check("MnPS3: 残差是 mÅ 量级 (纯精度问题, 不是缺对称)",
          r_mn.get("worst_mismatch_ang", 1.0) < 1e-2, "%s Å" % r_mn.get("worst_mismatch_ang"))
else:
    check("MnPS3 fixture 存在", False, POSCAR_MNPS3)

if os.path.exists(POSCAR_CROCL):
    lat, atoms, frac, _ = A.MO.parse_poscar(POSCAR_CROCL)
    r_cr = A.check(lat, atoms, frac, "Cr")
    check("CrOCl: 两套亚晶格", r_cr["n_sub"] == 2)
    check("CrOCl: 反演连接 -> 宇称禁止", r_cr["verdict"] == A.VERDICT_PARITY,
          str(sorted({o["kind"] for o in r_cr["ops"]})))
    check("CrOCl: 反演在场即使还有旋转也判禁止",
          "inversion" in {o["kind"] for o in r_cr["ops"]}
          and len({o["kind"] for o in r_cr["ops"]}) > 1,
          str(sorted({o["kind"] for o in r_cr["ops"]})))
    check("CrOCl: 说明点出「反演在场即已禁止」", "反演在场" in r_cr["note"], r_cr["note"])
else:
    check("CrOCl fixture 存在", False, POSCAR_CROCL)

if os.path.exists(POSCAR_FE3SN):
    lat, atoms, frac, _ = A.MO.parse_poscar(POSCAR_FE3SN)
    r_fe = A.check(lat, atoms, frac, "Fe")
    check("Fe3Sn: 3 套等价类 -> 不适用", r_fe["verdict"] == A.VERDICT_NA and r_fe["n_sub"] == 3)
else:
    check("Fe3Sn fixture 存在", False, POSCAR_FE3SN)

# ============================================================ (D) CLI
print("--- (D) CLI ---")
p_self = run_cli(os.path.join(SCRIPTS, "altermagnet_check.py"), "--selftest")
check("--selftest exit 0 且 ALL PASS",
      p_self.returncode == 0 and "ALL PASS" in p_self.stdout,
      "exit=%d" % p_self.returncode)

p_no = run_cli(os.path.join(SCRIPTS, "altermagnet_check.py"))
check("缺 --poscar -> exit 2", p_no.returncode == 2)

p_bad = run_cli(os.path.join(SCRIPTS, "altermagnet_check.py"), "--poscar", POSCAR_MNPS3,
                "--element", "Xx")
check("未知元素 -> exit 2 且提示可选元素",
      p_bad.returncode == 2 and "没有元素" in (p_bad.stdout + p_bad.stderr))

jp = os.path.join(TD, "mn.json")
p_json = run_cli(os.path.join(SCRIPTS, "altermagnet_check.py"), "--poscar", POSCAR_MNPS3,
                 "--element", "Mn", "--json", jp)
check("--json exit 0", p_json.returncode == 0)
obj = {}
if os.path.exists(jp):
    with open(jp, encoding="utf-8") as fh:
        obj = json.load(fh)
need = {"element", "n_sub", "sub_sizes", "n_mag", "compensated", "ops", "verdict",
        "note", "disclaimer"}
check("--json 键齐全", need <= set(obj), str(sorted(need - set(obj))))
check("--json 与 API 结论一致", obj.get("verdict") == A.VERDICT_CANDIDATE)
check("--json 带残差信息", "worst_mismatch_ang" in obj or "best_mismatch_ang" in obj)

p_spin = run_cli(os.path.join(SCRIPTS, "altermagnet_check.py"), "--poscar", POSCAR_MNPS3,
                 "--element", "Mn", "--spins", "+1,+1", "--json", os.path.join(TD, "sp.json"))
obj2 = {}
if os.path.exists(os.path.join(TD, "sp.json")):
    with open(os.path.join(TD, "sp.json"), encoding="utf-8") as fh:
        obj2 = json.load(fh)
check("--spins +1,+1 -> compensated False", obj2.get("compensated") is False)

p_tol = run_cli(os.path.join(SCRIPTS, "altermagnet_check.py"), "--poscar", POSCAR_MNPS3,
                "--element", "Mn", "--tol", "0.0000001")
check("--tol 极小 -> 报告「不成立」+ 最接近者",
      "不成立" in p_tol.stdout and "最接近者" in p_tol.stdout)

# ============================================================ (E) MC 侧接线
print("--- (E) MC 侧接线 ---")
check("parse_order_spec(auto)", M.parse_order_spec("auto") == ("auto", None))
check("parse_order_spec(sublattice:0,2)",
      M.parse_order_spec("sublattice:0,2") == ("sublattice", (0, 2)))
check("parse_order_spec(staggered)", M.parse_order_spec("staggered") == ("staggered", None))
check("parse_order_spec 未知名原样返回 (由调用方拦)",
      M.parse_order_spec("bogus") == ("bogus", None))
for bad in ("sublattice:", "sublattice:0,1,", "sublattice:x"):
    try:
        M.parse_order_spec(bad)
        ok = False
        msg = "没有报错"
    except ValueError as e:
        ok = True
        msg = str(e)[:40]
    check("parse_order_spec(%r) 报错" % bad, ok, msg)

sg = M.sublattice_signs(2, (4, 4), (0,))
check("sublattice_signs 长度 = N", sg.shape == (2 * 4 * 4,), str(sg.shape))
check("sublattice_signs 一半 +1 一半 -1",
      int((sg > 0).sum()) == sg.size // 2 and int((sg < 0).sum()) == sg.size // 2)
basis_rule = np.where(np.arange(sg.size) % 2 == 0, 1.0, -1.0)
check("sublattice_signs({0}) == 老的 basis 规则", np.array_equal(sg, basis_rule))
sg2 = M.sublattice_signs(3, (2, 2), (1, 2))
check("sublattice_signs 多子晶格掩码", int((sg2 > 0).sum()) == 8 and int((sg2 < 0).sum()) == 4)
for bad_mask, why in (((0, 1), "覆盖全部亚晶格"), ((5,), "越界")):
    try:
        M.sublattice_signs(2, (4, 4), bad_mask)
        ok = False
        msg = "没有报错"
    except ValueError as e:
        ok = True
        msg = str(e)[:32]
    check("sublattice_signs(%s) 报错 (%s)" % (list(bad_mask), why), ok, msg)
check("ORDER_HELP 提到 sublattice 与掩码写法",
      "sublattice" in M.ORDER_HELP and "<mask>" in M.ORDER_HELP)

# ---- 端到端: 老 staggered 与新 sublattice:0 必须逐位一致 (二亚晶格模型)
if os.path.exists(POSCAR_MNPS3):
    model = os.path.join(TD, "mnps3.json")
    r_model = run_cli(os.path.join(SCRIPTS, "poscar_to_mc3.py"), POSCAR_MNPS3,
                      "--element", "Mn", "--S", "2.5", "--J1", "-0.79", "--J2", "-0.04",
                      "--J3", "-0.23", "--MAE", "0.0", "-o", model)
    check("poscar_to_mc3 造 mc 模型", r_model.returncode == 0 and os.path.exists(model),
          r_model.stdout.strip().splitlines()[-1] if r_model.stdout.strip() else "")
    if os.path.exists(model):
        common = ["--size", "4", "4", "--tmin", "5", "--tmax", "40", "--nt", "3",
                  "--therm", "20", "--meas", "50", "--seeds", "1"]
        pa = run_cli(os.path.join(SCRIPTS, "mc_2d.py"), model, "--order", "staggered",
                     *common, "--out", os.path.join(TD, "a"))
        pb = run_cli(os.path.join(SCRIPTS, "mc_2d.py"), model, "--order", "sublattice:0",
                     *common, "--out", os.path.join(TD, "b"))
        check("mc_2d --order staggered exit 0", pa.returncode == 0,
              (pa.stdout + pa.stderr).strip().splitlines()[-1] if pa.returncode else "")
        check("mc_2d --order sublattice:0 exit 0", pb.returncode == 0,
              (pb.stdout + pb.stderr).strip().splitlines()[-1] if pb.returncode else "")
        fa, fb = os.path.join(TD, "a_mt.csv"), os.path.join(TD, "b_mt.csv")
        if os.path.exists(fa) and os.path.exists(fb):
            da = open(fa, "rb").read()
            db = open(fb, "rb").read()
            check("sublattice:0 与 staggered 逐位一致 (同一初态/种子)",
                  da == db, "%d vs %d 字节" % (len(da), len(db)))
            check("CSV 头仍是 4 列契约",
                  da.splitlines()[0].decode() == "T_Kelvin,order_param,Susceptibility,Binder_U4")
        else:
            check("mc_2d 产出 CSV", False, fa)
        check("sublattice 日志写明掩码与分组",
              "掩码" in pb.stdout and "Néel" in pb.stdout)
        check("老的 staggered 日志仍打图二着色说明", "图二着色" in pa.stdout)
        p_guard = run_cli(os.path.join(SCRIPTS, "mc_2d.py"), model, "--order", "sublattice:0,1",
                          "--size", "4", "4", "--nt", "2", "--out", os.path.join(TD, "c"))
        check("覆盖全部亚晶格 -> 非零退出且说明 B 组为空",
              p_guard.returncode != 0 and "B 组为空" in (p_guard.stdout + p_guard.stderr))
        p_fss = run_cli(os.path.join(SCRIPTS, "mc_fss.py"), model, "--order", "sublattice:0,1",
                        "--sizes", "4", "6", "--nt", "2", "--out", os.path.join(TD, "d"))
        check("mc_fss 复用同一解析 (掩码非法 -> 报错)",
              p_fss.returncode != 0 and "B 组为空" in (p_fss.stdout + p_fss.stderr))
        p_cv = run_cli(os.path.join(SCRIPTS, "mc_fss.py"), model, "--order", "cv",
                       "--sizes", "4", "6", "--nt", "2", "--out", os.path.join(TD, "e"))
        check("mc_fss 明确拒绝 cv (U4 交点不适用)",
              p_cv.returncode != 0 and "U4" in (p_cv.stdout + p_cv.stderr))

p_drv = run_cli(os.path.join(HERE, "run_2d_magnet.py"), "--help")
check("驱动 --help 列出 --mc-order", "--mc-order" in p_drv.stdout)
check("驱动 --help 说明掩码写法", "sublattice:<mask>" in p_drv.stdout)

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
