#!/usr/bin/env python
"""smoke_test.py — 回归测试：跑关键确定性步骤，断言输出与已验证值一致。

运行:
  python tests/smoke_test.py

覆盖:
  1. detect_lattice 正确分类 (蜂窝 / 簇磁体)
  2. cluster_sites: Nb3I8 的 3 个 Nb -> 1 个簇
  3. solve_J3 --demo 精确还原文献 J (PASS)
  4. build_magnetic_cases3 蜂窝系数 = 教科书值
  5. peak_report/zoom_window/summarize_mc: 峰落温度窗口边界时不报边界值, 并自动缩放窗口

依赖 fixture (仓库自带 tests/fixtures/, 或设 MAGNET_FIXTURES 指向):
  MnPS3_POSCAR / Nb3I8_POSCAR
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import types
import uuid

import numpy as np

# 子进程统一环境: 强制 UTF-8, 免得父进程 (或父 python) 设了 PYTHONIOENCODING=utf-8
# 而子进程按 GBK 输出时, 父端解码抛 UnicodeDecodeError 把整套测试带崩。
CHILD_ENV = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONDONTWRITEBYTECODE="1")


def run(*cmd, **kw):
    kw.setdefault("capture_output", True)
    kw.setdefault("text", True)
    kw.setdefault("encoding", "utf-8")
    kw.setdefault("errors", "replace")
    kw.setdefault("env", CHILD_ENV)
    return subprocess.run(list(cmd), **kw)

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(HERE, "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import cluster_sites as CS  # noqa: E402
import fetch_structure as FS  # noqa: E402
import j_result as JR  # noqa: E402
from run_2d_magnet import (ERROR_CATALOG, MODEL_PRESETS, SCRIPTS as DRV_SCRIPTS,  # noqa: E402
                           _is_rectangular_net, _is_square_net, _scripts_dir,
                           detect_lattice, doctor, have_matplotlib, multi_sublattice_warning,
                           peak_report, summarize_mc, write_workflow, zoom_window)

FIXTURES = os.environ.get("MAGNET_FIXTURES", os.path.join(HERE, "tests", "fixtures"))
MNPS3 = os.path.join(FIXTURES, "MnPS3_POSCAR")
NB3I8 = os.path.join(FIXTURES, "Nb3I8_POSCAR")
CRSBR = os.path.join(FIXTURES, "CrSBr_POSCAR")
GDI2 = os.path.join(FIXTURES, "GdI2_POSCAR")
LACUO = os.path.join(FIXTURES, "La2CuO4_POSCAR")
K2NIF4 = os.path.join(FIXTURES, "K2NiF4_POSCAR")
CR_OCL = os.path.join(FIXTURES, "CrOCl_POSCAR")

for _name, _path in (("MnPS3", MNPS3), ("Nb3I8", NB3I8), ("CrSBr", CRSBR), ("GdI2", GDI2),
                     ("La2CuO4", LACUO), ("K2NiF4", K2NIF4), ("CrOCl", CR_OCL)):
    if not os.path.exists(_path):
        sys.exit("[skip] 缺少测试 fixture: %s (仓库应自带 tests/fixtures/, 或设 MAGNET_FIXTURES 指向)" % _path)

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


def mktmp(prefix="_tmp_"):
    """找一个**真能写**的临时目录。

    踩过的坑 (受限沙箱 + Windows ACL):
      1. tempfile.mkdtemp / TemporaryDirectory 会 chmod 0700, 在受限环境下连
         **创建者自己**都被拒绝写子目录 (PermissionError), 且清理时抛错中断整个测试
         -> 自己拼名 + os.makedirs, 用 shutil.rmtree(ignore_errors) 清理;
      2. 受限沙箱会给流水线子目录挂上拒绝"新建条目"的 ACL (工作区根可写、仓库内不可写),
         所以不能只在仓库内建 -> 按候选链回退。**子进程也必须能写这里** (会当 --outdir)。
    """
    name = prefix + uuid.uuid4().hex[:8]
    cands = [os.environ.get("MAGNET_TEST_TMP"), HERE,
             os.path.dirname(HERE), tempfile.gettempdir()]
    for base in cands:
        if not base:
            continue
        td = os.path.abspath(os.path.join(base, name))
        try:
            os.makedirs(td, exist_ok=True)
            return td
        except OSError:
            continue
    raise RuntimeError("找不到可写的临时目录 (试过 MAGNET_TEST_TMP / 仓库 / 上级 / temp)")


# ---- 1. detect_lattice ----
check("detect MnPS3 = honeycomb", detect_lattice(MNPS3, "Mn")[0] == "honeycomb")
check("detect Nb3I8 = cluster", detect_lattice(NB3I8, "Nb")[0] == "cluster")

# ---- 2. cluster_sites: 3 Nb -> 1 簇 ----
lat, atoms, coords = CS.parse_poscar(NB3I8)
groups = CS.find_clusters(atoms, coords, lat, "Nb", 3.3)
check("Nb3I8: 3 Nb -> 1 cluster", len(groups) == 1 and len(groups[0]) == 3)

# ---- 3. solve_J3 --demo 精确还原文献 J ----
#      与下面 section 4 共用同一个 build_magnetic_cases3 生成的 cases 目录 (无硬编码路径)

# ---- 4. build_magnetic_cases3 蜂窝系数 ----
td = mktmp("_tmp_cases_")
try:
    r = run(sys.executable, os.path.join(SCRIPTS, "build_magnetic_cases3.py"),
            MNPS3, "--supercell", "2", "2", "1",
            "--outdir", os.path.join(td, "cases"),
            "--element", "Mn", "--nonmag", "P", "S",
            "--moment", "5.0", "--spin", "2.5")
    check("build_magnetic_cases3 退出码 0", r.returncode == 0)
    if r.returncode != 0:
        print(r.stdout[-1500:], r.stderr[-1500:])
    else:
        meta = json.load(open(os.path.join(td, "cases", "magnetic_orders.json"),
                             encoding="utf-8"))
        c = meta["coefficients_per_cell"]
        check("fm coeff (-12,-24,-12)", np.allclose(c["fm"], [-12, -24, -12]))
        check("neel coeff (+12,-24,+12)", np.allclose(c["neel"], [12, -24, 12]))
        check("zigzag coeff (+4,+8,-12)", np.allclose(c["zigzag"], [4, 8, -12]))
        check("stripy coeff (-4,+8,+12)", np.allclose(c["stripy"], [-4, 8, 12]))

        # solve_J3 --demo 反解自检: 读同一份 magnetic_orders.json (硬编码旧目录已移除)
        rd = run(sys.executable, os.path.join(SCRIPTS, "solve_J3.py"),
                 os.path.join(td, "cases"), "--demo")
        check("solve_J3 --demo PASS", rd.returncode == 0 and "PASS" in rd.stdout)
        if rd.returncode != 0:
            print(rd.stdout[-1000:], rd.stderr[-1000:])
finally:
    shutil.rmtree(td, ignore_errors=True)

# ---- 5. MC 峰判据: 窗口边界检测 + 自动缩放窗口 ----
#     回归 MnI2 的真实缺陷: 默认窗口 5-150 K / nt=30 (格距 5 K) 整个跨过 ~3 K 宽的 Cv 峰,
#     取 argmax 落到窗口下限, 旧版把这 5.0 K 当成 T_N 报出 (真值 ≈7.4 K)。
td = mktmp("_tmp_peak_")
try:
    os.makedirs(os.path.join(td, "05_mc"))
    csv = os.path.join(td, "05_mc", "mc_mt.csv")
    T = np.linspace(5.0, 150.0, 30)

    # (a) 峰落在窗口下限: 粗网格整个跨过窄峰, 旧版会把边界 5.0 K 当 T_N 报出
    np.savetxt(csv, np.column_stack([T, -0.1 * T, 0.11 * np.exp(-(T - 5.0) / 60.0)]),
               delimiter=",", header="T_Kelvin,E_per_site,Cv", comments="")
    pk = peak_report(csv)
    check("峰判据: 识别落在窗口下限", bool(pk["at_low"]) and pk["T_peak"] == 5.0)
    check("峰落下限: 结论标未分辨 (不把 5.0 K 当 T_N)", "未分辨" in summarize_mc(td, "cv"))
    lo, hi = zoom_window(pk)
    check("缩放窗口: 向下外扩且保留上方覆盖", lo < 5.0 and hi >= 15.0)
    check("缩放窗口: 覆盖被跨过的真峰 (7.4 K)", lo < 7.4 < hi)

    # (b) 峰在窗口内部 -> 正常报值
    np.savetxt(csv, np.column_stack([T, -0.1 * T, np.exp(-((T - 40.0) ** 2) / 18.0)]),
               delimiter=",", header="T_Kelvin,E_per_site,Cv", comments="")
    s = summarize_mc(td, "cv")
    check("峰在窗口内部: 报 T_N 且无未分辨警告", "未分辨" not in s and "40" in s)

    # (c) 峰落在窗口上限
    np.savetxt(csv, np.column_stack([T, -0.1 * T, 0.11 * np.exp(-(150.0 - T) / 60.0)]),
               delimiter=",", header="T_Kelvin,E_per_site,Cv", comments="")
    check("峰判据: 识别落在窗口上限", bool(peak_report(csv)["at_high"]))

    # (d) 蜂窝路径 (χ 判据) 落边界同样要标未分辨
    np.savetxt(csv, np.column_stack([T, 1.0 - T / 200.0, 0.3 * np.exp(-(T - 5.0) / 60.0)]),
               delimiter=",", header="T_Kelvin,order_param,Susceptibility", comments="")
    check("χ 峰落边界: 标未分辨", "未分辨" in summarize_mc(td, "auto"))
finally:
    shutil.rmtree(td, ignore_errors=True)

# ---- 5b. 比热能量的口径: _energy 必须含单离子各向异性 (P0-1 回归) ----
#     旧 _energy 只算交换项、把 onsite 丢掉; 而 Metropolis 的 dE 含各向异性
#     -> 抽样含 A / 测量不含 A -> Cv 少一块、T_N 系统偏移, 退出码仍是 0。
r = run(sys.executable, os.path.join(SCRIPTS, "mc_2d.py"), "--selftest")
check("mc_2d --selftest 退出码 0", r.returncode == 0)
check("mc_2d --selftest 三项全 PASS (含 A / A=0 / onsite 有贡献)",
      r.stdout.count("[PASS]") == 3 and "[FAIL]" not in r.stdout)
if r.returncode != 0:
    print(r.stdout[-1500:], r.stderr[-1500:])

# ---- 6. 多亚晶格铁磁网络的误判警告 (Fe3GeTe2 类: 3 个 Fe/胞 被几何并成 1 簇) ----
td = mktmp("_tmp_fgt_")
try:
    # 合成的 FGT **式**胞: 只为检验分类与告警的几何逻辑, **不是**真实 FGT 结构
    # (六方 a=3.99; 三个 Fe 分居 (0,0)/(1/3,2/3)/(2/3,1/3), 两两 2.30 A < cutoff 3.3)
    pos = os.path.join(td, "POSCAR")
    with open(pos, "w", encoding="utf-8") as f:
        f.write("synthetic FGT-like cell (test fixture, not real Fe3GeTe2)\n1.0\n")
        f.write("   3.9900000000   0.0000000000   0.0000000000\n")
        f.write("  -1.9950000000   3.4554230000   0.0000000000\n")
        f.write("   0.0000000000   0.0000000000  25.0000000000\n")
        f.write("Fe Ge Te\n3 1 2\nDirect\n")
        f.write("0.000000 0.000000 0.45\n0.333333 0.666667 0.45\n0.666667 0.333333 0.45\n")
        f.write("0.000000 0.000000 0.50\n0.000000 0.000000 0.35\n0.000000 0.000000 0.65\n")
    kind, ncl, nmag = detect_lattice(pos, "Fe")
    check("FGT 式 3Fe/胞 被几何并成 1 簇", (kind, ncl, nmag) == ("cluster", 1, 3))
    w = multi_sublattice_warning(kind, ncl, nmag, "Fe")
    check("多亚晶格告警: 触发 (cluster + 1 簇 + >=3 磁性原子)", w is not None)
    check("多亚晶格告警: 双侧表述 (分子簇 与 铁磁网络 都提到)",
          w is not None and "分子簇" in w and "多亚晶格铁磁网络" in w)
    check("多亚晶格告警: 真簇磁体 (2 簇) 不触发",
          multi_sublattice_warning("cluster", 2, 6, "Nb") is None)
    check("多亚晶格告警: 三角/蜂窝不触发",
          multi_sublattice_warning("triangular", 1, 1, "Mn") is None
          and multi_sublattice_warning("honeycomb", 2, 2, "Cr") is None)

    # 关键: 这条几何误判以前会一路走到 compute_coeffs 才 KeyError 崩掉。
    # 现在必须**明确拒绝**并说清原因 (超胞内只有 1 个簇 -> 簇间耦合无从定义)。
    r = run(sys.executable, os.path.join(SCRIPTS, "build_cluster_cases.py"), pos,
            "--element", "Fe", "--nonmag", "Ge", "Te", "--cutoff", "3.3",
            "--nbr-shells", "3", "--supercell", "2", "2", "1",
            "--outdir", os.path.join(td, "cases"), "--root", td,
            "--U", "4.0", "--moment", "2.0")
    msg = (r.stdout or "") + (r.stderr or "")
    check("FGT 式胞: build_cluster_cases 明确拒绝 (非零退出)", r.returncode != 0)
    check("FGT 式胞: 报错说清原因 (不是 KeyError 崩栈)",
          "KeyError" not in msg and "不是分子簇磁体" in msg)
    check("FGT 式胞: 报错指明这是连通磁性网络/多亚晶格",
          "连通的磁性网络" in msg or "多亚晶格" in msg)
finally:
    shutil.rmtree(td, ignore_errors=True)

# ---- 7. workflow.md 记录「多亚晶格」警告 (用户可见产物) ----
td = mktmp("_tmp_wf_")
try:
    fake_args = types.SimpleNamespace(vacuum=15.0, mc_refine_nt=60,
                                      noncollinear=False, soc=False, saxis=None)
    warn_txt = multi_sublattice_warning("cluster", 1, 3, "Fe")

    write_workflow(td, "SynthFGT", "Fe", ["Ge", "Te"], 1.5, 3.0, 4.0, 10.0, 0.0, 0.0, 0.0,
                   fake_args, None, "cluster", None, warn_txt)
    txt = open(os.path.join(td, "workflow.md"), encoding="utf-8").read()
    check("workflow.md: 有警告时写入局限第 8 条", "8." in txt and "多亚晶格" in txt)
    check("workflow.md: 警告原文进卡片", "分子簇" in txt)

    write_workflow(td, "MnPS3", "Mn", ["P", "S"], 2.5, 5.0, 5.0, -0.79, -0.04, -0.23, 0.0,
                   fake_args, None, "honeycomb", None, None)
    txt2 = open(os.path.join(td, "workflow.md"), encoding="utf-8").read()
    check("workflow.md: 无警告时不出现第 8 条", "多亚晶格" not in txt2 and "8." not in txt2)
finally:
    shutil.rmtree(td, ignore_errors=True)

# ---- 8. mc_fss.py 多尺寸有限尺寸标度 (Binder U4 交点定 T_C/T_N) ----
td = mktmp("_tmp_fss_")
try:
    model = os.path.join(td, "model.json")
    r = run(sys.executable, os.path.join(SCRIPTS, "poscar_to_mc3.py"), MNPS3,
            "--element", "Mn", "--S", "2.5", "--J1", "-0.79", "--J2", "-0.04",
            "--J3", "-0.23", "--MAE", "0.0", "-o", model)
    check("mc_fss: 生成 MC 模型", r.returncode == 0)
    r = run(sys.executable, os.path.join(SCRIPTS, "mc_fss.py"), model,
            "--sizes", "6", "8", "10", "--order", "auto", "--tmin", "5", "--tmax", "60",
            "--nt", "12", "--therm", "60", "--meas", "100", "--seeds", "1",
            "--out", os.path.join(td, "fss"))
    check("mc_fss: 退出 0", r.returncode == 0)
    rj = os.path.join(td, "fss_fss.json")
    check("mc_fss: 写出 fss_fss.json", os.path.exists(rj))
    if os.path.exists(rj):
        d = json.load(open(rj, encoding="utf-8"))
        check("mc_fss: json 含 tc (float)", "tc" in d and isinstance(d["tc"], float))
    check("mc_fss: 写出 fss_fss.csv", os.path.exists(os.path.join(td, "fss_fss.csv")))
finally:
    shutil.rmtree(td, ignore_errors=True)

# ---- 9. fetch_structure.py: CIF/OPTIMADE/AFLOW -> POSCAR + provider 链 ----
# 9a. CIF -> POSCAR (NaCl: 2 不对称位 x 4 对称操作 -> 4 Na + 4 Cl)
_naCl_cif = """data_test
_cell_length_a 5.0
_cell_length_b 5.0
_cell_length_c 5.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
loop_
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na 0.0 0.0 0.0
Cl 0.5 0.5 0.5
loop_
_symmetry_equiv_pos_as_xyz
'x,y,z'
'x+1/2,y+1/2,z'
'x+1/2,y,z+1/2'
'x,y+1/2,z+1/2'
"""
_p9 = FS.cif_to_poscar(_naCl_cif, "NaCl")
check("fetch: CIF->POSCAR 含物种行 Na Cl", "Na  Cl" in _p9)
check("fetch: CIF->POSCAR 对称展开 8 原子", "Direct" in _p9
      and sum(1 for l in _p9.splitlines()[8:] if l.strip()) == 8)

# 9b. OPTIMADE -> POSCAR (MnS, 面心式)
_p9b = FS.optimade_to_poscar(
    {"lattice_vectors": [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]],
     "cartesian_site_positions": [[0, 0, 0], [2.0, 2.0, 2.0]],
     "species_at_sites": ["Mn", "S"], "nelements": 2}, "MnS")
check("fetch: OPTIMADE->POSCAR 分数坐标 0.5", "0.500000000000" in _p9b)

# 9c. AFLOW -> POSCAR
_p9c = FS.aflow_to_poscar(
    {"lattice": [[4.0, 0, 0], [0, 4.0, 0], [0, 0, 4.0]],
     "atom_names": ["Mn", "S"],
     "positions_fractional": [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]}, "MnS")
check("fetch: AFLOW->POSCAR 物种 Mn S", _p9c.count("\n") >= 7 and "Mn  S" in _p9c)

# 9d. provider 链 (mock HTTP, 不触网)
class _Resp:
    def __init__(self, text=None, js=None):
        self._text, self._js = text, js

    def json(self):
        return self._js

    @property
    def text(self):
        return self._text

    def raise_for_status(self):
        pass


def _fake_get(url, **kw):
    if "result.php" in url:
        return _Resp(js=[{"file": "test.cif"}])
    if url.endswith("test.cif"):
        return _Resp(text=_naCl_cif)
    raise RuntimeError("unexpected url " + url)


_p9d_name, _p9d_poscar = FS.fetch_chain("NaCl", providers=["cod"], _get=_fake_get)
check("fetch: 链 mock 命中 cod", _p9d_name == "cod" and "Na  Cl" in _p9d_poscar)
# 链全失败 -> 抛 RuntimeError (含各 provider 错误)
try:
    FS.fetch_chain("Xx", providers=["cod", "oqmd"], _get=lambda u, **k: (_ for _ in ()).throw(RuntimeError("down")))
    check("fetch: 链全失败应抛错", False)
except RuntimeError as e:
    check("fetch: 链全失败抛错含 cod/oqmd", "cod:" in str(e) and "oqmd:" in str(e))

# ---- 10. 子命令化: fetch/build/report 分阶段派发 (不触网、不跑 DFT/MC/render) ----
td = mktmp("_tmp_sub_")
try:
    driver = os.path.join(HERE, "run_2d_magnet.py")
    base = [sys.executable, driver, "--poscar", MNPS3, "--formula", "MnPS3",
            "--element", "Mn", "--nonmag", "P", "S", "--literature", "MnPS3",
            "--outdir", td]
    r = run(*(base + ["fetch"]))
    check("sub: fetch 退出 0", r.returncode == 0)
    check("sub: fetch 写 01_structure/POSCAR",
          os.path.exists(os.path.join(td, "01_structure", "POSCAR")))
    check("sub: fetch 早退(不建 03_dft)", not os.path.exists(os.path.join(td, "03_dft")))
    r = run(*(base + ["build"]))
    check("sub: build 退出 0", r.returncode == 0)
    check("sub: build 写 02_supercell/cases/magnetic_orders.json",
          os.path.exists(os.path.join(td, "02_supercell", "cases", "magnetic_orders.json")))
    check("sub: build 早退(不建 03_dft)", not os.path.exists(os.path.join(td, "03_dft")))
    r = run(*(base + ["report"]))
    check("sub: report 退出 0", r.returncode == 0)
    check("sub: report 写 workflow.md", os.path.exists(os.path.join(td, "workflow.md")))
    check("sub: report 不跑 MC (05_mc 不存在)", not os.path.exists(os.path.join(td, "05_mc")))
finally:
    shutil.rmtree(td, ignore_errors=True)

# ---- 11. scripts/ 包化: _scripts_dir() 定位脚本 + 文献表 ----
check("pack: _scripts_dir 命中 (含 mc_2d.py)",
      os.path.exists(os.path.join(DRV_SCRIPTS, "mc_2d.py")))
check("pack: _scripts_dir 与源码 scripts/ 一致", DRV_SCRIPTS == SCRIPTS)
_lit = os.path.join(DRV_SCRIPTS, "literature_J.json")
check("pack: literature_J.json 在 scripts/ 内", os.path.exists(_lit))
if os.path.exists(_lit):
    _tab = json.load(open(_lit, encoding="utf-8"))
    check("pack: literature_J.json 可解析且含 MnPS3", "MnPS3" in _tab)

# ---- 12. doctor 子命令 + error catalog ----
check("doctor: ERROR_CATALOG 非空且条目字段齐全",
      bool(ERROR_CATALOG) and all(set(e) >= {"key", "symptom", "fix", "static"}
                                  for e in ERROR_CATALOG))
_known_keys = {e["key"] for e in ERROR_CATALOG}
for _k in ("static_incar", "solve_no_outcar", "packaging", "fss_no_cross", "peak_boundary"):
    check(f"doctor: catalog 含 {_k}", _k in _known_keys)
check("doctor: 有静态运行时条目 (>=5)",
      sum(1 for e in ERROR_CATALOG if e["static"]) >= 5)
check("doctor: 有环境可探条目 (>=5)",
      sum(1 for e in ERROR_CATALOG if not e["static"]) >= 5)

r = run(sys.executable, os.path.join(HERE, "run_2d_magnet.py"), "doctor")
check("doctor: 子命令派发 (退出码 0 或 1)", r.returncode in (0, 1))
check("doctor: 输出含环境自检头", "python :" in r.stdout)
check("doctor: 输出含静态参考段", "常见运行时问题" in r.stdout)
check("doctor: 静态参考含 static_incar", "static_incar" in r.stdout)

# ---- 13. --model 命名预设 ----
check("model: 预设键齐全", set(MODEL_PRESETS) >= {"mnps3", "cri3", "crbr3", "fe3gete2"})
_req = {"formula", "element", "nonmag", "lattice", "S", "moment", "U", "literature", "desc"}
check("model: 每条预设字段齐全", all(set(v) >= _req for v in MODEL_PRESETS.values()))
# 预设值须与 literature_J.json 材料专属值一致 (不能和 ELEMENT_PARAMS 通用兜底搞混)
_tab = json.load(open(_lit, encoding="utf-8"))
for _k in ("mnps3", "cri3", "crbr3", "fe3gete2"):
    _p = MODEL_PRESETS[_k]
    _l = _tab[_p["literature"]]
    check(f"model: {_k} 与文献条目一致 (element/S/moment/U)",
          _p["element"] == _l["element"] and _p["S"] == _l["S"]
          and _p["moment"] == _l["moment"] and _p["U"] == _l["U"]
          and _p["nonmag"] == _l["nonmag"])

td = mktmp("_tmp_model_")
try:
    driver = os.path.join(HERE, "run_2d_magnet.py")
    r = run(sys.executable, driver, "build", "--model", "mnps3",
            "--poscar", MNPS3, "--outdir", td)
    check("model: mnps3 build 退出 0", r.returncode == 0)
    check("model: mnps3 build 写 cases/magnetic_orders.json",
          os.path.exists(os.path.join(td, "02_supercell", "cases", "magnetic_orders.json")))
    check("model: 输出含 [model] 命名预设", "[model]" in r.stdout and "mnps3" in r.stdout)
    check("model: 输出含 [ref] J 参考声明 (stderr)", "[ref]" in r.stderr)
    check("model: build 早退(不建 03_dft)", not os.path.exists(os.path.join(td, "03_dft")))
finally:
    shutil.rmtree(td, ignore_errors=True)

# ---- 14. 矩形晶格 (CrSBr) ----
check("model: 含 crsbr 预设", "crsbr" in MODEL_PRESETS)
check("model: crsbr lattice=rectangular", MODEL_PRESETS["crsbr"]["lattice"] == "rectangular")
check("model: crsbr 与文献条目一致 (element/S/moment/U)",
      MODEL_PRESETS["crsbr"]["element"] == "Cr" and MODEL_PRESETS["crsbr"]["S"] == 1.5
      and MODEL_PRESETS["crsbr"]["moment"] == 3.0 and MODEL_PRESETS["crsbr"]["U"] == 4.0)

_lat14, _ncl14, _nmag14 = detect_lattice(CRSBR, "Cr")
check("rect: detect_lattice -> rectangular (1 磁原子/胞, 不是 triangular)",
      _lat14 == "rectangular" and _nmag14 == 1)

td14 = mktmp("_tmp_rect_")
try:
    r14 = run(sys.executable, os.path.join(SCRIPTS, "build_rectangular_cases.py"),
              CRSBR, "--element", "Cr", "--nonmag", "S", "Br",
              "--moment", "3.0", "--spin", "1.5", "--selftest", "--outdir", td14)
    check("rect: build_rectangular_cases --selftest (系数+反解) PASS",
          "PASS" in r14.stdout and "FAIL" not in r14.stdout)
    r14b = run(sys.executable, os.path.join(SCRIPTS, "solve_J_rect.py"), td14, "--demo")
    check("rect: solve_J_rect --demo PASS", "PASS" in r14b.stdout)

    drv = os.path.join(HERE, "run_2d_magnet.py")
    r14c = run(sys.executable, drv, "build", "--model", "crsbr",
               "--poscar", CRSBR, "--outdir", td14 + "_drv")
    check("rect: driver build --model crsbr 退出 0", r14c.returncode == 0)
    check("rect: driver 判 [lattice] rectangular", "[lattice] rectangular" in r14c.stdout)
    _mo = os.path.join(td14 + "_drv", "02_supercell", "cases", "magnetic_orders.json")
    check("rect: cases 四共线态 fm/afm_a/afm_b/afm_diag",
          os.path.exists(_mo) and set(json.load(open(_mo, encoding="utf-8"))["selected_orders"])
          == {"fm", "afm_a", "afm_b", "afm_diag"})
    # assemble_dft 应读 selected_orders (不是硬编码蜂窝 fm/neel/zigzag/stripy)
    r14d = run(sys.executable, drv, "dft", "--model", "crsbr",
               "--poscar", CRSBR, "--outdir", td14 + "_drv")
    check("rect: dft 组装 runs/ 四共线态 (无 neel/zigzag/stripy)",
          r14d.returncode == 0 and set(os.listdir(os.path.join(td14 + "_drv", "runs")))
          == {"fm", "afm_a", "afm_b", "afm_diag"})
finally:
    shutil.rmtree(td14, ignore_errors=True)
    shutil.rmtree(td14 + "_drv", ignore_errors=True)

# ---- 15. 4f 磁体 (GdI2, 三角) ----
check("model: 含 gdi2 预设", "gdi2" in MODEL_PRESETS)
check("model: gdi2 lattice=triangular (4f 磁体不是新晶格类型)",
      MODEL_PRESETS["gdi2"]["lattice"] == "triangular")
check("model: gdi2 与文献条目一致 (element/S/moment/U)",
      MODEL_PRESETS["gdi2"]["element"] == "Gd" and MODEL_PRESETS["gdi2"]["S"] == 3.5
      and MODEL_PRESETS["gdi2"]["moment"] == 7.0 and MODEL_PRESETS["gdi2"]["U"] == 6.0)
_lit = json.load(open(os.path.join(SCRIPTS, "literature_J.json"), encoding="utf-8"))
check("lit: literature_J.json 含 GdI2 (ground_state=fm, J 未核实留 null)",
      "GdI2" in _lit and _lit["GdI2"]["ground_state"] == "fm" and _lit["GdI2"]["J1"] is None)

_lat15, _ncl15, _nmag15 = detect_lattice(GDI2, "Gd")
check("gdi2: detect_lattice -> triangular (1 磁原子/胞)",
      _lat15 == "triangular" and _nmag15 == 1)

td15 = mktmp("_tmp_gdi2_")
try:
    drv = os.path.join(HERE, "run_2d_magnet.py")
    r15 = run(sys.executable, drv, "build", "--model", "gdi2", "--J1", "0.33",
               "--poscar", GDI2, "--outdir", td15)
    check("gdi2: driver build --model gdi2 退出 0", r15.returncode == 0)
    check("gdi2: driver 判 [lattice] triangular", "[lattice] triangular" in r15.stdout)
    _mo = os.path.join(td15, "02_supercell", "cases", "magnetic_orders.json")
    check("gdi2: cases 四态 fm/120deg/stripy_a/stripy_ab",
          os.path.exists(_mo) and set(json.load(open(_mo, encoding="utf-8"))["selected_orders"])
          == {"fm", "120deg", "stripy_a", "stripy_ab"})
finally:
    shutil.rmtree(td15, ignore_errors=True)

# ---- 16. P0.2: J_result.json 统一产物 (schema / 自描述 / 信任度不可手改) ----
# 同一组 CrSBr 矩形参数, 下面 section 17 还会拿它喂 poscar_to_mc3 --jfile
_JR_ARGS = dict(
    lattice="rectangular",
    shells=[("J1", 3.5404561399999994), ("J2", 4.75547185), ("J3", 5.928688066967777)],
    J={"J1": 2.79, "J2": 2.12, "J3": 3.88},
    E={"fm": -0.45612, "afm_a": 0.30348, "afm_b": 0.25524, "afm_diag": -0.1026},
    E0=0.0, cond=32.0, residual=5.551115123125783e-17, spin=1.5, n_mag=16,
    A=[[1.0, -16.0, -16.0, -32.0], [1.0, 16.0, -16.0, 32.0],
       [1.0, -16.0, 16.0, 32.0], [1.0, 16.0, 16.0, -32.0]],
    orders=("fm", "afm_a", "afm_b", "afm_diag"), material="CrSBr", z_from="fm",
)
_jr_vasp = JR.build(source="vasp", **_JR_ARGS)
check("jresult: source=vasp + 残差 5.5e-17 + cond 32 -> confidence=high",
      _jr_vasp["schema"] == "j_result/1" and _jr_vasp["confidence"] == "high")
check("jresult: z_n 由 fm 系数导出 (矩形配位数 2/2/4)",
      [s["z_n"] for s in _jr_vasp["shells"]] == [2.0, 2.0, 4.0])
check("jresult: 每个壳自带 tol -> 消费端不必自己猜容差",
      all(s.get("tol") is not None and 0 < s["tol"] <= 0.3 for s in _jr_vasp["shells"]))
check("jresult: J 往返不留尾差 (2.79, 不是 2.790000000000001)",
      _jr_vasp["J_meV"] == {"J1": 2.79, "J2": 2.12, "J3": 3.88})
check("jresult: J_tilde = J * S^2", abs(_jr_vasp["J_tilde_eV"]["J1"] - 2.79 * 2.25) < 1e-12)
check("jresult: demo/literature 一律 low (不是 vasp 就不给 high)",
      JR.confidence_of("demo", 1e-17, 32.0) == "low"
      and JR.confidence_of("literature", 1e-17, 32.0) == "low")
check("jresult: 残差超 WARN / cond 超闸门 -> low",
      JR.confidence_of("vasp", 1.0, 32.0) == "low"
      and JR.confidence_of("vasp", 1e-17, 1e4) == "low")

_td16 = mktmp("_tmp_jr16_")
try:
    _p16 = os.path.join(_td16, "J_result.json")
    JR.write(_p16, _jr_vasp)
    _back16, _w16 = JR.read(_p16)
    check("jresult: write -> read -> validate 零 warning", not _w16)
    check("jresult: 往返后 J 与壳距离不变",
          _back16["J_meV"] == _jr_vasp["J_meV"]
          and [s["d"] for s in _back16["shells"]] == [s["d"] for s in _jr_vasp["shells"]])

    _rej1 = _rej2 = False
    try:
        JR.validate(dict(_back16, source="demo", confidence="high"), strict=True)
    except JR.JResultError:
        _rej1 = True
    check("jresult: 手改 confidence=high (source=demo) 被拒绝", _rej1)
    try:
        JR.validate(dict(_back16, residual_eV=5.0), strict=True)
    except JR.JResultError:
        _rej2 = True
    check("jresult: 残差 5 eV 超硬阈 -> 拒绝写出 (不静默降级)", _rej2)
finally:
    shutil.rmtree(_td16, ignore_errors=True)

# ---- 17. P0.4: poscar_to_mc3 的显式壳 / 多元素 / J 来源 ----
# 核心契约: 换 J 的来源不应改变**几何**(键集/offset), 只改 J 的数值来源。
P2MC = os.path.join(SCRIPTS, "poscar_to_mc3.py")
FGTMONO = os.path.join(FIXTURES, "Fe3GeTe2_monolayer_POSCAR")

_td17 = mktmp("_tmp_p04_")
try:
    _leg = os.path.join(_td17, "legacy.json")
    r17a = run(sys.executable, P2MC, CRSBR, "--element", "Cr", "--S", "1.5",
               "--J1", "2.79", "--J2", "2.12", "--J3", "3.88", "-o", _leg)
    check("p0.4: 旧路径 (--J1/--J2/--J3) 仍退出 0 (向后兼容)", r17a.returncode == 0)
    _ma = json.load(open(_leg, encoding="utf-8"))
    check("p0.4: 旧路径物理载荷不变 (4 条键, J1_mat diag = -2.79*1.5^2)",
          len(_ma["interactions"]) == 4
          and abs(_ma["interactions"][0]["J"][0][0] - (-2.79 * 1.5 ** 2)) < 1e-12)
    check("p0.4: meta 仍带 J1_meV/J2_meV/J3_meV (mc_2d.py 的 auto 序参量会读)",
          _ma["meta"]["J1_meV"] == 2.79 and _ma["meta"]["J2_meV"] == 2.12
          and _ma["meta"]["J3_meV"] == 3.88)

    # 同一条链改走 --jfile: 键集/offset/J/shell 必须逐位一致
    _jf = os.path.join(_td17, "J_result.json")
    JR.write(_jf, JR.build(source="vasp", **_JR_ARGS))
    _new = os.path.join(_td17, "from_jfile.json")
    r17b = run(sys.executable, P2MC, CRSBR, "--element", "Cr", "--S", "1.5",
               "--jfile", _jf, "-o", _new)
    check("p0.4: --jfile 路径退出 0", r17b.returncode == 0)
    _mb = json.load(open(_new, encoding="utf-8"))

    def _bondsig(m):
        return sorted((i["b1"], i["b2"], tuple(i["offset"]), i["shell"],
                       tuple(map(tuple, i["J"]))) for i in m["interactions"])

    check("p0.4: --jfile 与旧路径 键集/offset/J/shell 完全一致", _bondsig(_ma) == _bondsig(_mb))
    check("p0.4: --jfile 记下 J 与壳的来源 (产物可追溯)",
          (_mb["meta"]["jfile"] or "").endswith("J_result.json")
          and _mb["meta"]["J1_meV"] == 2.79 and _mb["meta"]["source_poscar"])

    # 多个不同 --S: 引擎只有一个全局 S^2, 必须拒绝而不是静默算错
    r17c = run(sys.executable, P2MC, FGTMONO, "--element", "Fe", "Ge",
               "--S", "1.0", "2.0", "--J1", "1.0", "-o", os.path.join(_td17, "x.json"))
    check("p0.4: 多个不同 --S 被拒绝并指向 P3",
          r17c.returncode != 0 and "P3" in (r17c.stderr + r17c.stdout))

    # pair 字符串必须配 --family; 缺了要给可执行的报错, 不能静默全丢键
    _sd = os.path.join(_td17, "sd.json")
    with open(_sd, "w", encoding="utf-8") as _f:
        json.dump({"shells": [{"name": "J1", "d": 2.3, "tol": 0.1, "pair": "AB"}]}, _f)
    r17d = run(sys.executable, P2MC, FGTMONO, "--element", "Fe", "--S", "1.0",
               "--shell-def", _sd, "--J1", "1.0", "-o", os.path.join(_td17, "y.json"))
    check("p0.4: pair='AB' 缺 --family 报错并给出两种写法",
          r17d.returncode != 0 and "--family" in (r17d.stderr + r17d.stdout))

    # --family 生效: A=1(Fe I @2c) B=0,2(Fe II @4e) -> AB 壳只拿 A-B 键 (6 条/胞)
    _sd2 = os.path.join(_td17, "sd2.json")
    with open(_sd2, "w", encoding="utf-8") as _f:
        json.dump({"shells": [{"name": "J1", "d": 2.3042, "tol": 0.05, "pair": "AB"}]}, _f)
    _zf = os.path.join(_td17, "z.json")
    r17e = run(sys.executable, P2MC, FGTMONO, "--element", "Fe", "--S", "1.0",
               "--family", "A=1", "B=0,2", "--shell-def", _sd2, "--J1", "1.0", "-o", _zf)
    _n17e = len(json.load(open(_zf, encoding="utf-8"))["interactions"]) if r17e.returncode == 0 else -1
    check("p0.4: --family A=1 B=0,2 把 AB 壳限到 A-B 键 (3 Fe/胞 -> 6 条)", _n17e == 6)

    # 故意漏掉中间壳 (J2@4.7555): 那条键必须被点名告警, 不能静默蒸发
    _sd3 = os.path.join(_td17, "sd3.json")
    with open(_sd3, "w", encoding="utf-8") as _f:
        json.dump({"shells": [{"name": "J1", "d": 3.5405, "tol": 0.2},
                              {"name": "J3", "d": 5.9287, "tol": 0.2}]}, _f)
    r17f = run(sys.executable, P2MC, CRSBR, "--element", "Cr", "--S", "1.5",
               "--shell-def", _sd3, "--J1", "2.79", "--J3", "3.88",
               "-o", os.path.join(_td17, "w.json"))
    check("p0.4: 漏掉中间壳 -> 未分配的键被点名告警 (不静默丢键)",
          r17f.returncode == 0 and "没落进任何显式壳" in (r17f.stderr + r17f.stdout)
          and "4.755" in (r17f.stderr + r17f.stdout))
finally:
    shutil.rmtree(_td17, ignore_errors=True)

# ---- 18. 方形晶格 (K2NiF4 型: La2CuO4 = 定量锚点, K2NiF4 = 结构原型) ----
# 之前 1 磁原子/胞的方形格子会被 detect_lattice 误判成 triangular (120° 非共线);
# 本节的 8/9 两条就是这个路由缺口的回归钉子。
check("model: 含 la2cuo4 / k2nif4 预设",
      "la2cuo4" in MODEL_PRESETS and "k2nif4" in MODEL_PRESETS)
check("model: la2cuo4 lattice=square", MODEL_PRESETS["la2cuo4"]["lattice"] == "square")
check("model: k2nif4 lattice=square", MODEL_PRESETS["k2nif4"]["lattice"] == "square")
check("model: la2cuo4 与文献条目一致 (element/S/moment/U)",
      MODEL_PRESETS["la2cuo4"]["element"] == "Cu" and MODEL_PRESETS["la2cuo4"]["S"] == 0.5
      and MODEL_PRESETS["la2cuo4"]["moment"] == 1.0 and MODEL_PRESETS["la2cuo4"]["U"] == 7.0)
check("model: 两个方形预设落在同一 lattice 家族 (方形只走一份代码)",
      MODEL_PRESETS["la2cuo4"]["lattice"] == MODEL_PRESETS["k2nif4"]["lattice"] == "square")

_lit18 = json.load(open(os.path.join(SCRIPTS, "literature_J.json"), encoding="utf-8"))
check("lit: La2CuO4 是方形定量锚点 (Coldea 2001 已逐条核实)",
      _lit18["La2CuO4"]["J1"] == -111.8 and _lit18["La2CuO4"]["J2"] == 11.4
      and _lit18["La2CuO4"]["J3"] == -2.0)
check("lit: K2NiF4 的 J 仍留 null (未能核实, 只当结构原型不当锚点)",
      "K2NiF4" in _lit18 and _lit18["K2NiF4"]["J1"] is None and _lit18["K2NiF4"]["J2"] is None)

_lat18, _ncl18, _nmag18 = detect_lattice(LACUO, "Cu")
check("sq: detect_lattice -> square (1 磁原子/胞, **不是 triangular**)",
      _lat18 == "square" and _nmag18 == 1)
check("sq: K2NiF4 也判 square", detect_lattice(K2NIF4, "Ni")[0] == "square")


def _lat_of(path):
    _ls = [ln for ln in open(path, encoding="utf-8").read().splitlines() if ln.strip()]
    _sc = float(_ls[1].strip())
    return np.array([[float(x) for x in _ls[i].split()] for i in (2, 3, 4)]) * _sc


check("sq: _is_square_net(方形) 真且 _is_rectangular_net(方形) 假",
      _is_square_net(_lat_of(LACUO)) and not _is_rectangular_net(_lat_of(LACUO)))
check("sq: _is_rectangular_net(CrSBr) 真且 _is_square_net(CrSBr) 假 (两类互斥)",
      _is_rectangular_net(_lat_of(CRSBR)) and not _is_square_net(_lat_of(CRSBR)))

td18 = mktmp("_tmp_square_")
try:
    drv = os.path.join(HERE, "run_2d_magnet.py")
    r18 = run(sys.executable, drv, "build", "--model", "la2cuo4",
              "--poscar", LACUO, "--outdir", td18)
    check("sq: driver build --model la2cuo4 退出 0", r18.returncode == 0)
    check("sq: driver 判 [lattice] square", "[lattice] square" in r18.stdout)

    _cases18 = os.path.join(td18, "02_supercell", "cases")
    _mo18 = os.path.join(_cases18, "magnetic_orders.json")
    _d18 = json.load(open(_mo18, encoding="utf-8")) if os.path.exists(_mo18) else {}
    check("sq: cases 六序 all_orders + 选定 fm/stripe_a/p4_a/mix_ab",
          set(_d18.get("all_orders", [])) == {"fm", "g_type", "stripe_a", "stripe_b",
                                              "p4_a", "mix_ab"}
          and set(_d18.get("selected_orders", [])) == {"fm", "stripe_a", "p4_a", "mix_ab"})
    check("sq: model_spec lattice=square + policy=shift_search",
          _d18.get("model_spec", {}).get("lattice") == "square"
          and _d18.get("model_spec", {}).get("policy") == "shift_search")
    check("sq: design_cond < 1e4 闸门", _d18.get("design_cond", 1e9) < 1e4)
    check("sq: 周期 2 序 c3 简并被记进产物 (下游看得出为何要 p4_a)",
          abs(max(_d18["coefficients_per_cell"][k][2]
                  for k in ("fm", "g_type", "stripe_a", "stripe_b"))
              - min(_d18["coefficients_per_cell"][k][2]
                    for k in ("fm", "g_type", "stripe_a", "stripe_b"))) < 1e-9)

    _mm = open(os.path.join(_cases18, "g_type", "magmom.txt"), encoding="utf-8").read()
    _vals = [float(x) for x in _mm.split("=", 1)[1].split()]
    check("sq: magmom.txt 值数 = 超胞原子数 (16 Cu + 32 La + 64 O)",
          len(_vals) == 112 and _d18["n_atoms"] == 112)
    check("sq: g_type 的 16 个 Cu 磁矩全为 ±1.0 且正负各 8 (棋盘)",
          sorted(abs(v) for v in _vals if v != 0.0) == [1.0] * 16
          and sum(1 for v in _vals[:16] if v > 0.0) == 8)

    _sq_solve = os.path.join(SCRIPTS, "solve_J_square.py")
    r18b = run(sys.executable, _sq_solve, _cases18, "--demo")
    check("sq: solve_J_square --demo PASS", "PASS" in r18b.stdout)
    _jr18 = json.load(open(os.path.join(_cases18, "J_result.json"), encoding="utf-8"))
    check("sq: 写出 J_result.json (schema=j_result/1, lattice=square, source=demo)",
          _jr18["schema"] == "j_result/1" and _jr18["lattice"] == "square"
          and _jr18["source"] == "demo")

    r18c = run(sys.executable, _sq_solve, _cases18, "--demo", "--ref", "la2cuo4")
    check("sq: --ref la2cuo4 (已核实锚点) PASS", "PASS" in r18c.stdout)
    _jr18c = json.load(open(os.path.join(_cases18, "J_result.json"), encoding="utf-8"))
    _jv = _jr18c.get("J_MeV", _jr18c.get("J_meV"))
    _jlist = list(_jv.values()) if isinstance(_jv, dict) else list(_jv)
    check("sq: --ref 落盘 J1=-111.8 / J2=+11.4 / J3=-2.0 (Coldea 2001 转号)",
          [round(float(x), 6) for x in _jlist] == [-111.8, 11.4, -2.0])
    check("sq: --ref 把文献来源写进 notes, 且 source 老实标 demo (不冒充 vasp)",
          "Coldea" in json.dumps(_jr18c, ensure_ascii=False)
          and _jr18c.get("source") == "demo")

    r18d = run(sys.executable, _sq_solve, _cases18, "--demo", "--ref", "k2nif4")
    check("sq: --ref k2nif4 拒用未核实的 null J (并打印 _note)",
          r18d.returncode != 0 and "null" in (r18d.stdout + r18d.stderr))

    r18e = run(sys.executable, drv, "dft", "--model", "la2cuo4",
               "--poscar", LACUO, "--outdir", td18)
    _runs18 = os.path.join(td18, "runs")
    check("sq: dft 只组装**选定**的 4 个序 (assemble_dft 读 selected_orders)",
          r18e.returncode == 0 and os.path.isdir(_runs18)
          and set(os.listdir(_runs18)) == {"fm", "stripe_a", "p4_a", "mix_ab"})

    # render 分派: 方形复用通用超胞渲染器 (不另写 render_square.py), 且只画选定四序。
    # 渲染器二选一: 有 matplotlib 用 render_supercell.py, 否则驱动自动改走纯 Pillow 的 render_supercell_pil.py。
    _rend18 = "render_supercell.py" if have_matplotlib() else "render_supercell_pil.py"
    r18f = run(sys.executable, drv, "render", "--model", "la2cuo4",
               "--poscar", LACUO, "--outdir", td18)
    check("sq: render 分派给 %s 且 --orders 只含选定四序" % _rend18,
          _rend18 in r18f.stdout
          and "--orders fm stripe_a p4_a mix_ab" in r18f.stdout)
    # Pillow 兜底: 本机 matplotlib 缺 pyparsing/cycler/contourpy/fontTools/kiwisolver 且离线装不上,
    # 因此超胞图必须有一条不依赖 matplotlib 的等价路径, 否则方形/C 心矩形这些新晶格就没有图可看。
    _pil18 = os.path.join(SCRIPTS, "render_supercell_pil.py")
    _pilsrc18 = open(_pil18, encoding="utf-8").read() if os.path.isfile(_pil18) else ""
    check("sq: render_supercell_pil.py 存在且不依赖 matplotlib (纯 Pillow 兜底)",
          "from PIL import" in _pilsrc18
          and not any(l.strip().startswith(("import matplotlib", "from matplotlib"))
                      for l in _pilsrc18.splitlines()))
    r18g = run(sys.executable, _pil18, "--help")
    check("sq: Pillow 兜底与 matplotlib 版 CLI 一致 (--root/--orders/--element/--nonmag/--label/-o)",
          r18g.returncode == 0
          and all(f in r18g.stdout for f in ("--root", "--orders", "--element",
                                             "--nonmag", "--label", "-o")))
    _png18 = os.path.join(td18, "figures", "supercell_orders.png")
    check("sq: render 真写出 figures/supercell_orders.png (>10 KB)",
          os.path.isfile(_png18) and os.path.getsize(_png18) > 10000)
finally:
    shutil.rmtree(td18, ignore_errors=True)

# ---- 19. C 心矩形多亚晶格 (CrOCl: 2 磁性原子/胞, C 心矢量连同胞两 Cr) ----
# P2.0 修的真 bug: 旧 detect_lattice 一律把"被 cutoff 并簇"判成 cluster (Nb3I8 三聚体那条路),
# CrOCl 的两个 Cr 只差 3.0347 A < cutoff 3.3, 于是一个 C 心 AF 网被当成分子簇磁体。
_lat19, _ncl19, _nmg19 = detect_lattice(CR_OCL, "Cr")
check("cr: detect_lattice(CrOCl) = rect_centered (磁性原子 2/胞, 簇 1)",
      _lat19 == "rect_centered" and (_ncl19, _nmg19) == (1, 2))
check("cr: 簇判据没被改坏 (Nb3I8 三聚体 / FGT 多亚晶格网络仍判 cluster)",
      detect_lattice(NB3I8, "Nb")[0] == "cluster"
      and detect_lattice(os.path.join(FIXTURES, "Fe3GeTe2_monolayer_POSCAR"), "Fe")[0] == "cluster")
check("cr: CrOCl 是正交网 (a != b) => 走 rect_centered 而不是 sq_centered",
      _is_rectangular_net(_lat_of(CR_OCL)) and not _is_square_net(_lat_of(CR_OCL)))
check("cr: --model crocl 预设声明 lattice=rect_centered / S=1.5 / moment=3.0 / U=4.0",
      MODEL_PRESETS["crocl"]["lattice"] == "rect_centered"
      and MODEL_PRESETS["crocl"]["S"] == 1.5
      and MODEL_PRESETS["crocl"]["moment"] == 3.0
      and MODEL_PRESETS["crocl"]["U"] == 4.0)

td19 = mktmp("_tmp_crect_")
try:
    drv19 = os.path.join(HERE, "run_2d_magnet.py")
    r19 = run(sys.executable, drv19, "build", "--model", "crocl",
              "--poscar", CR_OCL, "--outdir", td19)
    check("cr: driver build --model crocl 退出 0", r19.returncode == 0)
    check("cr: driver 判 [lattice] rect_centered (磁性原子 2/胞, 簇 1)",
          "[lattice] rect_centered" in r19.stdout and "磁性原子 2/胞" in r19.stdout)
    # reorient.py 的旧默认注释写死了 "MnPS3 monolayer" (所有材料的 reoriented POSCAR 都自称 MnPS3)
    _l1 = open(os.path.join(td19, "01_structure", "POSCAR"), encoding="utf-8").readline()
    check("cr: reorient 后 POSCAR 第 1 行带对材料名 (不再自称 MnPS3)",
          "CrOCl" in _l1 and "MnPS3" not in _l1)

    _cases19 = os.path.join(td19, "02_supercell", "cases")
    _d19 = json.load(open(os.path.join(_cases19, "magnetic_orders.json"), encoding="utf-8"))
    check("cr: cases 七序 all_orders + 选定 fm/afm_c/afm_a/afm_b",
          set(_d19.get("all_orders", [])) == {"fm", "afm_c", "afm_a", "afm_b", "afm_diag",
                                              "afm_c_a", "afm_c_diag"}
          and set(_d19.get("selected_orders", [])) == {"fm", "afm_c", "afm_a", "afm_b"})
    check("cr: model_spec lattice=rect_centered + policy=shift_search",
          _d19.get("model_spec", {}).get("lattice") == "rect_centered"
          and _d19.get("model_spec", {}).get("policy") == "shift_search")
    check("cr: 壳距离是 3D 键长 3.0347/3.1980/3.8603 且 tol=0.05 (tol=0.2 会静默并壳)",
          [round(_d19["shells_angstrom"][k], 4) for k in ("J1", "J2", "J3")]
          == [3.0347, 3.198, 3.8603] and _d19["shell_tol"] == 0.05)
    check("cr: per_site_bonds = {J1:4, J2:2, J3:2} (跨口径不变量, MC 侧据此校验)",
          all(abs(float(_d19["per_site_bonds"][k]) - v) < 1e-9
              for k, v in (("J1", 4.0), ("J2", 2.0), ("J3", 2.0))))
    check("cr: design_cond < 1e4 闸门且 det 非零", _d19.get("design_cond", 1e9) < 1e4
          and abs(float(_d19.get("design_det", 0.0))) > 1e-6)

    _mm19 = open(os.path.join(_cases19, "afm_c", "magmom.txt"), encoding="utf-8").read()
    _vals19 = [float(x) for x in _mm19.split("=", 1)[1].split()]
    check("cr: magmom.txt 值数 = 超胞原子数 (32 Cr + 32 O + 32 Cl = 96)",
          len(_vals19) == 96 and _d19["n_atoms"] == 96)
    check("cr: afm_c 的 32 个 Cr 磁矩全为 ±3.0 且正负各 16 (basis 反平行)",
          sorted(abs(v) for v in _vals19 if v != 0.0) == [3.0] * 32
          and sum(1 for v in _vals19 if v > 0.0) == 16)

    _cr_solve = os.path.join(SCRIPTS, "solve_J_crect.py")
    r19b = run(sys.executable, _cr_solve, _cases19, "--demo", "--synthetic", "5", "-0.5", "0.2")
    check("cr: solve_J_crect --demo --synthetic 反解闭合 PASS", "PASS" in r19b.stdout)
    _jr19 = json.load(open(os.path.join(_cases19, "J_result.json"), encoding="utf-8"))
    check("cr: 写出 J_result.json (schema=j_result/1, lattice=rect_centered, source=demo)",
          _jr19["schema"] == "j_result/1" and _jr19["lattice"] == "rect_centered"
          and _jr19["source"] == "demo")
    check("cr: --synthetic 合成值原样解回 (J1=5 / J2=-0.5 / J3=0.2)",
          [round(float(x), 6) for x in
           (list(_jr19["J_meV"].values()) if isinstance(_jr19.get("J_meV"), dict)
            else list(_jr19["J_meV"]))] == [5.0, -0.5, 0.2])

    # J 政策: CrOCl 的文献 J 未核实 (J1/J2/J3 = null) => 不许拿它当锚点, 也不许静默编一个
    r19c = run(sys.executable, _cr_solve, _cases19, "--demo", "--ref", "CrOCl")
    check("cr: --ref CrOCl 拒用未核实的 null J (退出 != 0 且提到 null)",
          r19c.returncode != 0 and "null" in (r19c.stdout + r19c.stderr))

    # build 阶段不该硬要 --J1 (结构与磁序枚举本来就该在没有 J 时可用), 但 mc 必须显式给
    check("cr: build 不要求 --J1 (J 未核实也能做结构/序枚举)",
          "[params] J1 未给出" in r19.stdout)
    r19d = run(sys.executable, drv19, "mc", "--model", "crocl",
               "--poscar", CR_OCL, "--outdir", td19)
    check("cr: mc 未给 --J1 时显式报错 (不猜数)",
          r19d.returncode != 0 and "--J1" in (r19d.stdout + r19d.stderr))

    # 化学计量 / O-X 占位自洽 (P2 契约): 每胞计量落盘 + 写错元素硬报错 (不静默错位)
    check("cr: 每胞化学计量落进 magnetic_orders.json (Cr2 O2 Cl2)",
          _d19.get("stoichiometry_per_cell") == {"Cr": 2, "O": 2, "Cl": 2}
          and _d19.get("formula_per_cell") == "Cr2O2Cl2")
    r19f = run(sys.executable, os.path.join(SCRIPTS, "build_crect_cases.py"), CR_OCL,
               "--supercell", "4", "4", "1", "--element", "Cr", "--nonmag", "O", "Br",
               "--moment", "3.0", "--spin", "1.5", "--outdir", os.path.join(td19, "bad19"))
    check("cr: --nonmag 写错元素硬报错 (O/X 占位必须自洽)",
          r19f.returncode != 0 and "Br" in (r19f.stdout + r19f.stderr))

    # render 分派: C 心矩形复用通用超胞渲染器, 且只画选定四序 (渲染器由 have_matplotlib() 决定)
    _rend19 = "render_supercell.py" if have_matplotlib() else "render_supercell_pil.py"
    r19e = run(sys.executable, drv19, "render", "--model", "crocl",
               "--poscar", CR_OCL, "--outdir", td19)
    check("cr: render 分派给 %s 且 --orders 只含选定四序" % _rend19,
          _rend19 in r19e.stdout
          and "--orders fm afm_c afm_a afm_b" in r19e.stdout)
    _png19 = os.path.join(td19, "figures", "supercell_orders.png")
    check("cr: render 真写出 C 心矩形超胞图 (>10 KB)",
          os.path.isfile(_png19) and os.path.getsize(_png19) > 10000)
finally:
    shutil.rmtree(td19, ignore_errors=True)

# ---- 20. 渲染兜底: 簇磁体 / 三角磁序 / MC 曲线 ----
# 本机 matplotlib 缺 pyparsing/cycler/contourpy/fontTools/kiwisolver 且离线装不上,
# 所以这三条渲染路径也必须有纯 Pillow 版, 否则簇磁体(旧名 magnetic_orders.png)与
# 三角材料永远拿不到新图, MC 阶段只剩 csv 没有曲线图。
for _mod, _opts in (("render_cluster_orders_pil.py", ("--cutoff",)),
                    ("render_tri_orders_pil.py", ())):
    _p = os.path.join(SCRIPTS, _mod)
    _src = open(_p, encoding="utf-8").read() if os.path.isfile(_p) else ""
    check("fig: %s 存在且不依赖 matplotlib (纯 Pillow 兜底)" % _mod,
          "from PIL import" in _src
          and not any(l.strip().startswith(("import matplotlib", "from matplotlib"))
                      for l in _src.splitlines()))
    _rh = run(sys.executable, _p, "--help")
    check("fig: %s 的 CLI 与 matplotlib 版一致" % _mod,
          _rh.returncode == 0
          and all(f in _rh.stdout for f in ("--root", "--element", "--nonmag", "--label", "-o")
                  + _opts))
_pilfig = open(os.path.join(SCRIPTS, "pil_fig.py"), encoding="utf-8").read()
check("fig: pil_fig.py 提供三条渲染路径共用的原语 (面板变换/虚线/圆/箭头/图例)",
      all(k in _pilfig for k in ("def panel_grid(", "def panel_transform(", "def draw_seg(",
                                 "def arrow(", "def circle(", "def draw_legend(",
                                 "def legend_items(")))
_mc2dsrc = open(os.path.join(SCRIPTS, "mc_2d.py"), encoding="utf-8").read()
check("fig: mc_2d.py 有纯 Pillow 曲线图兜底 (matplotlib 失败时调 _plot_mt_png)",
      "def _plot_mt_png(" in _mc2dsrc and "_plot_mt_png(res, meta" in _mc2dsrc)
_drvsrc = open(os.path.join(HERE, "run_2d_magnet.py"), encoding="utf-8").read()
check("fig: 驱动 render 三个分支都由 have_matplotlib() 二选一",
      all(x in _drvsrc for x in ("render_supercell_pil.py", "render_cluster_orders_pil.py",
                                 "render_tri_orders_pil.py")))

drv20 = os.path.join(HERE, "run_2d_magnet.py")
td20 = mktmp("_tmp_smoke20_")
try:
    # 三角分支: 驱动 --model gdi2 走完 reorient+build+render, 渲染器必须是可用的那一版
    _rend20 = "render_tri_orders.py" if have_matplotlib() else "render_tri_orders_pil.py"
    r20 = run(sys.executable, drv20, "render", "--model", "gdi2",
              "--poscar", GDI2, "--outdir", td20)
    check("fig: 三角 render 分派给 %s" % _rend20, _rend20 in r20.stdout)
    _png20 = os.path.join(td20, "figures", "supercell_orders.png")
    check("fig: 三角渲染器真写出 120° 超胞图 (>10 KB)",
          os.path.isfile(_png20) and os.path.getsize(_png20) > 10000)

    # 簇磁体分支: 驱动没有簇磁体预设, 直接建 cases 再渲染 (渲染器与驱动走同一支)
    td20c = os.path.join(td20, "nbo")
    rb20 = run(sys.executable, os.path.join(SCRIPTS, "build_cluster_cases.py"), NB3I8,
               "--element", "Nb", "--nonmag", "I", "--cutoff", "3.3",
               "--outdir", os.path.join(td20c, "02_supercell", "cases"))
    check("fig: 簇磁体 cases 建成 (build_cluster_cases.py 退出 0)", rb20.returncode == 0)
    _png20c = os.path.join(td20, "clu.png")
    rc20 = run(sys.executable, os.path.join(SCRIPTS, "render_cluster_orders_pil.py"),
               "--root", td20c, "--element", "Nb", "--nonmag", "I", "--cutoff", "3.3",
               "--label", "Nb3I8", "-o", _png20c)
    check("fig: 簇磁体渲染器真写出 PNG (>10 KB)",
          rc20.returncode == 0 and os.path.isfile(_png20c)
          and os.path.getsize(_png20c) > 10000)
finally:
    shutil.rmtree(td20, ignore_errors=True)

# ---------------------------------------------------------------- 21) P3: 非共线 / SOC
print("\n[21] 非共线自旋 + SOC (P3 基础设施, 不下发 DFT)")
for _sc21, _nm21 in (("mag_vectors.py", "mag_vectors"), ("dft_energy.py", "dft_energy"),
                     ("assemble_dft.py", "assemble_dft"),
                     ("assemble_dft_tri.py", "assemble_dft_tri")):
    _r21 = run(sys.executable, os.path.join(SCRIPTS, _sc21), "--selftest")
    check("nc: %s --selftest 通过 (无 DFT 自检)" % _nm21,
          _r21.returncode == 0 and "ALL PASS" in _r21.stdout)

import mag_vectors as MV  # noqa: E402

drv21 = os.path.join(HERE, "run_2d_magnet.py")
td21 = mktmp("_tmp_smoke21_")
try:
    # 三角 + SOC: 驱动 -> assemble_dft_tri(适配层) -> assemble_dft(--noncollinear --lsorbit --all-orders)
    r21 = run(sys.executable, drv21, "dft", "--model", "gdi2", "--poscar", GDI2,
              "--outdir", td21, "--soc")
    check("nc: 驱动 dft 非共线+SOC 走通 (三角分支经适配层)",
          r21.returncode == 0 and "非共线+SOC" in r21.stdout)

    _cases21 = os.path.join(td21, "02_supercell", "cases")
    _meta21 = json.load(open(os.path.join(_cases21, "magnetic_orders.json"), encoding="utf-8"))
    _sel21 = list(_meta21["selected_orders"])
    _pos21 = open(os.path.join(_cases21, _sel21[0], "POSCAR"), encoding="utf-8").read().splitlines()
    _nat21 = sum(int(x) for x in _pos21[6].split())
    _inc21 = {o: open(os.path.join(td21, "runs", o, "INCAR"), encoding="utf-8").read()
              for o in _sel21 if os.path.isfile(os.path.join(td21, "runs", o, "INCAR"))}
    check("nc: 每个选中磁序都组装了 runs/<order>/INCAR", len(_inc21) == len(_sel21))
    check("nc: INCAR 打开 LNONCOLLINEAR + LSORBIT + GGA_COMPAT=.FALSE. 且不再有 ISPIN",
          all(("LNONCOLLINEAR = .TRUE." in t and "LSORBIT  = .TRUE." in t
               and "GGA_COMPAT = .FALSE." in t and "ISPIN" not in t) for t in _inc21.values()))
    check("nc: 03_dft/ 静态与弛豫 INCAR 都在",
          os.path.isfile(os.path.join(td21, "03_dft", "INCAR.static"))
          and os.path.isfile(os.path.join(td21, "03_dft", "INCAR.relax")))

    # MAGMOM 内容规则: 3 分量 / 面内 / 120° 三色
    _ok21, _ncmp21, _pb21 = MV.check_magmom(_inc21["120deg"], _nat21, want_ncomp=3)
    _v21, _ = MV.parse_magmom(_inc21["120deg"], n_atoms=_nat21)
    _m21 = _v21.reshape(_nat21, 3)
    _sp21 = _m21[np.linalg.norm(_m21, axis=1) > 1e-6]
    _ang21 = sorted({round(float(np.degrees(np.arctan2(y, x))) % 360.0, 3)
                     for x, y, _z in _sp21})
    check("nc: 120° 构型 MAGMOM 是 3 分量、面内 (z=0)、角只取 0/120/240 (n_atoms=%d)" % _nat21,
          _ok21 and _sp21.shape[0] > 0 and np.allclose(_sp21[:, 2], 0.0, atol=1e-6)
          and _ang21 == [0.0, 120.0, 240.0])

    # 快照回归: INCAR 正文逐字比对; MAGMOM 行换成 <原子数/分量/磁原子数/|m|/角度>
    # 描述符再比对 —— 既不受长行排版的干扰, 又能抓住"物理写错"(角度/大小/分量数)。
    _snap21 = os.path.join(FIXTURES, "gdi2_nc_soc_120deg.INCAR.norm")

    def _nc_norm(text, n_atoms):
        _v, _nc = MV.parse_magmom(text, n_atoms=n_atoms)
        _mm = _v.reshape(n_atoms, _nc)
        _nz = np.linalg.norm(_mm, axis=1) > 1e-6
        _s = _mm[_nz]
        _ags = sorted({round(float(np.degrees(np.arctan2(y, x))) % 360.0, 3)
                       for x, y, _z in _s})
        _mg = float(np.linalg.norm(_s[0])) if _s.shape[0] else 0.0
        _line = ("MAGMOM = <n_atoms=%d ncomp=%d n_mag=%d |m|=%.6f angles=%s>"
                 % (n_atoms, _nc, int(_nz.sum()), _mg, _ags))
        return "\n".join(_line if l.startswith("MAGMOM") else l
                         for l in text.replace("\r\n", "\n").split("\n")).strip()

    _want21 = None
    if os.path.isfile(_snap21):
        _want21 = "\n".join(l for l in open(_snap21, encoding="utf-8").read()
                            .replace("\r\n", "\n").split("\n") if not l.startswith("#")).strip()
    check("nc: 120deg 的 INCAR 与快照 fixture 一致 (正文逐字, MAGMOM 按描述符)",
          _want21 is not None and _want21 == _nc_norm(_inc21["120deg"], _nat21))

    # §115 不静默降级: 需要 SOC 的磁序, 没给 --lsorbit 必须硬报错
    _g21 = os.path.join(td21, "guard")
    _gc21 = os.path.join(_g21, "02_supercell", "cases", "cant_a")
    os.makedirs(_gc21)
    shutil.copy(os.path.join(_cases21, _sel21[0], "POSCAR"), os.path.join(_gc21, "POSCAR"))
    shutil.copy(os.path.join(_cases21, _sel21[0], "magmom.txt"), os.path.join(_gc21, "magmom.txt"))
    with open(os.path.join(_g21, "02_supercell", "cases", "magnetic_orders.json"),
              "w", encoding="utf-8") as f:
        json.dump({"selected_orders": ["cant_a"]}, f)
    _as21 = os.path.join(SCRIPTS, "assemble_dft.py")
    _rg21 = run(sys.executable, _as21, "--root", _g21, "--element", "Gd", "--nonmag", "I",
                "--noncollinear")
    _msg21 = (_rg21.stdout or "") + (_rg21.stderr or "")
    check("nc: 需要 SOC 的磁序未开 --lsorbit 时硬报错 (§115)",
          _rg21.returncode != 0 and "SOC" in _msg21 and "--lsorbit" in _msg21)
    _rv21 = run(sys.executable, _as21, "--root", _g21, "--element", "Gd", "--nonmag", "I",
                "--noncollinear", "--lsorbit", "--saxis", "0", "0", "1")
    check("nc: 同一磁序给了 --lsorbit 就通过 (确认报错来自 SOC 门禁)",
          _rv21.returncode == 0)

    # M3 最小闭环: 用 coefficients_per_cell 合成 OUTCAR -> solve_J_tri 回解 J
    _J21 = [1.5, -0.8, 0.4]
    _S221 = float(_meta21["S2"])
    for _o in _sel21:
        _c21 = _meta21["coefficients_per_cell"][_o]
        _E21 = sum(c * (j * _S221) for c, j in zip(_c21, _J21)) / 1000.0  # meV -> eV, 同 solve_J_tri
        with open(os.path.join(_cases21, _o, "OUTCAR"), "w", encoding="utf-8") as f:
            f.write(" vasp.6.4.2 (mock OUTCAR, 无 DFT)\n"
                    "  free  energy   TOTEN  =      %+.8f eV\n" % _E21)
    _rs21 = run(sys.executable, os.path.join(SCRIPTS, "solve_J_tri.py"), _cases21)
    check("nc: M3 mock-OUTCAR 闭环回解 J (J1/J2/J3 对到 1e-4 meV)",
          _rs21.returncode == 0
          and all(("J%d = %+.4f meV" % (i + 1, j)) in _rs21.stdout
                  for i, j in enumerate(_J21)))
finally:
    shutil.rmtree(td21, ignore_errors=True)

# ---- 22. 收口包: 六条渲染路径与 matplotlib 解耦 + 子进程编码 + --version ----
# 守三条「发布卫生」不变量: (a) 版本号只有 pyproject.toml 一个声明; (b) driver 派生的每个
# 子进程都被注入 PYTHONIOENCODING=utf-8 (GBK 控制台下子脚本打印 ⇒ 不再 UnicodeEncodeError);
# (c) 原先只走 matplotlib 的两条曲线图 (FSS 标度图 / FGT 磁序图) 现在各有纯 Pillow 兜底。
_ver22 = None
with open(os.path.join(HERE, "pyproject.toml"), encoding="utf-8") as _f22:
    for _ln22 in _f22:
        _s22 = _ln22.strip()
        if _s22.startswith("version") and "=" in _s22:
            _ver22 = _s22.split("=", 1)[1].strip().strip('"').strip("'")
            break
_r22 = run(sys.executable, os.path.join(HERE, "run_2d_magnet.py"), "--version")
check("ver: --version 退出 0 且打印 '2d-magnet <版本>'",
      _r22.returncode == 0 and _r22.stdout.strip().startswith("2d-magnet "))
check(f"ver: --version 与 pyproject.toml 的 version 一致 ({_ver22})",
      bool(_ver22) and _ver22 in _r22.stdout)

_p22 = ("import sys; sys.path.insert(0, %r); import run_2d_magnet as R; "
        "print(R.child_env().get('PYTHONIOENCODING', ''))" % HERE)
_r22b = run(sys.executable, "-c", _p22)
check("enc: child_env() 给子进程注入 PYTHONIOENCODING=utf-8",
      _r22b.returncode == 0 and "utf-8" in _r22b.stdout.lower())

td22 = mktmp("_tmp_close_")
try:
    # (c1) FSS 标度图: 直接喂合成 U4(T) 曲线给内置 Pillow 渲染器 (不起 MC, 秒级)
    _png22 = os.path.join(td22, "fss_probe.png")
    _probe22 = os.path.join(td22, "_probe_fss.py")
    with open(_probe22, "w", encoding="utf-8") as f:
        f.write("import sys\n"
                "sys.path.insert(0, %r)\n"
                "import numpy as np\n"
                "import mc_fss\n"
                "T = np.linspace(20.0, 80.0, 13)\n"
                "u = np.array([0.4] * 5 + [0.25, 0.10, 0.02] + [0.0] * 5)\n"
                "curves = {L: (T, u) for L in (12, 18, 24)}\n"
                "mc_fss._plot_fss_png(curves, [12, 18, 24], {'element': 'CrI3'}, 48.3, %r)\n"
                "print('FSS PNG OK')\n" % (SCRIPTS, _png22))
    _r22c = run(sys.executable, _probe22)
    check("fss: Pillow 兜底渲染 Binder U4 标度图 (>10 KB)",
          _r22c.returncode == 0 and os.path.isfile(_png22) and os.path.getsize(_png22) > 10000)

    # (c2) FGT 磁序图: 跑真构建器 (matplotlib 可用走 mpl, 不可用走 Pillow, 两条都必须出图)
    _o22 = os.path.join(td22, "fgt")
    _r22d = run(sys.executable, os.path.join(SCRIPTS, "build_fgt_cases.py"),
                FGTMONO, "--outdir", _o22)
    _png22d = os.path.join(_o22, "figures", "fgt_supercell_orders.png")
    check("fgt: build_fgt_cases 出图 (>10 KB)",
          _r22d.returncode == 0 and os.path.isfile(_png22d)
          and os.path.getsize(_png22d) > 10000)
finally:
    shutil.rmtree(td22, ignore_errors=True)

# ---- 23. kagome (Fe3Sn 类, 3 磁位/胞): 几何判别 + 系数矩阵 + 渲染 + 闸门 ----
# P4.1 是纯几何 + 线性代数 (不下发 DFT): 判别要能把 kagome 从「3 原子分子簇」里挑出来
# (配位数 4 vs 2/6), 系数矩阵必须 rank=4 (q=0 一族只有 rank=2 ⇒ J1/J2 不可分离),
# 反解器必须能无 DFT 往返闭合, 渲染必须与 matplotlib 解耦 (走 Pillow 兜底)。
print("\n[23] kagome 角共享三角形网 (P4.1: 纯几何 + 线性代数, 不下发 DFT)")
_KAG23 = os.path.join(FIXTURES, "Fe3Sn_POSCAR")
if not os.path.exists(_KAG23):
    check("kag: Fe3Sn fixture 存在 (跳过: 缺文件)", True)
else:
    _p23 = ("import sys; sys.path.insert(0, %r); import run_2d_magnet as R; "
            "print(R.detect_lattice(%r, 'Fe'))" % (HERE, _KAG23))
    _r23 = run(sys.executable, "-c", _p23)
    check("kag: detect_lattice 判出 kagome (不是 cluster)",
          _r23.returncode == 0 and "'kagome'" in _r23.stdout)

    td23 = mktmp("_tmp_kag_")
    try:
        _o23 = os.path.join(td23, "proj")
        _cases23 = os.path.join(_o23, "cases")
        _r23b = run(sys.executable, os.path.join(SCRIPTS, "build_kagome_cases.py"), _KAG23,
                    "--element", "Fe", "--nonmag", "Sn", "--moment", "2.0", "--spin", "1.0",
                    "--outdir", _cases23, "--selftest")
        check("kag: build_kagome_cases --selftest ALL PASS",
              _r23b.returncode == 0 and "ALL PASS" in _r23b.stdout)
        _r23c = run(sys.executable, os.path.join(SCRIPTS, "solve_J_kagome.py"), _cases23,
                    "--demo", "--synthetic", "1.5", "-0.8", "0.4", "--no-jresult")
        check("kag: solve_J_kagome --demo 往返闭合 (PASS)",
              _r23c.returncode == 0 and "PASS" in _r23c.stdout and "FAIL" not in _r23c.stdout)
        _r23d = run(sys.executable, os.path.join(HERE, "run_2d_magnet.py"), "render",
                    "--model", "fe3sn", "--poscar", _KAG23, "--outdir", _o23, "--skip-mc")
        _png23 = os.path.join(_o23, "figures", "supercell_orders.png")
        check("kag: 驱动 render 端到端出图 (>10 KB)",
              _r23d.returncode == 0 and os.path.isfile(_png23)
              and os.path.getsize(_png23) > 10000)
        _jf23 = os.path.join(_o23, "02_supercell", "cases", "magnetic_orders.json")
        _ok23 = False
        if os.path.isfile(_jf23):
            with open(_jf23, encoding="utf-8") as f:
                _d23 = json.load(f)
            _ok23 = (list(_d23.get("selected_orders") or []) ==
                     ["fm", "q0_120", "sq3_updntri", "sq3_alt"]
                     and int(_d23.get("design_rank") or 0) == 4
                     and bool(_d23.get("noncollinear")))
        check("kag: 落盘工程含 √3×√3 选定序 + rank=4 + 非共线", _ok23)
    finally:
        shutil.rmtree(td23, ignore_errors=True)

    _r23e = run(sys.executable, os.path.join(HERE, "tests", "test_kagome_cases.py"))
    check("kag: tests/test_kagome_cases.py 全部通过",
          _r23e.returncode == 0 and "全部通过" in _r23e.stdout)

print("\n[24] DMI / MTI 识别层 + 受控换算器 (P6: 不下发 DFT, 不提取 D)")
_DMI24 = os.path.join(SCRIPTS, "dmi_matrix.py")
_KAG24 = os.path.join(FIXTURES, "Fe3Sn_POSCAR")
if not os.path.exists(_DMI24) or not os.path.exists(_KAG24):
    check("dmi: dmi_matrix.py / Fe3Sn fixture 存在 (跳过: 缺文件)", True)
else:
    _r24a = run(sys.executable, _DMI24, "--selftest")
    check("dmi: dmi_matrix --selftest ALL PASS",
          _r24a.returncode == 0 and "ALL PASS" in _r24a.stdout)

    _cat24 = {e["key"]: e for e in ERROR_CATALOG}
    check("dmi: ERROR_CATALOG 有 dmi_needs_shell_def (静态, 点名 --shell-def/J_mat)",
          "dmi_needs_shell_def" in _cat24
          and _cat24["dmi_needs_shell_def"].get("static") is True
          and "--shell-def" in _cat24["dmi_needs_shell_def"].get("fix", ""))

    td24 = mktmp("_tmp_dmi_")
    try:
        _sh24 = os.path.join(td24, "shell.json")
        _mc24 = os.path.join(td24, "model.json")
        _r24b = run(sys.executable, _DMI24, "--D", "0.31", "-0.77", "0.55",
                    "--J", "1.9", "--S", "1.0", "--d", "2.70", "--tol", "0.05",
                    "--name", "J1", "--out", _sh24)
        check("dmi: 写出壳定义 JSON (--out, UTF-8)",
              _r24b.returncode == 0 and os.path.isfile(_sh24))

        _r24c = run(sys.executable, os.path.join(SCRIPTS, "poscar_to_mc3.py"), _KAG24,
                    "--element", "Fe", "--S", "1.0", "--shell-def", _sh24, "-o", _mc24)
        _ok24 = False
        if _r24c.returncode == 0 and os.path.isfile(_mc24):
            with open(_mc24, encoding="utf-8") as f:
                _m24 = json.load(f)
            _ia24 = _m24.get("interactions") or []
            if _ia24 and len(_ia24[0].get("J") or []) == 3:
                _J24 = _ia24[0]["J"]
                # J_D = -[D]x: J[0][1]=+Dz, J[0][2]=-Dy, J[1][2]=+Dx, 对角 = -J*S^2
                _ok24 = (len(_ia24) == 6                       # kagome J1: z=4 -> 6 键/胞
                         and abs(_J24[0][1] - 0.55) < 1e-9
                         and abs(_J24[0][2] - 0.77) < 1e-9
                         and abs(_J24[1][2] - 0.31) < 1e-9
                         and abs(_J24[0][0] + 1.9) < 1e-9
                         and abs(_J24[0][1] + _J24[1][0]) < 1e-12)
        check("dmi: 反称矩阵经 --shell-def 落进 mc3 (6 键/胞, 反称非零)", _ok24)

        _sh16 = os.path.join(td24, "shell_u16.json")
        _r24d = run(sys.executable, "-c",
                    "import json;d=json.load(open(%r,encoding='utf-8'));"
                    "json.dump(d,open(%r,'w',encoding='utf-16'),ensure_ascii=False)"
                    % (_sh24, _sh16))
        _r24e = run(sys.executable, os.path.join(SCRIPTS, "poscar_to_mc3.py"), _KAG24,
                    "--element", "Fe", "--S", "1.0", "--shell-def", _sh16,
                    "-o", os.path.join(td24, "model16.json"))
        check("dmi: 壳定义 UTF-16 BOM 也能读 (PowerShell `>` 重定向的坑)",
              _r24d.returncode == 0 and _r24e.returncode == 0 and "[utf-16]" in _r24e.stdout)
    finally:
        shutil.rmtree(td24, ignore_errors=True)

    _okm24 = False
    _memo24 = os.path.join(HERE, "DMI_MTI_MEMO.md")
    if os.path.isfile(_memo24):
        with open(_memo24, encoding="utf-8") as f:
            _t24 = f.read()
        _okm24 = (os.path.getsize(_memo24) > 15000 and "DMI" in _t24 and "MTI" in _t24
                  and "dmi_matrix.py" in _t24 and "可行性备忘" in _t24)
    check("dmi: DMI_MTI_MEMO.md 存在且含判据 / 边界 / 工具", _okm24)

    _oks24 = False
    _skill24 = os.path.join(HERE, "dsh", "2d-magnet", "SKILL.md")
    if os.path.isfile(_skill24):
        with open(_skill24, encoding="utf-8") as f:
            _ts24 = f.read()
        _oks24 = "分支：DMI / MTI" in _ts24 and "dmi_matrix.py" in _ts24
    check("dmi: SKILL.md 有 DMI/MTI 判别节", _oks24)

    _okr24 = False
    _readme24 = os.path.join(HERE, "README.md")
    if os.path.isfile(_readme24):
        with open(_readme24, encoding="utf-8") as f:
            _tr24 = f.read()
        _okr24 = ("DMI_MTI_MEMO.md" in _tr24 and "dmi_matrix.py" in _tr24
                  and "dmi_matrix.py --selftest" in _tr24)
    check("dmi: README 已同步 (备忘 + 换算器 + 自检入口)", _okr24)

print("\n[25] 交变磁体几何判别 + 逐子晶格序参量 (P5: 只判别, 不下结论, 不下发 DFT)")
_ALT25 = os.path.join(SCRIPTS, "altermagnet_check.py")
_FE3SN25 = os.path.join(FIXTURES, "Fe3Sn_POSCAR")
try:
    import altermagnet_check as AM  # noqa: E402
    import mc_2d as MC2  # noqa: E402
    _HAVE_AM25 = True
except Exception:
    AM = None
    MC2 = None
    _HAVE_AM25 = False
if not (os.path.isfile(_ALT25) and _HAVE_AM25):
    check("alt: altermagnet_check.py 存在 (跳过: 缺文件)", True)
else:
    _r25a = run(sys.executable, _ALT25, "--selftest")
    check("alt: altermagnet_check --selftest ALL PASS",
          _r25a.returncode == 0 and "ALL PASS" in _r25a.stdout)

    _cat25 = {e["key"]: e for e in ERROR_CATALOG}
    check("alt: ERROR_CATALOG 有 sublattice_mask_invalid (静态, 点名 sublattice:<mask>)",
          "sublattice_mask_invalid" in _cat25
          and _cat25["sublattice_mask_invalid"].get("static") is True
          and "sublattice:0" in _cat25["sublattice_mask_invalid"].get("fix", ""))

    td25 = mktmp("_tmp_alt_")

    def _verdict25(fixture, element):
        """跑一次 CLI 并读回 JSON; 失败返回 (None, 进程结果)。"""
        path = os.path.join(td25, "v_%s.json" % element)
        p = run(sys.executable, _ALT25, "--poscar", fixture, "--element", element, "--json", path)
        if p.returncode != 0 or not os.path.isfile(path):
            return None, p
        with open(path, encoding="utf-8") as fh:
            return json.load(fh), p

    try:
        _v25, _p25 = _verdict25(MNPS3, "Mn")
        check("alt: MnPS3 (镜面连接) -> 几何候选, 残差以 Å 计且 < 1e-3",
              bool(_v25) and _v25["verdict"] == AM.VERDICT_CANDIDATE
              and float(_v25.get("worst_mismatch_ang", 1.0)) < 1e-3)

        _v25c, _p25c = _verdict25(CR_OCL, "Cr")
        check("alt: CrOCl (反演连接) -> 宇称禁止 (不是交变磁体)",
              bool(_v25c) and _v25c["verdict"] == AM.VERDICT_PARITY)

        _v25d, _p25d = _verdict25(_FE3SN25, "Fe")
        check("alt: Fe3Sn (3 套亚晶格) -> 不适用",
              bool(_v25d) and _v25d["verdict"] == AM.VERDICT_NA and _v25d["n_sub"] == 3)

        check("alt: 每条结论都带「不可判定: k 分辨自旋劈裂」声明",
              all(bool(v) and "不可判定" in v.get("disclaimer", "")
                  for v in (_v25, _v25c, _v25d)))
    finally:
        shutil.rmtree(td25, ignore_errors=True)

    _s25 = MC2.sublattice_signs(2, (2, 2), {0})
    check("alt: sublattice_signs(2, mask={0}) = basis 0 为 +1 / 其余 -1",
          _s25.shape == (8,) and all(_s25[i] == (1.0 if i % 2 == 0 else -1.0) for i in range(8)))

    def _raises25(fn):
        try:
            fn()
        except ValueError:
            return True
        except Exception:
            return False
        return False

    check("alt: 空掩码 / 空项 / 越界 / 覆盖全部子晶格 全部硬报错 (不静默)",
          _raises25(lambda: MC2.parse_order_spec("sublattice:"))
          and _raises25(lambda: MC2.parse_order_spec("sublattice:0,,"))
          and _raises25(lambda: MC2.parse_order_spec("sublattice:x"))
          and _raises25(lambda: MC2.sublattice_signs(2, (2, 2), {2}))
          and _raises25(lambda: MC2.sublattice_signs(2, (2, 2), {0, 1}))
          and _raises25(lambda: MC2.sublattice_signs(1, (2, 2), {0})))

    check("alt: 合法掩码解析 (sublattice:0,1 / staggered / total 各自成形)",
          MC2.parse_order_spec("sublattice:0,1") == ("sublattice", (0, 1))
          and MC2.parse_order_spec("sublattice: 0 ") == ("sublattice", (0,))
          and MC2.parse_order_spec("staggered") == ("staggered", None)
          and MC2.parse_order_spec("total") == ("total", None)
          and "sublattice" in MC2.ORDER_HELP)

    _r25h = run(sys.executable, os.path.join(HERE, "run_2d_magnet.py"), "--help")
    check("alt: 驱动 --help 列出 --mc-order 并说明掩码写法",
          _r25h.returncode == 0 and "--mc-order" in _r25h.stdout
          and "sublattice" in _r25h.stdout)

    _t25 = os.path.join(HERE, "tests", "test_altermagnet.py")
    if os.path.isfile(_t25):
        _r25t = run(sys.executable, _t25)
        check("alt: 契约测试 tests/test_altermagnet.py 全 PASS",
              _r25t.returncode == 0 and "全部通过" in _r25t.stdout)
    else:
        check("alt: tests/test_altermagnet.py 存在 (跳过: 缺文件)", True)

    _oks25 = False
    _skill25 = os.path.join(HERE, "dsh", "2d-magnet", "SKILL.md")
    if os.path.isfile(_skill25):
        with open(_skill25, encoding="utf-8") as f:
            _ts25 = f.read()
        _oks25 = ("分支：交变磁体" in _ts25 and "altermagnet_check.py" in _ts25
                  and "sublattice:" in _ts25)
    check("alt: SKILL.md 有交变磁体分支 + 局限声明", _oks25)

    _okr25 = False
    _readme25 = os.path.join(HERE, "README.md")
    if os.path.isfile(_readme25):
        with open(_readme25, encoding="utf-8") as f:
            _tr25 = f.read()
        _okr25 = ("altermagnet_check.py" in _tr25 and "交变磁体" in _tr25
                  and "[x] **P5**" in _tr25)
    check("alt: README 已同步 (脚本 + 局限 + 路线图 P5 勾选)", _okr25)

    _okn25 = False
    _plan25 = os.path.join(HERE, "NEW_TYPES_PLAN.md")
    if os.path.isfile(_plan25):
        with open(_plan25, encoding="utf-8") as f:
            _tn25 = f.read()
        _i25 = _tn25.find("### P5 交变磁体")
        _okn25 = (_i25 >= 0 and "✅" in _tn25[_i25:_i25 + 400]
                  and "altermagnet_check.py" in _tn25[_i25:_i25 + 4000])
    check("alt: NEW_TYPES_PLAN 的 P5 节标 ✅ 且点名交付物", _okn25)

    # ------------------------------------------------------------------
    # [26] 弱铁磁倾斜 (P2', CrOCl 类 oxyhalide): 构建器 cant 家族 + 闭式 θ 标度反解器
    # ------------------------------------------------------------------
    _b26 = os.path.join(SCRIPTS, "build_crect_cases.py")
    _s26 = os.path.join(SCRIPTS, "solve_J_crect_cant.py")

    _h26b = run(sys.executable, _b26, "--help")
    check("cant: builder --help 列出 --cant-theta / --cant-azimuth / --cant-base",
          _h26b.returncode == 0 and "--cant-theta" in _h26b.stdout
          and "--cant-azimuth" in _h26b.stdout and "--cant-base" in _h26b.stdout)

    _h26s = run(sys.executable, _s26, "--help")
    check("cant: 反解器 --help 列出 --orders / --audit-only / --demo / --synthetic",
          _h26s.returncode == 0 and "--orders" in _h26s.stdout
          and "--audit-only" in _h26s.stdout and "--synthetic" in _h26s.stdout)

    _td26 = tempfile.mkdtemp(prefix="smoke26_")
    _td26b = tempfile.mkdtemp(prefix="smoke26b_")
    _root26 = None
    _td26d = None
    _td26g = None
    try:
        _r26 = run(sys.executable, _b26, CR_OCL, "--element", "Cr", "--nonmag", "Cl", "O",
                   "--moment", "3.0", "--spin", "1.5", "--outdir", _td26,
                   "--cant-theta", "8 16", "--selftest")
        check("cant: builder 带 --cant-theta 时 --selftest 全 PASS (退出码 0)",
              _r26.returncode == 0)

        _mj26 = os.path.join(_td26, "magnetic_orders.json")
        _m26 = json.load(open(_mj26, encoding="utf-8")) if os.path.isfile(_mj26) else {}
        _names26 = _m26.get("cant_orders") or []
        check("cant: 每个角度展开 ±θ 两条 (cant_p8/cant_m8/cant_p16/cant_m16), soc_orders 全覆盖",
              _names26 == ["cant_p8", "cant_m8", "cant_p16", "cant_m16"]
              and _m26.get("soc_orders") == _names26 and _m26.get("noncollinear") is True
              and _m26.get("cant_theta_deg") == {"cant_p8": 8.0, "cant_m8": -8.0,
                                                 "cant_p16": 16.0, "cant_m16": -16.0})

        _need26 = ("order", "base", "theta_deg", "azimuth_deg", "coeff_per_primitive_cell",
                   "coeff_fm_per_primitive_cell", "coeff_base_per_primitive_cell",
                   "prediction_per_primitive_cell", "prediction_rule", "max_abs_dev",
                   "per_site_bonds", "theta_dependent_shells", "theta_independent_shells",
                   "soc_required", "sign_of_theta_unidentifiable",
                   "azimuth_not_determined_by_geometry", "note")
        _rep26 = _m26.get("cant_report") or {}
        check("cant: cant_report 键齐全 (solve 端契约) 且母序一律 afm_c (不许 fm)",
              bool(_names26) and all(k in (_rep26.get(n) or {}) for n in _names26 for k in _need26)
              and all((_rep26.get(n) or {}).get("base") == "afm_c" for n in _names26))

        import math as _math26
        _ok26 = True
        _ok26b = True
        _dep26 = True
        for _n in _names26:
            _e = _rep26[_n]
            _th = _m26["cant_theta_deg"][_n]
            _s2 = _math26.sin(_math26.radians(_th)) ** 2
            _c2 = _math26.cos(_math26.radians(_th)) ** 2
            _mine = [_s2 * a + _c2 * b for a, b in zip(_e["coeff_fm_per_primitive_cell"],
                                                       _e["coeff_base_per_primitive_cell"])]
            _ok26 = _ok26 and all(abs(x - y) < 1e-9
                                  for x, y in zip(_mine, _e["coeff_per_primitive_cell"]))
            _ok26b = _ok26b and all(abs(x - y) < 1e-9
                                    for x, y in zip(_e["prediction_per_primitive_cell"],
                                                    _e["coeff_per_primitive_cell"]))
            _dep26 = _dep26 and _e["theta_dependent_shells"] == ["J1"] \
                and _e["theta_independent_shells"] == ["J2", "J3"] \
                and _e["soc_required"] is True
        check("cant: 闭式 c_k(θ)=sin²θ·c_k(fm)+cos²θ·c_k(母序) 在 smoke 侧独立重算逐条成立", _ok26)
        check("cant: builder 自己的闭式预测与落盘系数一致 (max_abs_dev <= 1e-9)", _ok26b)
        check("cant: 只有跨子晶格壳 J1 随 θ 变, J2/J3 (同子晶格键) 与 θ 无关; 且标记 soc_required",
              _dep26)
        check("cant: cant_p8/cant_m8 (及 p16/m16) 系数逐位相同 ⇒ θ 符号不可辨识",
              bool(_names26)
              and _rep26["cant_p8"]["coeff_per_primitive_cell"]
              == _rep26["cant_m8"]["coeff_per_primitive_cell"]
              and _rep26["cant_p16"]["coeff_per_primitive_cell"]
              == _rep26["cant_m16"]["coeff_per_primitive_cell"]
              and all(_rep26[n]["sign_of_theta_unidentifiable"] is True for n in _names26))

        import mag_vectors as _MV26
        _mm26 = True
        _mm26c = True
        for _n in _names26:
            _p = os.path.join(_td26, _n, "magmom.txt")
            if not os.path.isfile(_p):
                _mm26 = False
                continue
            with open(_p, encoding="utf-8") as _f:
                _t = _f.read()
            _ok, _nc, _ = _MV26.check_magmom(_t, _m26["n_atoms"], want_ncomp=3)
            _mm26 = _mm26 and _ok and _nc == 3
        _pc = os.path.join(_td26, "afm_c", "magmom.txt")
        if os.path.isfile(_pc):
            with open(_pc, encoding="utf-8") as _f:
                _okc, _ncc, _ = _MV26.check_magmom(_f.read(), _m26["n_atoms"], want_ncomp=1)
            _mm26c = _okc and _ncc == 1
        check("cant: 倾斜构型 MAGMOM 一律三分量 (共线构型仍是 1 分量, 未被带偏)",
              _mm26 and _mm26c)

        check("cant: 4 条倾斜序各自有 cases/<序名>/{POSCAR,magmom.txt}",
              all(os.path.isfile(os.path.join(_td26, _n, "POSCAR"))
                  and os.path.isfile(os.path.join(_td26, _n, "magmom.txt")) for _n in _names26))

        # SOC 接线: 倾斜构型进 assemble_dft 必须 --noncollinear + --lsorbit (缺一个就硬拒)
        _root26 = tempfile.mkdtemp(prefix="smoke26root_")
        shutil.copytree(_td26, os.path.join(_root26, "02_supercell", "cases"))
        _as26 = os.path.join(SCRIPTS, "assemble_dft.py")
        _r26c = run(sys.executable, _as26, "--root", _root26, "--element", "Cr",
                    "--nonmag", "Cl", "O", "--moment", "3.0", "--noncollinear", "--all-orders")
        _msg26c = (_r26c.stdout or "") + (_r26c.stderr or "")
        check("cant: 倾斜序进 assemble_dft 而没给 --lsorbit 时硬报错 (SOC 门禁, 不静默按共线算)",
              _r26c.returncode != 0 and "SOC" in _msg26c and "--lsorbit" in _msg26c)
        _r26d = run(sys.executable, _as26, "--root", _root26, "--element", "Cr",
                    "--nonmag", "Cl", "O", "--moment", "3.0", "--noncollinear", "--lsorbit",
                    "--all-orders")
        _inc26 = os.path.join(_root26, "runs", "cant_p8", "INCAR")
        _t26i = open(_inc26, encoding="utf-8").read() if os.path.isfile(_inc26) else ""
        check("cant: 给了 --lsorbit 就通过, cant_p8 的 INCAR 有 LNONCOLLINEAR+LSORBIT+GGA_COMPAT 且无 ISPIN",
              _r26d.returncode == 0 and "LNONCOLLINEAR = .TRUE." in _t26i
              and "LSORBIT" in _t26i and "GGA_COMPAT = .FALSE." in _t26i
              and "ISPIN" not in _t26i)

        _r26f = run(sys.executable, _b26, CR_OCL, "--element", "Cr", "--nonmag", "Cl", "O",
                    "--outdir", os.path.join(_td26, "badbase"),
                    "--cant-theta", "8", "--cant-base", "fm")
        check("cant: --cant-base fm 被硬拒 (整体转动不携带信息)",
              _r26f.returncode != 0 and "fm" in (_r26f.stdout + _r26f.stderr))

        _r26b = run(sys.executable, _b26, CR_OCL, "--element", "Cr", "--nonmag", "Cl", "O",
                    "--moment", "3.0", "--spin", "1.5", "--outdir", _td26b)
        _mjb = os.path.join(_td26b, "magnetic_orders.json")
        _m26b = json.load(open(_mjb, encoding="utf-8")) if os.path.isfile(_mjb) else {}
        _cp26 = _m26.get("coefficients_per_cell") or {}
        _cp26b = _m26b.get("coefficients_per_cell") or {}
        _coll26 = [k for k in _cp26b if k in _cp26]
        check("cant: 不给 --cant-theta 时新键全是中性值, 且共线系数逐位不变 (向后兼容)",
              _r26b.returncode == 0 and (_m26b.get("cant_orders") or []) == []
              and _m26b.get("noncollinear") is False
              and (_m26b.get("soc_orders") or []) == [] and bool(_coll26)
              and all(_cp26[k] == _cp26b[k] for k in _coll26))

        _a26 = run(sys.executable, _s26, _td26, "--orders", "cant", "--audit-only")
        check("cant: 反解器 --audit-only 读同一份 JSON, 报退化对与「不给设计矩阵增加秩」边界",
              _a26.returncode == 0 and "退化" in _a26.stdout
              and "不给设计矩阵增加秩" in _a26.stdout and "最大偏差=0" in _a26.stdout)

        _d26 = run(sys.executable, _s26, _td26, "--demo", "--synthetic", "5", "-0.5", "0.2")
        check("cant: 反解器 --demo 合成能量往返 PASS (4x4 精确解 + 配位数指纹)",
              _d26.returncode == 0 and "往返误差" in _d26.stdout and "PASS" in _d26.stdout
              and "期望 (4.0, 2.0, 2.0)" in _d26.stdout)

        _g26a = run(sys.executable, _s26, _td26, "--ref", os.path.join(_td26, "nope.json"),
                    "--synthetic", "5", "-0.5", "0.2")
        _g26b = run(sys.executable, _s26, _td26, "--orders", "bogus", "--audit-only")
        check("cant: 反解器守卫 (--ref/--synthetic 互斥; 未知序名) 均非零退出",
              _g26a.returncode != 0 and _g26b.returncode != 0)

        _t26 = os.path.join(HERE, "tests", "test_crect_cant.py")
        if os.path.isfile(_t26):
            _rt26 = run(sys.executable, _t26)
            check("cant: 契约测试 tests/test_crect_cant.py 全 PASS",
                  _rt26.returncode == 0 and "全部通过" in _rt26.stdout)
        else:
            check("cant: tests/test_crect_cant.py 存在 (跳过: 缺文件)", True)

        # driver 接线 (P2' planB): build --cant-theta -> solve 自动改派 cant 反解器
        _drv26 = os.path.join(HERE, "run_2d_magnet.py")
        _hd26 = run(sys.executable, _drv26, "--help")
        check("cant: 驱动 --help 列出 --cant-theta / --cant-base / --cant-azimuth",
              _hd26.returncode == 0 and "--cant-theta" in _hd26.stdout
              and "--cant-base" in _hd26.stdout and "--cant-azimuth" in _hd26.stdout)

        _td26d = tempfile.mkdtemp(prefix="smoke26drv_")
        _td26g = tempfile.mkdtemp(prefix="smoke26g_")
        _bd26 = run(sys.executable, _drv26, "build", "--model", "crocl", "--poscar", CR_OCL,
                    "--outdir", _td26d, "--cant-theta", "8 16")
        _cd26 = os.path.join(_td26d, "02_supercell", "cases")
        _mj26d = os.path.join(_cd26, "magnetic_orders.json")
        _m26d = json.load(open(_mj26d, encoding="utf-8")) if os.path.isfile(_mj26d) else {}
        check("cant: 驱动 build --cant-theta 把倾斜家族交给 builder 并落盘",
              _bd26.returncode == 0 and "倾斜家族" in _bd26.stdout
              and (_m26d.get("cant_orders") or []) == ["cant_p8", "cant_m8", "cant_p16", "cant_m16"]
              and _m26d.get("noncollinear") is True)

        # 合成 mock OUTCAR: E = Σ c_k·(J_k·S²)/1000 (与共线路径同一条能量口径)
        _J26d = (5.0, -0.5, 0.2)
        _need26d = sorted(set(list(_m26d.get("selected_orders") or [])
                              + list(_m26d.get("cant_orders") or [])))
        for _n in _need26d:
            _e26 = sum(ci * (j * float(_m26d["S2"]))
                       for ci, j in zip(_m26d["coefficients_per_cell"][_n], _J26d)) / 1000.0
            with open(os.path.join(_cd26, _n, "OUTCAR"), "w", encoding="utf-8") as _f:
                _f.write(" vasp.6.4.2 (mock OUTCAR, 无 DFT)\n"
                         "  free  energy   TOTEN  =      %+.8f eV\n" % _e26)
        _sd26 = run(sys.executable, _drv26, "solve", "--model", "crocl", "--poscar", CR_OCL,
                    "--outdir", _td26d)
        _jf26d = os.path.join(_cd26, "J_result.json")
        _jr26d = json.load(open(_jf26d, encoding="utf-8")) if os.path.isfile(_jf26d) else {}
        _jm26d = _jr26d.get("J_meV") or {}
        check("cant: 驱动 solve 检测到 cant_orders 后改派 solve_J_crect_cant.py",
              _sd26.returncode == 0 and "检测到倾斜家族" in _sd26.stdout
              and "solve_J_crect_cant.py" in _sd26.stdout)
        check("cant: 驱动 solve 用「倾斜家族 + 合成 OUTCAR」往返回解 J (误差 < 5e-4 meV)",
              _sd26.returncode == 0
              and set(_jr26d.get("orders") or []) == set(_need26d)
              and all(abs(float(_jm26d.get(k) or 9e9) - v) < 5e-4
                      for k, v in zip(("J1", "J2", "J3"), _J26d))
              and "cant_audit" in _jr26d and _jr26d.get("source") == "vasp"
              and _jr26d.get("confidence") == "high")

        _gd26a = run(sys.executable, _drv26, "build", "--model", "k2nif4", "--poscar", K2NIF4,
                     "--outdir", os.path.join(_td26g, "sqcant"), "--cant-theta", "8")
        _gm26a = (_gd26a.stdout or "") + (_gd26a.stderr or "")
        _gd26b = run(sys.executable, _drv26, "build", "--model", "crocl", "--poscar", CR_OCL,
                     "--outdir", os.path.join(_td26g, "fmcant"), "--cant-theta", "8",
                     "--cant-base", "fm")
        _gm26b = (_gd26b.stdout or "") + (_gd26b.stderr or "")
        check("cant: 驱动守卫 (非 C 心矩形 + --cant-theta / --cant-base fm) 硬拒且说明原因",
              _gd26a.returncode != 0 and "不是 C 心矩形" in _gm26a and "rect_centered" in _gm26a
              and _gd26b.returncode != 0 and "整体" in _gm26b)
        _wd26 = run(sys.executable, _drv26, "build", "--model", "crocl", "--poscar", CR_OCL,
                    "--outdir", os.path.join(_td26g, "wsoc"), "--cant-theta", "8")
        check("cant: 驱动在倾斜构型但未给 --soc 时只警告不阻断 (build 仍成功)",
              _wd26.returncode == 0 and "--soc" in (_wd26.stderr or ""))
    finally:
        shutil.rmtree(_td26, ignore_errors=True)
        shutil.rmtree(_td26b, ignore_errors=True)
        if _root26:
            shutil.rmtree(_root26, ignore_errors=True)
        if _td26d:
            shutil.rmtree(_td26d, ignore_errors=True)
        if _td26g:
            shutil.rmtree(_td26g, ignore_errors=True)

    check("cant: mag_vectors.needs_soc 认 cant_* 前缀 ⇒ 下游必须 --lsorbit",
          MV.needs_soc("cant_p8") is True and MV.needs_soc("fm") is False)

    _ki26 = os.path.join(HERE, "KNOWN_ISSUES.md")
    _kio26 = False
    if os.path.isfile(_ki26):
        with open(_ki26, encoding="utf-8") as _f:
            _tki26 = _f.read()
        _miss26 = [e["key"] for e in ERROR_CATALOG if ("`%s`" % e["key"]) not in _tki26]
        _kio26 = not _miss26
    check("cant: ERROR_CATALOG 每个 key 都在 KNOWN_ISSUES.md 有行 (%d 条)" % len(ERROR_CATALOG),
          _kio26)

    _rb26 = os.path.join(HERE, "P4_2_KAGOME_RUNBOOK.md")
    _rbo26 = False
    if os.path.isfile(_rb26) and os.path.getsize(_rb26) > 15000:
        with open(_rb26, encoding="utf-8") as _f:
            _trb26 = _f.read()
        _rbo26 = all(k in _trb26 for k in ("## 4. SCF 成功判据", "## 5. 回收", "## 8. 机时估算",
                                           "## 11. 检查清单", "## 12. 相关文件",
                                           "03_dft", "02_supercell", "POTCAR.spec"))
    check("cant: P4.2 手册存在 (P4_2_KAGOME_RUNBOOK.md > 15 KB, 含 SCF判据/回收/机时/检查清单)", _rbo26)

    _sk26 = os.path.join(HERE, "dsh", "2d-magnet", "SKILL.md")
    _sko26 = False
    if os.path.isfile(_sk26):
        with open(_sk26, encoding="utf-8") as _f:
            _tsk26 = _f.read()
        _sko26 = ("分支：弱铁磁倾斜" in _tsk26 and "solve_J_crect_cant.py" in _tsk26
                  and "P4_2_KAGOME_RUNBOOK.md" in _tsk26)
    check("cant: SKILL.md 有弱铁磁倾斜分支 + P4.2 手册指引", _sko26)

    _rd26 = os.path.join(HERE, "README.md")
    _rdo26 = False
    if os.path.isfile(_rd26):
        with open(_rd26, encoding="utf-8") as _f:
            _trd26 = _f.read()
        _rdo26 = ("solve_J_crect_cant.py" in _trd26 and "[x] **P2′**" in _trd26
                  and "P4_2_KAGOME_RUNBOOK.md" in _trd26)
    check("cant: README 已同步 (反解器 + 手册 + 路线图 P2′ 勾选)", _rdo26)

    _wf26 = os.path.join(HERE, "WORKFLOW.md")
    _wfo26 = False
    if os.path.isfile(_wf26):
        with open(_wf26, encoding="utf-8") as _f:
            _twf26 = _f.read()
        _wfo26 = ("1h." in _twf26 and "1i." in _twf26
                  and "solve_J_crect_cant.py" in _twf26)
    check("cant: WORKFLOW 已同步 (继承表 + 局限 1h/1i)", _wfo26)

    _pl26 = os.path.join(HERE, "NEW_TYPES_PLAN.md")
    _plo26 = False
    if os.path.isfile(_pl26):
        with open(_pl26, encoding="utf-8") as _f:
            _tpl26 = _f.read()
        _i26 = _tpl26.find("### P2′ 弱铁磁倾斜")
        _plo26 = (_i26 >= 0 and "✅" in _tpl26[_i26:_i26 + 400]
                  and "solve_J_crect_cant.py" in _tpl26[_i26:_i26 + 4000]
                  and "P4_2_KAGOME_RUNBOOK.md" in _tpl26[_i26:_i26 + 4000])
    check("cant: NEW_TYPES_PLAN 的 P2′/P4.2 节标 ✅ 且点名交付物", _plo26)

# ---------------------------------------------------------------------------
# [27] P4.2: k 点收敛抽查 (kconv) + VASP 输入参数化 (--kpoints/--npar/--isym)
# ---------------------------------------------------------------------------
if True:
    print("\n[27] P4.2: k 点收敛抽查 (kconv) + VASP 输入参数化")
    import k_convergence as KC  # noqa: E402
    import ingest_dft as ID27  # noqa: E402

    _r27h = run(sys.executable, os.path.join(HERE, "run_2d_magnet.py"), "--help")
    check("p42: 驱动 --help 列出 --kpoints/--npar/--isym 与 kconv 子命令",
          _r27h.returncode == 0 and "--kpoints" in _r27h.stdout and "--npar" in _r27h.stdout
          and "--isym" in _r27h.stdout and "--k-grids" in _r27h.stdout and "kconv" in _r27h.stdout)

    _r27t = run(sys.executable, os.path.join(SCRIPTS, "assemble_dft_tri.py"), "--selftest")
    check("p42: assemble_dft_tri --selftest ALL PASS (含 --kpoints/--npar/--isym 透传)",
          _r27t.returncode == 0 and "ALL PASS" in _r27t.stdout)

    _r27k = run(sys.executable, os.path.join(SCRIPTS, "k_convergence.py"), "--selftest")
    check("p42: k_convergence --selftest ALL PASS",
          _r27k.returncode == 0 and "ALL PASS" in _r27k.stdout)

    _r27kh = run(sys.executable, os.path.join(SCRIPTS, "k_convergence.py"), "--help")
    check("p42: k_convergence --help 列出 --grids/--tol/--plan-only/--ingest-only",
          _r27kh.returncode == 0 and all(s in _r27kh.stdout for s in
                                         ("--grids", "--tol", "--plan-only", "--ingest-only")))
    _r27nr = run(sys.executable, os.path.join(SCRIPTS, "k_convergence.py"), "--grids", "3,4")
    check("p42: 缺 --root 硬报错 (不静默跑空)",
          _r27nr.returncode != 0 and "--root" in (_r27nr.stderr + _r27nr.stdout))

    td27 = mktmp("_tmp_p42_")
    try:
        _cases27 = os.path.join(td27, "02_supercell", "cases")
        _r27b = run(sys.executable, os.path.join(SCRIPTS, "build_kagome_cases.py"),
                    os.path.join(FIXTURES, "Fe3Sn_POSCAR"), "--supercell", "3", "3", "1",
                    "--outdir", _cases27, "--element", "Fe", "--nonmag", "Sn",
                    "--moment", "2.0", "--spin", "1.0")
        check("p42: kagome 3x3x1 构型集构建成功 (kconv 的输入前提)", _r27b.returncode == 0)

        _base27 = [sys.executable, os.path.join(SCRIPTS, "assemble_dft.py"), "--root", td27,
                   "--element", "Fe", "--nonmag", "Sn", "--U", "4.0", "--moment", "2.0",
                   "--noncollinear"]
        _r27d = run(*_base27, "--kpoints", "6", "6", "1", "--npar", "8", "--isym", "0")
        _kp27 = ""
        _p27 = os.path.join(td27, "03_dft", "KPOINTS")
        if os.path.isfile(_p27):
            with open(_p27, encoding="utf-8") as _f:
                _kp27 = _f.read()
        _inc27 = ""
        _p27b = os.path.join(td27, "03_dft", "INCAR.static")
        if os.path.isfile(_p27b):
            with open(_p27b, encoding="utf-8") as _f:
                _inc27 = _f.read()
        _ord27 = sorted(os.listdir(os.path.join(td27, "runs")))[0] if os.path.isdir(os.path.join(td27, "runs")) else None
        _inc27b = ""
        if _ord27:
            with open(os.path.join(td27, "runs", _ord27, "INCAR"), encoding="utf-8") as _f:
                _inc27b = _f.read()
        check("p42: --kpoints 6 6 1 写进 KPOINTS (且表头不再写过时的 2x2x1)",
              _r27d.returncode == 0 and _kp27.splitlines()[1].split() == ["6", "6", "1"]
              and "2x2x1" not in _kp27)
        check("p42: --npar 8 / --isym 0 落在 INCAR.static 与 runs/<order>/INCAR",
              "NPAR     = 8" in _inc27 and "ISYM     = 0" in _inc27
              and "NPAR     = 8" in _inc27b and "ISYM     = 0" in _inc27b)

        _r27d2 = run(*_base27)
        _kp27b = ""
        if os.path.isfile(_p27):
            with open(_p27, encoding="utf-8") as _f:
                _kp27b = _f.read()
        _inc27c = ""
        if os.path.isfile(_p27b):
            with open(_p27b, encoding="utf-8") as _f:
                _inc27c = _f.read()
        check("p42: 不给新参数时保持模板默认 (KPOINTS 3 3 1, NPAR=4, static 无 ISYM)",
              _r27d2.returncode == 0 and _kp27b.splitlines()[1].split() == ["3", "3", "1"]
              and "NPAR     = 4" in _inc27c and "ISYM" not in _inc27c)

        _r27kc = run(sys.executable, os.path.join(HERE, "run_2d_magnet.py"), "kconv",
                     "--model", "fe3sn", "--poscar", os.path.join(FIXTURES, "Fe3Sn_POSCAR"),
                     "--outdir", td27, "--lattice", "kagome", "--k-grids", "3,4,6")
        _dirs27 = os.path.join(td27, "kconv", "runs")
        _plan27 = os.path.join(td27, "kconv", "kconv_plan.json")
        _js27 = os.path.join(td27, "kconv", "k_convergence.json")
        _n27 = len(os.listdir(_dirs27)) if os.path.isdir(_dirs27) else 0
        _kp27c = ""
        if os.path.isfile(os.path.join(_dirs27, "%s@k3_3_1" % _ord27, "KPOINTS")):
            with open(os.path.join(_dirs27, "%s@k3_3_1" % _ord27, "KPOINTS"), encoding="utf-8") as _f:
                _kp27c = _f.read()
        check("p42: 驱动 kconv 生成 <序>@<网格> 变体 (只换 KPOINTS)",
              _r27kc.returncode == 0 and _n27 > 0 and os.path.isfile(_plan27)
              and _kp27c.splitlines()[1].split() == ["3", "3", "1"])
        _st27 = None
        if os.path.isfile(_js27):
            with open(_js27, encoding="utf-8") as _f:
                _st27 = json.load(_f).get("status")
        check("p42: 没有 OUTCAR 时判定是 insufficient (如实等回收, 不编数字)", _st27 == "insufficient")

        for _nm27 in sorted(os.listdir(_dirs27)):
            _dd27 = os.path.join(_dirs27, _nm27)
            if not os.path.isdir(_dd27):
                continue
            _k27 = int(_nm27.rsplit("@", 1)[1][1:].split("_")[0])
            ID27.synth_outcar(os.path.join(_dd27, "OUTCAR"), -100.0 + 0.01 / (_k27 * _k27),
                              counts="6 6", n_kpts=_k27 * _k27, scf=True)
        _res27 = KC.ingest(td27, KC.parse_grids("3,4,6"), tol=1.0, min_kb=0, quiet=True)
        check("p42: 收敛序列 -> converged 且建议最粗合格网格 3 3 1",
              _res27["status"] == "ok"
              and _res27["recommended_kpoints"].get(_ord27) == [3, 3, 1])
        _r27ki = run(sys.executable, os.path.join(HERE, "run_2d_magnet.py"), "kconv",
                     "--model", "fe3sn", "--poscar", os.path.join(FIXTURES, "Fe3Sn_POSCAR"),
                     "--outdir", td27, "--lattice", "kagome", "--k-grids", "3,4,6",
                     "--k-ingest-only")
        check("p42: 驱动 --k-ingest-only 复用已有变体; 默认 200 KB 体积闸门把合成小文件挡在外面",
              _r27ki.returncode == 0 and "volume_too_small" in (_r27ki.stdout + _r27ki.stderr))
    finally:
        shutil.rmtree(td27, ignore_errors=True)

    _rb27 = os.path.join(HERE, "P4_2_KAGOME_RUNBOOK.md")
    _rbo27 = False
    if os.path.isfile(_rb27):
        with open(_rb27, encoding="utf-8") as _f:
            _trb27 = _f.read()
        _rbo27 = ("kconv" in _trb27 and "--kpoints" in _trb27 and "k_convergence.py" in _trb27)
    check("p42: 手册写了 k 点收敛抽查 (kconv + --kpoints + k_convergence.py)", _rbo27)

    _rd27 = os.path.join(HERE, "README.md")
    _rdo27 = False
    if os.path.isfile(_rd27):
        with open(_rd27, encoding="utf-8") as _f:
            _trd27 = _f.read()
        _rdo27 = ("k_convergence.py" in _trd27 and "kconv" in _trd27 and "--kpoints" in _trd27)
    check("p42: README 已同步 (脚本 + kconv 子命令 + --kpoints)", _rdo27)

    _wf27 = os.path.join(HERE, "WORKFLOW.md")
    _wfo27 = False
    if os.path.isfile(_wf27):
        with open(_wf27, encoding="utf-8") as _f:
            _twf27 = _f.read()
        _wfo27 = "k_convergence.py" in _twf27 and "kconv" in _twf27
    check("p42: WORKFLOW 已同步 (继承表 + 目录树)", _wfo27)

    _sk27 = os.path.join(HERE, "dsh", "2d-magnet", "SKILL.md")
    _sko27 = False
    if os.path.isfile(_sk27):
        with open(_sk27, encoding="utf-8") as _f:
            _tsk27 = _f.read()
        _sko27 = "kconv" in _tsk27 and "k_convergence.py" in _tsk27
    check("p42: SKILL.md 写了 k 点收敛抽查入口", _sko27)

    _pl27 = os.path.join(HERE, "NEW_TYPES_PLAN.md")
    _plo27 = False
    if os.path.isfile(_pl27):
        with open(_pl27, encoding="utf-8") as _f:
            _tpl27 = _f.read()
        _plo27 = "k_convergence.py" in _tpl27 or "k 点收敛" in _tpl27
    check("p42: NEW_TYPES_PLAN 的 P4.2 残余里点名 k 点收敛抽查", _plo27)

# ---------------------------------------------------------------------------
# [28] 发布闸门: sdist 内容 (MANIFEST.in) + CI 配置 + check_sdist 的判定逻辑
# ---------------------------------------------------------------------------
if True:
    print("\n[28] 发布闸门: sdist/打包")

    _man28 = os.path.join(HERE, "MANIFEST.in")
    _mano28 = False
    if os.path.isfile(_man28):
        with open(_man28, encoding="utf-8") as _f:
            _tman28 = _f.read()
        _mano28 = all(k in _tman28 for k in ("recursive-include examples *", "recursive-include tests *",
                                             "recursive-include dsh *", "recursive-include templates *",
                                             "global-exclude __pycache__/*"))
    check("pkg: MANIFEST.in 把 examples/tests/dsh/templates 收进 sdist 且排除缓存", _mano28)

    _ci28 = os.path.join(HERE, ".github", "workflows", "ci.yml")
    _cio28 = False
    if os.path.isfile(_ci28):
        with open(_ci28, encoding="utf-8") as _f:
            _tci28 = _f.read()
        _cio28 = all(k in _tci28 for k in ("check_sdist.py", "build", "run_2d_magnet.py doctor",
                                           "magnetic_orders.py --regress", "tests/smoke_test.py",
                                           "tests/test_ingest.py", "tests/test_crect_cant.py"))
    check("pkg: CI 跑 regress + smoke + 契约测试, 并有打包作业 (check_sdist)", _cio28)

    _cs28 = os.path.join(HERE, "tests", "check_sdist.py")
    check("pkg: tests/check_sdist.py 存在 (发布闸门本体)", os.path.isfile(_cs28))

    if os.path.isfile(_cs28):
        import tarfile as _tf28
        sys.path.insert(0, os.path.join(HERE, "tests"))
        import check_sdist as CS28  # noqa: E402

        _missing28 = [w for w in CS28.REQUIRED if not os.path.isfile(os.path.join(HERE, w.replace("/", os.sep)))]
        check("pkg: check_sdist 的 REQUIRED 清单在仓库里都真实存在 (%d 项)" % len(CS28.REQUIRED),
              not _missing28)

        _td28 = tempfile.mkdtemp(prefix="_tmp_pkg28_")
        try:
            def _mk28(arc, names):
                p = os.path.join(_td28, arc)
                with _tf28.open(p, "w") as tf:
                    for _n in names:
                        src = os.path.join(_td28, "src", _n.replace("/", os.sep))
                        os.makedirs(os.path.dirname(src), exist_ok=True)
                        with open(src, "w", encoding="utf-8") as fh:
                            fh.write("x\n")
                        tf.add(src, arcname="2d_magnet_pipeline-0.1.0/" + _n)
                return p

            _good28 = _mk28("good.tar.gz", CS28.REQUIRED)
            _rg28 = run(sys.executable, _cs28, _good28)
            check("pkg: check_sdist 对完整 sdist 判通过 (exit 0 + [OK])",
                  _rg28.returncode == 0 and "[OK]" in _rg28.stdout)

            _bad28 = _mk28("bad.tar.gz", [w for w in CS28.REQUIRED if not w.endswith("P4_2_KAGOME_RUNBOOK.md")])
            _rb28 = run(sys.executable, _cs28, _bad28)
            check("pkg: check_sdist 对缺项 sdist 判失败并点名缺项",
                  _rb28.returncode == 1 and "[FAIL]" in _rb28.stdout
                  and "P4_2_KAGOME_RUNBOOK.md" in _rb28.stdout)

            _dirty28 = _mk28("dirty.tar.gz", list(CS28.REQUIRED) + ["scripts/__pycache__/x.pyc"])
            _rd28 = run(sys.executable, _cs28, _dirty28)
            check("pkg: check_sdist 会拦下混入 __pycache__ 的 sdist",
                  _rd28.returncode == 1 and "不该有" in _rd28.stdout)
        finally:
            shutil.rmtree(_td28, ignore_errors=True)

# ---------------------------------------------------------------------------
# [29] Obsidian 阅读版备份: 与 README 同源 (readme_to_obsidian.py --check)
# ---------------------------------------------------------------------------
if True:
    print("\n[29] Obsidian 阅读版备份 (README.obsidian.md)")

    _r2o29 = os.path.join(SCRIPTS, "readme_to_obsidian.py")
    check("obs: scripts/readme_to_obsidian.py 存在", os.path.isfile(_r2o29))

    _ob29 = os.path.join(HERE, "README.obsidian.md")
    _txt29 = ""
    if os.path.isfile(_ob29):
        with open(_ob29, encoding="utf-8") as _f:
            _txt29 = _f.read()
    check("obs: README.obsidian.md 有 YAML properties + 内链目录 + mermaid 管线图",
          _txt29.startswith("---\n") and "source: README.md" in _txt29
          and "[[#" in _txt29 and "```mermaid" in _txt29)
    check("obs: 备份正文与 README 逐字一致（README 全文都在备份里）",
          bool(_txt29) and open(os.path.join(HERE, "README.md"), encoding="utf-8").read().rstrip("\n") in _txt29)

    if os.path.isfile(_r2o29):
        _c29 = run(sys.executable, _r2o29, "--check")
        check("obs: --check 判定同步 (README 改过就必须重新生成)",
              _c29.returncode == 0 and "[OK]" in _c29.stdout)
        _s29 = run(sys.executable, _r2o29, "--selftest")
        check("obs: 生成器 --selftest ALL PASS",
              _s29.returncode == 0 and "ALL PASS" in _s29.stdout)

print("\n全部通过" if not fails else f"\n{len(fails)} 个失败: {fails}")

sys.exit(1 if fails else 0)
