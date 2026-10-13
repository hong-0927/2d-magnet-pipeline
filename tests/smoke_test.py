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
from run_2d_magnet import (SCRIPTS as DRV_SCRIPTS, _scripts_dir, detect_lattice,  # noqa: E402
                           multi_sublattice_warning, peak_report, summarize_mc,
                           write_workflow, zoom_window)

FIXTURES = os.environ.get("MAGNET_FIXTURES", os.path.join(HERE, "tests", "fixtures"))
MNPS3 = os.path.join(FIXTURES, "MnPS3_POSCAR")
NB3I8 = os.path.join(FIXTURES, "Nb3I8_POSCAR")

for _name, _path in (("MnPS3", MNPS3), ("Nb3I8", NB3I8)):
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
    fake_args = types.SimpleNamespace(vacuum=15.0, mc_refine_nt=60)
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

print("\n全部通过" if not fails else f"\n{len(fails)} 个失败: {fails}")
sys.exit(1 if fails else 0)
