#!/usr/bin/env python
"""run_2d_magnet.py
二维磁性单层「结构 -> 磁耦合常数 -> Monte Carlo」通用流水线入口。

流程 (每一步调用 scripts/ 下对应脚本):
  1. fetch      : Materials Project 下载单层结构 (--formula / --mp-id), 或跳过 (--poscar)
  2. reorient   : 层法线对齐 z + 加真空
  3. supercell  : 蜂窝 2x2x1 超胞 + 四磁态 (FM/Neel/zigzag/stripy) + J1/J2/J3 系数
                 (簇磁体: build_cluster_cases.py 枚举 fm+波矢态 + J_intra/J1/J2/J3)
                 (三角: build_magnetic_cases_tri.py 枚举 fm+120°非共线+stripy/uud + J1/J2/J3)
  4. dft        : 组装 INCAR/KPOINTS/POTCAR.spec + runs/<order>/ (三角用非共线 INCAR)
  5. mc         : POSCAR -> MC 模型 (含 J3) -> Heisenberg MC -> T_C/T_N
  6. render     : 磁序俯视图 (蜂窝四磁序 / 簇磁体簇心视图 / 三角 120° 矢量箭头)
  7. workflow   : 写出 <outdir>/workflow.md
  8. ship-mc    : 复制 mc_2d.py + poscar_to_mc3.py + write_hamiltonian.py (簇/三角再加对应枚举+反解脚本) 到 <outdir>/05_mc/
  9. hamiltonian: 生成 <outdir>/hamiltonian.md (Heisenberg 哈密顿卡片)

用法 (最简: 用文献 J 表 + MP 下载):
  python run_2d_magnet.py --formula MnPS3 --literature MnPS3

  # 或手动给全部参数 (新材料):
  python run_2d_magnet.py --formula VI3 --element V --nonmag I --S 1.0 --moment 2.0 --U 3.0 \
      --J1 -1.0 --J2 -0.1 --J3 0.0 --MAE 0.2

参数优先级: 显式 CLI 参数 > --literature 表 > 元素表 ELEMENT_PARAMS (按磁性元素兜底 S/moment/U) > 默认。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np

# 控制台编码兜底: stdout 默认 errors="strict", 而文献注记/警告里含 '⇒' '≈' '×'
# 等在 GBK 控制台下无法编码的字符, print 会直接 UnicodeEncodeError 崩掉整条流水线
# (stderr 默认 backslashreplace 所以不会有这个问题)。这里只放宽错误处理,
# 不改变控制台自身的编码, 非 ASCII 字符最多显示成 '?' 而不是让命令失败。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors="replace")
    except Exception:
        pass

HERE = os.path.dirname(os.path.abspath(__file__))


def _scripts_dir():
    """定位 scripts/ 目录, 三种运行方式都能找到脚本:

    1. 非可编辑 `pip install .`: scripts 是 site-packages 里的真包, `import scripts` 命中;
    2. 可编辑安装: setuptools 把 scripts 映射回源码目录, 同样命中;
    3. 源码树直跑 (`python run_2d_magnet.py`): sys.path[0] = 脚本所在目录 = 仓库根, 命中仓库 scripts/。
    以 `mc_2d.py` 作标记校验, 避免 `scripts` 这个通用名撞上环境里无关的同名包;
    校验失败 (或删了 scripts/__init__.py 且不在仓库里跑) 才回退 __file__ 同目录。
    """
    try:
        import scripts as _pkg
        _d = os.path.dirname(os.path.abspath(_pkg.__file__))
        if os.path.exists(os.path.join(_d, "mc_2d.py")):
            return _d
    except Exception:
        pass
    return os.path.join(HERE, "scripts")


SCRIPTS = _scripts_dir()
LIT = os.path.join(SCRIPTS, "literature_J.json")

QUIET = False

# 常用磁性离子的经验默认参数: S (自旋), moment (磁矩初猜 μB), U (DFT+U 加在 d 轨道, eV)。
# 仅在「显式 CLI 未给」且「--literature 条目未提供」时作为兜底; 任一显式值都优先。
# U 是经验典型值(非材料本征量), 详见 README「元素参数表」。Fe 故意不列: Fe2+/Fe3+/巡游
# 有效自旋三种取值差异太大 (S=2/2.5/1), 兜底反而会静默给出错误物理, 须显式给 --S/--moment。
ELEMENT_PARAMS = {
    "Mn": {"S": 2.5, "moment": 5.0, "U": 4.0},   # Mn2+ 高自旋 3d5
    "Cr": {"S": 1.5, "moment": 3.0, "U": 3.0},   # Cr3+ 3d3
    "Ni": {"S": 1.0, "moment": 2.0, "U": 5.0},   # Ni2+ 3d8
    "Co": {"S": 1.5, "moment": 3.0, "U": 3.3},   # Co2+ 高自旋 3d7
    "V":  {"S": 1.0, "moment": 2.0, "U": 3.0},   # V3+ 3d2
    "Gd": {"S": 3.5, "moment": 7.0, "U": 6.0},   # Gd2+ 4f7 (S=7/2); 4f 磁矩, DFT+U 需 LDAUL=3 (与 3d 元素不同)
}

# 命名模型预设: `--model <key>` 一键设定 {化学式, 磁性元素, 非磁元素, 晶格类型, S, moment, U}。
# 与 ELEMENT_PARAMS 的兜底不同, 这里的值来自文献条目 (literature_J.json), 是**材料专属**的
# (如 MnPS3 的 U=5.0 取自 Sivadas 2015, 而非元素表里 Mn 的通用 4.0)。
# 关联的 `literature` 键只作**参考**: 文献 J 用于数量级 sanity-check 与 workflow.md 标注来源,
# 不冒充本项目的计算值 —— 真实 J 须跑 VASP (03_dft/runs 导出后手动提交) + `solve` 反解,
# 或用 `--J1/--J2/--J3` 显式给定。VASP 侧始终只导出输入 (INCAR/KPOINTS/POTCAR.spec/POSCAR),
# 不替用户提交作业。
MODEL_PRESETS = {
    "mnps3": {
        "formula": "MnPS3", "element": "Mn", "nonmag": ["P", "S"],
        "lattice": "honeycomb", "S": 2.5, "moment": 5.0, "U": 5.0,
        "literature": "MnPS3",
        "desc": "蜂窝 Néel 反铁磁, Mn2+ S=5/2 (Sivadas 2015)",
    },
    "cri3": {
        "formula": "CrI3", "element": "Cr", "nonmag": ["I"],
        "lattice": "honeycomb", "S": 1.5, "moment": 3.0, "U": 3.0,
        "literature": "CrI3",
        "desc": "蜂窝铁磁, Cr3+ S=3/2 (文献 J 待核实)",
    },
    "crbr3": {
        "formula": "CrBr3", "element": "Cr", "nonmag": ["Br"],
        "lattice": "honeycomb", "S": 1.5, "moment": 3.0, "U": 3.0,
        "literature": "CrBr3",
        "desc": "蜂窝铁磁, Cr3+ S=3/2 (文献 J 待核实)",
    },
    "fe3gete2": {
        "formula": "Fe3GeTe2", "element": "Fe", "nonmag": ["Ge", "Te"],
        "lattice": "cluster", "S": 1.0, "moment": 2.0, "U": 4.0,
        "literature": "Fe3GeTe2",
        "desc": "⚠ 巡游铁磁金属 (多亚晶格), 几何判据误判成 cluster; "
                "Fe 的 S/moment 需自洽 (元素表不列 Fe)",
    },
    "crsbr": {
        "formula": "CrSBr", "element": "Cr", "nonmag": ["S", "Br"],
        "lattice": "rectangular", "S": 1.5, "moment": 3.0, "U": 4.0,
        "literature": "CrSBr",
        "desc": "矩形铁磁, Cr3+ S=3/2 (Yang&Wu PRB 104,144416 (2021): 距离壳升序 "
                "J1(a)=2.79/J2(b)=2.12/J3(对角)=3.88 meV, 单层 T_C≈175 K)",
    },
    "gdi2": {
        "formula": "GdI2", "element": "Gd", "nonmag": ["I"],
        "lattice": "triangular", "S": 3.5, "moment": 7.0, "U": 6.0,
        "literature": "GdI2",
        "desc": "三角铁磁, Gd2+ 4f7 S=7/2 (4f-5d 巡游中介交换, 单层 T_C≈241 K; J 需 DFT 反解)",
    },
    "la2cuo4": {
        "formula": "La2CuO4", "element": "Cu", "nonmag": ["La", "O"],
        "lattice": "square", "S": 0.5, "moment": 1.0, "U": 7.0,
        "literature": "La2CuO4",
        "desc": "方形反铁磁 G 型, Cu2+ S=1/2 —— 方形格子的**定量锚点** "
                "(Coldea et al., PRL 86, 5377 (2001) 295 K 自由拟合: J=111.8, J'=-11.4 meV; "
                "转成流水线约定 J1=-111.8/J2=+11.4/J3=-2.0)。"
                "注意 T_N=325 K 是 3D 有序温度, 单层 MC 只能当量级检查",
    },
    "k2nif4": {
        "formula": "K2NiF4", "element": "Ni", "nonmag": ["K", "F"],
        "lattice": "square", "S": 1.0, "moment": 2.0, "U": 4.0,
        "literature": "K2NiF4",
        "desc": "方形反铁磁 G 型, Ni2+ S=1 —— 方形格子的**结构原型** (理想胞就是真实胞, "
                "无 Jahn-Teller 畸变/八面体倾转, 壳距离 a/a√2/2a 严格成立)。"
                "J 未核实 (文献表 J1/J2/J3=null) ⇒ 只用于几何/规格层检查, "
                "跑 MC 必须显式 --J1/--J2/--J3; 定量锚点见 la2cuo4",
    },
    "crocl": {
        "formula": "CrOCl", "element": "Cr", "nonmag": ["O", "Cl"],
        "lattice": "rect_centered", "S": 1.5, "moment": 3.0, "U": 4.0,
        "literature": "CrOCl",
        "desc": "⚠ C 心矩形 (2 磁位/胞, 正交 a≠b, Pmmn No.59; a=3.198/b=3.860 A) —— "
                "**多亚晶格**支持的理由: 最近邻 J1(3.0347 A) 是 C 心矢量, 连的是同一胞内的两个 Cr, "
                "没有 afm_c 这个 basis 序就无法把 J1 从 E0 里分出来。"
                "文献 J 未核实 (literature_J.json 里 J1/J2/J3=null) ⇒ 跑 MC 必须显式 "
                "--J1/--J2/--J3; solve --demo 也需要 --synthetic 或已核实的 --ref",
    },
    "fe3sn": {
        "formula": "Fe3Sn", "element": "Fe", "nonmag": ["Sn"],
        "lattice": "kagome", "S": 1.0, "moment": 2.0, "U": 4.0,
        "literature": None,
        "desc": "kagome 角共享三角形网 (3 磁位/胞, 六角 a=b=5.4 A; Sn 在六角中心) —— "
                "最近邻 a/2 每个磁位点有 4 个 (角共享三角形网)。"
                "**注意别把 J1/J1′ 说成「同距离两条不等价键」**: 在**理想** kagome 里 a/2 只有一条键、"
                "只有一个 J1(两种「键型」是同一根键按所在三角形的**取向**说的, 不是两个距离)。"
                "所谓「J1 是 J1 与 J1′ 的平均」是说**真实**(呼吸 kagome) 材料的两类三角边长不同、"
                "而本模型的几何看不见这个差别 ⇒ 反解出的 J1 是平均值。"
                "⚠ q=0 一族 (fm/q0_120/q0_ud/q0_duu) 对 J1/J2 **简并** (rank=2, cond≈2.3e34) ⇒ "
                "必须用含 √3×√3 的序组合 (builder 的 selected); rank<4 时 solve 端直接报错。"
                "S=1.0 是占位值 (Fe3Sn 的第一性 S/moment 未定, 元素表不列 Fe)。"
                "**文献条目已建 (literature_J.json 的 Fe3Sn), 但 J1/J2/J3 仍为 null** —— 已核实来源只给了 "
                "「它是铁磁体 + 有序温度 + 各向异性」, 没给壳分辨的 J ⇒ 跑 MC 仍必须显式 --J1/--J2/--J3; "
                "该条目带 T_C_K=705 (居里, 不是奈尔!) 与结构/各向异性锚点, 供 P4.2 验收比对。",
    },
}


def check_presets():
    """预设自检: 每个 --model 预设若声明了 literature 锚点, 必须与文献表逐项一致。

    返回问题字符串列表 (空 = 通过)。文献表缺失/损坏不在这里报 (doctor 的 lit_missing 负责),
    这样 fetch/build 等不依赖文献表的阶段不会被预设自检拖死。
    比较字段: formula / element / nonmag(顺序无关) / S / moment / U。
    """
    try:
        with open(LIT, encoding="utf-8") as f:
            tab = json.load(f)
    except Exception:
        return []
    probs = []
    for key, preset in sorted(MODEL_PRESETS.items()):
        lit_key = preset.get("literature")
        if not lit_key:
            continue
        ent = tab.get(lit_key)
        if ent is None:
            probs.append(f"预设 {key!r} 指向文献条目 {lit_key!r}, 但 literature_J.json 里没有该条")
            continue
        for fld in ("formula", "element", "nonmag", "S", "moment", "U"):
            a, b = preset.get(fld), ent.get(fld)
            if fld == "nonmag":
                a = None if a is None else sorted(a)
                b = None if b is None else sorted(b)
            if a != b:
                probs.append(f"预设 {key!r} 与文献条目 {lit_key!r} 的 {fld} 不一致: {a!r} != {b!r}")
    return probs


def log(*a, **k):
    if not QUIET:
        print(*a, **k)


def have_matplotlib():
    """matplotlib 是否真的能 import。

    本机 (Windows, Python 3.12) 装了 matplotlib 3.11.1 但缺 pyparsing / cycler / contourpy /
    fontTools / kiwisolver 且无网络、pip/uv 缓存里也没有这些 wheel => import 必失败。
    超胞图因此在 render_supercell.py 之外多了一条纯 Pillow 的等价路径 render_supercell_pil.py。
    """
    try:
        import matplotlib  # noqa: F401
        return True
    except Exception:
        return False


def child_env():
    """所有子进程统一的 stdio 编码环境: PYTHONIOENCODING=utf-8。

    为什么必须由父进程注入: 子进程默认继承控制台编码 (Windows 中文环境 = cp936),
    脚本 print 出的非 GBK 字符 (如 '⇒' / '≈' / '²') 会让**子进程**以
    `UnicodeEncodeError: 'gbk' codec can't encode character ...` 直接崩; QUIET 模式下
    父进程还用 cp936 解码子进程输出, 同样可能解码失败。历史教训: 文献注记里的 '⇒'
    让 `dft` 阶段在 `run_2d_magnet.py` 打印处炸过一次。
    原先只有 tests/ 下的脚本自己设了, driver 路径没设 —— 而 README / KNOWN_ISSUES
    一直写着"driver 子进程已内置", 这里补上, 让文档与代码一致。
    """
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    return env


def run(pyfile, args, cwd=None):
    """跑 scripts/ 下的脚本。QUIET 时捕获 stdout/stderr, 只在失败或警告时回显。"""
    cmd = [sys.executable, os.path.join(SCRIPTS, pyfile)] + args
    log("\n$", " ".join(cmd))
    if QUIET:
        r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                           encoding="utf-8", errors="replace", env=child_env())
        if r.returncode != 0:
            print(f"[FAIL] {pyfile} 退出码 {r.returncode}\n"
                  f"$ {' '.join(cmd)}\n{r.stdout}{r.stderr}", file=sys.stderr)
            sys.exit(r.returncode)
        if r.stderr and r.stderr.strip():
            print(r.stderr.rstrip(), file=sys.stderr)
    else:
        r = subprocess.run(cmd, cwd=cwd, env=child_env())
        if r.returncode != 0:
            print(f"[FAIL] {pyfile} 退出码 {r.returncode}", file=sys.stderr)
            sys.exit(r.returncode)


def try_mp_fetch(outdir, formula, mp_id):
    """跑 mp_fetch_2d.py (MP 下载 + 剥单层); 成功返回选中 POSCAR 路径, 失败返回 None (不退出)。"""
    fetch_dir = os.path.join(outdir, "_fetch_tmp")
    os.makedirs(fetch_dir, exist_ok=True)
    # mp_fetch_2d.py 固定写脚本同目录 structures/ + candidates.csv; 复制到临时目录运行
    shutil.copy(os.path.join(SCRIPTS, "mp_fetch_2d.py"), fetch_dir)
    fa = [os.path.join(fetch_dir, "mp_fetch_2d.py")]
    if mp_id:
        fa += ["fetch", "--ids", mp_id]
    else:
        fa += ["fetch", "--formula", formula]
    log("\n$", " ".join([sys.executable] + fa))
    kw = {"capture_output": True, "text": True, "encoding": "utf-8",
          "errors": "replace"} if QUIET else {}
    try:
        r = subprocess.run([sys.executable] + fa, cwd=fetch_dir, env=child_env(), **kw)
    except Exception as e:  # noqa: BLE001
        print(f"[fetch] MP 运行异常: {e}", file=sys.stderr)
        return None
    if r.returncode != 0:
        print(f"[fetch] MP 退出码 {r.returncode}: {(r.stderr or r.stdout or '')[-300:]}", file=sys.stderr)
        return None
    try:
        return pick_structure(os.path.join(fetch_dir, "structures"), formula, mp_id)
    except SystemExit as e:
        print(f"[fetch] MP 选结构失败: {e}", file=sys.stderr)
        return None


# ---- error catalog: 已知失败签名 -> 症状 -> 修复 (doctor 联动 + KNOWN_ISSUES.md) ----
# static=True 是「运行时才出现、环境自检探不到」的条目, doctor 作为静态参考列出;
# static=False 是可被环境自检命中的条目, 由 doctor 的检查项按 key 联动。
ERROR_CATALOG = [
    {"key": "py39", "static": False,
     "symptom": "语法/类型报错, 或 numpy/numba 版本不兼容",
     "fix": "升级到 Python >= 3.9 (如 py -3.12 或 conda)"},
    {"key": "missing_dep", "static": False,
     "symptom": "ImportError: No module named 'xxx'",
     "fix": "pip install -r requirements.txt"},
    {"key": "packaging", "static": False,
     "symptom": "render 阶段: ModuleNotFoundError: No module named 'packaging' "
               "(matplotlib 顶层 import 连带需要它)",
     "fix": "pip install packaging"},
    {"key": "fig_backend", "static": False,
     "symptom": "doctor 报「缺少依赖 matplotlib」或「缺少 Pillow」; render 阶段 "
               "`[warn] 画图跳过: No module named 'pyparsing'`, figures/ 空着没有 PNG",
     "fix": "6 条渲染路径都有 matplotlib / 纯 Pillow 两条等价实现, driver 或脚本自动二选一 "
            "(超胞 render_supercell_pil.py、簇磁体 render_cluster_orders_pil.py、三角 "
            "render_tri_orders_pil.py、MC 曲线 mc_2d.py 的 _plot_mt_png、FSS 标度图 "
            "mc_fss.py 的 _plot_fss_png、FGT 图 build_fgt_cases.py 的 _plot_fgt_png) "
            "⇒ 装 Pillow 即可出图。要修 matplotlib 需补齐 "
            "pyparsing/cycler/contourpy/fontTools/kiwisolver (matplotlib 在 import 期就要它们)"},
    {"key": "scripts_missing", "static": False,
     "symptom": "driver 报缺脚本 / 子阶段找不到脚本",
     "fix": "确认 scripts/ 目录完整 (重新 git clone 或解压分发包)"},
    {"key": "lit_missing", "static": False,
     "symptom": "--literature <key> 报找不到文献条目 / literature_J.json 缺失",
     "fix": "确认 scripts/literature_J.json 存在且键名正确"},
    {"key": "preset_mismatch", "static": False,
     "symptom": "--model 预设自检报「预设 X 与文献条目 Y 的 <字段> 不一致」/ 预设指向不存在的文献条",
     "fix": "把 MODEL_PRESETS[X] 与 literature_J.json[Y] 的 formula/element/nonmag/S/moment/U 对齐"
            " (预设是材料专属值, 不跟随文献表漂移); 回归见 tests/test_crect_cases.py 的预设自检两例"},
    {"key": "mp_key", "static": False,
     "symptom": "fetch 阶段 (--formula 下载) 报 401/403 或空结果",
     "fix": "设 MP_API_KEY; 或 --poscar 用本地结构, 或 --provider cod 走公开库"},
    {"key": "cwd_ro", "static": False,
     "symptom": "PermissionError: 建 --outdir / 写文件失败",
     "fix": "cd 到可写目录, 或 --outdir 指到可写路径 (如系统临时目录)"},

    {"key": "static_incar", "static": True,
     "symptom": "runs/<order>/INCAR 无 IBRION/NSW -> 静默假弛豫, OUTCAR 其实是静态总能",
     "fix": "用 03_dft/INCAR.relax 做弛豫; 或 relax_backfill 回填后用位移检查核对"},
    {"key": "atom_mismatch", "static": True,
     "symptom": "relax_backfill 报「原子数/物种不符」并拒绝",
     "fix": "核对 CONTCAR 与 POSCAR.pristine 的物种与原子数一致后再回填"},
    {"key": "gbk", "static": True,
     "symptom": "子进程输出 UnicodeDecodeError (GBK 编码)",
     "fix": "设 PYTHONIOENCODING=utf-8 再跑 (driver 子进程已内置; 手动跑脚本需自设)"},
    {"key": "solve_no_outcar", "static": True,
     "symptom": "solve_J3 报「需要 fm/neel/zigzag/stripy 四态 OUTCAR」",
     "fix": "先跑 DFT 生成各磁态 OUTCAR; 或加 --demo 自检反解"},
    {"key": "fss_no_cross", "static": True,
     "symptom": "FSS 报「未分辨」/ fss_fss.json 里 tc=null (尺寸间 U4 无干净交点)",
     "fix": "加大 --fss-sizes (默认 12 18 24 30) 或 --mc-seeds 压噪"},
    {"key": "peak_boundary", "static": True,
     "symptom": "T_C/T_N 标「未分辨」(峰落在温度窗口首/末格点)",
     "fix": "已内置自动缩放重扫; 或调 --mc-tmax / --mc-nt 覆盖真峰"},
    {"key": "soc_required", "static": True,
     "symptom": "assemble_dft.py 报「磁序 <order> 需要 SOC 但没给 --lsorbit」/「magmom.txt 是 3 分量"
                "但当前是共线模式」/「--lsorbit 需要 --noncollinear」",
     "fix": "driver 侧加 --noncollinear (弱铁磁倾斜/DMI/kagome 再加 --soc, 可配 --saxis); "
            "手动跑脚本时给 --noncollinear --lsorbit。这是 P3 的硬门禁 (§115 不静默降级): "
            "缺 SOC 时倾斜/螺旋构型物理上不成立, 不会按共线近似给你一个错的 J"},
    {"key": "sq_centered", "static": True,
     "symptom": "driver 报「检测到「方形 + 多磁原子/胞」(sq_centered)」并退出",
     "fix": "该网 (a≈b 且胞内多个磁原子被 cutoff 并簇) 本仓库未支持; 别硬套 square "
            "(会丢掉两个 basis 的反平行自由度 -> 错的 J)。要支持需新写方形版 basis 序 + 显式壳距 + "
            "独立反解脚本"},
    {"key": "crect_singular", "static": True,
     "symptom": "build_crect_cases.py 报 det<1e-6 / cond>1e4 / 「母序与其 basis 翻倍版同选」",
     "fix": "去掉退化组合 (afm_c_a 等于 afm_a, afm_c_diag 等于 afm_diag), DFT 只提交 "
            "fm/afm_c/afm_a/afm_b; 壳容差保持 0.05 (tol=0.2 会把 J1=3.0347 与 J2=3.1980 并壳且不报错)"},
    {"key": "mc_shell_mismatch", "static": True,
     "symptom": "driver 报「[check] MC 每格点键数 ... != 规格层 per_site_bonds」并 SystemExit",
     "fix": "说明 J1/J2/J3 很可能被绑到了错的壳上 (距离自动排序与规格层不一致): 核对 "
            "build_*_cases.py 的壳定义与 poscar_to_mc3 的口径 (CrOCl 应 4/2/2, La2CuO4 应 4/4/4)"},
    {"key": "energy_convention", "static": True,
     "symptom": "solve_* 报「no 'free energy TOTEN' in <OUTCAR>」, 或 dft_energy.mixed_convention "
                "报「两种能量口径混用」",
     "fix": "P3 起所有 solve_* 只认 OUTCAR 的 free energy TOTEN (取最后一次); 缺它说明计算没跑完或 "
            "输出被裁剪 -> 补跑/换完整 OUTCAR (完全没有 OUTCAR 时才走 E.txt 兜底)。不要用 "
            "energy(sigma->0) 凑数: 与 free energy TOTEN 混用会让 J 差一个常数项 (阈值 1e-3 eV)"},
    {"key": "kagome_degenerate", "static": True,
     "symptom": "solve_J_kagome 报「设计矩阵 rank=... < 4 ⇒ 这几个磁序**无法唯一确定** J1/J2/J3」, "
                "或报「设计矩阵病态: cond=... (>1e4)」; 也可能报「selected_orders=... 与本脚本期望的 "
                "不一致」",
     "fix": "kagome 的 q=0 一族 (fm/q0_120/q0_ud/q0_duu) 满足 c1 == c2 且 c3 相同 ⇒ rank=2 "
            "(cond≈2.3e34), J1/J2 在那个子块里不可分离。必须用含 √3×√3 的组合 "
            "(builder 的 selected = fm/q0_120/sq3_updntri/sq3_alt, rank=4, cond=98.51), "
            "且超胞 n0/n1 是 3 的倍数。换序后必须重算 rank/cond —— "
            "「cond 没超阈」不等于「三个壳都分得开」"},
    {"key": "dmi_needs_shell_def", "static": True,
     "symptom": "想跑含 DMI（反称交换）的格点模型: `poscar_to_mc3 --jfile J_result.json` 只用标量 J "
                "⇒ 反称项**全部丢失**、模型退化成纯 Heisenberg，而输出看起来一切正常（不报错）",
     "fix": "DMI 只有 3×3 `J_mat` 装得下，而 `J_result.json` 目前只有标量 `J_meV` ⇒ 必须走 "
            "`--shell-def` 并显式给 `J_mat`: 用 `scripts/dmi_matrix.py --D ... --J ... --S ... --out shell.json` "
            "(J_D = -[D]× 对应 H_DM = +D·(S_i×S_j), 矩阵已吸收 S²), 再 "
            "`poscar_to_mc3.py ... --shell-def shell.json`; 两个脚本的 `--S` 必须一致。"
            "D 本身**必须来自文献或你自己的 DFT** —— 本仓库不提取 D (见 DMI_MTI_MEMO.md §3.2/§3.3), "
            "且引擎无外场、无拓扑荷 ⇒ 只能看零场基态能, 不许报 skyrmion 相图 / 拓扑荷"},
    {"key": "sublattice_mask_invalid", "static": True,
     "symptom": "`mc_2d` / `mc_fss` / 驱动 `--mc-order` 报「--order sublattice:<mask> 需要掩码」/「掩码有空项」/"
                "「子晶格编号越界」/「掩码 [...] 覆盖了全部 N 个子晶格 ⇒ B 组为空」",
     "fix": "`sublattice:<mask>` 的 mask 是**逗号分隔的子晶格编号**（不写就是空掩码）: "
            "`sublattice:0` = 老的 `staggered`（双亚晶格体系逐位等价），`sublattice:0,1` 只在 num_basis >= 3 时合法"
            "（两套亚晶格模型给出 sublattice:0,1 会把 B 组清空 ⇒ mbx/cb 除零，脚本硬报错而不静默）。"
            "编号是 `i % num_basis` 的 basis 序号（来自 POSCAR 里的原子顺序），不是超胞格点号；"
            "掩码**至少留一个**子晶格在 B 组。基础模型只有 2 个 basis 时不要试图分组三角/kagome 的三个子晶格"},
    {"key": "cant_theta_needs_crect", "static": True,
     "symptom": "驱动报「给了 --cant-theta 但 lattice=...，不是 C 心矩形 (rect_centered)」或"
                "「给了 --cant-theta 却用 --cant-base fm」",
     "fix": "倾斜家族只在 C 心矩形网（CrOCl 类，2 磁原子/胞）上定义：`--cant-theta \"8 16\"` 配 "
            "`--model crocl`（或 `--lattice rect_centered`）。母序取 `afm_c`（默认）/`afm_a`/`afm_b` —— "
            "`fm` 被拒，因为全同向母序上的 ±θ 只是一次整体自旋旋转、系数行逐位等于 fm 行（不带信息）。"
            "倾斜构型要 `--soc` 装配（缺它 assemble_dft 会硬拒，不会静默按共线算）；"
            "θ 是**输入参数**，不是自洽解出来的基态角，φ 由几何定不出来（见 README 局限 17）"},
]


def pkg_version():
    """包版本号: 声明处只有一个 —— pyproject.toml 的 [project].version。

    优先读已安装元数据 (pip install 之后它是权威来源), 拿不到就解析与本文件同目录的
    pyproject.toml (源码树直接跑的情况); 两条都失败返回 "unknown" —— 不在这里写死一个
    版本号, 否则就有第二个"权威声明"了 (README 明确说版本只在 pyproject.toml)。
    """
    try:
        from importlib.metadata import version as _pkg_v
        return _pkg_v("2d-magnet-pipeline")
    except Exception:  # noqa: BLE001  (未安装 / 元数据缺失)
        pass
    try:
        import tomllib
        p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pyproject.toml")
        with open(p, "rb") as f:
            return tomllib.load(f)["project"]["version"]
    except Exception:  # noqa: BLE001  (Python < 3.11 无 tomllib / 文件缺失)
        return "unknown"


def doctor():
    """环境自检 (doctor 子命令): python / 依赖 / 脚本 / 文献表 / MP key / CWD 写权限,
    并把命中的已知问题映射到 ERROR_CATALOG 给「症状 -> 修复」。

    返回 True 表示硬性项全部通过; 缺依赖/缺脚本/不可写返回 False (exit 1), 便于 CI `&&` 短路。
    """
    import importlib
    ok = True
    hits = []  # 命中的 catalog key (可能重复, 打印前会去重)

    def line(good, msg, key=None):
        nonlocal ok
        if not good:
            ok = False
            if key:
                hits.append(key)
        print(("  [OK]   " if good else "  [WARN] ") + msg)

    print(f"python : {sys.version.split()[0]}   ({sys.executable})")
    print(f"pipeline: 2d-magnet {pkg_version()}   ({os.path.dirname(os.path.abspath(__file__))})")
    line(sys.version_info >= (3, 9),
         f"Python >= 3.9 (当前 {sys.version.split()[0]})", key="py39")
    for pkg in ("numpy", "numba", "matplotlib", "requests"):
        try:
            m = importlib.import_module(pkg)
            line(True, f"{pkg} {getattr(m, '__version__', '?')}")
        except Exception:
            if pkg == "matplotlib":
                line(False, "缺少依赖 matplotlib -> 全部 6 条渲染路径自动改走纯 Pillow 兜底 "
                            "(render_supercell_pil.py / render_cluster_orders_pil.py / "
                            "render_tri_orders_pil.py / mc_2d.py 的 _plot_mt_png / "
                            "mc_fss.py 的 _plot_fss_png / build_fgt_cases.py 的 _plot_fgt_png)",
                     key="fig_backend")
            else:
                line(False, f"缺少依赖 {pkg} -> pip install -r requirements.txt",
                     key="missing_dep")

    # Pillow 专项: matplotlib 缺席时它是全部渲染路径的唯一出路 (本机就是这样)
    try:
        from PIL import __version__ as _pil_ver
        line(True, f"Pillow {_pil_ver} (matplotlib 缺席时 6 条渲染路径的兜底)")
    except Exception:
        line(False, "缺少 Pillow -> 若 matplotlib 也不可用则所有图都出不来 "
                    "(pip install Pillow)",
             key="fig_backend")

    # packaging 专项: matplotlib 顶层 import 会连带需要它, 缺了 render 直接崩
    try:
        importlib.import_module("packaging")
        line(True, "packaging (matplotlib 顶层 import 依赖)")
    except Exception:
        line(False, "缺少 packaging -> render 阶段会崩, pip install packaging",
             key="packaging")

    need = ["reorient.py", "mc_2d.py", "poscar_to_mc3.py", "write_hamiltonian.py",
            "build_magnetic_cases3.py", "build_cluster_cases.py",
            "build_magnetic_cases_tri.py", "mp_fetch_2d.py"]
    miss = [s for s in need if not os.path.exists(os.path.join(SCRIPTS, s))]
    line(not miss, "scripts/ 关键脚本齐全" if not miss else f"缺脚本: {miss}",
         key="scripts_missing")

    if os.path.exists(LIT):
        try:
            with open(LIT, encoding="utf-8") as f:
                tab = json.load(f)
            line(True, f"literature_J.json: {len(tab)} 条 ({', '.join(tab.keys())})")
        except Exception as e:
            line(False, f"literature_J.json 解析失败: {e}", key="lit_missing")
    else:
        line(False, f"literature_J.json 不存在: {LIT}", key="lit_missing")

    # 预设自检: --model 预设是材料专属值, 声明了 literature 锚点的必须与文献表逐项一致。
    problems = check_presets() if os.path.exists(LIT) else None
    if problems:
        for p in problems:
            line(False, p, key="preset_mismatch")
    elif problems == []:
        n_lit = sum(1 for v in MODEL_PRESETS.values() if v.get("literature"))
        line(True, f"--model 预设自检: {len(MODEL_PRESETS)} 个预设 (含 {n_lit} 个文献锚点) "
                   f"与 literature_J.json 一致")

    # MP key 与 fixtures 是条件项 (非硬失败): 缺 MP key 只挡 fetch, 缺 fixtures 只挡回归测试。
    key = os.environ.get("MP_API_KEY")
    if key:
        line(True, "MP_API_KEY 已设置")
    else:
        print("  [INFO] MP_API_KEY 未设置 (仅 --poscar 可用; 下载结构需 --mp-key / --provider cod)")
        hits.append("mp_key")  # 条件项: 命中 catalog 但不置 ok=False

    fx = os.path.join(HERE, "tests", "fixtures")
    if os.path.isdir(fx):
        line(True, "tests/fixtures/ 存在 (回归测试可用)")
    else:
        print("  [INFO] tests/fixtures/ 缺失 (回归测试会跳过)")

    probe = os.path.join(os.getcwd(), "._doctor_probe")
    try:
        with open(probe, "w", encoding="utf-8") as f:
            f.write("x")
        os.remove(probe)
        line(True, f"CWD 可写 ({os.getcwd()})")
    except Exception as e:
        line(False, f"CWD 不可写 ({os.getcwd()}): {e}", key="cwd_ro")

    # ---- 命中条目: 症状 -> 修复 (去重) ----
    by_key = {e["key"]: e for e in ERROR_CATALOG}
    shown = []
    for k in hits:
        if k in by_key and k not in shown:
            shown.append(k)
    if shown:
        print("\n--- 命中 error catalog (症状 -> 修复) ---")
        for k in shown:
            e = by_key[k]
            print(f"  * {e['symptom']}")
            print(f"    -> {e['fix']}")

    # ---- 静态参考条目 (运行时才出现, 环境自检探不到) ----
    print("\n--- 常见运行时问题 (静态参考, 详见 KNOWN_ISSUES.md) ---")
    for e in ERROR_CATALOG:
        if e["static"]:
            print(f"  * [{e['key']}] {e['symptom']}")
            print(f"    -> {e['fix']}")

    return ok


def peak_report(csv):
    """定位峰判据 (Cv 或 χ, 即第 3 列) 在温度网格上的位置。

    返回 dict(T_peak, y_peak, i, n, at_low, at_high, dT, T_lo, T_hi); 读不到返回 None。
    **峰落在首/末格点 = 窗口截断, 不能当作 T_C/T_N 报出** (见 summarize_mc / refine_mc_peak)。
    """
    if not os.path.exists(csv):
        return None
    try:
        d = np.loadtxt(csv, delimiter=",", skiprows=1)
        if d.ndim != 2 or d.shape[1] < 3 or len(d) < 3:
            return None
        T, y = d[:, 0], d[:, 2]
        i = int(np.argmax(y))
        dT = float(np.median(np.abs(np.diff(T))))
        return {"T_peak": float(T[i]), "y_peak": float(y[i]), "i": i, "n": len(T),
                "at_low": i == 0, "at_high": i == len(T) - 1,
                "dT": dT, "T_lo": float(T[0]), "T_hi": float(T[-1])}
    except Exception:
        return None


def zoom_window(pk, floor=0.1):
    """峰落在窗口边界 -> 给出下一轮「向该侧外扩 + 自适应加密」的温度窗口。

    注意低侧情形: 真峰可能只是在 [T_lo, T_lo+dT] 之间被粗网格**跨过**了 (不是位于 T_lo 之下),
    所以外扩的同时必须保留窗口上方 4 个旧格距, 否则又会把真峰排除在新窗口之外。
    """
    T_lo, T_hi, dT = pk["T_lo"], pk["T_hi"], pk["dT"]
    if dT <= 0:
        return None
    if pk["at_low"]:
        lo = max(floor, T_lo / 4.0)
        hi = min(T_hi, T_lo + 4.0 * dT)
    else:
        lo = max(floor, T_hi - 4.0 * dT)
        hi = min(T_hi * 4.0, T_hi + 40.0 * dT)
    if hi - lo <= 0:
        return None
    return lo, hi


def verify_mc_shells(outdir, model):
    """校验 MC 模型的每格点键数与规格层一致 —— 防止 J1/J2/J3 被绑到错的壳上。

    两边用的是**不同的几何口径**: 规格层 (magnetic_orders.json) 走 3D 键长, 而
    poscar_to_mc3.py 把 POSCAR 投到 2D 再搜索。同一个 C 心键, 前者 3.0347 A (含
    1.71 A 面外起伏), 后者 2.506 A。距离没法对, 但**每格点键数**是跨口径不变量,
    而 poscar_to_mc3.py 是按距离自动排序把 J1/J2/J3 绑到壳上的 —— 只要那个次序与
    规格层的壳序不同, J2 与 J3 就会静默互换 (反解/拟合照样跑完, 只是物理错了)。
    这里把「靠巧合」变成「被验证」。
    """
    spec_p = os.path.join(outdir, "02_supercell", "cases", "magnetic_orders.json")
    if not (os.path.exists(spec_p) and os.path.exists(model)):
        return
    with open(spec_p, encoding="utf-8-sig") as f:
        spec = json.load(f)
    ps = spec.get("per_site_bonds")
    if not ps:
        return
    with open(model, encoding="utf-8-sig") as f:
        m = json.load(f)
    want = [float(ps.get(n, 0.0)) for n in ("J1", "J2", "J3")]
    sh = m.get("shells", {}) or {}
    got = [float(sh.get(k) or 0.0)
           for k in ("nn_per_site", "nnn_per_site", "nnnn_per_site")]
    if any(abs(a - b) > 1e-6 for a, b in zip(want, got)):
        raise SystemExit(
            f"MC 模型的每格点键数 {got} 与规格层 {want} 不一致 -> J1/J2/J3 可能被绑到了"
            f"错的壳上 (距离自动排序的次序与规格层不同)。\n"
            f"  规格层壳 (3D 键长): {spec.get('shells_angstrom')}\n"
            f"  MC 模型壳 (2D 投影): {[s.get('d') for s in (sh.get('all') or [])]}\n"
            f"  修法: 给 poscar_to_mc3.py 一个显式 --shell-def, 或修正规格层的壳序。")
    log(f"[check] MC 每格点键数 {got} == 规格层 {want} (J1/J2/J3 绑定正确)")


def refine_mc_peak(outdir, model, mc_order, size, nt, max_rounds, enabled=True, seeds=1):
    """峰判据落在温度窗口边界时自动缩放窗口 + 加密重扫 (退火), 直到峰落进窗口内部。

    这是修一个真实缺陷: 默认窗口 5–150 K / nt=30 → 格距 5 K, 会**整个跨过**只有几 K 宽的
    Cv 峰; 取 argmax 就落到窗口边界, 被当成 T_C/T_N 报出 (MnI2 曾因此把 7.4 K 报成 5.0 K)。
    首轮结果保留为 `mc_mt_survey.csv/.png` 作溯源, 规范文件名 `mc_mt.csv/.png` 留给重扫结果。
    返回 dict(survey, final, rounds, resolved) 或 None。
    """
    csv = os.path.join(outdir, "05_mc", "mc_mt.csv")
    survey = peak_report(csv)
    if survey is None:
        return None
    if not (survey["at_low"] or survey["at_high"]):
        return {"survey": survey, "final": survey, "rounds": [], "resolved": True}
    if not enabled or max_rounds <= 0:
        return {"survey": survey, "final": survey, "rounds": [], "resolved": False}

    keep = os.path.join(outdir, "05_mc", "mc_mt_survey.csv")
    if not os.path.exists(keep):
        shutil.copyfile(csv, keep)
        png = os.path.splitext(csv)[0] + ".png"
        if os.path.exists(png):
            shutil.copyfile(png, os.path.splitext(keep)[0] + ".png")

    rounds, pk = [], survey
    for r in range(1, max_rounds + 1):
        w = zoom_window(pk)
        if w is None:
            break
        lo, hi = w
        if rounds and lo == rounds[-1]["window"][0] and hi == rounds[-1]["window"][1]:
            break  # 窗口不再变化, 再跑也不会改善
        side = "下限" if pk["at_low"] else "上限"
        log(f"[mc] 峰判据落在窗口{side} {pk['T_peak']:.2f} K (格距 {pk['dT']:.2f} K) -> "
            f"第 {r} 轮缩放重扫 [{lo:.2f}, {hi:.2f}] K, nt={nt}, 退火")
        run("mc_2d.py", [model, "--order", mc_order, "--size"] + [str(x) for x in size] +
            ["--tmin", f"{lo:.4f}", "--tmax", f"{hi:.4f}", "--nt", str(nt),
             "--seeds", str(seeds),
             "--cool", "--out", os.path.join(outdir, "05_mc", "mc")], cwd=outdir)
        pk2 = peak_report(csv)
        if pk2 is None:
            break
        rounds.append({"round": r, "window": [lo, hi], "nt": nt, "peak": pk2})
        pk = pk2
        if not (pk2["at_low"] or pk2["at_high"]):
            return {"survey": survey, "final": pk2, "rounds": rounds, "resolved": True}
    return {"survey": survey, "final": pk, "rounds": rounds, "resolved": False}


def summarize_mc(outdir, mc_order):
    """从 05_mc/mc_mt.csv 读 T_C/T_N (用于结尾摘要, QUIET 下也能看到结论)。

    峰判据落在温度窗口边界时**不报该边界值** —— 那是窗口截断/网格太粗的假象, 显式声明未分辨。
    """
    csv = os.path.join(outdir, "05_mc", "mc_mt.csv")
    pk = peak_report(csv)
    if pk is None:
        return None
    side = "下限" if pk["at_low"] else "上限"
    win = f"窗口 [{pk['T_lo']:.1f},{pk['T_hi']:.1f}] K, 格距 {pk['dT']:.2f} K"
    if mc_order == "cv":
        if pk["at_low"] or pk["at_high"]:
            return (f"比热 Cv 峰落在窗口{side} {pk['T_peak']:.1f} K → T_N **未分辨** ({win})")
        return f"比热 Cv 峰 T_N ≈ {pk['T_peak']:.2f} K ({win})"
    try:
        d = np.loadtxt(csv, delimiter=",", skiprows=1)
        T = d[:, 0]
        i_chi = int(np.argmax(d[:, 2]))
        i_dm = int(np.argmax(-np.gradient(d[:, 1], T)))
        edge = lambda i: " (落在窗口边界, 未分辨)" if i in (0, len(T) - 1) else ""
        return (f"χ 峰 T ≈ {T[i_chi]:.2f} K{edge(i_chi)}, "
                f"|dm/dT| 最大 T ≈ {T[i_dm]:.2f} K{edge(i_dm)} ({win})")
    except Exception:
        return None


def load_literature(key):
    with open(LIT, encoding="utf-8") as f:
        tab = json.load(f)
    if key not in tab:
        raise SystemExit(f"literature_J.json 里没有 {key!r}; 可用: {list(tab.keys())}")
    return tab[key]


def pick_structure(struct_dir, formula, mp_id=None):
    """从 candidates.csv 选最稳定 (hull 最小) 的单层 POSCAR。"""
    import csv
    csv_path = os.path.join(os.path.dirname(struct_dir), "candidates.csv")
    rows = list(csv.DictReader(open(csv_path, encoding="utf-8-sig")))
    if not rows:
        raise SystemExit(
            "candidates.csv 为空: MP 没有返回任何候选结构。\n"
            "  最常见原因: **当前 MP 数据库里没有这个化学式** (2026-10 实测: Fe3GeTe2 在 MP 中"
            "不存在, 连整个 Fe-Ge-Te 三元空间都是 0 条; 对照 CrI3 / Li-Fe-O 正常返回)。\n"
            "  办法: 用 --poscar 自带单层 POSCAR (含真空), 或用 --mp-id 指定 MP 里已有的条目。")
    if mp_id:
        rows = [r for r in rows if r["material_id"] == mp_id] or rows
    rows = sorted(rows, key=lambda r: (r.get("hull") or 9.0))
    best = rows[0]
    log(f"[pick] 选 {best['material_id']} ({best['formula']}) hull={best.get('hull')}")
    if len(rows) > 1:
        log(f"[pick] 共 {len(rows)} 个候选, 其余: "
            + ", ".join(r['material_id'] for r in rows[1:])
            + " (可用 --mp-id 指定)")
    return os.path.join(struct_dir, best["poscar"])


def _is_rectangular_net(lat, tol_cos=0.05, tol_ratio=0.05):
    """2D 面内晶格是否正交且 a≠b (矩形晶格, 如 CrSBr: a⊥b, |a|≠|b|)。

    原子数判据 (1→triangular / 2→honeycomb) 只看个数不看几何, 会把矩形单离子磁体
    误判成三角/蜂窝; 这里先验正交+不等长, 把矩形晶格救出来。
    """
    a = lat[0]
    b = lat[1]
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-8 or nb < 1e-8:
        return False
    cosang = abs(np.dot(a, b)) / (na * nb)
    ratio = abs(1.0 - na / nb)
    return cosang < tol_cos and ratio > tol_ratio


def _is_square_net(lat, tol_cos=0.05, tol_ratio=0.05):
    """2D 面内晶格是否正交且 a≈b (方形晶格 ab 面, 如 K2NiF4 / La2CuO4)。

    与 _is_rectangular_net 互补: 两者用同一个 tol_ratio 把 a≈b / a≠b 分开, 互斥。
    不加这一档的话, 方形单离子磁体 (1 磁性原子/胞) 会被下面的 n_mag == 1 判成
    triangular 而去解 120° 非共线三子晶格 —— 几何与磁序全错。
    """
    a = lat[0]
    b = lat[1]
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na < 1e-8 or nb < 1e-8:
        return False
    cosang = abs(np.dot(a, b)) / (na * nb)
    ratio = abs(1.0 - na / nb)
    return cosang < tol_cos and ratio <= tol_ratio


def _coordination_z(coords, mag, lat, tol=0.05, nshift=1):
    """返回 (d_min, z_list): 每个磁性位点在**最近磁-磁距离**上的真实配位数 z。

    必须显式枚举**周期像**: kagome 的原胞里每个位点在 a/2 上只有 2 个**胞内**伙伴
    (另两条键连的是同一批伙伴的周期像), 只数胞内原子会得到 z=2 —— 与"3 原子分子簇"
    无法区分。枚举 ±nshift 个平移后, kagome 得到 z=4 (角共享三角形), Fe3GeTe2 型的
    金属近邻网得到 z=6, Nb3 三聚体得到 z=2。z 是**与超胞无关**的纯几何量。
    """
    shifts = [np.array([i, j, k], dtype=float)
              for i in range(-nshift, nshift + 1)
              for j in range(-nshift, nshift + 1)
              for k in range(-nshift, nshift + 1)]
    n = len(mag)
    d_min = None
    for i in range(n):
        for j in range(n):
            for s in shifts:
                if i == j and not s.any():
                    continue
                d = coords[mag[j]] - coords[mag[i]] + s
                r = float(np.linalg.norm(d @ lat))
                d_min = r if d_min is None else min(d_min, r)
    zs = []
    for i in range(n):
        z = 0
        for j in range(n):
            for s in shifts:
                if i == j and not s.any():
                    continue
                d = coords[mag[j]] - coords[mag[i]] + s
                if abs(float(np.linalg.norm(d @ lat)) - d_min) <= 1e-6 + tol * d_min:
                    z += 1
        zs.append(z)
    return float(d_min), zs


def detect_lattice(poscar_path, element, cutoff=3.3):
    """读 reorient 后 POSCAR, 数磁性原子/胞 + 簇数, 返回 (lattice, n_clusters, n_mag)。

    lattice ∈ {"honeycomb", "triangular", "rectangular", "square", "rect_centered",
              "kagome", "sq_centered", "cluster", "other"}:
      cluster       = 磁性原子成**有限分子簇** (簇数 < 原子数, 且晶格不是正交网)
      kagome        = 3 个磁性原子/胞, 且每个位点在最近距离 a/2 上恰好 **4** 个近邻
                      (角共享三角形网; Fe3Sn 类)。判据用配位数而不是"六角 120°":
                      Fe3GeTe2 同样是 3 原子/胞 + 六角, 但 Fe 有 6 个同距离近邻 ⇒ 不是 kagome。
      rect_centered = 正交网 a≠b, 但磁性原子被 cutoff 并簇 ⇒ 多磁原子/胞的 C 心矩形网
                      (CrOCl 类 oxyhalide: 2 个 Cr/胞, 最近 Cr-Cr 3.035 A)
      sq_centered   = 同上但 a≈b (暂不支持, 显式报错而不是静默走方形路径)
      rectangular   = 正交且 a≠b 的矩形晶格 (CrSBr 类, 非受挫)
      square        = 正交且 a≈b 的方形晶格 (K2NiF4 / La2CuO4 类, 非受挫)
      triangular    = 1 磁性原子/胞
      honeycomb     = 2 磁性原子/胞 (不分簇)
    """
    lines = [l for l in open(poscar_path, encoding="utf-8").read().splitlines() if l.strip()]
    lat = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)])
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    idx = 7
    if lines[idx].strip().lower()[:1] == "s":
        idx += 1
    mode = lines[idx].strip().lower()
    idx += 1
    atoms = []
    for sp, n in zip(species, counts):
        atoms += [sp] * n
    coords = np.array([[float(x) for x in lines[idx + i].split()[:3]]
                       for i in range(len(atoms))])
    if mode[:1] in "ck":
        coords = coords @ np.linalg.inv(lat)
    coords -= np.floor(coords)
    mag = [i for i, a in enumerate(atoms) if a == element]
    n_mag = len(mag)
    if n_mag == 0:
        raise SystemExit(f"POSCAR 里找不到磁性元素 {element}")

    parent = list(range(n_mag))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n_mag):
        for j in range(i + 1, n_mag):
            d = coords[mag[j]] - coords[mag[i]]
            d -= np.floor(d + 0.5)
            if np.linalg.norm(d @ lat) < cutoff:
                ra, rb = find(i), find(j)
                if ra != rb:
                    parent[rb] = ra
    n_clusters = len({find(i) for i in range(n_mag)})

    # kagome 网 (角共享三角形): 3 个磁性原子/胞, 每人在 a/2 上恰好 4 个近邻。
    # 必须放在 `n_clusters < n_mag` **之前**: kagome 的最近邻 a/2 (~2.7 A) 远小于
    # 默认 cutoff=3.3 A ⇒ 三个原子会被并成 1 簇, 旧版一律路由到「分子簇磁体」
    # (Fe3GeTe2 也落在那里)。配位数 4 vs 6 是唯一可靠的几何判据。
    if n_mag == 3:
        _dmin, _zs = _coordination_z(coords, mag, lat)
        if all(z == 4 for z in _zs):
            return "kagome", n_clusters, n_mag

    if n_clusters < n_mag:
        # 磁性原子被 cutoff 并成了「簇」。这里其实有两种体系:
        #   (a) 真分子簇磁体 (Nb3I8 的 Nb3 三聚体): 簇是**有限分子**, 晶格本身不是正交网;
        #   (b) 「多磁原子/胞」的层内无限网 (CrOCl 单层: 2 个 Cr/胞, 最近 Cr-Cr=3.035 A < cutoff):
        #       两个原子被并成 1 簇, 但它是贯穿整层的网格, 不是分子。
        # 旧版一律判 cluster, 于是 CrOCl 会被送进 Nb3I8 的三聚体模型 —— 而且
        # build_rectangular_cases 的 tol=0.2 还会把 3.035 A 与 3.198 A 两个物理上
        # 不同的壳**静默合并**成一个 J1。所以这里按几何分开:
        # 晶格是正交网就说明是层内网格 (centered), 否则才交给簇路径。
        if _is_rectangular_net(lat):
            return "rect_centered", n_clusters, n_mag
        if _is_square_net(lat):
            return "sq_centered", n_clusters, n_mag
        return "cluster", n_clusters, n_mag
    if _is_rectangular_net(lat):
        return "rectangular", n_clusters, n_mag
    if _is_square_net(lat):
        return "square", n_clusters, n_mag
    if n_mag == 1:
        return "triangular", 1, 1
    if n_mag == 2:
        return "honeycomb", 2, 2
    return "other", n_clusters, n_mag


def multi_sublattice_warning(lattice, n_clusters, n_mag, element):
    """「分子簇磁体」与「多亚晶格铁磁网络」几何上无法区分 -> 返回警告文本或 None。

    Nb3I8 型分子簇 (Nb3 三聚体, 簇内强键合 / 簇间弱) 与 Fe3GeTe2 型巡游铁磁金属
    都有「磁性近邻 < cutoff」的特征: FGT 的 3d 金属层内 Fe-Fe ~2.6-3.0 A, 同样会被
    并成 1 簇, 于是整胞塌成 1 个有效磁位、套用「分子簇/单重态」模型 —— 对铁磁金属
    没有物理意义 (簇内 J 不是可测量, T_N 也不可信)。
    这里**只警告不改路由**: 体系是不是分子簇是化学事实, 几何判据判不出来。
    """
    if lattice == "cluster" and n_clusters == 1 and n_mag >= 3:
        return (f"{element} {n_mag} 个/胞 被并成 1 簇 → 走「分子簇磁体」模型。"
                f"注意几何判据**分不出**下面两种体系 (签名会完全相同: Nb3I8 与 Fe3GeTe2 "
                f"都是 cluster / 1 簇 / 3 个磁性原子): "
                f"(a) 真分子簇 (Nb3I8 型: 簇内强键合、簇间弱) → 本分支适用; "
                f"(b) 多亚晶格铁磁网络 (Fe3GeTe2/Fe3GaTe2 型: 连通的金属网络) → 本分支**不适用**: "
                f"超胞里会连成一整个簇, build_cluster_cases.py 会直接报错退出, "
                f"簇内 J 无物理意义、T_N 不可信, 需另建双亚晶格模型 (见 SKILL.md「多亚晶格铁磁网络」)")
    return None


def _print_summary(outdir, lattice, S, J1, J2, J3, MAE, mc_order, mc_refine, args):
    """打印 MC/FSS 结果摘要 (供 mc 子命令与全流水线共用, 避免重复)。"""
    print(f"      lattice={lattice} | J1={J1} J2={J2} J3={J3} MAE={MAE} meV | S={S}")
    s = summarize_mc(outdir, mc_order)
    if s:
        if mc_refine and mc_refine["rounds"]:
            s += f" [峰落边界, 已自动缩放重扫 {len(mc_refine['rounds'])} 轮]"
        print(f"      MC: {s}")
    if args.fss:
        _fj = os.path.join(outdir, "05_mc", "fss_fss.json")
        if os.path.exists(_fj):
            try:
                _r = json.load(open(_fj, encoding="utf-8"))
                if _r.get("tc") is not None:
                    print(f"      FSS: Binder U4 交点 T_C/T_N ≈ {_r['tc']:.2f} K ± {_r['sd']:.2f} K "
                          f"(尺寸 {_r['sizes']}, {_r['n_crossings']} 个交点)"
                          + (f"; collapse {_r['tc_collapse']:.2f} K" if _r.get("tc_collapse") else ""))
                else:
                    print("      FSS: 未分辨 (无干净交点, 见 05_mc/fss_fss.csv/.png)")
            except Exception as e:
                print(f"      FSS: 读取结果失败 ({e})", file=sys.stderr)
        else:
            print("      FSS: 未产生结果文件 (见 05_mc/fss_fss.csv / .png)", file=sys.stderr)


def has_cant_orders(cases):
    """cases/magnetic_orders.json 里是否有倾斜家族 (cant_orders 非空)。

    solve 端据此在 crect 路径上改派 solve_J_crect_cant.py —— 不去猜「用户是不是想反解倾斜」,
    只看落盘元数据; 读不到/没有该键一律当 False (老产物向后兼容)。
    """
    path = os.path.join(cases, "magnetic_orders.json")
    if not os.path.isfile(path):
        return False
    try:
        with open(path, encoding="utf-8") as f:
            meta = json.load(f)
    except (OSError, ValueError):
        return False
    return bool(meta.get("cant_orders"))


def main():
    ap = argparse.ArgumentParser(description="二维磁体 结构->J->MC 通用流水线")
    # 子命令 (缺省=全流水线; 只跑单个阶段, 复用已有产物)
    ap.add_argument("command", nargs="?", default=None,
                    choices=["fetch", "build", "dft", "ingest", "kconv", "mc", "solve", "render", "report", "doctor"],
                    help="只跑单个阶段 (fetch/build/dft/ingest/kconv/mc/solve/render/report); "
                         "kconv=P4.2 k 点收敛抽查 (生成 KPOINTS 变体 + 回收判定); "
                         "doctor=环境自检+error catalog; 缺省跑全流水线")
    # 结构来源
    ap.add_argument("--formula", help="化学式, 如 MnPS3 (触发 MP 下载)")
    ap.add_argument("--mp-id", help="指定 MP material_id (覆盖自动挑选)")
    ap.add_argument("--poscar", help="已有单层 POSCAR (跳过下载)")
    ap.add_argument("--mp-key", default=None, help="MP_API_KEY (默认读环境变量)")
    ap.add_argument("--provider", choices=["auto", "mp", "cod", "oqmd", "aflow", "jarvis"],
                    default="auto", help="结构来源 provider (auto=MP 优先, 失败回退 cod→oqmd→aflow→jarvis)")
    # 材料/方法参数
    ap.add_argument("--element", help="磁性元素, 如 Mn")
    ap.add_argument("--nonmag", nargs="+", help="非磁元素, 如 P S")
    ap.add_argument("--S", type=float, help="自旋大小 (Mn2+=2.5, Cr3+=1.5)")
    ap.add_argument("--moment", type=float, help="磁矩初猜 muB (Mn2+=5.0, Cr3+=3.0)")
    ap.add_argument("--U", type=float, help="DFT+U 加在磁性元素 d 轨道 (eV)")
    # J 值
    ap.add_argument("--model", type=str.lower, choices=sorted(MODEL_PRESETS), default=None,
                    help="命名模型预设: 一键设定 formula/element/nonmag/lattice/S/moment/U "
                         "(关联文献 J 仅作参考)。可选: "
                         + "; ".join(f"{k}={v['desc']}" for k, v in sorted(MODEL_PRESETS.items())))
    ap.add_argument("--literature", help="literature_J.json 里的键, 如 MnPS3")
    ap.add_argument("--J1", type=float, default=None)
    ap.add_argument("--J2", type=float, default=None)
    ap.add_argument("--J3", type=float, default=None)
    ap.add_argument("--MAE", type=float, default=None)
    # 输出 / 开关
    ap.add_argument("--outdir", default=None, help="项目目录 (默认 ./<formula>_pipeline)")
    ap.add_argument("--vacuum", type=float, default=15.0)
    ap.add_argument("--supercell", type=int, nargs=3, default=[2, 2, 1])
    ap.add_argument("--lattice",
                    choices=["auto", "honeycomb", "triangular", "cluster", "rectangular", "square",
                             "rect_centered", "kagome"],
                    default="auto",
                    help="磁性格点类型 (默认自动判别; rectangular=正交 a≠b 且 1 磁位/胞; "
                         "square=正交 a≈b 且 1 磁位/胞; rect_centered=正交 a≠b 且多磁位/胞的 C 心矩形网, 如 CrOCl; "
                         "kagome=3 磁位/胞的角共享三角形网, 如 Fe3Sn)")
    ap.add_argument("--cluster-cutoff", type=float, default=3.3,
                    help="簇识别成键阈值 A (簇磁体用)")
    ap.add_argument("--cluster-nbr-shells", type=int, default=3,
                    help="簇间壳数 (1=只 J1; 3=J1/J2/J3, 簇磁体用)")
    ap.add_argument("--cluster-supercell", type=int, nargs=3, default=[5, 5, 1],
                    help="簇磁体能量映射超胞 (默认 5 5 1)")
    ap.add_argument("--tri-supercell", type=int, nargs=3, default=[5, 5, 1],
                    help="三角磁体能量映射超胞 (默认 5 5 1, J3=2a 需 smax>=2)")
    ap.add_argument("--rect-supercell", type=int, nargs=3, default=[4, 4, 1],
                    help="矩形磁体能量映射超胞 (默认 4 4 1, 偶数保证条纹态纯净)")
    ap.add_argument("--square-supercell", type=int, nargs=3, default=[4, 4, 1],
                    help="方形磁体能量映射超胞 (默认 4 4 1; a/b 必须相等且是 4 的倍数——"
                         "打破 J3/E0 简并的 p4 序要求周期 4)")
    ap.add_argument("--crect-supercell", type=int, nargs=3, default=[4, 4, 1],
                    help="C 心矩形磁体能量映射超胞 (默认 4 4 1; n0/n1 必须是偶数——"
                         "afm_a/afm_b 沿 a/b 的 (-1)^n 条纹要求偶数才能纯净)")
    ap.add_argument("--kagome-supercell", type=int, nargs=3, default=[3, 3, 1],
                    help="kagome 能量映射超胞 (默认 3 3 1; n0/n1 必须是 3 的倍数——"
                         "√3×√3 序的定义相位要求周期 3, 少了它 J1/J2 在 q=0 简并下不可分离)")
    ap.add_argument("--mc-size", type=int, nargs=2, default=[30, 30])
    ap.add_argument("--mc-tmin", type=float, default=5.0)
    ap.add_argument("--mc-tmax", type=float, default=150.0)
    ap.add_argument("--mc-nt", type=int, default=30)
    ap.add_argument("--mc-refine-nt", type=int, default=60,
                    help="峰判据落在窗口边界时, 缩放重扫的格点数 (默认 60)")
    ap.add_argument("--mc-refine-max", type=int, default=3,
                    help="峰判据落在窗口边界时, 最多自动缩放重扫轮数 (默认 3)")
    ap.add_argument("--no-mc-refine", action="store_true",
                    help="关闭「峰落边界自动缩放重扫」(只声明未分辨, 不重跑)")
    ap.add_argument("--mc-seeds", type=int, default=1,
                    help="MC 独立种子数 (误差棒): 峰位 ± 种子间标准差 (默认 1 = 关闭)")
    ap.add_argument("--mc-order", default=None,
                    help="强制 MC 序参量 (默认按 lattice 自动选): auto | staggered | total | "
                         "cv | sublattice:<mask>。交变磁体 (P5) 用 sublattice:<mask> 指定 A 组"
                         "子晶格号, 例: --mc-order sublattice:0 或 sublattice:0,1")
    ap.add_argument("--fss", action="store_true",
                    help="额外跑多尺寸有限尺寸标度 (Binder U4 交点自动定 T_C/T_N)")
    ap.add_argument("--fss-sizes", type=int, nargs="+", default=[12, 18, 24, 30],
                    help="FSS 的线性尺寸列表 (默认 12 18 24 30)")
    # ingest (P4.2 收回通道): 只对 ingest 子命令生效
    ap.add_argument("--runs", default=None,
                    help="ingest: OUTCAR 回收目录 (默认 <outdir>/runs; 也可指向从集群拷回的任意目录)")
    ap.add_argument("--ingest-dry-run", action="store_true",
                    help="ingest: 只校验不落盘 (先看闸门过不过)")
    ap.add_argument("--ingest-force", action="store_true",
                    help="ingest: 允许覆盖 cases/<order>/OUTCAR 里已有的旧文件")
    ap.add_argument("--ingest-min-kb", type=float, default=None,
                    help="ingest: OUTCAR 最小体积 KB (默认按 ingest_dft 的判据; 明显偏小=被 kill 的空壳)")
    ap.add_argument("--ingest-max-spread", type=float, default=None,
                    help="ingest: 四态能量口径最大允许离散 eV (默认 ingest_dft 的判据)")
    ap.add_argument("--ingest-min-gap", type=float, default=None,
                    help="ingest: 四态之间要求的最小能量差 eV (防四个序其实是同一个构型)")
    ap.add_argument("--ingest-lenient", nargs="+", default=None, metavar="CHECK",
                    help="ingest: 把指定硬闸门降级为警告 (如 ionic_steps); 慎用, 会削弱 J 的可信度")
    ap.add_argument("--skip-dft", action="store_true")
    ap.add_argument("--skip-mc", action="store_true")
    ap.add_argument("--noncollinear", action="store_true",
                    help="DFT 输入用非共线 (LNONCOLLINEAR + 每原子 3 分量 MAGMOM); "
                         "三角晶格本来就默认非共线")
    ap.add_argument("--soc", action="store_true",
                    help="非共线 + SOC (LSORBIT + GGA_COMPAT=.FALSE.); 弱铁磁倾斜/DMI/kagome 必需")
    ap.add_argument("--saxis", nargs=3, type=float, default=None,
                    help="自旋量子化轴 (SAXIS), 如 --saxis 0 0 1")
    # ---- P4.2: VASP 输入的可调项（默认与 03_dft 模板一致，只有显式给了才透传）----
    ap.add_argument("--kpoints", type=int, nargs=3, default=None,
                    help="显式 k 网格 (默认用模板的 3 3 1; 上真实机时前建议先做一次 k 收敛抽查, "
                         "见 P4_2_KAGOME_RUNBOOK.md §1)")
    ap.add_argument("--npar", type=int, default=None, help="显式 NPAR (默认 4)")
    ap.add_argument("--isym", type=int, default=None,
                    help="显式 ISYM (0 = 关对称性); 默认 static 不写、relax 用模板的 0")
    # ---- P4.2: k 点收敛抽查 (子命令 kconv) ----
    ap.add_argument("--k-grids", default="3,4,6,8",
                    help="kconv: k 网格列表, 如 '3,4,6,8' (面内 k, kz=1) 或 '3 3 1;4 4 1'")
    ap.add_argument("--k-tol", type=float, default=1.0,
                    help="kconv: 判据 |ΔE| <= 该值 (meV/原子, 默认 1.0)")
    ap.add_argument("--k-orders", nargs="+", default=None, help="kconv: 只抽查这些序 (默认全部)")
    ap.add_argument("--k-plan-only", action="store_true", help="kconv: 只生成网格变体, 不判定")
    ap.add_argument("--k-ingest-only", action="store_true", help="kconv: 只判定 (不重新生成)")
    # ---- P2′: 弱铁磁倾斜 (+θ/−θ 家族) ----
    ap.add_argument("--cant-theta", default=None, metavar="\"8 16\"",
                    help="(只对 C 心矩形 / CrOCl 类) 额外生成 ±θ 弱铁磁倾斜家族。角度列表用空格或逗号分隔, "
                         "要求 |θ| < 90; 构型写 3 分量 MAGMOM ⇒ 装配需要 --soc; 反解自动改派 "
                         "solve_J_crect_cant.py (P2′)。θ 是**输入参数**, 不是自洽解出来的基态角")
    ap.add_argument("--cant-base", default="afm_c",
                    help="倾斜的共线母序 (afm_c | afm_a | afm_b; 默认 afm_c; fm 被拒: 全同向母序上叠 ±θ "
                         "只是一次整体自旋旋转, 系数行恒等于 fm 行、不带信息)")
    ap.add_argument("--cant-azimuth", type=float, default=0.0,
                    help="倾斜方位角 φ (度, 落在 [0,360)); 注意 φ **由几何定不出来** ⇒ 这是取值不是解")
    ap.add_argument("--quiet", action="store_true",
                    help="静默模式: 抑制子进程输出, 只回显错误/警告 + 结尾摘要")
    ap.add_argument("--version", action="version", version=f"2d-magnet {pkg_version()}",
                    help="打印版本号 (取自 pyproject.toml) 并退出")
    args = ap.parse_args()

    if args.command == "doctor":
        sys.exit(0 if doctor() else 1)

    global QUIET
    QUIET = args.quiet

    # ---- 解析参数 (文献表/模型预设 < 显式 CLI) ----
    model = None
    if args.model:
        model = MODEL_PRESETS[args.model]
        log(f"[model] 命名预设 {args.model!r}: {model['desc']}")
        # 预设自带文献 J 作参考; 显式 --literature 可覆盖
        if not args.literature:
            args.literature = model.get("literature")

    lit = None
    if args.literature:
        lit = load_literature(args.literature)
        log(f"[lit] 使用文献值 ({args.literature}): {lit.get('source')}")

    def pick(a, b, d=None):
        """解析优先级: 显式 CLI > --model 预设 > 文献表非 null 值 > 默认 d。

        注意: 文献条目里的 `null` 表示「未核实/未给出」, 必须落回默认值 —— 不能把它
        当字符串传下去 (FGT 条就只核实了 J1, J2/J3 是 null; 旧写法会让 J2 变成 'None'
        再喂给 mc_2d.py)。J1 的默认是 None: 在需要 J 的阶段 (mc/report/全流水线) 会触发
        「缺少必要参数」, 但在只用几何的阶段 (build/dft/solve/render) 允许缺省 —— 这样文献 J
        未核实的材料仍能进主链做结构/磁序枚举, 不必编一个 J 出来。
        """
        v = a if a is not None else (model.get(b) if model else None)
        if v is None:
            v = lit.get(b) if lit else None
        return d if v is None else v

    formula = args.formula or (model.get("formula") if model else None) \
        or (lit.get("formula") if lit else None)
    element = pick(args.element, "element")
    nonmag = args.nonmag or (model.get("nonmag") if model else None) \
        or (lit.get("nonmag") if lit else [])
    S = pick(args.S, "S")
    moment = pick(args.moment, "moment")
    U = pick(args.U, "U")
    J1 = pick(args.J1, "J1")
    J2 = pick(args.J2, "J2", 0.0)
    J3 = pick(args.J3, "J3", 0.0)
    MAE = pick(args.MAE, "MAE", 0.0)

    # J 落回文献条目时明确声明「仅参考」, 免得把文献值当成本项目结果。
    # 只对「确实取自文献非 null 条目」的 J 报警; 文献 null -> 落默认 0.0 的不算 (那是默认, 不是文献)。
    if lit:
        ref_J = [k for k, given in (("J1", args.J1), ("J2", args.J2), ("J3", args.J3))
                 if given is None and lit.get(k) is not None]
        if ref_J:
            print(f"[ref] {', '.join(ref_J)} 来自文献参考 ({args.literature}), "
                  f"非本项目 DFT 计算值; 请跑 VASP (03_dft/runs 导出后手动提交) + `solve` 反解真实 J, "
                  f"或用 --J1/--J2/--J3 覆盖", file=sys.stderr)

    # 元素参数表兜底: 显式 CLI > 文献表 > 元素表 (ELEMENT_PARAMS)。落表时打 stderr 提醒,
    # 免得静默给出默认物理参数 (S/moment/U 都可被 --S/--moment/--U 覆盖)。
    ep = ELEMENT_PARAMS.get(element) if element else None
    if ep:
        filled = []
        if S is None and "S" in ep:
            S = ep["S"]; filled.append("S")
        if moment is None and "moment" in ep:
            moment = ep["moment"]; filled.append("moment")
        if U is None and "U" in ep:
            U = ep["U"]; filled.append("U")
        if filled:
            print(f"[params] {element} 的 {', '.join(filled)} 未显式给出 -> 元素表默认 "
                  f"S={S}, moment={moment}, U={U} (用 --S/--moment/--U 覆盖)",
                  file=sys.stderr)
    elif element and any(v is None for v in (S, moment, U)):
        # 不在元素表里 (Fe 是刻意的): 这里**不猜值**, 直接说清原因与下一步 —— 否则用户看到的是
        # 后面「--S/--moment/--U 缺失」的症状, 而不是「本表故意不列该元素」的原因。
        _miss = ", ".join(k for k, v in (("--S", S), ("--moment", moment), ("--U", U)) if v is None)
        print(f"[params] {element} 不在内置元素参数表 (Mn / Cr / Ni / Co / V / Gd) 里 —— 该表**故意**不列 "
              f"Fe 等变价/巡游元素的有效自旋 (猜一个值会静默给出错的物理)。请显式给 {_miss}; "
              f"或用 --model <预设> / --literature <条目> 取材料专属值。", file=sys.stderr)

    # J 只在**真正需要它**的阶段才必须: mc / report / 全流水线。
    # build/dft/solve/render 只用几何与规格层 —— 这一点很重要: 文献 J 未核实的材料
    # (K2NiF4 / CrOCl / CrI3 / CrBr3 / GdI2) 的**结构与磁序枚举本来就该在没有 J 的情况下可用**,
    # 旧版一律硬要 --J1, 等于逼用户编一个 J 出来 (与「不猜数」的 J 政策自相矛盾)。
    need_j = args.command in (None, "mc", "report")
    need = [("--element", element), ("--S", S), ("--moment", moment), ("--U", U)]
    if need_j:
        need.append(("--J1", J1))
    missing = [k for k, v in need if v is None]
    if missing:
        raise SystemExit(f"缺少必要参数 {missing}; 用 --literature <key> 或显式给出")
    if J1 is None:
        # 走 stdout 而不是 stderr: 这是正常路径的说明, 不是错误 (写到 stderr 会让调用方
        # /PowerShell 把它当 NativeCommandError, 误判成失败)。
        print(f"[params] J1 未给出 (该材料文献 J 未核实) -> 本阶段 {args.command!r} 只用几何/规格层; "
              f"跑 mc/report 时必须显式给 --J1/--J2/--J3")

    if not formula and not args.poscar and not args.mp_id:
        raise SystemExit("需要 --formula / --mp-id / --poscar 之一")

    if args.outdir:
        outdir = os.path.abspath(args.outdir)
    else:
        outdir = os.path.abspath((formula or "material") + "_pipeline")
    os.makedirs(outdir, exist_ok=True)
    log(f"[run] 项目目录: {outdir}")

    cmd = args.command
    # ---- P4.2: k 点收敛抽查 (只读已有 runs/, 生成 KPOINTS 变体 + 回收后判定) ----
    if cmd == "kconv":
        _a = ["--root", outdir, "--grids", args.k_grids, "--tol", str(args.k_tol)]
        if args.k_orders:
            _a += ["--orders"] + list(args.k_orders)
        if args.k_plan_only:
            _a += ["--plan-only"]
        if args.k_ingest_only:
            _a += ["--ingest-only"]
        run("k_convergence.py", _a, cwd=outdir)
        return
    # 子命令语义: 只跑对应阶段(及其数据依赖), 复用已有产物; 各阶段边界处早退
    run_dft = cmd in (None, "dft")
    run_mc = cmd in (None, "mc")
    run_render = cmd in (None, "render")

    # ---- 1. fetch ----
    os.makedirs(os.path.join(outdir, "01_structure"), exist_ok=True)
    poscar_mp = os.path.join(outdir, "01_structure", "POSCAR.mp")
    if args.poscar:
        src_poscar = os.path.abspath(args.poscar)
        log(f"[fetch] 跳过下载, 用 {src_poscar}")
        shutil.copy(src_poscar, poscar_mp)
    elif cmd is not None and os.path.exists(poscar_mp):
        log("[fetch] 子命令复用已下载的 01_structure/POSCAR.mp")
    else:
        if args.mp_key:
            os.environ["MP_API_KEY"] = args.mp_key
        provider = args.provider
        got = False

        # 1) Materials Project (需 key; 会剥单层)。--provider mp/auto 时先试。
        if provider in ("auto", "mp"):
            if "MP_API_KEY" in os.environ:
                src_poscar = try_mp_fetch(outdir, formula, args.mp_id)
                if src_poscar:
                    shutil.copy(src_poscar, poscar_mp)
                    got = True
                    log("[fetch] Materials Project 成功")
                else:
                    log("[fetch] MP 失败, 回退其他 provider")
            else:
                log("[fetch] 无 MP_API_KEY, 跳过 MP (用 --mp-key, 或 --provider cod/oqmd/aflow/jarvis)")
            if provider == "mp" and not got:
                raise SystemExit("MP provider 失败 (需 --mp-key 且 MP 里存在该结构)")

        # 2) 其他 provider 链 (COD/OQMD/AFLOW/JARVIS, 无 key)。
        if not got:
            chain = ["cod", "oqmd", "aflow", "jarvis"] if provider == "auto" else [provider]
            for p in chain:
                if p == "mp":
                    continue
                log(f"[fetch] 试 {p} ...")
                r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fetch_structure.py"),
                                    formula, "--provider", p, "--out", poscar_mp],
                                   capture_output=True, text=True, encoding="utf-8",
                                   errors="replace", env=child_env())
                if r.returncode == 0 and os.path.exists(poscar_mp):
                    log(f"[fetch] {p} 成功")
                    got = True
                    break
                print(f"[fetch] {p} 失败: {(r.stderr or r.stdout or '').strip()[-200:]}", file=sys.stderr)
            if not got:
                raise SystemExit("所有结构 provider 均失败; 用 --poscar 自带单层 POSCAR, "
                                 "或 --provider mp 配 --mp-key")

    # ---- 2. reorient ----
    # POSCAR 第 1 行是注释 (VASP 忽略), 但它是产物的一部分、会被人读到, 所以要带对材料名。
    _preset = MODEL_PRESETS.get(args.model) if args.model else None
    mat_label = (_preset or {}).get("formula") if isinstance(_preset, dict) else None
    if not mat_label:
        mat_label = os.path.basename(os.path.abspath(outdir)).replace("_pipeline", "")
    run("reorient.py", [os.path.join(outdir, "01_structure", "POSCAR.mp"),
                        "--vacuum", str(args.vacuum),
                        "--comment", f"{mat_label} monolayer (reoriented, vacuum)",
                        "-o", os.path.join(outdir, "01_structure", "POSCAR")], cwd=outdir)
    poscar_re = os.path.join(outdir, "01_structure", "POSCAR")

    if cmd == "fetch":
        print(f"\n[run] fetch 完成: {poscar_re} (reorient 后, 真空 {args.vacuum} Å)")
        return

    # ---- 2.5 检测磁性格点类型 ----
    lattice, n_clusters, n_mag = detect_lattice(poscar_re, element, args.cluster_cutoff)
    detected = lattice
    if args.lattice != "auto":
        lattice = args.lattice
    elif model and model.get("lattice"):
        lattice = model["lattice"]
        if lattice != detected:
            log(f"[model] 预设晶格类型 {lattice} (自动判别是 {detected})")
    log(f"[lattice] {lattice} (磁性原子 {n_mag}/胞, 簇 {n_clusters})")
    if lattice == "sq_centered":
        raise SystemExit(
            "检测到「方形 + 多磁原子/胞」(sq_centered): 磁性原子被 cutoff 并簇, 但晶格 a≈b。\n"
            "  这与 C 心矩形 (rect_centered) 是两类不同的网, 本仓库目前**只支持 rect_centered**。\n"
            "  静默套用方形路径 (1 磁位/胞) 会把两个 basis 的反平行自由度丢掉、给出错的 J ——\n"
            "  所以这里显式拒绝而不是降级。要支持它需要: 方形版的 basis 序 + 显式壳距 + 独立的反解脚本。")
    lat_warn = multi_sublattice_warning(lattice, n_clusters, n_mag, element)
    if lat_warn:
        # 走 stderr: --quiet 下 log() 被抑制, 但这条必须让用户看见
        print(f"[warn] {lat_warn}", file=sys.stderr)

    # ---- P2′: 倾斜家族的参数守卫 (不静默忽略) ----
    if args.cant_theta:
        if lattice != "rect_centered":
            raise SystemExit(
                f"给了 --cant-theta 但 lattice={lattice}, 不是 C 心矩形 (rect_centered)。\n"
                "  倾斜家族只在 C 心矩形网 (CrOCl 类, 2 磁原子/胞) 上定义: 叠在母序上的 ±θ 是**胞内两个\n"
                "  basis** 的反平行 z 分量 + 同向面内分量。静默忽略这个参数会让你以为跑了倾斜、实际只\n"
                "  跑了共线四序 —— 所以这里显式拒绝。要用请加 --lattice rect_centered (或 --model crocl)。")
        if args.cant_base == "fm":
            raise SystemExit(
                "给了 --cant-theta 却用 --cant-base fm: 全同向母序上叠 ±θ 只是一次**整体**自旋旋转,\n"
                "  任意键的点积恒为 +1 ⇒ 倾斜行的系数与 fm 行逐位相同 (θ 标度自检还会假通过)。\n"
                "  母序取 afm_c (默认) / afm_a / afm_b。")
        if not args.soc:
            # 不硬拒: build/solve 两个阶段不需要 SOC, 只有 dft 装配需要。但必须让用户看见,
            # 否则 dft 阶段会被 assemble_dft 的 SOC 门禁拒掉 (那才是硬错)。
            print("[warn] --cant-theta 给了但没给 --soc: 倾斜构型装配需要 LSORBIT, "
                  "dft 阶段会被 assemble_dft 的 SOC 门禁硬拒 (不会静默按共线算)", file=sys.stderr)

    # ---- solve 子命令: 对已算完的 OUTCAR 反解 J (需先跑过 dft 并手动跑 VASP) ----
    if cmd == "solve":
        cases = os.path.join(outdir, "02_supercell", "cases")
        if not os.path.exists(os.path.join(cases, "magnetic_orders.json")):
            raise SystemExit("solve 需要 02_supercell/cases (先跑 build 或全流水线, 再跑 VASP 填 OUTCAR)")
        if lattice == "cluster":
            run("solve_cluster_J.py", [cases], cwd=outdir)
        elif lattice == "triangular":
            run("solve_J_tri.py", [cases], cwd=outdir)
        elif lattice == "rectangular":
            run("solve_J_rect.py", [cases], cwd=outdir)
        elif lattice == "square":
            run("solve_J_square.py", [cases], cwd=outdir)
        elif lattice == "rect_centered":
            if has_cant_orders(cases):
                log("[solve] 检测到倾斜家族 (cant_orders 非空) ⇒ 改派 solve_J_crect_cant.py "
                    "(--orders cant: 倾斜行 + 张成所需的 fm/母序行)")
                run("solve_J_crect_cant.py", [cases, "--orders", "cant"], cwd=outdir)
            else:
                run("solve_J_crect.py", [cases], cwd=outdir)
        elif lattice == "kagome":
            run("solve_J_kagome.py", [cases], cwd=outdir)
        else:
            run("solve_J3.py", [cases], cwd=outdir)
        print(f"\n[run] solve 完成 (lattice={lattice}, cases={cases})")
        return

    # ---- ingest 子命令: 把集群跑完的 OUTCAR 收回 cases/ (P4.2, 唯一收回通道) ----
    # 这一步以前是手工 cp, 是整条链上最危险的一步 (未收敛/混口径/错序的 OUTCAR 会被
    # 静默当成好数据)。现在走 ingest_dft.py 的硬闸门, 不通过就一件都不写。
    if cmd == "ingest":
        cases = os.path.join(outdir, "02_supercell", "cases")
        if not os.path.exists(os.path.join(cases, "magnetic_orders.json")):
            raise SystemExit("ingest 需要 02_supercell/cases (先跑 build 或全流水线)")
        runs_dir = args.runs or os.path.join(outdir, "runs")
        if not os.path.isdir(runs_dir):
            raise SystemExit(f"ingest 找不到回收目录: {runs_dir} "
                             f"(用 --runs 指定, 或先把集群产物放回那里)")
        ig_args = [runs_dir, cases]
        for flag, val in (("--min-kb", args.ingest_min_kb),
                          ("--max-spread", args.ingest_max_spread),
                          ("--min-energy-gap", args.ingest_min_gap)):
            if val is not None:
                ig_args += [flag, str(val)]
        if args.ingest_dry_run:
            ig_args.append("--dry-run")
        if args.ingest_force:
            ig_args.append("--force")
        for name in (args.ingest_lenient or []):
            ig_args += ["--lenient", name]
        log(f"[ingest] 收回通道: {runs_dir} -> {cases}")
        run("ingest_dft.py", ig_args, cwd=outdir)
        print(f"\n[run] ingest 完成 (cases={cases}); 接着跑: "
              f"python run_2d_magnet.py solve --outdir {outdir}")
        return

    mc_poscar = poscar_re
    mc_order = "auto"
    do_dft = False
    do_render = False
    render_kind = None

    if lattice == "cluster":
        # 簇磁体: 磁序枚举 + 能量映射 (J_intra/J1/J2/J3) -> 02_supercell/cases + 03_dft/runs
        run("build_cluster_cases.py",
            [poscar_re, "--element", element, "--nonmag"] + nonmag +
            ["--cutoff", str(args.cluster_cutoff),
             "--nbr-shells", str(args.cluster_nbr_shells),
             "--supercell"] + [str(x) for x in args.cluster_supercell] +
            ["--outdir", os.path.join(outdir, "02_supercell", "cases"),
             "--root", outdir, "--U", str(U), "--moment", str(moment)], cwd=outdir)
        # 磁性簇 -> 有效位 (供 MC)
        poscar_clusters = os.path.join(outdir, "01_structure", "POSCAR.clusters.vasp")
        run("cluster_sites.py", [poscar_re, "--element", element,
                                 "--cutoff", str(args.cluster_cutoff),
                                 "-o", poscar_clusters], cwd=outdir)
        mc_poscar = poscar_clusters
        mc_order = "cv"
        do_render = True
        render_kind = "cluster"
    elif lattice == "triangular":
        # 三角单离子: 磁序枚举 (含 120° 非共线) + J1/J2/J3 映射 -> 02_supercell/cases + 03_dft/runs
        run("build_magnetic_cases_tri.py",
            [poscar_re, "--supercell"] + [str(x) for x in args.tri_supercell] +
            ["--outdir", os.path.join(outdir, "02_supercell", "cases"),
             "--element", element, "--nonmag"] + nonmag +
            ["--moment", str(moment), "--spin", str(S)], cwd=outdir)
        mc_order = "cv"
        do_dft = True
        do_render = True
        render_kind = "triangular"
    elif lattice == "rectangular":
        # 矩形单离子磁体 (正交 a≠b, 如 CrSBr): 四共线磁态 + J1/J2/J3 映射
        run("build_rectangular_cases.py",
            [poscar_re, "--supercell"] + [str(x) for x in args.rect_supercell] +
            ["--outdir", os.path.join(outdir, "02_supercell", "cases"),
             "--element", element, "--nonmag"] + nonmag +
            ["--moment", str(moment), "--spin", str(S)], cwd=outdir)
        mc_poscar = poscar_re
        mc_order = "auto"
        do_dft = True
        do_render = True
        render_kind = "rectangular"
    elif lattice == "square":
        # 方形单离子磁体 (正交 a≈b, 如 K2NiF4 / La2CuO4): 四共线磁态 + J1/J2/J3 映射。
        # 必须含周期 4 的 p4 序, 否则 J3 与 E0 不可分离 (见 build_square_cases.py docstring)。
        run("build_square_cases.py",
            [poscar_re, "--supercell"] + [str(x) for x in args.square_supercell] +
            ["--outdir", os.path.join(outdir, "02_supercell", "cases"),
             "--element", element, "--nonmag"] + nonmag +
            ["--moment", str(moment), "--spin", str(S)], cwd=outdir)
        mc_poscar = poscar_re
        mc_order = "auto"
        do_dft = True
        do_render = True
        render_kind = "square"
    elif lattice == "rect_centered":
        # C 心矩形多亚晶格 (正交 a≠b, 2 磁原子/胞, 如 CrOCl): 四共线磁态 fm/afm_c/afm_a/afm_b。
        # afm_c (胞内两个 basis 反平行) 是唯一能把 J1 —— 即 C 心矢量 3.0347 A, 连的是
        # **同一胞内的两个 Cr** —— 从 E0 里分出来的序; 少了它 J1 与 E0 在设计矩阵里合并。
        # 壳距离必须**显式声明**: J1 与 J2 只差 0.163 A, 广义 magnetic_shells(tol=0.2)
        # 会把两个物理上不同的壳静默合并 (反解照样收敛、残差照样是 0)。
        run("build_crect_cases.py",
            [poscar_re, "--supercell"] + [str(x) for x in args.crect_supercell] +
            ["--outdir", os.path.join(outdir, "02_supercell", "cases"),
             "--element", element, "--nonmag"] + nonmag +
            ["--moment", str(moment), "--spin", str(S)] +
            (["--cant-theta", args.cant_theta, "--cant-base", args.cant_base,
              "--cant-azimuth", str(args.cant_azimuth)] if args.cant_theta else []), cwd=outdir)
        if args.cant_theta:
            log(f"[cant] 倾斜家族: θ = {args.cant_theta.strip()} (母序 {args.cant_base}, φ = "
                f"{args.cant_azimuth}°); θ 是输入参数 —— 基态 θ 需 SOC 总能量极小, 本流水线不算")
        mc_poscar = poscar_re
        mc_order = "auto"
        do_dft = True
        do_render = True
        render_kind = "rect_centered"
    elif lattice == "kagome":
        # kagome 角共享三角形网 (六角 a=b, 3 磁原子/胞, 如 Fe3Sn): 七个磁态 + 三个壳。
        # 壳 1: J1 = a/2, 角共享三角形内的近邻 —— 但**两条键不等价** (J1 与 J1', 同距离
        #       不同亚晶格对), 本 3 壳模型取两者的**平均** (要分开需降对称超胞, P4.2);
        # 壳 2: J2 = a√3/2 跨六角环; 壳 3: J3 = a 同亚晶格对角。
        # ⚠ q=0 一族 (fm/q0_120/q0_ud/q0_duu) 对 J1/J2 **简并** (rank=2, cond≈2.3e34):
        #   只用 q=0 序反解会得到一个"看似正常"的 J1/J2 —— 其实只是 (J1+J2)/2 的任意分配。
        #   builder 的 selected = (fm, q0_120, sq3_updntri, sq3_alt) 含 √3×√3 序 ⇒ rank=4
        #   (cond=98.51), 并把 rank/cond 落盘供 solve 端复核; solve 端 rank<4 直接报错。
        run("build_kagome_cases.py",
            [poscar_re, "--supercell"] + [str(x) for x in args.kagome_supercell] +
            ["--outdir", os.path.join(outdir, "02_supercell", "cases"),
             "--element", element, "--nonmag"] + nonmag +
            ["--moment", str(moment), "--spin", str(S)], cwd=outdir)
        mc_poscar = poscar_re
        mc_order = "auto"
        do_dft = True
        do_render = True
        render_kind = "kagome"
    elif lattice == "honeycomb":
        # 紧凑胞标准化 + 四磁态 + 系数
        run("standardize_cell.py", [poscar_re, "--element", element,
                                    "-o", os.path.join(outdir, "01_structure", "POSCAR.std")],
            cwd=outdir)
        poscar_std = os.path.join(outdir, "01_structure", "POSCAR.std")
        run("build_magnetic_cases3.py",
            [poscar_std, "--supercell"] + [str(x) for x in args.supercell] +
            ["--outdir", os.path.join(outdir, "02_supercell", "cases"),
             "--element", element, "--nonmag"] + nonmag +
            ["--moment", str(moment), "--spin", str(S)], cwd=outdir)
        mc_poscar = poscar_std
        mc_order = "auto"
        do_dft = True
        do_render = True
        render_kind = "honeycomb"
    else:
        raise SystemExit(f"无法处理晶格类型 {lattice!r} ({n_mag} 磁性原子/胞, {n_clusters} 簇); "
                         f"用 --lattice 指定或路由到 comp-chem-workflow 人工处理")

    if cmd == "build":
        print(f"\n[run] build 完成: 02_supercell/cases (lattice={lattice})")
        return

    # ---- DFT 输入 ----
    # P3.16: 自旋模式由 --noncollinear / --soc 决定; 三角晶格恒为非共线（适配层自带该开关）。
    _soc_flags = (["--lsorbit"] if args.soc else []) + \
                 (["--saxis"] + [str(x) for x in args.saxis] if args.saxis else [])
    # P4.2: k 网格 / NPAR / ISYM 只有显式给了才往装配脚本传（缺省保持模板默认，快照不变）
    _kflags = []
    if args.kpoints:
        _kflags += ["--kpoints"] + [str(x) for x in args.kpoints]
    if args.npar is not None:
        _kflags += ["--npar", str(args.npar)]
    if args.isym is not None:
        _kflags += ["--isym", str(args.isym)]
    if run_dft and do_dft and not args.skip_dft:
        if render_kind == "triangular":
            run("assemble_dft_tri.py",
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--U", str(U), "--moment", str(moment)] + _soc_flags + _kflags, cwd=outdir)
        else:
            _nc = ["--noncollinear"] if (args.noncollinear or args.soc or args.cant_theta) else []
            run("assemble_dft.py",
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--U", str(U), "--moment", str(moment)] + _nc + _soc_flags + _kflags, cwd=outdir)

    if run_dft and do_dft and not args.skip_dft:
        _dft_what = ("非共线+SOC" if args.soc else
                     ("非共线" if (args.noncollinear or render_kind == "triangular") else "共线"))
        log(f"[dft] VASP 输入已导出 (03_dft/ + runs/<order>/, 自旋模式: {_dft_what}"
            + (f", k 网格 {' '.join(str(x) for x in args.kpoints)}" if args.kpoints else "")
            + (f", NPAR {args.npar}" if args.npar is not None else "")
            + (f", ISYM {args.isym}" if args.isym is not None else "")
            + "); 请手动在集群提交, 本流水线不代提交")

    if cmd == "dft":
        print(f"\n[run] dft 完成: 03_dft/ + runs/<order>/ 已导出 (lattice={lattice}, 未提交 VASP)")
        return

    # ---- MC ----
    mc_refine = None
    if run_mc and not args.skip_mc:
        # 交变磁体 (P5) 等非默认分组: 允许用 --mc-order 覆盖按 lattice 的默认序参量
        if args.mc_order and args.mc_order != mc_order:
            log(f"[mc] 序参量被 --mc-order 覆盖: {mc_order} -> {args.mc_order}")
            mc_order = args.mc_order
        os.makedirs(os.path.join(outdir, "05_mc"), exist_ok=True)
        model = os.path.join(outdir, "05_mc", "mc_model.json")
        run("poscar_to_mc3.py",
            [mc_poscar,
             "--element", element, "--S", str(S),
             "--J1", str(J1), "--J2", str(J2), "--J3", str(J3), "--MAE", str(MAE),
             "-o", model], cwd=outdir)
        verify_mc_shells(outdir, model)
        run("mc_2d.py", [model, "--order", mc_order,
                         "--size"] + [str(x) for x in args.mc_size] +
            ["--tmin", str(args.mc_tmin), "--tmax", str(args.mc_tmax), "--nt", str(args.mc_nt),
             "--seeds", str(args.mc_seeds),
             "--cool", "--out", os.path.join(outdir, "05_mc", "mc")], cwd=outdir)
        # 峰判据落窗口边界 (网格跨过窄峰) -> 自动缩放 + 加密重扫, 免得把边界值当 T_C/T_N
        mc_refine = refine_mc_peak(outdir, model, mc_order, args.mc_size,
                                   args.mc_refine_nt, args.mc_refine_max,
                                   enabled=not args.no_mc_refine, seeds=args.mc_seeds)
        # 可选: 多尺寸有限尺寸标度 (Binder U4 交点) 自动定 T_C/T_N, 比单一尺寸峰判据更稳
        if args.fss:
            fss_out = os.path.join(outdir, "05_mc", "fss")
            run("mc_fss.py", [model, "--order", mc_order,
                              "--sizes"] + [str(x) for x in args.fss_sizes] +
                ["--tmin", str(args.mc_tmin), "--tmax", str(args.mc_tmax),
                 "--nt", str(args.mc_nt), "--seeds", str(args.mc_seeds),
                 "--cool", "--out", fss_out], cwd=outdir)
        # 附带 MC 脚本到 05_mc/, 使每条流水线可脱离主仓独立重跑 MC
        ship = ["mc_2d.py", "mc_fss.py", "poscar_to_mc3.py", "write_hamiltonian.py"]
        if lattice == "cluster":
            ship += ["build_cluster_cases.py", "solve_cluster_J.py"]
        elif lattice == "triangular":
            ship += ["build_magnetic_cases_tri.py", "solve_J_tri.py"]
        elif lattice == "rectangular":
            ship += ["build_rectangular_cases.py", "solve_J_rect.py"]
        elif lattice == "square":
            ship += ["build_square_cases.py", "solve_J_square.py"]
        elif lattice == "rect_centered":
            ship += ["build_crect_cases.py", "solve_J_crect.py", "solve_J_crect_cant.py"]
        elif lattice == "kagome":
            ship += ["build_kagome_cases.py", "solve_J_kagome.py"]
        for name in ship:
            shutil.copy(os.path.join(SCRIPTS, name), os.path.join(outdir, "05_mc", name))
        log("[run] 附带 MC 脚本 -> " + ", ".join(f"05_mc/{n}" for n in ship))
        # 生成海森堡哈密顿卡片 (顶层, 与 workflow.md 并列)
        run("write_hamiltonian.py", [model, "--label", formula or element,
                                     "-o", os.path.join(outdir, "hamiltonian.md")], cwd=outdir)

    if cmd == "mc":
        print(f"\n[run] mc 完成: 05_mc/mc_mt.csv/.png (lattice={lattice})")
        _print_summary(outdir, lattice, S, J1, J2, J3, MAE, mc_order, mc_refine, args)
        return

    # ---- render 超胞/磁序图 (蜂窝 / 簇磁体 / 三角) ----
    if run_render and do_render:
        # matplotlib 不可用时 (本机缺依赖) 三条渲染路径都自动改走纯 Pillow 的等价脚本:
        # 同样的 CLI / 输入 / 输出, 几何与配色复用各自 matplotlib 版 (均已改为延迟 import)。
        _mpl = have_matplotlib()
        rend_super = "render_supercell.py" if _mpl else "render_supercell_pil.py"
        rend_clu = "render_cluster_orders.py" if _mpl else "render_cluster_orders_pil.py"
        rend_tri = "render_tri_orders.py" if _mpl else "render_tri_orders_pil.py"
        if not _mpl:
            log("[warn] matplotlib 不可用 (缺依赖) -> 磁序图改用 Pillow 兜底: "
                "render_supercell_pil.py / render_cluster_orders_pil.py / render_tri_orders_pil.py")
        if render_kind == "cluster":
            run(rend_clu,
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--cutoff", str(args.cluster_cutoff), "--label", formula or element,
                 "-o", os.path.join(outdir, "figures", "supercell_orders.png")], cwd=outdir)
        elif render_kind == "triangular":
            run(rend_tri,
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--label", formula or element,
                 "-o", os.path.join(outdir, "figures", "supercell_orders.png")], cwd=outdir)
        elif render_kind == "rectangular":
            run(rend_super,
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--label", formula or element, "--orders", "fm", "afm_a", "afm_b", "afm_diag",
                 "-o", os.path.join(outdir, "figures", "supercell_orders.png")],
                cwd=outdir)
        elif render_kind == "square":
            # 复用通用 render_supercell.py: 它从 magnetic_orders.json 读壳距离与 coefficients_per_mag,
            # 会把 4 个方形磁序的 (J1,J2,J3) 系数标在副标题上; 单独再写 render_square.py 是重复代码
            run(rend_super,
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--label", formula or element, "--orders", "fm", "stripe_a", "p4_a", "mix_ab",
                 "-o", os.path.join(outdir, "figures", "supercell_orders.png")],
                cwd=outdir)
        elif render_kind == "rect_centered":
            # 复用通用 render_supercell.py (同方形): 它从 magnetic_orders.json 读壳距离与
            # coefficients_per_mag, 会把 4 个磁序的 (J1,J2,J3) 系数标在副标题上
            run(rend_super,
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--label", formula or element, "--orders", "fm", "afm_c", "afm_a", "afm_b",
                 "-o", os.path.join(outdir, "figures", "supercell_orders.png")],
                cwd=outdir)
        elif render_kind == "kagome":
            # 复用 render_tri_orders.py: 它是仓库里唯一已支持**3 分量非共线 magmom** 的
            # 渲染器 (每原子 reshape(-1,3) 后画面内箭头), 且壳距离从 magnetic_orders.json
            # 的 shells_angstrom 读、序列表直接从 cases/ 目录扫 —— 没有三角专属硬编码。
            # render_supercell.py 只认共线 (mom[k] > 0), 给它 kagome 的 108 个值会画出乱码。
            # 所以「P4.1 新增 render_kagome.py」由复用等效实现 (与方形/C 心矩形复用
            # render_supercell.py 同一口径, 避免第二份 3 分量箭头代码)。
            run(rend_tri,
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--label", formula or element,
                 # 复用三角渲染器时必须覆盖标题: 默认标题写的是「三角晶格」, 对 kagome 是错的
                 "--title", f"{formula or element} kagome (Fe3Sn 型) 磁序俯视图（含 120° 非共线）",
                 "-o", os.path.join(outdir, "figures", "supercell_orders.png")],
                cwd=outdir)
        else:
            run(rend_super,
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--label", formula or element,
                 "-o", os.path.join(outdir, "figures", "supercell_orders.png")],
                cwd=outdir)

    if cmd == "render":
        print(f"\n[run] render 完成: figures/supercell_orders.png (lattice={lattice})")
        return

    # ---- 7. workflow.md ----
    write_workflow(outdir, formula, element, nonmag, S, moment, U, J1, J2, J3, MAE,
                   args, lit, lattice, mc_refine, lat_warn)

    if cmd == "report":
        print(f"\n[run] report 完成: workflow.md (lattice={lattice})")
        return

    print("\n[run] 完成。项目目录:", outdir)
    _print_summary(outdir, lattice, S, J1, J2, J3, MAE, mc_order, mc_refine, args)
    print(f"      产出: workflow.md, hamiltonian.md, 01_structure/, 02_supercell/, 03_dft/, 05_mc/, figures/")


def write_workflow(outdir, formula, element, nonmag, S, moment, U, J1, J2, J3, MAE,
                   args, lit, lattice, mc_refine=None, lat_warn=None):
    if lattice == "triangular":
        gs = "120°/螺旋 (三角受挫)" if J1 < 0 else "铁磁"
    elif lattice == "rect_centered":
        # J1 连的是胞内两个 basis ⇒ J1<0 时基态就是 C 心 Néel (afm_c), 这正是必须支持的序
        gs = "C 心 Néel 反铁磁 (afm_c, 胞内两个 basis 反平行)" if J1 < 0 else "铁磁"
    elif lattice == "kagome":
        # 角共享三角形 ⇒ 几何受挫: J1<0 时经典基态是 q=0 的 120° 态 (与 √3×√3 态近简并)
        gs = "q=0 120° 非共线 (kagome 受挫; √3×√3 态近简并)" if J1 < 0 else "铁磁"
    else:
        gs = "Néel 反铁磁" if J1 < 0 else "铁磁"
    src = lit.get("source") if lit else "手动给定 (未核实文献)"
    # MC 峰分辨记录 (窗口/格距/是否自动缩放重扫) —— 防止把窗口边界值当成 T_C/T_N
    peak_lbl = "Cv 峰" if lattice in ("triangular", "cluster", "kagome") else "χ 峰"
    if mc_refine is None:
        mc_note = "- MC 峰分辨: 未跑 MC (`--skip-mc`)"
    elif not mc_refine["rounds"]:
        s = mc_refine["final"]
        mc_note = (f"- MC 峰分辨: 峰在首轮窗口内部 (窗口 [{s['T_lo']:.1f}, {s['T_hi']:.1f}] K, "
                   f"格距 {s['dT']:.2f} K), {peak_lbl} T ≈ {s['T_peak']:.2f} K")
    elif mc_refine["resolved"]:
        s0, sn = mc_refine["survey"], mc_refine["final"]
        win = " -> ".join(f"[{r['window'][0]:.2f}, {r['window'][1]:.2f}]" for r in mc_refine["rounds"])
        mc_note = (f"- MC 峰分辨: 首轮窗口 [{s0['T_lo']:.1f}, {s0['T_hi']:.1f}] K "
                   f"(格距 {s0['dT']:.2f} K) 的峰值落在 {s0['T_peak']:.2f} K 边界 → 自动缩放重扫 "
                   f"{len(mc_refine['rounds'])} 轮 {win} K (nt={args.mc_refine_nt}), "
                   f"{peak_lbl} T ≈ {sn['T_peak']:.2f} K; 首轮存 `05_mc/mc_mt_survey.csv/.png`")
    else:
        s = mc_refine["final"]
        mc_note = (f"- MC 峰分辨: ⚠ {peak_lbl}仍落在窗口边界 {s['T_peak']:.2f} K → "
                   f"T_C/T_N **未分辨** (已重扫 {len(mc_refine['rounds'])} 轮); "
                   f"请手动收窄 `--mc-tmin/--mc-tmax` 或加密 `--mc-nt`")
    _dft_tag = ("非共线+SOC (LNONCOLLINEAR + LSORBIT)" if args.soc else
                "非共线 (LNONCOLLINEAR)" if args.noncollinear else "共线 (ISPIN=2)")
    _dft_line = f"- DFT 输入 ({_dft_tag}): `03_dft/` + `runs/<order>/` (未运行, 待 POTCAR)\n"
    _dft_tag_tri = ("非共线+SOC (LNONCOLLINEAR + LSORBIT)" if args.soc else "非共线 (LNONCOLLINEAR)")
    _dft_line_tri = f"- DFT 输入 ({_dft_tag_tri}): `03_dft/` + `runs/<order>/` (未运行, 待 POTCAR)\n"
    if lattice == "honeycomb":
        model_desc = ("蜂窝 Heisenberg 模型 H = -J1 Σ_<ij> S·S - J2 Σ_<<ij>> S·S "
                      "- J3 Σ_<<<ij>>> S·S (J>0 铁磁)")
        mc_desc = "序参量 " + ("交错(反铁磁)" if J1 < 0 else "总磁化(铁磁)")
        struct_desc = ("- 超胞/磁态: `02_supercell/cases/{{fm,neel,zigzag,stripy}}` "
                       "+ `magnetic_orders.json`\n- DFT 输入: `03_dft/` + `runs/<order>/` (未运行, 待 POTCAR)")
    elif lattice == "rectangular":
        model_desc = ("矩形晶格单离子 Heisenberg 模型 H = -J1 Σ_<ij> S·S - J2 Σ_<<ij>> S·S "
                      "- J3 Σ_<<<ij>>> S·S (J1=a 轴, J2=b 轴, J3=对角; J>0 铁磁)")
        mc_desc = "序参量 " + ("总磁化(铁磁)" if J1 > 0 else "交错(反铁磁)")
        struct_desc = ("- 磁序枚举: `02_supercell/cases/{{fm,afm_a,afm_b,afm_diag}}` "
                       "+ `magnetic_orders.json`\n"
                       + _dft_line +
                       "- 反解: `05_mc/solve_J_rect.py` (OUTCAR -> J1/J2/J3)\n"
                       "- 图: `figures/supercell_orders.png` (四共线磁态)")
    elif lattice == "square":
        model_desc = ("方形晶格单离子 Heisenberg 模型 H = -J1 Σ_<ij> S·S - J2 Σ_<<ij>> S·S "
                      "- J3 Σ_<<<ij>>> S·S (J1=a 近邻, J2=对角 a√2, J3=2a 同轴向; J>0 铁磁)")
        mc_desc = "序参量 " + ("总磁化(铁磁)" if J1 > 0 else "交错(反铁磁)")
        struct_desc = ("- 磁序枚举: `02_supercell/cases/{{fm,g_type,stripe_a,stripe_b,p4_a,mix_ab}}` "
                       "+ `magnetic_orders.json` (selected = fm / stripe_a / p4_a / mix_ab)\n"
                       "- ⚠ 必须含周期 4 的 `p4_a`: 一切周期 2 的序 (fm/g_type/stripe) 给出的 J3 "
                       "系数完全相同, 少了它 J3 与 E0 在能量里无法分离 (与 kagome 的 q=0 简并同构)\n"
                       + _dft_line +
                       "- 反解: `05_mc/solve_J_square.py` (OUTCAR -> J1/J2/J3)\n"
                       "- 图: `figures/supercell_orders.png` (四共线磁态)")
    elif lattice == "rect_centered":
        model_desc = ("C 心矩形晶格**多亚晶格** Heisenberg 模型 H = -J1 Σ_<ij>C心 S·S "
                      "- J2 Σ_<ij>沿a S·S - J3 Σ_<ij>沿b S·S "
                      "(2 磁原子/胞; J1 = C 心矢量 x4, 连**同一胞内的两个**磁原子; "
                      "J2 = 沿 a x2; J3 = 沿 b x2; J>0 铁磁, meV)")
        mc_desc = "序参量 " + ("总磁化(铁磁)" if J1 > 0 else "交错(反铁磁)")
        struct_desc = ("- 磁序枚举: `02_supercell/cases/"
                       "{{fm,afm_c,afm_a,afm_b,afm_diag,afm_c_a,afm_c_diag}}` "
                       "+ `magnetic_orders.json` (selected = fm / afm_c / afm_a / afm_b)\n"
                       "- ⚠ 必须含 `afm_c` (胞内两个 basis 反平行): 它是唯一能把 J1 从 E0 里"
                       "分出来的序 —— 少了它 J1 与 E0 在设计矩阵里合并\n"
                       "- ⚠ 壳距离必须显式声明: J1=3.0347 A 与 J2=3.1980 A 只差 0.163 A, "
                       "广义 `magnetic_shells(tol=0.2)` 会把两个物理上不同的交换路径**静默合并**\n"
                       + _dft_line +
                       "- 反解: `05_mc/solve_J_crect.py` (OUTCAR -> J1/J2/J3)\n"
                       "- 图: `figures/supercell_orders.png` (四共线磁态)")
    elif lattice == "kagome":
        model_desc = ("kagome 角共享三角形网 **3 亚晶格** Heisenberg 模型 "
                      "H = -J1 Σ_<ij>三角形内 S·S - J2 Σ_<<ij>>跨六角 S·S "
                      "- J3 Σ_<<<ij>>>同亚晶格 S·S (3 磁原子/胞; J1 = a/2 x4 [J1 与 J1' "
                      "同距离不同亚晶格对, 本模型取其平均], J2 = a√3/2 x4, J3 = a x6; "
                      "J>0 铁磁, meV)")
        mc_desc = "序参量 比热 Cv (受挫/非二分晶格)"
        struct_desc = ("- 磁序枚举: `02_supercell/cases/{{fm,q0_120,q0_ud,q0_duu,sq3_120,"
                       "sq3_updntri,sq3_alt}}` + `magnetic_orders.json` "
                       "(selected = fm / q0_120 / sq3_updntri / sq3_alt, rank=4, cond=98.51)\n"
                       "- ⚠ **必须含 √3×√3 序**: q=0 一族 (fm/q0_120/q0_ud/q0_duu) 满足 "
                       "c1 == c2 且 c3 相同 ⇒ 其子块 rank=2 (cond≈2.3e34), 只用 q=0 序反解会得到"
                       "一个「看似正常」的 J1/J2 —— 其实只是 (J1+J2)/2 的任意分配。"
                       "`sq3_alt` 的 c1 != c2 把 J1/J2 两列分开, `sq3_updntri` 的 c3 与 q=0 不同"
                       "把 J3 与 E0 分开。solve 端 rank<4 直接报错 (不静默给 J)\n"
                       "- ⚠ 反解出的 J1 是三角形两条不等价最近邻键 (J1/J1') 的**平均**; "
                       "把二者分开需要降对称超胞 (P4.2)\n"
                       "- 超胞的 n0/n1 必须是 3 的倍数 (√3×√3 的定义相位要求周期 3)\n"
                       + _dft_line +
                       "- 反解: `05_mc/solve_J_kagome.py` (OUTCAR -> J1/J2/J3 + `J_result.json`)\n"
                       "- 图: `figures/supercell_orders.png` (七磁态俯视, 含 120° 矢量箭头)")
    elif lattice == "cluster":
        model_desc = ("簇磁体: 磁性原子成簇 (Nb3 等), 每簇 = 1 个有效自旋, 簇中心成"
                      "三角/非二分晶格; 低能可能是单重态, 经典 Heisenberg 为近似")
        mc_desc = "序参量 比热 Cv (受挫/非二分晶格)"
        struct_desc = ("- 簇识别: `01_structure/POSCAR.clusters.vasp` (磁性簇 → 有效位)\n"
                       "- 磁序枚举: `02_supercell/cases/` (fm + 波矢态) + `magnetic_orders.json`"
                       " (模型 E = E0 + J_intra*C_intra + Σ Jk*Ck)\n"
                       "- DFT 输入: `03_dft/` + `runs/<order>/` (未运行, 待 POTCAR)\n"
                       "- 反解: `05_mc/solve_cluster_J.py` (OUTCAR -> E0/J_intra/J1..Jk)\n"
                       "- 图: `figures/supercell_orders.png` (簇心磁序视图)")
    else:  # triangular
        model_desc = "三角晶格单离子 Heisenberg 模型 (1 磁性原子/胞), 含非共线 120°"
        mc_desc = "序参量 比热 Cv (受挫/非二分晶格)"
        struct_desc = ("- 磁序枚举: `02_supercell/cases/` (fm + 120° 非共线 + 共线 stripy/uud) "
                       "+ `magnetic_orders.json`\n"
                       + _dft_line_tri +
                       "- 反解: `05_mc/solve_J_tri.py` (OUTCAR -> J1/J2/J3)\n"
                       "- 图: `figures/supercell_orders.png` (含 120° 矢量箭头)")

    lat_note = ""
    if lat_warn:
        lat_note = ("8. **本体系被判为「分子簇磁体」, 但可能是多亚晶格铁磁网络**: "
                    + lat_warn + "\n")
    elif lattice == "square":
        # 方形两条专属局限, 不写进通用清单 (会污染其它晶格)
        lat_note = ("8. **方形单层的层间耦合被丢掉**: K2NiF4 / La2CuO4 都是 3D 氧化物 (非 vdW 可解理, "
                    "层间靠 K–F / La–O 键连接), 切成单层后层间 J 被显式丢弃 ⇒ 二维 MC 的临界温度与体相 "
                    "T_N (K2NiF4 97.1 K / La2CuO4 325 K) 只能当**量级检查**, 不能当验证。\n"
                    "9. **J3 依赖周期 4 的 `p4_a` 磁序**: 少了它, 一切周期 2 的序给出相同 J3 系数, "
                    "J3 与 E0 在设计矩阵里合并且 cond(A) 仍可能很小 —— 所以「cond 没超阈」不等于"
                    "「三个壳都分得开」, 必须同时检查 selected 里有没有破简并的序。\n")
    elif lattice == "rect_centered":
        # C 心矩形三条专属局限, 不写进通用清单 (会污染其它晶格)
        lat_note = ("8. **C 心矩形是「多亚晶格」, 不是单离子晶格**: 每胞 2 个磁性原子, 最近邻 J1 "
                    "(3.0347 A) 是 C 心矢量, 连的是**同一胞内的两个**原子。推论: (a) 磁序集合必须含 "
                    "basis 自由度 (`afm_c`), 否则 J1 与 E0 在设计矩阵里合并; (b) `detect_lattice` "
                    "只按几何并簇时会把 2 个原子并成 1 个「簇」并路由到 Nb3I8 三聚体模型 —— "
                    "所以簇分支必须先按正交网判据分开 (P2 已修); (c) 反解出的 J3 依赖 `afm_a`/`afm_b` "
                    "这一对沿 a/b 的序 —— a≠b 使它们不简并, 因此**不需要**方形那样的周期 4 序。\n"
                    "9. **CrOCl 的文献 J 未核实**: `literature_J.json` 的 J1/J2/J3 = null ⇒ 跑 MC 必须"
                    "显式给 `--J1/--J2/--J3`; `solve --demo` 也必须给 `--synthetic J1 J2 J3` 或一个"
                    "已核实的 `--ref` (用 0 顶上会被显式拒绝)。物性口径 (层内铁磁 + 层间反平行的 "
                    "A 型反铁磁) 与体相 T_N 均待正文核实; 2D 模型显式丢掉层间耦合, 所以 MC 临界温度"
                    "**不能**对上体相 T_N。\n"
                    "10. **J1/J2 只差 0.163 A**: 这是 oxyhalide C 心网的结构特征 (C 心键 3.0347 A vs "
                    "a 轴 3.1980 A)。任何用「容差自动归并距离壳」的旧工具都会把两个交换路径算成一个, "
                    "而反解仍然收敛、残差仍然是 0 —— 唯一的防线是产物里的显式壳距离与 tol。\n")
    elif lattice == "kagome":
        # kagome 三条专属局限 (与方形局限 9 同构: 「cond 没超阈」≠「三个壳都分得开」)
        lat_note = ("8. **q=0 一族对 J1/J2 简并**: fm / q0_120 / q0_ud / q0_duu 的系数列满足 "
                    "c1 == c2 且 c3 相同 ⇒ 该子块 rank=2 (cond≈2.3e34, 数值上完全奇异)。"
                    "只用 q=0 序反解**不会报错**, 会给出一个看似正常的 J1/J2 —— 它们只是 "
                    "(J1+J2)/2 的任意分配。必须含 √3×√3 序: `sq3_alt` 的 c1 != c2 分开 J1/J2 两列, "
                    "`sq3_updntri` 的 c3 与 q=0 一族不同把 J3 与 E0 分开。"
                    "`solve_J_kagome.py` 对 rank<4 与 cond>1e4 都直接报错, 不静默给 J。\n"
                    "9. **J1 与 J1' 同距离, 本模型取其平均**: 角共享三角形的两条最近邻键都是 a/2, "
                    "但连接的是不同的亚晶格对。3 壳方案把二者并成一条 J1 ⇒ 反解值与 MC 的 T_N 只对 "
                    "(J1+J1')/2 的组合敏感; 分开二者需要降对称 (√3×√3 或 2×2) 的超胞与更多磁序 (P4.2 待做)。"
                    "同理 `poscar_to_mc3 --jfile` 必须用 `pair=` 显式声明壳, 旧的自动路径会把同距离的两条并成一条。\n"
                    "10. **Fe3Sn 的文献 J 未核实 + S 是占位值**: `literature_J.json` 里有 Fe3Sn 条目(铁磁基态 + "
                    "T_C = 705 K + 结构/各向异性锚点), 但 J1/J2/J3 按设计仍是 null (原文自承交换常数未确立) ⇒ "
                    "跑 MC 必须显式给 `--J1/--J2/--J3`; `solve --demo` 必须给 `--synthetic` 或已核实的 `--ref`。"
                    "预设里的 S=1.0 / moment=2.0 / U=4.0 是占位值 (元素表不列 Fe, 需自洽计算), "
                    "所以本晶格目前是**几何/线性代数层可用, 物理结论待 DFT**。\n")

    lines = f"""# {formula} 单层 磁耦合 (J) → Monte Carlo 通用流水线

## 模型与方法

| 项 | 值 |
|---|---|
| 体系 | {formula} 单层, 磁性格点类型 = **{lattice}** |
| 磁性元素 / 非磁 | {element} / {', '.join(nonmag) if nonmag else '-'} |
| 自旋 S / 磁矩 | S={S} / {moment} μB |
| 泛函 | PBE + DFT+U (U={U} eV 加 {element} d) |
| 基态 | {gs} |
| 交换常数 (meV) | J1={J1}, J2={J2}, J3={J3}, MAE={MAE} |

模型: {model_desc}

J 值来源: {src}

## 结果

- 结构: `01_structure/POSCAR` (reorient 后, 真空 {args.vacuum} Å)
{struct_desc}
- MC: `05_mc/mc_mt.csv/.png`, {mc_desc}
{mc_note}
- 哈密顿卡片: `hamiltonian.md` (顶层, 符号式 H + 壳摘要 + 逐键明细)
- MC 脚本 (自包含, 可独立重跑): `05_mc/mc_2d.py` (引擎) + `05_mc/poscar_to_mc3.py` (POSCAR -> 模型) + `05_mc/write_hamiltonian.py`

## 局限 (通用)

1. 磁性格点类型 {lattice}; 蜂窝/簇/三角/矩形/方形/C 心矩形均已支持磁序枚举 + DFT 输入 + 反解脚本 (三角用非共线 120°); J 用文献/给定, 未跑 DFT 前非本项目计算值。
2. J 取自文献/给定; 未跑 DFT 前非本项目计算值。
3. 纯 Heisenberg MAE 默认为 0 → 二维无真正相变 (Mermin-Wagner), 临界温度为有限尺寸估计。
4. 簇磁体可能是单重态 (非磁有序), 经典 Heisenberg-MC 的 T_N **无物理意义**, 仅作机制演示 (exploratory)。
5. mp_fetch_2d.py 的 exfoliate 对菱面体层状/已存成紧胞的材料有已知缺口。
6. k 网格 3x3x1、U 值为可调假设。
7. 峰判据 (Cv/χ) 只有落在温度窗口**内部**才算分辨: 峰落在首/末格点 = 窗口截断或网格跨过窄峰, 该边界值**不是** T_C/T_N。driver 检测到这种情形会自动缩放窗口 + 加密重扫 (可调 `--mc-refine-nt`/`--mc-refine-max`, `--no-mc-refine` 关闭); 仍落在边界则结论显式标「未分辨」。
{lat_note}"""
    with open(os.path.join(outdir, "workflow.md"), "w", encoding="utf-8") as f:
        f.write(lines)
    log("[run] 写出 workflow.md")


if __name__ == "__main__":
    main()
