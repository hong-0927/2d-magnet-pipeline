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
}


def log(*a, **k):
    if not QUIET:
        print(*a, **k)


def run(pyfile, args, cwd=None):
    """跑 scripts/ 下的脚本。QUIET 时捕获 stdout/stderr, 只在失败或警告时回显。"""
    cmd = [sys.executable, os.path.join(SCRIPTS, pyfile)] + args
    log("\n$", " ".join(cmd))
    if QUIET:
        r = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"[FAIL] {pyfile} 退出码 {r.returncode}\n"
                  f"$ {' '.join(cmd)}\n{r.stdout}{r.stderr}", file=sys.stderr)
            sys.exit(r.returncode)
        if r.stderr and r.stderr.strip():
            print(r.stderr.rstrip(), file=sys.stderr)
    else:
        r = subprocess.run(cmd, cwd=cwd)
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
    kw = {"capture_output": True, "text": True} if QUIET else {}
    try:
        r = subprocess.run([sys.executable] + fa, cwd=fetch_dir, **kw)
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
    {"key": "scripts_missing", "static": False,
     "symptom": "driver 报缺脚本 / 子阶段找不到脚本",
     "fix": "确认 scripts/ 目录完整 (重新 git clone 或解压分发包)"},
    {"key": "lit_missing", "static": False,
     "symptom": "--literature <key> 报找不到文献条目 / literature_J.json 缺失",
     "fix": "确认 scripts/literature_J.json 存在且键名正确"},
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
]


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
    line(sys.version_info >= (3, 9),
         f"Python >= 3.9 (当前 {sys.version.split()[0]})", key="py39")
    for pkg in ("numpy", "numba", "matplotlib", "requests"):
        try:
            m = importlib.import_module(pkg)
            line(True, f"{pkg} {getattr(m, '__version__', '?')}")
        except Exception:
            line(False, f"缺少依赖 {pkg} -> pip install -r requirements.txt",
                 key="missing_dep")

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


def detect_lattice(poscar_path, element, cutoff=3.3):
    """读 reorient 后 POSCAR, 数磁性原子/胞 + 簇数, 返回 (lattice, n_clusters, n_mag)。

    lattice ∈ {"honeycomb", "triangular", "rectangular", "cluster", "other"}:
      cluster     = 磁性原子成簇 (簇数 < 原子数), 簇 = 有效磁位
      rectangular = 正交且 a≠b 的矩形晶格 (CrSBr 类, 非受挫)
      triangular  = 1 磁性原子/胞
      honeycomb   = 2 磁性原子/胞 (不分簇)
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

    if n_clusters < n_mag:
        return "cluster", n_clusters, n_mag
    if _is_rectangular_net(lat):
        return "rectangular", n_clusters, n_mag
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


def main():
    ap = argparse.ArgumentParser(description="二维磁体 结构->J->MC 通用流水线")
    # 子命令 (缺省=全流水线; 只跑单个阶段, 复用已有产物)
    ap.add_argument("command", nargs="?", default=None,
                    choices=["fetch", "build", "dft", "mc", "solve", "render", "report", "doctor"],
                    help="只跑单个阶段 (fetch/build/dft/mc/solve/render/report); doctor=环境自检+error catalog; 缺省跑全流水线")
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
    ap.add_argument("--lattice", choices=["auto", "honeycomb", "triangular", "cluster", "rectangular"],
                    default="auto", help="磁性格点类型 (默认自动判别; rectangular=正交 a≠b 矩形)")
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
    ap.add_argument("--fss", action="store_true",
                    help="额外跑多尺寸有限尺寸标度 (Binder U4 交点自动定 T_C/T_N)")
    ap.add_argument("--fss-sizes", type=int, nargs="+", default=[12, 18, 24, 30],
                    help="FSS 的线性尺寸列表 (默认 12 18 24 30)")
    ap.add_argument("--skip-dft", action="store_true")
    ap.add_argument("--skip-mc", action="store_true")
    ap.add_argument("--quiet", action="store_true",
                    help="静默模式: 抑制子进程输出, 只回显错误/警告 + 结尾摘要")
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
        再喂给 mc_2d.py)。J1 的默认是 None, 所以 null 的 J1 仍会触发「缺少必要参数」。
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

    missing = [k for k, v in [("--element", element), ("--S", S), ("--moment", moment),
                              ("--U", U), ("--J1", J1)] if v is None]
    if missing:
        raise SystemExit(f"缺少必要参数 {missing}; 用 --literature <key> 或显式给出")

    if not formula and not args.poscar and not args.mp_id:
        raise SystemExit("需要 --formula / --mp-id / --poscar 之一")

    if args.outdir:
        outdir = os.path.abspath(args.outdir)
    else:
        outdir = os.path.abspath((formula or "material") + "_pipeline")
    os.makedirs(outdir, exist_ok=True)
    log(f"[run] 项目目录: {outdir}")

    cmd = args.command
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
                                   capture_output=True, text=True)
                if r.returncode == 0 and os.path.exists(poscar_mp):
                    log(f"[fetch] {p} 成功")
                    got = True
                    break
                print(f"[fetch] {p} 失败: {(r.stderr or r.stdout or '').strip()[-200:]}", file=sys.stderr)
            if not got:
                raise SystemExit("所有结构 provider 均失败; 用 --poscar 自带单层 POSCAR, "
                                 "或 --provider mp 配 --mp-key")

    # ---- 2. reorient ----
    run("reorient.py", [os.path.join(outdir, "01_structure", "POSCAR.mp"),
                        "--vacuum", str(args.vacuum),
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
    lat_warn = multi_sublattice_warning(lattice, n_clusters, n_mag, element)
    if lat_warn:
        # 走 stderr: --quiet 下 log() 被抑制, 但这条必须让用户看见
        print(f"[warn] {lat_warn}", file=sys.stderr)

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
        else:
            run("solve_J3.py", [cases], cwd=outdir)
        print(f"\n[run] solve 完成 (lattice={lattice}, cases={cases})")
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
    if run_dft and do_dft and not args.skip_dft:
        if render_kind == "triangular":
            run("assemble_dft_tri.py",
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--U", str(U), "--moment", str(moment)], cwd=outdir)
        else:
            run("assemble_dft.py",
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--species"] + ([element] + nonmag) +
                ["--U", str(U), "--moment", str(moment)], cwd=outdir)

    if run_dft and do_dft and not args.skip_dft:
        log("[dft] VASP 输入已导出 (03_dft/ + runs/<order>/); 请手动在集群提交, 本流水线不代提交")

    if cmd == "dft":
        print(f"\n[run] dft 完成: 03_dft/ + runs/<order>/ 已导出 (lattice={lattice}, 未提交 VASP)")
        return

    # ---- MC ----
    mc_refine = None
    if run_mc and not args.skip_mc:
        os.makedirs(os.path.join(outdir, "05_mc"), exist_ok=True)
        model = os.path.join(outdir, "05_mc", "mc_model.json")
        run("poscar_to_mc3.py",
            [mc_poscar,
             "--element", element, "--S", str(S),
             "--J1", str(J1), "--J2", str(J2), "--J3", str(J3), "--MAE", str(MAE),
             "-o", model], cwd=outdir)
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
        if render_kind == "cluster":
            run("render_cluster_orders.py",
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--cutoff", str(args.cluster_cutoff), "--label", formula or element,
                 "-o", os.path.join(outdir, "figures", "supercell_orders.png")], cwd=outdir)
        elif render_kind == "triangular":
            run("render_tri_orders.py",
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--label", formula or element,
                 "-o", os.path.join(outdir, "figures", "supercell_orders.png")], cwd=outdir)
        elif render_kind == "rectangular":
            run("render_supercell.py",
                ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
                ["--label", formula or element, "--orders", "fm", "afm_a", "afm_b", "afm_diag",
                 "-o", os.path.join(outdir, "figures", "supercell_orders.png")],
                cwd=outdir)
        else:
            run("render_supercell.py",
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
    else:
        gs = "Néel 反铁磁" if J1 < 0 else "铁磁"
    src = lit.get("source") if lit else "手动给定 (未核实文献)"
    # MC 峰分辨记录 (窗口/格距/是否自动缩放重扫) —— 防止把窗口边界值当成 T_C/T_N
    peak_lbl = "Cv 峰" if lattice in ("triangular", "cluster") else "χ 峰"
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
                       "- DFT 输入 (共线): `03_dft/` + `runs/<order>/` (未运行, 待 POTCAR)\n"
                       "- 反解: `05_mc/solve_J_rect.py` (OUTCAR -> J1/J2/J3)\n"
                       "- 图: `figures/supercell_orders.png` (四共线磁态)")
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
                       "- DFT 输入 (非共线): `03_dft/` + `runs/<order>/` (未运行, 待 POTCAR)\n"
                       "- 反解: `05_mc/solve_J_tri.py` (OUTCAR -> J1/J2/J3)\n"
                       "- 图: `figures/supercell_orders.png` (含 120° 矢量箭头)")

    lat_note = ""
    if lat_warn:
        lat_note = ("8. **本体系被判为「分子簇磁体」, 但可能是多亚晶格铁磁网络**: "
                    + lat_warn + "\n")

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

1. 磁性格点类型 {lattice}; 蜂窝/簇/三角均已支持磁序枚举 + DFT 输入 + 反解脚本 (三角用非共线 120°); J 用文献/给定, 未跑 DFT 前非本项目计算值。
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
