#!/usr/bin/env python
"""test_fgt_e2e.py — FGT 双亚晶格反解链路的端到端测试 (mock DFT, 不跑 VASP)

覆盖:
  1. solve_fgt_J.py: 用已知 J 正演出各磁态能量 (4 个写 OUTCAR 的 TOTEN 行, 1 个只写 E.txt,
     两条读能量路径都走到), 反解必须精确还原 J1..J4;
  2. fgt_verify_states.py: 磁矩符号与目标 sigma 一致 -> PASS (退出码 0);
     故意把一个 Fe 的磁矩翻号 -> FAIL (退出码 1), 证明"静默跑到别的态"能被抓住;
  3. fgt_relax_backfill.py: 用假的 CONTCAR 回填, 5 个磁态 POSCAR 必须都换成弛豫后坐标,
     且原坐标存进 POSCAR.pristine; 物种/原子数不符时必须拒绝。

运行:
  python tests/test_fgt_e2e.py

依赖已验证产物:
  C:/Users/吴/Desktop/项目/Fe3GeTe2_pipeline/02_supercell/cases  (magnetic_orders.json + 5 态)
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import uuid

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(HERE, "scripts")
FGT_CASES = r"C:\Users\吴\Desktop\项目\Fe3GeTe2_pipeline\02_supercell\cases"

fails = []
# 子进程统一用 UTF-8 输出、父进程统一按 UTF-8 解码: 免得在中文 Windows 上
# 父环境是 GBK 而子进程被设成 UTF-8 时抛 UnicodeDecodeError。
CHILD_ENV = dict(os.environ, PYTHONIOENCODING="utf-8", PYTHONDONTWRITEBYTECODE="1")


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


def mktmp(prefix="_tmp_"):
    """找一个**真能写**的临时目录。

    踩过的坑:
      1. tempfile.mkdtemp/TemporaryDirectory 会 chmod 0700, 某些受限环境下连创建者自己
         都被拒绝写子目录 -> 自己拼名 + os.makedirs, 用 shutil.rmtree(ignore_errors) 清理;
      2. 受限沙箱会给流水线子目录挂上拒绝"新建条目"的 ACL (仓库根可写、`2d_magnet_pipeline\\`
         内不可写), 所以不能只在仓库内建 -> 按候选链回退。子进程也要能写这里。
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


def fake_outcar(energy_ev, moments):
    """最小可读 OUTCAR: 一行 TOTEN + 一个 magnetization (x) 块。"""
    L = [f"  free  energy   TOTEN  =      {energy_ev:.8f} eV", "",
         " magnetization (x)", "",
         "# of ion       s       p       d       tot",
         "-" * 42]
    for i, m in enumerate(moments, start=1):
        L.append(f"{i:5d}   0.002   0.001   {m - 0.10:7.3f} {m:8.3f}")
    L += ["-" * 42, "tot         0.054   0.027   " + f"{sum(moments):8.3f}", ""]
    return "\n".join(L) + "\n"


def target_moments(magmom_path):
    """magmom.txt -> 每个原子的目标磁矩 (含 0.0 的非磁性原子), 按 POSCAR 顺序。"""
    out = []
    for l in open(magmom_path, encoding="utf-8"):
        l = l.split("#")[0].strip()
        if not l:
            continue
        p = l.split()
        if len(p) >= 4 and p[1].isdigit():
            out.append(float(p[-1]))
    return out


if not os.path.isdir(FGT_CASES):
    print(f"缺依赖目录: {FGT_CASES}")
    sys.exit(1)

meta_src = json.load(open(os.path.join(FGT_CASES, "magnetic_orders.json"), encoding="utf-8"))
rows = meta_src["coefficient_matrix"]["rows"]
labels = meta_src["coefficient_matrix"]["cols"][1:]
M = np.array(meta_src["coefficient_matrix"]["M"], float)

# 已知"真值" J (meV, S=1): 正演出能量, 反解应精确还原
TRUTH = {"J1": 113.0, "J2": -5.5, "J3": -10.9, "J4": 2.25}
E0_TRUE_EV = -52.0

td = mktmp("_tmp_fgte2e_")
try:
    cases = os.path.join(td, "cases")
    runs = os.path.join(td, "runs")
    shutil.copytree(FGT_CASES, cases)
    for o in rows:                      # 外部目录若已回填过, 先从 pristine 恢复干净坐标
        p = os.path.join(cases, o, "POSCAR.pristine")
        if os.path.exists(p):
            shutil.copyfile(p, os.path.join(cases, o, "POSCAR"))
            os.remove(p)

    truth_vec = np.array([TRUTH[l] for l in labels], float)
    E_vec = M[:, 0] * E0_TRUE_EV + M[:, 1:] @ (truth_vec / 1000.0)

    # ---- 1. 合成各磁态能量: 前 4 个走 OUTCAR, 最后一个只给 E.txt ----
    moments_by_order = {}
    for k, r in enumerate(rows):
        d = os.path.join(runs, r)
        os.makedirs(d, exist_ok=True)
        mm = target_moments(os.path.join(cases, r, "magmom.txt"))
        moments_by_order[r] = mm
        if k < len(rows) - 1:
            with open(os.path.join(d, "OUTCAR"), "w", encoding="utf-8") as f:
                f.write(fake_outcar(E_vec[k], mm))
        else:
            with open(os.path.join(d, "E.txt"), "w", encoding="utf-8") as f:
                f.write(f"{E_vec[k]:.8f}\n")
    check("mock: 5 个磁态能量就位 (4 OUTCAR + 1 E.txt)",
          all(os.path.exists(os.path.join(runs, r, "OUTCAR")) or
              os.path.exists(os.path.join(runs, r, "E.txt")) for r in rows))

    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "solve_fgt_J.py"), cases,
                        "--runs", runs, "--spin", "1.0"], capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=CHILD_ENV)
    check("solve_fgt_J 退出码 0", r.returncode == 0)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:])
    else:
        js = json.load(open(os.path.join(td, "J_solved.json"), encoding="utf-8"))
        for l in labels:
            check(f"反解还原 {l} = {TRUTH[l]}", abs(js["J_meV"][l] - TRUTH[l]) < 1e-6)
        check("反解 E0 还原", abs(js["E0_meV_per_cell"] - E0_TRUE_EV * 1000.0) < 1e-6)
        check("反解残差 ~ 0", js["residual_eV"] < 1e-9)
        check("反解两条读能量路径都被用到",
              any("OUTCAR" in s for s in js["sources"]) and any("E.txt" in s for s in js["sources"]))

    # ---- 2. 磁矩核对: 一致 -> PASS ----
    #      verify_states 只看 OUTCAR, 所以给"只有 E.txt"那个态也补上 OUTCAR
    with open(os.path.join(runs, rows[-1], "OUTCAR"), "w", encoding="utf-8") as f:
        f.write(fake_outcar(E_vec[-1], moments_by_order[rows[-1]]))
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fgt_verify_states.py"), cases,
                        "--runs", runs], capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=CHILD_ENV)
    check("verify_states: 符号一致时退出码 0", r.returncode == 0)
    check("verify_states: 输出 VERDICT: PASS", "VERDICT: PASS" in r.stdout)
    check("verify_states: 5 态全部 [OK]",
          r.stdout.count("[OK]") == len(rows) and "[MISMATCH]" not in r.stdout)

    # ---- 2b. 把一个 Fe 的磁矩翻号 -> 必须抓住 ----
    bad_order = "afm_b_vert"
    mm = list(moments_by_order[bad_order])
    flip = next(i for i, m in enumerate(mm) if abs(m) > 1e-9)
    mm[flip] = -mm[flip]
    with open(os.path.join(runs, bad_order, "OUTCAR"), "w", encoding="utf-8") as f:
        f.write(fake_outcar(E_vec[rows.index(bad_order)], mm))
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fgt_verify_states.py"), cases,
                        "--runs", runs], capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=CHILD_ENV)
    check("verify_states: 翻号后退出码 1", r.returncode == 1)
    check("verify_states: 输出 VERDICT: FAIL", "VERDICT: FAIL" in r.stdout)
    check("verify_states: 点名不一致的原子", "[MISMATCH]" in r.stdout and "Fe#" in r.stdout)
    with open(os.path.join(runs, bad_order, "OUTCAR"), "w", encoding="utf-8") as f:
        f.write(fake_outcar(E_vec[rows.index(bad_order)], moments_by_order[bad_order]))

    # ---- 2c. 缺一个磁态 -> 默认失败, --allow-missing 放行 ----
    miss = os.path.join(runs, rows[-1], "OUTCAR")
    os.rename(miss, miss + ".bak")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fgt_verify_states.py"), cases,
                        "--runs", runs], capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=CHILD_ENV)
    check("verify_states: 缺 OUTCAR 默认判失败",
          r.returncode == 1 and "[MISSING]" in r.stdout)
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fgt_verify_states.py"), cases,
                        "--runs", runs, "--allow-missing"], capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=CHILD_ENV)
    check("verify_states: --allow-missing 放行", r.returncode == 0 and "VERDICT: PASS" in r.stdout)
    os.rename(miss + ".bak", miss)

    # ---- 3. CONTCAR 回填 ----
    src = os.path.join(cases, "fm", "POSCAR")
    orig = open(src, encoding="utf-8").read().splitlines()
    cont = list(orig)
    p = cont[8].split()
    cont[8] = f"{float(p[0]) + 0.037:.12f}  {float(p[1]):.12f}  {float(p[2]):.12f}"
    cont[0] = "Fe3GeTe2 relaxed FM (CONTCAR)"
    os.makedirs(os.path.join(runs, "fm"), exist_ok=True)
    with open(os.path.join(runs, "fm", "CONTCAR"), "w", encoding="utf-8") as f:
        f.write("\n".join(cont) + "\n")

    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fgt_relax_backfill.py"), cases],
                       capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=CHILD_ENV)
    check("backfill 退出码 0", r.returncode == 0)
    if r.returncode != 0:
        print(r.stdout[-2000:], r.stderr[-2000:])
    else:
        ok_all, ok_pri, ok_hdr = True, True, True
        for o in rows:
            cp = os.path.join(cases, o, "POSCAR")
            rp = os.path.join(runs, o, "POSCAR")
            body = open(cp, encoding="utf-8").read().splitlines()
            ok_hdr &= body[0].startswith("Fe3GeTe2") and "(relaxed)" in body[0]
            ok_all &= body[8] == cont[8] and os.path.exists(rp) and \
                open(rp, encoding="utf-8").read().splitlines()[8] == cont[8]
            pr = cp + ".pristine"
            ok_pri &= os.path.exists(pr) and \
                open(pr, encoding="utf-8").read().splitlines()[8] == orig[8]
        check("backfill: 5 个磁态 POSCAR 都换成弛豫坐标", ok_all)
        check("backfill: 原坐标存入 POSCAR.pristine", ok_pri)
        check("backfill: 首行保留原注释 + (relaxed)", ok_hdr)

        # 幂等: 再跑一次结果不变, 且 pristine 不被覆盖
        r2 = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fgt_relax_backfill.py"), cases],
                            capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=CHILD_ENV)
        body = open(os.path.join(cases, "fm", "POSCAR"), encoding="utf-8").read().splitlines()
        pr = open(os.path.join(cases, "fm", "POSCAR.pristine"), encoding="utf-8").read().splitlines()
        check("backfill: 可重复运行 (幂等, pristine 不被覆盖)",
              r2.returncode == 0 and body[8] == cont[8] and pr[8] == orig[8] and
              body[0].count("(relaxed)") == 1)
        # 位移检查: 真弛豫过 -> 打印位移, 且**不能**报"逐位相同"
        check("backfill: 打印最大坐标位移且不误报未弛豫",
              "最大坐标位移" in r.stdout and "逐位相同" not in r.stdout)

    # ---- 3b. 物种/原子数不符的 CONTCAR -> 必须拒绝 ----
    bad = list(cont)
    bad[6] = "26 9 18"
    with open(os.path.join(runs, "fm", "CONTCAR"), "w", encoding="utf-8") as f:
        f.write("\n".join(bad) + "\n")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fgt_relax_backfill.py"), cases],
                       capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=CHILD_ENV)
    check("backfill: 原子数不符时拒绝 (非零退出)",
          r.returncode != 0 and "26, 9, 18" in (r.stdout + r.stderr))

    # ---- 3c. CONTCAR 与原始 POSCAR 逐位相同 (= 静态自洽 IBRION=-1 的产物) ----
    #      VASP 照样写 CONTCAR, 于是"回填"等于没做: 必须吼出来而不是静默成功。
    static_cont = list(orig)
    with open(os.path.join(runs, "fm", "CONTCAR"), "w", encoding="utf-8") as f:
        f.write("\n".join(static_cont) + "\n")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fgt_relax_backfill.py"), cases],
                       capture_output=True, text=True,
                       encoding="utf-8", errors="replace", env=CHILD_ENV)
    check("backfill: 未弛豫 (CONTCAR 逐位相同) 时警告静态 INCAR",
          r.returncode == 0 and "[warn]" in r.stdout and "IBRION" in r.stdout and
          "0.0000" in r.stdout)
finally:
    shutil.rmtree(td, ignore_errors=True)

print("\n全部通过" if not fails else f"\n{len(fails)} 个失败: {fails}")
sys.exit(1 if fails else 0)
