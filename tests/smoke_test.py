#!/usr/bin/env python
"""smoke_test.py — 回归测试：跑关键确定性步骤，断言输出与已验证值一致。

运行:
  python tests/smoke_test.py

覆盖:
  1. detect_lattice 正确分类 (蜂窝 / 簇磁体)
  2. cluster_sites: Nb3I8 的 3 个 Nb -> 1 个簇
  3. solve_J3 --demo 精确还原文献 J (PASS)
  4. build_magnetic_cases3 蜂窝系数 = 教科书值

依赖已验证结构:
  C:/Users/吴/Desktop/项目/MnPS3_pipeline/01_structure/POSCAR
  C:/Users/吴/Desktop/项目/Nb3I8_pipeline/01_structure/POSCAR
"""
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(HERE, "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import cluster_sites as CS  # noqa: E402
from run_2d_magnet import detect_lattice  # noqa: E402

MNPS3 = r"C:\Users\吴\Desktop\项目\MnPS3_pipeline\01_structure\POSCAR"
NB3I8 = r"C:\Users\吴\Desktop\项目\Nb3I8_pipeline\01_structure\POSCAR"

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


# ---- 1. detect_lattice ----
check("detect MnPS3 = honeycomb", detect_lattice(MNPS3, "Mn")[0] == "honeycomb")
check("detect Nb3I8 = cluster", detect_lattice(NB3I8, "Nb")[0] == "cluster")

# ---- 2. cluster_sites: 3 Nb -> 1 簇 ----
lat, atoms, coords = CS.parse_poscar(NB3I8)
groups = CS.find_clusters(atoms, coords, lat, "Nb", 3.3)
check("Nb3I8: 3 Nb -> 1 cluster", len(groups) == 1 and len(groups[0]) == 3)

# ---- 3. solve_J3 --demo ----
r = subprocess.run(
    [sys.executable, os.path.join(SCRIPTS, "solve_J3.py"),
     os.path.join(r"C:\Users\吴\Desktop\项目\MnPS3_pipeline\02_supercell\cases"), "--demo"],
    capture_output=True, text=True)
check("solve_J3 --demo PASS", r.returncode == 0 and "PASS" in r.stdout)

# ---- 4. build_magnetic_cases3 蜂窝系数 ----
with tempfile.TemporaryDirectory() as td:
    subprocess.run(
        [sys.executable, os.path.join(SCRIPTS, "build_magnetic_cases3.py"),
         MNPS3, "--supercell", "2", "2", "1",
         "--outdir", os.path.join(td, "cases"),
         "--element", "Mn", "--nonmag", "P", "S",
         "--moment", "5.0", "--spin", "2.5"],
        capture_output=True)
    meta = json.load(open(os.path.join(td, "cases", "magnetic_orders.json"), encoding="utf-8"))
    c = meta["coefficients_per_cell"]
    check("fm coeff (-12,-24,-12)", np.allclose(c["fm"], [-12, -24, -12]))
    check("neel coeff (+12,-24,+12)", np.allclose(c["neel"], [12, -24, 12]))
    check("zigzag coeff (+4,+8,-12)", np.allclose(c["zigzag"], [4, 8, -12]))
    check("stripy coeff (-4,+8,+12)", np.allclose(c["stripy"], [-4, 8, 12]))

print("\n全部通过" if not fails else f"\n{len(fails)} 个失败: {fails}")
sys.exit(1 if fails else 0)
