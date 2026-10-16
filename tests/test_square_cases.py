#!/usr/bin/env python
"""test_square_cases.py — 方形晶格磁构型回归 (≥20 例, 自带 __main__ 直跑入口)。

设计要点
--------
方形格子 (K2NiF4 型, 1 磁原子/胞) 的系数矩阵有一个**结构性陷阱**:

  J3 = 2a 的分数位移在 n0 方向恰好是**半周期** (4 宽超胞时 (0.5, 0, 0))。
  最小镜像会把 +2a 与 -2a 两个像折叠成同一位点, 每胞只数到 1 根 2a 键
  (无限晶格正确值 2 根) ⇒ c3 少算一半, J3 被系统性低估 2 倍。

所以本测试用三种**互相独立**的算法交叉锚定同一批系数:

  (1) build_square_cases.py 自带的几何实现 (显式平移搜索, classify_bonds nshift=1)
  (2) P0.1 共享规格层 magnetic_orders.py (显式 policy="shift_search", nshift=1)
  (3) 4x4x1 的手推解析理想值 (per mag: 每胞 2 根 NN / 2 根对角 / 2 根 2a 键)

三者必须逐位一致; 任一处漂移立刻炸出来, 而不是等某个产物悄悄错掉。

运行:
  python tests/test_square_cases.py        (不需要 pytest)

环境变量:
  MAGNET_FIXTURES   fixture 目录 (默认 <repo>/tests/fixtures)
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(HERE, "scripts")
FIXTURES = os.environ.get("MAGNET_FIXTURES", os.path.join(HERE, "tests", "fixtures"))
POSCAR_LA = os.path.join(FIXTURES, "La2CuO4_POSCAR")
POSCAR_K = os.path.join(FIXTURES, "K2NiF4_POSCAR")

for _p in (POSCAR_LA, POSCAR_K):
    if not os.path.exists(_p):
        sys.exit("[skip] 缺少 fixture: %s (仓库应自带 tests/fixtures/, 或设 MAGNET_FIXTURES)"
                 % _p)

sys.path.insert(0, SCRIPTS)

import build_square_cases as B  # noqa: E402
import magnetic_orders as MO  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + (("   " + extra) if extra else ""))
    if not cond:
        fails.append(name)


# ---------------------------------------------------------------- 两种独立实现

def builder_coeffs(poscar, element, nonmag, sc):
    """(1) build_square_cases 自带的几何实现 —— 独立基准。"""
    lat, atoms, coords, _mode = B.parse_poscar(poscar)
    lat_sc, atoms_sc, coords_sc, cell_idx = B.make_supercell(lat, atoms, coords, list(sc))
    mag = [i for i, a in enumerate(atoms_sc) if a == element]
    shells = B.magnetic_shells(mag, coords_sc, lat_sc, 3)
    n0 = cell_idx[mag][:, 0]
    n1 = cell_idx[mag][:, 1]
    sig = B.sigmas_of(n0, n1)
    sums = {k: B.classify_bonds(mag, coords_sc, lat_sc, s, shells) for k, s in sig.items()}
    coeff = {k: [-0.5 * s[0], -0.5 * s[1], -0.5 * s[2]] for k, s in sums.items()}
    return coeff, shells, len(mag)


def spec_coeffs(poscar, element, nonmag, sc, shells, spin):
    """(2) P0.1 共享规格层 —— 与 (1) 完全独立的第二套实现。"""
    spec = MO.Spec(
        key="sq_test", label="方形 (test)",
        structure=MO.StructureSpec(poscar=os.path.abspath(poscar), element=element,
                                   nonmag=tuple(nonmag), supercell=tuple(sc)),
        shells=tuple(MO.Shell("J%d" % (i + 1), d) for i, d in enumerate(shells)),
        orders=tuple(MO.Order(nm, (lambda c, nm=nm: B.sigmas_of(c.n0, c.n1)[nm]))
                     for nm in B.ORDERS),
        selected=tuple(B.SELECTED), spin=spin,
        # 必须显式给 shift_search: 规格层默认 min_image, 而 2a 恰是半周期,
        # min_image 会把 ±2a 两个像折叠成同一位点 (c3 少算一半)。
        policy="shift_search", nshift=1,
    )
    res, _ctx = MO.compute(spec)
    return {nm: [float(x) for x in res[nm].coeffs] for nm in B.ORDERS}


def note_no_silent_drift(a, b):
    """两套实现的差: 必须为 0; 但把「错一半」这个经典症状单独标出来。"""
    d = np.asarray(a) - np.asarray(b)
    m = float(np.max(np.abs(d)))
    tag = ""
    if m > 1e-9 and np.allclose(d, -0.5 * np.asarray(b), atol=1e-6):
        tag = "  <== 典型症状: min_image 把 ±2a 折叠, 2a 键少算一半"
    return m, tag


# ------------------------------------------------- (A) 超胞 x 磁序 逐位交叉校验
# 6 个磁序 x 2 个超胞 = 12 例。8x8x1 用来证明这不是 4 宽超胞的巧合。
CACHE = {}
for _sc in ((4, 4, 1), (8, 8, 1)):
    _cb, _sh, _n = builder_coeffs(POSCAR_LA, "Cu", ("La", "O"), _sc)
    _cs = spec_coeffs(POSCAR_LA, "Cu", ("La", "O"), _sc, _sh, 0.5)
    CACHE[_sc] = (_cb, _sh, _n)
    for _nm in B.ORDERS:
        _err, _tag = note_no_silent_drift(_cb[_nm], _cs[_nm])
        check("sq %dx%dx1 %-9s 构建器 == 规格层 (逐位)" % (_sc[0], _sc[1], _nm),
              _err <= 1e-9,
              "per mag %+.3f,%+.3f,%+.3f  max|d|=%.3g%s"
              % (_cb[_nm][0] / _n, _cb[_nm][1] / _n, _cb[_nm][2] / _n, _err, _tag))

CB4, SH4, N4 = CACHE[(4, 4, 1)]

# --------------------------------------- (B) 4x4x1 手推解析理想值 (第三重锚定)
# 4x4 周期超胞下, 每胞: 2 根 NN(±x 各一被 i<j 去重成 2) / 2 根对角 / 2 根 2a。
# 于是 c_k(per mag) = -(1/2) * sum_{shell k} sigma_i sigma_j / n_mag, 手推得下表。
IDEAL_44 = {
    "fm":       (-2.0, -2.0, -2.0),
    "g_type":   (+2.0, -2.0, -2.0),
    "stripe_a": (0.0, +2.0, -2.0),
    "stripe_b": (0.0, +2.0, -2.0),
    "p4_a":     (-1.0, 0.0, 0.0),
    "mix_ab":   (+1.0, 0.0, 0.0),
}
for _nm, _ideal in IDEAL_44.items():
    _got = tuple(c / N4 for c in CB4[_nm])
    _err = max(abs(g - t) for g, t in zip(_got, _ideal))
    check("sq 4x4x1 %-9s 手推解析理想值 (per mag)" % _nm, _err < 1e-9,
          "got %+.3f,%+.3f,%+.3f  want %+.3f,%+.3f,%+.3f  err=%.3g"
          % (_got + _ideal + (_err,)))

# ------------------------------------------- (C) 结构性简并: 为什么必须用周期 4
_c3_p2 = [CB4[k][2] for k in ("fm", "g_type", "stripe_a", "stripe_b")]
check("sq 4x4x1 周期 2 序 (fm/g_type/stripe_a/stripe_b) 的 c3 完全简并",
      max(_c3_p2) - min(_c3_p2) < 1e-9,
      "c3 = %s  -> J3 只能靠 p4_a/mix_ab 分离" % (["%.3f" % x for x in _c3_p2],))
_ds = max(abs(a - b) for a, b in zip(CB4["stripe_a"], CB4["stripe_b"]))
check("sq 4x4x1 stripe_a == stripe_b (|a|=|b| 下两者必简并)", _ds < 1e-9,
      "max|d|=%.3g" % _ds)
_c3_p4 = [CB4[k][2] for k in ("p4_a", "mix_ab")]
check("sq 4x4x1 周期 4 序的 c3 与周期 2 序**不同** (简并被打破)",
      all(abs(x - _c3_p2[0]) > 1e-9 for x in _c3_p4),
      "p4_a/mix_ab c3 = %s vs 周期 2 的 %.3f" % (["%.3f" % x for x in _c3_p4], _c3_p2[0]))

# ---------------------------------------------------------- (D) 合成 J 反解闭合
A = np.array([[1.0] + CB4[k] for k in B.SELECTED])
for _i, J_true in enumerate(((5.0, -0.5, 0.2), (-111.8, 11.4, -2.0), (0.0, 3.3, -7.7))):
    Jt = [j * 0.5 ** 2 for j in J_true]
    E = np.array([sum(CB4[k][m] * Jt[m] for m in range(3))
                  for k in B.SELECTED]) / 1000.0
    sol = np.linalg.solve(A, E)
    got = [sol[1 + m] / 0.5 ** 2 * 1000.0 for m in range(3)]
    check("sq 4x4x1 合成 J%d=%s 反解闭合" % (_i + 1, J_true),
          all(abs(g - t) < 1e-6 for g, t in zip(got, J_true)),
          "-> %s" % (["%.6f" % g for g in got],))

# -------------------------------------------------------------- (E) 条件数闸门
_cond = float(np.linalg.cond(A))
check("sq 4x4x1 cond(A) < 1e4 闸门", _cond < 1e4, "cond(A)=%.4f" % _cond)

# --------------------------- (F) K2NiF4 同构: 系数只由几何决定, 与 a 无关
CBK, SHK, NK = builder_coeffs(POSCAR_K, "Ni", ("K", "F"), (4, 4, 1))
CSK = spec_coeffs(POSCAR_K, "Ni", ("K", "F"), (4, 4, 1), SHK, 1.0)
_errK = max(note_no_silent_drift(CBK[nm], CSK[nm])[0] for nm in B.ORDERS)
check("sq K2NiF4 4x4x1 构建器 == 规格层 (逐位)", _errK <= 1e-9, "max|d|=%.3g" % _errK)
_same = all(abs(CBK[nm][m] - CB4[nm][m]) < 1e-9 for nm in B.ORDERS for m in range(3))
check("sq K2NiF4 与 La2CuO4 的系数逐位相同 (纯几何设计: 系数与 a 无关)", _same)

print("\n全部通过" if not fails else "\n%d 个失败: %s" % (len(fails), fails))
sys.exit(1 if fails else 0)
