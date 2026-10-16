#!/usr/bin/env python
"""test_crect_cases.py — C 心矩形多亚晶格 (CrOCl 类) builder 的契约测试

覆盖:
  1. **两套独立实现逐位一致**: build_crect_cases 自带几何 vs P0.1 规格层;
  2. 7 个磁序在 4x4x1 下的手推解析理想值 (第三重锚定);
  3. **basis 简并是陷阱**: J1 = C 心矢量 (1/2,1/2) 连的是**同一原胞内**的两个 Cr,
     所以 basis 符号对一切 (n0,n1)-奇偶序完全不可见 (afm_c_a == afm_a, 逐位相同)。
     真陷阱是: 一个序与它的 basis 翻倍版**同时选** => 设计矩阵两行完全相同 => det 精确为 0。
     afm_c 的价值在于它是唯一翻转 basis 乘积的选定序, 把 J1 列从奇偶序的 0 变成 ±64;
     但**纯奇偶序组其实也能解** ((fm,afm_a,afm_b,afm_diag): det=-262144, cond=59.17),
     所以 afm_c 不是数学必需 (本测试把这两条都钉住, 免得理由被夸大);
  4. 壳容差陷阱: tol=0.2 会把 3.0347 与 3.1980 两个物理不同的壳**静默并成一个**
     (间距只有 0.1632) —— 默认 0.05 才是对的;
  5. per_site_bonds = -2*c(fm)/n = (4,2,2), 且与超胞大小无关 (跨口径不变量,
     驱动层用它校验 poscar_to_mc3 的 2D 投影没把 J 绑到错壳上);
  6. **两把尺子**: J1 的 3D 键长 3.0347 A (含 ~1.71 A 层起伏) 与面内投影
     2.5064 A 都是事实; 选壳必须用 3D, 面内值不得出现在壳表里;
  7. 守卫: 胞内磁性原子数 != 2 / 超胞前两维 < 2 必须显式报错;
  8. 落盘 fixture 与现场重算逐位一致 (产物可复现);
  9. `--model` 预设自检 (driver 的 `check_presets`): 预设与 `literature_J.json` 逐项一致,
     并反向验证「把预设改脏后自检真的会报」——免得自检永远是空列表还不被发现。

运行:
  python tests/test_crect_cases.py
  pytest tests/test_crect_cases.py -q
"""
import contextlib
import io
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(HERE)
SCRIPTS = os.path.join(HERE, "scripts")
FIXTURES = os.environ.get("MAGNET_FIXTURES", os.path.join(HERE, "tests", "fixtures"))
POSCAR_CR = os.path.join(FIXTURES, "CrOCl_POSCAR")
POSCAR_LA = os.path.join(FIXTURES, "La2CuO4_POSCAR")

if not os.path.exists(POSCAR_CR):
    sys.exit("[skip] 缺少 fixture: %s (仓库应自带 tests/fixtures/, 或设 MAGNET_FIXTURES)"
             % POSCAR_CR)

sys.path.insert(0, SCRIPTS)

import build_crect_cases as B  # noqa: E402
import magnetic_orders as MO  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + (("   " + extra) if extra else ""))
    if not cond:
        fails.append(name)


# ---------------------------------------------------------------- 两套独立实现

def _geometry(poscar, element, sc, shell_tol=0.05):
    """(1) build_crect_cases 自带的几何实现 —— 独立基准。"""
    lat, atoms, coords, _mode = B.parse_poscar(poscar)
    n_prim = len(atoms)
    cr_prim = [i for i, a in enumerate(atoms) if a == element]
    lat_sc, atoms_sc, coords_sc, cell_idx = B.make_supercell(lat, atoms, coords, list(sc))
    mag = [i for i, a in enumerate(atoms_sc) if a == element]
    n = len(mag)
    prim_of_mag = np.array([i % n_prim for i in mag], dtype=int)
    basis = np.where(prim_of_mag == cr_prim[0], 1.0, -1.0)
    shells = B.magnetic_shells(mag, coords_sc, lat_sc, max_shells=3, tol=shell_tol)
    n0 = cell_idx[mag][:, 0]
    n1 = cell_idx[mag][:, 1]
    sig = B.sigmas_of(n0, n1, basis)
    sums = {k: B.classify_bonds(mag, coords_sc, lat_sc, s, shells, tol=shell_tol)
            for k, s in sig.items()}
    coeff = {k: [-0.5 * s[0], -0.5 * s[1], -0.5 * s[2]] for k, s in sums.items()}
    return coeff, shells, n, basis


def builder_coeffs(poscar, element, sc, shell_tol=0.05):
    coeff, shells, n, _basis = _geometry(poscar, element, sc, shell_tol)
    return coeff, shells, n


def spec_coeffs(poscar, element, nonmag, sc, shells, spin, basis_prim, shell_tol=0.05):
    """(2) P0.1 共享规格层 —— 与 (1) 完全独立的第二套实现。"""
    spec = MO.Spec(
        key="crect_test", label="C 心矩形 (test)",
        structure=MO.StructureSpec(poscar=os.path.abspath(poscar), element=element,
                                   nonmag=tuple(nonmag), supercell=tuple(sc)),
        shells=tuple(MO.Shell("J%d" % (i + 1), d, tol=shell_tol)
                     for i, d in enumerate(shells)),
        orders=tuple(MO.Order(nm, (lambda c, nm=nm: B.sigmas_of(
            c.n0, c.n1, np.where(c.prim_idx[c.mag_idx] == basis_prim, 1.0, -1.0))[nm]))
            for nm in B.ORDERS),
        selected=tuple(B.SELECTED), spin=spin,
        # 必须显式给 shift_search: 规格层默认 min_image, C 心矢量是半周期。
        policy="shift_search", nshift=1,
    )
    res, _ctx = MO.compute(spec)
    return {nm: [float(x) for x in res[nm].coeffs] for nm in B.ORDERS}


def note_no_silent_drift(a, b):
    d = np.asarray(a) - np.asarray(b)
    m = float(np.max(np.abs(d)))
    tag = ""
    if m > 1e-9 and np.allclose(d, -0.5 * np.asarray(b), atol=1e-6):
        tag = "  <== 典型症状: 半周期平移像被折叠, 该壳的键少算一半"
    return m, tag


# ------------------------------------------------- (A) 超胞 x 磁序 逐位交叉校验
CACHE = {}
for _sc in ((4, 4, 1), (8, 8, 1)):
    _cb, _sh, _n = builder_coeffs(POSCAR_CR, "Cr", _sc)
    _cs = spec_coeffs(POSCAR_CR, "Cr", ("O", "Cl"), _sc, _sh, 1.5, 0)
    CACHE[_sc] = (_cb, _sh, _n)
    for _nm in B.ORDERS:
        _err, _tag = note_no_silent_drift(_cb[_nm], _cs[_nm])
        check("cr %dx%dx1 %-11s 构建器 == 规格层 (逐位)" % (_sc[0], _sc[1], _nm),
              _err <= 1e-9,
              "per mag %+.3f,%+.3f,%+.3f  max|d|=%.3g%s"
              % (_cb[_nm][0] / _n, _cb[_nm][1] / _n, _cb[_nm][2] / _n, _err, _tag))

CB4, SH4, N4 = CACHE[(4, 4, 1)]

# --------------------------------------- (B) 4x4x1 手推解析理想值 (第三重锚定)
# 每胞 J1 壳 4 根键 (两个 Cr 各自的 (1/2,1/2) 与 (1/2,-1/2), i<j 去重后 4)、
# J2 壳 2 根 (a 方向) 、J3 壳 2 根 (b 方向)。
# c_k(per mag) = -(1/2) * sum_{shell k} sigma_i sigma_j / n_mag, 手推得下表。
IDEAL_44 = {
    "fm":           (-2.0, -1.0, -1.0),
    "afm_c":        (+2.0, -1.0, -1.0),
    "afm_a":        (0.0, +1.0, -1.0),
    "afm_b":        (0.0, -1.0, +1.0),
    "afm_diag":     (0.0, +1.0, +1.0),
    "afm_c_a":      (0.0, +1.0, -1.0),
    "afm_c_diag":   (0.0, +1.0, +1.0),
}
for _nm, _ideal in IDEAL_44.items():
    _got = tuple(c / N4 for c in CB4[_nm])
    _err = max(abs(g - t) for g, t in zip(_got, _ideal))
    check("cr 4x4x1 %-11s 手推解析理想值 (per mag)" % _nm, _err < 1e-9,
          "got %+.3f,%+.3f,%+.3f  want %+.3f,%+.3f,%+.3f  err=%.3g"
          % (_got + _ideal + (_err,)))

# --------------------------- (C) basis 简并: C 心键连的是同一胞内的两个 basis
_da = max(abs(a - b) for a, b in zip(CB4["afm_c_a"], CB4["afm_a"]))
_dd = max(abs(a - b) for a, b in zip(CB4["afm_c_diag"], CB4["afm_diag"]))
check("cr 4x4x1 basis 符号对奇偶序不可见 (afm_c_a == afm_a)", _da < 1e-9, "max|d|=%.3g" % _da)
check("cr 4x4x1 basis 符号对奇偶序不可见 (afm_c_diag == afm_diag)", _dd < 1e-9,
      "max|d|=%.3g" % _dd)
# fm 与 afm_c 只差 c1 (J1 壳), c2/c3 完全相同 —— 这正是"少了 afm_c 就解不出 J1"的原因
_f, _a = CB4["fm"], CB4["afm_c"]
check("cr 4x4x1 fm 与 afm_c 只在 c1 上不同 (c2/c3 逐位相同)",
      abs(_f[1] - _a[1]) < 1e-9 and abs(_f[2] - _a[2]) < 1e-9 and abs(_f[0] - _a[0]) > 1e-9,
      "fm c1=%+.3f vs afm_c c1=%+.3f" % (_f[0], _a[0]))

# fm 只与 afm_c 在 c1 上不同 (c2/c3 逐位相同) —— 这是 J1 符号的**直接**可见后果。
# 但注意: 反过来说"少了 afm_c 就解不出 J1"是**错的** (纯奇偶序组 det=-262144, cond=59.17,
# 因为 J1 连的是 (n0-1,n1) 这样的跨胞 C 心键, 奇偶序照样给出 c1 = 0 这一行的信息)。
# 真正的奇异来自"序 + 它的 basis 翻倍版"同时入选。
for _pair, _nm in ((("fm", "afm_a", "afm_c_a", "afm_b"), "afm_a + afm_c_a"),
                   (("fm", "afm_diag", "afm_c_diag", "afm_b"), "afm_diag + afm_c_diag")):
    _Ab = np.array([[1.0] + CB4[k] for k in _pair])
    _detb = float(np.linalg.det(_Ab))
    check("cr 同时选 %s (逐位同序) 的设计矩阵**奇异**" % _nm, abs(_detb) < 1e-6,
          "det=%.3g" % _detb)
check("cr fm/afm_c 的 c1 反号且绝对值相等 (翻转 basis 乘积的唯一可见后果)",
      abs(_f[0] + _a[0]) < 1e-9, "c1: fm=%+.1f  afm_c=%+.1f" % (_f[0], _a[0]))
check("cr 选定的四个序里只有 fm/afm_c 的 c1 != 0 (另两个恒为 0)",
      sum(1 for k in B.SELECTED if abs(CB4[k][0]) > 1e-9) == 2,
      "c1 = %s" % ([(k, "%.1f" % CB4[k][0]) for k in B.SELECTED],))
_Ap = np.array([[1.0] + CB4[k] for k in ("fm", "afm_a", "afm_b", "afm_diag")])
check("cr 纯奇偶序组 (fm,afm_a,afm_b,afm_diag) **也可解** (afm_c 非数学必需)",
      abs(float(np.linalg.det(_Ap))) > 1e-6 and float(np.linalg.cond(_Ap)) < 1e4,
      "det=%.1f cond=%.4f" % (float(np.linalg.det(_Ap)), float(np.linalg.cond(_Ap))))
_A = np.array([[1.0] + CB4[k] for k in B.SELECTED])
_det = float(np.linalg.det(_A))
_cond = float(np.linalg.cond(_A))
check("cr 选定四个序的设计矩阵非奇异", abs(_det) > 1e-6, "det=%.1f" % _det)
check("cr 4x4x1 cond(A) < 1e4 闸门", _cond < 1e4, "cond(A)=%.4f" % _cond)

# ---------------------------------------------------------- (D) 合成 J 反解闭合
for _i, J_true in enumerate(((5.0, -0.5, 0.2), (12.0, 3.3, -7.7), (-8.0, 1.0, 1.0))):
    Jt = [j * 1.5 ** 2 for j in J_true]
    E = np.array([sum(CB4[k][m] * Jt[m] for m in range(3))
                  for k in B.SELECTED]) / 1000.0
    sol = np.linalg.solve(_A, E)
    got = [sol[1 + m] / 1.5 ** 2 * 1000.0 for m in range(3)]
    check("cr 4x4x1 合成 J%d=%s 反解闭合" % (_i + 1, J_true),
          all(abs(g - t) < 1e-6 for g, t in zip(got, J_true)),
          "-> %s" % (["%.6f" % g for g in got],))

# --------------------------------------------------------------- (E) 壳容差陷阱
_lat, _atoms, _coords, _m = B.parse_poscar(POSCAR_CR)
_lsc, _asc, _csc, _ci = B.make_supercell(_lat, _atoms, _coords, [4, 4, 1])
_mag = [i for i, a in enumerate(_asc) if a == "Cr"]
_sh005 = B.magnetic_shells(_mag, _csc, _lsc, max_shells=3, tol=0.05)
_sh020 = B.magnetic_shells(_mag, _csc, _lsc, max_shells=3, tol=0.2)
check("cr J1/J2 是两个物理不同的壳 (3.0347 vs 3.1980, 间距 0.1632 < 0.2)",
      abs(_sh005[0] - 3.03470713) < 1e-6 and abs(_sh005[1] - 3.19795066) < 1e-6
      and 0.0 < (_sh005[1] - _sh005[0]) < 0.2,
      "tol=0.05 -> %s" % (["%.4f" % d for d in _sh005],))
check("cr tol=0.2 会把 J1/J2 **静默合并** (第二个壳跳到 3.8603 = J3)",
      abs(_sh020[1] - 3.86026201) < 1e-6,
      "tol=0.2 -> %s" % (["%.4f" % d for d in _sh020],))
_cb020, _sh020b, _n020 = builder_coeffs(POSCAR_CR, "Cr", (4, 4, 1), shell_tol=0.2)
_drift = max(abs(a - b) for k in B.ORDERS for a, b in zip(_cb020[k], CB4[k]))
check("cr tol=0.2 会静默改变系数 (错绑 J2/J3) —— 所以默认必须是 0.05",
      _drift > 1e-9, "max|d| vs tol=0.05 = %.3g" % _drift)

# ------------------------------------------- (F) 每格点键数 (跨口径不变量)
_per_site = tuple(-2.0 * c / N4 for c in CB4["fm"])
check("cr per_site_bonds = (4, 2, 2)", all(abs(g - t) < 1e-9 for g, t in
                                          zip(_per_site, (4.0, 2.0, 2.0))),
      "-> (%.1f, %.1f, %.1f)" % _per_site)
CB8, SH8, N8 = CACHE[(8, 8, 1)]
_per_site8 = tuple(-2.0 * c / N8 for c in CB8["fm"])
check("cr per_site_bonds 与超胞大小无关 (4x4x1 == 8x8x1)",
      all(abs(a - b) < 1e-9 for a, b in zip(_per_site, _per_site8)),
      "4x4 -> %s  8x8 -> %s" % (["%.1f" % x for x in _per_site],
                                ["%.1f" % x for x in _per_site8]))
_sh8 = CACHE[(8, 8, 1)][1]
check("cr 壳距离也与超胞大小无关",
      all(abs(a - b) < 1e-9 for a, b in zip(SH4, _sh8)),
      "%s vs %s" % (["%.4f" % d for d in SH4], ["%.4f" % d for d in _sh8]))

# --------------------------------------------------------- (G) 两把尺子 (3D vs 面内)
_d_frac = np.array([0.5, 0.5, 0.457453116420 - 0.542571911159])
_cart = _d_frac @ _lat
_d3 = float(np.linalg.norm(_cart))
_d2 = float(np.linalg.norm(_cart[:2]))
_dz = abs(float(_cart[2]))
check("cr J1 的 3D 键长 = 3.0347 A (含层起伏)", abs(_d3 - SH4[0]) < 1e-9,
      "|d|_3D=%.4f  vs 壳表 J1=%.4f" % (_d3, SH4[0]))
check("cr 同一对 Cr 的**面内投影** = 2.5064 A (两把尺子差 %.3f A)" % _dz,
      abs(_d2 - 2.5064) < 5e-4, "|d|_2D=%.4f  |dz|=%.4f" % (_d2, _dz))
check("cr 面内值 2.5064 **不得**出现在壳表里 (选壳必须用 3D 口径)",
      all(abs(_d2 - s) > 1e-3 for s in SH4), "壳表 %s" % (["%.4f" % d for d in SH4],))

# ------------------------------------------------------------- (H) 守卫 / 报错
def _run_main(argv):
    """跑 build_crect_cases.main(), 返回 SystemExit.code (0 = 没退出)。"""
    old = sys.argv
    sys.argv = ["build_crect_cases.py"] + argv
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            B.main()
    except SystemExit as exc:
        return exc.code
    finally:
        sys.argv = old
    return 0


if os.path.exists(POSCAR_LA):
    _code = _run_main([POSCAR_LA, "--element", "Cu", "--nonmag", "La", "O", "--selftest"])
    check("cr 胞内 1 个磁性原子时显式报错 (提示改用 build_rectangular_cases.py)",
          isinstance(_code, str) and "2 个磁性原子/胞" in _code, "code=%r" % (_code,))
else:
    check("cr 胞内 1 个磁性原子时显式报错 (跳过: 缺 La2CuO4 fixture)", True)
_code = _run_main([POSCAR_CR, "--element", "Cr", "--supercell", "1", "4", "1", "--selftest"])
check("cr 超胞前两维 < 2 时显式报错 (C 心键是半周期)",
      isinstance(_code, str) and "必须 >= 2" in _code, "code=%r" % (_code,))

# ------------------------------------------------- (I) 落盘 fixture 可复现
_fx = os.path.join(ROOT, "CrOCl_pipeline", "02_supercell", "cases", "magnetic_orders.json")
if os.path.exists(_fx):
    with open(_fx, encoding="utf-8") as _fh:
        _d = json.load(_fh)
    _got = {o["name"]: o["coeffs"] for o in _d["order_spec"]}
    _worst = max(max(abs(a - b) for a, b in zip(_got[nm], CB4[nm])) for nm in B.ORDERS)
    check("cr 落盘 magnetic_orders.json 与现场重算逐位一致", _worst <= 1e-9,
          "max|d|=%.3g  (%s)" % (_worst, os.path.relpath(_fx, ROOT)))
    _ps = _d.get("per_site_bonds") or {}
    _ok = all(abs(float(_ps.get("J%d" % (i + 1), 1e9)) - t) < 1e-9
              for i, t in enumerate((4.0, 2.0, 2.0)))
    check("cr 落盘 per_site_bonds = {J1:4, J2:2, J3:2} (驱动层据此校验 MC 模型)",
          _ok, "per_site_bonds=%s" % (_ps,))
else:
    check("cr 落盘 fixture 存在 (跳过: %s 不在)" % _fx, True)

# ------------------------------------------- (J) --model 预设自检 (driver 侧)
import importlib.util as _ilu
_dspec = _ilu.spec_from_file_location("_drv", os.path.join(HERE, "run_2d_magnet.py"))
_drv = _ilu.module_from_spec(_dspec)
_dspec.loader.exec_module(_drv)
_probs = _drv.check_presets()
check("cr --model 预设自检: 预设的文献锚点全部与 literature_J.json 一致",
      _probs == [], "problems=%s" % (_probs,))
# 反向验证: 把预设改脏后自检必须报出来 (否则这条自检可能永远返回空列表而不被发现)
_p = _drv.MODEL_PRESETS["crocl"]
_old_u = _p["U"]
try:
    _p["U"] = _old_u + 1.0
    _bad = _drv.check_presets()
finally:
    _p["U"] = _old_u
check("cr --model 预设自检能抓漂移 (改脏 U 后必须报出 crocl)",
      any("crocl" in s and "U" in s for s in _bad), "problems=%s" % (_bad,))

print("\n全部通过" if not fails else "\n%d 个失败: %s" % (len(fails), fails))
sys.exit(1 if fails else 0)
