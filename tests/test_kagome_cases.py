#!/usr/bin/env python
"""test_kagome_cases.py — kagome (Fe3Sn 类, 3 磁位/胞) builder + 反解器的契约测试

覆盖:
  1. **两套独立实现逐位一致**: build_kagome_cases 自带几何 (classify_bonds, 显式平移)
     vs P0.1 规格层 (magnetic_orders.compute, policy="shift_search"), 两个超胞各 7 个序;
  2. 3x3x1 下 7 个磁序的**手推解析理想值** (per mag, 第三重锚定);
  3. **q=0 简并是核心陷阱**: q=0 一族 (fm/q0_120/q0_ud/q0_duu) 的 c1 == c2 逐位相同,
     所以那一族的设计矩阵 rank 只有 2 —— J1/J2 在 q=0 子块里**不可分离**; 只有加入
     √3×√3 调制 (sq3_*) 才把 rank 提到 4。本测试把「q=0 无解」与「selected 可解」
     两条都钉住, 免得"必须上 3 倍胞"被说成玄学;
  4. √3×√3 的三个序给出的信息互不相同 (sq3_updntri 的 c3 ≠ q0 族; sq3_alt 的 c1 ≠ c2):
     这是 J1/J2/J3 能被分离的直接来源; 并钉住 √3×√3 标签 (n0+n1)%3 的良定义性;
  5. per_site_bonds = (4, 4, 6) 且**与超胞大小无关** (跨口径不变量, 驱动层/poscar_to_mc3
     据此校验 J 没被绑到错壳上);
  6. 守卫: 超胞非 3 的倍数 / 胞内非 3 个磁性原子 / 非 kagome 的 3 原子胞 (Fe3GeTe2)
     必须显式报错; 反解器遇到 rank<4 或 cond 病态也必须报错 (**不静默出 J**);
  7. **无 DFT 的闭环**: 用系数矩阵合成 4 个磁态能量写进 mock OUTCAR, 反解器必须把
     J1/J2/J3 逐位还原 (M3/mock 口径, 不下发 DFT);
  8. 落盘 fixture 可复现 + 驱动侧集成 (detect_lattice 认识 kagome 且不误伤
     Fe3GeTe2/Nb3I8; ERROR_CATALOG 有 kagome_degenerate; --model fe3sn 预设自检)。

运行:
  python tests/test_kagome_cases.py
  pytest tests/test_kagome_cases.py -q
"""
import contextlib
import io
import json
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.dirname(HERE)
SCRIPTS = os.path.join(HERE, "scripts")
FIXTURES = os.environ.get("MAGNET_FIXTURES", os.path.join(HERE, "tests", "fixtures"))
POSCAR_FE3SN = os.path.join(FIXTURES, "Fe3Sn_POSCAR")
POSCAR_FGT = os.path.join(FIXTURES, "Fe3GeTe2_monolayer_POSCAR")
POSCAR_NB3I8 = os.path.join(FIXTURES, "Nb3I8_POSCAR")
POSCAR_LA = os.path.join(FIXTURES, "La2CuO4_POSCAR")

if not os.path.exists(POSCAR_FE3SN):
    sys.exit("[skip] 缺少 fixture: %s (仓库应自带 tests/fixtures/, 或设 MAGNET_FIXTURES)"
             % POSCAR_FE3SN)

sys.path.insert(0, SCRIPTS)

# 控制台编码兜底: 本机默认 cp936, 而断言的 `extra` 里会出现 ℓ / ⇒ 等非 GBK 字符
# (builder/solver 的报错原文就带 ℓ)。只放宽 errors, 不改编码 —— 与 run_2d_magnet.py 顶部同款。
for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors="replace")
    except Exception:
        pass

import build_kagome_cases as B  # noqa: E402
import ingest_dft as IG  # noqa: E402
import magnetic_orders as MO  # noqa: E402
import solve_J_kagome as S  # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + (("   " + extra) if extra else ""))
    if not cond:
        fails.append(name)


# ---------------------------------------------------------------- 两套独立实现

def _geometry(poscar, element, sc, shell_tol=0.2):
    """(1) build_kagome_cases 自带的几何实现 —— 独立基准。"""
    lat, atoms, coords, _mode = B.parse_poscar(poscar)
    sc = list(sc)
    sc_f = np.array(sc, dtype=float)
    lat_sc, atoms_sc, coords_sc, cell_idx = B.make_supercell(lat, atoms, coords, sc)
    _idxs, sub_map, _n_sub = B.sublattice_labels(atoms, coords, element)
    mag = [i for i, a in enumerate(atoms_sc) if a == element]
    sub = np.array([sub_map[tuple(np.round((coords_sc[i] * sc_f) % 1.0, 6))] for i in mag],
                   dtype=int)
    n0 = cell_idx[mag][:, 0].astype(int)
    n1 = cell_idx[mag][:, 1].astype(int)
    shells = B.magnetic_shells(mag, coords_sc, lat_sc, max_shells=3, tol=shell_tol)
    spins = B.spins_of(n0, n1, sub)
    coeff = {}
    for nm in B.ORDERS:
        sums = B.classify_bonds(mag, coords_sc, lat_sc, spins[nm], shells)
        coeff[nm] = [-0.5 * s for s in sums]
    return coeff, shells, len(mag)


def spec_coeffs(poscar, element, nonmag, sc, shells, spin=1.0, shell_tol=0.2):
    """(2) P0.1 共享规格层 —— 与 (1) 完全独立的第二套实现。

    序的构造用**规格层自己的**相位助手 (sub_label 由 magnetic_orders 从分数坐标推),
    所以这一节同时钉住「两套亚晶格标签推导同序」这件事。
    """
    spec = MO.Spec(
        key="kagome_test", label="kagome (test)",
        structure=MO.StructureSpec(poscar=os.path.abspath(poscar), element=element,
                                   nonmag=tuple(nonmag), supercell=tuple(sc),
                                   sublattice_names=("A", "B", "C")),
        shells=tuple(MO.Shell("J%d" % (i + 1), d, tol=shell_tol)
                     for i, d in enumerate(shells)),
        orders=tuple(MO.Order(nm, (lambda c, nm=nm: B.spins_of(c.n0, c.n1, c.sub_label)[nm]))
                     for nm in B.ORDERS),
        selected=tuple(B.SELECTED), spin=spin,
        # 必须显式给 shift_search: kagome 的 √3×√3 序依赖胞相位, 且 min_image 在
        # 小超胞会把远壳的像折叠到近处 (与 build_* 的显式平移搜索同口径)。
        policy="shift_search", nshift=1,
    )
    res, _ctx = MO.compute(spec)
    return {nm: [float(x) for x in res[nm].coeffs] for nm in B.ORDERS}


# ------------------------------------------------- (A) 超胞 x 磁序 逐位交叉校验
CACHE = {}
for _sc in ((3, 3, 1), (6, 6, 1)):
    _cb, _sh, _n = _geometry(POSCAR_FE3SN, "Fe", _sc)
    _cs = spec_coeffs(POSCAR_FE3SN, "Fe", ("Sn",), _sc, _sh)
    CACHE[_sc] = (_cb, _sh, _n)
    for _nm in B.ORDERS:
        _err = max(abs(a - b) for a, b in zip(_cb[_nm], _cs[_nm]))
        check("kag %dx%dx1 %-12s 构建器 == 规格层 (逐位)" % (_sc[0], _sc[1], _nm),
              _err <= 1e-9,
              "per mag %+.4f,%+.4f,%+.4f  max|d|=%.3g"
              % (_cb[_nm][0] / _n, _cb[_nm][1] / _n, _cb[_nm][2] / _n, _err))

CB3, SH3, N3 = CACHE[(3, 3, 1)]

# --------------------------------------- (B) 3x3x1 手推解析理想值 (第三重锚定)
# kagome: 每胞 3 个磁性位, 前三个壳是 a/2 (J1, 每格点 4 根键)、a√3/2 (J2, 4 根)、
# a (J3, 6 根)。c_k(per mag) = -(1/2) * Σ_shell σ_i σ_j / n_mag, 手推得下表。
for _nm, _ideal in B.IDEAL_33.items():
    _got = tuple(c / N3 for c in CB3[_nm])
    _err = max(abs(g - t) for g, t in zip(_got, _ideal))
    check("kag 3x3x1 %-12s 手推解析理想值 (per mag)" % _nm, _err < 1e-9,
          "got %+.4f,%+.4f,%+.4f  want %+.4f,%+.4f,%+.4f  err=%.3g"
          % (_got + _ideal + (_err,)))

# --------------------------------------------------------- (C) q=0 简并 (核心陷阱)
_q0 = B.Q0_FAMILY
_d12 = max(abs(CB3[k][0] - CB3[k][1]) for k in _q0)
check("kag q=0 一族 c1 == c2 逐位相同 (J1 与 J2 在 q=0 子块里不可分辨)", _d12 < 1e-9,
      "max|c1-c2| = %.3g  (c=%s)" % (_d12, [(k, [round(x, 3) for x in CB3[k]]) for k in _q0]))
_du = max(abs(a - b) for a, b in zip(CB3["q0_ud"], CB3["q0_duu"]))
check("kag q0_ud 与 q0_duu 逐位相同 (纯亚晶格翻倍序对系数不可见)", _du < 1e-9,
      "max|d| = %.3g" % _du)
_Aq0 = np.array([[1.0] + CB3[k] for k in _q0])
_rq0 = int(np.linalg.matrix_rank(_Aq0))
check("kag **q=0 一族的设计矩阵 rank == 2** (四个方程只有两个独立)", _rq0 == 2,
      "rank=%d  cond=%.3g" % (_rq0, float(np.linalg.cond(_Aq0))))
_raised = False
try:
    np.linalg.solve(_Aq0, np.array([1.0, 2.0, 3.0, 4.0]))
except np.linalg.LinAlgError:
    _raised = True
check("kag q=0 一族直接 np.linalg.solve 抛 LinAlgError (不会给出假的 J)", _raised)
_A_sel = np.array([[1.0] + CB3[k] for k in B.SELECTED])
_rk = int(np.linalg.matrix_rank(_A_sel))
_det = float(np.linalg.det(_A_sel))
_cond = float(np.linalg.cond(_A_sel))
check("kag 选定四个序 (含 √3×√3) 的设计矩阵 rank == 4", _rk == 4, "rank=%d" % _rk)
check("kag 选定四个序的设计矩阵非奇异", abs(_det) > 1e-6, "det=%.1f" % _det)
check("kag 选定四个序 cond(A) < 1e4 闸门 (反解噪声放大可控)", _cond < 1e4,
      "cond(A)=%.4f" % _cond)
check("kag 选定四个序里有 √3×√3 调制 (去掉 sq3_* 就退回 rank=2)",
      all(k in B.SELECTED for k in ("sq3_updntri", "sq3_alt"))
      and int(np.linalg.matrix_rank(np.array([[1.0] + CB3[k]
                                              for k in B.SELECTED
                                              if not k.startswith("sq3")]))) == 2,
      "selected=%s" % (B.SELECTED,))
# 备选组合: 不用 fm, 全用 √3×√3 也能解 (rank 4); 全 7 个序同样 rank 4
for _alt, _lbl in (((("q0_120", "sq3_120", "sq3_updntri", "sq3_alt")), "q0_120+三个 sq3"),
                   ((tuple(B.ORDERS)), "全 7 个序")):
    _Aa = np.array([[1.0] + CB3[k] for k in _alt])
    check("kag 备选组合 %s 也可解 (rank 4 且 cond < 1e4)" % _lbl,
          int(np.linalg.matrix_rank(_Aa)) == 4 and float(np.linalg.cond(_Aa)) < 1e4,
          "rank=%d cond=%.4f" % (int(np.linalg.matrix_rank(_Aa)), float(np.linalg.cond(_Aa))))
# √3×√3 给出的**新信息**: sq3_updntri 的 c3 与 q0 族不同 (J3 的可分辨来源),
# sq3_alt 的 c1 != c2 (J1/J2 的可分辨来源)。
check("kag sq3_updntri 的 c3 与 q0 族不同 (√3×√3 才能分辨 J3 的第二类信息)",
      all(abs(CB3["sq3_updntri"][2] - CB3[k][2]) > 1e-9 for k in _q0),
      "sq3_updntri c3=%+.3f  q0_120 c3=%+.3f" % (CB3["sq3_updntri"][2], CB3["q0_120"][2]))
check("kag sq3_alt 的 c1 != c2 (J1 与 J2 分离的直接来源)",
      abs(CB3["sq3_alt"][0] - CB3["sq3_alt"][1]) > 1e-9,
      "c1=%+.4f c2=%+.4f" % (CB3["sq3_alt"][0], CB3["sq3_alt"][1]))
check("kag sq3_120 的 c1 == c2 == 0 但 c3 > 0 (纯 120° 调制对 J1/J2 完全不可见)",
      abs(CB3["sq3_120"][0]) < 1e-9 and abs(CB3["sq3_120"][1]) < 1e-9
      and CB3["sq3_120"][2] > 0.0,
      "c=%s" % ([round(x, 4) for x in CB3["sq3_120"]],))
# √3×√3 标签的良定义性: ℓ = (n0+n1) % 3 在 √3×√3 的两个格矢 (a1-a2, a1+2a2) 下不变
_lab = np.array([(n0 + n1) % 3 for n0 in range(-6, 7) for n1 in range(-6, 7)])
_shift1 = np.array([((n0 + 1) + (n1 - 1)) % 3 for n0 in range(-6, 7) for n1 in range(-6, 7)])
_shift2 = np.array([((n0 + 1) + (n1 + 2)) % 3 for n0 in range(-6, 7) for n1 in range(-6, 7)])
check("kag √3×√3 标签 (n0+n1)%3 在 (a1-a2) 与 (a1+2a2) 平移下不变",
      np.array_equal(_lab, _shift1) and np.array_equal(_lab, _shift2))
check("kag 规格层的 sqrt3_label 与 builder 的 (n0+n1)%3 同序",
      all(MO.sqrt3_label(np.array([a]), np.array([b]))[0] == (a + b) % 3
          for a, b in ((0, 0), (1, 0), (0, 2), (-1, 4), (5, -3))))

# ---------------------------------------------------------- (D) 合成 J 反解闭合
for _i, J_true in enumerate(((1.5, -0.8, 0.4), (12.0, 3.3, -7.7), (-8.0, 1.0, 1.0))):
    Jt = [j * 1.0 ** 2 for j in J_true]
    E = np.array([sum(CB3[k][m] * Jt[m] for m in range(3)) for k in B.SELECTED]) / 1000.0
    sol = np.linalg.solve(_A_sel, E)
    got = [sol[1 + m] / 1.0 ** 2 * 1000.0 for m in range(3)]
    check("kag 3x3x1 合成 J%d=%s 反解闭合" % (_i + 1, J_true),
          all(abs(g - t) < 1e-6 for g, t in zip(got, J_true)),
          "-> %s" % (["%.6f" % g for g in got],))

# ------------------------------------------- (E) 每格点键数 (跨口径不变量)
_per_site = tuple(-2.0 * c / N3 for c in CB3["fm"])
check("kag per_site_bonds = (4, 4, 6)", all(abs(g - t) < 1e-9 for g, t in
                                            zip(_per_site, (4.0, 4.0, 6.0))),
      "-> (%.1f, %.1f, %.1f)" % _per_site)
check("kag per_site_bonds 与 builder 的 EXPECT_Z 自洽",
      all(abs(g - t) < 1e-9 for g, t in zip(_per_site, B.EXPECT_Z)),
      "per_site=%s  EXPECT_Z=%s" % (["%.1f" % x for x in _per_site], B.EXPECT_Z))
CB6, SH6, N6 = CACHE[(6, 6, 1)]
_per_site6 = tuple(-2.0 * c / N6 for c in CB6["fm"])
check("kag per_site_bonds 与超胞大小无关 (3x3x1 == 6x6x1)",
      all(abs(a - b) < 1e-9 for a, b in zip(_per_site, _per_site6)),
      "3x3 -> %s  6x6 -> %s" % (["%.1f" % x for x in _per_site],
                                ["%.1f" % x for x in _per_site6]))
check("kag 壳距离也与超胞大小无关",
      all(abs(a - b) < 1e-9 for a, b in zip(SH3, SH6)),
      "%s vs %s" % (["%.4f" % d for d in SH3], ["%.4f" % d for d in SH6]))
check("kag 三个壳 = a/2, a√3/2, a (含 Sn 在六角中心不影响磁壳)",
      abs(SH3[0] - 2.7) < 1e-6 and abs(SH3[1] - 2.7 * np.sqrt(3)) < 1e-6
      and abs(SH3[2] - 5.4) < 1e-6,
      "壳 = %s" % (["%.4f" % d for d in SH3],))

# ------------------------------------------------------------- (F) 守卫 / 报错
def _run_main(argv, mod, name):
    """跑模块的 main(), 返回 SystemExit.code (0 = 没退出)。"""
    old = sys.argv
    sys.argv = [name] + argv
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            mod.main()
    except SystemExit as exc:
        return exc.code
    finally:
        sys.argv = old
    return 0


_code = _run_main([POSCAR_FE3SN, "--element", "Fe", "--supercell", "2", "2", "1"],
                  B, "build_kagome_cases.py")
check("kag 超胞非 3 的倍数时显式报错 (√3×√3 需要 3 倍胞)",
      isinstance(_code, str) and "3 的倍数" in _code, "code=%r" % (_code,))
if os.path.exists(POSCAR_LA):
    _code = _run_main([POSCAR_LA, "--element", "Cu", "--nonmag", "La", "O"],
                      B, "build_kagome_cases.py")
    check("kag 胞内 1 个磁性原子时显式报错 (提示改用别的 builder)",
          isinstance(_code, str) and "3 个/胞" in _code, "code=%r" % (_code,))
else:
    check("kag 胞内 1 个磁性原子时显式报错 (跳过: 缺 La2CuO4 fixture)", True)
if os.path.exists(POSCAR_FGT):
    _code = _run_main([POSCAR_FGT, "--element", "Fe", "--nonmag", "Ge", "Te",
                       "--supercell", "3", "3", "1"], B, "build_kagome_cases.py")
    check("kag 3 个磁性原子但配位数非 (4,4,6) 时显式报错 (Fe3GeTe2 不是 kagome)",
          isinstance(_code, str) and "配位数" in _code, "code=%r" % (_code,))
else:
    check("kag 非 kagome 的 3 原子胞显式报错 (跳过: 缺 Fe3GeTe2 fixture)", True)

# 反解器闸门: 用 rank 只有 2 的序集合作一份脏 fixture, 反解必须报错而不是出 J
_TMP = tempfile.mkdtemp(prefix="_kag_test_")
try:
    _cases = os.path.join(_TMP, "cases")
    _code = _run_main([POSCAR_FE3SN, "--element", "Fe", "--nonmag", "Sn", "--moment", "2.0",
                       "--spin", "1.0", "--supercell", "3", "3", "1", "--outdir", _cases],
                      B, "build_kagome_cases.py")
    _jf = os.path.join(_cases, "magnetic_orders.json")
    check("kag 构建器能在临时目录落盘工程 (无 DFT)", os.path.exists(_jf), "%s" % _jf)
    with open(_jf, encoding="utf-8") as _fh:
        _data = json.load(_fh)
    _coeff = {o["name"]: o["coeffs"] for o in _data["order_spec"]}
    _S2 = float(_data["S2"])
    # (7) mock OUTCAR 闭环: E = Σ c_k · (J_k S^2) / 1000 (meV -> eV), 4 个选定序
    # 关键: 写进 **runs/** 而不是直接写进 cases/ —— 必须经 ingest_dft.py 的硬闸门收回。
    # 反解器现在强制要求 cases/ingest_manifest.json (收回凭证), 手工 cp 的 OUTCAR 会被拒。
    _J_TRUE = [1.5, -0.8, 0.4]
    _runs = os.path.join(_TMP, "runs")
    # 合成 OUTCAR 的晶格矢量必须与 cases/<order>/POSCAR 一致, 否则闸门 8 (结构一致性) 会拒
    _pl = [ln for ln in open(os.path.join(_cases, "fm", "POSCAR"), encoding="utf-8")
           .read().splitlines() if ln.strip()]
    _sc = float(_pl[1].split()[0])
    _lat = [[float(x) * _sc for x in _pl[2 + i].split()[:3]] for i in range(3)]
    for _nm in B.SELECTED:
        _e = sum(_coeff[_nm][m] * (_J_TRUE[m] * _S2) for m in range(3)) / 1000.0
        os.makedirs(os.path.join(_runs, _nm), exist_ok=True)
        IG.synth_outcar(os.path.join(_runs, _nm, "OUTCAR"), _e, steps=1,
                        n_ions=36, counts="27 9", lattice=_lat)
    _ok, _prob, _warn, _st, _mpath = IG.ingest(
        _runs, _cases, orders=list(B.SELECTED), min_kb=0, force=False, quiet=True)
    check("kag mock-OUTCAR 闭环: runs/ 经 ingest 硬闸门收回 cases/", _ok and not _prob,
          "ok=%r problems=%s" % (_ok, _prob))
    check("kag mock-OUTCAR 闭环: 收回后留下 ingest_manifest.json",
          os.path.exists(os.path.join(_cases, "ingest_manifest.json")),
          "%s" % os.path.join(_cases, "ingest_manifest.json"))
    _buf = io.StringIO()
    _old = sys.argv
    sys.argv = ["solve_J_kagome.py", _cases, "--no-jresult"]
    try:
        with contextlib.redirect_stdout(_buf):
            S.main()
        _out = _buf.getvalue()
        _code2 = 0
    except SystemExit as _exc:
        _out = _buf.getvalue()
        _code2 = _exc.code
    finally:
        sys.argv = _old
    check("kag mock-OUTCAR 闭环: 反解器退出 0", _code2 == 0, "code=%r" % (_code2,))
    check("kag mock-OUTCAR 闭环: J1/J2/J3 逐位还原 (容差 1e-4 meV)",
          all(("J%d = %+.4f meV" % (i + 1, _J_TRUE[i])) in _out for i in range(3)),
          "expect %s" % (["J%d = %+.4f meV" % (i + 1, _J_TRUE[i]) for i in range(3)],))
    check("kag 反解器读到了收回凭证 (打印 manifest 路径)", "收回凭证" in _out)

    # 凭证闸门: 把 manifest 拿掉, 反解器必须**拒绝**出 J (而不是静默放行)
    _noprov = os.path.join(_TMP, "noprov")
    shutil.copytree(_cases, _noprov)
    os.remove(os.path.join(_noprov, "ingest_manifest.json"))
    _buf3 = io.StringIO()
    sys.argv = ["solve_J_kagome.py", _noprov, "--no-jresult"]
    try:
        with contextlib.redirect_stdout(_buf3):
            S.main()
        _code3 = 0
    except SystemExit as _exc:
        _code3 = _exc.code
    finally:
        sys.argv = _old
    check("kag 缺收回凭证时反解器拒绝出 J (provenance 闸门)", _code3 not in (0, None),
          "code=%r" % (_code3,))

    # 脏 fixture: 把 selected_orders 换成 q=0 一族 (rank 2) ⇒ 反解器必须报错
    _dirty = os.path.join(_TMP, "dirty")
    shutil.copytree(_cases, _dirty)
    with open(os.path.join(_dirty, "magnetic_orders.json"), encoding="utf-8") as _fh:
        _d2 = json.load(_fh)
    _d2["selected_orders"] = list(B.Q0_FAMILY)
    with open(os.path.join(_dirty, "magnetic_orders.json"), "w", encoding="utf-8") as _fh:
        json.dump(_d2, _fh, ensure_ascii=False, indent=1)
    _code3 = _run_main([_dirty, "--no-jresult"], S, "solve_J_kagome.py")
    check("kag 反解器在 rank<4 的序集合上显式报错 (不静默出 J)",
          isinstance(_code3, str) and ("rank" in _code3 or "无法唯一确定" in _code3),
          "code=%r" % (_code3,))
finally:
    shutil.rmtree(_TMP, ignore_errors=True)

# ------------------------------------------------- (G) 落盘 fixture 可复现
_fx = os.path.join(ROOT, "Fe3Sn_pipeline", "02_supercell", "cases", "magnetic_orders.json")
if os.path.exists(_fx):
    with open(_fx, encoding="utf-8") as _fh:
        _d = json.load(_fh)
    _got = {o["name"]: o["coeffs"] for o in _d["order_spec"]}
    _worst = max(max(abs(a - b) for a, b in zip(_got[nm], CB3[nm]))
                 for nm in B.ORDERS if nm in _got)
    check("kag 落盘 magnetic_orders.json 与现场重算逐位一致", _worst <= 1e-9,
          "max|d|=%.3g  (%s)" % (_worst, os.path.relpath(_fx, ROOT)))
    _ps = _d.get("per_site_bonds") or {}
    _ok = all(abs(float(_ps.get("J%d" % (i + 1), 1e9)) - t) < 1e-9
              for i, t in enumerate((4.0, 4.0, 6.0)))
    check("kag 落盘 per_site_bonds = {J1:4, J2:4, J3:6}", _ok, "per_site_bonds=%s" % (_ps,))
    check("kag 落盘 selected_orders == builder 的 SELECTED",
          tuple(_d.get("selected_orders") or ()) == tuple(B.SELECTED),
          "selected=%s" % (_d.get("selected_orders"),))
    check("kag 落盘默认全共线 (spin_mode=collinear, noncollinear=False, 无 120° 序)",
          (not bool(_d.get("noncollinear")))
          and str(_d.get("spin_mode") or "").startswith("collinear")
          and not any("120" in str(o) for o in (_d.get("all_orders") or ()))
          and not any(os.path.isdir(os.path.join(os.path.dirname(_fx), o))
                      and "120" in o for o in ("q0_120", "sq3_120")),
          "noncollinear=%r spin_mode=%r all=%s"
          % (_d.get("noncollinear"), _d.get("spin_mode"), _d.get("all_orders")))
else:
    check("kag 落盘 fixture 存在 (跳过: %s 不在)" % _fx, True)

# ------------------------------------------- (H) 驱动侧集成 (detect_lattice / 预设)
import importlib.util as _ilu  # noqa: E402

_dspec = _ilu.spec_from_file_location("_drv", os.path.join(HERE, "run_2d_magnet.py"))
_drv = _ilu.module_from_spec(_dspec)
_dspec.loader.exec_module(_drv)

_lat_fe3sn = _drv.detect_lattice(POSCAR_FE3SN, "Fe")
check("drv detect_lattice(Fe3Sn) == ('kagome', 1, 3)", _lat_fe3sn == ("kagome", 1, 3),
      "-> %s" % (_lat_fe3sn,))
for _p, _el, _tag in ((POSCAR_FGT, "Fe", "Fe3GeTe2"), (POSCAR_NB3I8, "Nb", "Nb3I8")):
    if not os.path.exists(_p):
        check("drv detect_lattice 不误伤 %s (跳过: 缺 fixture)" % _tag, True)
        continue
    _l = _drv.detect_lattice(_p, _el)
    check("drv detect_lattice 不误伤 %s (仍是 cluster)" % _tag, _l[0] != "kagome",
          "-> %s" % (_l,))
# 配位数判据本身: kagome z=(4,4,4); Fe3GeTe2 z=(4,6,4) (混合配位); Nb3I8 z=(2,2,2) (三聚体)
def _coord(poscar, element):
    lat, atoms, coords, _m = B.parse_poscar(poscar)
    mag = [i for i, a in enumerate(atoms) if a == element]
    return _drv._coordination_z(coords, mag, lat)


_d, _z = _coord(POSCAR_FE3SN, "Fe")
check("drv 配位数判据: Fe3Sn z == (4,4,4) 且 d_min = 2.7000 A",
      list(_z) == [4, 4, 4] and abs(_d - 2.7) < 1e-6, "d_min=%.4f z=%s" % (_d, _z))
if os.path.exists(POSCAR_FGT):
    _d2, _z2 = _coord(POSCAR_FGT, "Fe")
    check("drv 配位数判据: Fe3GeTe2 z == (4,6,4) (混合配位 ⇒ 不是 kagome)",
          list(_z2) == [4, 6, 4], "d_min=%.4f z=%s" % (_d2, _z2))
if os.path.exists(POSCAR_NB3I8):
    _d3, _z3 = _coord(POSCAR_NB3I8, "Nb")
    check("drv 配位数判据: Nb3I8 z == (2,2,2) (三聚体 ⇒ 不是 kagome)",
          list(_z3) == [2, 2, 2], "d_min=%.4f z=%s" % (_d3, _z3))
check("drv LATTICE_OF 把 kagome_fe3sn 映射到 'kagome'",
      MO.LATTICE_OF.get("kagome_fe3sn") == "kagome",
      "-> %r" % (MO.LATTICE_OF.get("kagome_fe3sn"),))
_cat = {e["key"]: e for e in _drv.ERROR_CATALOG}
check("drv ERROR_CATALOG 有 kagome_degenerate (静态, 含 cond/rank 症状)",
      "kagome_degenerate" in _cat and bool(_cat["kagome_degenerate"].get("static"))
      and "rank" in _cat["kagome_degenerate"]["symptom"],
      "keys=%d" % len(_cat))
_probs = _drv.check_presets()
check("drv --model 预设自检: 全部预设与 literature_J.json 一致", _probs == [],
      "problems=%s" % (_probs,))
check("drv --model fe3sn 预设声明 lattice=kagome 且**不谎称**有文献 J",
      _drv.MODEL_PRESETS.get("fe3sn", {}).get("lattice") == "kagome"
      and _drv.MODEL_PRESETS.get("fe3sn", {}).get("literature") is None,
      "preset=%s" % ({k: _drv.MODEL_PRESETS.get("fe3sn", {}).get(k)
                      for k in ("lattice", "literature", "S")},))
# 反向验证: 自检不是空列表恒真
_p = _drv.MODEL_PRESETS["cri3"]
_old_u = _p["U"]
try:
    _p["U"] = _old_u + 1.0
    _bad = _drv.check_presets()
finally:
    _p["U"] = _old_u
check("drv --model 预设自检能抓漂移 (改脏 cri3 的 U 后必须报出)",
      any("cri3" in s and "U" in s for s in _bad), "problems=%s" % (_bad,))
# 接线存在性: build/solve 分派里真的调到了 kagome 的两个脚本
_src = open(os.path.join(HERE, "run_2d_magnet.py"), encoding="utf-8").read()
check("drv 源码里有 kagome 的 build/solve 分派",
      "build_kagome_cases.py" in _src and "solve_J_kagome.py" in _src
      and '"kagome"' in _src)
_buf_h = io.StringIO()
_old_argv = sys.argv
try:
    sys.argv = ["run_2d_magnet.py", "--help"]
    with contextlib.redirect_stdout(_buf_h):
        try:
            _drv.main()
        except SystemExit:
            pass
finally:
    sys.argv = _old_argv
_hlp = _buf_h.getvalue()
_lseg = _hlp.split("--lattice")[1].split("]")[0] if "--lattice" in _hlp else ""
check("drv --help 的 --lattice choices 含全部取值 (含 kagome)",
      all(x in _lseg for x in ("auto", "honeycomb", "triangular", "cluster",
                               "rectangular", "square", "rect_centered", "kagome")))

print("\n全部通过" if not fails else "\n%d 个失败: %s" % (len(fails), fails))
sys.exit(1 if fails else 0)
