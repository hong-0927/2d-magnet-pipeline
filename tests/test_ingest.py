#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""`scripts/ingest_dft.py` 的契约测试 —— P4.2 收回通道 (无 VASP 依赖)。

为什么单独一个文件: `ingest_dft.py --selftest` 覆盖的是"闸门本身对不对"; 这里覆盖的是
**契约**: 哪些闸门是硬的 (不给 --lenient 时不许放行)、原子性 (任何硬错都不许落半个盘)、
`--lenient` 真的生效、manifest 的结构能被下游 (`solve_J_kagome.py` 的凭证闸门) 读懂。

关键: 全部用**合成 OUTCAR** (IG.synth_outcar), 不碰真实 DFT。

    python tests/test_ingest.py
"""
import json
import os
import shutil
import sys
import tempfile
import contextlib
import io

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(HERE, "scripts")
FIXTURES = os.environ.get("MAGNET_FIXTURES", os.path.join(HERE, "tests", "fixtures"))
POSCAR_FE3SN = os.path.join(FIXTURES, "Fe3Sn_POSCAR")

if not os.path.exists(POSCAR_FE3SN):
    sys.exit("[skip] 缺少 fixture: %s (仓库应自带 tests/fixtures/, 或设 MAGNET_FIXTURES)"
             % POSCAR_FE3SN)

sys.path.insert(0, SCRIPTS)

for _stream in (sys.stdout, sys.stderr):
    try:
        _stream.reconfigure(errors="replace")
    except Exception:
        pass

import build_kagome_cases as B        # noqa: E402
import ingest_dft as IG               # noqa: E402
import solve_J_kagome as S            # noqa: E402

fails = []


def check(name, cond, extra=""):
    print(("PASS " if cond else "FAIL ") + name + (("   " + extra) if extra else ""))
    if not cond:
        fails.append(name)


ORDERS = list(B.SELECTED)          # ['fm', 'q0_ud', 'sq3_updntri', 'sq3_alt']
TMP = tempfile.mkdtemp(prefix="_ingest_test_")


def _mk_cases(tag):
    """建一套干净的 cases/ (真几何, 由构造器生成 —— 不是手搓的 fixture)。"""
    d = os.path.join(TMP, tag, "cases")
    old = sys.argv
    sys.argv = ["build_kagome_cases.py", POSCAR_FE3SN, "--element", "Fe", "--nonmag", "Sn",
                "--moment", "2.0", "--spin", "1.0", "--supercell", "3", "3", "1",
                "--outdir", d]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            B.main()
    finally:
        sys.argv = old
    if not os.path.exists(os.path.join(d, "magnetic_orders.json")):
        raise SystemExit("build_kagome_cases 没有产出 magnetic_orders.json: %s" % d)
    return d


CASES = _mk_cases("main")


def _poscar_lat(cases, order="fm"):
    """从 cases/<order>/POSCAR 取 (scale 后的) 晶格矢量, 供合成 OUTCAR 对齐。"""
    p = os.path.join(cases, order, "POSCAR")
    ln = [x for x in open(p, encoding="utf-8").read().splitlines() if x.strip()]
    sc = float(ln[1].split()[0])
    return [[float(v) * sc for v in ln[2 + i].split()[:3]] for i in range(3)]


LAT = _poscar_lat(CASES)
COUNTS = "27 9"          # Fe3Sn 3x3x1: 27 Fe + 9 Sn
N_IONS = 36
J_TRUE = [1.5, -0.8, 0.4]

_meta = json.load(open(os.path.join(CASES, "magnetic_orders.json"), encoding="utf-8"))
_COEFF = _meta["coefficients_per_cell"]
_S2 = float(_meta["S2"])


def _energies():
    return {nm: sum(_COEFF[nm][m] * (J_TRUE[m] * _S2) for m in range(3)) / 1000.0
            for nm in ORDERS}


E = _energies()

# ===========================================================================
# (1) 闸门策略表本身: 硬/软的分工不能被悄悄改
# ===========================================================================
names = [n for n, _h, _d in IG.CHECKS]
check("ingest: CHECKS 表里有 11 条闸门 (5 硬 + 2 软 + 4 硬)",
      len(IG.CHECKS) == 11, "got %d: %s" % (len(IG.CHECKS), names))
check("ingest: 闸门名不重复", len(set(names)) == len(names), "%s" % names)

HARD_WANT = {"files_present", "size", "timing", "scf", "ionic_steps", "convention",
             "structure", "not_shadowed", "spread"}
SOFT_WANT = {"accuracy", "moment"}
hard_got = {n for n, h, _d in IG.CHECKS if h}
soft_got = set(names) - hard_got
check("ingest: 9 条硬闸门 (收敛/口径/结构/防旧文件/离散度都在里面)",
      hard_got == HARD_WANT, "hard=%s" % sorted(hard_got))
check("ingest: `accuracy` 与 `moment` 是软的 (静态 NSW=1 本就不打印 accuracy)",
      soft_got == SOFT_WANT, "soft=%s" % sorted(soft_got))
check("ingest: check_policy 对已知闸门回真值",
      IG.check_policy("scf") is True and IG.check_policy("accuracy") is False)
try:
    IG.check_policy("no_such_check")
    _kerr = False
except KeyError:
    _kerr = True
check("ingest: check_policy 对未知闸门 raise KeyError (不静默当软闸门)", _kerr)


# ===========================================================================
# (2) parse_outcar: 合成 OUTCAR 的解析口径
# ===========================================================================
_p = os.path.join(TMP, "one", "OUTCAR")
os.makedirs(os.path.dirname(_p), exist_ok=True)
IG.synth_outcar(_p, -100.12345678, steps=1, n_ions=N_IONS, counts=COUNTS, lattice=LAT,
                lsorbit=False)
d = IG.parse_outcar(_p, min_kb=0)
check("parse: free energy TOTEN 读准 (8 位小数)",
      abs(d["energy_eV"] - (-100.12345678)) < 1e-9, "got %r" % (d["energy_eV"],))
check("parse: n_ionic_steps == 1", d["n_ionic_steps"] == 1, "got %r" % (d["n_ionic_steps"],))
check("parse: SCF 收敛标记被识别", d["scf_converged"] is True)
check("parse: timing 段被识别", d["has_timing"] is True)
check("parse: ions per type 合计 = 36", d["n_ions"] == N_IONS, "got %r" % (d["n_ions"],))
check("parse: 元素块计数 = [27, 9]", d["ion_counts"] == [27, 9], "got %r" % (d["ion_counts"],))
check("parse: 非共线/无 SOC 被识别 (LNONCOLLINEAR=T, LSORBIT=F)",
      d["noncollinear"] is True and d["lsorbit"] is False,
      "nc=%r soc=%r" % (d["noncollinear"], d["lsorbit"]))
# SOC 开关也要能区分 (P4.2 的 SOC 抽查靠它)
_p_soc = os.path.join(TMP, "one_soc", "OUTCAR")
os.makedirs(os.path.dirname(_p_soc), exist_ok=True)
IG.synth_outcar(_p_soc, -100.0, steps=1, n_ions=N_IONS, counts=COUNTS, lattice=LAT,
                lsorbit=True)
d_soc = IG.parse_outcar(_p_soc, min_kb=0)
check("parse: LSORBIT=T 被识别", d_soc["lsorbit"] is True, "got %r" % (d_soc["lsorbit"],))
check("parse: 读出 direct lattice vectors (结构闸门 8 需要)",
      d["lattice"] is not None
      and max(abs(d["lattice"][i][j] - LAT[i][j]) for i in range(3) for j in range(3)) < 1e-6,
      "got %r" % (d["lattice"],))

# 两个离子步 ⇒ 松弛过, 不是静态
_p2 = os.path.join(TMP, "two", "OUTCAR")
os.makedirs(os.path.dirname(_p2), exist_ok=True)
IG.synth_outcar(_p2, -100.0, steps=2, n_ions=N_IONS, counts=COUNTS, lattice=LAT)
d2 = IG.parse_outcar(_p2, min_kb=0)
check("parse: 两个 TOTEN ⇒ n_ionic_steps == 2 (松弛过能被看出来)",
      d2["n_ionic_steps"] == 2, "got %r" % (d2["n_ionic_steps"],))


# ===========================================================================
# (3) 好数据: 四态齐全 ⇒ 收回成功, manifest 结构与下游契约一致
# ===========================================================================
def _write_runs(tag, mutate=None, orders=None):
    """写一套 runs/<order>/OUTCAR; mutate(order, path, energy) 可逐序改写。"""
    rd = os.path.join(TMP, tag, "runs")
    for nm in (orders or ORDERS):
        os.makedirs(os.path.join(rd, nm), exist_ok=True)
        p = os.path.join(rd, nm, "OUTCAR")
        e = E.get(nm, -1.0)
        kw = {}
        if mutate:
            mutate(nm, p, e, kw)
        IG.synth_outcar(p, e, steps=kw.get("steps", 1), n_ions=kw.get("n_ions", N_IONS),
                        counts=kw.get("counts", COUNTS), lattice=kw.get("lattice", LAT),
                        scf=kw.get("scf", True), timing=kw.get("timing", True))
    return rd


_runs = _write_runs("good")
_cases = os.path.join(TMP, "good", "cases")
shutil.copytree(CASES, _cases)
ok, prob, warn, st, mpath = IG.ingest(_runs, _cases, orders=ORDERS, min_kb=0, quiet=True)
check("good: 四态齐全 ⇒ 收回成功", ok and not prob, "ok=%r problems=%s" % (ok, prob))
check("good: 四个 OUTCAR 都落到 cases/<order>/", 
      all(os.path.exists(os.path.join(_cases, o, "OUTCAR")) for o in ORDERS))
check("good: 每个 OUTCAR 带 .sha256 旁证",
      all(os.path.exists(os.path.join(_cases, o, "OUTCAR.sha256")) for o in ORDERS))
man = json.load(open(mpath, encoding="utf-8"))
check("good: manifest schema = dft_ingest/v1", man["schema"] == "dft_ingest/v1",
      "got %r" % (man["schema"],))
check("good: manifest 记了四个序的能量, 且与合成值一致",
      all(abs(man["orders"][o]["energy_eV"] - E[o]) < 1e-9 for o in ORDERS))
check("good: manifest 的每序 sha256 是 64 位十六进制",
      all(len(f.get("sha256") or "") == 64 for f in man["files"]),
      "%s" % [len(f.get("sha256") or "") for f in man["files"]])
check("good: manifest 声明「不证明物理正确性」(免责条款不许丢)",
      "不证明" in (man.get("note") or ""), "%r" % (man.get("note"),))
check("good: 写出 ingest_report.md", os.path.exists(os.path.join(_cases, "ingest_report.md")))

# 下游契约: solve_J_kagome 的凭证闸门必须能读懂这份 manifest
prov, pprob = S.read_ingest_manifest(_cases)
check("下游契约: solve_J_kagome.read_ingest_manifest 读得懂 (prov 非 None)",
      prov is not None and pprob is None, "problem=%r" % (pprob,))
check("下游契约: 凭证里四个序齐 + 记了口径",
      prov and set(prov["orders"]) == set(ORDERS) and prov["energy_convention"],
      "%s" % (sorted(prov["orders"]) if prov else None,))
check("下游契约: 凭证里带逐序 sha256",
      prov and all(len(v or "") == 64 for v in prov["outcar_sha256"].values()),
      "%s" % (prov.get("outcar_sha256") if prov else None,))
# 缺一个序的 manifest ⇒ 下游必须拒 (而不是拿三个序凑)
_cases_short = os.path.join(TMP, "shortprov", "cases")
shutil.copytree(_cases, _cases_short)
_man_short = json.load(open(os.path.join(_cases_short, "ingest_manifest.json"),
                            encoding="utf-8"))
_man_short["orders"].pop("sq3_alt")
json.dump(_man_short, open(os.path.join(_cases_short, "ingest_manifest.json"), "w",
                           encoding="utf-8"), ensure_ascii=False)
prov2, pprob2 = S.read_ingest_manifest(_cases_short)
check("下游契约: manifest 缺一个序 ⇒ 明确报出缺哪个", prov2 is None and "sq3_alt" in (pprob2 or ""),
      "problem=%r" % (pprob2,))


# ===========================================================================
# (4) 逐条硬闸门: 各自的坏数据都必须被拒, 且**一个盘都不许落**
# ===========================================================================
def _reject_case(tag, mutate=None, orders=None, **kw):
    """跑一次 ingest, 返回 (ok, problems, 目标 cases 里 OUTCAR 的个数)。"""
    rd = _write_runs(tag, mutate=mutate, orders=orders)
    cd = os.path.join(TMP, tag, "cases")
    shutil.copytree(CASES, cd)
    ok, prob, _warn, _st, _mp = IG.ingest(rd, cd, orders=ORDERS, min_kb=kw.pop("min_kb", 0),
                                          quiet=True, **kw)
    n_out = sum(1 for o in ORDERS if os.path.exists(os.path.join(cd, o, "OUTCAR")))
    return ok, prob, n_out, cd


def _mut_scf(order, path, e, kw):
    if order == "q0_ud":
        kw["scf"] = False


_ok, _prob, _n, _ = _reject_case("bad_scf", mutate=_mut_scf)
check("闸门 scf: 有一个序 SCF 未收敛 ⇒ 拒绝", not _ok and _prob)
check("闸门 scf: 拒绝时 problems 指向 scf 闸门",
      any(c == "scf" for _o, c, _m in _prob), "%s" % [(o, c) for o, c, _m in _prob])
check("闸门 scf: **原子性** —— 拒绝后 cases/ 里一个 OUTCAR 都没有",
      _n == 0, "落盘 %d 个" % _n)


def _mut_timing(order, path, e, kw):
    if order == "sq3_alt":
        kw["timing"] = False


_ok, _prob, _n, _ = _reject_case("bad_timing", mutate=_mut_timing)
check("闸门 timing: 缺 timing 段 ⇒ 拒绝", not _ok and any(c == "timing" for _o, c, _m in _prob))
check("闸门 timing: 原子性 (一个都没落)", _n == 0, "落盘 %d 个" % _n)


def _mut_steps(order, path, e, kw):
    if order == "fm":
        kw["steps"] = 2


_ok, _prob, _n, _ = _reject_case("bad_steps", mutate=_mut_steps)
check("闸门 ionic_steps: 两离子步 (松弛过) ⇒ 拒绝",
      not _ok and any(c == "ionic_steps" for _o, c, _m in _prob),
      "%s" % [(o, c) for o, c, _m in _prob])


def _mut_struct(order, path, e, kw):
    if order == "sq3_updntri":
        kw["n_ions"] = 35
        kw["counts"] = "26 9"


_ok, _prob, _n, _ = _reject_case("bad_struct", mutate=_mut_struct)
check("闸门 structure: 原子数/元素块不一致 ⇒ 拒绝",
      not _ok and any(c == "structure" for _o, c, _m in _prob),
      "%s" % [(o, c) for o, c, _m in _prob])


def _mut_lattice(order, path, e, kw):
    if order == "q0_ud":
        bad = [list(v) for v in LAT]
        bad[0][0] += 0.05                    # 0.05 Å > 1e-4 容差
        kw["lattice"] = bad


_ok, _prob, _n, _ = _reject_case("bad_lattice", mutate=_mut_lattice)
check("闸门 structure: 晶格矢量差 0.05 Å (> 1e-4) ⇒ 拒绝",
      not _ok and any(c == "structure" for _o, c, _m in _prob),
      "%s" % [(o, c) for o, c, _m in _prob])


def _mut_nogap(order, path, e, kw):
    pass


_rd_nogap = os.path.join(TMP, "bad_nogap", "runs")
for _nm in ORDERS:
    os.makedirs(os.path.join(_rd_nogap, _nm), exist_ok=True)
    IG.synth_outcar(os.path.join(_rd_nogap, _nm, "OUTCAR"), -100.0, steps=1, n_ions=N_IONS,
                    counts=COUNTS, lattice=LAT)
_cd_nogap = os.path.join(TMP, "bad_nogap", "cases")
shutil.copytree(CASES, _cd_nogap)
_ok, _prob, _w, _st, _mp = IG.ingest(_rd_nogap, _cd_nogap, orders=ORDERS, min_kb=0, quiet=True)
_n = sum(1 for _o in ORDERS if os.path.exists(os.path.join(_cd_nogap, _o, "OUTCAR")))
check("闸门 spread/能量不退化: 四态能量全同 ⇒ 拒绝 (磁序没生效)",
      not _ok and any(c == "spread" for _o, c, _m in _prob),
      "%s" % [(o, c) for o, c, _m in _prob])
check("闸门 spread: 原子性 (一个都没落)", _n == 0, "落盘 %d 个" % _n)


_rd_size = os.path.join(TMP, "bad_size", "runs")
for _nm in ORDERS:
    os.makedirs(os.path.join(_rd_size, _nm), exist_ok=True)
    IG.synth_outcar(os.path.join(_rd_size, _nm, "OUTCAR"), E[_nm], steps=1, n_ions=N_IONS,
                    counts=COUNTS, lattice=LAT)
_cd_size = os.path.join(TMP, "bad_size", "cases")
shutil.copytree(CASES, _cd_size)
_ok, _prob, _w, _st, _mp = IG.ingest(_rd_size, _cd_size, orders=ORDERS, min_kb=200, quiet=True)
_n = sum(1 for _o in ORDERS if os.path.exists(os.path.join(_cd_size, _o, "OUTCAR")))
check("闸门 size: --min-kb 200 但文件只有几 KB ⇒ 拒绝 (空壳)",
      not _ok and any(c == "size" for _o, c, _m in _prob),
      "%s" % [(o, c) for o, c, _m in _prob])


# 缺一个序: 这是**结构级**硬错 —— ingest 直接 SystemExit (在跑闸门之前), 不返回元组
_rd_miss = _write_runs("missing", orders=[o for o in ORDERS if o != "sq3_alt"])
_cd_miss = os.path.join(TMP, "missing", "cases")
shutil.copytree(CASES, _cd_miss)
_sexit = None
try:
    with contextlib.redirect_stdout(io.StringIO()):
        IG.ingest(_rd_miss, _cd_miss, orders=ORDERS, min_kb=0, quiet=True)
except SystemExit as _e:
    _sexit = str(_e)
check("闸门 files_present: 缺一个序 ⇒ SystemExit (结构级硬错, 不是可降级闸门)",
      _sexit is not None, "SystemExit=%r" % (_sexit,))
check("闸门 files_present: 报错点名缺的那个序", "sq3_alt" in (_sexit or ""),
      "%r" % (_sexit,))
check("闸门 files_present: 缺序时一个盘都不落 (连 gate 都不进)",
      not any(os.path.exists(os.path.join(_cd_miss, o, "OUTCAR")) for o in ORDERS))


# ===========================================================================
# (5) --lenient 真的生效 (曾经有一个 bug: policy_override 被全量预填, --lenient 形同虚设)
# ===========================================================================
_ok, _prob, _n, cd = _reject_case("lenient", mutate=_mut_steps)
check("lenient 基线: 两离子步默认被拒", not _ok and _n == 0)

_rd = _write_runs("lenient2", mutate=_mut_steps)
_cd = os.path.join(TMP, "lenient2", "cases")
shutil.copytree(CASES, _cd)
ok, prob, warn, st, mp = IG.ingest(_rd, _cd, orders=ORDERS, min_kb=0, quiet=True,
                                   lenient=["ionic_steps"])
check("lenient: --lenient ionic_steps ⇒ 不再拒 (硬降软真的生效)", ok and not prob,
      "ok=%r problems=%s" % (ok, prob))
check("lenient: 降级的那条进了 warnings (留痕, 不是静默放过)",
      any(c == "ionic_steps" for _o, c, _m in warn), "%s" % [(o, c) for o, c, _m in warn])
man_l = json.load(open(mp, encoding="utf-8"))
check("lenient: manifest 的 checks 里该闸门记为 hard=False (下游能查到)",
      any(c["name"] == "ionic_steps" and not c["hard"] for c in man_l["checks"]),
      "%s" % [(c["name"], c["hard"]) for c in man_l["checks"]])
prov_l, _p = S.read_ingest_manifest(_cd)
check("lenient: 凭证里记下降级闸门名 (J_result.json 的 notes 会带上它)",
      prov_l and "ionic_steps" in prov_l["lenient"],
      "%s" % (prov_l.get("lenient") if prov_l else None,))
check("lenient: 凭证的 lenient 只列「非硬」闸门, 且含 accuracy/moment (它们天生是软的)",
      prov_l and set(prov_l["lenient"]) == {"ionic_steps", "accuracy", "moment"},
      "%s" % (prov_l.get("lenient") if prov_l else None,))
# 反例: 不降级时 lenient 只剩那两条天生软的闸门
prov_h, _p = S.read_ingest_manifest(_cases)
check("lenient: 没降级时凭证里只剩 accuracy/moment (硬闸门都在岗)",
      prov_h and set(prov_h["lenient"]) == {"accuracy", "moment"},
      "%s" % (prov_h.get("lenient") if prov_h else None,))


# ===========================================================================
# (6) 防旧文件充数 (not_shadowed) + --force + --dry-run
# ===========================================================================
ok1, _p1, _w1, _s1, _m1 = IG.ingest(_runs, _cases, orders=ORDERS, min_kb=0, quiet=True,
                                     force=True)
_rd_same = _write_runs("shadow", orders=ORDERS)
_cd_sh = os.path.join(TMP, "shadow", "cases")
shutil.copytree(_cases, _cd_sh)          # 先有 good 的 OUTCAR
ok_s, prob_s, _w, _s, _m = IG.ingest(_rd_same, _cd_sh, orders=ORDERS, min_kb=0, quiet=True)
check("not_shadowed: 陈旧 OUTCAR (sha256 与 cases 里已有的一样) ⇒ 拒绝",
      not ok_s and any(c == "not_shadowed" for _o, c, _m in prob_s),
      "%s" % [(o, c) for o, c, _m in prob_s])
ok_sf, prob_sf, _w, _s, _m = IG.ingest(_rd_same, _cd_sh, orders=ORDERS, min_kb=0, quiet=True,
                                       force=True)
check("not_shadowed: --force 放行", ok_sf and not prob_sf, "problems=%s" % (prob_sf,))

_rd_d = _write_runs("dry", orders=ORDERS)
_cd_d = os.path.join(TMP, "dry", "cases")
shutil.copytree(CASES, _cd_d)
ok_d, prob_d, _w, _s, _m = IG.ingest(_rd_d, _cd_d, orders=ORDERS, min_kb=0, quiet=True,
                                     dry_run=True)
check("dry-run: 好数据下 dry_run 返回 ok", ok_d and not prob_d, "problems=%s" % (prob_d,))
check("dry-run: **不落盘** (cases/<order>/OUTCAR 不存在)",
      not any(os.path.exists(os.path.join(_cd_d, o, "OUTCAR")) for o in ORDERS))
check("dry-run: 不写 manifest", not os.path.exists(os.path.join(_cd_d, "ingest_manifest.json")))
_rd_bad = _write_runs("drybad", mutate=_mut_scf)
_cd_bad = os.path.join(TMP, "drybad", "cases")
shutil.copytree(CASES, _cd_bad)
ok_db, prob_db, _w, _s, _m = IG.ingest(_rd_bad, _cd_bad, orders=ORDERS, min_kb=0, quiet=True,
                                      dry_run=True)
check("dry-run: 坏数据下 dry_run 也报出来 (先看闸门再决定收不收)",
      not ok_db and any(c == "scf" for _o, c, _m in prob_db),
      "%s" % [(o, c) for o, c, _m in prob_db])
check("dry-run: 坏数据下也**不落盘**",
      not any(os.path.exists(os.path.join(_cd_bad, o, "OUTCAR")) for o in ORDERS))

# ===========================================================================
print()
if fails:
    print("%d 个失败: %s" % (len(fails), fails))
    sys.exit(1)
print("全部通过")
