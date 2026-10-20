#!/usr/bin/env python
"""mag_vectors.py —— 非共线 / 共线自旋矢量的**唯一生成器**（P3.15）。

把原先散落在三处的角度规则收拢成「一张表 + 一组纯函数」：

  * `build_magnetic_cases_tri.py` 的 `spins_120()` / `spins_collinear()`
  * `magnetic_orders.py` 的 `planar120()`
  * `assemble_dft_tri.py` 的 3 分量 MAGMOM 模板

规格层（`magnetic_orders.py`）、builder（`build_*_cases.py`）与 DFT 组装层
（`assemble_dft.py` / `assemble_dft_tri.py`）从此共用同一实现：加一类序只需要
在下面的表里加一行，不再是复制一份角度代码。

## 角度规则表

| 序 | 自 z 轴倾角 | 面内方位 | 返回 | 需要 SOC |
|---|---|---|---|---|
| `collinear` | 0 / 180° | — | `(0, 0, sigma)` | 否 |
| `planar120` | 90° | `((n0+n1) % period) * 120°` | `(cos, sin, 0)` | 否 |
| `cant` | `±theta` | `azimuth` | `(sin t cos p, sin t sin p, ±cos t)` | **是** |

上表在代码里就是 `ANGLE_TABLE` + `spin_*` 三个函数。

## 约定（与主链其它部分一致，别改）

1. 返回的矢量一律是**单位矢量**（`|S| = 1`）；乘 `moment`(μB) 才是 MAGMOM 的物理值。
2. 矢量按 **POSCAR 原子顺序**给出（非磁原子为 0）；只算磁性原子时用 n_mag 视图。
3. **弱铁磁倾斜（`cant`）的两个子晶格「面内分量同向、z 分量反平行」** —— 这才是净
   铁磁矩的来源（每对净矩 = `2 sin(theta)`），也是 DMI / 弱倾斜的物理图像。
4. 需要 SOC 的序（`needs_soc`）在组装 DFT 输入时**必须**开 `LSORBIT`，否则直接报错，
   不允许按共线近似继续（NEW_TYPES_PLAN §115「不静默降级」）。

## 用法

  python mag_vectors.py --selftest
"""
import argparse
import math
import os

import numpy as np

# ---------------------------------------------------------------------------
# 需要 SOC (LSORBIT) 才能算的序：弱铁磁倾斜 / 自旋螺旋 / 任何非共线小量
# ---------------------------------------------------------------------------
SOC_ORDERS = ("cant", "cant_a", "cant_b", "cant_diag", "spiral", "tilt")

# 子晶格相位/角度表（纯数据，供文档与 selftest 交叉校验）
ANGLE_TABLE = {
    "collinear": "S = (0, 0, sigma)                  sigma = ±1",
    "planar120": "color = (n0+n1) % 3                S = (cos(120°·color), sin(120°·color), 0)",
    "cant": "S = (sin t·cos p, sin t·sin p, sigma·cos t)   两子晶格面内同向、z 反平行",
}


def needs_soc(order):
    """该磁序是否必须有 SOC（`LSORBIT`）才能算。"""
    name = str(order).strip().lower()
    return name in SOC_ORDERS or name.startswith(("cant", "tilt", "spiral"))


# ---------------------------------------------------------------------------
# 自旋矢量生成
# ---------------------------------------------------------------------------
def unit(v):
    """归一化为单位矢量。"""
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    return v / n if n > 0 else v


def spin_collinear(signs):
    """共线序：`(0, 0, sigma)`，signs 为 ±1 数组。返回 (n, 3)。"""
    signs = np.asarray(signs, dtype=float)
    s = np.zeros((signs.size, 3))
    s[:, 2] = signs
    return s


def spin_planar(theta_deg):
    """面内自旋：`(cos θ, sin θ, 0)`。θ 可以是标量或数组。返回 (n, 3)。"""
    th = np.radians(np.atleast_1d(np.asarray(theta_deg, dtype=float)))
    s = np.zeros((th.size, 3))
    s[:, 0] = np.cos(th)
    s[:, 1] = np.sin(th)
    return s


def spin_120deg(colors, period=3):
    """三角/蜂窝的 120° 共面非共线序：θ = （颜色 · 360°/period）。

    颜色通常取 `(n0 + n1) % period`。**不需要 SOC**（受挫晶格的几何阻挫基态）。
    """
    colors = np.asarray(colors, dtype=float)
    return spin_planar(colors * (360.0 / float(period)))


def spin_cant(signs, theta_deg, azimuth_deg=0.0):
    """弱铁磁倾斜（canted AFM）：自 z 轴倾斜 `theta_deg`，面内方位 `azimuth_deg`。

    约定（见模块 docstring 第 3 条）：**两个子晶格的面内分量同向、z 分量反平行** ——
    `sigma = +1` 给 `(sin t cos p, sin t sin p, +cos t)`，`sigma = -1` 给
    `(sin t cos p, sin t sin p, -cos t)`，于是每对净矩 = `2 sin t`（沿面内方位）。
    `theta = 0` 退化为纯共线 `(0, 0, ±1)`；`theta = 90°` 退化为面内 90° 序。
    **需要 SOC**（DMI / 磁各向异性驱动），`needs_soc("cant")` 为 True。
    """
    signs = np.asarray(signs, dtype=float)
    t = math.radians(float(theta_deg))
    p = math.radians(float(azimuth_deg))
    s = np.zeros((signs.size, 3))
    s[:, 0] = math.sin(t) * math.cos(p)
    s[:, 1] = math.sin(t) * math.sin(p)
    s[:, 2] = signs * math.cos(t)
    return s


def spin_all(mag_spins, mag_idx, n_atoms):
    """把 n_mag 个单位自旋按 `mag_idx` 放进 POSCAR 原子顺序的 (n_atoms, 3) 数组。"""
    mag_spins = np.asarray(mag_spins, dtype=float).reshape(len(mag_idx), 3)
    out = np.zeros((int(n_atoms), 3))
    for k, i in enumerate(mag_idx):
        out[int(i)] = mag_spins[k]
    return out


# ---------------------------------------------------------------------------
# MAGMOM 文本：生成 / 解析 / 校验
# ---------------------------------------------------------------------------
def magmom_values(mag_spins, moment, mag_idx=None, n_atoms=None):
    """单位自旋 -> MAGMOM 物理值（乘 moment）。返回 (n_atoms, 3) 或 (n_mag, 3)。"""
    v = np.asarray(mag_spins, dtype=float) * float(moment)
    if mag_idx is None:
        return v
    if n_atoms is None:
        n_atoms = int(max(mag_idx)) + 1
    return spin_all(np.asarray(mag_spins, dtype=float), mag_idx, n_atoms) * float(moment)


def render_magmom(values, prefix="MAGMOM   = "):
    """(n, 3) 值 -> VASP 的 `MAGMOM = ...` 行（每原子 3 分量，POSCAR 顺序）。

    值一个都不省略（不写 `N*val` 缩写）—— 非共线 MAGMOM 用缩写极易被 VASP 按
    「值的个数」错解，显式展开最安全。
    """
    vals = np.asarray(values, dtype=float).reshape(-1, 3)
    body = " ".join(f"{x:.6f} {y:.6f} {z:.6f}" for x, y, z in vals)
    return f"{prefix}{body}"


def expand_repeats(tokens):
    """展开 VASP 的 `N*val` 缩写（也接受 `N*vx vy vz` 这种非标准写法）。"""
    out = []
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        if "*" in tok:
            n_s, _, v_s = tok.partition("*")
            try:
                n = int(n_s)
                v = float(v_s)
            except ValueError:
                raise ValueError(f"无法解析 MAGMOM 缩写 {tok!r}")
            out.extend([v] * n)
        else:
            out.append(float(tok))
        i += 1
    return out


def parse_magmom(text, n_atoms=None):
    """解析 `MAGMOM = ...` 行 -> (values_flat, n_per_atom)。

    也接受**整份 INCAR 文本**（自动挑出含 `MAGMOM` 的那一行）—— 所以
    `parse_magmom(open("runs/fm/INCAR").read(), n_atoms=N)` 可以直接用。

    `n_per_atom` 是 1（共线）或 3（非共线）；无法整除时返回 (values, 0)。
    值个数能同时被 3 和 1 整除时（例如 3 个值：1 个非共线原子 vs 3 个共线原子）
    **必须**给 `n_atoms` 才能定解；不给则偏向 3 分量（本模块的语境是非共线）。
    """
    lines = [ln.strip() for ln in str(text).splitlines() if ln.strip()]
    line = next((ln for ln in lines if ln.upper().startswith("MAGMOM")), None)
    if line is None:
        line = " ".join(lines) if lines else str(text)
    body = line.split("=", 1)[1] if "=" in line else line
    body = body.split("!", 1)[0]
    vals = expand_repeats(body.split())
    n = len(vals)
    if n_atoms:
        for k in (3, 1):
            if n == int(n_atoms) * k:
                return np.array(vals), k
        return np.array(vals), 0
    for k in (3, 1):
        if n % k == 0:
            return np.array(vals), k
    return np.array(vals), 0


def check_magmom(text, n_atoms, moment=None, tol=1e-6, want_ncomp=None, n_mag=None):
    """校验 MAGMOM 行。返回 (ok, n_per_atom, problems[])。

    查：值个数 == n_atoms * ncomp；分量数符合 `want_ncomp`；有磁矩的原子个数 == `n_mag`
    （给了 `n_mag` 时）；每个有磁矩原子的 |m| == `moment`（给了 `moment` 时）。
    共线行用在非共线上下文里，会被 `want_ncomp=3` 直接逮住。
    """
    problems = []
    vals, ncomp = parse_magmom(text, n_atoms=n_atoms)
    if ncomp == 0:
        problems.append(f"MAGMOM 值个数 {vals.size} 既不是 n_atoms({n_atoms}) 的 1 倍也不是 3 倍")
        return False, 0, problems
    if want_ncomp and ncomp != want_ncomp:
        problems.append(f"MAGMOM 是 {ncomp} 分量, 期望 {want_ncomp} 分量")
    m = vals.reshape(n_atoms, ncomp)
    mag = np.linalg.norm(m, axis=1) if ncomp == 3 else np.abs(m[:, 0])
    nz = mag > tol
    if n_mag is not None and int(nz.sum()) != int(n_mag):
        problems.append(f"有磁矩的原子数 {int(nz.sum())} != n_mag {int(n_mag)}")
    if moment is not None:
        if not nz.any():
            problems.append(f"全部原子磁矩为 0, moment={moment} 没落实")
        else:
            bad = mag[nz][np.abs(mag[nz] - abs(float(moment))) > 1e-4]
            if bad.size:
                problems.append(f"{bad.size} 个原子的 |m| 不等于 moment={moment}: 例 {bad[:3].tolist()}")
    return not problems, ncomp, problems


# ---------------------------------------------------------------------------
# selftest（无 DFT）
# ---------------------------------------------------------------------------
def _selftest():
    checks = []

    def ck(name, cond, extra=""):
        checks.append(bool(cond))
        print(f"  [{'ok ' if cond else 'FAIL'}] {name}" + (f"  {extra}" if extra else ""))

    s = spin_collinear([1, -1, 1])
    ck("collinear: 单位矢量 (0,0,±1)", np.allclose(np.abs(s[:, 2]), 1.0)
       and np.allclose(s[:, :2], 0.0))

    c3 = spin_120deg(np.array([0, 1, 2]))
    ang = np.degrees(np.arctan2(c3[:, 1], c3[:, 0])) % 360.0
    ck("planar120: 角度恰为 0/120/240 度", np.allclose(np.sort(ang), [0.0, 120.0, 240.0], atol=1e-9),
       f"got {np.round(np.sort(ang), 6).tolist()}")
    ck("planar120: 三子晶格闭合 sum S = 0 (受挫三角)",
       np.allclose(c3.sum(axis=0), 0.0, atol=1e-12))
    ck("planar120: |S| = 1 且 Sz = 0", np.allclose(np.linalg.norm(c3, axis=1), 1.0)
       and np.allclose(c3[:, 2], 0.0))

    c = spin_cant(np.array([1, -1]), 8.0, azimuth_deg=30.0)
    t = math.radians(8.0)
    ck("cant: |S| = 1", np.allclose(np.linalg.norm(c, axis=1), 1.0))
    ck("cant: z 分量 = ±cos(theta)", np.allclose(c[:, 2], [math.cos(t), -math.cos(t)]))
    ck("cant: 两子晶格**面内同向** (净铁磁矩的来源)",
       np.allclose(c[0, :2], c[1, :2]) and abs(np.linalg.norm(c[0, :2]) - math.sin(t)) < 1e-12,
       f"面内模 {np.linalg.norm(c[0, :2]):.6f} = sin(theta)")
    ck("cant: 每对净矩 = 2 sin(theta)",
       abs(np.linalg.norm(c[0] + c[1]) - 2 * math.sin(t)) < 1e-12)
    ck("cant: theta=0 退化为纯共线", np.allclose(spin_cant([1, -1], 0.0), spin_collinear([1, -1])))

    vals = magmom_values(spin_120deg(np.array([0, 1, 2])), 3.0)
    txt = render_magmom(vals)
    back, ncomp = parse_magmom(txt)
    ck("MAGMOM 文本往返 (render -> parse) 逐位一致",
       ncomp == 3 and np.allclose(back.reshape(-1, 3), vals, atol=1e-6))
    ok, nc, probs = check_magmom(txt, 3, moment=3.0, want_ncomp=3, n_mag=3)
    ck("check_magmom: 合法非共线行通过", ok and nc == 3, str(probs))
    ok1, nc1, probs1 = check_magmom("MAGMOM = 3*5.0", 3, moment=5.0, want_ncomp=3, n_mag=3)
    ck("check_magmom: 逮住『共线行 (1 分量) 用于非共线上下文』",
       (not ok1) and nc1 == 1, str(probs1)[:60])
    ok2, _, probs2 = check_magmom(render_magmom(np.zeros((3, 3))), 3, moment=3.0, n_mag=3)
    ck("check_magmom: 逮住全 0 行 (moment 没落实)", not ok2, str(probs2)[:60])
    v_bad = np.array([[3.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 3.0]])
    ok3, _, probs3 = check_magmom(render_magmom(v_bad), 3, moment=3.0, n_mag=3)
    ck("check_magmom: 逮住单个原子 |m| 偏离 moment", not ok3, str(probs3)[:60])

    idx = [2, 4]
    allv = spin_all(spin_collinear([1, -1]), idx, 6)
    ck("spin_all: 按 mag_idx 落位, 其余为 0",
       np.allclose(allv[2], [0, 0, 1]) and np.allclose(allv[4], [0, 0, -1])
       and np.allclose(np.delete(allv, idx, axis=0), 0.0))

    ck("needs_soc: cant/spiral 为真, fm/120deg/stripy 为假",
       needs_soc("cant_a") and needs_soc("spiral") and needs_soc("tilt_b")
       and not needs_soc("fm") and not needs_soc("120deg") and not needs_soc("stripy_ab"))

    # 与规格层 cross-check: magnetic_orders.planar120 必须与本模块一致 (同一物理, 两处实现)
    try:
        import magnetic_orders as MO
        cell = np.array([[0, 0], [1, 0], [0, 1], [1, 1], [2, 0], [2, 1], [3, 3]], dtype=float)
        ck("与 magnetic_orders.planar120 逐位一致 (同一角规则的单一来源)",
           np.allclose(MO.planar120(cell, 3), spin_120deg((cell[:, 0] + cell[:, 1]) % 3)))
    except Exception as e:  # pragma: no cover
        ck(f"与 magnetic_orders.planar120 交叉校验 (导入失败: {e})", False)

    # 与构造器 cross-check: 现场建一套三角构型 (GdI2 夹具), 核对共线标量 vs 非共线 3 分量。
    # 不读外部项目目录: 构造器是唯一来源, 样本也由它现场生成 (P6: 默认全共线)。
    import json
    import shutil
    import subprocess
    import sys
    import tempfile

    here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.dirname(here)
    poscar = os.path.join(repo, "tests", "fixtures", "GdI2_POSCAR")
    tri = os.path.join(here, "build_magnetic_cases_tri.py")
    tmp = tempfile.mkdtemp(prefix="_mv_tri_")

    def _build(tag, extra=()):
        out = os.path.join(tmp, tag)
        cmd = [sys.executable, tri, poscar, "--supercell", "5", "5", "1",
               "--element", "Gd", "--nonmag", "I", "--moment", "7.0", "--spin", "3.5",
               "--outdir", out] + list(extra)
        r = subprocess.run(cmd, capture_output=True, text=True)
        return out if r.returncode == 0 else None

    def _read(out):
        meta = json.load(open(os.path.join(out, "magnetic_orders.json"), encoding="utf-8"))
        rows = {}
        for name in meta["selected_orders"]:
            p = os.path.join(out, name, "magmom.txt")
            text = open(p, encoding="utf-8").read()
            toks = text.split("=", 1)[-1].split()
            rows[name] = (text, len(toks))
        return meta, rows

    try:
        out_c = _build("collinear")
        out_n = _build("noncollinear", ["--include-noncollinear"])
        if out_c is None or out_n is None:
            ck("构造器交叉校验: 现场建三角构型 (构造器失败)", False)
        else:
            meta_c, rows_c = _read(out_c)
            _meta_n, rows_n = _read(out_n)
            n_atoms = int(meta_c["n_atoms"])
            bad = [f"{k}: {n_tok} 个值" for k, (_t, n_tok) in rows_c.items()
                   if n_tok != n_atoms]
            ck(f"构造器交叉校验: 默认四态全共线标量 MAGMOM (n_atoms={n_atoms})",
               not bad, "; ".join(bad)[:80])
            text120, n_tok120 = rows_n.get("120deg", ("", 0))
            ck("构造器交叉校验: 只有 --include-noncollinear 才给 120deg, 且为 3 分量",
               n_tok120 == 3 * n_atoms, f"值个数={n_tok120} (应为 3*{n_atoms})")
            if text120 and n_tok120 == 3 * n_atoms:
                v = np.array([float(x) for x in text120.split("=", 1)[-1].split()]) * 7.0
                m = v.reshape(-1, 3)
                m = m[np.linalg.norm(m, axis=1) > 1e-6]
                a = np.degrees(np.arctan2(m[:, 1], m[:, 0]))
                cols = sorted({int(round(x / 120.0)) % 3 for x in a})
                ck("构造器交叉校验: 120deg 恰好只有 0/120/240 三色",
                   cols == [0, 1, 2], f"cols={cols}")
            else:
                ck("构造器交叉校验: 120deg 三色 (缺少 120deg 或分量不对, 跳过)", False)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    n_ok = sum(checks)
    print(f"\n--selftest: {'ALL PASS' if n_ok == len(checks) else 'FAIL'} ({n_ok}/{len(checks)} 项)")
    return 0 if n_ok == len(checks) else 1


def main():
    ap = argparse.ArgumentParser(description="非共线/共线自旋矢量生成器 (P3.15)")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return _selftest()
    ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
