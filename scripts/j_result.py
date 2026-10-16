#!/usr/bin/env python
"""j_result.py — P0.2: J 反解的统一产物 `J_result.json` 的 schema / 写出 / 读取 / 校验。

**为什么要它**: 现在 5 条 solve_*_J.py 各自 print 结果就完事, 下游 (poscar_to_mc3 / mc_2d)
只能靠人肉抄数。抄错一个是常数级的静默错误 —— 能量看着对, 物理错了。统一产物让
「J 及其来源、拟合质量」变成机器可读、可校验的一等公民。

**契约 (schema = "j_result/1")**:
    {
      "schema": "j_result/1", "created": "<ISO>",
      "material": "CrSBr" | null, "lattice": "rectangular",
      "source": "vasp" | "demo" | "literature",     # 数据来源, 不是拟合方法
      "confidence": "high" | "medium" | "low",      # 由 source/residual/cond 导出, 不可手填
      "spin": 1.5, "S2": 2.25, "n_mag": 16, "supercell": [4,4,1],
      "orders": ["fm", "afm_a", "afm_b", "afm_diag"],
      "A": [[1.0, -16.0, -16.0, -32.0], ...],       # 设计矩阵, 首列定 E0
      "shells": [
        {"name": "J1", "d": 3.5404561399999994,
         "sum_coeff": [-16.0, 16.0, -16.0, 16.0],   # 该壳在各磁序下的系数 (A 的第 k+1 列)
         "z_n": 2.0,                                # 每磁位点的键数 (配位数), 由 fm 行导出
         "pair": null}                              # 亚晶格过滤 (FGT 的 J3=AA/J4=BB 靠它区分)
      ],
      "E_eV": {"fm": -0.123, ...},                  # 每超胞总能 (eV)
      "E0_eV": 0.0123, "J_meV": {"J1": 2.79, ...}, "J_tilde_eV": {"J1": 6.28, ...},
      "cond": 2.34,                                 # cond(A); 病态闸门见 check_design
      "residual_eV": 1.2e-16,                       # |A·x - b|_inf; 拟合是否自洽
      "notes": ["..."]
    }

**"不静默降级" 规则** (这是本模块存在的理由, 不是形式主义):
    residual > RESIDUAL_FAIL  -> write() 直接 raise, 除非显式 force=True;
    residual > RESIDUAL_WARN  -> confidence 强制降到 "low" 并在 stdout 大声警告;
    source != "vasp"          -> confidence 封顶 "low" (demo/literature 不是本体系的结果)。

用法:
    from j_result import build, write, read, load_shells
    obj = build(lattice="rectangular", shells=[...], J={"J1": 2.79}, E={...}, E0=..., 
                cond=..., residual=..., spin=1.5, n_mag=16, orders=[...], A=A)
    write("cases/J_result.json", obj)
"""
import json
import math
import os
from datetime import datetime

SCHEMA = "j_result/1"

SOURCES = ("vasp", "demo", "literature")
CONFIDENCE = ("high", "medium", "low")

# 残差闸门 (eV / 超胞): 反解自洽性 |A·x - b|_inf
RESIDUAL_WARN = 1e-4      # 超过 -> 降级 "low" + 警告
RESIDUAL_FAIL = 1e-2      # 超过 -> 拒绝写出 (除非 force)
# 条件数闸门 (与 magnetic_orders.check_design 同源)
COND_WARN = 1e3


class JResultError(ValueError):
    """产物不合法 / 质量不达标。**故意用异常**: 宁可中断, 也不静默出 J。"""


# --------------------------------------------------------------------------- 构造
def confidence_of(source: str, residual: float, cond: float) -> str:
    """由 来源 + 残差 + 条件数 导出可信度。手填没有意义, 所以只在这里算。"""
    if source != "vasp":
        return "low"
    if not (isinstance(residual, (int, float)) and math.isfinite(residual)):
        return "low"
    if residual > RESIDUAL_WARN or (math.isfinite(cond) and cond > COND_WARN):
        return "low"
    return "high"


def build(*, lattice, shells, J, E, E0, cond, residual, spin, n_mag,
          orders=(), A=None, S2=None, supercell=None, material=None,
          source="vasp", z_from=None, notes=None, **extra):
    """组装 J_result 字典。shells: [{name, d, sum_coeff?, z_n?, pair?}] 或 (name, d) 元组。

    z_from: 用来反推配位数的**铁磁参照磁序名** (元胞内键数 = -c_k(fm))。默认取 orders[0]"。
    """
    if source not in SOURCES:
        raise JResultError(f"source 必须是 {SOURCES} 之一, 收到 {source!r}")
    orders = list(orders)
    S2 = float(spin) ** 2 if S2 is None else float(S2)

    sh = []
    for s in shells:
        if isinstance(s, dict):
            d = dict(s)
        elif isinstance(s, (tuple, list)):
            d = {"name": str(s[0]), "d": float(s[1])}
            if len(s) > 2:
                d["pair"] = s[2]
        else:
            raise JResultError(f"shell 无法解析: {s!r}")
        d.setdefault("pair", None)
        d.setdefault("tol", None)      # 距离容差; None = 由壳间距推出 (见下)
        sh.append(d)

    names = [s["name"] for s in sh]
    if len(set(names)) != len(names):
        raise JResultError(f"壳名重复: {names}")
    for nm in J:
        if nm not in names:
            raise JResultError(f"J 里有壳 {nm!r}, 但 shells 里没有 (shells={names})")

    # tol 自描述: 没显式给就按相邻壳间距的一半推 (壳之间不重叠), 上限 0.3。
    # 消费端 (poscar_to_mc3) 靠它把键分到最近的显式壳, 写进产物才不会靠猜。
    ds = [float(s["d"]) for s in sh]
    for i, s in enumerate(sh):
        if s.get("tol") is None:
            gaps = [abs(ds[i] - ds[j]) for j in range(len(ds)) if j != i]
            s["tol"] = round(min(0.3, 0.5 * (min(gaps) if gaps else 0.6)), 6)

    # sum_coeff: 设计矩阵里对应的列 (如果有 A)
    if A is not None:
        Am = _as_matrix(A)
        if Am.shape[1] != len(sh) + 1:
            raise JResultError(f"A 的形状 {Am.shape} 与 shells 数 {len(sh)} 不匹配 (应为 N x {len(sh)+1})")
        for i, s in enumerate(sh):
            s.setdefault("sum_coeff", Am[:, i + 1].tolist())

    # z_n: 元胞内**无向键数** = -c_k(fm) (铁磁态 σ_iσ_j ≡ 1), 每磁位点键数 = 2·n_bonds/n_mag
    fm = z_from or (orders[0] if orders else None)
    for s in sh:
        if "z_n" in s and s["z_n"] is not None:
            continue
        sc = s.get("sum_coeff")
        row = None
        if sc is not None and fm in orders:
            row = float(sc[orders.index(fm)])
        s["z_n"] = (round(2.0 * (-row) / n_mag, 6) if row is not None and n_mag else None)

    obj = {
        "schema": SCHEMA,
        "created": datetime.now().isoformat(timespec="seconds"),
        "material": material,
        "lattice": lattice,
        "source": source,
        "confidence": confidence_of(source, residual, cond),
        "spin": float(spin), "S2": S2, "n_mag": int(n_mag),
        "supercell": list(supercell) if supercell is not None else None,
        "orders": orders,
        "A": _as_matrix(A).tolist() if A is not None else None,
        "shells": sh,
        "E_eV": {k: round(float(v), 12) for k, v in E.items()},
        "E0_eV": round(float(E0), 12),
        # J 取 10 位小数: 物理解算会经 J_tilde = J*S^2 往返, 不圆整会留下
        # 2.79 -> 2.790000000000001 这类尾差, 让消费端算出的 J_mat 与原始输入差一个 ulp。
        "J_meV": {k: round(float(v), 10) for k, v in J.items()},
        "J_tilde_eV": {k: round(float(v) * S2, 10) for k, v in J.items()},
        "cond": float(cond),
        "residual_eV": float(residual),
        "notes": list(notes or []),
    }
    obj.update(extra)
    return obj


def _as_matrix(A):
    import numpy as np
    return np.asarray(A, dtype=float)


# --------------------------------------------------------------------------- 校验
def validate(obj, *, strict=True):
    """返回 warnings 列表; strict 时对硬伤 raise JResultError。"""
    w = []
    if not isinstance(obj, dict):
        raise JResultError(f"产物必须是 dict, 收到 {type(obj).__name__}")

    if obj.get("schema") != SCHEMA:
        msg = f"schema 不明: {obj.get('schema')!r} (本代码认 {SCHEMA!r})"
        if strict:
            raise JResultError(msg)
        w.append(msg)

    for key in ("lattice", "source", "confidence", "spin", "n_mag", "shells", "J_meV",
                "cond", "residual_eV"):
        if key not in obj:
            msg = f"缺字段 {key!r}"
            if strict:
                raise JResultError(msg)
            w.append(msg)
    if w and not strict:
        return w

    if obj["source"] not in SOURCES:
        raise JResultError(f"source 非法: {obj['source']!r}")
    if obj["confidence"] not in CONFIDENCE:
        raise JResultError(f"confidence 非法: {obj['confidence']!r}")

    # confidence 必须能被 (source, residual, cond) 复现 —— 防止有人手改成 high
    want = confidence_of(obj["source"], obj["residual_eV"], obj["cond"])
    if obj["confidence"] != want:
        msg = (f"confidence 自相矛盾: 记录 {obj['confidence']!r}, 但按 "
               f"source={obj['source']!r}/residual={obj['residual_eV']:.3g}/"
               f"cond={obj['cond']:.3g} 应为 {want!r}")
        if strict:
            raise JResultError(msg)
        w.append(msg)

    if not obj["shells"]:
        raise JResultError("shells 为空")
    have = {s.get("name") for s in obj["shells"]}
    missing = [k for k in obj["J_meV"] if k not in have]
    if missing:
        raise JResultError(f"J 里的壳 {missing} 不在 shells {sorted(have)} 里")

    if strict:
        r = obj["residual_eV"]
        if not math.isfinite(r) or r > RESIDUAL_FAIL:
            raise JResultError(
                f"残差 {r:.6g} eV 超过硬阈值 {RESIDUAL_FAIL:g} eV —— 反解没收敛/设计矩阵错, "
                f"拒绝写出 (确实要留档请显式 force=True)")
    if obj["residual_eV"] > RESIDUAL_WARN:
        w.append(f"残差 {obj['residual_eV']:.3g} eV 偏大 (> {RESIDUAL_WARN:g}), 可信度已降为 low")
    if obj["confidence"] == "low":
        w.append("confidence=low: 不要直接拿这组 J 去跑 MC 出物理结论")
    return w


# --------------------------------------------------------------------------- 读写
def write(path, obj, *, force=False):
    """校验后写出。返回 warnings。force=True 才允许写出残差不合格的产物。"""
    w = validate(obj, strict=not force)
    if force:
        w += validate(obj, strict=False)
    d = os.path.dirname(os.path.abspath(path))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)
        f.write("\n")
    return w


def read(path, *, strict=True):
    with open(path, encoding="utf-8-sig") as f:
        obj = json.load(f)
    w = validate(obj, strict=strict)
    return obj, w


def format_summary(obj) -> str:
    """一行行的可读摘要 (给 CLI 用)。"""
    lines = [f"J_result [{obj.get('schema')}] {obj.get('material') or '-'} "
             f"lattice={obj.get('lattice')} source={obj.get('source')} "
             f"confidence={obj.get('confidence')}"]
    lines.append(f"  spin={obj.get('spin')} S2={obj.get('S2')} n_mag={obj.get('n_mag')} "
                 f"supercell={obj.get('supercell')}")
    lines.append(f"  cond(A)={obj.get('cond'):.4g}  residual={obj.get('residual_eV'):.3g} eV  "
                 f"E0={obj.get('E0_eV'):+.6f} eV")
    for s in obj.get("shells", []):
        j = obj.get("J_meV", {}).get(s["name"])
        jt = obj.get("J_tilde_eV", {}).get(s["name"])
        pair = f" pair={s['pair']}" if s.get("pair") else ""
        lines.append(f"  {s['name']:<4s} d={s['d']:.4f} A  z_n={s.get('z_n')}  "
                     f"J={j:+.4f} meV  J~={jt:+.6f} eV{pair}")
    return "\n".join(lines)


def load_shells(path_or_json):
    """从 J_result.json (P0.2) **或** magnetic_orders.json (旧的) 或内联 JSON 里取壳定义。

    返回 [{name, d, tol, pair, J_meV?}] —— 供 poscar_to_mc3 --shell-def/--jfile 使用。
    这是 P0.4 的桥: MC 模型必须用**同一套壳**, 否则 J 与键集对不上。
    """
    if isinstance(path_or_json, (list, dict)):
        data = path_or_json
        src = "<inline>"
    else:
        src = path_or_json
        text = path_or_json.strip()
        if text.startswith("{") or text.startswith("["):
            data = json.loads(text)
        else:
            with open(path_or_json, encoding="utf-8-sig") as f:
                data = json.load(f)

    out = []
    if data.get("schema") == SCHEMA or "shells" in data and isinstance(data["shells"], list):
        for s in data["shells"]:
            out.append(dict(name=s["name"], d=float(s["d"]), tol=float(s.get("tol", 0.3)),
                            pair=s.get("pair"),
                            J_meV=(data.get("J_meV") or {}).get(s["name"])))
        if not out:
            raise JResultError(f"{src}: shells 为空")
        return out
    raise JResultError(f"{src}: 认不出壳定义 (既不是 j_result/1 也不是含 shells[] 的 JSON)")


def shells_from_orders_json(path, tol=0.3):
    """把旧的 `magnetic_orders.json` 的 `shells_angstrom` 适配成壳定义列表 (兼容旧产物)。"""
    with open(path, encoding="utf-8-sig") as f:
        data = json.load(f)
    sa = data.get("shells_angstrom")
    if isinstance(sa, dict):
        return [dict(name=k, d=float(v), tol=tol, pair=None, J_meV=None) for k, v in sa.items()]
    if isinstance(sa, (list, tuple)):
        return [dict(name=f"J{i+1}", d=float(v), tol=tol, pair=None, J_meV=None)
                for i, v in enumerate(sa)]
    raise JResultError(f"{path}: 没有 shells_angstrom, 无法适配")
