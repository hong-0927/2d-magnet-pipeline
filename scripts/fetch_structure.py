#!/usr/bin/env python3
"""结构 provider 链: COD / OQMD / AFLOW / JARVIS (Materials Project 由 mp_fetch_2d.py 单独处理)。

给一个化学式, 依次尝试各 provider 下载结构并转成 POSCAR(纯文本)。三种解析器统一收敛到
POSCAR: CIF(COD)、OPTIMADE(OQMD / JARVIS)、AFLUX(AFLOW)。HTTP 请求隔离在 `_http_get`
(便于测试 monkeypatch); 转换器为纯函数(不触网), 可独立单测。

用法:
  python fetch_structure.py MnPS3 --provider cod   --out POSCAR.mp
  python fetch_structure.py MnPS3 --provider auto  --out POSCAR.mp   # cod -> oqmd -> aflow -> jarvis

无 API key 即可用(MP 除外)。说明: 各 provider 的公开 REST 端点可能随时间调整, 解析均为
best-effort; 失败会被链捕获并尝试下一家。
"""
import os
import re
import sys


# ------------------------------------------------------------------ 基础工具
def formula_elements(formula):
    """'MnPS3' -> ['Mn', 'P', 'S'] (去重, 保持出现顺序, 支持双字母元素)。"""
    toks = re.findall(r"[A-Z][a-z]?", formula)
    seen = []
    for t in toks:
        if t and t not in seen:
            seen.append(t)
    return seen or [formula]


def _frac(s):
    """'1/2' -> 0.5, '2' -> 2.0, '0.5' -> 0.5。"""
    s = s.strip()
    if "/" in s:
        a, _, b = s.partition("/")
        return float(a) / float(b)
    return float(s)


# ------------------------------------------------------------------ CIF 解析
def _cif_tokens(text):
    """CIF 词法: 产出 ('tag', name) / ('loop', None) / ('value', val)。"""
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if c in " \t\r\n":
            i += 1
            continue
        if c == "#":  # 注释到行尾
            j = text.find("\n", i)
            i = n if j == -1 else j + 1
            continue
        if c == "_" or text[i:i + 5] == "loop_":
            m = re.match(r"(_\S+|loop_)", text[i:])
            tok = m.group(1)
            i += len(tok)
            yield ("loop", None) if tok == "loop_" else ("tag", tok)
            continue
        if c == ";" and (i == 0 or text[i - 1] in "\r\n"):
            j = text.find("\n;", i + 1)
            if j == -1:
                yield ("value", text[i + 1:].strip())
                i = n
            else:
                yield ("value", text[i + 1:j].strip())
                i = j + 2
            continue
        if c in "\"'":
            q = c
            j = text.find(q, i + 1)
            if j == -1:
                yield ("value", text[i + 1:].strip())
                i = n
            else:
                yield ("value", text[i + 1:j].strip())
                i = j + 1
            continue
        m = re.match(r"(\S+)", text[i:])
        yield ("value", m.group(1))
        i += len(m.group(1))


def parse_cif(text):
    """最小 CIF 读取: 标量 tag -> str/None; loop 列 tag -> list[str] (按行对齐)。"""
    toks = list(_cif_tokens(text))
    data = {}
    i, n = 0, len(toks)
    while i < n:
        kind, val = toks[i]
        if kind == "loop":
            tags = []
            i += 1
            while i < n and toks[i][0] == "tag":
                tags.append(toks[i][1])
                i += 1
            vals = []
            while i < n and toks[i][0] == "value":
                vals.append(toks[i][1])
                i += 1
            ncol = len(tags)
            if ncol:
                nrow = len(vals) // ncol
                for ci, tag in enumerate(tags):
                    data[tag] = [vals[r * ncol + ci] for r in range(nrow)]
            continue
        if kind == "tag":
            tag = val
            if i + 1 < n and toks[i + 1][0] == "value":
                data[tag] = toks[i + 1][1]
                i += 2
            else:
                data[tag] = None
                i += 1
            continue
        i += 1  # 游离 value, 跳过
    return data


def _cell_to_lattice(a, b, c, alpha, beta, gamma):
    import math
    al, be, ga = map(math.radians, (alpha, beta, gamma))
    ca, cb, cg = math.cos(al), math.cos(be), math.cos(ga)
    sg = math.sin(ga)
    ax = [a, 0.0, 0.0]
    bx = [b * cg, b * sg, 0.0]
    c1, c2 = c * cb, c * (ca - cb * cg) / sg
    c3 = math.sqrt(max(0.0, c * c - c1 * c1 - c2 * c2))
    return [ax, bx, [c1, c2, c3]]


def _parse_component(expr):
    """'1/2+x-y' -> (0.5, {'x':1.0,'y':-1.0,'z':0.0})。"""
    expr = expr.replace(" ", "")
    const = 0.0
    coeffs = {"x": 0.0, "y": 0.0, "z": 0.0}
    for tok in re.findall(r"[+-]?[^+-]+", expr):
        sign = 1.0
        if tok.startswith("-"):
            sign, tok = -1.0, tok[1:]
        elif tok.startswith("+"):
            tok = tok[1:]
        if not tok:
            continue
        matched = next((v for v in ("x", "y", "z") if v in tok), None)
        if matched:
            co = tok.replace(matched, "")
            co = 1.0 if co in ("", "+") else (-1.0 if co == "-" else _frac(co))
            coeffs[matched] += sign * co
        else:
            const += sign * _frac(tok)
    return const, coeffs


def _sym_matrix(op):
    """CIF 对称操作串 -> (M(3x3), t(3)); 解析失败返回 None。"""
    parts = [p for p in op.split(",") if p.strip()]
    if len(parts) != 3:
        return None
    M = [[0.0] * 3 for _ in range(3)]
    t = [0.0] * 3
    for i, part in enumerate(parts):
        const, coeffs = _parse_component(part)
        t[i] = const
        for j, v in enumerate(("x", "y", "z")):
            M[i][j] = coeffs[v]
    return M, t


def _expand_sym(sites, mats):
    out, seen = [], set()
    for spec, xyz in sites:
        for M, t in mats:
            new = [M[i][0] * xyz[0] + M[i][1] * xyz[1] + M[i][2] * xyz[2] + t[i] for i in range(3)]
            new = [round(v % 1.0, 5) for v in new]
            key = (spec, new[0], new[1], new[2])
            if key not in seen:
                seen.add(key)
                out.append((spec, new))
    return out


def _write_poscar(species_coords, lattice, title="structure"):
    """species_coords: [(spec, [fx,fy,fz])] (分数坐标); lattice: 3x3 行向量 a,b,c。"""
    order = []
    for spec, _ in species_coords:
        if spec not in order:
            order.append(spec)
    counts = [sum(1 for s, _ in species_coords if s == spec) for spec in order]
    lines = [title, "1.0"]
    lines += ["  ".join(f"{x:.12f}" for x in v) for v in lattice]
    lines.append("  ".join(order))
    lines.append("  ".join(str(c) for c in counts))
    lines.append("Direct")
    for spec in order:
        for s, xyz in species_coords:
            if s == spec:
                lines.append("  ".join(f"{x:.12f}" for x in xyz))
    return "\n".join(lines) + "\n"


def cif_to_poscar(cif_text, title="structure"):
    """CIF 文本 -> POSCAR 文本。自动按对称操作展开不对称单元。"""
    d = parse_cif(cif_text)
    need = ["_cell_length_a", "_cell_length_b", "_cell_length_c",
            "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma",
            "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"]
    missing = [k for k in need if k not in d]
    if missing:
        raise ValueError(f"CIF 缺字段: {missing}")
    lattice = _cell_to_lattice(float(d["_cell_length_a"]), float(d["_cell_length_b"]),
                               float(d["_cell_length_c"]), float(d["_cell_angle_alpha"]),
                               float(d["_cell_angle_beta"]), float(d["_cell_angle_gamma"]))
    sym = d.get("_atom_site_type_symbol") or d.get("_atom_site_label")
    fx, fy, fz = d["_atom_site_fract_x"], d["_atom_site_fract_y"], d["_atom_site_fract_z"]
    sites = [(sym[i], [float(fx[i]), float(fy[i]), float(fz[i])]) for i in range(len(sym))]
    ops = d.get("_symmetry_equiv_pos_as_xyz")
    if ops:
        mats = [m for m in (_sym_matrix(o) for o in ops) if m]
        if mats:
            sites = _expand_sym(sites, mats)
    return _write_poscar(sites, lattice, title)


# ------------------------------------------------------------------ OPTIMADE
def optimade_to_poscar(attrs, title="structure"):
    """OPTIMADE structures entry.attributes -> POSCAR 文本。"""
    import numpy as np
    lat = attrs["lattice_vectors"]
    cart = attrs["cartesian_site_positions"]
    species = attrs["species_at_sites"]
    L = np.array(lat, dtype=float)
    inv = np.linalg.inv(L)
    frac = (np.array(cart, dtype=float) @ inv) % 1.0
    sites = [(species[i], frac[i].tolist()) for i in range(len(species))]
    return _write_poscar(sites, lat, title)


# ------------------------------------------------------------------ AFLOW
def aflow_to_poscar(entry, title="structure"):
    """aflowlib entry JSON -> POSCAR 文本。"""
    import numpy as np
    lat = entry["lattice"]
    species = entry.get("atom_names") or entry.get("atom_species")
    if not species:
        raise ValueError("AFLOW entry 缺 atom_names")
    if "positions_fractional" in entry:
        frac = entry["positions_fractional"]
    elif "positions_cartesian" in entry:
        L = np.array(lat, dtype=float)
        frac = ((np.array(entry["positions_cartesian"], dtype=float) @ np.linalg.inv(L)) % 1.0).tolist()
    else:
        raise ValueError("AFLOW entry 缺 positions")
    sites = [(species[i], frac[i]) for i in range(len(species))]
    return _write_poscar(sites, lat, title)


# ------------------------------------------------------------------ HTTP + provider
def _http_get(url, **kw):
    import requests
    kw.setdefault("timeout", 60)
    r = requests.get(url, **kw)
    r.raise_for_status()
    return r


def fetch_cod(formula, _get=_http_get):
    """COD: result.php 搜索 -> 取第一条 CIF -> cif_to_poscar。"""
    r = _get("https://www.crystallography.net/cod/result.php",
             params={"formula": formula, "format": "json"})
    data = r.json()
    entries = data if isinstance(data, list) else (data.get("data") or data.get("results") or [])
    cif_url = None
    for e in entries:
        f = None
        if isinstance(e, dict):
            f = e.get("file") or e.get("url") or e.get("href")
        elif isinstance(e, list) and e:
            f = e[0] if not str(e[0]).startswith("http") else e[0]
        if f:
            cif_url = f if str(f).startswith("http") else "https://www.crystallography.net/cod/" + str(f)
            break
    if not cif_url:
        raise RuntimeError(f"COD 无 {formula} 候选")
    return cif_to_poscar(_get(cif_url).text, title=formula)


OPTIMADE_BASES = {
    "oqmd": "https://oqmd.org/optimade/v1/structures",
    "jarvis": "https://jarvis.nist.gov/optimade/jarvisdft/v1/structures",
}


def fetch_optimade(base, formula, _get=_http_get):
    els = formula_elements(formula)
    filt = "elements HAS ALL " + ", ".join(f'"{e}"' for e in els)
    r = _get(base, params={"filter": filt, "page_limit": 5,
                           "response_fields": ("lattice_vectors,cartesian_site_positions,"
                                               "species_at_sites,nelements,chemical_formula_reduced")})
    data = r.json().get("data") or []
    for entry in data:
        attrs = entry.get("attributes", entry)
        if attrs.get("nelements") == len(els):
            return optimade_to_poscar(attrs, title=formula)
    if data:
        return optimade_to_poscar(data[0].get("attributes", data[0]), title=formula)
    raise RuntimeError(f"{base} 无 {formula} 候选")


def fetch_oqmd(formula, _get=_http_get):
    return fetch_optimade(OPTIMADE_BASES["oqmd"], formula, _get=_get)


def fetch_jarvis(formula, _get=_http_get):
    return fetch_optimade(OPTIMADE_BASES["jarvis"], formula, _get=_get)


def fetch_aflow(formula, _get=_http_get):
    """AFLOW: AFLUX 搜索 -> aurl 取几何 JSON -> aflow_to_poscar。"""
    spec = ",".join(formula_elements(formula))
    r = _get("https://aflow.org/API/aflux/",
             params={"species": f"({spec})", "size": "4", "format": "json", "paging": "0"})
    data = r.json().get("data") or []
    last = None
    for d in data:
        aurl = (d or {}).get("aurl") or (d or {}).get("auid")
        if not aurl:
            continue
        try:
            entry = _get("https://aflow.org/" + aurl + "/", params={"format": "json"}).json()
            return aflow_to_poscar(entry, title=formula)
        except Exception as e:  # noqa: BLE001 逐条尝试
            last = e
    raise RuntimeError(f"AFLOW 解析失败 ({formula}): {last}")


PROVIDERS = {
    "cod": fetch_cod,
    "oqmd": fetch_oqmd,
    "aflow": fetch_aflow,
    "jarvis": fetch_jarvis,
}
CHAIN = ["cod", "oqmd", "aflow", "jarvis"]


def fetch_chain(formula, providers=None, _get=_http_get):
    """依次尝试 providers(默认 CHAIN), 返回 (provider_name, poscar_text)。全失败抛 RuntimeError。"""
    providers = providers or CHAIN
    errs = []
    for p in providers:
        fn = PROVIDERS[p]
        try:
            return p, fn(formula, _get=_get)
        except Exception as e:  # noqa: BLE001
            errs.append(f"{p}: {e}")
    raise RuntimeError(" / ".join(errs))


# ------------------------------------------------------------------ CLI
def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="结构 provider 链 (COD/OQMD/AFLOW/JARVIS)")
    ap.add_argument("formula")
    ap.add_argument("--provider", choices=["auto"] + CHAIN, default="auto")
    ap.add_argument("--out", required=True, help="输出 POSCAR 路径")
    args = ap.parse_args(argv)
    providers = CHAIN if args.provider == "auto" else [args.provider]
    name, poscar = fetch_chain(args.formula, providers)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(poscar)
    print(f"[fetch] {name} -> {args.out} ({len(poscar.splitlines()) - 8} 原子)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
