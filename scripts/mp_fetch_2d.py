#!/usr/bin/env python
"""
mp_fetch_2d.py
访问 Materials Project 官方 REST API, 找二维磁性单层候选, 剥离单层 + 加真空, 写成 VASP POSCAR。

依赖: 只有 requests / numpy (本机已装)。
认证: 环境变量 MP_API_KEY (materialsproject.org/api 获取)。

用法:
  # 0) 只看候选 (不下载结构)
  python mp_fetch_2d.py list --max-hull 0.05 --nsites-max 30

  # 1) 按化学式精确取 (如 CrI3)
  python mp_fetch_2d.py fetch --formula CrI3 --formula VI3 --formula CrCl3

  # 2) 批量: 磁性 + 近稳定 + 小胞, 自动剥单层
  python mp_fetch_2d.py fetch --auto --max-hull 0.05 --nsites-max 30 --limit 60

  # 3) 指定 material_id
  python mp_fetch_2d.py fetch --ids mp-27638 mp-1018793

  # 4) 自检 (不联网): 把已有单层叠成双层再剥回来
  python mp_fetch_2d.py selftest

输出: structures/<material_id>_<formula>.vasp  (POSCAR) + candidates.csv
"""

import argparse
import csv
import json
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(BASE, "structures")
API = "https://api.materialsproject.org/materials/summary/"
FIELDS = ("material_id,formula_pretty,chemsys,is_magnetic,ordering,total_magnetization,"
          "num_magnetic_sites,types_of_magnetic_species,nsites,energy_above_hull,is_stable,"
          "symmetry,structure")


# ------------------------------------------------------------------ POSCAR / 结构
def write_poscar(path, lattice, frac, symbols, comment=""):
    order = []
    for s in symbols:
        if s not in order:
            order.append(s)
    counts = [symbols.count(s) for s in order]
    out = [comment, "1.0"]
    for v in lattice:
        out.append("  %18.10f %18.10f %18.10f" % tuple(v))
    out.append("  " + "  ".join(order))
    out.append("  " + "  ".join(str(c) for c in counts))
    out.append("Direct")
    for s in order:
        for i, si in enumerate(symbols):
            if si == s:
                out.append("  %18.10f %18.10f %18.10f" % tuple(frac[i]))
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


def read_poscar(path):
    lines = [l.rstrip("\n") for l in open(path)]
    scale = float(lines[1].split()[0])
    L = np.array([[float(x) for x in lines[i].split()[:3]] for i in range(2, 5)]) * scale
    tok = lines[5].split()
    if all(_isfloat(t) for t in tok):
        symbols, idx = None, 6
        counts = [int(t) for t in tok]
    else:
        symbols = tok
        counts = [int(t) for t in lines[6].split()]
        idx = 7
    if lines[idx].strip()[:1] in "sS":
        idx += 1
    mode = lines[idx].strip()[:1]
    idx += 1
    nat = sum(counts)
    coords = np.array([[float(x) for x in lines[idx + i].split()[:3]] for i in range(nat)])
    if mode in "cCkK":
        coords = coords @ np.linalg.inv(L)
    if symbols is None:
        symbols = ["X"] * nat
    else:
        symbols = [s for s, c in zip(symbols, counts) for _ in range(c)]
    return L, coords - np.floor(coords), symbols


def _isfloat(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def structure_arrays(st):
    L = np.array(st["lattice"]["matrix"], dtype=float)
    frac = np.array([s["abc"] for s in st["sites"]], dtype=float)
    symbols = [s["species"][0]["element"] for s in st["sites"]]
    return L, frac, symbols


# ------------------------------------------------------------------ 层状分析 / 剥单层
def stacking_axis(L, frac):
    """返回(层堆叠方向轴, 最大层间间隙 Å)。"""
    best = (0, -1.0)
    for ax in range(3):
        p = np.sort(frac[:, ax] % 1.0)
        gaps = np.diff(np.concatenate([p, [p[0] + 1.0]]))
        g = float(gaps.max()) * float(np.linalg.norm(L[ax]))
        if g > best[1]:
            best = (ax, g)
    return best


def split_layers(L, frac, ax, gap_floor=2.0, rel=0.7):
    """沿 ax 轴把原子分成若干层 (周期)。层间 vdW 间隙判据: gap > max(floor, rel*最大间隙),
    避免把层内原子间距 (通常 < 2 A) 误判成层间间隙。"""
    clen = float(np.linalg.norm(L[ax]))
    proj = (frac[:, ax] % 1.0) * clen
    order = np.argsort(proj)
    ps = proj[order]
    n = len(ps)
    if n <= 1:
        return [list(order)]
    gaps = np.diff(np.concatenate([ps, [ps[0] + clen]]))
    thr = max(gap_floor, rel * float(gaps.max()))
    cuts = [i for i in range(n) if gaps[i] > thr]
    if not cuts:
        return [list(order)]
    groups = []
    for a, b in zip(cuts, cuts[1:] + [cuts[0]]):
        idxs, k = [], (a + 1) % n
        while True:
            idxs.append(int(order[k]))
            if k == b:
                break
            k = (k + 1) % n
        groups.append(idxs)
    return groups


def reduced_ratio(symbols):
    """元素比的约简形式, 用于校验剥离层与体相是否同一化学计量。"""
    from collections import Counter
    from math import gcd
    c = Counter(symbols)
    g = 0
    for v in c.values():
        g = gcd(g, v)
    g = g or 1
    return tuple(sorted((k, v // g) for k, v in c.items()))


def exfoliate(L, frac, symbols, mag_elements=(), vacuum=15.0, gap_floor=2.0):
    """把层状块体剥成单层并加真空。返回 (L2, frac2, symbols2, info)。
    info['stoich_ok'] = 剥离层元素比是否与体相一致 (不一致说明这不是干净的 vdW 层)。"""
    ax, gap = stacking_axis(L, frac)
    layers = split_layers(L, frac, ax, gap_floor)
    if len(layers) <= 1:
        return L, frac, symbols, dict(exfoliated=False, n_layers=1, stacking_axis=ax,
                                      stoich_ok=True)

    # 选含磁性元素最多的层; 并列时取原子数最多的
    def score(g):
        mag = sum(1 for i in g if symbols[i] in mag_elements) if mag_elements else 0
        return (mag, len(g))
    layer = max(layers, key=score)
    keep = sorted(layer)

    # 层厚
    axis_dir = L[ax] / np.linalg.norm(L[ax])
    proj = frac[keep] @ L @ axis_dir
    thick = float(proj.max() - proj.min())
    c_new = thick + vacuum

    L2 = L.copy()
    L2[ax] = axis_dir * c_new
    if L2[ax, 2] < 0:            # 规整: 让堆叠轴指向 +z
        L2[ax] = -L2[ax]

    # 重新表达分数坐标: 层中心放到 0.5
    cart = frac[keep] @ L
    center = cart.mean(axis=0)
    cart = cart - center
    cart = cart + 0.5 * L2[ax]          # 层中心沿新 c 到 0.5
    frac2 = cart @ np.linalg.inv(L2)
    frac2[:, ax] %= 1.0
    for a in range(3):
        if a != ax:
            frac2[:, a] %= 1.0
    sym2 = [symbols[i] for i in keep]
    stoich_ok = reduced_ratio(sym2) == reduced_ratio(symbols)
    return L2, frac2, sym2, dict(exfoliated=True, n_layers=len(layers),
                                 stacking_axis=ax, gap_A=gap, thickness_A=thick,
                                 vacuum_A=vacuum, stoich_ok=stoich_ok)


def layered_metrics(L, frac, gap_floor=2.0):
    ax, gap = stacking_axis(L, frac)
    layers = split_layers(L, frac, ax, gap_floor)
    return ax, gap, len(layers)


# ------------------------------------------------------------------ scan
TMS = ["V", "Cr", "Mn", "Fe", "Co", "Ni"]
XS = ["F", "Cl", "Br", "I", "S", "Se", "Te"]


def scan(key, args):
    """按 3d 磁性元素 x 卤素/硫属 扫 chemsys, 保留层状且磁性、近稳定的, 按化学式去重。"""
    seen = {}
    for tm in args.tms:
        for x in args.xs:
            cs = f"{tm}-{x}"
            try:
                docs = query(key, {"chemsys": cs,
                                   "energy_above_hull_max": args.max_hull,
                                   "nsites_max": args.nsites_max}, 100)
            except Exception as e:
                print(f"  {cs}: skip ({type(e).__name__}: {str(e)[:80]})")
                continue
            n_lay = 0
            for d in docs:
                st = d.get("structure")
                if not st:
                    continue
                L, frac, sym = structure_arrays(st)
                ax, gap, nl = layered_metrics(L, frac, args.gap_floor)
                m = d.get("total_magnetization") or 0.0
                mag = d.get("types_of_magnetic_species") or []
                if gap < args.min_gap or (abs(m) < 0.1 and not mag):
                    continue
                n_lay += 1
                keyf = d.get("formula_pretty")
                hull = d.get("energy_above_hull")
                hull = 9.0 if hull is None else hull
                if keyf not in seen or hull < seen[keyf][0]:
                    seen[keyf] = (hull, d)
            print(f"  {cs:6s}: {len(docs):3d} 条, 层状磁性 {n_lay:2d}, 去重累计 {len(seen)}")
    return [d for _, d in sorted(seen.values(), key=lambda t: t[0])]


# ------------------------------------------------------------------ MP API
def api_get(key, params, tries=5):
    import time
    import requests
    last = None
    for k in range(tries):
        try:
            r = requests.get(API, headers={"X-API-KEY": key, "accept": "application/json"},
                             params=params, timeout=90)
            if r.status_code == 200:
                return r.json()
            if r.status_code in (429, 500, 502, 503, 504):
                last = f"{r.status_code}: {r.text[:200]}"
                time.sleep(2.0 * (k + 1))
                continue
            raise SystemExit(f"MP API {r.status_code}: {r.text[:400]}")
        except requests.exceptions.RequestException as e:
            last = str(e)
            time.sleep(2.0 * (k + 1))
    raise SystemExit(f"MP API 请求失败 (已重试 {tries} 次): {last}")


def query(key, params, limit):
    out, skip = [], 0
    while len(out) < limit:
        p = dict(params)
        p["_fields"] = FIELDS
        p["_limit"] = min(200, limit - len(out))
        p["_skip"] = skip
        js = api_get(key, p)
        data = js.get("data") or []
        out.extend(data)
        if len(data) < p["_limit"]:
            break
        skip += p["_limit"]
    return out


def emit(docs, key, vacuum, gap_floor, args):
    os.makedirs(OUTDIR, exist_ok=True)
    rows = []
    for d in docs:
        st = d.get("structure")
        if not st:
            continue
        L, frac, symbols = structure_arrays(st)
        mag_el = [e for e in (d.get("types_of_magnetic_species") or [])]
        if not mag_el:
            mag_el = [e for e in ("V", "Cr", "Mn", "Fe", "Co", "Ni") if e in symbols]
        L2, frac2, sym2, info = exfoliate(L, frac, symbols, mag_el, vacuum, gap_floor)
        mid = str(d["material_id"])
        formula = d.get("formula_pretty", "")
        if not info.get("stoich_ok", True):
            print(f"  skip {mid:12s} {formula:10s} 剥离层元素比与体相不符 "
                  f"({reduced_ratio(sym2)} != {reduced_ratio(symbols)}), 判为非 vdW 层")
            continue
        fn = os.path.join(OUTDIR, f"{mid}_{formula}.vasp")
        comment = (f"{formula} ({mid}) mono | layers={info.get('n_layers')} "
                   f"vac={vacuum}A | ordering={d.get('ordering')} "
                   f"m={d.get('total_magnetization')}")
        write_poscar(fn, L2, frac2, sym2, comment)
        rows.append(dict(material_id=mid, formula=formula, chemsys=d.get("chemsys"),
                         ordering=d.get("ordering"), mag_sites=d.get("num_magnetic_sites"),
                         mag_species=";".join(mag_el), nsites=d.get("nsites"),
                         hull=d.get("energy_above_hull"), n_layers=info.get("n_layers"),
                         stacking_axis=info.get("stacking_axis"), poscar=os.path.basename(fn)))
        print(f"  {mid:12s} {formula:10s} order={str(d.get('ordering')):8s} "
              f"m={d.get('total_magnetization')} hull={d.get('energy_above_hull')} "
              f"layers={info.get('n_layers')} -> {os.path.basename(fn)}")
    csv_path = os.path.join(BASE, "candidates.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["material_id"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n写出 {len(rows)} 个 POSCAR 到 {OUTDIR}")
    print(f"候选清单: {csv_path}")


# ------------------------------------------------------------------ selftest
def selftest():
    src = os.path.join(BASE, "..", "CrI3_J", "POSCAR")
    if not os.path.exists(src):
        src = os.path.join(BASE, "structures")
        print("需要 CrI3 单层 POSCAR 做自检; 未找到, 跳过"); return
    L, frac, symbols = read_poscar(src)
    print(f"读入单层: {len(symbols)} 原子, c = {np.linalg.norm(L[2]):.2f} A")
    # 叠成双层: 把单层复制一份平移半个 c
    Lb = L.copy()
    frac_b = np.vstack([frac, frac + np.array([0, 0, 0.5])])
    sym_b = symbols + symbols
    print(f"人工双层: {len(sym_b)} 原子")
    L2, frac2, sym2, info = exfoliate(Lb, frac_b, sym_b, mag_elements=("Cr",))
    print("剥离结果:", {k: (round(v, 3) if isinstance(v, float) else v) for k, v in info.items()})
    print(f"剥回原子数: {len(sym2)} (期望 {len(symbols)})")
    print(f"新 c = {np.linalg.norm(L2[2]):.2f} A")


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description="Materials Project -> 二维磁性单层 POSCAR")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p1 = sub.add_parser("list")
    p2 = sub.add_parser("fetch")
    for p in (p1, p2):
        p.add_argument("--max-hull", type=float, default=0.05)
        p.add_argument("--nsites-max", type=int, default=30)
        p.add_argument("--limit", type=int, default=60)
        p.add_argument("--vacuum", type=float, default=15.0)
        p.add_argument("--gap-floor", type=float, default=2.0)
        p.add_argument("--min-gap", type=float, default=3.0,
                       help="层间 vdW 间隙下限 (A), 用于筛层状材料")
    p2.add_argument("--formula", action="append", default=[])
    p2.add_argument("--ids", nargs="*", default=[])
    p2.add_argument("--auto", action="store_true")

    p3 = sub.add_parser("scan")
    p3.add_argument("--tms", nargs="*", default=TMS)
    p3.add_argument("--xs", nargs="*", default=XS)
    p3.add_argument("--max-hull", type=float, default=0.05)
    p3.add_argument("--nsites-max", type=int, default=24)
    p3.add_argument("--vacuum", type=float, default=15.0)
    p3.add_argument("--gap-floor", type=float, default=2.0)
    p3.add_argument("--min-gap", type=float, default=3.0)

    sub.add_parser("selftest")

    args = ap.parse_args()

    if args.cmd == "selftest":
        selftest(); return

    key = os.environ.get("MP_API_KEY")
    if not key:
        raise SystemExit("未设置 MP_API_KEY 环境变量")

    if args.cmd == "list":
        docs = query(key, {"num_magnetic_sites_min": 1, "energy_above_hull_max": args.max_hull,
                           "nsites_max": args.nsites_max, "nelements_min": 2, "nelements_max": 3},
                     args.limit)
        rows = []
        for d in docs:
            st = d.get("structure")
            if not st:
                continue
            L, frac, sym = structure_arrays(st)
            ax, gap, nl = layered_metrics(L, frac, args.gap_floor)
            rows.append((gap, nl, d))
        rows.sort(key=lambda r: -(r[0]))
        shown = 0
        print(f"{'material_id':12s} {'formula':12s} {'chemsys':12s} {'order':6s} "
              f"{'m':>8s} {'hull':>7s} {'gapA':>6s} layers")
        for gap, nl, d in rows:
            if gap < args.min_gap:
                continue
            print(f"  {d['material_id']:12s} {str(d.get('formula_pretty','')):12s} "
                  f"{str(d.get('chemsys')):12s} {str(d.get('ordering')):6s} "
                  f"{d.get('total_magnetization'):8.3f} {d.get('energy_above_hull'):7.4f} "
                  f"{gap:6.2f} {nl}")
            shown += 1
        print(f"\n{len(rows)} 个磁性候选, 其中 {shown} 个层状间隙 >= {args.min_gap} A (二维候选)")
        return

    # scan: 3d 磁性元素 x 卤素/硫属 系统扫描
    if args.cmd == "scan":
        docs = scan(key, args)
        print(f"\n扫描得到 {len(docs)} 个去重后的层状磁性候选, 开始剥离单层...")
        emit(docs, key, args.vacuum, args.gap_floor, args)
        return

    # fetch
    if args.formula:
        docs = []
        for f in args.formula:
            docs += query(key, {"formula": f}, 50)
    elif args.ids:
        docs = query(key, {"material_ids": ",".join(args.ids)}, len(args.ids) + 5)
    elif args.auto:
        docs = query(key, {"num_magnetic_sites_min": 1, "energy_above_hull_max": args.max_hull,
                           "nsites_max": args.nsites_max, "nelements_min": 2,
                           "nelements_max": 3}, args.limit)
    else:
        raise SystemExit("fetch 需要 --formula / --ids / --auto 之一")

    # 去重
    seen, uniq = set(), []
    for d in docs:
        if d["material_id"] not in seen:
            seen.add(d["material_id"]); uniq.append(d)
    print(f"取到 {len(uniq)} 个结构, 开始剥离单层...")
    emit(uniq, key, args.vacuum, args.gap_floor, args)


if __name__ == "__main__":
    main()
