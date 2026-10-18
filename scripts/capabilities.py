#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""capabilities.json 的校验器与展示器 —— 机器可读的能力边界 (复审 P1-5)。

为什么需要这个文件
------------------
这个仓库的能力边界原本只写在 README 的「已知局限」和 SKILL.md 的散文里。对**人**够用,
对**agent** 不够: agent 会读 README, 但更倾向于读代码/跑命令; 于是「没有真实 DFT 能量
就不要报 T_N」这类约束, 只能靠它自己记得住。`capabilities.json` 把边界变成一份可被
程序读取、可被 CI 校验的**数据**:

    python scripts/capabilities.py --check     # 校验 (stdout 逐条 [OK]/[FAIL])
    python scripts/capabilities.py --show      # 人读的能力表
    python scripts/capabilities.py --check --json other.json

`--check` 拦什么 (这些错都**不报错但很危险**)
---------------------------------------------
1. 引用了不存在的脚本 (能力声明与代码漂移);
2. 某条能力把 `physical_verified` 写成 `true` —— 本版本没有真实 DFT 结果进仓库,
   把它改 true 就等于对外宣称「这是材料物理结论」。要翻成 true, 必须:
   真实 OUTCAR + 实验/文献对照 + 把 examples/ 的 source 从 demo 换掉, 三者缺一不可;
3. `wrong_why`: level 是 coefficient/estimate/input 却没有 `why` (读者不知道该怎么用);
4. `--lattice auto` 之外少一个晶格 / 多一个没实现的晶格 (与 run_2d_magnet.py 的
   choices 对齐);
5. `release.physical_verified` 与 `pyproject.toml` 的 version 对不上 (版本漂移);
6. `gates` 里的脚本不存在, 或 `fails_closed` 不是 bool。

设计约束
--------
- 只读: 不修改任何文件; 不做网络访问。
- 纯标准库: 与 check_sdist.py 一致, 不引入新依赖。
- 任何一条 FAIL 都 `return 1`; 全过才 `return 0`。
"""

import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DEFAULT_JSON = os.path.join(ROOT, "capabilities.json")

SCHEMA = "capabilities/1"
LEVELS = ("coefficient", "estimate", "input", "not_extracted", "not_computed")
# 本版本一律 false: 没有真实 DFT 结果与实验对照进仓库
MUST_BE_FALSE = ("physical_verified",)
# run_2d_magnet.py 的 --lattice choices (去掉 auto)
EXPECTED_LATTICES = ("honeycomb", "triangular", "cluster", "rectangular",
                     "square", "rect_centered", "kagome")


def _version_from_pyproject():
    p = os.path.join(ROOT, "pyproject.toml")
    try:
        with open(p, encoding="utf-8") as f:
            txt = f.read()
    except OSError:
        return None
    m = re.search(r'(?m)^\s*version\s*=\s*"([^"]+)"', txt)
    return m.group(1) if m else None


def _exists(rel):
    return os.path.exists(os.path.join(ROOT, rel.replace("/", os.sep)))


def check(path=DEFAULT_JSON, verbose=True):
    """返回 (失败条数, 明细列表)。"""
    fails = []

    def ok(msg):
        if verbose:
            print(f"  [OK] {msg}")

    def bad(msg):
        fails.append(msg)
        if verbose:
            print(f"  [FAIL] {msg}")

    try:
        with open(path, encoding="utf-8") as f:
            cap = json.load(f)
    except OSError as e:
        bad(f"读不到 {path}: {e}")
        return len(fails), fails
    except json.JSONDecodeError as e:
        bad(f"JSON 解析失败 ({path}): {e}")
        return len(fails), fails
    ok(f"JSON 可解析: {os.path.relpath(path, ROOT)}")

    if cap.get("schema") != SCHEMA:
        bad(f"schema 应为 {SCHEMA!r}, 实得 {cap.get('schema')!r}")
    else:
        ok(f"schema = {SCHEMA}")
    if not str(cap.get("note") or "").strip():
        bad("缺 note (说明这份表该怎么用)")
    else:
        ok("note 非空")

    # ---- 晶格 ----
    lats = cap.get("lattices") or {}
    missing = [k for k in EXPECTED_LATTICES if k not in lats]
    extra = [k for k in lats if k not in EXPECTED_LATTICES]
    if missing:
        bad(f"lattices 缺 {missing} (与 run_2d_magnet.py 的 --lattice choices 对齐)")
    if extra:
        bad(f"lattices 多出 {extra} (driver 没有这条分支)")
    if not missing and not extra:
        ok(f"lattices 与 driver 的 {len(EXPECTED_LATTICES)} 类晶格一一对齐")
    for name, d in sorted(lats.items()):
        tag = f"lattices.{name}"
        for key in ("builder", "solver"):
            rel = d.get(key)
            if not rel or not _exists(str(rel)):
                bad(f"{tag}.{key} 指向的脚本不存在: {rel!r}")
        for rel in d.get("solver_extra") or []:
            if not _exists(str(rel)):
                bad(f"{tag}.solver_extra 指向的脚本不存在: {rel!r}")
        if d.get("physical_verified") is not False:
            bad(f"{tag}.physical_verified 必须显式 false (本版本没有真实 DFT 结果)")
        if not str(d.get("why") or "").strip():
            bad(f"{tag}.why 为空 (读者不知道该晶格的能力边界在哪)")
        shells = d.get("shells")
        if not isinstance(shells, int) or shells < 1:
            bad(f"{tag}.shells 应为 >=1 的整数, 实得 {shells!r}")
    ok(f"lattices 共 {len(lats)} 条: builder/solver 存在, physical_verified 全 false")

    # ---- 特殊链 ----
    for name, d in sorted((cap.get("chains") or {}).items()):
        tag = f"chains.{name}"
        for key in ("builder", "solver"):
            rel = d.get(key)
            if not rel or not _exists(str(rel)):
                bad(f"{tag}.{key} 指向的脚本不存在: {rel!r}")
        if d.get("physical_verified") is not False:
            bad(f"{tag}.physical_verified 必须显式 false")
    ok(f"chains 共 {len(cap.get('chains') or {})} 条")

    # ---- 物理量 ----
    qs = cap.get("quantities") or {}
    if not qs:
        bad("quantities 为空")
    for name, d in sorted(qs.items()):
        tag = f"quantities.{name}"
        lvl = d.get("level")
        if lvl not in LEVELS:
            bad(f"{tag}.level={lvl!r} 不在 {LIST_LEVELS}")
        if lvl != "not_computed" and lvl != "not_extracted" and not str(d.get("why") or "").strip():
            bad(f"{tag}.why 为空 (level={lvl} 必须解释怎么用/怎么错)")
        if lvl in ("not_computed", "not_extracted") and not str(d.get("why") or "").strip():
            bad(f"{tag}.why 为空 (说明为什么不算/不提取)")
        for key in MUST_BE_FALSE:
            if d.get(key) is not False:
                bad(f"{tag}.{key} 必须显式 false (改 true 等于宣称这是材料物理结论)")
    ok(f"quantities 共 {len(qs)} 条: level 合法, physical_verified 全 false")

    # ---- 不支持 ----
    ns = cap.get("not_supported") or []
    if not ns:
        bad("not_supported 为空 (至少列出结构弛豫/真实 DFT 能量/screening)")
    for i, d in enumerate(ns):
        if not str(d.get("item") or "").strip():
            bad(f"not_supported[{i}] 缺 item")
        if not str(d.get("why") or "").strip():
            bad(f"not_supported[{i}] ({d.get('item')!r}) 缺 why")
    ok(f"not_supported 共 {len(ns)} 条, 每条都有 item + why")

    # ---- 闸门 ----
    gs = cap.get("gates") or []
    if not gs:
        bad("gates 为空 (verify_states / ingest / S² / 壳秩 至少各一条)")
    for i, d in enumerate(gs):
        tag = f"gates[{i}]"
        if not str(d.get("name") or "").strip():
            bad(f"{tag} 缺 name")
        rel = d.get("script")
        if not rel or not _exists(str(rel)):
            bad(f"{tag}.script 指向的脚本不存在: {rel!r}")
        if not isinstance(d.get("fails_closed"), bool):
            bad(f"{tag}.fails_closed 应为 bool, 实得 {d.get('fails_closed')!r}")
    ok(f"gates 共 {len(gs)} 条: 脚本存在, fails_closed 是 bool")

    # ---- release 与版本 ----
    rel = cap.get("release") or {}
    if rel.get("physical_verified") is not False:
        bad("release.physical_verified 必须显式 false")
    pv = _version_from_pyproject()
    if pv and rel.get("version") != pv:
        bad(f"release.version={rel.get('version')!r} 与 pyproject.toml 的 {pv!r} 不一致")
    elif pv:
        ok(f"release.version = {pv} (与 pyproject.toml 一致)")

    if not fails:
        ok(f"全部通过 ({len(EXPECTED_LATTICES)} 类晶格 / {len(qs)} 个物理量 / "
           f"{len(ns)} 条不支持 / {len(gs)} 道闸门)")
    return len(fails), fails


LIST_LEVELS = " / ".join(LEVELS)


def show(path=DEFAULT_JSON):
    with open(path, encoding="utf-8") as f:
        cap = json.load(f)
    print(f"能力边界 (schema {cap.get('schema')}, 版本 {(cap.get('release') or {}).get('version')})")
    print("=" * 78)
    print(f"{'晶格':<14}{'壳':<4}{'非共线':<8}{'物理实证':<10}builder")
    for name, d in sorted((cap.get("lattices") or {}).items()):
        print(f"{name:<14}{d.get('shells', '?'):<4}{str(d.get('noncollinear')):<8}"
              f"{str(d.get('physical_verified')):<10}{d.get('builder')}")
    print("\n物理量:")
    for name, d in sorted((cap.get("quantities") or {}).items()):
        print(f"  {name:<26}level={d.get('level'):<14}physical_verified={d.get('physical_verified')}")
        why = str(d.get("why") or "")
        print(f"    {why[:150]}{'…' if len(why) > 150 else ''}")
    print("\n不支持:")
    for d in cap.get("not_supported") or []:
        print(f"  - {d.get('item')}: {str(d.get('why'))[:110]}")
    print("\n闸门 (fail-close):")
    for d in cap.get("gates") or []:
        print(f"  - {d.get('name')} [{d.get('script')}] fails_closed={d.get('fails_closed')}")
    print("\n结论: 本版本所有物理量都是系数层/估算层, 见 release:",
          (cap.get("release") or {}).get("why"))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="校验/展示 capabilities.json (机器可读的能力边界)",
        epilog="例: python scripts/capabilities.py --check")
    ap.add_argument("--json", default=DEFAULT_JSON, help="capabilities.json 路径")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true", help="逐条校验 (默认)")
    g.add_argument("--show", action="store_true", help="打印人读的能力表")
    ap.add_argument("--quiet", action="store_true", help="只打印 FAIL")
    args = ap.parse_args(argv)

    if args.show:
        return show(args.json)
    if not os.path.exists(args.json):
        print(f"[FAIL] 找不到 {args.json}", file=sys.stderr)
        return 1
    print(f"capabilities.json 校验: {os.path.relpath(args.json, ROOT)}")
    n, fails = check(args.json, verbose=not args.quiet)
    if fails:
        print(f"\n[FAIL] {n} 项不通过:", file=sys.stderr)
        for m in fails:
            print(f"  - {m}", file=sys.stderr)
        return 1
    print("\n[OK] capabilities.json 自洽, 与代码/版本/闸门一致")
    return 0


if __name__ == "__main__":
    sys.exit(main())
