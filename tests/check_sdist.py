#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""校验 sdist 里该有的东西都在（发布闸门，配合 .github/workflows/ci.yml 的 package 作业）。

用法::

    python tests/check_sdist.py dist/2d_magnet_pipeline-0.1.0.tar.gz
    python tests/check_sdist.py "dist\\*.tar.gz"      # 通配符由 shell 展开（Windows PowerShell 不展开 ⇒ 本脚本自己兜底）

为什么需要它：`wheel` 只装 `run_2d_magnet.py` + `scripts/`（见 pyproject.toml），
文档、手册、agent skill、样例产物、测试夹具、扩展骨架全靠 `MANIFEST.in` 进 sdist。
少了 MANIFEST.in 时，`pip install` 出来的包能跑，但 clone 里有的手册/样例在 sdist 里消失 —— 这类
「clone 有、sdist 没有」的发布事故不会被任何测试抓到，所以单独做成一道闸门。
"""
from __future__ import annotations

import glob
import os
import sys
import tarfile

# 必须是 sdist 里存在的条目（用 endswith 匹配，容忍顶层的 <name>-<version>/ 前缀）
REQUIRED = [
    "README.md",
    "LICENSE",
    "requirements.txt",
    "pyproject.toml",
    "MANIFEST.in",
    "KNOWN_ISSUES.md",
    "WORKFLOW.md",
    "NEW_TYPES_PLAN.md",
    "P4_2_KAGOME_RUNBOOK.md",
    "DMI_MTI_MEMO.md",
    "CHANGELOG.md",
    "CITATION.cff",
    "RELEASE_NOTES_v0.1.0.md",
    "run_2d_magnet.py",
    "capabilities.json",
    "scripts/capabilities.py",
    "scripts/make_release.py",
    "scripts/literature_J.json",
    "scripts/magnetic_orders.py",
    "scripts/verify_states.py",
    "scripts/readme_to_obsidian.py",
    "README.zh.md",
    "README.obsidian.md",
    "README.zh.obsidian.md",
    "tests/smoke_test.py",
    "tests/fixtures/MnPS3_POSCAR",
    "dsh/2d-magnet/SKILL.md",
    "templates/material_preset.json",
    "examples/README.md",
    "examples/la2cuo4_square/README.md",
    "examples/kagome_fe3sn_dryrun/02_supercell/cases/ingest_report.md",
]

# 明确**不该**出现的（本地脏物 / 缓存 / 探针 / 笔记软件工作区）
FORBIDDEN_SUFFIXES = ("__pycache__/", ".pyc", ".pyo", ".egg-info/", "_tmp_", "_probe_", ".git/", ".obsidian/", ".trash/")


def resolve(spec: str) -> str:
    """把可能的通配符/目录展开成唯一的 tar.gz 路径。"""
    if os.path.isdir(spec):
        hits = sorted(glob.glob(os.path.join(spec, "*.tar.gz")))
    elif any(ch in spec for ch in "*?["):
        hits = sorted(glob.glob(spec))
    else:
        hits = [spec]
    hits = [h for h in hits if h.endswith(".tar.gz")]
    if not hits:
        raise SystemExit("[error] 找不到 sdist: %s（先跑 `python -m build`）" % spec)
    if len(hits) > 1:
        raise SystemExit("[error] 匹配到多个 sdist, 请只留一个: %s" % hits)
    return hits[0]


def main(argv) -> int:
    if len(argv) != 2:
        raise SystemExit("用法: python tests/check_sdist.py <sdist.tar.gz|dist 目录|通配符>")
    path = resolve(argv[1])
    with tarfile.open(path) as tf:
        names = [n for n in tf.getnames() if not n.endswith("/")]
    print("sdist: %s" % path)
    print("文件数: %d" % len(names))

    miss = [w for w in REQUIRED if not any(n.endswith(w) for n in names)]
    dirty = sorted({n for n in names if n.endswith(FORBIDDEN_SUFFIXES) or "/__pycache__/" in n})

    for w in REQUIRED:
        print("  %s %s" % ("[ok]  " if w not in miss else "[缺]  ", w))
    if miss:
        print("\n[FAIL] 缺少 %d 项 —— 检查 MANIFEST.in（文档/样例/skill/测试夹具是不是没进去）" % len(miss))
    if dirty:
        print("\n[FAIL] sdist 里混入了不该有的东西:")
        for d in dirty[:20]:
            print("   ", d)
    if miss or dirty:
        return 1
    print("\n[OK] sdist 内容完整（文档 / 手册 / skill / 样例 / 测试夹具 / 骨架 都在），且无缓存与临时物")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
