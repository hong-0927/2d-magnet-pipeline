#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""把 `README.md` 转成 Obsidian 友好的备份 `README.obsidian.md`。

为什么单独做一个脚本而不是手抄一份：
手抄的备份一定会和 README 漂移（README 是本仓最常改的文件）。这里做的是**确定性**转换——
正文与 `README.md` **逐字一致**，只在前面加 Obsidian 才吃得到的东西（properties / 内链目录 /
mermaid 管线图 / callout），并把 `README.md` 的 sha256 写进 frontmatter；`--check` 用来当闸门
（README 改了但忘了重新生成 ⇒ 非零退出，smoke §29 会跑它）。

用法::

    python scripts/readme_to_obsidian.py                 # 生成 README.obsidian.md
    python scripts/readme_to_obsidian.py --check         # 只校验是否同步（CI/smoke 用）
    python scripts/readme_to_obsidian.py --stdout        # 打到标准输出，不落盘
    python scripts/readme_to_obsidian.py --selftest

在 Obsidian 里怎么用：**把 `2d_magnet_pipeline` 整个目录当 vault 打开**（Open folder as vault），
再打开 `README.obsidian.md`。这样正文里的相对链接（`KNOWN_ISSUES.md` / `WORKFLOW.md` / …
以及 `examples/`）都能解析；`dsh/2d-magnet/SKILL.md` 也在同一个 vault 里，可以直接跳。
`.md` 之外的脚本/夹具在 Obsidian 里是纯文本条目，看代码建议用编辑器。
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "README.md")
OUT = os.path.join(ROOT, "README.obsidian.md")

FRONT_KEYS = ("title", "aliases", "tags", "source", "source_sha256", "generated_by")

PIPELINE_MERMAID = """```mermaid
flowchart TD
    F["fetch<br/>mp_fetch_2d.py / fetch_structure.py"] --> B["build<br/>build_*_cases.py → POSCAR + magmom.txt"]
    B --> D["dft<br/>assemble_dft[_tri].py → INCAR / KPOINTS / POTCAR.spec + runs/&lt;序&gt;/"]
    D --> V(["外部机时：VASP（本仓库只导出输入，不代提交作业）"])
    V --> I["ingest<br/>ingest_dft.py：11 条硬闸门 → cases/"]
    I --> S["solve<br/>solve_J_*.py → J_result.json"]
    S --> M["mc<br/>mc_2d.py / mc_fss.py → T_C / T_N + Binder U4"]
    M --> R["render / report<br/>PNG + workflow.md"]
    D --> K["kconv<br/>k_convergence.py：只换 KPOINTS"]
    K --> I
    B -.->|不需要 VASP 的旁路| M
    classDef ext fill:#fde,stroke:#c39;
    class V ext;
```"""


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def headings(md: str, level: int = 2):
    """取 `level` 级标题（跳过 ``` 围栏内的内容——README 的目录树注释里满是 `#`）。"""
    out = []
    fence = False
    for line in md.splitlines():
        s = line.strip()
        if s.startswith("```") or s.startswith("~~~"):
            fence = not fence
            continue
        if fence:
            continue
        m = re.match(r"^(#{%d})\s+(.*?)\s*$" % level, line)
        if m:
            out.append(m.group(2))
    return out


def toc(entries):
    return "\n".join("- [[#%s]]" % e for e in entries)


def build(readme: str) -> str:
    """README 文本 → Obsidian 版全文（确定性：同一输入必得同一输出）。"""
    h2 = headings(readme, 2)
    head = "\n".join([
        "---",
        "title: 2d-magnet-pipeline — README（Obsidian 阅读版）",
        'aliases: ["2D 磁性格点流水线 README", "2d-magnet README"]',
        "tags: [materials-science, magnetism, monte-carlo, dft, vasp, python]",
        "source: README.md",
        "source_sha256: %s" % sha256_text(readme),
        "generated_by: scripts/readme_to_obsidian.py",
        "---",
        "",
        "> [!info] 这份文件是什么",
        "> 它是 [`README.md`](README.md) 的**阅读版备份**，由 `python scripts/readme_to_obsidian.py` 生成，"
        "正文与 README **逐字一致**（只在前面加了 properties、内链目录、管线图和这段说明）。",
        "> README 改过之后请重跑生成器；`python scripts/readme_to_obsidian.py --check` 是同步闸门"
        "（smoke §29 会跑它）。",
        "",
        "> [!tip] 在 Obsidian 里怎么用",
        "> 把 `2d_magnet_pipeline` **整个目录当 vault 打开**（Open folder as vault）——这样 "
        "`KNOWN_ISSUES.md` / `WORKFLOW.md` / `dsh/2d-magnet/SKILL.md` / `examples/` 这些相对链接都能跳。",
        "> 下面的 mermaid 图与 callout 都是 Obsidian 原生渲染，不需要插件；要交汇报可以用 "
        "Obsidian 的 Export to PDF。",
        "",
        "## 目录",
        "",
        toc(h2),
        "",
        "## 管线一览（主链）",
        "",
        PIPELINE_MERMAID,
        "",
        "---",
        "",
    ])
    return head + readme.rstrip("\n") + "\n"


def read(path: str) -> str:
    with open(path, encoding="utf-8") as f:
        return f.read()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="README.md → README.obsidian.md（Obsidian 阅读版备份）")
    ap.add_argument("--src", default=SRC)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--check", action="store_true", help="只校验同步（不同步则 exit 1）")
    ap.add_argument("--stdout", action="store_true", help="打到标准输出，不写文件")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)

    if a.selftest:
        return _selftest()

    if not os.path.isfile(a.src):
        print("[error] 找不到源文件: %s" % a.src)
        return 2
    text = read(a.src)
    want = build(text)

    if a.stdout:
        sys.stdout.write(want)
        return 0

    if a.check:
        if not os.path.isfile(a.out):
            print("[FAIL] 备份不存在: %s（跑 python scripts/readme_to_obsidian.py）" % a.out)
            return 1
        have = read(a.out)
        if have != want:
            print("[FAIL] 备份与 README.md 不同步: %s" % os.path.basename(a.out))
            print("       README 的 sha256 = %s" % sha256_text(text))
            m = re.search(r"source_sha256: ([0-9a-f]{64})", have)
            print("       备份里记录的     = %s" % (m.group(1) if m else "(缺)"))
            print("       修复: python scripts/readme_to_obsidian.py")
            return 1
        print("[OK] 备份与 README.md 同步（%d 字节, sha256 %s）"
              % (len(have.encode("utf-8")), sha256_text(text)[:12]))
        return 0

    with open(a.out, "w", encoding="utf-8", newline="\n") as f:
        f.write(want)
    print("[OK] 已生成 %s（%d 字节, README sha256 %s）"
          % (a.out, len(want.encode("utf-8")), sha256_text(text)[:12]))
    return 0


def _selftest() -> int:
    fails = []

    def ck(name, cond):
        print("  %s %s" % ("PASS" if cond else "FAIL", name))
        if not cond:
            fails.append(name)

    sample = ("# T\n\n## 一 A\n\ntext\n\n```\n## 不该进目录\n```\n\n"
              "## 二 B（括号）\n\nmore\n\n### 三级不进目录\n")
    h = headings(sample, 2)
    ck("headings: 只取二级标题", h == ["一 A", "二 B（括号）"])
    ck("headings: 跳过围栏内的 ## 注释", "不该进目录" not in h)
    ck("toc: 生成 Obsidian 内链", toc(h).splitlines()[0] == "- [[#一 A]]")

    out = build(sample)
    ck("build: frontmatter 键齐全", out.startswith("---\n") and all(("%s:" % k) in out for k in FRONT_KEYS))
    ck("build: 正文逐字保留", sample.rstrip("\n") in out)
    ck("build: 记录 README 的 sha256", sha256_text(sample) in out)
    ck("build: 带 mermaid 管线图", "```mermaid" in out)
    ck("build: 确定性（两次构建一致）", build(sample) == build(sample))
    ck("build: 不同输入 ⇒ 不同 sha256", sha256_text(sample) != sha256_text(sample + "\n"))

    # 真实 README（若在仓库里就顺带查一遍，保证脚本对真文件不炸）
    if os.path.isfile(SRC):
        real = read(SRC)
        ro = build(real)
        ck("real: README 的二级标题都进了目录（%d 条）" % len(headings(real, 2)),
           all(("[[#%s]]" % x) in ro for x in headings(real, 2)))
        ck("real: 正文没被改写", real.rstrip("\n") in ro)

    print("\n--selftest: %s" % ("ALL PASS" if not fails else "%d FAIL: %s" % (len(fails), fails)))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
