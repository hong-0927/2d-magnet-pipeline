"""Build the uploadable agent-skill bundle.

与 ``make_release.py``（仓库包：整棵源码树原样）不同，本脚本产出的是**能直接上传/发布的 skill 包**：

    dist/2d-magnet-skill-<version>.zip
    dist/2d-magnet-skill-<version>.zip.sha256

包结构（解压后这一层目录就是 SKILL.md 里说的 ``<repo>``）::

    2d-magnet-skill-<version>/
        SKILL.md                  ← skill 入口（frontmatter 保留；路径措辞改成「包根」）
        bin/2d-magnet.cmd         ← 封装脚本（ROOT 改成 %~dp0..）
        run_2d_magnet.py
        scripts/  knowledge/  tests/  examples/  templates/
        capabilities.json  requirements.txt  pyproject.toml  install.*  README*.md ...

三条硬约束：

1. **SKILL.md 必须在包根**：多数 harness 的技能目录是 ``<skills>/<name>/SKILL.md``，
   只有入口在包根，解压即被识别。
2. **包里必须自带流水线**：SKILL.md 驱动的是 ``run_2d_magnet.py`` + ``scripts/``；
   只拷 ``skill/2d-magnet/`` 会得到一个指向不存在路径的说明文件。仓库里保留 ``skill/`` 原样，
   包内改成根布局，因此两处的路径措辞不同（本脚本负责改写，不手改仓库那份）。
3. **确定性**：文件清单来自 git（已跟踪 + 未忽略），时间戳取 HEAD 提交时间，
   同一棵树重复构建字节一致。

用法::

    python scripts/make_skill_bundle.py . dist
"""

import io
import os
import sys
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from make_release import commit_time, release_files, sha256, version_of  # noqa: E402

# 这些不是 skill 包该带的东西：仓库专属配置、CI、仓库内的 skill 目录（包内改成根布局）、构建产物。
SKIP_TOP = (".github/", "dist/", "skill/", "__pycache__/", ".git/")
SKIP_FILES = {".gitignore", ".gitattributes", "MANIFEST.in"}

# SKILL.md 里必须改写到位、否则解压后路径全断的措辞。
SKILL_REWRITES = (
    ("（本 SKILL.md 所在目录的\n上两级）", "（即本 SKILL.md 所在的包根目录）"),
    ("`<repo>\\skill\\bin\\2d-magnet.cmd`", "`<repo>\\bin\\2d-magnet.cmd`"),
    ("（它按 `%~dp0..\\..` 自己定位", "（它按 `%~dp0..` 自己定位"),
)
CMD_REWRITES = (
    ('set "ROOT=%~dp0..\\.."', 'set "ROOT=%~dp0.."'),
    ("Pipeline root = two levels up from this file (..\\..\\ = <repo>).",
     "Pipeline root = one level up from this file (..\\ = <repo>, skill bundle root)."),
)


def rewrite(text, rules, label):
    for old, new in rules:
        if old not in text:
            raise SystemExit("[error] %s 里找不到要改写的片段: %r" % (label, old))
        text = text.replace(old, new)
    return text


def wants(rel):
    if any(rel.startswith(p) for p in SKIP_TOP):
        return False
    return rel not in SKIP_FILES


def main():
    root = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
    outdir = os.path.abspath(sys.argv[2] if len(sys.argv) > 2 else os.path.join(root, "dist"))
    ver = version_of(root)
    prefix = "2d-magnet-skill-%s/" % ver

    names = [n for n in release_files(root) if wants(n.replace(os.sep, "/"))]

    skill_src = os.path.join(root, "skill", "2d-magnet", "SKILL.md")
    cmd_src = os.path.join(root, "skill", "bin", "2d-magnet.cmd")
    with io.open(skill_src, encoding="utf-8") as fh:
        skill = rewrite(fh.read(), SKILL_REWRITES, "SKILL.md")
    with io.open(cmd_src, encoding="utf-8") as fh:
        cmd = rewrite(fh.read(), CMD_REWRITES, "2d-magnet.cmd")
    if "\nname: 2d-magnet" not in skill:
        raise SystemExit("[error] SKILL.md 缺 frontmatter name: 2d-magnet")

    import time as _t
    lt = _t.gmtime(commit_time(root))[:6]
    lt = (max(lt[0], 1980),) + lt[1:]

    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, "2d-magnet-skill-%s.zip" % ver)
    entries = [(prefix + "SKILL.md", skill.encode("utf-8")),
               (prefix + "bin/2d-magnet.cmd", cmd.encode("utf-8"))]
    for rel in names:
        with open(os.path.join(root, rel), "rb") as fh:
            entries.append((prefix + rel.replace(os.sep, "/"), fh.read()))

    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for arcname, data in entries:
            info = zipfile.ZipInfo(arcname, date_time=lt)
            info.external_attr = (0o644 & 0xFFFF) << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            zf.writestr(info, data)

    digest = sha256(out)
    sidecar = out + ".sha256"
    io.open(sidecar, "w", encoding="utf-8", newline="\n").write(
        "# sha256 (2d-magnet skill %s, %d entries)\n%s  %s\n"
        % (ver, len(entries), digest, os.path.basename(out)))

    print("[OK] %s（%d 字节，%d 项）" % (os.path.basename(out), os.path.getsize(out), len(entries)))
    print("[OK] %s" % os.path.basename(sidecar))
    print("     SKILL.md -> %sSKILL.md" % prefix)
    return 0


if __name__ == "__main__":
    sys.exit(main())
