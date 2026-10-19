"""Build a deterministic release archive from the working tree (tracked + untracked-but-not-ignored files).

用法::

    python scripts/make_release.py . dist

产物: ``dist/2d-magnet-pipeline-<version>.tar.gz`` / ``.zip`` / ``SHA256SUMS``。

- 名单来自 ``git ls-files --cached --others --exclude-standard``（**已跟踪 + 未跟踪但没被 .gitignore 忽略**的文件；
  未跟踪的 ``*.pptx`` / ``dist/`` / ``_tmp_*`` / ``__pycache__`` 都被 .gitignore 挡住，不会进包）
  —— 这样「先打包验证、后提交」的流程也能出完整包，提交之后同一棵树打出来的名单不变
- 内容取**工作树当前内容**（因此包含尚未提交的改动 —— 想发「某个提交的精确内容」，先提交/打 tag 再跑）
- 元数据固定（uid/gid=0、mtime 取 HEAD 提交时间、gzip 不写时间戳），同一棵树重复构建**字节一致**
- 打包后用 ``python tests/check_sdist.py dist/*.tar.gz``（发布闸门）与
  ``python tests/smoke_test.py``（解包后跑一次）验证
"""

import hashlib
import io
import os
import subprocess
import sys
import tarfile
import zipfile


def version_of(root):
    for line in io.open(os.path.join(root, "pyproject.toml"), encoding="utf-8"):
        if line.strip().startswith("version"):
            return line.split("=", 1)[1].strip().strip('"\'')
    raise SystemExit("[error] pyproject.toml 里没有 version")


def git(root, *args):
    return subprocess.run(["git", "-C", root] + list(args),
                          capture_output=True, check=True).stdout


def commit_time(root):
    return int(git(root, "log", "-1", "--format=%ct").decode().strip())


def release_files(root):
    raw = git(root, "ls-files", "-z", "--cached", "--others", "--exclude-standard")
    names = [n.decode("utf-8") for n in raw.split(b"\0") if n]
    # 已 git rm 但还没提交的文件仍在索引里; 工作树里没有就跳过, 否则打不出包。
    kept = []
    skipped = []
    for n in names:
        if os.path.isfile(os.path.join(root, n)):
            kept.append(n)
        else:
            skipped.append(n)
    if skipped:
        print("[warn] 跳过索引里有、工作树里没有的文件 %d 个: %s"
              % (len(skipped), ", ".join(skipped[:5])))
    return sorted(kept)


def build_tar(root, out, prefix, mtime, names):
    # 用 mtime=0 的 gzip 包装: 默认 GzipFile 会把「当前时间」写进 gzip 头,
    # 那样同一棵树在不同秒构建会得到不同字节, 与「可复现构建」的承诺矛盾。
    import gzip as _gz
    with open(out, "wb") as raw:
        with _gz.GzipFile(fileobj=raw, mode="wb", compresslevel=9, mtime=0) as gz:
            with tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as tf:
                for rel in names:
                    path = os.path.join(root, rel)
                    data = open(path, "rb").read()
                    info = tarfile.TarInfo(prefix + rel.replace(os.sep, "/"))
                    info.size = len(data)
                    info.mtime = mtime
                    info.mode = 0o644
                    info.uid = info.gid = 0
                    info.uname = info.gname = "root"
                    tf.addfile(info, io.BytesIO(data))


def build_zip(root, out, prefix, mtime, names):
    import time as _t
    lt = _t.gmtime(mtime)[:6]  # zip 的 date_time 最小年份是 1980
    lt = (max(lt[0], 1980),) + lt[1:]
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for rel in names:
            data = open(os.path.join(root, rel), "rb").read()
            info = zipfile.ZipInfo(prefix + rel.replace(os.sep, "/"), date_time=lt)
            info.external_attr = (0o644 & 0xFFFF) << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            zf.writestr(info, data)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    root, outdir = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])
    os.makedirs(outdir, exist_ok=True)
    ver = version_of(root)
    prefix = "2d-magnet-pipeline-%s/" % ver
    mtime = commit_time(root)
    names = release_files(root)
    tar_path = os.path.join(outdir, "2d-magnet-pipeline-%s.tar.gz" % ver)
    zip_path = os.path.join(outdir, "2d-magnet-pipeline-%s.zip" % ver)
    build_tar(root, tar_path, prefix, mtime, names)
    build_zip(root, zip_path, prefix, mtime, names)
    lines = ["# sha256  (2d-magnet-pipeline %s, %d 个文件, 内容取自工作树)" % (ver, len(names))]
    for p in (tar_path, zip_path):
        lines.append("%s  %s" % (sha256(p), os.path.basename(p)))
    sums = os.path.join(outdir, "SHA256SUMS")
    io.open(sums, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    print("[OK] %s（%d 字节）" % (os.path.basename(tar_path), os.path.getsize(tar_path)))
    print("[OK] %s（%d 字节）" % (os.path.basename(zip_path), os.path.getsize(zip_path)))
    print("[OK] SHA256SUMS")
    print("\n".join(lines[1:]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
