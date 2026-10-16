#!/usr/bin/env python
"""dft_handoff.py —— P4.2: 「本机 <-> 外部 DFT 机时」的交接包与提交脚本。

## 为什么需要它

本机**没有 VASP、没有机时**(硬约束)。P4.2 的物理层只能这样走:

    本机 build/dft  ──交接包──>  集群跑 VASP  ──回收──>  本机 ingest/solve/MC

P4.2 之前, 中间两步是「手工 tar + 手写作业脚本 + 手工 cp 回来」——这是整条链上
最不可复现的一段。本脚本把它变成三个可重放的命令:

    python scripts/dft_handoff.py pack   --root <proj> --out <handoff_dir>
    # 把 handoff_dir 上传到集群, 在里面 bash submit.sh
    python scripts/dft_handoff.py ingest --root <proj> --from <放回本机的目录>

## pack 产出什么

    <handoff_dir>/
      README_HANDOFF.md      <- 交接说明 (人读: 要什么、怎么交、怎么自查)
      MANIFEST.json          <- 每个文件的 sha256 + 尺寸 (集群侧可核对没传坏)
      POTCAR.spec            <- 元素顺序 + 建议 PAW (**POTCAR 本身有 license, 不随包**)
      CHECKLIST.md           <- 提交前/回收前的逐条勾选
      submit.sh              <- 自动识别 SLURM/PBS/SGE/LSF/裸机 的提交脚本
      jobscript.template     <- 作业脚本模板 (占位符由 submit.sh 填)
      runs/<order>/          <- INCAR / KPOINTS / POSCAR (VASP 输入, 原样打包)

## 关键规矩

  * **不打包 POTCAR** (license 不允许随包分发) —— 集群侧自己 `cat` 出来。
  * **不代提交** —— 本脚本只生成脚本, 提交动作由人在集群上做。
  * 打包时对 runs/ 做**完整性预检**: 缺 INCAR/KPOINTS/POSCAR 的序直接报错,
    免得人跑到集群才发现少文件。

## 用法

  python dft_handoff.py pack --root <proj> --out <dir> [--scheduler auto|slurm|pbs|sge|lsf|none]
  python dft_handoff.py ingest --root <proj> --from <dir> [--dry-run] [--force]
  python dft_handoff.py --selftest
"""
import argparse
import hashlib
import json
import os
import shutil
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ingest_dft as IG  # noqa: E402

VERSION = "1.0"

NEEDED = ("INCAR", "KPOINTS", "POSCAR")

# 调度器 -> (作业脚本模板, 提交命令)。模板里的 %(...)s 占位符由 submit.sh 用 envsubst
# 风格的 sed 填 —— 不用 heredoc, 免得踩不同 shell 的兼容性坑。
SCHEDULERS = {
    "slurm": {
        "shebang": "#!/bin/bash",
        "header": [
            "#SBATCH --job-name=%(name)s",
            "#SBATCH --nodes=1",
            "#SBATCH --ntasks=%(ncpu)s",
            "#SBATCH --time=%(walltime)s",
            "#SBATCH --output=slurm-%%j.out",
            "#SBATCH --error=slurm-%%j.err",
            "#SBATCH --partition=%(partition)s",
        ],
        "submit": "sbatch jobscript.sh",
        "note": "SLURM。walltime 形如 24:00:00; 若集群用 --cpus-per-task 代替 --ntasks, "
                "改 submit.sh 里的 HEADER 数组。",
    },
    "pbs": {
        "shebang": "#!/bin/bash",
        "header": [
            "#PBS -N %(name)s",
            "#PBS -l nodes=1:ppn=%(ncpu)s",
            "#PBS -l walltime=%(walltime)s",
            "#PBS -j oe",
            "#PBS -q %(partition)s",
        ],
        "submit": "qsub jobscript.sh",
        "note": "PBS/Torque。若集群是 OpenPBS 且 ppn 不被识别, 改成 "
                "`#PBS -l select=1:ncpus=%(ncpu)s`。",
    },
    "sge": {
        "shebang": "#!/bin/bash",
        "header": [
            "#$ -N %(name)s",
            "#$ -pe smp %(ncpu)s",
            "#$ -l h_rt=%(walltime)s",
            "#$ -cwd",
            "#$ -j y",
            "#$ -q %(partition)s",
        ],
        "submit": "qsub jobscript.sh",
        "note": "SGE。h_rt 用 H:M:S。",
    },
    "lsf": {
        "shebang": "#!/bin/bash",
        "header": [
            "#BSUB -J %(name)s",
            "#BSUB -n %(ncpu)s",
            "#BSUB -W %(walltime_hhmm)s",
            "#BSUB -o lsf-%%J.out",
            "#BSUB -q %(partition)s",
        ],
        "submit": "bsub < jobscript.sh",
        "note": "LSF。",
    },
    "none": {
        "shebang": "#!/bin/bash",
        "header": [],
        "submit": "bash jobscript.sh",
        "note": "裸机/交互式 (无调度器)。逐个序串行跑, 慢但最不容易出配置问题。",
    },
}
DEFAULT_WALLTIME = {"slurm": "24:00:00", "pbs": "24:00:00", "sge": "24:00:00",
                    "lsf": "24:00", "none": ""}
DEFAULT_PARTITION = {"slurm": "compute", "pbs": "workq", "sge": "all.q",
                     "lsf": "normal", "none": ""}


# ---------------------------------------------------------------------------
def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def read_meta(root):
    """读 02_supercell/cases/magnetic_orders.json -> (orders, meta)。"""
    p = os.path.join(root, "02_supercell", "cases", "magnetic_orders.json")
    if not os.path.isfile(p):
        raise SystemExit(f"缺 {p} —— 先跑 `run_2d_magnet.py build` (或全流水线)")
    with open(p, encoding="utf-8-sig") as f:
        meta = json.load(f)
    orders = list(meta.get("selected_orders") or [])
    if not orders:
        raise SystemExit(f"{p} 没有 selected_orders —— 无法确定要提交哪些磁态")
    return orders, meta


def preflight_runs(root, orders):
    """打包前的完整性预检: 每个序必须有 INCAR/KPOINTS/POSCAR。缺一件就报错。

    这一步存在的原因: 人跑到集群才发现少文件, 一个作业排 6 小时。
    """
    runs = os.path.join(root, "runs")
    problems = []
    if not os.path.isdir(runs):
        raise SystemExit(f"缺 {runs} —— 先跑 `run_2d_magnet.py dft` 导出 VASP 输入")
    for o in orders:
        d = os.path.join(runs, o)
        if not os.path.isdir(d):
            problems.append(f"  [缺目录] runs/{o}/")
            continue
        for nm in NEEDED:
            if not os.path.isfile(os.path.join(d, nm)):
                problems.append(f"  [缺文件] runs/{o}/{nm}")
    if problems:
        raise SystemExit("runs/ 不完整, 拒绝打包 (免得人跑到集群才发现):\n"
                         + "\n".join(problems)
                         + "\n  修法: 重跑 `python run_2d_magnet.py dft --outdir <proj>`")
    return runs


# ---------------------------------------------------------------------------
# POTCAR.spec 校验 (打包时把"元素顺序"再核一遍)
# ---------------------------------------------------------------------------
def potcar_species(root, order):
    """从 runs/<order>/POSCAR 取元素顺序 (POTCAR 必须按这个顺序 cat)。"""
    p = os.path.join(root, "runs", order, "POSCAR")
    lines = [ln for ln in open(p, encoding="utf-8").read().splitlines() if ln.strip()]
    return lines[5].split(), [int(x) for x in lines[6].split()]


# ---------------------------------------------------------------------------
# 作业脚本
# ---------------------------------------------------------------------------
JOBSCRIPT = """%(shebang)s
%(header)s
#
# P4.2 kagome DFT 作业脚本 —— 由 dft_handoff.py 生成, 请按集群实际情况改。
# 用法 (在解包后的目录里):  bash submit.sh
#
set -euo pipefail

# ---- 集群相关: 这三行**必须**按你的环境改 ------------------------------------
VASP_BIN="${VASP_BIN:-vasp_std}"     # 可执行名: vasp_std / vasp_ncl / mpirun -np N vasp_ncl
POTCAR_DIR="${POTCAR_DIR:-.}"        # 你的 POTCAR 库存根 (下面会 cat $POTCAR_DIR/<El>/POTCAR)
NCPU="${NCPU:-%(ncpu)s}"

ORDERS="%(orders)s"

echo "=== P4.2 kagome DFT: $(date -u +%%FT%%TZ) ==="
echo "VASP_BIN=$VASP_BIN  POTCAR_DIR=$POTCAR_DIR  NCPU=$NCPU"

for o in $ORDERS; do
  echo
  echo "--- [$o] 组装 POTCAR ---"
  if [ ! -d "$o" ]; then echo "!! 缺目录 $o"; exit 1; fi
  cp KPOINTS "$o/KPOINTS"
  # POTCAR 顺序必须与 POSCAR 的元素块顺序一致 (见 POTCAR.spec / README_HANDOFF.md)
  : > "$o/POTCAR"
  for el in %(species)s; do
    if [ ! -f "$POTCAR_DIR/$el/POTCAR" ]; then
      echo "!! 找不到 $POTCAR_DIR/$el/POTCAR —— 请设 POTCAR_DIR 或确认 POTCAR 库布局"
      exit 2
    fi
    cat "$POTCAR_DIR/$el/POTCAR" >> "$o/POTCAR"
  done
  echo "    POTCAR 大小 $(wc -c < "$o/POTCAR") 字节 (应为数个 10^5, 太小说明 cat 错了)"

  echo "--- [$o] 跑 VASP ---"
  cd "$o"
  cp ../INCAR ./INCAR
  if [ "$VASP_BIN" = "vasp_ncl" ] || [ "$VASP_BIN" = "vasp_gam" ]; then
    :   # 非共线必须用 vasp_ncl: 见 README_HANDOFF.md 的自旋模式一节
  fi
  # shellcheck disable=SC2086
  $VASP_BIN > vasp.out 2>&1 || { echo "!! [$o] VASP 退出码非零, 看 $o/vasp.out"; cd ..; exit 3; }

  # 就地留证: 记录这次用的二进制/输入指纹, 供 ingest 的 provenance 用
  {
    echo "{"
    echo "  \\"order\\": \\"$o\\","
    echo "  \\"vasp_bin\\": \\"$VASP_BIN\\","
    echo "  \\"host\\": \\"$(hostname)\\","
    echo "  \\"finished_utc\\": \\"$(date -u +%%FT%%TZ)\\","
    echo "  \\"incar_sha256\\": \\"$(sha256sum INCAR 2>/dev/null | cut -d' ' -f1 || echo '')\\","
    echo "  \\"outcar_bytes\\": $(wc -c < OUTCAR 2>/dev/null || echo 0)"
    echo "}"
  } > run_meta.json
  cd ..
  echo "--- [$o] 完成 (OUTCAR $(wc -c < "$o/OUTCAR" 2>/dev/null || echo 0) 字节)"
done

echo
echo "=== 全部结束: $(date -u +%%FT%%TZ) ==="
echo "回收前请先自查:"
echo "  1) 每个序有 OUTCAR 且不是空壳"
echo "  2) OUTCAR 里有 'General timing and accounting'"
echo "  3) 有 'aborting loop because EDIFF is reached' (SCF 收敛)"
echo "  4) 本机侧会再跑 ingest_dft.py 的 10 条硬闸门, 不合格会被拒收"
"""

SUBMIT_SH = """#!/bin/bash
# submit.sh —— 由 dft_handoff.py 生成。在**这个目录**里执行:  bash submit.sh
# 作用: 把 jobscript.template 里的占位符填成实际值 -> jobscript.sh, 然后按调度器提交。
set -euo pipefail
cd "$(dirname "$0")"

SCHEDULER="${SCHEDULER:-%(scheduler)s}"
NCPU="${NCPU:-%(ncpu)s}"
WALLTIME="${WALLTIME:-%(walltime)s}"
PARTITION="${PARTITION:-%(partition)s}"

case "$SCHEDULER" in
  slurm) SHEBANG="#!/bin/bash"; HEADER=(%(slurm_header)s); SUBMIT="sbatch jobscript.sh" ;;
  pbs)   SHEBANG="#!/bin/bash"; HEADER=(%(pbs_header)s); SUBMIT="qsub jobscript.sh" ;;
  sge)   SHEBANG="#!/bin/bash"; HEADER=(%(sge_header)s); SUBMIT="qsub jobscript.sh" ;;
  lsf)   SHEBANG="#!/bin/bash"; HEADER=(%(lsf_header)s); SUBMIT="bsub < jobscript.sh" ;;
  none)  SHEBANG="#!/bin/bash"; HEADER=(); SUBMIT="bash jobscript.sh" ;;
  *) echo "!! 未知 SCHEDULER=$SCHEDULER (可选 slurm/pbs/sge/lsf/none)"; exit 2 ;;
esac

{
  echo "$SHEBANG"
  for h in "${HEADER[@]:-}"; do [ -n "$h" ] && echo "$h"; done
  # 模板正文 (跳过模板自己的 shebang 与 header 区)
  sed -n '/^# ---- 正文开始/,$p' jobscript.template
} > jobscript.sh
chmod +x jobscript.sh

echo "已生成 jobscript.sh (SCHEDULER=$SCHEDULER NCPU=$NCPU WALLTIME=$WALLTIME PARTITION=$PARTITION)"
echo "提交命令: $SUBMIT"
echo
echo "自查: 少 1 个序 => 本机 ingest 会拒收 (设计矩阵需要四个独立能量)。"
echo "跑完把整个目录 (含各 <order>/OUTCAR) 拷回本机, 然后:"
echo "  python scripts/dft_handoff.py ingest --root <proj> --from <这个目录>"
"""

# 模板正文的分隔标记 (submit.sh 会 sed 出这一段)
BODY_MARK = "# ---- 正文开始"


def render_jobscript(orders, species, scheduler, ncpu, walltime, partition):
    sch = SCHEDULERS[scheduler]
    fmt = dict(name="kagome_dft", ncpu=ncpu, walltime=walltime, partition=partition,
               walltime_hhmm=walltime.replace(":", "")[:4] or "2400",
               orders=" ".join(orders), species=" ".join(species),
               shebang=sch["shebang"])
    head = "\n".join(h % fmt for h in sch["header"])
    body = JOBSCRIPT % dict(fmt, header=BODY_MARK)
    return body.replace(BODY_MARK + "\n", BODY_MARK + "\n"), head


# ---------------------------------------------------------------------------
# README / CHECKLIST
# ---------------------------------------------------------------------------
def render_readme(root, out, orders, species, counts, meta, scheduler, want_noncollinear):
    nc = meta.get("n_mag")
    S2 = meta.get("S2")
    L = []
    L.append("# P4.2 kagome (Fe3Sn 型) DFT 交接包")
    L.append("")
    L.append(f"- 生成: {datetime.now(timezone.utc).isoformat(timespec='seconds')} (UTC) "
             f"by `dft_handoff.py` v{VERSION}")
    L.append(f"- 源项目: `{os.path.abspath(root)}`")
    L.append(f"- 调度器模板: `{scheduler}`")
    L.append(f"- 磁态 ({len(orders)} 个): `{orders}`")
    L.append(f"- 元素顺序 (POTCAR 必须按这个顺序 cat): `{species}` = `{counts}`")
    L.append(f"- 每胞磁原子数 `n_mag` = {nc}; `S^2` = {S2}")
    L.append("")
    L.append("## 这个包要什么")
    L.append("")
    L.append("**只要一件东西: 每个磁态一个 `OUTCAR`**, 且满足下面全部条件。")
    L.append("")
    L.append("| # | 条件 | 怎么自查 |")
    L.append("|---|---|---|")
    L.append("| 1 | **静态**计算: 每个 OUTCAR 里 "
             "`free  energy   TOTEN` **只出现一次** | `grep -c 'free  energy   TOTEN' OUTCAR` = 1 |")
    L.append("| 2 | SCF **收敛**: 有 `aborting loop because EDIFF is reached` | `grep 'aborting loop' OUTCAR` |")
    L.append("| 3 | 跑完: 有 `General timing and accounting` | `grep -c 'General timing' OUTCAR` |")
    L.append("| 4 | 四个序**同一个结构、同一套数值设置** "
             "(ENCUT/K 点/POTCAR 版本/自旋模式全同) | 看 `MANIFEST.json` + 本文件的设置节 |")
    L.append("| 5 | OUTCAR 不能是空壳 (至少几百 KB) | `ls -l */OUTCAR` |")
    L.append("")
    L.append("⚠ **不要**只交 `OSZICAR` / `vasp.out` / 能量汇总表 —— 收回通道只认 OUTCAR, "
             "因为要核对收敛标志与自旋模式。")
    L.append("")
    L.append("## 自旋模式 (必须与 INCAR 一致)")
    L.append("")
    if want_noncollinear:
        L.append("**非共线** (`LNONCOLLINEAR = .TRUE.`)。kagome 的 120° 与 √3×√3 序是"
                 "**面内非共线**, 共线设置算不出来。")
        L.append("")
        L.append("- 用 **`vasp_ncl`** 而不是 `vasp_std` (VASP 对 LNONCOLLINEAR 有专门的二进制)。")
        L.append("- 每个 INCAR 的 `MAGMOM` 有 **3N 个分量** (N = 原子数); "
                 "共线设置下只会读前 N 个 —— 这是最隐蔽的一类错误。")
        L.append("- 本包是否开 SOC 由原项目的 `--soc` 决定: 看 INCAR 里有没有 `LSORBIT`。")
    else:
        L.append("共线 (`ISPIN = 2`)。")
    L.append("")
    L.append("## POTCAR")
    L.append("")
    L.append("**包里没有 POTCAR** (license 不允许随包分发)。请在集群上自己拼:")
    L.append("")
    L.append("```bash")
    L.append("# 顺序必须严格按 POSCAR 的元素块顺序")
    L.append(f"cat {' '.join('$POTCAR_DIR/' + s + '/POTCAR' for s in species)} > POTCAR")
    L.append("# 或者用仓库给出的建议 (见 POTCAR.spec)")
    L.append("```")
    L.append("")
    L.append("推荐 PAW 版本见 `POTCAR.spec`。**四个序必须用同一套 POTCAR** —— 换版本会让能量"
             "整体平移, 而反解 J 吃的是序间能量**差**。")
    L.append("")
    L.append("## 怎么跑")
    L.append("")
    L.append("```bash")
    L.append("tar xzf <这个包>.tar.gz && cd <解包目录>")
    L.append("export POTCAR_DIR=/path/to/potpaw_PBE      # 你的 POTCAR 库")
    L.append("export VASP_BIN=vasp_ncl                   # 非共线用 ncl")
    L.append("export NCPU=64")
    L.append("bash submit.sh                             # 生成并提交 jobscript.sh")
    L.append("# 或者不想用调度器:  SCHEDULER=none bash submit.sh")
    L.append("```")
    L.append("")
    L.append("## 交回来什么")
    L.append("")
    L.append("把**整个目录**(含每个 `<order>/OUTCAR`)打包寄回, 或者只寄 4 个 OUTCAR:"
             "保持 `<order>/OUTCAR` 的相对路径即可。")
    L.append("")
    L.append("```bash")
    L.append("# 本机侧 (唯一的收回通道, 会跑 10 条硬闸门)")
    L.append("python scripts/dft_handoff.py ingest --root <proj> --from <放回来的目录>")
    L.append("# 只想先看看合不合格, 不写入:")
    L.append("python scripts/dft_handoff.py ingest --root <proj> --from <目录> --dry-run")
    L.append("```")
    L.append("")
    L.append("收回通道会**拒收**这些情况 (每一条都在本机上跑过 selftest):")
    L.append("")
    L.append("- 缺任何一个磁态 / 空壳 OUTCAR / 没跑完 / SCF 未收敛")
    L.append("- 有多个离子步 (= 结构被松弛过, 与 `cases/` 的系数不对应)")
    L.append("- 能量口径混用 / 读不出 `free energy TOTEN`")
    L.append("- 四个序的结构指纹或自旋模式/ENCUT/K 点不一致")
    L.append("- 某个序与 `fm` 的能量差 ≈ 0 (= 该磁序没生效, MAGMOM 没进 INCAR)")
    L.append("- OUTCAR 与 `cases/` 里已有的**逐位相同** (= 旧结果充新, 需要 `--force` 才放行)")
    L.append("")
    L.append("## 这个包**不**解决什么 (P4.2 的边界)")
    L.append("")
    L.append("- 它不产生 J, 不产生 `T_N`/`T_C`。那些要等 OUTCAR 回收 + `solve` 反解 + MC。")
    L.append("- `J1` 是 kagome 三角形内两条**不等价**最近邻键 (`J1` 与 `J1'`) 的**平均**: "
             "3 壳模型的固有近似。要分开它们需要降对称超胞 (P4.2 未做)。")
    L.append("- 本模型**没有层间耦合与磁各向异性** ⇒ 二维各向同性 MC 的临界温度**不能**"
             "直接当 `T_N`/`T_C` 与文献比 (Mermin-Wagner: 二维各向同性 Heisenberg 在有限"
             "温度没有长程序, MC 的峰是有限尺寸伪临界点)。量级对照见 runbook。")
    return "\n".join(L) + "\n"


def render_checklist(orders, scheduler):
    L = ["# 提交 / 回收 逐条勾选", ""]
    L.append("## 提交前 (本机已自动预检的部分打 x)")
    L.append("")
    L.append("- [x] `runs/` 四个序的 INCAR/KPOINTS/POSCAR 齐全 (pack 时已硬检)")
    L.append("- [ ] `POTCAR` 已在集群按 POSCAR 元素顺序拼好, 四个序**同一套**")
    L.append("- [ ] 调度器模板按集群改好 (`submit.sh` 的 HEADER / `POTCAR_DIR` / `VASP_BIN`)")
    L.append("- [ ] 非共线体系用的是 `vasp_ncl`")
    L.append("- [ ] 四个序用**同一个** ENCUT / K 点 / POTCAR 版本")
    L.append("- [ ] 作业时长够 (静态单点: 108 磁原子 3x3x1 超胞, 见 runbook §9 机时估算)")
    L.append("")
    L.append("## 回收前 (每一条都要过, 否则本机 ingest 会拒收)")
    L.append("")
    for o in orders:
        L.append(f"- [ ] `{o}/OUTCAR` 存在且 >= 200 KB, "
                 f"`grep -c 'free  energy   TOTEN' {o}/OUTCAR` = 1")
    L.append("- [ ] 每个 OUTCAR 都有 `aborting loop because EDIFF is reached`")
    L.append("- [ ] 每个 OUTCAR 都有 `General timing and accounting`")
    L.append("- [ ] 没有任何一个序是重跑过的旧文件 (sha256 会与本机 cases/ 对比)")
    L.append("")
    L.append("## 回收命令")
    L.append("")
    L.append("```bash")
    L.append("python scripts/dft_handoff.py ingest --root <proj> --from <目录> --dry-run  # 先自查")
    L.append("python scripts/dft_handoff.py ingest --root <proj> --from <目录>            # 正式收回")
    L.append("python run_2d_magnet.py solve --outdir <proj>                               # 反解 J")
    L.append("```")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
def pack(root, out, scheduler="auto", ncpu=32, walltime=None, partition=None,
         verbose=True):
    """组装交接包。返回 dict 摘要。"""
    say = print if verbose else (lambda *a, **k: None)
    root = os.path.abspath(root)
    orders, meta = read_meta(root)
    runs = preflight_runs(root, orders)
    species, counts = potcar_species(root, orders[0])

    # 自旋模式: 由 INCAR 判定 (LNONCOLLINEAR 在不在)
    inc0 = os.path.join(runs, orders[0], "INCAR")
    inc_txt = open(inc0, encoding="utf-8").read()
    nl = [ln for ln in inc_txt.splitlines() if ln.strip() and not ln.strip().startswith("!")]
    want_nc = any(ln.split("=")[0].strip().upper() == "LNONCOLLINEAR" for ln in nl)
    has_soc = any(ln.split("=")[0].strip().upper() == "LSORBIT" for ln in nl)

    if scheduler == "auto":
        scheduler = "slurm"      # 集群最常见; README 里说明了怎么改
    if scheduler not in SCHEDULERS:
        raise SystemExit(f"未知 --scheduler {scheduler!r}; 可选 {sorted(SCHEDULERS)}")
    walltime = walltime or DEFAULT_WALLTIME[scheduler]
    partition = partition or DEFAULT_PARTITION[scheduler]

    if os.path.isdir(out):
        say(f"[pack] 目标已存在, 只覆盖文件: {out}")
    os.makedirs(out, exist_ok=True)

    # runs/<order>/ -> <out>/<order>/
    for o in orders:
        dst = os.path.join(out, o)
        os.makedirs(dst, exist_ok=True)
        for nm in NEEDED + ("run_meta.json",):
            src = os.path.join(runs, o, nm)
            if os.path.isfile(src):
                shutil.copy2(src, os.path.join(dst, nm))
        # INCAR 放到包根 (submit.sh 会 cp 进去), 每个序的原件也留着
    shutil.copy2(inc0, os.path.join(out, "INCAR"))
    shutil.copy2(os.path.join(runs, orders[0], "KPOINTS"), os.path.join(out, "KPOINTS"))

    spec = os.path.join(root, "03_dft", "POTCAR.spec")
    if os.path.isfile(spec):
        shutil.copy2(spec, os.path.join(out, "POTCAR.spec"))

    body, _head = render_jobscript(orders, species, scheduler, ncpu, walltime, partition)
    with open(os.path.join(out, "jobscript.template"), "w", encoding="utf-8", newline="\n") as f:
        f.write(body)

    def _q(s):
        return '"' + s.replace('"', '\\"') + '"'

    fmt = dict(name="kagome_dft", ncpu=ncpu, walltime=walltime, partition=partition,
               walltime_hhmm=walltime.replace(":", "")[:4] or "2400",
               orders=" ".join(orders), species=" ".join(species))

    def hdr(sch):
        return " ".join(_q(h % fmt) for h in SCHEDULERS[sch]["header"])

    submit = SUBMIT_SH % dict(scheduler=scheduler, ncpu=ncpu, walltime=walltime,
                              partition=partition, slurm_header=hdr("slurm"),
                              pbs_header=hdr("pbs"), sge_header=hdr("sge"),
                              lsf_header=hdr("lsf"))
    with open(os.path.join(out, "submit.sh"), "w", encoding="utf-8", newline="\n") as f:
        f.write(submit)

    with open(os.path.join(out, "README_HANDOFF.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write(render_readme(root, out, orders, species, counts, meta, scheduler, want_nc))
    with open(os.path.join(out, "CHECKLIST.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write(render_checklist(orders, scheduler))

    # MANIFEST.json: 每个文件的 sha256 (集群侧可核对没传坏)
    files = []
    for dirpath, dirnames, filenames in os.walk(out):
        dirnames[:] = [d for d in dirnames if d != "__pycache__"]
        for fn in sorted(filenames):
            if fn == "MANIFEST.json":
                continue
            p = os.path.join(dirpath, fn)
            files.append({"path": os.path.relpath(p, out).replace(os.sep, "/"),
                          "bytes": os.path.getsize(p), "sha256": sha256_file(p)})
    man = {
        "schema": "dft_handoff/v1", "handoff_version": VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_root": root, "orders": orders,
        "species": species, "counts": counts,
        "n_mag": meta.get("n_mag"), "S2": meta.get("S2"),
        "supercell": meta.get("supercell"),
        "spin_mode": ("noncollinear" + ("+soc" if has_soc else "")) if want_nc else "collinear",
        "scheduler": scheduler, "ncpu": ncpu, "walltime": walltime, "partition": partition,
        "coefficients_per_cell": meta.get("coefficients_per_cell"),
        "expects_back": [f"{o}/OUTCAR" for o in orders],
        "not_included": ["POTCAR (license: 集群侧自备)"],
        "files": files,
        "note": ("本包只含 VASP 输入。交付物 = 每个磁态一个通过收敛/口径/结构闸门的 OUTCAR; "
                 "回收走 `dft_handoff.py ingest` (10 条硬闸门)。"),
    }
    with open(os.path.join(out, "MANIFEST.json"), "w", encoding="utf-8") as f:
        json.dump(man, f, indent=2, ensure_ascii=False)

    say(f"[pack] 交接包已生成: {out}")
    say(f"       磁态 {len(orders)} 个: {orders}")
    say(f"       元素顺序 {species} = {counts};  自旋模式 {man['spin_mode']}")
    say(f"       调度器 {scheduler} (NCPU={ncpu} WALLTIME={walltime} PARTITION={partition})")
    say(f"       文件 {len(files)} 个, MANIFEST.json 记了 sha256")
    say(f"       集群上: export POTCAR_DIR=... VASP_BIN=... ; bash submit.sh")
    return man


# ---------------------------------------------------------------------------
def ingest_from(root, src, dry_run=False, force=False, lenient=None, min_kb=200,
                verbose=True):
    """从回收目录收 OUTCAR -> cases/。只认 <order>/OUTCAR 的相对路径。"""
    root = os.path.abspath(root)
    src = os.path.abspath(src)
    orders, meta = read_meta(root)
    cases = os.path.join(root, "02_supercell", "cases")
    say = print if verbose else (lambda *a, **k: None)

    found, missing = {}, []
    for o in orders:
        p = os.path.join(src, o, "OUTCAR")
        if os.path.isfile(p):
            found[o] = p
        else:
            missing.append(o)
    if missing:
        raise SystemExit(
            f"回收目录里缺这些磁态的 OUTCAR: {missing}\n"
            f"  期望: {os.path.join(src, '<order>', 'OUTCAR')}\n"
            f"  找到了: {sorted(found)}")

    # 走 ingest_dft 的同一套闸门: 用临时 runs/ 视图, 保证"只有一条收回通道"
    import tempfile
    td = tempfile.mkdtemp(prefix="handoff_ingest_")
    try:
        for o, p in found.items():
            d = os.path.join(td, "runs", o)
            os.makedirs(d, exist_ok=True)
            shutil.copy2(p, os.path.join(d, "OUTCAR"))
            meta_p = os.path.join(src, o, "run_meta.json")
            if os.path.isfile(meta_p):
                shutil.copy2(meta_p, os.path.join(d, "run_meta.json"))
        ok, problems, warns, stats, mp = IG.ingest(
            os.path.join(td, "runs"), cases, orders=orders, dry_run=dry_run,
            force=force, lenient=lenient, min_kb=min_kb, material=meta.get("material"),
            progress=(lambda m: say(m)))
    finally:
        shutil.rmtree(td, ignore_errors=True)

    if warns:
        say(f"\n[ingest] {len(warns)} 条警告:")
        for o, c, m in warns:
            say(f"  [warn] [{c}] {('(' + o + ') ') if o else ''}{m}")
    if not ok:
        say(f"\n[ingest] 拒收 ({len(problems)} 条硬问题):")
        for o, c, m in problems:
            say(f"  [FAIL] [{c}] {('(' + o + ') ') if o else ''}{m}")
        say("\n  这些都在本机 selftest 里跑过 (见 ingest_dft.py --selftest)。"
            "确认某条是误判才用 --lenient <检查名> 放行。")
        return 1
    if not dry_run:
        # provenance: 把 run_meta 索引进 manifest
        if mp:
            with open(mp, encoding="utf-8") as f:
                man = json.load(f)
            prov = {}
            for o in orders:
                m2 = os.path.join(src, o, "run_meta.json")
                if os.path.isfile(m2):
                    try:
                        prov[o] = json.load(open(m2, encoding="utf-8"))
                    except (OSError, ValueError):
                        prov[o] = {"_error": "run_meta.json 解析失败"}
            if prov:
                man["run_meta"] = prov
                with open(mp, "w", encoding="utf-8") as f:
                    json.dump(man, f, indent=2, ensure_ascii=False)
                say(f"[ingest] 已记录 {len(prov)} 个序的 run_meta (集群侧 provenance)")
        say(f"\n[ingest] 已收回 -> {cases}")
        say(f"  下一步: python run_2d_magnet.py solve --outdir {root}")
    return 0


# ---------------------------------------------------------------------------
def _selftest():
    import tempfile
    checks = []

    def ck(name, cond, extra=""):
        checks.append(bool(cond))
        print(f"  [{'ok ' if cond else 'FAIL'}] {name}" + (f"  {extra}" if extra else ""))

    td = tempfile.mkdtemp(prefix="handoff_")
    try:
        root = os.path.join(td, "Fe3Sn_pipeline")
        cases = os.path.join(root, "02_supercell", "cases")
        runs = os.path.join(root, "runs")
        orders = ["fm", "q0_120", "sq3_updntri", "sq3_alt"]
        species, counts = ["Fe", "Sn"], [108, 36]
        poscar = ("synth\n1.0\n9.3 0 0\n-4.65 8.05 0\n0 0 24\nFe Sn\n108 36\nDirect\n"
                  + "".join("0.1 0.2 0.3\n" for _ in range(sum(counts))))
        for o in orders:
            os.makedirs(os.path.join(cases, o), exist_ok=True)
            os.makedirs(os.path.join(runs, o), exist_ok=True)
            with open(os.path.join(cases, o, "POSCAR"), "w", encoding="utf-8") as f:
                f.write(poscar)
            with open(os.path.join(runs, o, "POSCAR"), "w", encoding="utf-8") as f:
                f.write(poscar)
            with open(os.path.join(runs, o, "INCAR"), "w", encoding="utf-8") as f:
                f.write("ENCUT = 520\nLNONCOLLINEAR = .TRUE.\nMAGMOM = 1 2 3\n")
            with open(os.path.join(runs, o, "KPOINTS"), "w", encoding="utf-8") as f:
                f.write("k\n1\nGamma\n3 3 1\n0 0 0\n")
        with open(os.path.join(cases, "magnetic_orders.json"), "w", encoding="utf-8") as f:
            json.dump({"selected_orders": orders, "n_mag": 108, "S2": 1.0,
                       "supercell": [3, 3, 1], "material": "Fe3Sn"}, f)
        os.makedirs(os.path.join(root, "03_dft"), exist_ok=True)
        with open(os.path.join(root, "03_dft", "POTCAR.spec"), "w", encoding="utf-8") as f:
            f.write("# POTCAR 组装说明\n元素顺序 **必须与 POSCAR 一致**:  Fe Sn\n")

        out = os.path.join(td, "handoff_out")
        man = pack(root, out, scheduler="slurm", ncpu=64, verbose=False)
        ck("pack: 四个序目录都在",
           all(os.path.isdir(os.path.join(out, o)) for o in orders))
        ck("pack: 每个序带 INCAR/KPOINTS/POSCAR",
           all(all(os.path.isfile(os.path.join(out, o, nm)) for nm in NEEDED)
               for o in orders))
        ck("pack: 产出 README/CHECKLIST/MANIFEST/submit.sh/jobscript.template",
           all(os.path.isfile(os.path.join(out, n)) for n in
               ("README_HANDOFF.md", "CHECKLIST.md", "MANIFEST.json", "submit.sh",
                "jobscript.template")))
        ck("pack: 没有打包 POTCAR (license)",
           not any(fn.upper().startswith("POTCAR") and fn != "POTCAR.spec"
                   for _d, _s, fs in os.walk(out) for fn in fs))
        ck("pack: MANIFEST 记了 sha256 且列了期望回收物",
           man["files"] and man["expects_back"] == [f"{o}/OUTCAR" for o in orders])
        ck("pack: 自旋模式从 INCAR 判定为非共线",
           man["spin_mode"] == "noncollinear", man["spin_mode"])
        js = open(os.path.join(out, "jobscript.template"), encoding="utf-8").read()
        ck("jobscript: 含 POTCAR 按元素顺序 cat", "for el in Fe Sn" in js)
        ck("jobscript: 含四个序的循环", "ORDERS=\"fm q0_120 sq3_updntri sq3_alt\"" in js)
        ck("jobscript: 提醒 vasp_ncl", "vasp_ncl" in js)
        ss = open(os.path.join(out, "submit.sh"), encoding="utf-8").read()
        ck("submit.sh: 五个调度器分支齐全",
           all(s in ss for s in ("slurm)", "pbs)", "sge)", "lsf)", "none)")))
        ck("submit.sh: slurm 模板带 #SBATCH",
           "#SBATCH --job-name=kagome_dft" in ss)

        # --- 回收: 好数据 -----------------------------------------------------
        recv = os.path.join(td, "received")
        S2, J = 1.0, (10.0, -2.0, 0.5)
        c = {"fm": [-54, -54, -81], "q0_120": [27, 27, -81],
             "sq3_updntri": [-6, -6, 27], "sq3_alt": [-6, 18, 27]}
        E = {o: round(sum(c[o][i] * J[i] * S2 for i in range(3)) / 1000.0, 8) for o in orders}
        for o in orders:
            os.makedirs(os.path.join(recv, o), exist_ok=True)
            IG.synth_outcar(os.path.join(recv, o, "OUTCAR"), E[o] + 100.0, steps=1,
                            accuracy=True, tot_moment=(150.0 if o == "fm" else 0.0))
            with open(os.path.join(recv, o, "run_meta.json"), "w", encoding="utf-8") as f:
                json.dump({"order": o, "host": "cluster01", "vasp_bin": "vasp_ncl"}, f)
        rc = ingest_from(root, recv, min_kb=0, verbose=False)
        ck("ingest: 好数据收回成功", rc == 0, f"rc={rc}")
        import dft_energy as DE
        got = {o: DE.read_toten(os.path.join(cases, o, "OUTCAR")) for o in orders}
        ck("ingest: 能量读回一致", max(abs(got[o] - E[o] - 100.0) for o in orders) < 1e-9)
        mp = os.path.join(cases, "ingest_manifest.json")
        mj = json.load(open(mp, encoding="utf-8"))
        ck("ingest: manifest 记了集群 provenance",
           (mj.get("run_meta") or {}).get("fm", {}).get("host") == "cluster01")

        # --- 回收: 缺一个序 ---------------------------------------------------
        os.remove(os.path.join(recv, "sq3_alt", "OUTCAR"))
        try:
            ingest_from(root, recv, min_kb=0, verbose=False)
            ck("ingest: 缺一个序报错", False)
        except SystemExit as e:
            ck("ingest: 缺一个序报错", "sq3_alt" in str(e), str(e).splitlines()[0][:56])

        # --- 回收: 旧文件充数 (与 cases 里逐位相同) ---------------------------
        IG.synth_outcar(os.path.join(recv, "sq3_alt", "OUTCAR"), E["sq3_alt"] + 100.0,
                        steps=1, accuracy=True, tot_moment=0.0)
        rc = ingest_from(root, recv, min_kb=0, verbose=False)
        ck("ingest: 旧 OUTCAR 充数被拒", rc == 1, f"rc={rc}")
        rc = ingest_from(root, recv, min_kb=0, force=True, verbose=False)
        ck("ingest: --force 放行", rc == 0, f"rc={rc}")

        # --- pack 预检: 缺文件必须报错 ---------------------------------------
        os.remove(os.path.join(runs, "fm", "INCAR"))
        try:
            pack(root, os.path.join(td, "out2"), verbose=False)
            ck("pack: runs 缺文件时拒绝打包", False)
        except SystemExit as e:
            ck("pack: runs 缺文件时拒绝打包", "INCAR" in str(e), str(e).splitlines()[0][:56])
    finally:
        shutil.rmtree(td, ignore_errors=True)

    n_ok = sum(checks)
    print(f"\n--selftest: {'ALL PASS' if n_ok == len(checks) else 'FAIL'} ({n_ok}/{len(checks)} 项)")
    return 0 if n_ok == len(checks) else 1


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="P4.2 kagome DFT 交接包 (pack) 与结果回收 (ingest)")
    sub = ap.add_subparsers(dest="cmd")

    p = sub.add_parser("pack", help="把 runs/ 打成可上传的交接包")
    p.add_argument("--root", required=True, help="项目目录 (<formula>_pipeline)")
    p.add_argument("--out", required=True, help="交接包输出目录")
    p.add_argument("--scheduler", default="auto",
                   choices=["auto"] + sorted(SCHEDULERS))
    p.add_argument("--ncpu", type=int, default=32)
    p.add_argument("--walltime", default=None)
    p.add_argument("--partition", default=None)
    p.add_argument("--quiet", action="store_true")

    q = sub.add_parser("ingest", help="把回收目录里的 OUTCAR 收进 cases/ (走硬闸门)")
    q.add_argument("--root", required=True)
    q.add_argument("--from", dest="src", required=True, help="回收目录 (含 <order>/OUTCAR)")
    q.add_argument("--dry-run", action="store_true")
    q.add_argument("--force", action="store_true")
    q.add_argument("--lenient", nargs="*", default=None)
    q.add_argument("--min-kb", type=int, default=200)
    q.add_argument("--quiet", action="store_true")

    ap.add_argument("--selftest", action="store_true")
    _argv = list(sys.argv[1:] if argv is None else argv)
    if "--selftest" in _argv:
        return _selftest()
    args = ap.parse_args(_argv)

    if args.cmd == "pack":
        pack(args.root, args.out, scheduler=args.scheduler, ncpu=args.ncpu,
             walltime=args.walltime, partition=args.partition, verbose=not args.quiet)
        return 0
    if args.cmd == "ingest":
        return ingest_from(args.root, args.src, dry_run=args.dry_run, force=args.force,
                           lenient=args.lenient, min_kb=args.min_kb,
                           verbose=not args.quiet)
    ap.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
