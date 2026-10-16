#!/usr/bin/env python
"""ingest_dft.py —— 「runs/ -> cases/」的唯一收回通道（P4.2）。

## 为什么需要这个脚本

`assemble_dft.py` 把每个磁态组装成 `runs/<order>/{INCAR,KPOINTS,POSCAR}` 交付给外部
DFT; 算完之后需要把 4 个 `<order>/OUTCAR` 搬回 `cases/<order>/OUTCAR`, `solve` 才能反解 J。
在 P4.2 之前这一步是**手工 cp**——而这正是整条链上最危险的一步:

  * 把还没收敛的 OUTCAR 收回来 -> J 被误差污染, 但**看起来完全正常**;
  * 把另一台机器/另一次设置的 OUTCAR 混进来 -> `free energy TOTEN` 的口径不同
    (`mixed_convention`), J 直接错;
  * 某个序的 OUTCAR 其实是**另一个磁序**的 -> 设计矩阵的行与能量错配, J 乱码;
  * 之前跑过的旧 OUTCAR 留在 cases/ 里 -> 你以为 4 个序齐了, 其实其中一个没重算。

本脚本把这些都变成**硬闸门**: 先全量校验, 全部通过才动 `cases/`; 任何一条不过就一件不写
(原子性), 并打印"哪一条、为什么、怎么修"。

## 硬闸门 (默认 ALL 打开)

| # | 检查 | 不过会怎样 |
|---|---|---|
| 1 | 四态 OUTCAR 齐全 (`magnetic_orders.json.selected_orders`) | 缺哪个报哪个 |
| 2 | 文件大小 >= `--min-kb` | 拒绝: 很可能是被 kill 的空壳 |
| 3 | `General timing and accounting` 段存在 | 拒绝: VASP 没跑到正常结束 |
| 4 | SCF 收敛 (`aborting loop because EDIFF is reached`) | 拒绝: 未收敛的能量不是能量 |
| 5 | **恰好一个**离子步 (静态计算) | 拒绝: 松弛过的结构与 `cases/` 的 POSCAR/系数不对应 |
| 6 | `reached required accuracy` 或 `reached required accuracy - stopping structural energy minimisation` | 警告 (静态 `NSW=1` 不打印这句, 属正常) |
| 7 | 能量口径统一 (`dft_energy.mixed_convention`) | 拒绝: 混口径的能量不可相减 |
| 8 | 结构指纹与参考序一致 (原子数/晶格/组成) | 拒绝: 不同结构算出的能量不可比 |
| 9 | 与 `cases/` 现有 OUTCAR 的 sha256 不同 (除非 `--force`) | 拒绝: 你以为重算了, 其实是旧的 |
| 10 | 能量离散度 < `--max-spread` eV/胞 | 拒绝: 多半是 INCAR 不同/磁矩没生效/几何不同 |

闸门 5/8 的**目的**: 反解 J 的前提是"唯一变量是磁序"。价格是外部的 OUTCAR 必须是
**同一个结构、同一次数值设置**下、**静态**算出来的。本脚本无法验证"数值设置"本身
(ENCUT/K 点/POTCAR 版本), 只能把能看到的指纹记进 `ingest_manifest.json` 供人工核对。

## 能量口径 (P3.17)

一律走 `dft_energy.py`: `free energy TOTEN`、取**最后一次**出现、缺了就报错——
不静默改用 `energy(sigma->0)`。

## 产物

  * `cases/<order>/OUTCAR`            —— 校验过的能量来源 (solve 读这个)
  * `cases/<order>/OUTCAR.sha256`     —— 便于日后核对"是不是同一份"
  * `cases/ingest_manifest.json`      —— 每个序的 sha256 / 能量 / 收敛标志 / 自旋模式 /
    磁矩诊断 / 警告, 加全局 `spread_eV`, `source_runs`, 脚本版本
  * `cases/ingest_report.md`          —— 人能读的一页 (交接收回后先看这个)

## 用法

  # 校验 + 搬进 cases/ (全部通过才写)
  python ingest_dft.py --runs <proj>/runs --cases <proj>/02_supercell/cases

  # 只校验不写 (退出码 0/1, 交接方先自查)
  python ingest_dft.py --runs <proj>/runs --cases <...>/cases --dry-run

  # 无 DFT 自检 (合成 OUTCAR, 同时验证"坏 OUTCAR 会被拒")
  python ingest_dft.py --selftest
"""
import argparse
import glob
import hashlib
import json
import os
import re
import shutil
import sys
import tempfile
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dft_energy as DE  # noqa: E402

VERSION = "1.0"

# ---- OUTCAR 里的确定性标志 (都是 VASP 会打印的固定字符串, 不是猜的) ----------
RE_TIMING = re.compile(r"General timing and accounting")
RE_SCF_OK = re.compile(r"aborting loop because EDIFF is reached", re.I)
RE_ACCURACY = re.compile(r"reached required accuracy", re.I)
RE_NIONS = re.compile(r"ions per type\s*=\s*([0-9\s]+)")
RE_NKPT = re.compile(r"NKPTS\s*=\s*(\d+)")
RE_ENCUT = re.compile(r"ENCUT\s*=\s*([-+0-9.Ee]+)")
RE_LNONCOL = re.compile(r"LNONCOLLINEAR\s*=\s*(T|F|\.TRUE\.|\.FALSE\.)", re.I)
RE_LSORBIT = re.compile(r"LSORBIT\s*=\s*(T|F|\.TRUE\.\.?|\.FALSE\.)", re.I)
RE_ISPIN = re.compile(r"ISPIN\s*=\s*(\d)")
RE_VASPVER = re.compile(r"vasp\.([0-9.]+)")
RE_MAG_BLOCK = re.compile(r"magnetization \(x\)")
RE_POSCAR_LATT = re.compile(r"direct lattice vectors", re.I)

# 每磁位点矩的**诊断**阈值 (仅用于"提醒人看一眼", 不作为拒绝依据)
DIAG_FM_MIN = 0.5      # 铁磁序: 每磁位点净矩低于此值 -> 提醒核对 MAGMOM
DIAG_AFM_MAX = 0.5     # 反铁磁/非共线序: 每磁位点净矩高于此值 -> 提醒核对 MAGMOM


# ---------------------------------------------------------------------------
# 文本工具
# ---------------------------------------------------------------------------
def _read_text(path, limit_bytes=None):
    """读 OUTCAR (ASCII/latin-1 兜底: 有些 VASP 构建会打印非 UTF-8 字节)。

    limit_bytes 用于只读头部/尾部——大 OUTCAR (几百 MB) 不该整个进内存。
    """
    if limit_bytes is None:
        with open(path, "rb") as f:
            raw = f.read()
    else:
        with open(path, "rb") as f:
            raw = f.read(limit_bytes)
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def _read_tail(path, nbytes=65536):
    size = os.path.getsize(path)
    with open(path, "rb") as f:
        if size > nbytes:
            f.seek(size - nbytes)
        raw = f.read()
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1")


def sha256_file(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _find_outcar(d):
    """在一个序目录里找 OUTCAR。优先 `OUTCAR`; 其次容忍 `OUTCAR.gz` 之外的常见改名。"""
    for nm in ("OUTCAR", "OUTCAR.txt", "outcar", "OUTCAR.static"):
        p = os.path.join(d, nm)
        if os.path.isfile(p):
            return p
    hits = [p for p in glob.glob(os.path.join(d, "OUTCAR*")) if os.path.isfile(p)]
    return sorted(hits)[0] if hits else None


# ---------------------------------------------------------------------------
# OUTCAR 解析
# ---------------------------------------------------------------------------
def parse_outcar(path, min_kb=200):
    """把一个 OUTCAR 解析成结构化 dict。**只解析, 不判断** (判断在 validate 里)。

    返回 dict:
      n_bytes, size_ok, has_timing, scf_converged, accuracy_reached, n_ionic_steps,
      energy_eV (经 dft_energy, 缺失 -> None 且 energy_error 记原因),
      sigma0, without_entropy, mixed_convention, n_ions, n_kpts, encut,
      noncollinear, lsorbit, ispin, vasp_version, has_mag_block,
      tot_moment_raw, moment_error
    """
    size = os.path.getsize(path)
    out = {
        "path": os.path.abspath(path),
        "n_bytes": size,
        "size_kb": round(size / 1024.0, 1),
        "size_ok": size >= min_kb * 1024,
        "has_timing": False,
        "scf_converged": False,
        "accuracy_reached": False,
        "n_ionic_steps": None,
        "energy_eV": None,
        "energy_error": None,
        "sigma0": None,
        "without_entropy": None,
        "mixed_convention": False,
        "n_ions": None,
        "n_kpts": None,
        "encut": None,
        "noncollinear": None,
        "lsorbit": None,
        "ispin": None,
        "vasp_version": None,
        "has_mag_block": False,
        "tot_moment_raw": None,
        "moment_error": None,
        "lattice": None,
    }

    head = _read_text(path, limit_bytes=8 << 20)
    tail = _read_tail(path)

    out["has_timing"] = bool(RE_TIMING.search(head) or RE_TIMING.search(tail))
    out["scf_converged"] = bool(RE_SCF_OK.search(head) or RE_SCF_OK.search(tail))
    out["accuracy_reached"] = bool(RE_ACCURACY.search(head) or RE_ACCURACY.search(tail))
    # 能量: 唯一口径入口 (dft_energy 是 P3.17 的唯一实现, 这里只转发)
    d = DE.read_energy_detail(path)
    out["n_ionic_steps"] = d.get("n_ionic_steps")
    try:
        out["energy_eV"] = DE.read_toten(path)
    except Exception as e:                                   # noqa: BLE001
        out["energy_error"] = str(e)
    out["sigma0"] = d.get("sigma0")
    out["without_entropy"] = d.get("without_entropy")
    out["mixed_convention"] = DE.mixed_convention(path)

    m = RE_NIONS.search(head)
    if m:
        counts = [int(x) for x in m.group(1).split()]
        out["n_ions"] = sum(counts)
        out["ion_counts"] = counts
    m = RE_NKPT.search(head)
    if m:
        out["n_kpts"] = int(m.group(1))
    m = RE_ENCUT.search(head)
    if m:
        try:
            out["encut"] = float(m.group(1))
        except ValueError:
            pass
    m = RE_LNONCOL.search(head)
    if m:
        out["noncollinear"] = m.group(1).upper().startswith("T")
    m = RE_LSORBIT.search(head)
    if m:
        out["lsorbit"] = m.group(1).upper().startswith("T")
    m = RE_ISPIN.search(head)
    if m:
        out["ispin"] = int(m.group(1))
    m = RE_VASPVER.search(head)
    if m:
        out["vasp_version"] = m.group(1)

    # 结构指纹: OUTCAR 自己的晶格矢量 (与 cases/<order>/POSCAR 逐位比对, 闸门 8)
    fp = outcar_fingerprint(head)
    if fp:
        out["lattice"] = fp["lattice"]

    # 磁矩诊断: OUTCAR 尾部 "magnetization (x)" 表中的 totals
    out["has_mag_block"] = bool(RE_MAG_BLOCK.search(tail))
    if out["has_mag_block"]:
        try:
            out["tot_moment_raw"] = _parse_total_moment_raw(tail)
        except Exception as e:                               # noqa: BLE001
            out["moment_error"] = str(e)
    return out


def _parse_total_moment_raw(tail):
    """从 OUTCAR 的 `magnetization (x)` 表里取最后一行的 `tot` 列。

    VASP 这张表形如:
        # of ion       s       p       d       ...
        -----------------------------------------
            1        0.001   0.002   1.998
        -----------------------------------------
         tot          0.001   0.002   1.998

    返回 float 或 3-tuple (非共线有 x/y/z 三张表时取每张的 tot)。
    这里**只做诊断**, 解析不到就返回 None —— 绝不用它做物理判断。
    """
    idxs = [m.start() for m in RE_MAG_BLOCK.finditer(tail)]
    if not idxs:
        return None
    totals = []
    for i in idxs:
        seg = tail[i:i + 400000]
        stop = seg.find("magnetization (y)")
        stop2 = seg.find("magnetization (z)")
        for s in (stop, stop2):
            if s > 0:
                seg = seg[:s]
        vals = None
        for line in seg.splitlines():
            t = line.strip()
            if t.startswith("tot"):
                parts = t.split()
                try:
                    vals = [float(x) for x in parts[1:] if x not in ("-",)]
                except ValueError:
                    vals = None
        totals.append(vals)
    return totals


# ---------------------------------------------------------------------------
# 结构指纹 (闸门 8)
# ---------------------------------------------------------------------------
def supercell_fingerprint(case_dir):
    """从 `cases/<order>/POSCAR` 取结构指纹: 原子数 / 晶格矢量 / 元素组成。

    反解 J 的前提是"唯一变量是磁序" ⇒ 四个序的**结构必须逐位相同**。这里给的是
    cases 侧 (即 DFT 输入的) 指纹; ingest 只核对 OUTCAR 里的原子数与元素块,
    晶格矢量需要 OUTCAR 侧解析 (`direct lattice vectors`) —— 见 `outcar_fingerprint`。
    """
    p = os.path.join(case_dir, "POSCAR")
    if not os.path.isfile(p):
        return None
    lines = [ln for ln in _read_text(p).splitlines() if ln.strip()]
    try:
        species = lines[5].split()
        counts = [int(x) for x in lines[6].split()]
        scale = float(lines[1].split()[0])
        lat = [[float(x) * scale for x in lines[2 + i].split()[:3]] for i in range(3)]
    except (IndexError, ValueError):
        return None
    return {"n_atoms": sum(counts), "species": species, "counts": counts,
            "n_mag_atoms": None, "lattice": lat}


def outcar_fingerprint(text_head):
    """从 OUTCAR 头部取 `direct lattice vectors` 的 3x3 (Å) 与 `ions per type`。"""
    m = RE_POSCAR_LATT.search(text_head)
    if not m:
        return None
    seg = text_head[m.end():m.end() + 2000]
    vecs = []
    for line in seg.splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            xyz = [float(x) for x in parts[:3]]
        except ValueError:
            continue
        vecs.append(xyz)
        if len(vecs) == 3:
            break
    if len(vecs) != 3:
        return None
    return {"lattice": vecs}


# ---------------------------------------------------------------------------
# 校验
# ---------------------------------------------------------------------------
# hard=True 的检查不过 -> 拒绝写入; hard=False 只记警告
CHECKS = [
    ("files_present", True, "四态 OUTCAR 齐全"),
    ("size", True, "OUTCAR 体积不小于阈值"),
    ("timing", True, "VASP 正常跑完 (有 timing 段)"),
    ("scf", True, "SCF 收敛 (EDIFF reached)"),
    ("ionic_steps", True, "恰好一个离子步 (静态)"),
    ("convention", True, "能量口径统一 (free energy TOTEN)"),
    ("structure", True, "与参考序结构一致"),
    ("not_shadowed", True, "不是 cases/ 里的旧 OUTCAR (sha256 必须变)"),
    ("spread", True, "能量离散度在阈值内"),
    ("accuracy", False, "打印了 reached required accuracy (静态可不打印)"),
    ("moment", False, "每磁位点矩与磁序类型相符 (诊断)"),
]


def check_policy(name):
    for nm, hard, _desc in CHECKS:
        if nm == name:
            return hard
    raise KeyError(name)


def validate(entries, ref_fp, max_spread, min_energy_gap, policy_override, quiet=False):
    """对解析好的 entries 跑全部检查。

    entries: {order: parsed dict}
    返回 (problems, warns, stats), problems 里每条是 (order_or_None, check_name, msg)
    """
    problems = []
    warns = []
    stats = {}

    def add(order, check, msg):
        hard = policy_override.get(check, check_policy(check))
        (problems if hard else warns).append((order, check, msg))

    fms = [e for o, e in entries.items() if e.get("energy_eV") is not None]
    for order, e in sorted(entries.items()):
        if e.get("energy_error"):
            add(order, "convention", f"读不出 `free energy TOTEN`: {e['energy_error']}")
        if not e["size_ok"]:
            add(order, "size",
                f"只有 {e['size_kb']} KB (阈值 {max(1, expected_min_kb)} KB) —— "
                "很可能是被 kill 掉的空壳")
        if not e["has_timing"]:
            add(order, "timing", "没有 `General timing and accounting` 段 —— VASP 没跑到正常结束")
        if not e["scf_converged"]:
            add(order, "scf", "没有 `aborting loop because EDIFF is reached` —— SCF 未收敛, "
                              "这个能量不是能量")
        if e["n_ionic_steps"] != 1:
            add(order, "ionic_steps",
                f"有 {e['n_ionic_steps']} 个 `free  energy   TOTEN` 出现次数 (期望 1) —— "
                "静态计算必须只跑一个离子步; 松弛过的结构与 cases/ 里的系数不对应")
        if e["mixed_convention"]:
            add(order, "convention", "同一文件内混用了不同能量口径 (free energy vs energy(sigma->0))")
        if not e["accuracy_reached"]:
            add(order, "accuracy", "没打印 `reached required accuracy` "
                                   "(静态 NSW=1 一般不打印, 属正常; 提醒确认一下)")

    # 结构一致性: OUTCAR 里的原子数必须与参考 POSCAR 的原子数一致
    if ref_fp and ref_fp.get("n_atoms"):
        for order, e in sorted(entries.items()):
            if e.get("n_ions") is not None and e["n_ions"] != ref_fp["n_atoms"]:
                add(order, "structure",
                    f"OUTCAR 的 ions per type 合计 = {e['n_ions']}, 但 cases 的 POSCAR 有 "
                    f"{ref_fp['n_atoms']} 个原子 —— 不是同一个结构, 能量不可比")
            if e.get("ion_counts") and ref_fp.get("counts") and \
                    e["ion_counts"] != ref_fp["counts"]:
                add(order, "structure",
                    f"元素块计数 {e['ion_counts']} != POSCAR 的 {ref_fp['counts']}")
            # 晶格矢量逐位比对 (容差 1e-4 Å): 四个序必须算的是**同一个胞**。
            # 只比 OUTCAR 的 `direct lattice vectors` 与 cases/<order>/POSCAR 的前三行。
            rl = ref_fp.get("lattice")
            ol = e.get("lattice")
            if rl and ol:
                worst = max(abs(rl[i][j] - ol[i][j]) for i in range(3) for j in range(3))
                stats.setdefault("lattice_max_dev", {})[order] = worst
                if worst > 1e-4:
                    add(order, "structure",
                        f"OUTCAR 的晶格矢量与 cases 的 POSCAR 最大偏差 {worst:.3g} Å "
                        f"(> 1e-4) —— 算的不是同一个胞, 能量不可比")
            elif rl and ol is None:
                add(order, "structure",
                    "OUTCAR 里读不到 `direct lattice vectors` —— 无法核对胞, 结构一致性存疑")

    # 自旋模式一致性 (四态必须同一种: 要么都非共线, 要么都共线)
    nc = {o: e.get("noncollinear") for o, e in entries.items()}
    if len({v for v in nc.values() if v is not None}) > 1:
        add(None, "structure",
            f"四个序的自旋模式不一致 (LNONCOLLINEAR: {nc}) —— 必须同一次设置下算出来")

    stats["noncollinear"] = nc
    stats["lsorbit"] = {o: e.get("lsorbit") for o, e in entries.items()}
    stats["encut"] = {o: e.get("encut") for o, e in entries.items()}
    stats["n_kpts"] = {o: e.get("n_kpts") for o, e in entries.items()}
    stats["vasp_version"] = {o: e.get("vasp_version") for o, e in entries.items()}

    # 数值设置一致性 (只警告: OUTCAR 里不一定都打印全)
    for key, label in (("encut", "ENCUT"), ("n_kpts", "NKPTS"), ("vasp_version", "VASP 版本"),
                       ("lsorbit", "LSORBIT")):
        vals = {o: v for o, v in stats[key].items() if v is not None}
        if len(set(map(str, vals.values()))) > 1:
            add(None, "structure", f"四个序的 {label} 不一致: {vals} —— "
                                   "必须同一次数值设置下算出来")

    # 能量离散度 (闸门 10): 三个非 fm 序相对 fm 的能量差
    if fms:
        e_fm = entries.get("fm", {}).get("energy_eV")
        stats["spread_eV"] = round(max(e["energy_eV"] for e in fms)
                                   - min(e["energy_eV"] for e in fms), 6)
        if stats["spread_eV"] > max_spread:
            add(None, "spread",
                f"能量离散度 {stats['spread_eV']:.4f} eV/胞 > 阈值 {max_spread} —— "
                "多半是 INCAR 不同 / 磁矩没生效 / 几何不同。先核对, 不要直接反解 J。")
        if e_fm is not None:
            gaps = {o: round(e["energy_eV"] - e_fm, 6) for o, e in entries.items()}
            stats["gap_vs_fm_eV"] = gaps
            tiny = [o for o, g in gaps.items() if o != "fm" and abs(g) < min_energy_gap]
            if tiny:
                add(None, "spread",
                    f"序 {tiny} 与 fm 的能量差 < {min_energy_gap} eV —— 该磁序很可能没有生效 "
                    "(MAGMOM 没进 INCAR / 自旋没磁化)。这正是 HANDOFF 里最常见的失败模式。")
            fm_row = [e for o, e in entries.items() if o == "fm"]
            n_mag = (ref_fp or {}).get("n_mag_atoms")
            if fm_row and n_mag:
                tm = fm_row[0].get("tot_moment_raw")
                stats["fm_total_moment_raw"] = tm
                v = _first_tot(tm)
                if v is not None and abs(v) / n_mag < DIAG_FM_MIN:
                    add("fm", "moment",
                        f"fm 的每磁位点净矩 ≈ {abs(v) / n_mag:.3f} μB < {DIAG_FM_MIN} —— "
                        "铁磁序应当有明显净矩, 请核对 MAGMOM/是否真的磁化")
            for o, e in sorted(entries.items()):
                if o == "fm" or not n_mag:
                    continue
                v = _first_tot(e.get("tot_moment_raw"))
                if v is not None and abs(v) / n_mag > DIAG_AFM_MAX:
                    add(o, "moment",
                        f"每磁位点净矩 ≈ {abs(v) / n_mag:.3f} μB > {DIAG_AFM_MAX} —— "
                        f"{o} 是反铁磁/非共线序, 净矩应当接近 0, 请核对")
    return problems, warns, stats


def _first_tot(tm):
    """tot_moment_raw 可能是 float / list-of-lists (x,y,z 三张表) -> 取第一个标量和。"""
    if tm is None:
        return None
    if isinstance(tm, (int, float)):
        return float(tm)
    for row in tm:
        if row:
            try:
                return float(sum(row))
            except (TypeError, ValueError):
                continue
    return None


expected_min_kb = 200   # 模块级, 供 validate 的提示信息引用 (main 里可覆盖)


# ---------------------------------------------------------------------------
# 写入 (原子)
# ---------------------------------------------------------------------------
def write_outputs(cases_dir, entries, sources, stats, problems, warns, args,
                  manifest_extra=None):
    """把校验过的 OUTCAR 搬进 cases/, 并写 manifest + 报告。**只在 problems 为空时调用。**"""
    os.makedirs(cases_dir, exist_ok=True)
    moved = []
    for order, e in sorted(entries.items()):
        dst_dir = os.path.join(cases_dir, order)
        os.makedirs(dst_dir, exist_ok=True)
        dst = os.path.join(dst_dir, "OUTCAR")
        src = sources[order]
        if os.path.abspath(src) != os.path.abspath(dst):
            shutil.copy2(src, dst)
        h = sha256_file(dst)
        with open(dst + ".sha256", "w", encoding="utf-8") as f:
            f.write(f"{h}  OUTCAR\n")
        moved.append({"order": order, "outcar": os.path.relpath(dst, cases_dir),
                      "sha256": h, "bytes": os.path.getsize(dst)})

    manifest = {
        "schema": "dft_ingest/v1",
        "ingest_version": VERSION,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "energy_convention": DE.CONVENTION,
        "source_runs": os.path.abspath(args.runs) if args.runs else None,
        "cases_dir": os.path.abspath(cases_dir),
        "energy_tolerance_eV": DE.MIX_TOL,
        "checks": [{"name": n, "hard": bool(policy_get(args, n)), "desc": d}
                   for n, _h, d in CHECKS],
        "files": moved,
        "orders": {o: {k: e.get(k) for k in
                       ("energy_eV", "sigma0", "without_entropy", "n_ionic_steps",
                        "scf_converged", "accuracy_reached", "has_timing", "size_kb",
                        "n_ions", "n_kpts", "encut", "noncollinear", "lsorbit", "ispin",
                        "vasp_version", "has_mag_block", "tot_moment_raw")}
                   for o, e in sorted(entries.items())},
        "stats": stats,
        "warnings": [{"order": o, "check": c, "msg": m} for o, c, m in warns],
        "problems": [],
        "note": ("本 manifest 只证明「OUTCAR 通过了格式/收敛/结构/口径闸门」, "
                 "**不证明** J 与 T_N 的物理正确性 —— 那需要真实 DFT 结果 + 文献锚点核对。"),
    }
    if manifest_extra:
        manifest.update(manifest_extra)
    mpath = os.path.join(cases_dir, "ingest_manifest.json")
    with open(mpath, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    with open(os.path.join(cases_dir, "ingest_report.md"), "w", encoding="utf-8") as f:
        f.write(render_report(manifest))
    return mpath, moved


def policy_get(args, name):
    """检查项是否仍然硬 (--lenient <name> 可把某项降为警告)。"""
    if name in (args.lenient or []):
        return False
    return check_policy(name)


def render_report(man):
    L = ["# DFT 回收报告 (ingest_dft.py)", ""]
    L.append(f"- 生成时间 (UTC): `{man['generated_utc']}`")
    L.append(f"- 来源 runs: `{man['source_runs']}`")
    L.append(f"- 目标 cases: `{man['cases_dir']}`")
    L.append(f"- 能量口径: `{man['energy_convention']}` (混口径容差 {man['energy_tolerance_eV']} eV)")
    L.append("- **本报告只证明 OUTCAR 过了格式/收敛/结构/口径闸门, 不证明 J/T_N 的物理正确性。**")
    L.append("")
    L.append("## 四个磁态")
    L.append("")
    L.append("| 序 | 能量 (eV) | ionic steps | SCF 收敛 | 与 fm 之差 (meV) | sha256 前 12 位 |")
    L.append("|---|---|---|---|---|---|")
    gaps = (man.get("stats") or {}).get("gap_vs_fm_eV") or {}
    for fo in man["files"]:
        o = fo["order"]
        e = man["orders"][o]
        g = gaps.get(o)
        L.append(f"| {o} | {e['energy_eV']!r} | {e['n_ionic_steps']} "
                 f"| {'是' if e['scf_converged'] else '**否**'} "
                 f"| {('%.3f' % (g * 1000)) if g is not None else '-'} "
                 f"| `{fo['sha256'][:12]}` |")
    st = man.get("stats") or {}
    L.append("")
    L.append(f"- 能量离散度: **{st.get('spread_eV')} eV/胞** (阈值见 manifest `checks`)")
    L.append(f"- LNONCOLLINEAR: `{st.get('noncollinear')}`   LSORBIT: `{st.get('lsorbit')}`")
    L.append(f"- ENCUT: `{st.get('encut')}`   NKPTS: `{st.get('n_kpts')}`   "
             f"VASP: `{st.get('vasp_version')}`")
    L.append("")
    if man.get("warnings"):
        L.append("## 警告 (不拦写入, 但请核对)")
        L.append("")
        for w in man["warnings"]:
            L.append(f"- `[{w['check']}]` {('(' + w['order'] + ') ') if w['order'] else ''}{w['msg']}")
        L.append("")
    L.append("## 下一步")
    L.append("")
    L.append("```")
    L.append("# 反解 J (kagome)")
    L.append("python run_2d_magnet.py solve --outdir <proj>")
    L.append("# 或直接")
    L.append("python scripts/solve_J_kagome.py <proj>/02_supercell/cases")
    L.append("```")
    return "\n".join(L) + "\n"


# ---------------------------------------------------------------------------
# 合成 OUTCAR (自检用; 与 tests 里的 mock 同构)
# ---------------------------------------------------------------------------
def synth_outcar(path, free_eV, *, sigma0=None, noent=None, steps=1, n_ions=None,
                 counts="108 36", n_kpts=6, encut=520.0, noncollinear=True,
                 lsorbit=True, vasp="6.4.2", timing=True, scf=True, accuracy=False,
                 tot_moment=None, lattice=None):
    """写一个**最小但结构正确**的 OUTCAR: 只含 validate 会去查的那些标志。

    这是自检/单测用的合成文件, **与真实 VASP 输出不等价** —— 只为验证闸门逻辑。
    """
    L = []
    L.append(f" vasp.6.4.2  (build date)  complex")
    L.append(" running on 64 total cores")
    L.append("")
    L.append(" POSCAR: synth")
    if counts:
        L.append(f"    ions per type = {counts}")
    L.append(f"   NKPTS = {n_kpts}   k-points in BZ")
    L.append(f"   ENCUT  =  {encut:.1f} eV")
    L.append(f"   ISPIN  =  2" if not noncollinear else "   ISPIN  =  1")
    L.append(f"   LNONCOLLINEAR = {'T' if noncollinear else 'F'}")
    L.append(f"   LSORBIT = {'T' if lsorbit else 'F'}")
    L.append("")
    L.append(" direct lattice vectors                 reciprocal lattice vectors")
    for v in (lattice or [[9.3, 0.0, 0.0], [-4.65, 8.05, 0.0], [0.0, 0.0, 24.0]]):
        L.append("  " + "  ".join(f"{x:10.6f}" for x in v))
    L.append("")
    for i in range(int(steps)):
        if i > 0:
            L.append(f" ----------------------------------- Iteration {i} -----------------------------------")
        if scf:
            L.append("  DAVID: aborting loop because EDIFF is reached")
        if accuracy:
            L.append("  reached required accuracy")
        L.append(f"  free  energy   TOTEN  =      {free_eV:+.8f} eV")
        L.append(f"  energy  without entropy=     {(noent if noent is not None else free_eV):+.8f}  "
                 f"energy(sigma->0) =     {(sigma0 if sigma0 is not None else free_eV):+.8f}")
    L.append("")
    if tot_moment is not None:
        L.append(" magnetization (x)")
        L.append("")
        L.append(" # of ion       s       p       d       tot")
        L.append(" ------------------------------------------")
        for i in range(6):
            L.append(f"     {i + 1}      0.001   0.002   {(tot_moment / 6.0):8.3f}")
        L.append(" ------------------------------------------")
        L.append(f"    tot         {0.001 * 6:.3f}   {0.002 * 6:.3f}   {tot_moment:.3f}")
        L.append("")
        L.append(" magnetization (y)")
        L.append("")
        L.append(" # of ion       s       p       d       tot")
        L.append(" ------------------------------------------")
        L.append(f"    tot         {0.0:.3f}   {0.0:.3f}   {0.0:.3f}")
        L.append("")
    if timing:
        L.append(" General timing and accounting informations for this run:")
        L.append(" ----------------------------------------------------------------")
        L.append("      Total CPU time used (sec):      123.456")
        L.append("")
    L.append(" Voluntary context switches:  0")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    return path


def _selftest_only_hard(problems):
    return [p for p in problems if p[1] in [n for n, h, _d in CHECKS if h]]


# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
def collect(runs_dir, orders, cases_dir, quiet=False):
    """从 runs_dir 收 OUTCAR。返回 (entries, sources, missing)。"""
    entries, sources, missing = {}, {}, []
    for o in orders:
        d = os.path.join(runs_dir, o)
        p = _find_outcar(d) if os.path.isdir(d) else None
        if p is None:
            missing.append(o)
            continue
        entries[o] = parse_outcar(p, min_kb=expected_min_kb)
        sources[o] = p
    return entries, sources, missing


def ingest(runs_dir, cases_dir, orders=None, max_spread=2.0, min_energy_gap=1e-4,
           dry_run=False, force=False, lenient=None, min_kb=200, quiet=False,
           material=None, progress=None):
    """runs/ -> cases/ 的完整通道。返回 (ok: bool, problems, warns, stats, manifest_path)。

    原子性: 先全量校验, 有问题就 **一件不写** (dry_run 时永远不写)。
    """
    global expected_min_kb
    expected_min_kb = int(min_kb)
    say = progress or (lambda m: None)

    meta_path = os.path.join(cases_dir, "magnetic_orders.json")
    meta = {}
    if os.path.isfile(meta_path):
        with open(meta_path, encoding="utf-8-sig") as f:
            meta = json.load(f)
    if orders is None:
        orders = list(meta.get("selected_orders") or [])
        if not orders:
            orders = sorted(d for d in os.listdir(cases_dir)
                            if os.path.isdir(os.path.join(cases_dir, d))) \
                if os.path.isdir(cases_dir) else []
        if not orders:
            raise SystemExit(f"无法确定磁序: {cases_dir} 没有 magnetic_orders.json, "
                             "也没有子目录。先跑 build。")
    orders = list(orders)

    say(f"[ingest] 磁序 {orders}  <-  runs: {os.path.abspath(runs_dir)}")
    entries, sources, missing = collect(runs_dir, orders, cases_dir, quiet)
    if missing:
        raise SystemExit(
            f"runs 里缺这些磁序的 OUTCAR: {missing}\n"
            f"  期望路径: {os.path.join(os.path.abspath(runs_dir), '<order>', 'OUTCAR')}\n"
            f"  四个序必须都算完才能反解 J (设计矩阵 rank 需要四个独立能量)。")

    # 参考结构指纹: 取第一个序的 POSCAR
    ref_fp = None
    for o in orders:
        ref_fp = supercell_fingerprint(os.path.join(cases_dir, o))
        if ref_fp:
            break
    if ref_fp is not None:
        n_mag = meta.get("n_mag")
        if n_mag:
            ref_fp["n_mag_atoms"] = int(n_mag)
    else:
        say("[ingest] [warn] cases/ 里没有可读的 POSCAR -> 跳过结构指纹核对")

    class _A:
        """write_outputs / policy_get 需要的最小参数视图 (避免在这里再开一份 argparse)。"""
        pass

    a = _A()
    a.lenient = list(lenient or [])
    a.runs = runs_dir

    # policy_override 只放**被 --lenient 改动过的**项; 其余走 CHECKS 里的默认 hard 值。
    # (曾经的 bug: 把全部检查都塞进 override, 于是 --lenient 永远不生效。)
    override = {n: False for n in a.lenient}
    problems, warns, stats = validate(entries, ref_fp, max_spread, min_energy_gap,
                                      override, quiet=quiet)

    # 闸门 9: 不许用 cases/ 里的旧 OUTCAR 充数
    if not force:
        old = {}
        for o in orders:
            op = os.path.join(cases_dir, o, "OUTCAR")
            if os.path.isfile(op):
                old[o] = sha256_file(op)
        shadowed = [o for o in orders if o in old and old[o] == sha256_file(sources[o])]
        if shadowed:
            problems.append((None, "not_shadowed",
                             f"序 {shadowed} 的 OUTCAR 与 cases/ 里已有的**逐位相同** "
                             "(sha256 一致) —— 你以为重算了, 其实是旧的。"
                             "确认无误后加 --force, 或重新跑 VASP。"))

    ok = not problems
    if not ok:
        say(f"[ingest] 拒绝写入: {len(problems)} 条硬问题")
        return False, problems, warns, stats, None
    if dry_run:
        say(f"[ingest] --dry-run: 全部闸门通过 ({len(warns)} 条警告), 未写入任何文件")
        return True, problems, warns, stats, None

    mpath, _moved = write_outputs(cases_dir, entries, sources, stats, problems, warns, a,
                                  manifest_extra=({"material": material} if material else None))
    say(f"[ingest] 已写入 {len(entries)} 个 OUTCAR -> {os.path.abspath(cases_dir)}")
    say(f"[ingest] manifest: {os.path.abspath(mpath)}")
    return True, problems, warns, stats, mpath


# ---------------------------------------------------------------------------
# 自检 (无 DFT)
# ---------------------------------------------------------------------------
def _selftest():
    checks = []

    def ck(name, cond, extra=""):
        checks.append(bool(cond))
        print(f"  [{'ok ' if cond else 'FAIL'}] {name}" + (f"  {extra}" if extra else ""))

    td = tempfile.mkdtemp(prefix="ingest_")
    try:
        proj = os.path.join(td, "Fe3Sn_pipeline")
        cases = os.path.join(proj, "02_supercell", "cases")
        runs = os.path.join(proj, "runs")
        orders = ["fm", "q0_120", "sq3_updntri", "sq3_alt"]
        for o in orders:
            os.makedirs(os.path.join(cases, o), exist_ok=True)
            os.makedirs(os.path.join(runs, o), exist_ok=True)
            with open(os.path.join(cases, o, "POSCAR"), "w", encoding="utf-8") as f:
                f.write("synth\n1.0\n9.3 0 0\n-4.65 8.05 0\n0 0 24\nFe Sn\n108 36\nDirect\n"
                        + "".join("0.1 0.2 0.3\n" for _ in range(144)))
        with open(os.path.join(cases, "magnetic_orders.json"), "w", encoding="utf-8") as f:
            json.dump({"selected_orders": orders, "n_mag": 108, "spin": 1.0,
                       "coefficients_per_cell": {
                           "fm": [-54, -54, -81], "q0_120": [27, 27, -81],
                           "sq3_updntri": [-6, -6, 27], "sq3_alt": [-6, 18, 27]}}, f)

        # 参考能量: 用 fm 行系数 + 一组合成 J 造出"自洽"的四个能量
        S2, J = 1.0, (10.0, -2.0, 0.5)
        c = {"fm": [-54, -54, -81], "q0_120": [27, 27, -81],
             "sq3_updntri": [-6, -6, 27], "sq3_alt": [-6, 18, 27]}
        E = {o: (sum(c[o][i] * J[i] * S2 for i in range(3))) / 1000.0 for o in orders}
        E = {o: round(v, 8) for o, v in E.items()}
        tots = {"fm": 150.0, "q0_120": 0.0, "sq3_updntri": 0.0, "sq3_alt": 0.0}
        for o in orders:
            synth_outcar(os.path.join(runs, o, "OUTCAR"), E[o] + 100.0, sigma0=E[o] + 100.0,
                         steps=1, accuracy=True, tot_moment=tots[o])

        # ---- 好数据: 四个序都过 ---------------------------------------------
        ok, probs, warns, stats, mp = ingest(runs, cases, orders=orders, min_kb=0)
        ck("合数据: 全部闸门通过", ok, f"problems={len(probs)} warns={len(warns)}")
        ck("合数据: 写出 manifest", mp is not None and os.path.isfile(mp))
        ck("合数据: 4 个 OUTCAR 就位",
           all(os.path.isfile(os.path.join(cases, o, "OUTCAR")) for o in orders))
        ck("合数据: 每个 OUTCAR 配 .sha256",
           all(os.path.isfile(os.path.join(cases, o, "OUTCAR.sha256")) for o in orders))
        ck("合数据: 能量离散度被记录", stats.get("spread_eV") is not None,
           f"spread={stats.get('spread_eV')}")
        ck("合数据: 报告文件写出", os.path.isfile(os.path.join(cases, "ingest_report.md")))

        # 能量必须能被 solve 口径读回来
        got = {o: DE.read_toten(os.path.join(cases, o, "OUTCAR")) for o in orders}
        ck("能量口径: 读回值与合成值一致 (1e-9)",
           max(abs(got[o] - E[o] - 100.0) for o in orders) < 1e-9)

        # ---- 坏数据 1: 未收敛 -------------------------------------------------
        synth_outcar(os.path.join(runs, "fm", "OUTCAR"), E["fm"] + 100.0, scf=False,
                     steps=1, accuracy=True)
        ok, probs, _w, _s, _m = ingest(runs, cases, orders=orders, dry_run=True,
                                       min_kb=0, force=True)
        ck("坏数据 SCF 未收敛: 拒绝", not ok and any(c == "scf" for _o, c, _msg in probs),
           f"报 {len(probs)} 条")

        # ---- 坏数据 2: 两个离子步 (松弛过的结构) ------------------------------
        synth_outcar(os.path.join(runs, "fm", "OUTCAR"), E["fm"] + 100.0, steps=2,
                     accuracy=True)
        ok, probs, _w, _s, _m = ingest(runs, cases, orders=orders, dry_run=True,
                                       min_kb=0, force=True)
        ck("坏数据 两个离子步: 拒绝", not ok and any(c == "ionic_steps" for _o, c, _msg in probs),
           f"报 {len(probs)} 条")

        # ---- 坏数据 3: 能量口径 (只有 energy(sigma->0)) ------------------------
        p = os.path.join(runs, "fm", "OUTCAR")
        with open(p, "w", encoding="utf-8") as f:
            f.write(" vasp.6.4.2\n ions per type = 108 36\n NKPTS = 6\n ENCUT = 520.0\n"
                    " LNONCOLLINEAR = T\n LSORBIT = T\n"
                    " General timing and accounting\n"
                    "  energy(sigma->0) =  -100.12345678\n")
        ok, probs, _w, _s, _m = ingest(runs, cases, orders=orders, dry_run=True,
                                       min_kb=0, force=True)
        ck("坏数据 缺 free energy TOTEN: 拒绝",
           not ok and any(c == "convention" for _o, c, _msg in probs), f"报 {len(probs)} 条")

        # ---- 坏数据 4: 旧 OUTCAR 充数 (闸门 9) --------------------------------
        synth_outcar(os.path.join(runs, "fm", "OUTCAR"), E["fm"] + 100.0, steps=1,
                     accuracy=True)
        # cases/fm/OUTCAR 是第一次 ingest 留下的 -> 与 runs 的逐位相同
        ok, probs, _w, _s, _m = ingest(runs, cases, orders=orders, dry_run=True, min_kb=0)
        ck("坏数据 旧 OUTCAR 充数: 拒绝 (闸门 9)",
           not ok and any(c == "not_shadowed" for _o, c, _msg in probs), f"报 {len(probs)} 条")
        ok2, _p2, _w2, _s2, _m2 = ingest(runs, cases, orders=orders, dry_run=True,
                                         min_kb=0, force=True)
        ck("--force 可以放行闸门 9", ok2)

        # ---- 坏数据 5: 缺一个序 ----------------------------------------------
        os.remove(os.path.join(runs, "sq3_alt", "OUTCAR"))
        try:
            ingest(runs, cases, orders=orders, dry_run=True, min_kb=0, force=True)
            ck("坏数据 缺一个序: 报错", False)
        except SystemExit as e:
            ck("坏数据 缺一个序: 报错", "sq3_alt" in str(e), str(e).splitlines()[0][:64])
        synth_outcar(os.path.join(runs, "sq3_alt", "OUTCAR"), E["sq3_alt"] + 100.0,
                     steps=1, accuracy=True, tot_moment=0.0)

        # ---- 坏数据 6: 磁序没生效 (与 fm 能量几乎相同) ------------------------
        synth_outcar(os.path.join(runs, "sq3_alt", "OUTCAR"), E["fm"] + 100.0, steps=1,
                     accuracy=True, tot_moment=0.0)
        ok, probs, _w, _s, _m = ingest(runs, cases, orders=orders, dry_run=True,
                                       min_kb=0, force=True, min_energy_gap=1e-4)
        ck("坏数据 磁序没生效 (gap≈0): 拒绝",
           not ok and any(c == "spread" for _o, c, _msg in probs), f"报 {len(probs)} 条")
        synth_outcar(os.path.join(runs, "sq3_alt", "OUTCAR"), E["sq3_alt"] + 100.0,
                     steps=1, accuracy=True, tot_moment=0.0)

        # ---- 坏数据 7: 结构不同 (原子数不一致) --------------------------------
        synth_outcar(os.path.join(runs, "fm", "OUTCAR"), E["fm"] + 100.0, steps=1,
                     accuracy=True, counts="54 18")
        ok, probs, _w, _s, _m = ingest(runs, cases, orders=orders, dry_run=True,
                                       min_kb=0, force=True)
        ck("坏数据 结构不同: 拒绝",
           not ok and any(c == "structure" for _o, c, _msg in probs), f"报 {len(probs)} 条")

        # ---- 坏数据 8: 空壳文件 ------------------------------------------------
        with open(os.path.join(runs, "fm", "OUTCAR"), "w", encoding="utf-8") as f:
            f.write("x")
        ok, probs, _w, _s, _m = ingest(runs, cases, orders=orders, dry_run=True,
                                       min_kb=200, force=True)
        ck("坏数据 空壳文件: 拒绝 (size)",
           not ok and any(c == "size" for _o, c, _msg in probs), f"报 {len(probs)} 条")

        # ---- --lenient 能把硬闸门降级 ----------------------------------------
        synth_outcar(os.path.join(runs, "fm", "OUTCAR"), E["fm"] + 100.0, steps=1,
                     accuracy=False)
        ok0, p0, _w, _s, _m = ingest(runs, cases, orders=orders, dry_run=True,
                                     min_kb=0, force=True)
        _ = ok0
        ck("accuracy 默认只是警告 (不算硬错)",
           not any(c == "accuracy" for _o, c, _msg in p0))
        synth_outcar(os.path.join(runs, "fm", "OUTCAR"), E["fm"] + 100.0, steps=2)
        ok1, p1, _w, _s, _m = ingest(runs, cases, orders=orders, dry_run=True, min_kb=0,
                                     force=True, lenient=["ionic_steps"])
        ck("--lenient ionic_steps: 该闸门降为警告",
           ok1 and not any(c == "ionic_steps" for _o, c, _msg in p1), f"problems={len(p1)}")
    finally:
        shutil.rmtree(td, ignore_errors=True)

    n_ok = sum(checks)
    print(f"\n--selftest: {'ALL PASS' if n_ok == len(checks) else 'FAIL'} ({n_ok}/{len(checks)} 项)")
    return 0 if n_ok == len(checks) else 1


# ---------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(
        description="DFT 结果回收闸门 (runs/ -> cases/); 全部校验通过才写入")
    ap.add_argument("--runs", help="交付出去的 runs/ 目录 (含 <order>/OUTCAR)")
    ap.add_argument("--cases", help="目标 02_supercell/cases")
    ap.add_argument("--orders", nargs="+", default=None,
                    help="磁序列表 (默认取 cases/magnetic_orders.json 的 selected_orders)")
    ap.add_argument("--dry-run", action="store_true", help="只校验不写")
    ap.add_argument("--force", action="store_true",
                    help="放行闸门 9 (允许与 cases 里已有 OUTCAR 逐位相同)")
    ap.add_argument("--lenient", nargs="*", default=None,
                    help="把这些硬检查降为警告: " + " ".join(n for n, h, _d in CHECKS if h))
    ap.add_argument("--min-kb", type=int, default=200, help="OUTCAR 最小体积 KB (默认 200)")
    ap.add_argument("--max-spread", type=float, default=2.0,
                    help="四态能量离散度上限 eV/胞 (默认 2.0)")
    ap.add_argument("--min-energy-gap", type=float, default=1e-4,
                    help="非 fm 序相对 fm 的最小能量差 eV (默认 1e-4); 低于此值判为磁序没生效")
    ap.add_argument("--material", default=None, help="材料名 (写进 manifest)")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    _argv = list(sys.argv[1:] if argv is None else argv)
    if "--selftest" in _argv:
        return _selftest()
    args = ap.parse_args(_argv)

    if not args.runs or not args.cases:
        raise SystemExit("需要 --runs 与 --cases (或用 --selftest)")
    for nm in args.lenient or []:
        if nm not in [n for n, _h, _d in CHECKS]:
            raise SystemExit(f"--lenient 不认识 {nm!r}; 可选: {[n for n, _h, _d in CHECKS]}")

    try:
        ok, problems, warns, stats, mp = ingest(
            args.runs, args.cases, orders=args.orders, max_spread=args.max_spread,
            min_energy_gap=args.min_energy_gap, dry_run=args.dry_run, force=args.force,
            lenient=args.lenient, min_kb=args.min_kb, quiet=args.quiet,
            material=args.material,
            progress=(lambda m: None) if args.quiet else (lambda m: print(m)))
    except SystemExit as e:
        print(f"[ingest] 拒绝: {e}", file=sys.stderr)
        return 2

    if warns and not args.quiet:
        print(f"\n[ingest] {len(warns)} 条警告:")
        for o, c, m in warns:
            print(f"  [warn] [{c}] {('(' + o + ') ') if o else ''}{m}")
    if not ok:
        print(f"\n[ingest] 拒绝写入 ({len(problems)} 条硬问题):", file=sys.stderr)
        for o, c, m in problems:
            print(f"  [FAIL] [{c}] {('(' + o + ') ') if o else ''}{m}", file=sys.stderr)
        print("\n  修法: 核对上面的每一条; 确认某条是误判才用 `--lenient <检查名>` 放行。",
              file=sys.stderr)
        return 1

    if not args.quiet:
        st = stats or {}
        print(f"\n[ingest] 全部闸门通过。能量离散度 = {st.get('spread_eV')} eV/胞")
        gaps = st.get("gap_vs_fm_eV") or {}
        for o in sorted(gaps):
            print(f"  {o:>12s}: E = {gaps[o] * 1000:+.3f} meV (相对 fm)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
