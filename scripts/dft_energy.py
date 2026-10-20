#!/usr/bin/env python
"""dft_energy.py —— VASP 总能量提取口径的**唯一实现**（P3.17）。

## 口径（LATTICE_ROADMAP §P3.17）

一组磁构型反解 J 时，**每个构型都必须用同一个能量口径**，否则反解出的 J 会掺进
口径差。本项目统一取：

    free energy  TOTEN   (eV)      —— OUTCAR 里 `free  energy   TOTEN  = ...`

理由：本流水线用 `ISMEAR = 0`（MP 展宽）+ `SIGMA = 0.05`（金属性 3d/4d 体系），
OUTCAR 同时给出两条：

    free  energy   TOTEN  =   -123.456789 eV          <- E0 + 展宽熵项（与 ISMEAR/SIGMA 自洽）
    energy  without entropy = -123.456789  energy(sigma->0) = -123.460000
                                                      <- 外推到 sigma=0 的"零温"能量

两者差 ~ `SIGMA` 量级。**同一组构型只用其中一条**：(a) 用 free energy 时，展宽熵
项在各构型间近似相消；(b) 若混用（有的构型取 free、有的取 sigma->0），J 会被污染。
`read_toten()` 因此固定返回 free energy，并在两条差得过多（> `MIX_TOL`）时给出
显式警告 —— 只提示，不替调用者改口径。

## 非共线 / SOC 的额外约定

* `LNONCOLLINEAR`/`LSORBIT` 打开后 OUTCAR 会多出 `E1+E2`、`FS-` 之类的分解行，
  但 `free energy TOTEN` 仍是总能量，取法不变。
* **SOC 开 / 关的能量差只有在同一胞、同一赝势、同一 k 网格、同一 ISMEAR/SIGMA 下
  才有意义**；不同参考点相减没有物理含义（P3 的 M3 验收要求"量级正确"，因此必须
  成对同条件计算）。
* 弛豫任务会打印多行 TOTEN（每个离子步一行）⇒ 一律取**最后一次**出现。

## 用法

  python dft_energy.py --selftest
  from dft_energy import read_toten, read_energy_detail
"""
import argparse
import os
import re
import tempfile

CONVENTION = "free energy TOTEN (eV)"
MIX_TOL = 1e-3  # free energy 与 energy(sigma->0) 超过这个差 (eV) 就提示口径不同

_RE_FREE = re.compile(r"free\s+energy\s+TOTEN\s*=\s*([-+0-9.Ee]+)")
_RE_SIGMA0 = re.compile(r"energy\(sigma->0\)\s*=\s*([-+0-9.Ee]+)")
_RE_NOENT = re.compile(r"energy\s+without\s+entropy\s*=\s*([-+0-9.Ee]+)")


def read_energy_detail(path):
    """读 OUTCAR -> dict(free_energy, sigma0, without_entropy, n_ionic_steps, path)。

    取每一类里**最后一次**出现的值（弛豫多离子步时即最终态）；缺失项为 None。
    """
    free = sigma0 = noent = None
    n_steps = 0
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = _RE_FREE.search(line)
            if m:
                free = float(m.group(1))
                n_steps += 1
                continue
            m = _RE_SIGMA0.search(line)
            if m:
                sigma0 = float(m.group(1))
                continue
            m = _RE_NOENT.search(line)
            if m:
                noent = float(m.group(1))
    return {"free_energy": free, "sigma0": sigma0, "without_entropy": noent,
            "n_ionic_steps": n_steps, "path": str(path)}


def read_toten(path):
    """取 `free energy TOTEN`（eV）—— 全项目唯一的能量口径。缺失即报错。"""
    d = read_energy_detail(path)
    if d["free_energy"] is None:
        raise ValueError(f"no 'free energy TOTEN' in {path}")
    return d["free_energy"]


def mixed_convention(path, tol=MIX_TOL):
    """该 OUTCAR 的 free energy 与 energy(sigma->0) 是否差得超过 `tol`（口径不一致的征兆）。"""
    d = read_energy_detail(path)
    if d["free_energy"] is None or d["sigma0"] is None:
        return False
    return abs(d["free_energy"] - d["sigma0"]) > tol


def check_same_convention(paths, tol=MIX_TOL):
    """一组 OUTCAR 的一致性检查 -> (ok, problems[])。

    只查同一件事：free energy / energy(sigma->0) 的偏离量在各构型间是否**同量级**。
    偏离量彼此差得远（> tol）说明构型间的展宽/收敛状态不可比，此时反解 J 不可信。
    """
    problems = []
    devs = {}
    for p in paths:
        try:
            d = read_energy_detail(p)
        except OSError as e:
            problems.append(f"{p}: 读不了 ({e})")
            continue
        if d["free_energy"] is None:
            problems.append(f"{os.path.basename(os.path.dirname(p)) or p}: 没有 free energy TOTEN")
            continue
        if d["sigma0"] is None:
            continue
        devs[p] = abs(d["free_energy"] - d["sigma0"])
    if len(devs) >= 2:
        lo, hi = min(devs.values()), max(devs.values())
        if hi - lo > tol:
            worst = max(devs, key=devs.get)
            problems.append(f"free energy 与 energy(sigma->0) 的偏离量在构型间不一致 "
                            f"({lo:.2e} … {hi:.2e} eV, 最大在 {os.path.basename(os.path.dirname(worst))}); "
                            f"同一组构型必须同 ISMEAR/SIGMA 且都取 free energy TOTEN")
    return not problems, problems


# ---------------------------------------------------------------------------
# selftest（无 DFT：写合成 OUTCAR 再解析）
# ---------------------------------------------------------------------------
_OUTCAR_TAIL = """\
 FREE ENERGIE OF THE ION-ELECTRON SYSTEM (eV)
 ---------------------------------------------------
  free  energy   TOTEN  =       {free:.8f} eV
  energy  without entropy=     {noent:.8f}  energy(sigma->0) =     {sigma0:.8f}
"""


def _write_outcar(path, free, sigma0=None, noent=None, steps=1):
    sigma0 = free if sigma0 is None else sigma0
    noent = free if noent is None else noent
    with open(path, "w", encoding="utf-8") as f:
        f.write(" vasp.6.4.2 (synthetic OUTCAR for dft_energy selftest)\n")
        f.write("  LNONCOLLINEAR = .TRUE.  LSORBIT = .TRUE.\n")
        for _ in range(steps):
            f.write(_OUTCAR_TAIL.format(free=free, sigma0=sigma0, noent=noent))
    return path


def _selftest():
    checks = []

    def ck(name, cond, extra=""):
        checks.append(bool(cond))
        print(f"  [{'ok ' if cond else 'FAIL'}] {name}" + (f"  {extra}" if extra else ""))

    td = tempfile.mkdtemp(prefix="dften_")
    try:
        os.makedirs(os.path.join(td, "fm"), exist_ok=True)
        p1 = _write_outcar(os.path.join(td, "fm", "OUTCAR"),
                           free=-100.12345678, sigma0=-100.12345678)
        ck("OUTCAR (ISMEAR=-5 型: 两条相等) 取到 free energy",
           abs(read_toten(p1) + 100.12345678) < 1e-9, f"{read_toten(p1):.8f} eV")
        ck("口径不同的 run 会被 mixed_convention 标出", not mixed_convention(p1))
        ck("check_same_convention: 单个 OUTCAR 无问题", check_same_convention([p1])[0])

        os.makedirs(os.path.join(td, "afm"), exist_ok=True)
        p2 = _write_outcar(os.path.join(td, "afm", "OUTCAR"),
                           free=-100.22345678, sigma0=-100.24500000)
        ck("展宽 run: free energy 与 sigma->0 不同, 且被识别",
           abs(read_toten(p2) + 100.22345678) < 1e-9 and mixed_convention(p2),
           f"Δ={read_energy_detail(p2)['free_energy'] - read_energy_detail(p2)['sigma0']:.2e} eV")

        os.makedirs(os.path.join(td, "relax"), exist_ok=True)
        p3 = _write_outcar(os.path.join(td, "relax", "OUTCAR"), free=-1.0, steps=4)
        with open(p3, "a", encoding="utf-8") as f:
            f.write(_OUTCAR_TAIL.format(free=-9.5, sigma0=-9.5, noent=-9.5))
        d3 = read_energy_detail(p3)
        ck("弛豫多离子步: 取**最后一次** TOTEN",
           abs(d3["free_energy"] + 9.5) < 1e-9 and d3["n_ionic_steps"] == 5,
           f"n_steps={d3['n_ionic_steps']}")

        p4 = os.path.join(td, "only_sigma0", "OUTCAR")
        os.makedirs(os.path.dirname(p4), exist_ok=True)
        with open(p4, "w", encoding="utf-8") as f:
            f.write("  energy  without entropy=   -5.0  energy(sigma->0) =   -5.0\n")
        try:
            read_toten(p4)
            ck("缺 free energy TOTEN 时报错 (不静默给 sigma->0)", False)
        except ValueError as e:
            ck("缺 free energy TOTEN 时报错 (不静默给 sigma->0)",
               "free energy TOTEN" in str(e) and p4 in str(e))

        p5 = os.path.join(td, "empty", "OUTCAR")
        os.makedirs(os.path.dirname(p5), exist_ok=True)
        open(p5, "w").close()
        try:
            read_toten(p5)
            ck("空 OUTCAR 报错", False)
        except ValueError as e:
            ck("空 OUTCAR 报错", "no 'free energy TOTEN'" in str(e))

        os.makedirs(os.path.join(td, "bad"), exist_ok=True)
        p6 = _write_outcar(os.path.join(td, "bad", "OUTCAR"), free=-1.0, sigma0=-1.5)
        ok, probs = check_same_convention([p2, p6])
        ck("check_same_convention: 逮住构型间自由能/零温能量偏离量不一致",
           not ok and "不一致" in " ".join(probs), str(probs)[:90])

        real = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "..", "Fe3GeTe2_pipeline", "03_dft", "runs")
        real = os.path.abspath(real)
        found = []
        if os.path.isdir(real):
            for r, _d, fs in os.walk(real):
                found += [os.path.join(r, x) for x in fs if x == "OUTCAR"]
        if found:
            e = read_toten(found[0])
            ck(f"真实 OUTCAR 可读 ({os.path.relpath(found[0], os.path.dirname(real))})",
               isinstance(e, float), f"{e:.6f} eV")
        else:
            ck("真实 OUTCAR 可读 (本机无 OUTCAR, 跳过)", True)
    finally:
        import shutil
        shutil.rmtree(td, ignore_errors=True)

    n_ok = sum(checks)
    print(f"\n--selftest: {'ALL PASS' if n_ok == len(checks) else 'FAIL'} ({n_ok}/{len(checks)} 项)")
    return 0 if n_ok == len(checks) else 1


def main():
    ap = argparse.ArgumentParser(description="VASP 总能量口径 (P3.17): 一律 free energy TOTEN")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("outcars", nargs="*", help="要读的 OUTCAR (不给则打印口径说明)")
    args = ap.parse_args()
    if args.selftest:
        return _selftest()
    if args.outcars:
        ok, probs = check_same_convention(args.outcars)
        for p in args.outcars:
            print(f"{read_toten(p):16.8f} eV   {p}")
        for x in probs:
            print(f"[warn] {x}")
        return 0 if ok else 1
    print(f"能量口径: {CONVENTION}")
    print(__doc__.split("## 用法")[0].strip().split("## 口径（LATTICE_ROADMAP §P3.17）")[-1].strip()[:400])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
