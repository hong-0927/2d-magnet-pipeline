#!/usr/bin/env python
"""assemble_dft_tri.py
三角晶格的非共线 VASP 输入组装。

P3.16 起本脚本是 **`assemble_dft.py` 的薄适配层**（NEW_TYPES_PLAN §106: 消灭复制粘贴的
builder 模式）。历史上它自带一份 3 分量 MAGMOM 模板和自己的 INCAR 文本，与
`assemble_dft.py` 是复制粘贴关系 —— 现在:

  * 角度/自旋生成规则  -> `mag_vectors.py`（单一来源）
  * INCAR / MAGMOM 渲染 -> `assemble_dft.py --noncollinear [--lsorbit]`
  * 本脚本只剩「三角 = 非共线 + 全体磁序同一设置（可比）+ 可选 SOC」这条策略差异

非共线 (含 120 度): `LNONCOLLINEAR = .TRUE.`, MAGMOM 每原子 3 分量。
共线态 (fm/stripy/zigzag) 也用同一非共线设置跑（自旋沿 z），保证整组可比。
弱铁磁倾斜/DMI 需要 `--lsorbit`（并会让 `assemble_dft.py` 检查磁序是否真需要 SOC）。

用法:
  python assemble_dft_tri.py --root <proj> --element Mn --nonmag I --U 4.0 --moment 5.0
  python assemble_dft_tri.py --root <proj> --element Mn --nonmag I --U 4.0 --moment 5.0 --lsorbit
  python assemble_dft_tri.py --selftest
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
TARGET = os.path.join(HERE, "assemble_dft.py")


def build_cmd(a):
    """把三角侧的参数翻译成 `assemble_dft.py` 的命令行（纯函数, 便于自检）。"""
    cmd = [sys.executable, TARGET, "--root", a.root, "--element", a.element,
           "--noncollinear", "--all-orders",
           "--U", str(a.U), "--moment", str(a.moment), "--encut", str(a.encut),
           "--cases", a.cases]
    if a.nonmag:
        cmd += ["--nonmag"] + list(a.nonmag)
    if a.lsorbit:
        cmd.append("--lsorbit")
    if a.saxis:
        cmd += ["--saxis"] + [str(x) for x in a.saxis]
    # P4.2: k 网格 / NPAR / ISYM 也可显式覆盖（缺省一律不传，保持 assemble_dft 的模板默认）
    if getattr(a, "kpoints", None):
        cmd += ["--kpoints"] + [str(x) for x in a.kpoints]
    if getattr(a, "npar", None) is not None:
        cmd += ["--npar", str(a.npar)]
    if getattr(a, "isym", None) is not None:
        cmd += ["--isym", str(a.isym)]
    return cmd


def main(argv=None):
    ap = argparse.ArgumentParser(description="三角晶格非共线 VASP 输入（转发给 assemble_dft.py）")
    ap.add_argument("--root", required=True)
    ap.add_argument("--element", required=True)
    ap.add_argument("--nonmag", nargs="+", default=[])
    ap.add_argument("--U", type=float, default=4.0)
    ap.add_argument("--moment", type=float, default=5.0)
    ap.add_argument("--encut", type=float, default=520)
    ap.add_argument("--kpoints", type=int, nargs=3, default=None,
                    help="显式 k 网格 (默认用 assemble_dft 的 3 3 1)")
    ap.add_argument("--npar", type=int, default=None, help="显式 NPAR (默认 4)")
    ap.add_argument("--isym", type=int, default=None, help="显式 ISYM (0 = 关对称性)")
    ap.add_argument("--cases", default="02_supercell/cases")
    ap.add_argument("--lsorbit", action="store_true", help="加 SOC (弱铁磁倾斜/DMI)")
    ap.add_argument("--saxis", nargs=3, type=float, default=None)
    ap.add_argument("--print-cmd", action="store_true", help="只打印转发命令")
    ap.add_argument("--selftest", action="store_true")
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--selftest" in argv:
        return _selftest()
    args = ap.parse_args(argv)

    cmd = build_cmd(args)
    if args.print_cmd:
        print(" ".join(cmd))
        return 0
    if not os.path.exists(TARGET):
        raise SystemExit(f"找不到 {TARGET}")
    return subprocess.call(cmd)


# ---------------------------------------------------------------------------
def _selftest():
    checks = []

    def ck(name, cond, extra=""):
        checks.append(bool(cond))
        print(f"  [{'ok ' if cond else 'FAIL'}] {name}" + (f"  {extra}" if extra else ""))

    class A:
        pass

    ck("assemble_dft.py 存在", os.path.exists(TARGET))

    a = A()
    a.root, a.element, a.nonmag, a.U, a.moment, a.encut = "/tmp/x", "Mn", ["I"], 4.0, 5.0, 520
    a.cases, a.lsorbit, a.saxis = "02_supercell/cases", False, None
    cmd = build_cmd(a)
    ck("转发命令含 --noncollinear 与 --all-orders（三角=非共线+全体可比）",
       "--noncollinear" in cmd and "--all-orders" in cmd)
    ck("转发命令带元素/非磁/参数", cmd[cmd.index("--element") + 1] == "Mn"
       and cmd[cmd.index("--nonmag") + 1:] == ["I"] and cmd[cmd.index("--U") + 1] == "4.0")
    ck("未给 --lsorbit 时命令里没有 --lsorbit", "--lsorbit" not in cmd)

    a.lsorbit, a.saxis = True, [0.0, 0.0, 1.0]
    cmd2 = build_cmd(a)
    ck("--lsorbit 透传, 且 --saxis 三个分量都在",
       "--lsorbit" in cmd2 and cmd2[cmd2.index("--saxis") + 1:] == ["0.0", "0.0", "1.0"])

    a.kpoints, a.npar, a.isym = [6, 6, 1], 8, 0
    cmd3 = build_cmd(a)
    ck("--kpoints/--npar/--isym 透传 (P4.2 参数化)",
       cmd3[cmd3.index("--kpoints") + 1:] == ["6", "6", "1", "--npar", "8", "--isym", "0"])
    a.kpoints, a.npar, a.isym = None, None, None
    ck("缺省时不传 k 网格/NPAR/ISYM (保持 assemble_dft 的模板默认)",
       "--kpoints" not in build_cmd(a) and "--npar" not in build_cmd(a)
       and "--isym" not in build_cmd(a))

    # 端到端: 造一个三角式 cases 目录, 真跑一次转发
    td = tempfile.mkdtemp(prefix="asmtri_")
    try:
        cases = os.path.join(td, "02_supercell", "cases")
        poscar = ("probe\n1.0\n4.17 0 0\n-2.085 3.611 0\n0 0 20\nMn I\n3 3\nDirect\n"
                  + "".join(f"{0.05 * i:.6f} {0.05 * i:.6f} 0.5\n" for i in range(6)))
        spin_120 = [[5.0, 0.0, 0.0], [-2.5, 4.330127, 0.0], [-2.5, -4.330127, 0.0]]
        for name, spins in (("fm", [[0.0, 0.0, 5.0]] * 3), ("120deg", spin_120)):
            d = os.path.join(cases, name)
            os.makedirs(d, exist_ok=True)
            open(os.path.join(d, "POSCAR"), "w", encoding="utf-8").write(poscar)
            body = " ".join(f"{x:.6f} {y:.6f} {z:.6f}" for x, y, z in spins)
            open(os.path.join(d, "magmom.txt"), "w", encoding="utf-8").write(
                f"MAGMOM   = {body} " + " ".join(["0.0 0.0 0.0"] * 3))
        json.dump({"selected_orders": ["fm"]}, open(os.path.join(cases, "magnetic_orders.json"), "w"))

        rc = main(["--root", td, "--element", "Mn", "--nonmag", "I",
                   "--U", "4.0", "--moment", "5.0"])
        runs = sorted(os.listdir(os.path.join(td, "runs")))
        inc = open(os.path.join(td, "runs", "120deg", "INCAR"), encoding="utf-8").read()
        ck("端到端: 退出码 0 且 03_dft/INCAR.static 写出",
           rc == 0 and os.path.exists(os.path.join(td, "03_dft", "INCAR.static")))
        ck("端到端: --all-orders 让 runs/ 含 ['120deg','fm']（不只 selected_orders）",
           runs == ["120deg", "fm"], f"runs={runs}")
        ck("端到端: INCAR 非共线 (有 LNONCOLLINEAR, 无 ISPIN)",
           "LNONCOLLINEAR = .TRUE." in inc and "ISPIN" not in inc)
        _ml = next(ln for ln in inc.splitlines() if ln.strip().upper().startswith("MAGMOM"))
        _nv = len(_ml.split("=", 1)[1].split())
        ck("端到端: MAGMOM 是 3 分量/原子 (6 原子 -> 18 个值)", _nv == 18, f"nvals={_nv}")
        ck("端到端: 未开 SOC 时 INCAR 里没有 LSORBIT", "LSORBIT" not in inc)
    finally:
        shutil.rmtree(td, ignore_errors=True)

    n_ok = sum(checks)
    print(f"\n--selftest: {'ALL PASS' if n_ok == len(checks) else 'FAIL'} ({n_ok}/{len(checks)} 项)")
    return 0 if n_ok == len(checks) else 1


if __name__ == "__main__":
    raise SystemExit(main() or 0)
