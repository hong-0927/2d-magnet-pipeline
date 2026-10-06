#!/usr/bin/env python
"""run_2d_magnet.py
二维磁性单层「结构 -> 磁耦合常数 -> Monte Carlo」通用流水线入口。

流程 (每一步调用 scripts/ 下对应脚本):
  1. fetch      : Materials Project 下载单层结构 (--formula / --mp-id), 或跳过 (--poscar)
  2. reorient   : 层法线对齐 z + 加真空
  3. supercell  : 2x2x1 超胞 + 四磁态 (FM/Neel/zigzag/stripy) + J1/J2/J3 系数
  4. dft        : 组装 INCAR/KPOINTS/POTCAR.spec + runs/<order>/
  5. mc         : POSCAR -> MC 模型 (含 J3) -> Heisenberg MC -> T_C/T_N
  6. render     : 四磁序俯视图
  7. workflow   : 写出 <outdir>/workflow.md

用法 (最简: 用文献 J 表 + MP 下载):
  python run_2d_magnet.py --formula MnPS3 --literature MnPS3

  # 或手动给全部参数 (新材料):
  python run_2d_magnet.py --formula VI3 --element V --nonmag I --S 1.0 --moment 2.0 --U 3.0 \
      --J1 -1.0 --J2 -0.1 --J3 0.0 --MAE 0.2

参数优先级: 显式 CLI 参数 > --literature 表 > 默认。
"""
import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(HERE, "scripts")
LIT = os.path.join(HERE, "literature_J.json")


def run(pyfile, args, cwd=None):
    cmd = [sys.executable, os.path.join(SCRIPTS, pyfile)] + args
    print("\n$", " ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd=cwd)
    if r.returncode != 0:
        print(f"[FAIL] {pyfile} 退出码 {r.returncode}", file=sys.stderr)
        sys.exit(r.returncode)


def load_literature(key):
    with open(LIT, encoding="utf-8") as f:
        tab = json.load(f)
    if key not in tab:
        raise SystemExit(f"literature_J.json 里没有 {key!r}; 可用: {list(tab.keys())}")
    return tab[key]


def pick_structure(struct_dir, formula, mp_id=None):
    """从 candidates.csv 选最稳定 (hull 最小) 的单层 POSCAR。"""
    import csv
    csv_path = os.path.join(os.path.dirname(struct_dir), "candidates.csv")
    rows = list(csv.DictReader(open(csv_path, encoding="utf-8-sig")))
    if not rows:
        raise SystemExit("candidates.csv 为空")
    if mp_id:
        rows = [r for r in rows if r["material_id"] == mp_id] or rows
    rows = sorted(rows, key=lambda r: (r.get("hull") or 9.0))
    best = rows[0]
    print(f"[pick] 选 {best['material_id']} ({best['formula']}) hull={best.get('hull')}")
    if len(rows) > 1:
        print(f"[pick] 共 {len(rows)} 个候选, 其余: "
              + ", ".join(r['material_id'] for r in rows[1:])
              + " (可用 --mp-id 指定)")
    return os.path.join(struct_dir, best["poscar"])


def detect_lattice(poscar_path, element, cutoff=3.3):
    """读 reorient 后 POSCAR, 数磁性原子/胞 + 簇数, 返回 (lattice, n_clusters, n_mag)。

    lattice ∈ {"honeycomb", "triangular", "cluster", "other"}:
      cluster    = 磁性原子成簇 (簇数 < 原子数), 簇 = 有效磁位
      triangular = 1 磁性原子/胞
      honeycomb  = 2 磁性原子/胞 (不分簇)
    """
    lines = [l for l in open(poscar_path, encoding="utf-8").read().splitlines() if l.strip()]
    lat = np.array([[float(x) for x in lines[i].split()] for i in range(2, 5)])
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    idx = 7
    if lines[idx].strip().lower()[:1] == "s":
        idx += 1
    mode = lines[idx].strip().lower()
    idx += 1
    atoms = []
    for sp, n in zip(species, counts):
        atoms += [sp] * n
    coords = np.array([[float(x) for x in lines[idx + i].split()[:3]]
                       for i in range(len(atoms))])
    if mode[:1] in "ck":
        coords = coords @ np.linalg.inv(lat)
    coords -= np.floor(coords)
    mag = [i for i, a in enumerate(atoms) if a == element]
    n_mag = len(mag)
    if n_mag == 0:
        raise SystemExit(f"POSCAR 里找不到磁性元素 {element}")

    parent = list(range(n_mag))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(n_mag):
        for j in range(i + 1, n_mag):
            d = coords[mag[j]] - coords[mag[i]]
            d -= np.floor(d + 0.5)
            if np.linalg.norm(d @ lat) < cutoff:
                ra, rb = find(i), find(j)
                if ra != rb:
                    parent[rb] = ra
    n_clusters = len({find(i) for i in range(n_mag)})

    if n_clusters < n_mag:
        return "cluster", n_clusters, n_mag
    if n_mag == 1:
        return "triangular", 1, 1
    if n_mag == 2:
        return "honeycomb", 2, 2
    return "other", n_clusters, n_mag


def main():
    ap = argparse.ArgumentParser(description="二维磁体 结构->J->MC 通用流水线")
    # 结构来源
    ap.add_argument("--formula", help="化学式, 如 MnPS3 (触发 MP 下载)")
    ap.add_argument("--mp-id", help="指定 MP material_id (覆盖自动挑选)")
    ap.add_argument("--poscar", help="已有单层 POSCAR (跳过下载)")
    ap.add_argument("--mp-key", default=None, help="MP_API_KEY (默认读环境变量)")
    # 材料/方法参数
    ap.add_argument("--element", help="磁性元素, 如 Mn")
    ap.add_argument("--nonmag", nargs="+", help="非磁元素, 如 P S")
    ap.add_argument("--S", type=float, help="自旋大小 (Mn2+=2.5, Cr3+=1.5)")
    ap.add_argument("--moment", type=float, help="磁矩初猜 muB (Mn2+=5.0, Cr3+=3.0)")
    ap.add_argument("--U", type=float, help="DFT+U 加在磁性元素 d 轨道 (eV)")
    # J 值
    ap.add_argument("--literature", help="literature_J.json 里的键, 如 MnPS3")
    ap.add_argument("--J1", type=float, default=None)
    ap.add_argument("--J2", type=float, default=None)
    ap.add_argument("--J3", type=float, default=None)
    ap.add_argument("--MAE", type=float, default=None)
    # 输出 / 开关
    ap.add_argument("--outdir", default=None, help="项目目录 (默认 ./<formula>_pipeline)")
    ap.add_argument("--vacuum", type=float, default=15.0)
    ap.add_argument("--supercell", type=int, nargs=3, default=[2, 2, 1])
    ap.add_argument("--lattice", choices=["auto", "honeycomb", "triangular", "cluster"],
                    default="auto", help="磁性格点类型 (默认自动判别)")
    ap.add_argument("--cluster-cutoff", type=float, default=3.3,
                    help="簇识别成键阈值 A (簇磁体用)")
    ap.add_argument("--mc-size", type=int, nargs=2, default=[30, 30])
    ap.add_argument("--mc-tmin", type=float, default=5.0)
    ap.add_argument("--mc-tmax", type=float, default=150.0)
    ap.add_argument("--mc-nt", type=int, default=30)
    ap.add_argument("--skip-dft", action="store_true")
    ap.add_argument("--skip-mc", action="store_true")
    args = ap.parse_args()

    # ---- 解析参数 (文献表 < 显式 CLI) ----
    lit = None
    if args.literature:
        lit = load_literature(args.literature)
        print(f"[lit] 使用文献值 ({args.literature}): {lit.get('source')}")

    def pick(a, b, d=None):
        return a if a is not None else (lit.get(b) if lit else d)

    formula = args.formula or (lit.get("formula") if lit else None)
    element = pick(args.element, "element")
    nonmag = args.nonmag or (lit.get("nonmag") if lit else [])
    S = pick(args.S, "S")
    moment = pick(args.moment, "moment")
    U = pick(args.U, "U")
    J1 = pick(args.J1, "J1")
    J2 = pick(args.J2, "J2", 0.0)
    J3 = pick(args.J3, "J3", 0.0)
    MAE = pick(args.MAE, "MAE", 0.0)

    missing = [k for k, v in [("--element", element), ("--S", S), ("--moment", moment),
                              ("--U", U), ("--J1", J1)] if v is None]
    if missing:
        raise SystemExit(f"缺少必要参数 {missing}; 用 --literature <key> 或显式给出")

    if not formula and not args.poscar and not args.mp_id:
        raise SystemExit("需要 --formula / --mp-id / --poscar 之一")

    if args.outdir:
        outdir = os.path.abspath(args.outdir)
    else:
        outdir = os.path.abspath((formula or "material") + "_pipeline")
    os.makedirs(outdir, exist_ok=True)
    print(f"[run] 项目目录: {outdir}")

    # ---- 1. fetch ----
    os.makedirs(os.path.join(outdir, "01_structure"), exist_ok=True)
    if args.poscar:
        src_poscar = os.path.abspath(args.poscar)
        print(f"[fetch] 跳过下载, 用 {src_poscar}")
    else:
        if args.mp_key:
            os.environ["MP_API_KEY"] = args.mp_key
        if "MP_API_KEY" not in os.environ:
            raise SystemExit("未设置 MP_API_KEY; 用 --mp-key 或环境变量")
        fetch_dir = os.path.join(outdir, "_fetch_tmp")
        os.makedirs(fetch_dir, exist_ok=True)
        # mp_fetch_2d.py 固定写脚本同目录 structures/ + candidates.csv; 复制到临时目录运行
        shutil.copy(os.path.join(SCRIPTS, "mp_fetch_2d.py"), fetch_dir)
        fa = [os.path.join(fetch_dir, "mp_fetch_2d.py")]
        if args.mp_id:
            fa += ["fetch", "--ids", args.mp_id]
        else:
            fa += ["fetch", "--formula", formula]
        print("\n$", " ".join([sys.executable] + fa), flush=True)
        r = subprocess.run([sys.executable] + fa, cwd=fetch_dir)
        if r.returncode != 0:
            sys.exit(r.returncode)
        struct_dir = os.path.join(fetch_dir, "structures")
        src_poscar = pick_structure(struct_dir, formula, args.mp_id)
    shutil.copy(src_poscar, os.path.join(outdir, "01_structure", "POSCAR.mp"))

    # ---- 2. reorient ----
    run("reorient.py", [os.path.join(outdir, "01_structure", "POSCAR.mp"),
                        "--vacuum", str(args.vacuum),
                        "-o", os.path.join(outdir, "01_structure", "POSCAR")], cwd=outdir)
    poscar_re = os.path.join(outdir, "01_structure", "POSCAR")

    # ---- 2.5 检测磁性格点类型 ----
    lattice, n_clusters, n_mag = detect_lattice(poscar_re, element, args.cluster_cutoff)
    if args.lattice != "auto":
        lattice = args.lattice
    print(f"[lattice] {lattice} (磁性原子 {n_mag}/胞, 簇 {n_clusters})")

    mc_poscar = poscar_re
    mc_order = "auto"
    do_dft = False
    do_render = False

    if lattice == "cluster":
        # 簇磁体: 磁性簇 -> 有效位, 三角/非二分晶格用比热 Cv 定 T_N
        poscar_clusters = os.path.join(outdir, "01_structure", "POSCAR.clusters.vasp")
        run("cluster_sites.py", [poscar_re, "--element", element,
                                 "--cutoff", str(args.cluster_cutoff),
                                 "-o", poscar_clusters], cwd=outdir)
        mc_poscar = poscar_clusters
        mc_order = "cv"
    elif lattice == "triangular":
        # 三角单离子: 1 位/胞, 用比热 Cv
        mc_order = "cv"
    elif lattice == "honeycomb":
        # 紧凑胞标准化 + 四磁态 + 系数
        run("standardize_cell.py", [poscar_re, "--element", element,
                                    "-o", os.path.join(outdir, "01_structure", "POSCAR.std")],
            cwd=outdir)
        poscar_std = os.path.join(outdir, "01_structure", "POSCAR.std")
        run("build_magnetic_cases3.py",
            [poscar_std, "--supercell"] + [str(x) for x in args.supercell] +
            ["--outdir", os.path.join(outdir, "02_supercell", "cases"),
             "--element", element, "--nonmag"] + nonmag +
            ["--moment", str(moment), "--spin", str(S)], cwd=outdir)
        mc_poscar = poscar_std
        mc_order = "auto"
        do_dft = True
        do_render = True
    else:
        raise SystemExit(f"无法处理晶格类型 {lattice!r} ({n_mag} 磁性原子/胞, {n_clusters} 簇); "
                         f"用 --lattice 指定或路由到 comp-chem-workflow 人工处理")

    # ---- DFT 输入 (仅蜂窝有完整四磁态反解) ----
    if do_dft and not args.skip_dft:
        run("assemble_dft.py",
            ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
            ["--species"] + ([element] + nonmag) +
            ["--U", str(U), "--moment", str(moment)], cwd=outdir)

    # ---- MC ----
    if not args.skip_mc:
        os.makedirs(os.path.join(outdir, "05_mc"), exist_ok=True)
        model = os.path.join(outdir, "05_mc", "mc_model.json")
        run("poscar_to_mc3.py",
            [mc_poscar,
             "--element", element, "--S", str(S),
             "--J1", str(J1), "--J2", str(J2), "--J3", str(J3), "--MAE", str(MAE),
             "-o", model], cwd=outdir)
        run("mc_2d.py", [model, "--order", mc_order,
                         "--size"] + [str(x) for x in args.mc_size] +
            ["--tmin", str(args.mc_tmin), "--tmax", str(args.mc_tmax), "--nt", str(args.mc_nt),
             "--cool", "--out", os.path.join(outdir, "05_mc", "mc")], cwd=outdir)

    # ---- render 超胞图 (仅蜂窝) ----
    if do_render:
        run("render_supercell.py",
            ["--root", outdir, "--element", element, "--nonmag"] + nonmag +
            ["--label", formula or element,
             "-o", os.path.join(outdir, "figures", "supercell_orders.png")],
            cwd=outdir)

    # ---- 7. workflow.md ----
    write_workflow(outdir, formula, element, nonmag, S, moment, U, J1, J2, J3, MAE,
                   args, lit, lattice)

    print("\n[run] 完成。项目目录:", outdir)
    print(f"      workflow.md, 01_structure/, 02_supercell/, 03_dft/, 05_mc/, figures/")


def write_workflow(outdir, formula, element, nonmag, S, moment, U, J1, J2, J3, MAE,
                   args, lit, lattice):
    gs = "Néel 反铁磁" if J1 < 0 else "铁磁"
    src = lit.get("source") if lit else "手动给定 (未核实文献)"
    if lattice == "honeycomb":
        model_desc = ("蜂窝 Heisenberg 模型 H = -J1 Σ_<ij> S·S - J2 Σ_<<ij>> S·S "
                      "- J3 Σ_<<<ij>>> S·S (J>0 铁磁)")
        mc_desc = "序参量 " + ("交错(反铁磁)" if J1 < 0 else "总磁化(铁磁)")
        struct_desc = ("- 超胞/磁态: `02_supercell/cases/{{fm,neel,zigzag,stripy}}` "
                       "+ `magnetic_orders.json`\n- DFT 输入: `03_dft/` + `runs/<order>/` (未运行, 待 POTCAR)")
    elif lattice == "cluster":
        model_desc = ("簇磁体: 磁性原子成簇 (Nb3 等), 每簇 = 1 个有效自旋, 簇中心成"
                      "三角/非二分晶格; 低能可能是单重态, 经典 Heisenberg 为近似")
        mc_desc = "序参量 比热 Cv (受挫/非二分晶格)"
        struct_desc = ("- 簇识别: `01_structure/POSCAR.clusters.vasp` (磁性簇 → 有效位)\n"
                       "- 注: DFT 反解 J 的三角/簇磁序枚举未实现, J 需文献/给定")
    else:  # triangular
        model_desc = "三角晶格单离子 Heisenberg 模型 (1 磁性原子/胞)"
        mc_desc = "序参量 比热 Cv (受挫/非二分晶格)"
        struct_desc = ("- 注: DFT 反解 J 的三角磁序枚举未实现, J 需文献/给定")

    lines = f"""# {formula} 单层 磁耦合 (J) → Monte Carlo 通用流水线

## 模型与方法

| 项 | 值 |
|---|---|
| 体系 | {formula} 单层, 磁性格点类型 = **{lattice}** |
| 磁性元素 / 非磁 | {element} / {', '.join(nonmag) if nonmag else '-'} |
| 自旋 S / 磁矩 | S={S} / {moment} μB |
| 泛函 | PBE + DFT+U (U={U} eV 加 {element} d) |
| 基态 | {gs} |
| 交换常数 (meV) | J1={J1}, J2={J2}, J3={J3}, MAE={MAE} |

模型: {model_desc}

J 值来源: {src}

## 结果

- 结构: `01_structure/POSCAR` (reorient 后, 真空 {args.vacuum} Å)
{struct_desc}
- MC: `05_mc/mc_mt.csv/.png`, {mc_desc}

## 局限 (通用)

1. 磁性格点类型 {lattice}; 非蜂窝/簇磁体的 DFT 反解 J 未实现, J 用文献/给定。
2. J 取自文献/给定; 未跑 DFT 前非本项目计算值。
3. 纯 Heisenberg MAE 默认为 0 → 二维无真正相变 (Mermin-Wagner), 临界温度为有限尺寸估计。
4. 簇磁体可能是单重态 (非磁有序), 经典 Heisenberg-MC 结果需谨慎解读。
5. mp_fetch_2d.py 的 exfoliate 对菱面体层状/已存成紧胞的材料有已知缺口。
6. k 网格 3x3x1、U 值为可调假设。
"""
    with open(os.path.join(outdir, "workflow.md"), "w", encoding="utf-8") as f:
        f.write(lines)
    print("[run] 写出 workflow.md")


if __name__ == "__main__":
    main()
