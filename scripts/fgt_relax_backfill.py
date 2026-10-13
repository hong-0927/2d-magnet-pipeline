#!/usr/bin/env python
"""fgt_relax_backfill.py — 把 FM 弛豫得到的 CONTCAR 回填成 5 个磁态的 POSCAR

流程:
    1) 先只跑 runs/fm/ (POSCAR + INCAR.relax) 得到 CONTCAR;
    2) 本脚本把 CONTCAR 写成 cases/<order>/POSCAR 与 runs/<order>/POSCAR;
    3) 5 个静态计算于是都在**同一个弛豫后结构**上做 —— 反解要求各磁态除自旋向外完全可比,
       哪怕 0.01 A 的结构差都会混进能量差里, 直接污染 J。

    INCAR 里的逐原子 MAGMOM 不动 (它只跟自旋模式有关, 与坐标无关)。

用法:
    python fgt_relax_backfill.py <cases_dir> [--runs <runs_dir>] [--contcar <路径>] [--dry-run]
    (第一个参数也可以直接给流水线 outdir, 会自动下钻到 02_supercell/cases)
默认 CONTCAR = <runs>/fm/CONTCAR

校验 (任一不过就拒绝写):
    - 5 个 cases/<order>/POSCAR 必须原本一致 (只允许首行注释不同);
    - CONTCAR 的物种顺序与原子数必须与它们一致;
    - CONTCAR 的坐标行数必须够。
同时会打印**最大坐标位移**; 位移 ~0 时警告 —— 静态自洽 (IBRION=-1) 也会写 CONTCAR,
但坐标与 POSCAR 逐位相同, 那种情况下"回填"等于什么都没做, 整个反解前提就废了。
同时在 cases/<order>/ 留一份 POSCAR.pristine (原始未弛豫), 可反复回填/对照;
再次运行时位移检查以 POSCAR.pristine 为基准 (所以幂等重跑不会误报)。
"""
import argparse
import os
import shutil


def read_poscar(path):
    lines = open(path, encoding="utf-8", errors="ignore").read().splitlines()
    if len(lines) < 8:
        raise SystemExit(f"{path} 太短, 不像 POSCAR")
    species = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    coord = lines[7].strip()
    if not (coord.lower().startswith("d") or coord.lower().startswith("c")):
        raise SystemExit(f"{path} 第 8 行不是 Direct/Cartesian: {coord!r}")
    return lines, species, counts


def lattice_of(lines):
    """返回 3x3 晶格矩阵 (已乘缩放系数)。"""
    s = float(lines[1].split()[0])
    return [[float(x) * s for x in lines[2 + i].split()[:3]] for i in range(3)]


def coords_cart(lines, lat, n):
    """前 n 个原子的 Cartesian 坐标 (A)。'Direct' 按晶格换算, 'Cartesian' 直接取 (含缩放)。"""
    s = float(lines[1].split()[0])
    direct = lines[7].strip().lower().startswith("d")
    out = []
    for i in range(8, 8 + n):
        x = [float(t) for t in lines[i].split()[:3]]
        if direct:
            out.append([sum(x[k] * lat[k][j] for k in range(3)) for j in range(3)])
        else:
            out.append([t * s for t in x])
    return out


def resolve_cases(p):
    """接受 cases 目录, 也接受流水线 outdir (自动下钻 02_supercell/cases)。"""
    p = os.path.abspath(p)
    if os.path.exists(os.path.join(p, "magnetic_orders.json")):
        return p
    c = os.path.join(p, "02_supercell", "cases")
    if os.path.exists(os.path.join(c, "magnetic_orders.json")):
        return c
    raise SystemExit(
        f"{p} 里既没有 magnetic_orders.json, 也没有 02_supercell/cases。\n"
        f"  请指向 cases 目录 (或流水线 outdir)。")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cases", help="02_supercell/cases 目录 (也接受流水线 outdir, 会自动下钻)")
    ap.add_argument("--runs", default=None, help="runs 目录 (默认 cases 的同级 runs)")
    ap.add_argument("--contcar", default=None, help="弛豫得到的 CONTCAR (默认 <runs>/fm/CONTCAR)")
    ap.add_argument("--dry-run", action="store_true", help="只检查, 不写文件")
    args = ap.parse_args()

    cases = resolve_cases(args.cases)
    runs = args.runs or os.path.join(os.path.dirname(cases), "runs")
    cont = args.contcar or os.path.join(runs, "fm", "CONTCAR")
    if not os.path.exists(cont):
        raise SystemExit(
            f"找不到 CONTCAR: {cont}\n"
            f"  先跑 runs/fm/ (POSCAR + INCAR.relax), 把 VASP 产出的 CONTCAR 放到那里,\n"
            f"  或用 --contcar 指定别处的路径。")

    orders = sorted(d for d in os.listdir(cases) if os.path.isdir(os.path.join(cases, d)))
    if not orders:
        raise SystemExit(f"{cases} 下没有磁态子目录")

    # 1) 各磁态 POSCAR 必须原本一致 (只允许首行不同)
    ref = None
    for o in orders:
        p = os.path.join(cases, o, "POSCAR")
        if not os.path.exists(p):
            raise SystemExit(f"缺 {p}")
        body = open(p, encoding="utf-8", errors="ignore").read().splitlines()[1:]
        if ref is None:
            ref, ref_o = body, o
        elif body != ref:
            raise SystemExit(
                f"{o}/POSCAR 与 {ref_o}/POSCAR 不一致 (除首行外)。\n"
                f"  5 个磁态必须共用同一坐标, 否则能量差里混进了结构差, 反解无意义。")
    # 位移/晶格检查的基准: 优先用原始未弛豫坐标 (pristine), 这样幂等重跑不会把"已回填"
    # 误判成"没弛豫"。
    base_path = os.path.join(cases, ref_o, "POSCAR")
    base_use = base_path + ".pristine" if os.path.exists(base_path + ".pristine") else base_path
    base_lines = open(base_use, encoding="utf-8", errors="ignore").read().splitlines()

    # 2) CONTCAR 与 case POSCAR 的物种/原子数一致性
    clines, cspe, ccnt = read_poscar(cont)
    _, rspe, rcnt = read_poscar(base_use)
    if cspe != rspe or ccnt != rcnt:
        raise SystemExit(
            f"CONTCAR 的物种/原子数与超胞不符:\n"
            f"  CONTCAR : {cspe} {ccnt}  (共 {sum(ccnt)})\n"
            f"  超胞    : {rspe} {rcnt}  (共 {sum(rcnt)})\n"
            f"  CONTCAR 必须是**同一个超胞**弛豫出来的 (molecule/supercell 别搞混)。")
    nat = sum(ccnt)
    body = clines[8:8 + nat]
    if len(body) < nat:
        raise SystemExit(f"CONTCAR 坐标行只有 {len(body)} 行, 少于 {nat} 个原子")

    # 3) 晶格变化提示 (ISIF=2 时不应变)
    lc, lr = lattice_of(clines), lattice_of(base_lines)
    worst = 0.0
    for i in range(3):
        for j in range(3):
            d = abs(lc[i][j] - lr[i][j])
            worst = max(worst, d / max(abs(lr[i][j]), 1e-6))
    if worst > 1e-3:
        print(f"[warn] 晶格相对变化最大 {worst:.2e} —— 若用了 ISIF=2 不该发生; "
              f"5 个磁态都用这个新晶格, 但请确认是有意为之 (例如 ISIF=3 全弛豫)。")

    # 4) 位移检查: 静态自洽 (IBRION=-1) 同样会输出 CONTCAR, 但坐标与输入 POSCAR 逐位相同。
    #    那种情况下 5 个磁态其实还是未弛豫结构 —— 不报错但整条反解前提是错的, 必须吼出来。
    maxd = 0.0
    cc, rc = coords_cart(clines, lc, nat), coords_cart(base_lines, lr, nat)
    for i in range(nat):
        for j in range(3):
            maxd = max(maxd, abs(cc[i][j] - rc[i][j]))
    if maxd < 1e-6:
        print("[warn] CONTCAR 与 POSCAR 逐位相同 (最大坐标位移 0.0000 A)")
        print("[warn] 这通常是**静态自洽** (IBRION=-1) 的产物, 结构并没有被弛豫 —— 回填等于没回填。")
        print("[warn] 请确认 runs/fm/ 跑的是 03_dft/INCAR.relax (IBRION=2, NSW=60),")
        print("[warn] 而不是 runs/fm/INCAR (那是静态的, 默认 IBRION=-1)。")
        print("[warn] 若你本就用的是已弛豫好的结构, 可忽略本警告。")
    else:
        print(f"最大坐标位移: {maxd:.4f} A  (弛豫前 -> 弛豫后)")

    print(f"CONTCAR : {cont}  (共 {nat} 原子, {' '.join(cspe)})")
    print(f"目标磁态: {', '.join(orders)}")
    print(f"{'目录':<10}{'POSCAR.pristine':<18}{'POSCAR'}")
    for o in orders:
        src = os.path.join(cases, o, "POSCAR")
        pristine = src + ".pristine"
        has_pristine = os.path.exists(pristine)
        print(f"{o:<10}{'已存在' if has_pristine else '将新建':<18}{'将覆盖' if not args.dry_run else '(dry-run)'}")
        if args.dry_run:
            continue
        if not has_pristine:
            shutil.copyfile(src, pristine)          # 存档原始未弛豫坐标
        own0 = open(src, encoding="utf-8", errors="ignore").read().splitlines()[0]
        new = [own0.split(" (relaxed)")[0] + " (relaxed)"] + clines[1:]
        with open(src, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(new) + "\n")
        rp = os.path.join(runs, o, "POSCAR")
        if os.path.isdir(os.path.dirname(rp)):
            with open(rp, "w", encoding="utf-8", newline="\n") as f:
                f.write("\n".join(new) + "\n")

    if args.dry_run:
        print("\n(dry-run: 没有写任何文件)")
    else:
        print("\n回填完成。5 个磁态现在共用同一弛豫后结构, INCAR 里逐原子 MAGMOM 未动。")
        print("下一步: 每个 runs/<order>/ 跑静态自洽 (INCAR.static), 然后")
        print(f"        python fgt_verify_states.py \"{cases}\"   # 先核对收敛到目标 sigma")
        print(f"        python solve_fgt_J.py \"{cases}\" --spin 1.0")


if __name__ == "__main__":
    main()
