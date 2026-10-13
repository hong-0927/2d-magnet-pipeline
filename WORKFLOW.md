# 二维磁性单层通用流水线 — MP → 磁耦合 J → Monte Carlo

给一个二维磁性材料的化学式，自动完成：**Materials Project 下载结构 → 层对齐 → 超胞/四磁态 → 磁耦合常数 J(J1/J2/J3) 能量映射 → 经典 Heisenberg Monte Carlo 临界温度 (T_C / T_N)**。

> 从 `MnPS3_pipeline` 抽象而来（它是已验证实例）；继承了 `CrBr3_pipeline` 的参数化脚本与 `mp_2d` 的下载脚本。之前的 CrI3/CrBr3 流水线只支持 J1/J2 + 铁磁总磁化，本框架扩展了 **J3** 与 **反铁磁交错序参量**，并把 MC 引擎重写成自包含版本（不再依赖外部 `MnSTe_repro` 仓库）。

---

## 1. 一句话运行

```powershell
cd <repo>

# 已收录文献 J 的材料（如 MnPS3）—— 最简：
python run_2d_magnet.py --formula MnPS3 --literature MnPS3 --mp-key <MP_API_KEY>

# 新材料（手动给参数）：
python run_2d_magnet.py --formula VI3 --element V --nonmag I --S 1.0 --moment 2.0 --U 3.0 `
    --J1 -1.0 --J2 -0.1 --J3 0.0 --MAE 0.2 --mp-key <KEY>

# 已有 POSCAR（跳过下载）：
python run_2d_magnet.py --poscar my.vasp --element Cr --nonmag I --S 1.5 --moment 3.0 --U 3.0 `
    --J1 2.4 --J2 0.6

# 无 MP key，改用公开 provider 链（COD/OQMD/AFLOW/JARVIS，返回体相）：
python run_2d_magnet.py --formula MnPS3 --literature MnPS3 --provider cod
```

输出到 `<公式>_pipeline/`（`--outdir` 可改）：`workflow.md`、`hamiltonian.md`（顶层海森堡哈密顿卡片）、`01_structure/`、`02_supercell/`、`03_dft/`、`05_mc/`（含自包含的 `mc_2d.py` + `poscar_to_mc3.py` + `write_hamiltonian.py`）、`figures/`。

---

## 2. 输入参数

| 参数 | 含义 | 例 (MnPS3) |
|---|---|---|
| `--formula` | MP 化学式 | `MnPS3` |
| `--element` | 磁性元素 | `Mn` |
| `--nonmag` | 非磁元素（可多个，空格分隔） | `P S` |
| `--S` | 自旋大小 | `2.5`（Mn²⁺ 5/2） |
| `--moment` | 磁矩初猜 μB | `5.0` |
| `--U` | DFT+U 加磁性元素 d 轨道 (eV) | `5.0` |
| `--J1/--J2/--J3` | 交换常数 meV，约定 `H = -J S·S`（J>0 铁磁） | `-0.79/-0.04/-0.23` |
| `--MAE` | 单离子各向异性 meV（默认 0） | `0.0` |

**J 符号约定**：`H = -J Σ S·S`，**J>0 = 铁磁，J<0 = 反铁磁**。多数文献用 `H = +J S·S`（J>0=反铁磁），转过来**要变号**。

---

## 3. J 值从哪来（关键）

1. **优先查** `literature_J.json`（当前只有 MnPS3 是**已核实**文献值）。
2. 否则**联网查文献**（Sivadas / Chittari / Wildes 等），或让用户给，**绝不猜**。
3. 或**真跑 DFT**：`03_dft` + `runs/<order>` 已备好，四磁态静态 + FM 弛豫后
   `scripts/solve_J3.py <cases> --demo` 自检 → 用真实能量反解。

新材料写进 `literature_J.json`（字段 `formula/element/nonmag/S/moment/U/J1/J2/J3/MAE/ground_state/source`）。**J 未核实的条目 `J1` 留 `null` 并写 `_note`**。

---

## 4. 内部步骤（driver 自动调用）

```
fetch      mp_fetch_2d.py           MP 下载 + 剥单层（取 hull 最小的候选）
           fetch_structure.py       结构 provider 链: COD/OQMD/AFLOW/JARVIS (--provider, 无 key, CIF/OPTIMADE/AFLUX -> POSCAR)
reorient   reorient.py              层法线对齐 z + 加真空（修 MP 单斜/紧胞缺口）
standardize standardize_cell.py     紧凑胞标准化（把胞内两磁性原子变成最近邻；MP 原胞可能非紧凑，如 CrI3 的 8.05 Å）
supercell  build_magnetic_cases3.py 2×2×1 超胞 + 四磁态 + J1/J2/J3 系数
          build_cluster_cases.py   簇磁体：fm+波矢态 + J_intra/J1/J2/J3
          build_magnetic_cases_tri.py 三角：fm + 120°非共线 + stripy/uud + J1/J2/J3
dft        assemble_dft.py / assemble_dft_tri.py  生成 INCAR/KPOINTS/POTCAR.spec + runs/<order>/ (三角用非共线)
mc         poscar_to_mc3.py + mc_2d.py  建模型 + Heisenberg MC（自写 fast kernel）
           之后 driver 自己再判一次峰位: 峰落在温度窗口首/末格点 (= 窗口截断或网格跨过窄峰)
           → 自动缩放窗口 + 加密退火重扫, 首轮留 `mc_mt_survey.csv` (避免把边界值当 T_C/T_N)
fss        mc_fss.py (--fss 时)  多尺寸 U₄(T) 有限尺寸标度: 两两交点 + collapse 交叉验证 → T_C/T_N
render     render_supercell.py / render_cluster_orders.py / render_tri_orders.py  磁序俯视图
ship-mc    (复制) mc_2d.py + poscar_to_mc3.py + write_hamiltonian.py (簇/三角再加对应枚举+反解脚本) -> <outdir>/05_mc/
hamiltonian write_hamiltonian.py    生成 <outdir>/hamiltonian.md (Heisenberg 哈密顿卡片, 顶层)
```

脚本均在 `scripts/`，可单独调用。`mc_2d.py` 自包含（numpy+numba+matplotlib）。

driver 另支持位置参 `<command>` 分阶段运行：`fetch` / `build` / `dft` / `mc` / `solve` / `render` / `report`。各子命令在阶段边界早退、复用已有产物（如 `mc` 不产 DFT 输入、`report` 不触发 MC/图）；省略 `<command>` 即完整流水线。

---

## 5. 目录结构

```
2d_magnet_pipeline/
├── run_2d_magnet.py     # 单入口 driver
├── WORKFLOW.md          # 本文件
└── scripts/
    ├── __init__.py      # 使 scripts/ 可作包安装 (非可编辑 pip install . 也能定位脚本)
    ├── literature_J.json # 已核实文献 J 表（MnPS3 已验证）
    ├── mp_fetch_2d.py   / reorient.py / standardize_cell.py
    ├── build_magnetic_cases3.py / solve_J3.py
    ├── poscar_to_mc3.py / mc_2d.py (--order cv 比热判据)
    ├── write_hamiltonian.py # 模型 JSON -> hamiltonian.md (海森堡哈密顿卡片)
    ├── cluster_sites.py # 簇磁体：磁性簇→有效位（Nb3I8 等）
    ├── build_cluster_cases.py # 簇磁体：超胞+磁序枚举+J_intra/J1/J2/J3 映射 (含 DFT runs)
    ├── solve_cluster_J.py # 簇磁体：OUTCAR -> E0/J_intra/J1..Jk
    ├── render_cluster_orders.py # 簇磁体：簇心磁序俯视图
    ├── build_magnetic_cases_tri.py # 三角：超胞 + 磁序枚举(含120°非共线) + J1/J2/J3 映射
    ├── assemble_dft_tri.py # 三角：非共线 INCAR (LNONCOLLINEAR) + runs
    ├── solve_J_tri.py # 三角：OUTCAR -> J1/J2/J3
    ├── render_tri_orders.py # 三角：磁序俯视图 (120° 矢量箭头)
    ├── build_fgt_cases.py # **双亚晶格铁磁网络 (Fe3GeTe2/Fe3GaTe2 类)**：超胞(默认 3x3, 含别名自检)
    │                      #   + 磁序枚举 + 系数矩阵 + 逐原子 MAGMOM 的 INCAR/KPOINTS/runs + 俯视图 + 哈密顿卡片
    ├── solve_fgt_J.py     # **双亚晶格反解**：各磁态总能量 -> J1..J4 + E0（`--demo` 自检, 应 PASS）
    ├── assemble_dft.py  / render_supercell.py
```

---

## 6. 从哪些项目继承
| 通用框架脚本 | 来源 | 改动 |
|---|---|---|
| `mp_fetch_2d.py` | `mp_2d/` | 原样复用 |
| `build_magnetic_cases3.py` | CrBr3 `build_magnetic_cases(2).py` | J1/J2→J1/J2/J3、多非磁元素 |
| `solve_J3.py` | CrBr3 `solve_J.py` | J1/J2→J1/J2/J3 |
| `render_supercell.py` | CrBr3 `render_supercell.py` | 加多非磁元素 + J1/J2/J3 壳分色/系数标注 |
| `assemble_dft.py` | CrBr3 `assemble_runs.py` 重写 | 参数化 U/LDAUL/物种 |
| `poscar_to_mc3.py` | MnPS3（源外部 `MnSTe_repro` 仓库） | 加 J3 + 修 bug |
| `mc_2d.py` | 新写 | 自包含 MC + 铁磁/反铁磁序参量 |
| `write_hamiltonian.py` | 新写 | 模型 JSON -> `hamiltonian.md` (海森堡哈密顿卡片) |
| `reorient.py` | 新写 | 修 MP 层法线/真空缺口 |
| `standardize_cell.py` | CrBr3 `standardize_cell.py` | 原样复用：把胞内两磁性原子标准化成最近邻（MP 原胞可能非紧凑，如 CrI3 胞内 8.05 Å = J3 距离，会破坏 M 点枚举） |
| `cluster_sites.py` | 新写 | 簇磁体：磁性簇识别 → 簇质心作为有效磁位（Nb3I8 的 Nb3→1 位） |
| `build_cluster_cases.py` | 新写 | 簇磁体：超胞 + 磁序枚举 + J_intra/J1/J2/J3 映射（含 runs/） |
| `solve_cluster_J.py` | 新写 | 簇磁体：OUTCAR -> E0/J_intra/J1..Jk（对齐 build_cluster_cases --selftest） |
| `render_cluster_orders.py` | 新写 | 簇磁体：簇心磁序俯视图 |
| `build_magnetic_cases_tri.py` | 新写 | 三角：超胞 + 磁序枚举(含 120° 非共线) + J1/J2/J3 映射 |
| `assemble_dft_tri.py` | 新写 | 三角：非共线 INCAR (LNONCOLLINEAR) + runs |
| `solve_J_tri.py` | 新写 | 三角：OUTCAR -> J1/J2/J3 |
| `render_tri_orders.py` | 新写 | 三角：磁序俯视图 (120° 矢量箭头) |
| `mc_2d.py --order cv` | 新写 | 比热 Cv 判据，受挫/非二分晶格（三角）通用定 T_N |
| 峰边界检测 + 自动缩放重扫 | 新写（本仓） | 内置在 `run_2d_magnet.py`：峰落窗口首/末格点 → 缩放加密退火重扫，首轮留 `mc_mt_survey.csv`；`mc_2d.py` 自身也会对边界峰打警告 |
| `build_fgt_cases.py` / `solve_fgt_J.py` | 新写（本仓, 2026-10） | **双亚晶格铁磁网络**（Fe3GeTe2 类）：A(2c)+B(4e) 建模、磁序枚举、系数矩阵(M=[1\|C], cond≈178)、逐原子 MAGMOM、俯视图、哈密顿卡片、总能量反解 J1..J4（`--demo` PASS） |

---

## 7. 已知局限（结论里务必声明）

1. **蜂窝/三角/簇磁体**均支持磁序枚举 + DFT 输入 + 反解脚本 + 图（三角用非共线 120°）；**簇磁体**若是**单重态/dimer（如 Nb3X8）**，经典 Heisenberg 的 T_N **无物理意义**（仅机制演示），需两套 J 另建模，见 skill `references/cluster-magnets.md`。
1b. **多亚晶格铁磁网络（Fe3GeTe2/Fe3GaTe2/Fe5GeTe2 类）不走通用分支**：几何并簇会把它误判成「分子簇」（Nb3I8 与 FGT 签名相同，都是 cluster/1 簇/3 原子），`build_cluster_cases.py` 会明确报错退出。请用 `build_fgt_cases.py` + `solve_fgt_J.py`。**超胞不能太小**：脚本内建 `r_max < 面内边长/2` 判据（2×2 对 3.99 Å 的键已坍缩）。
2. **J 取文献/给定**，未跑 DFT 前非本项目计算值。
3. **MAE 默认 0** → 二维无真正相变（Mermin-Wagner），临界温度为有限尺寸伪相变估计。
4. **`mp_fetch_2d.py` 的 `exfoliate` 有已知缺口**：菱面体层状（层法线沿 [111]）和 MP 已存成单层紧胞会失效；`reorient.py` 修了「层未对齐 + 真空不足」，但多层体相仍需人工剥单层。
5. **k 网格 3×3×1、U 值**为默认/可调假设。
6. **zigzag/stripy 命名**可能和文献互换（M 点两成员标注约定），只影响名字不影响反解。
7. **MC 峰判据必须落在温度窗口内部**才算分辨：峰落在首/末格点 = 窗口截断或网格跨过窄峰，
   该边界值**不是** T_C/T_N。driver 检测到会自动缩放窗口 + 加密退火重扫
   （`--mc-refine-nt` 默认 60 / `--mc-refine-max` 默认 3 / `--no-mc-refine` 关闭），
   首轮原始扫描保留为 `mc_mt_survey.csv/.png`；重扫后仍落边界则结论显式标「未分辨」。
   典型例：MnI2 默认窗口 5–150 K、nt=30（格距 5 K）整个跨过 ~3 K 宽的 Cv 峰，
   曾把 T_N ≈ 7.3 K 报成窗口下限 5.0 K。

---

## 8. Agent 接手检查清单

收到「去 MP 找 XX 结构算磁耦合 + MC」：

1. **先分类磁性格点**：蜂窝（原胞 2 磁原子）/ 三角（1 磁原子）/ 簇磁体（磁性原子成簇）；确认价态（S、μB；簇磁体用簇净自旋）。
2. 查 `literature_J.json`；没有就联网查 J（**转号！**），查不到就问用户或走 DFT。簇磁体需簇间 J。
3. 跑 `run_2d_magnet.py --formula XX --literature ...`（或手动参数）；三角/受挫晶格用 `--order cv`。
4. 读 MC 输出 T_C/T_N（峰判据须落在窗口**内部**；落边界会被自动重扫，仍落边界则标「未分辨」），
   对照实验体相量级（单层偏低）。
5. 在 `workflow.md` 声明：J 来源、MAE、有限尺寸、单层 vs 体相、U/k 假设。
6. **交付 py 文件（ship-mc）**：走 `run_2d_magnet.py` 时由 driver 自动完成（`05_mc/` 内
   `mc_2d.py` + `poscar_to_mc3.py` + `write_hamiltonian.py`，按分支再加枚举+反解脚本，见第 1 节产物清单）。
   **非标准分支（多亚晶格网络，如 Fe3GeTe2）不走 driver ⇒ ship-mc 不会发生，必须手动补齐**
   （首次交付极易漏掉，只剩数据+图而无法脱离主仓重跑），并跑一次
   `python 05_mc/run_rest.py "<05_mc 绝对路径>" 8` 验证自包含：已完成的 run 全 `skip`、退出码 0、
   **目录文件数不变**。参考实例 `Fe3GeTe2_pipeline/05_mc/`（15 个 py）。
