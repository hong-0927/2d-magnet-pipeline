# 二维磁性单层通用流水线 — MP → 磁耦合 J → Monte Carlo

给一个二维磁性材料的化学式，自动完成：**Materials Project 下载结构 → 层对齐 → 超胞/四磁态 → 磁耦合常数 J(J1/J2/J3) 能量映射 → 经典 Heisenberg Monte Carlo 临界温度 (T_C / T_N)**。

> 从 `MnPS3_pipeline` 抽象而来（它是已验证实例）；继承了 `CrBr3_pipeline` 的参数化脚本与 `mp_2d` 的下载脚本。之前的 CrI3/CrBr3 流水线只支持 J1/J2 + 铁磁总磁化，本框架扩展了 **J3** 与 **反铁磁交错序参量**，并把 MC 引擎重写成自包含版本（不再依赖 `Desktop\MnSTe_repro`）。

---

## 1. 一句话运行

```powershell
cd C:\Users\吴\Desktop\项目\2d_magnet_pipeline

# 已收录文献 J 的材料（如 MnPS3）—— 最简：
python run_2d_magnet.py --formula MnPS3 --literature MnPS3 --mp-key <MP_API_KEY>

# 新材料（手动给参数）：
python run_2d_magnet.py --formula VI3 --element V --nonmag I --S 1.0 --moment 2.0 --U 3.0 `
    --J1 -1.0 --J2 -0.1 --J3 0.0 --MAE 0.2 --mp-key <KEY>

# 已有 POSCAR（跳过下载）：
python run_2d_magnet.py --poscar my.vasp --element Cr --nonmag I --S 1.5 --moment 3.0 --U 3.0 `
    --J1 2.4 --J2 0.6
```

输出到 `<公式>_pipeline/`（`--outdir` 可改）：`workflow.md`、`01_structure/`、`02_supercell/`、`03_dft/`、`05_mc/`、`figures/`。

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
reorient   reorient.py              层法线对齐 z + 加真空（修 MP 单斜/紧胞缺口）
standardize standardize_cell.py     紧凑胞标准化（把胞内两磁性原子变成最近邻；MP 原胞可能非紧凑，如 CrI3 的 8.05 Å）
supercell  build_magnetic_cases3.py 2×2×1 超胞 + 四磁态 + J1/J2/J3 系数
dft        assemble_dft.py          生成 INCAR/KPOINTS/POTCAR.spec + runs/<order>/
mc         poscar_to_mc3.py + mc_2d.py  建模型 + Heisenberg MC（自写 fast kernel）
render     render_supercell.py      四磁序俯视图
```

脚本均在 `scripts/`，可单独调用。`mc_2d.py` 自包含（numpy+numba+matplotlib）。

---

## 5. 目录结构

```
2d_magnet_pipeline/
├── run_2d_magnet.py     # 单入口 driver
├── literature_J.json    # 已核实文献 J 表（MnPS3 已验证）
├── WORKFLOW.md          # 本文件
└── scripts/
    ├── mp_fetch_2d.py   / reorient.py / standardize_cell.py
    ├── build_magnetic_cases3.py / solve_J3.py
    ├── poscar_to_mc3.py / mc_2d.py (--order cv 比热判据)
    ├── cluster_sites.py # 簇磁体：磁性簇→有效位（Nb3I8 等）
    ├── assemble_dft.py  / render_supercell.py
```

---

## 6. 从哪些项目继承
| 通用框架脚本 | 来源 | 改动 |
|---|---|---|
| `mp_fetch_2d.py` | `mp_2d/` | 原样复用 |
| `build_magnetic_cases3.py` | CrBr3 `build_magnetic_cases(2).py` | J1/J2→J1/J2/J3、多非磁元素 |
| `solve_J3.py` | CrBr3 `solve_J.py` | J1/J2→J1/J2/J3 |
| `render_supercell.py` | CrBr3 `render_supercell.py` | 加多非磁元素 |
| `assemble_dft.py` | CrBr3 `assemble_runs.py` 重写 | 参数化 U/LDAUL/物种 |
| `poscar_to_mc3.py` | MnPS3（源 `Desktop\MnSTe_repro`） | 加 J3 + 修 bug |
| `mc_2d.py` | 新写 | 自包含 MC + 铁磁/反铁磁序参量 |
| `reorient.py` | 新写 | 修 MP 层法线/真空缺口 |
| `standardize_cell.py` | CrBr3 `standardize_cell.py` | 原样复用：把胞内两磁性原子标准化成最近邻（MP 原胞可能非紧凑，如 CrI3 胞内 8.05 Å = J3 距离，会破坏 M 点枚举） |
| `cluster_sites.py` | 新写 | 簇磁体：磁性簇识别 → 簇质心作为有效磁位（Nb3I8 的 Nb3→1 位） |
| `mc_2d.py --order cv` | 新写 | 比热 Cv 判据，受挫/非二分晶格（三角）通用定 T_N |

---

## 7. 已知局限（结论里务必声明）

1. **蜂窝单离子磁体**全链支持；**三角单离子** MC 可用（`--order cv`），但 DFT 反解 J 的三角磁序枚举未实现；**簇磁体/单重态磁体（如 Nb3X8）不适用**——簇间是 dimer/单重态物理，需两套 J（J_intra/J_inter）另建模，见 skill `references/cluster-magnets.md`。
2. **J 取文献/给定**，未跑 DFT 前非本项目计算值。
3. **MAE 默认 0** → 二维无真正相变（Mermin-Wagner），临界温度为有限尺寸伪相变估计。
4. **`mp_fetch_2d.py` 的 `exfoliate` 有已知缺口**：菱面体层状（层法线沿 [111]）和 MP 已存成单层紧胞会失效；`reorient.py` 修了「层未对齐 + 真空不足」，但多层体相仍需人工剥单层。
5. **k 网格 3×3×1、U 值**为默认/可调假设。
6. **zigzag/stripy 命名**可能和文献互换（M 点两成员标注约定），只影响名字不影响反解。

---

## 8. Agent 接手检查清单

收到「去 MP 找 XX 结构算磁耦合 + MC」：

1. **先分类磁性格点**：蜂窝（原胞 2 磁原子）/ 三角（1 磁原子）/ 簇磁体（磁性原子成簇）；确认价态（S、μB；簇磁体用簇净自旋）。
2. 查 `literature_J.json`；没有就联网查 J（**转号！**），查不到就问用户或走 DFT。簇磁体需簇间 J。
3. 跑 `run_2d_magnet.py --formula XX --literature ...`（或手动参数）；三角/受挫晶格用 `--order cv`。
4. 读 MC 输出 T_C/T_N，对照实验体相量级（单层偏低）。
5. 在 `workflow.md` 声明：J 来源、MAE、有限尺寸、单层 vs 体相、U/k 假设。
