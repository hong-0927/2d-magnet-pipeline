---
name: 2d-magnet
description: >-
  Compute magnetic exchange coupling constants J (J1/J2/J3) and classical
  Heisenberg Monte Carlo critical temperature (T_C / T_N) for a 2D magnetic
  monolayer, starting from a Materials Project structure or a POSCAR. Use
  whenever the user wants 磁耦合常数 / 磁交换 / J1 J2 J3 / 磁性 Monte Carlo /
  居里温度 Curie / Néel 温度 / 二维磁体 / 铁磁 / 反铁磁 / 超胞磁序 / 海森堡哈密顿.
  Drives run_2d_magnet.py (--quiet) in C:\Users\吴\Desktop\项目\2d_magnet_pipeline.
---

# 2D 磁体：结构 → 磁耦合 J → Monte Carlo

单入口是 `run_2d_magnet.py`（已封装，勿手写脚本）。流水线根目录
`C:\Users\吴\Desktop\项目\2d_magnet_pipeline\`，脚本在 `scripts\`，说明在 `WORKFLOW.md`。

**命令封装**（Windows cmd，自动带 `--quiet`）：
`C:\Users\吴\Desktop\项目\2d_magnet_pipeline\dsh\bin\2d-magnet.cmd`
等价于 `python "<root>\run_2d_magnet.py" --quiet <args>`。Python 需有 numpy/numba/matplotlib
（Anaconda `C:\ProgramData\anaconda3\python.exe` 与 dsh 自带 Python 均满足）。

## 跑之前先定三件事

1. **磁性格点类型**：数磁性原子/原胞（`detect_lattice` 自动判）：
   - 蜂窝（原胞 2 磁原子）→ 全链；三角（1 磁原子）→ 含 120° 非共线；簇（成簇）→ 簇模型。
   - ⚠ **判据只有几何**（近邻 < `--cluster-cutoff`，默认 3.3 Å）。**多亚晶格铁磁网络**
     （Fe₃GeTe₂/Fe₃GaTe₂/Fe₅GeTe₂ 类）会被**误判成「簇」**，见下文专节。
   - 自旋 S / 磁矩 μB：Mn²⁺ S=5/2、Cr³⁺ S=3/2、Fe²⁺ S=2；簇用簇净自旋。
     **FGT 类不要套 Fe²⁺ 默认值**（不等价 Fe 位，磁矩/J 都须来自文献或自算）。
2. **J 从哪来**：文献（`literature_J.json`，目前仅 MnPS3 核实）、用户给定、或跑 DFT 反解。
3. **MP_API_KEY**：`--formula` 下载需要（环境变量 `MP_API_KEY` 或 `--mp-key`）；已有 POSCAR 用 `--poscar` 跳过。

## 命令

```powershell
cd C:\Users\吴\Desktop\项目\2d_magnet_pipeline

# 文献 J（已核实）
python run_2d_magnet.py --quiet --formula MnPS3 --literature MnPS3

# 新材料（手动参数）
python run_2d_magnet.py --quiet --formula VI3 --element V --nonmag I --S 1.0 --moment 2.0 --U 3.0 `
    --J1 -1.0 --J2 -0.1 --J3 0.0 --MAE 0.2 --mp-key <KEY>

# 已有 POSCAR（跳过下载）
python run_2d_magnet.py --quiet --poscar my.vasp --element Cr --nonmag I --S 1.5 --moment 3.0 --U 3.0 `
    --J1 2.4 --J2 0.6
```

## 特例：多亚晶格铁磁网络（Fe₃GeTe₂ / Fe₃GaTe₂ / Fe₅GeTe₂ 类）

**这是自动分类会误判的一类**（2026-10 由 Fe₃GeTe₂ 暴露，已加代码警告 + 回归测试）。

**为什么会误判**：`detect_lattice` 只有几何判据 —— 磁性原子近邻 < cutoff（默认 3.3 Å）就并簇。
Fe₃GeTe₂（FGT）单层胞内有 **3 个 Fe**，且有两个**不等价位** Fe(1)/Fe(2)，层内 Fe–Fe ≈ 2.6–3.0 Å：

- `detect_lattice` → **`cluster`**（3 个 Fe 并成 1 簇）；
- 于是走「分子簇磁体」分支：`cluster_sites.py` 把**整胞塌成 1 个有效磁位**（Fe₃ 网络 → 一个巨自旋），
  `build_cluster_cases.py` 给的是「簇内 `J_intra` + 簇间 J1/J2/J3」模型。

**但这与 Nb3I8 的分子簇不是一回事**：FGT 是**巡游铁磁金属**，Fe₃ 层是延展的金属网络，
"簇内 J" 不是可测量物理量；分子簇那套「低温可能是单重态、T_N 无物理意义」的结论**不适用**，
反过来把它的 T_N 当作 FGT 的居里温度也**不成立**。

**已加的保护**（`multi_sublattice_warning()`）：当 `cluster` + 1 簇 + ≥3 磁性原子时，
driver 会往 **stderr 打警告**（`--quiet` 下同样可见），并在 `workflow.md` 的「局限」里记一条
第 8 项；回归测试 `tests/smoke_test.py` 第 6 节覆盖（含「真簇磁体 2 簇不触发」的反例）。
**警告只提醒、不擅自改路由** —— 体系是不是分子簇是化学事实，几何判据判不出来。

**FGT 类怎么用这个流水线（2026-10 实测后的结论）**：

1. **正解：走专用的双亚晶格反解工具链**（本仓新增，已自检通过；**全链路 5 步，只有第 2 步要 VASP**）：
   ```powershell
   # 1) 超胞 + 磁序 + 系数矩阵 + DFT 输入 + 图 + 哈密顿卡片 一次生成
   python scripts/build_fgt_cases.py <单层POSCAR> --supercell 3 3 1 --spin 1.0 --moment 2.0 --U 4.0
   # 2) 跑 FM 弛豫 拿 CONTCAR, 回填成 5 个磁态共用结构:
   #    ⚠ 必须用 03_dft/INCAR.relax (== runs/fm/INCAR.relax, IBRION=2/NSW=60);
   #      runs/<order>/INCAR 全是**静态**的 (IBRION 默认 -1), 直接跑它不会弛豫,
   #      但 VASP 照样写 CONTCAR(与 POSCAR 逐位相同) —— 不报错却让整个反解前提失效。
   python scripts/fgt_relax_backfill.py <outdir>\02_supercell\cases   # 会打印最大坐标位移
   # 3) 跑 5 个磁态的静态自洽 (INCAR.static) —— 唯一需要 VASP 的一步
   # 4) 反解前先核对每个磁态真的收敛到目标 sigma (逐原子比符号):
   python scripts/fgt_verify_states.py <outdir>\02_supercell\cases
   # 5) 反解 J1..J4 (+E0):
   python scripts/solve_fgt_J.py <outdir>\02_supercell\cases --runs <outdir>\runs --spin 1.0
   # 不需要 DFT 的自检 (含端到端 mock, 合成能量 -> 反解 -> 校验还原):
   python scripts/solve_fgt_J.py <outdir>\02_supercell\cases --demo     # 应打印 DEMO PASS
   python tests/test_fgt_e2e.py                                        # 应打印 全部通过
   ```
   `fgt_verify_states.py` 是**防静默错误的关键闸门**：能量映射假设每态真的收敛到指定 σ，
   但 VASP 可能从初始 MAGMOM 跑到别的态（尤其 `afm_b_vert` 的 2.55 Å 强键合二聚体反向态），
   那样代进系数矩阵的能量是错的、而且**不报错**。它逐原子比 `magnetization (x)` 的 `tot`
   与 `magmom.txt`，出 `VERDICT: FAIL` 并退出码 1。
   `fgt_relax_backfill.py` 同理是**防"假弛豫"的闸门**：它打印最大坐标位移，位移 ≈0 时警告
   `CONTCAR 与 POSCAR 逐位相同 ... IBRION` —— 那是静态自洽的产物，结构没动过。
   模型 = A(Fe I@2c, 1/胞) + B(Fe II@4e, 2/胞)，壳层从结构实测：
   **J1 A–B 2.63 Å (6 键/胞)、J2 B–B 2.55 Å (1, 胞内竖直对)、J3 A–A 3.99 Å (3)、J4 B–B 3.99 Å (6)**；
   5 个磁态 (fm / afm_ab / afm_b_vert / afm_a / afm_b_plane) 与系数矩阵 **cond ≈ 178**，可逆性好。
   ⚠ **超胞不能太小**：脚本内建几何判据 `r_max < 面内边长/2`，否则周期像会把不同的键折到同一对
   原子上（**2×2 对 3.99 Å 的 J3/J4 已经坍缩**），反解出的 J 会系统性偏大 → 用 3×3 及以上。
   INCAR 里是**逐原子 MAGMOM**（54 原子），不是 `N*2.0` 那种写法。

**第二条路线（完全不用 VASP）：文献 J + 双亚晶格 MC** —— 已跑通，报告与全部产物在
`Fe3GeTe2_pipeline/05_mc/`（**自包含**：`mc_2d.py` 引擎 + 11 个脚本 + 4 个模型 JSON + 扫描数据 + 4 张图
+ `README_mc_literature_J.md`）。整套目录可整体拷走，不依赖 `2d_magnet_pipeline`：

```powershell
cd C:\Users\吴\Desktop\项目\Fe3GeTe2_pipeline\05_mc
# 1) 建模型: A=Fe I@2c / B=Fe II@4e x2 自动识别, 按「距离+配对」把 Kim 2026 的 J1..J4
#    映射到实测壳层, 共 28 根键/胞; 自检 Σ(-J_mat) = +185.28 meV/胞 => 全 FM 基态 -61.76 meV/格点
python build_fgt_mc_model.py ..\01_structure\POSCAR -o fgt_mc_model_mae.json --mae 1.05
#    --mae 0 看各向异性影响; --scale-J 0.5/2.0 做文献 J 的敏感性括号
# 2) 扫描 (8 路并发; 引擎默认用同目录 mc_2d.py; 已完成的 run 自动跳过)
python run_rest.py . 8
# 3) 分析: 抛物线插值定峰 + 标窗口边界 / 种子平均 + 误差棒
python analyze_sweeps.py .
python analyze_seeds.py  .
```

**结果：T_C ≈ 110 K**（四条独立估计 104.5 / 110.5 / 111.2 / 114.6 K，L=24 与 L=36 无系统漂移），
对照单层实验 **130 K**、文献双亚晶格 MC **152 K** ⇒ 低 15% / 28%；同套 J 的平均场估计 **518 K**
⇒ MC/MF ≈ 0.21（2D 磁体的正常压降）。**Kim 的约定与本流水线相同（`H=−J Σ S·S`，J>0=FM）⇒ 不需要转号**；
Pushkarev 那套必须转号且主导耦合差 1.4–4.6 倍。

**这条路线最容易踩的坑（都实测过）**：

- **粗窗 + 单次运行定不出 T_C**：Cv 峰高只有 ~0.10 而 50 K 处已 0.084（峰是"平台"）；
  50–250 K / 17 点的粗窗下峰位散布 **57.9–140.1 K**，细窗（80–200 K / 25 点）+ 2–8 个随机种子
  （`--seed`）才收敛到 104–115 K。**可靠协议 = 降温扫描 `--cool` + 细窗 + 长测量（therm 4000 / meas 12000）+ 多种子。**
- **升温扫描在低温端欠平衡**：L=36 的升温扫描把 Cv 峰顶到窗口下限 50 K（假峰）；峰落窗口边界时引擎会报「未分辨」。
- **主导误差是文献 J 本身**：整套 J ×0.5 / ×2 ⇒ T_C 52 / 224 K（近似线性），
  ⇒ **T_C ≈ 110 K 是「Kim 参数集下的 MC 结果」，不是 FGT 的 T_C 定论**。
- **MC 内核单线程**：`mc_2d.py` 一个 run 只吃一个核，必须并发跑（本机 24 核，16 个 run 并发 4 min vs 串行 30 min）。
- 引擎侧顺带修了 3 处（`2d_magnet_pipeline/scripts/mc_2d.py`）：`--seed` 没暴露到 CLI、
  `--order auto` 对多壳层模型 `KeyError: 'J1_meV'`、模型 JSON 按 GBK 解码会崩（`Invalid \escape`）。

2. **流水线自带的三个分支都不要用**：`honeycomb`(2 位/胞)/`triangular`(1 位/胞) 与 FGT 的
   Fe 网络不对应；`cluster` 分支会把 FGT 判成"分子簇"，`build_cluster_cases.py` 会直接报错退出
   （超胞内只有 1 个簇）。这些是**设计行为**，不是 bug。
3. **输入侧**：`--element Fe --nonmag Ge Te`；S / μB / J / MAE 必须来自文献或自算 DFT。
   `literature_J.json` 的 `Fe3GeTe2` 条已填 **Pushkarev 单层（统一 S=1）+ 转号**
   （J1 = +113、J3 = −10.9、MAE = +0.35、S = 1、μB = 2；J2 该文未印出 → `null` → 落 0.0）。
   **可直接 `--literature Fe3GeTe2` 取用**（已实测参数解析正确），但那只对**通用分支**有意义；
   走上面的双亚晶格反解时 J 是**待求量**，不要预填。
4. **结构性前提**：`--poscar` 给单层 POSCAR（含真空）。**`--formula Fe3GeTe2` 不可用**
   —— MP 里查不到该材料（见「已核实的 FGT 数据」末两条）。

### 已核实的 FGT 数据（2026-10，含出处；用前读符号警告）

- **结构**：空格群 **P6₃/mmc (No. 194)** —— 注意**不是** R-3m（常见误记）。
  Wyckoff：**Fe I = 2c、Ge I = 2d、Fe II = 4e、Te I = 4f** ⇒ 每层胞 **3 个 Fe**
  （不等价位 Fe I / Fe II）。a = 3.991 Å、c = 16.33 Å（XRD 293 K，Deiseroth 2006）；
  单层 DFT a = 4.010 Å。Fe–Fe：J1 2.46 Å、J2 2.63 Å（Kim 2026 表 1）。
- **Tc 参考**：单层实验 ≈ **130 K**（Fei et al., Nat. Mater. 17, 778 (2018)）；
  体相实验 ≈ 200–230 K；双亚晶格 MC：单层 **152 K** / 体相 205 K（arXiv:2009.05701）。
- **⚠ 交换常数的符号约定互相冲突，禁止混用**（三篇量级差 1.5–5 倍）：
  | 来源 | 约定 | 模型 | J1 / 其他 |
  |---|---|---|---|
  | Pushkarev, JMMM 588, 171456 (2023) arXiv:2303.01888（**单层**） | `H=+J S·S`，**J<0 = FM** | 双亚晶格折算为**统一 S=1** | J1 = −113、J3^I = +10.9 meV、A = 0.35 meV/Fe |
  | Kim, PRR 8, 023333 (2026) arXiv:2504.17998 | `H=−J S·S`，**J>0 = FM**（**与本流水线同约定**） | 双亚晶格 | J1 = +78.84、J2 = +24.61、J3 = −2.46、J3′ = −4.56、J4 = −2.13 meV；E_MCA = 3.15 meV/f.u. |
  | Shen, arXiv:2009.05701 | `H=+J S·S`，**J<0 = FM** | 双亚晶格 (Fe1/Fe5) | J(Fe1–Fe2) = −23.82、J(Fe1–Fe5) = −15.52、J(Fe5–Fe5) = +2.50 meV |
  → **转号规则**：本流水线 `H = −J Σ S·S`（J>0 = FM）。`H=+J S·S` 的文献**必须变号**。
  `literature_J.json` 里的 `Fe3GeTe2` 条用的是 **Pushkarev 单层 + 已转号**
  （J1 = **+113**、J3 = **−10.9**、MAE = **+0.35**、S = 1、μB = 2）；该文的 J2 未印出，
  留 `null` → 落回流水线默认 0.0，**结论里必须声明**。
- **单层 POSCAR 怎么来（已跑通，可复现）**：MP 里**没有**这个材料（见下），所以走「自己剥层」：
  1. 拿 AFLOW 原型 CIF：`https://aflow.org/p/A3BC2_hP12_194_ce_d_f-001/aflow.cif`
     （数值坐标只有 CIF 里有；原型页只给符号 `z3/z4`）；
  2. 关键两个数：**Fe II@4e z₃ = 0.6718**、**Te I@4f z₄ = 0.40982**（a = 3.991 Å, c = 16.330 Å, P6₃/mmc）；
  3. 用 24 个对称操作展开体相胞（→ 12 原子 = 6 Fe + 2 Ge + 4 Te），取 z ∈ [0.05, 0.45) 的**一层**
     （→ 6 原子 = 3 Fe + 1 Ge + 2 Te，层厚 Te–Te = 5.220 Å），c 设为 层厚 + 真空（15 Å）；
  4. 部分占据理想化：CIF 里 Fe1 occ 0.83 / Ge 0.99 → 取 1.0（模型惯例）。
  **自校验（务必做）**：算出的 Fe(2c)–Fe(4e) 距离应为 **≈2.63 Å**（文献 J2 键长）；
  得到 2.634 Å ✓。跑通的实例与其出处见 `Fe3GeTe2_pipeline/01_structure/`
  （`POSCAR.fgt_monolayer.vasp` + `STRUCTURE_SOURCE.md` + `aflow_prototype_Fe3GeTe2.cif` + `build_poscar.py`）。
- **⚠ MP 里查不到 FGT**（2026-10 实测，`MP_API_KEY` 有效时的对照实验）：
  `formula=Fe3GeTe2` → 0 条、`chemsys=Fe-Ge-Te`（整个三元空间）→ **0 条**、姊妹化合物
  `Fe3GaTe2`/`Fe-Ga-Te` 同样 0 条；而 `CrI3` → 1、`chemsys=Cr-I/Fe-Ge/Fe-Te/Ge-Te` → 4/17/7/8。
  ⇒ **`--formula Fe3GeTe2` 这条路不可用**，必须 `--poscar` 自带结构（driver 对空候选的报错已写明此点）。
- **⚠ 实际跑进去会被拒（这是设计行为，不是 bug）**：FGT 的 3 个 Fe/胞在 5×5 超胞里被 3.3 Å 判据
  **连成 1 个 75 原子的簇**（不是三聚体），`build_cluster_cases.py` 现在会**明确报错退出**并指向本节
  （以前是一路走到 `compute_coeffs` 才 `KeyError: (4,4,0)` 崩栈，看不出原因）。
  也就是说：**流水线自带的三个通用分支对 FGT 没有可用路径**，必须走本节的两条专用路线
  （第 1 条 DFT 双亚晶格反解、上面那条文献 J + MC；两者共用同一套壳层定义与键表）。注意
  `Nb3I8 → (cluster,1,3)` 与 `FGT → (cluster,1,3)` **签名完全相同** —— 几何判据分不出
  「真分子簇」与「连通磁性网络」，所以那条警告是对**双侧**说的。

## 约定（务必遵守）

- **J 符号**：`H = -J Σ S·S`，**J>0=铁磁，J<0=反铁磁**。文献多用 `H = +J S·S`，转过来**要变号**。
- **MAE 默认 0** → 二维无真正相变（Mermin-Wagner），T_C/T_N 是有限尺寸伪相变估计。
- 三角螺旋磁体（MnI2/NiI2）的基态序与 T_N 是 J2/J3 竞争的后果，120° 枚举只用于反解 J。
- 簇磁体若低温是单重态/dimer（Nb3X8），经典 MC 的 T_N **无物理意义**，仅机制演示。

## 读结果

- 每条 `<公式>_pipeline/`：`workflow.md`（摘要）、顶层 `hamiltonian.md`（海森堡哈密顿卡片）、
  `01_structure/`、`02_supercell/cases/`（磁序+`magnetic_orders.json`）、`03_dft/`+`runs/`、
  `figures/supercell_orders.png`（磁序俯视图）、`05_mc/`（`mc_mt.csv/.png` + 自包含脚本）。
- `--quiet` 结尾只打印 `lattice / J / T_C(T_N) / 产出清单`；`mc_mt.csv` 里是序参量 vs T。
  **T_C/T_N 后面会带温度窗口和格距**——先看这个再信数值。
- **峰判据要落在窗口内部**：峰落在首/末格点 = 窗口截断或网格跨过窄峰，那**不是** T_C/T_N。
  driver 会自动缩放窗口+加密退火重扫（首轮留 `mc_mt_survey.csv`，可调 `--mc-refine-nt`/
  `--mc-refine-max`，`--no-mc-refine` 关闭）；重扫后仍落边界会直接写「未分辨」。
  默认窗口 5–150 K / nt=30（格距 5 K）对 T_C ≲ 30 K 的材料**太粗**，必要时显式收窄+加密。
- 反解 J 需真跑 DFT 后：蜂窝 `solve_J3.py`、三角 `solve_J_tri.py`、簇 `solve_cluster_J.py`（`--demo` 自检）；
  多亚晶格网络用 `fgt_relax_backfill.py`（回填弛豫结构）→ 静态自洽 → `fgt_verify_states.py`
  （核对每态收敛到目标 σ）→ `solve_fgt_J.py`。

## 交付：py 文件必须随产物一起 ship（**非标准分支要手动做**）

约定出处：`WORKFLOW.md:26`（产物清单里 `05_mc/` **含自包含的 `mc_2d.py` + `poscar_to_mc3.py` +
`write_hamiltonian.py`**）、`WORKFLOW.md:72`（ship-mc 步骤）、`WORKFLOW.md:76`（"`mc_2d.py` 自包含"）。

- **走 `run_2d_magnet.py` 的通用路线这一步是自动的**：driver 在 `run_2d_magnet.py:522-530` 复制脚本，
  原注释「附带 MC 脚本到 05_mc/, 使每条流水线可脱离主仓独立重跑 MC」；先复制通用三项，再按分支追加
  枚举+反解脚本（簇 `build_cluster_cases.py`+`solve_cluster_J.py`；三角 `build_magnetic_cases_tri.py`+`solve_J_tri.py`）。
- ⚠ **非标准分支（多亚晶格网络，如 FGT）不走 driver ⇒ ship-mc 不会发生，必须手动补齐**，
  否则交付只有数据+图，脱离主仓就重跑不了（**本会话首次交付就漏了这一步**）。

交付前自检清单：

1. `05_mc/` 内有 **`mc_2d.py`**（自包含引擎，与 `2d_magnet_pipeline/scripts/mc_2d.py` 同版本，
   可用 md5 对比）；按约定再放 `poscar_to_mc3.py`、`write_hamiltonian.py`（通用建模/卡片生成器，
   即便本路线没用上也放，标清哪份是路线专用）。
2. **分支的枚举 + 反解脚本**：双亚晶格 = `build_fgt_cases.py` + `solve_fgt_J.py`
   （对应三角的 `build_magnetic_cases_tri.py` + `solve_J_tri.py`）。
3. **本路线自己的**建模/驱动/分析/核验脚本（双亚晶格：`build_fgt_mc_model.py`、`run_*.py`、
   `analyze_*.py`、`shell_census.py`、`check_bonds.py`）。
4. 目录名用 **`05_mc/`**、顶层有 `workflow.md`（目标/交付/方法/结果/局限）与 `hamiltonian.md`。
5. **验证自包含**：`python 05_mc/run_rest.py <05_mc 的绝对路径> 8` ⇒ 已完成的 run 全 `skip`、
   `exit=0`、**文件数不变**（证明只依赖目录内的 `mc_2d.py`，不指回主仓）。
6. 参考实例：`Fe3GeTe2_pipeline/05_mc/`（15 个 py + 模型 JSON + 原始 CSV/日志 + 图 + README）。

## 局限（结论里务必声明）

J 来源、MAE、有限尺寸、单层 vs 体相、U/k 网格假设；三角/簇的 J 未跑 DFT 前非本项目计算值。
另：**报 T_C/T_N 时必须带温度窗口/格距，并说明峰是否真落在窗口内部**（落边界 = 未分辨）。
另：**多亚晶格铁磁网络（Fe₃GeTe₂/Fe₃GaTe₂/Fe₅GeTe₂）不能用通用分支**：它们会被几何并簇
误判成「分子簇」，警告会发（双侧，因 Nb3I8 与 FGT 签名相同），随后 `build_cluster_cases.py`
会明确报错退出（超胞内只有 1 个簇）。**要用专用的双亚晶格反解工具链**
（`build_fgt_cases.py` + `fgt_relax_backfill.py` + `fgt_verify_states.py` + `solve_fgt_J.py`，
见上文专节）；不要拿文献的 per-Fe J 去凑单胞塌缩模型。
另：**反解出 J 之前必须过 `fgt_verify_states.py`** —— 磁态没收敛到指定 σ 时能量映射是错的，
且不会自己报错（该脚本自检与 mock 端到端测试都在 `tests/test_fgt_e2e.py` 里覆盖）。
