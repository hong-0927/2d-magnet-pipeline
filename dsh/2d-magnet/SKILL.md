---
name: 2d-magnet
description: >-
  Compute magnetic exchange coupling constants J (J1/J2/J3) and classical
  Heisenberg Monte Carlo critical temperature (T_C / T_N) for a 2D magnetic
  monolayer, starting from a Materials Project structure or a POSCAR. Use
  whenever the user wants 磁耦合常数 / 磁交换 / J1 J2 J3 / 磁性 Monte Carlo /
  居里温度 Curie / Néel 温度 / 二维磁体 / 铁磁 / 反铁磁 / 超胞磁序 / 海森堡哈密顿.
  Drives run_2d_magnet.py (--quiet) from the pipeline root (paths below are repo-relative).
---

# 2D 磁体：结构 → 磁耦合 J → Monte Carlo

单入口是 `run_2d_magnet.py`（已封装，勿手写脚本）。流水线根目录 `<repo>`（本 SKILL.md 所在目录的
上两级），脚本在 `<repo>\scripts\`，说明在 `WORKFLOW.md`。

**运行前提（先确认，否则相对路径全断）**：下面所有示例都假定**当前工作目录就是 `<repo>`**
（`cd /d <repo>`）—— `tests/fixtures/...`、`scripts/...`、`run_2d_magnet.py` 一律是相对 `<repo>` 的路径。
从别处调用时，要么先 `cd`，要么用封装脚本 `<repo>\dsh\bin\2d-magnet.cmd`（它按 `%~dp0..\..` 自己定位
`<repo>`，**不依赖 CWD**）。**一条自检入口**（不跑 DFT、也不依赖 CWD）：`python "<repo>\tests\smoke_test.py"`
—— §1–§21 主回归，通过时最后一行是「全部通过」且退出码 0；只想核对某一层就用各脚本的 `--selftest` /
`--demo`（清单见 README 的「脚本级自检」节）。

**命令封装**（Windows cmd，自动带 `--quiet`）：
`<repo>\dsh\bin\2d-magnet.cmd`
等价于 `python "<repo>\run_2d_magnet.py" --quiet <args>`。Python 需 numpy/numba/matplotlib/requests
（`pip install -r requirements.txt` 即可；`.cmd` 会按 `%MAGNET_PYTHON%` → `python` → `py` 顺序找解释器）。

**画图不怕 matplotlib 坏**：matplotlib 在 `import` 期就需要 `pyparsing`/`cycler`/`contourpy`/`fontTools`/`kiwisolver`，缺任何一个都会让 `render`/`mc`/`--fss` 阶段只能打 `[warn] 画图跳过` 而拿不到 PNG。本仓为**每条**渲染路径都配了纯 Pillow 的等价实现（`render_supercell_pil.py` / `render_cluster_orders_pil.py` / `render_tri_orders_pil.py`，同 CLI/输入/输出/配色；曲线类则内置于各自脚本：`mc_2d.py` 的 `_plot_mt_png`、`mc_fss.py` 的 `_plot_fss_png`、`build_fgt_cases.py` 的 `_plot_fgt_png`），共用原语 `scripts/pil_fig.py`；driver 用 `have_matplotlib()` 现场二选一，**装了 Pillow 就能出图**（macOS/Linux 上同样适用）。六条路径全部与 matplotlib 解耦，没有"只能等 matplotlib 修好"的图了。

## 跑之前先定三件事

1. **磁性格点类型**（`detect_lattice` 自动判：先查正交几何，再数磁性原子/原胞）：
   - 蜂窝（原胞 2 磁原子）→ 全链；三角（1 磁原子）→ 含 120° 非共线；簇（成簇）→ 簇模型；
   - **矩形（正交 a≠b，1 磁原子/胞，如 CrSBr/CrSCl）→ `rectangular`**（4 共线态 fm/afm_a/afm_b/afm_diag，见专节）；
   - **方形（正交 a=b，1 磁原子/胞，如 La₂CuO₄/K₂NiF₄）→ `square`**（6 共线枚举、DFT 只提交 4 个 fm/stripe_a/p4_a/mix_ab，见专节）；
   - **C 心矩形（正交 a≠b，2 磁原子/胞，两者由 (½,½) C 心矢量相连，如 CrOCl）→ `rect_centered`**（见专节；**oxyhalide 走这条，不是 `rectangular`**）；
   - **六角 kagome（六角 a=b，3 磁原子/胞，每个磁原子在最近邻 a/2 上恰好 4 个近邻，如 Fe₃Sn）→ `kagome`**（见专节；**只到系数层 P4.1**，J/T_N 待 P4.2）；
   - ⚠ **仍会误判的**：**多亚晶格铁磁网络**（Fe₃GeTe₂/Fe₃GaTe₂/Fe₅GeTe₂ 类）会被**误判成「簇」**（见专节）；
     方形几何 + 2 个磁原子/胞会判成 **`sq_centered`** 并**显式报错**（本仓未支持，但不会静默套方形路径）。
   - 自旋 S / 磁矩 μB：Mn²⁺ S=5/2、Cr³⁺ S=3/2、Fe²⁺ S=2、Gd²⁺ S=7/2（4f7，见局限第 10 条）；簇用簇净自旋。
     **FGT 类不要套 Fe²⁺ 默认值**（不等价 Fe 位，磁矩/J 都须来自文献或自算）。
2. **J 从哪来**：文献（`literature_J.json`；**已核实**：MnPS3 / CrSBr / Fe3GeTe2 / La₂CuO₄；CrI₃ / CrBr₃ / GdI₂ / K₂NiF₄ / **CrOCl** 留 `null`）、用户给定、或跑 DFT 反解。
3. **MP_API_KEY**：`--formula` 下载需要（环境变量 `MP_API_KEY` 或 `--mp-key`）；已有 POSCAR 用 `--poscar` 跳过。

## 遇到没见过 / 相似但没见过的结构：自主分类并归纳完善（自我学习协议）

流水线内置 **8** 类晶格模型（蜂窝 / 三角 / 簇 / 矩形 / **方形** / **C 心矩形** / 双亚晶格 FGT / **kagome**），二维磁体远不止这些。
拿到新结构时**不要直接报错或硬套**，按下述流程自主判定并把它补进本 skill，让下一次能直接处理。

**第一步 · 判定「没见过」还是「相似但没见过」**

- **相似但没见过**：磁性亚晶格类型已在 8 类里（只是换了元素 / 交换机制 / 磁态）→ 走「扩展」，不开新分支。
- **没见过**：磁性亚晶格不属于任何已有类型（新对称 / 磁位不等价 / 多亚晶格网络等）→ 走「新开分支」。

**第二步 · 相似但没见过 → 扩展（成本低，先走这条）**

1. 换元素 → 加 `ELEMENT_PARAMS`（S/moment/U）+ `MODEL_PRESETS`（`--model` 一键）+ `literature_J.json` 条目（**J 留 null 待核实，不自作权威**）。
2. 换机制（如 3d 超交换 → 4f-5d 巡游中介）→ 在 `literature_J.json` `_note` 与 `SKILL.md` 局限里记机制差异，强调 **J 必须 DFT 反解**。
3. 加 fixture + smoke 断言 + README/WORKFLOW 文档 → 回归 → 提交。
   - 例：**GdI₂**（三角晶格 + 4f 磁体，非新晶格）→ 加 Gd 元素表 + `--model gdi2` + 文献条目，走三角分支。

**第三步 · 真没见过 → 新开分支（完整流程，参考 CrSBr）**

1. 判几何：`detect_lattice` 现有判据（矩形 + 数磁原子数）会漏新对称，先设计新的几何判据。
2. 枚举磁态 + 能量映射系数（每态系数表 + `--selftest` 前向合成→反解自检）。
3. 写新 `build_*_cases.py` + `solve_J_*.py`（对齐命名约定），加 fixture + smoke 节。
4. 接线 `detect_lattice` + `main()` 的 build/solve/dft/ship/render/write_workflow 各阶段。
5. 把「为何误判 / 正解怎么走 / J 命名坑」写进 `SKILL.md` 专节 + 局限 → 回归 → 提交。
   - 例：**CrSBr**（矩形晶格，被误判成蜂窝）→ 新开 `rectangular` 分支 + `--model crsbr`。

**硬约束（始终遵守）**

- 几何判据判不出的化学事实（如 FGT 多亚晶格网络）→ **只警告、不擅自改路由**，判断留给用户。
- 新材料的 **J 只作参考**，正式 MC 输入的 J 必须来自用户自己的 DFT 反解；**VASP 只导出输入、不代提交**。
- 每完成一个分支/扩展都要 **回归（smoke + e2e）+ 提交**，不留半成品。

## 命令

```powershell
cd <repo>

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
cd <outdir>\05_mc
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

## 分支：矩形晶格铁磁（CrSBr / CrSCl 类）——已实现

**这是自动分类会误判的第二类**（2026-10 由 CrSBr 暴露）。CrSBr 是**矩形（正交）晶格的 2D 铁磁半导体**：
层内 Cr 排成矩形格子（a ≠ b），交换沿**三个方向各向异性** —— 与蜂窝（六角、2 位/胞）、三角（六角、1 位/胞）都不同。

**为什么会误判**：`detect_lattice` 只数磁性原子数（1/胞→`triangular`、2/胞→`honeycomb`），**不检查几何形状**：

- 体相 CrSBr（2 层/胞 = 2 Cr/胞）→ **误判 `honeycomb`**，走 `build_magnetic_cases3.py`，
  生成蜂窝的 4 态（fm/neel/zigzag/stripy），J 按**距离排序**（J1=3.54 Å、J2=4.76 Å、J3=5.93 Å）；
- 单层 CrSBr（1 Cr/胞）→ 会**误判 `triangular`**（矩形 a≠b，不是三角的 a=b 六角）。

**但 CrSBr 的真实模型是三项各向异性 FM 交换，且交换不按距离排序**（对角线键**最长却是最强**）：

| 键 | 方向 | 键长 | J (meV) | 超交换路径 |
|---|---|---|---|---|
| J1 | ab 对角线 | ≈5.9 Å（**最长**） | **+3.88（最强）** | S–S 共边（近 90° Cr–S–Cr，强 FM） |
| J2 | a 轴 | 3.54 Å | +2.79 | S–Br 共边 |
| J3 | b 轴 | 4.73 Å | +2.12 | S 共角 |

> 所以「按距离给 J1/J2/J3」的蜂窝约定**直接错位**：把最强的对角耦合（J1=3.88）当成最远邻 J3。

**CrSBr 关键参数（已核实 2026-10，出处 arXiv:2108.02418 / PRB 104, 144416 (2021)）**：

- **结构**：Pmmn (No. 59) 正交。单层 DFT a = 3.54 Å、b = 4.73 Å（矩形）；实验体相 a = 3.509 / b = 4.742 / c = 7.953 Å；
  MP 里是 **mp-22998**（体相 a=3.540 / b=4.755 / c=8.394 Å）。单层胞 1 Cr + 1 S + 1 Br（3 原子）；体相 2 层。
- **磁性**：Cr³⁺（d³，S = 3/2），自旋矩 3.06 μB（整胞 3.0 μB/fu），U = 4 eV（GGA+U，J_H = 0.9 eV）。
  **三轴磁各向异性**：易轴 b、中间 a（+16 μeV/Cr）、难轴 c（+34 μeV/Cr）；实验饱和场 b/a/c = 0.58/1.17/2.00 T。
- **层间耦合极弱**：层间 AF 与 FM 仅差 0.04 meV/fu（比层内小 3 个量级）→ 体相是 **A 型反铁磁 T_N = 132 K**（层内仍 FM），
  **单层是 2D 铁磁 T_C ≈ 175 K（理论）/ 146 K（实验）**。→ **体相与单层磁基态不同，报结果必须写明是哪个**。
- **J 符号**：文献用 `H = −J Σ S·S`（J>0=FM），**与本流水线同约定，不需要转号**；E_FM = (−2J1 − J2 − J3)S²。

**怎么用（2026-10 起已实现为独立分支 `rectangular`）**：`detect_lattice` 现在先查矩形几何（a⊥b 且 a≠b），
单层 CrSBr（1 Cr/胞）会正确判成 `rectangular`，走 `build_rectangular_cases.py`（4 共线态）+ `solve_J_rect.py`（反解）：

1. **一键预设**：`--model crsbr` 自动注入 Cr/S/Br、S=1.5 / moment=3.0 / U=4.0 与文献 J（仅参考）。
   `build` 生成 4 个**共线态**：`fm`（全↑）/ `afm_a`（沿 a 反铁磁）/ `afm_b`（沿 b）/ `afm_diag`（对角）。
   `solve` 子命令自动派发 `solve_J_rect.py`，从这 4 态 OUTCAR 总能反解 J1/J2/J3。
2. **J 命名（注意！）**：本流水线 J 一律**按距离壳升序**给：`J1`=最近邻 a 轴 2.79、`J2`=b 轴 2.12、`J3`=对角 3.88 meV；
   论文原始命名是 J1=对角 3.88 / J2=a 2.79 / J3=b 2.12（**对角最长却最强，不按距离排序**），
   所以正确写法是 `--J1 2.79 --J2 2.12 --J3 3.88`，**别照论文直接抄 `--J1 3.88`**。
3. **MC**：`poscar_to_mc3.py` 读真实几何，矩形会自动生成三项各向异性 J 的 MC 模型（`literature_J.json` 已按距离重排）。

**输入侧**：`--model crsbr`（或手动 `--element Cr --nonmag S Br --S 1.5 --moment 3.0 --U 4.0 --J1 2.79 --J2 2.12 --J3 3.88`）；
`--poscar` 给单层 POSCAR（含真空；`tests/fixtures/CrSBr_POSCAR` 已内置）。

## 分支：方形晶格（K₂NiF₄ / K₂CuF₄ / La₂CuO₄ 类）——已实现

**判据**：正交 **a=b** 且原胞 1 个磁原子 → `detect_lattice` 返回 `square`（**必须**有 `_is_square_net` 判据，否则会被误路由到 `triangular` 去解 120°）。
`--model la2cuo4`（Cu, S=1/2, J 已核实于 Coldea PRL 86, 5377 (2001)）/ `--model k2nif4`（Ni, J 未核实）。
`build_square_cases.py` 枚举 **6** 个共线态（fm/g_type/stripe_a/stripe_b/p4_a/mix_ab），**DFT 只提交选定的 4 个** `fm/stripe_a/p4_a/mix_ab`；
J1=a / J2=a√2 / J3=2a，cond(A)≈40.41。`solve_J_square.py` 反解（`--demo` / `--ref la2cuo4` 均应 PASS）。

两个**必须记住的坑**：
1. **J3 与 e0 不可分离**：周期 2 的序一律给出相同 c3 ⇒ 不引入**周期 4** 的 `p4_a` 就解不出 J3（与 kagome 的 q=0 简并同构，成因不同）。
2. **`policy="shift_search"` 是必需的**：J3 = 2a 在 4 宽超胞里是半周期，最小镜像会把 ±2a 折叠到同一像 ⇒ **J3 被低估 2 倍**。

**局限**：这类是 **3D 氧化物、不可解理**，单层切面丢掉层间超交换 ⇒ MC 的 T_C/T_N 只能当量级检查
（La₂CuO₄ 实测 T_N = 325 K 是 3D 驱动；本项目单层 2D 有限尺寸 χ 峰 **266.67 K**，比值 0.82）。

## 分支：C 心矩形 / oxyhalide（CrOCl 类，2 磁原子/胞）——已实现

**判据**：正交 **a≠b** 且原胞 **2** 个磁原子，且两者由 **(½,½) C 心矢量**相连 → `detect_lattice` 返回 `rect_centered`
（旧版会因两 Cr 只差 3.03 Å < 簇判据 3.3 Å 而误返回 `cluster`，走 Nb₃I₈ 分子簇路径）。
`--model crocl`（Cr / nonmag `O Cl` / S=1.5 / moment=3.0 / U=4.0；**J 与 T_N 均未核实**）。
`build_crect_cases.py` 枚举 **7** 个共线态（fm/afm_c/afm_a/afm_b/afm_diag/afm_c_a/afm_c_diag），**DFT 只提交 `fm/afm_c/afm_a/afm_b`**；
J1 = C 心 **3D 键长** 3.0347 Å / J2 = a 3.1980 / J3 = b 3.8603，cond(A)=64.02。`solve_J_crect.py` 反解。

**为什么不能拿 `build_rectangular_cases.py` 顶**：那个脚本的 `--selftest` 理想值表假设 **1 磁原子/胞**，
喂 2 磁原子/胞的材料会报 `系数最大偏差=2.000 FAIL`（能跑但反解口径是错的）。

三个**静默陷阱**（都会"跑得通但物理错"）：
1. **壳容差**：J1=3.0347 与 J2=3.1980 只差 **0.1632 Å** ⇒ `magnetic_shells(tol=0.2)` 会把两壳**并成一个且不报错**。必须 `tol=0.05`。
2. **basis 退化**：C 心键连的是**同一胞内**两个 basis ⇒ basis 符号对一切 (n0,n1)-奇偶序**完全不可见**
   （`afm_c_a == afm_a`、`afm_c_diag == afm_diag` 逐位相同）；母序与其 basis 翻倍版**同选必使 det=0**（脚本硬拒）。
3. **两把尺子**：`poscar_to_mc3.py` 在**面内 2D 投影**里搜键（J1 读到 2.5064 Å），`build_crect_cases.py` 用 **3D 键长**（3.0347 Å，含 1.71 Å 层起伏）。
   两者都对但**不可混用** —— driver 用跨口径不变量 `per_site_bonds`(4/2/2) 在 MC 前校验 J 绑定，不一致直接退出
   （见 `[check] MC 每格点键数 ... == 规格层 ...`）。

**局限**：CrOCl 的基态是**弱铁磁倾斜（需 SOC/DMI）**，本仓库共线 Heisenberg 层**不覆盖** —— 倾斜已落地为
**系数层 + 反解器**（P2′，见下一节），但 θ 是输入参数、一次倾斜 SCF 定不出唯一 (J, D)；
单层切面丢掉层间耦合；**J/T_N 未核实** ⇒ `mc`/`report` 必须显式给 `--J1/--J2/--J3`，`--demo` 必须配 `--ref` 或 `--synthetic`。

## 分支：弱铁磁倾斜（CrOCl 类 oxyhalide 的 P2′）——**系数层 + 反解器；θ 不是自洽解**

**什么时候用**：基态是共线 Heisenberg 不覆盖的**弱铁磁倾斜**（CrOCl 类；判据与 C 心矩形同一套几何 `rect_centered`，
区别只在**磁序**：共线族解释不了总能差，要在母序上叠一个 ±θ 倾斜）。

**工具**：`scripts/solve_J_crect_cant.py`（同一个 cases 目录，第二个反解器；`build_crect_cases.py` 提供倾斜家族）

```powershell
# ① 走 driver（推荐；driver 会自己把倾斜家族交给 builder，solve 端也会自动改派反解器）
python run_2d_magnet.py build --model crocl --poscar <POSCAR> --outdir <out> --cant-theta "8 16"
python run_2d_magnet.py solve --model crocl --poscar <POSCAR> --outdir <out>   # 见到 cant_orders 自动改派
python run_2d_magnet.py dft   --model crocl --poscar <POSCAR> --outdir <out> --cant-theta "8 16" --soc
# ② 脚本级（driver 之外，用于单看审计或合成闭环）
python scripts/build_crect_cases.py <POSCAR> --element Cr --nonmag Cl O --moment 3.0 --spin 1.5 `
    --cant-theta "8 16" --outdir <cases>                                       # 先出倾斜家族（3 分量 magmom.txt）
python scripts/solve_J_crect_cant.py <cases> --orders cant --audit-only        # 只看系数/自由度审计（不需要 OUTCAR）
python scripts/solve_J_crect_cant.py <cases> --demo --synthetic 5 -0.5 0.2     # 无 DFT 闭环自检
python scripts/solve_J_crect_cant.py <cases>                                   # 真值路径（需 <cases>/<order>/OUTCAR）
```
契约测试 `tests/test_crect_cant.py`（44 项）+ smoke §26。

**物理内核是一条闭式关系**（整条分支的地基）：同一母序倾斜 θ 后的系数 =
**共线母序的系数 + `sin²θ ×` 一个与 θ 无关的基向量**（`canted_coeff`，已与独立实现逐位对拍）。
正因为是闭式 ⇒ **θ 是输入参数**：脚本算的是「给定 θ 时系数是多少」，**不是**「基态 θ 是多少」。

**三条边界（写结论时必须带上）**：
1. **θ 不自洽**：真实 θ 需 SOC 总能量的极小 —— 本仓库不做，也没有真实 SOC 能量。
2. **一次倾斜 SCF 定不出唯一 (J, D)**：单壳闭式只有一条冗余自由度 ⇒ `--orders cant` 取的是
   「倾斜家族 **+ 张成所需的 fm/母序行**」，而不是「倾斜行自己」（只给倾斜行会退化成 rank 2）；
   `--audit-only` 打印每行系数与自由度。
3. **SOC 是硬的**：倾斜序属 `mag_vectors.SOC_ORDERS` ⇒ 装配必须 `--soc`，否则**硬拒**（绝不静默按共线算）；
   CrOCl 的文献 J **未核实** ⇒ 真值路径必须 `--ref` 或 `--synthetic`（拒绝猜锚点）。
4. **易轴是假设、路由只是透传**：倾斜按「面内分量 + z 分量」写 MAGMOM（`SAXIS` 缺省 `0 0 1`），
   各向异性大小与易轴方向**都没算过**，`--cant-azimuth` 是取值不是解（φ 由几何定不出来）；
   `run_2d_magnet.py` 现在**已接通**（`build --cant-theta "8 16"` 透传给 builder，`solve` 见到
   `cases/magnetic_orders.json` 里的 `cant_orders` 就自动改派 `solve_J_crect_cant.py --orders cant`，
   非 C 心矩形 / `--cant-base fm` / 有倾斜序却没 `--soc` 三条守卫是硬的）——但**路由不等于物理解**：
   θ 仍是你给的输入，`solve` 只是反解同一组 J 并做自洽校验（倾斜行不增加设计矩阵的秩）。

**不许写**：CrOCl 的 J/D 数值、基态倾斜角、由倾斜推出的 T_N、以及任何「本流水线算了 DMI」的说法。
真值路径**从未跑过真实 DFT**，验收只到「合成能量闭环 + 与独立实现逐位对拍」。

## 分支：六角 kagome（Fe₃Sn 类，3 磁原子/胞）——**只到系数层（P4.1）**

**判据**：六角 **a=b**（γ=120°）且原胞 **3** 个磁原子，且每个磁原子在最近邻距离 **a/2** 上恰好 **4** 个近邻
→ `detect_lattice` 返回 `kagome`。配位数判据是**必需**的：**Fe₃GeTe₂ 同为「六角 + 3 磁原子/胞」，但同一距离上有 6 个近邻**
（金属网而非角共享三角形）⇒ 两者只靠配位数区分；判别在**并簇之前**做，所以 Fe₃Sn 不会再落到 `cluster`（旧版会）。
`--model fe3sn`（Fe / nonmag `Sn` / **S=1.0、moment=2.0、a=5.4 Å 是理想化占位** / U=4.0；**Fe₃Sn 文献 J 未核实**，预设 `literature=None`）。

`build_kagome_cases.py` 枚举 **7** 个非共线磁态（`fm` / q=0 三种 `q0_120`、`q0_ud`、`q0_duu` / √3×√3 三种 `sq3_120`、`sq3_updntri`、`sq3_alt`），
**DFT 只提交选定的 4 个** `fm/q0_120/sq3_updntri/sq3_alt`；J1 = a/2 = 2.70 Å / J2 = a√3/2 = 4.6765 / J3 = a = 5.40，cond(A)=98.51、rank=4。
`solve_J_kagome.py` 反解（`--demo --synthetic 3 -0.4 0.15` 应 PASS），产物走 `j_result.py` 的 `J_result.json`。

> **物理层（真实非共线 DFT → J → T_N）的手册与工具在 `P4_2_KAGOME_RUNBOOK.md`**（仓库根）：入口条件（真实机时 + 已核实的文献锚点 + 参数表填满）、装配与已知坑（`KPOINTS` 默认 `3 3 1` 可用 `--kpoints 6 6 1` 覆盖、`NPAR` 默认 4 可用 `--npar`、static 默认不写 `ISYM` 可用 `--isym 0`、`POTCAR` **不生成**、理想化胞未弛豫）、SCF 四条判据与失败协议、回收通道**已自动化**（`scripts/ingest_dft.py` 的 11 条硬闸门 + `scripts/dft_handoff.py` 的交接包，见手册 §5.0/§5.1）、**k 点收敛抽查**（`python run_2d_magnet.py kconv --outdir <proj> --k-grids 3,4,6,8`，或直调 `scripts/k_convergence.py`：生成 `<序>@<k网格>` 变体、**只换 KPOINTS**、OUTCAR 放回后判「总能量收敛到哪个网格够了」；**总能量收敛 ≠ J 收敛**，它只是必要条件）、反解可信度、机时**纯估算**、以及与原计划的偏差（含 J1/J1′ 分离被明确列为做不到）。**手册与工具都不是进度**：本仓库没有任何真实 kagome J/T_N。

三个**必须记住的坑**：
1. **q=0 一族天然简并**：`fm/q0_120/q0_ud/q0_duu` 满足 **c1 == c2** 且 c3 相同 ⇒ 子块 **rank=2（cond≈2.3e34）**，
   J1 与 J2 在那里**不可分离** —— 这正是 `NEW_TYPES_PLAN.md` 风险行警告的那条。必须上含 **√3×√3** 的序（相位 ℓ=(n0+n1) mod 3，
   在 √3×√3 平移下不变）。反过来注意：**「cond 没超阈」不等于「三个壳都分得开」**，换序后要重算 rank/cond。
2. **超胞 n0/n1 必须是 3 的倍数**（√3×√3 相位周期 3，默认 `--kagome-supercell 3 3 1`）；配位数必须是 (4,4,6)，否则 builder 硬报错。
3. **J1 与 J1′ 被合并**：kagome 有两条**同距离、不同亚晶格对**的最近邻交换（J1 与 J1′），3 壳模型取**平均**并成一条 `J1`
   ⇒ `J_result.json` 里的 `J1` 是两者平均，**不许说成「解出了 J1/J1′」**。要分离需更细的构型（属 P4.2+）。

**局限**：**P4.1 只交付系数层**（纯几何 + 线性代数，`tests/test_kagome_cases.py` 68 项全绿）——**没有任何真实非共线 DFT 能量**，
所以**不许报 J 或 T_N**；真 J 与 T_N 属 **P4.2**。渲染复用三角渲染器（`render_tri_orders*.py`，非共线 3 分量；`render_supercell.py` 只认共线会画错），
`mc`/`report` 必须显式给 `--J1/--J2/--J3`；`literature_J.json` 里**没有** Fe₃Sn 条目。

## 分支：非共线 + SOC（自旋倾斜 / 螺旋；kagome 的 DFT 前置）——能力已实现，物理未验证

**怎么开**（`dft` 阶段）：
```bash
python run_2d_magnet.py dft --model crocl --poscar ... --noncollinear        # 3 分量 MAGMOM (LNONCOLLINEAR)
python run_2d_magnet.py dft --model crocl --poscar ... --soc --saxis 0 0 1   # 再加 LSORBIT + GGA_COMPAT=.FALSE.
```
- 默认 = 共线 `ISPIN=2`；`--noncollinear` = `LNONCOLLINEAR=.TRUE.` + 3 分量 MAGMOM；`--soc` = 再加
  `LSORBIT=.TRUE.` + `GGA_COMPAT=.FALSE.`（可配 `--saxis x y z`）。驱动日志与 `workflow.md` 都会标出「非共线 / 非共线+SOC」。
- **非共线/SOC 下必须没有 `ISPIN`**（与 `LNONCOLLINEAR` 互斥）；`GGA_COMPAT` 不关掉 SOC 会被 GGA 兼容模式抵消。
- 角度表在 `scripts/mag_vectors.py`：**已实现**共线（`spin_collinear`）/ 面内 120°（`spin_120deg`）/ ±θ 倾斜（`spin_cant`）；
  `SOC_ORDERS` 只是 **SOC 门禁的名字表**（含尚无生成器的 `spiral`/`tilt`/`cant_a|cant_b|cant_diag`）——加新磁序要写生成器
  并把名字加进 `SOC_ORDERS`（否则 `--lsorbit` 门禁不生效），`ANGLE_TABLE` 是给人看的角度约定表、不参与分派，
  再让 `build_*_cases.py` 出 3 分量 `magmom.txt`。
- 三角晶格走 `assemble_dft_tri.py`，它只是**薄适配层**（转调 `assemble_dft.py --noncollinear --all-orders`），
  两条路径共用同一份 INCAR 渲染，不会分叉（`--print-cmd` 可看它实际执行的命令）。

**门禁（不静默降级）**：需要 SOC 的磁序（`cant*`/`spiral`/`tilt`）没开 `--soc` ⇒ **硬报错**；
3 分量 `magmom.txt` 配共线模式 ⇒ 硬报错；`--lsorbit` 缺 `--noncollinear` ⇒ 硬报错。理由：缺 SOC 时这些构型
物理上根本不成立，静默降级会产出「看着正常但物理错」的 J。

**能量口径（P3 起统一）**：所有 `solve_*` 只认 OUTCAR 的 **`free energy  TOTEN`**、取**最后一次**出现；
缺这一行就报错，**不会**改用 `energy(sigma->0)` —— 两种口径混用会让反解出的 J 差一个常数项
（`scripts/dft_energy.py` 的 `mixed_convention()` 专抓这种混用，阈值 1e-3 eV）。E.txt 兜底只在没有 OUTCAR 时启用。

**收敛判据**：`EDIFF=1E-6`、`NELM=200`、`ENCUT=520`、`ISMEAR=0 + SIGMA=0.05`（金属自行改 `ISMEAR=-5`）；
非共线 + SOC 必须 `GGA_COMPAT=.FALSE.` 且**删掉 `ISPIN`**。

**局限（写结论时必须声明）**：本仓库**从未真跑过非共线 DFT**。这套能力的验收只到「INCAR/MAGMOM 内容 +
快照回归 + mock-OUTCAR 闭环回解 J」（`tests/smoke_test.py` §21）——所以 `--soc` 的含义是**把 VASP 输入写对**，
**不是**「SOC 物理已验证」：不许据此报 T_C/T_N 或 J。DMI 的反对称项在 MC 侧已就位，但 J 提取侧仍需自旋螺旋法；
kagome 的 q=0 简并已按上面「分支：六角 kagome」处理（P4.1 给系数层 + 硬门禁，不许静默出 J），
但 kagome 的真实非共线 DFT 与 J/T_N 仍属 **P4.2**。

## 分支：DMI / MTI（P6）——**只做判别，不进主链**

遇到「弱铁磁倾斜 / 螺旋 / 手性 / skyrmion / 拓扑绝缘」字样，先读 **`DMI_MTI_MEMO.md`**（判据 + 拒绝模板 + 可行性备忘）。三句话版：

- **DMI**（Dzyaloshinskii–Moriya）：判据 = 重元素 + 强 SOC + **反演对称破缺**（界面/表面/二元层状）；实验线索 = 弱 FM 倾斜、螺旋/摆线、skyrmion、交换偏置。**本仓不提取 D**（要 SOC 四态法或自旋螺旋法的第一性计算 + 3 分量 D 矢量）。
- **MTI**（磁性拓扑绝缘体）：判据 = 层状 vdW + 拓扑能带（Z₂ / 表面态）+ 层间 A 型 AFM。**本仓不做能带拓扑**，也不做层分辨 A 型 AFM 与层间 `J_c`（识别线索可用 `mp_fetch_2d.py` 的 `split_layers` / `layered_metrics`）。
- **能说 / 不能说**：可以说「该材料属 DMI/MTI 候选类，本流水线未覆盖其关键物理」；**不许**说「本流水线算了 DMI」「该材料是 DMI 磁体 / 拓扑绝缘体」，也不许报 D 分量、`|D|/J`、skyrmion 半径、拓扑荷 `Q`、B–T 相图、层间 `J_c`、Z₂、T_N。

**唯一可执行的例外** —— 把**文献 D 值**塞进格点模型跑 MC（≠ 算了 DMI）：

```bat
python scripts\dmi_matrix.py --D 0.31 -0.77 0.55 --J 1.9 --S 1.5 --d 2.70 --tol 0.05 --out shell.json
python scripts\poscar_to_mc3.py <POSCAR> --element Fe --S 1.5 --shell-def shell.json -o model.json
```

`J_D = -[D]×`（对应 `H_DM = +D·(S_i×S_j)`；**反号就变成帮倒忙**），矩阵已把 `S²` 吸收进去，
反向键由引擎用 `J.T` 处理，所以反称部分会被正确消费。两个坑：给两个脚本的 `--S` 必须一致；别用 shell 的
`>` 重定向写 JSON（PowerShell 会写成 UTF-16LE，用 `--out`）。只能看**零场基态能** —— 引擎无外场、无拓扑荷，
**做不出 skyrmion 的 B–T 相图**。结论里必须写明「D 来自文献 ⟨出处⟩」。

## 分支：交变磁体（altermagnet，MnTe / CrSb / RuO₂ 类，P5）——**只做几何判别**

遇到「交变磁体 / altermagnet / 自旋劈裂 / 零净磁矩但 k 分辨劈裂」字样，先跑判别器：

```bat
python scripts\altermagnet_check.py <POSCAR> --element Mn --nonmag Te     REM 只给几何结论
python scripts\altermagnet_check.py <POSCAR> --element Mn --spins 1 -1    REM 可选: 附带补偿性判据
```

判据是**可判定的几何问题** —— 两套磁性亚晶格的连接操作是什么：

- **平移**连接 ⇒ 常规共线反铁磁（自旋劈裂被对称性禁止，`(gT)²` 给 −1）。
- **反演**连接（且两套亚晶格共面）⇒ 宇称禁止劈裂（`(gI)²` 给 −1），**也不是**交变磁体。
- **旋转 / 镜面**连接 ⇒ **几何候选**（`(gR)² = +1`，劈裂被允许）。

**判据到此为止**：是否真有 k 分辨自旋劈裂**要能带计算**，本仓不做能带 ⇒ 结论里只能写「几何候选 / 常规共线反铁磁 /
宇称禁止」，**不许**写「这是交变磁体」「劈裂多大」「带隙里有劈裂」。三个已知坑：
① 判别用的容差以 **Å** 计（不是分数坐标元组比较；坐标抖 3e-6 也能认出对称操作）；
② **镜面连接也有可能是常规反铁磁**（例如 MnPS3 是蜂窝 Néel 序 + 存在镜面），所以「候选」不等于交变磁体；
③ 磁性原子必须恰好两套亚晶格（3 套 → 不适用）。

序参量侧：`mc_2d.py --order sublattice:<mask>`（driver `--mc-order`）算**逐子晶格 Néel 矢量**的标量 |L| 与
Binder U4 —— `sublattice:0` 就是原来的交错序参量（双亚晶格体系逐位等价），`num_basis>=3` 时可以取任意子晶格组合
（老的 `--order staggered` 只会把 basis 1..n−1 并成一组，现在会打 warn）。**它只给 |L|，不分辨两套亚晶格的
Néel 矢量方向** ⇒ 交变磁体的序参量判据仍属未覆盖。掩码必须留至少一个子晶格在 B 组（否则除零，脚本硬报错）。

真 fixture 的判别结论：**MnPS3 几何候选**（镜面连接，残差 0.012 mÅ）、**CrOCl 宇称禁止**（反演连接）、
**Fe3Sn 不适用**（3 套亚晶格）。契约测试 `tests/test_altermagnet.py`（83 项）。

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
另：**矩形晶格铁磁（CrSBr/CrSCl 类）已有独立分支**（`detect_lattice` 先查 a⊥b 且 a≠b → `rectangular`），
走 `build_rectangular_cases.py`（4 共线态）+ `solve_J_rect.py`；注意 J 按**距离壳升序**给（J1=a 轴 2.79 / J2=b 轴 2.12 / J3=对角 3.88），
与论文原始命名（J1=对角）不同，别照抄。
另：**4f 磁体（Gd²⁺ S=7/2，如 GdI₂）不是新晶格类型，走既有三角分支**。GdI₂ 单层是 2H-TaS₂ 型，
三角 Gd 亚晶格（单层胞 1 Gd），`detect_lattice` 正确判 `triangular`。但交换机制与 3d 超交换不同 ——
层内铁磁交换经 Gd 的 5d 巡游电子中介（double-exchange 类似），比纯超交换强，单层 T_C≈241 K（GdI₂, MC）。
Heisenberg 模型仍适用（S=7/2 近乎纯自旋），但 **J 必须来自 DFT 反解**（4f 自动加 `LDAUL=3`，组装脚本按元素
判 d/f 轨道；`--U 6.0` 只是经验假设值）；已加 Gd 元素表兜底（S=3.5/moment=7.0/U=6.0）+ `--model gdi2` 预设 + `literature_J.json`
GdI₂ 条目（J 留 null 待核实）。
另：**方形晶格（K₂NiF₄/K₂CuF₄/La₂CuO₄ 类）是 3D 氧化物、不可解理**：单层切面丢掉层间超交换 ⇒
MC 的 T_C/T_N 只能当量级检查（La₂CuO₄ 实测 T_N = 325 K 是 3D 驱动，本项目单层 χ 峰 266.67 K）。
方形**必须**用 `policy="shift_search"` 且**必须**含周期 4 的 `p4_a`，否则 J3 低估 2 倍 / 与 e0 不可分离。
另：**C 心矩形（CrOCl 类，`rect_centered`）三个静默陷阱**：① 两 Cr 只差 3.03 Å < 簇判据 3.3 Å（旧版会当分子簇）；
② J1=3.0347 与 J2=3.1980 只差 0.1632 Å ⇒ `tol=0.2` 并壳且不报错（必须 `tol=0.05`）；
③ basis 符号对 (n0,n1)-奇偶序不可见，母序与 basis 翻倍版同选必使 det=0。J1 的 3D 键长与 `poscar_to_mc3` 的 2D 投影
（2.5064 Å）都对但不可混用 ⇒ 靠 `per_site_bonds`(4/2/2) 跨口径校验。**CrOCl 的基态是弱铁磁倾斜（需 SOC/DMI），
本仓库共线 Heisenberg 不覆盖**，且其 **J/T_N 均未核实** ⇒ `mc`/`report` 必须显式给 `--J1/--J2/--J3`。
另：**交变磁体（altermagnet）只做几何判别（P5）**：`altermagnet_check.py` 只能给「平移 ⇒ 常规共线反铁磁 /
反演 ⇒ 宇称禁止 / 旋转·镜面 ⇒ 几何候选」，**是否真有 k 分辨自旋劈裂要能带计算** ⇒ 不许写「这是交变磁体」。
MC 侧 `--order sublattice:<mask>` 只给**标量 |L|** 的 U4 交点，不分辨两套亚晶格的 Néel 矢量方向。
判据容差以 **Å** 计，且**镜面连接也可能是常规反铁磁**（MnPS3 即如此）⇒「候选」≠ 交变磁体。
另：**弱铁磁倾斜（CrOCl 类，P2′）只有系数层 + 反解器**（`solve_J_crect_cant.py` + 44 项契约测试）：
闭式 `sin²θ` 标度可解析验证，但 **θ 是输入参数不是基态解**、**一次倾斜 SCF 定不出唯一 (J, D)**、
装配**必须** `--soc`、CrOCl 文献 J 未核实 ⇒ 真值路径必须 `--ref`/`--synthetic`。
⇒ 不许报 CrOCl 的 J/D、基态倾斜角或由它推的 T_N；验收只到合成能量闭环。
另：**kagome 的物理层（P4.2）只有操作手册**（`P4_2_KAGOME_RUNBOOK.md`），**没有任何真实 J/T_N**；
手册 §8 的机时是纯估算（先用一个 FM 参考点替换成实测，否则不许当预算），§9 明确「J1/J1′ 分离」做不到。
