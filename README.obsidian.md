---
title: 2d-magnet-pipeline — README（Obsidian 阅读版）
aliases: ["2D 磁性格点流水线 README", "2d-magnet README"]
tags: [materials-science, magnetism, monte-carlo, dft, vasp, python]
source: README.md
source_sha256: 035790a452b2016507a17e152696ffe2e5b334929bf34edf64631403c886015f
generated_by: scripts/readme_to_obsidian.py
---

> [!info] 这份文件是什么
> 它是 [`README.md`](README.md) 的**阅读版备份**，由 `python scripts/readme_to_obsidian.py` 生成，正文与 README **逐字一致**（只在前面加了 properties、内链目录、管线图和这段说明）。
> README 改过之后请重跑生成器；`python scripts/readme_to_obsidian.py --check` 是同步闸门（smoke §29 会跑它）。

> [!tip] 在 Obsidian 里怎么用
> 把 `2d_magnet_pipeline` **整个目录当 vault 打开**（Open folder as vault）——这样 `KNOWN_ISSUES.md` / `WORKFLOW.md` / `dsh/2d-magnet/SKILL.md` / `examples/` 这些相对链接都能跳。
> 下面的 mermaid 图与 callout 都是 Obsidian 原生渲染，不需要插件；要交汇报可以用 Obsidian 的 Export to PDF。

## 目录

- [[#该读哪份（文档索引）]]
- [[#快速开始]]
- [[#安装]]
- [[#CLI 参考（节选，完整见 `python run_2d_magnet.py --help`）]]
- [[#J 产物与规格层]]
- [[#目录结构]]
- [[#测试]]
- [[#已知局限]]
- [[#路线图]]
- [[#发布状态]]
- [[#License]]

## 管线一览（主链）

```mermaid
flowchart TD
    F["fetch<br/>mp_fetch_2d.py / fetch_structure.py"] --> B["build<br/>build_*_cases.py → POSCAR + magmom.txt"]
    B --> D["dft<br/>assemble_dft[_tri].py → INCAR / KPOINTS / POTCAR.spec + runs/&lt;序&gt;/"]
    D --> V(["外部机时：VASP（本仓库只导出输入，不代提交作业）"])
    V --> I["ingest<br/>ingest_dft.py：11 条硬闸门 → cases/"]
    I --> S["solve<br/>solve_J_*.py → J_result.json"]
    S --> M["mc<br/>mc_2d.py / mc_fss.py → T_C / T_N + Binder U4"]
    M --> R["render / report<br/>PNG + workflow.md"]
    D --> K["kconv<br/>k_convergence.py：只换 KPOINTS"]
    K --> I
    B -.->|不需要 VASP 的旁路| M
    classDef ext fill:#fde,stroke:#c39;
    class V ext;
```

---
# 2d-magnet-pipeline

给一个二维磁性单层，自动完成：**结构 → 磁交换耦合常数 J → 经典 Heisenberg Monte Carlo 临界温度 T_C / T_N**。

- 支持磁性格点：**蜂窝**（J1/J2/J3）、**三角**（120° 非共线）、**分子簇**、**双亚晶格铁磁网络**（Fe₃GeTe₂ 类，走 DFT 反解）、**矩形**（CrSBr 类，正交 a≠b）、**C 心矩形**（CrOCl 类，正交 a≠b 且原胞 **2** 个磁原子，两者由 (½,½) C 心矢量相连）、**方形**（La₂CuO₄ / K₂NiF₄ 类，正交 a=b）、**六角 kagome**（Fe₃Sn 类，角共享三角形网、原胞 **3** 个磁原子，最近邻 a/2 上 4 配位）。
- 两种 J 来源：**文献/给定值 → 正向 Monte Carlo**；或 **DFT 多磁态能量映射 → 反解 J**（后者需要 VASP）。
- 八种晶格共用一套**规格层 + 统一产物**（`scripts/magnetic_orders.py` / `J_result.json`）：J 是怎么算出来的可追溯、可回归比对。见「J 产物与规格层」。
- 附 DeepSeek Harness 的 agent skill（`dsh/2d-magnet/SKILL.md`）。

## 该读哪份（文档索引）

七份文档各管一段，按这个顺序读不会绕路：

| 顺序 | 文档 | 读它为了 |
|---|---|---|
| 1 | 本 README | 装、跑、产物在哪、**已知局限 18 条**——先读局限，才不会把「工具就绪」当成「结果就绪」 |
| 2 | [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) | **症状 → 修复表**（24 条 error key；`doctor` 命中的就是它）。跑挂了第一站来这里 |
| 3 | [`WORKFLOW.md`](WORKFLOW.md) | 每个脚本的出处、阶段交接、局限 1a–1i（比 README 更细的边界） |
| 4 | [`dsh/2d-magnet/SKILL.md`](dsh/2d-magnet/SKILL.md) | agent 技能本体：跑之前先定三件事、逐晶格专节、**文献符号约定冲突表** |
| 5 | [`NEW_TYPES_PLAN.md`](NEW_TYPES_PLAN.md) | 各磁性格点的判定依据与工期、**明确不做什么** |
| 6 | [`P4_2_KAGOME_RUNBOOK.md`](P4_2_KAGOME_RUNBOOK.md) | 真要上机时用：kagome 的装配 / 提交 / 回收 / 判定（**手册不是结果**） |
| 7 | [`DMI_MTI_MEMO.md`](DMI_MTI_MEMO.md) | DMI / MTI 为什么**不进主链**：能说与不能说的清单 |

> 只想先看产物长什么样：直接看 [`examples/`](examples/)（合成/文献 J 的样例产物，不需要 VASP、也不需要装依赖）。
>
> 想用 [Obsidian](https://obsidian.md) 读（内链目录 / mermaid 管线图 / callout / 导出 PDF）：
> 打开 Obsidian → 左下角 vault 图标 → **Open another vault → Open folder as vault** → 选 `2d_magnet_pipeline` **目录本身**
> → 在左侧文件列表里点 [`README.obsidian.md`](README.obsidian.md)。它是本 README 的**生成物备份**（正文逐字一致，由
> `scripts/readme_to_obsidian.py` 生成，`--check` 保证不漂移）。
>
> ⚠️ **Obsidian 只能打开「文件夹（vault）」，不能打开孤立的单个 `.md` 文件**——把 `README.obsidian.md` 拖到桌面双击是打不开的
> （Windows 上 `.md` 默认也没有关联程序），必须让它待在 vault 目录里。以本目录为 vault 时 Obsidian 会生成 `.obsidian/` 配置目录，
> 它已在 [`.gitignore`](.gitignore) 里忽略、也不会进 sdist（`MANIFEST.in` + `tests/check_sdist.py` 双重拦截）⇒ **发布物不受影响**。

---

## 快速开始

```powershell
git clone <this repo>
cd 2d_magnet_pipeline
pip install -r requirements.txt     # numpy / numba / matplotlib / requests / Pillow（兜底渲染器）

# 1) 命名模型预设：一键设定 MnPS3 的 formula/element/nonmag/lattice/S/moment/U（文献 J 仅作参考）
python run_2d_magnet.py --model mnps3 --poscar tests/fixtures/MnPS3_POSCAR

# 2) 文献 J 的正向计算（需 Materials Project key 下载结构）
python run_2d_magnet.py --formula MnPS3 --literature MnPS3 --mp-key <MP_API_KEY>

# 3) 已有 POSCAR + 显式 J（跳过下载；J 需自备）
python run_2d_magnet.py --poscar tests/fixtures/MnPS3_POSCAR `
  --element Mn --nonmag P S --S 2.5 --moment 5.0 --U 4.0 --J1 -0.3 --J2 0.1 --J3 0.02

# 4) 双亚晶格反解（FGT 类，需 VASP；详见 dsh/2d-magnet/SKILL.md 或 hamiltonian.md 用法）
python scripts/build_fgt_cases.py <monolayer_POSCAR> --outdir Fe3GeTe2_pipeline
```

> `--formula` 下载结构需要 Materials Project API key（环境变量 `MP_API_KEY` 或 `--mp-key`）。
> 已有 POSCAR 用 `--poscar` 跳过下载。`--J1` 必填（J>0 铁磁 / J<0 反铁磁），除非用 `--literature` / `--model`（此时 J 取文献参考值，仅作数量级 sanity-check，非本项目计算）。

## 安装

- **Python 3.9+**（`numba` 兼容性），依赖见 [requirements.txt](requirements.txt)。
- 可选：`pip install .`（或 `pip install -e .`）得到 `2d-magnet` 命令。`scripts/` 已包化（含 `literature_J.json`），非可编辑安装后入口也能通过 `import scripts` 定位脚本；否则直接用 `python run_2d_magnet.py`。`pip install .` **不含 Pillow**（它是可选兜底渲染器，matplotlib 装不全时超胞图/簇图/三角图/MC 曲线图自动改走纯 Pillow，见已知局限 12）——要一起装就打 `pip install ".[render-fallback]"`。
- Windows 封装：`dsh/bin/2d-magnet.cmd`（自动带 `--quiet`；解释器查找顺序 `%MAGNET_PYTHON%` → `python` → `py`）。

### 分发、升级与卸载

- **分发形态以本仓库为准**：`git clone` / 复制目录即可用，`python run_2d_magnet.py` 是唯一必需入口；`pip install .` 是可选的便利层（`pyproject.toml` 给 `2d-magnet` 控制台入口，`MANIFEST.in` 给 sdist 的内容清单，闸门是 `tests/check_sdist.py` + CI 的 `package` 作业，见「发布状态」）。因此**「仓库内容完整」等于「分发包完整」**——发布前必须确认新增脚本/测试/fixture/样例都已 `git add`。
- **Windows 封装只为源码树设计**：`dsh/bin/2d-magnet.cmd` 用 `..\..` 反推 `<repo>`，只在裸检出/解压目录里有效；`pip install .` 之后请改用控制台命令 `2d-magnet`（等价于 `python run_2d_magnet.py --quiet`）。两种方式都不要混用。
- **升级**：`git pull`（或覆盖解压新版目录）后强制重跑一次验证——`python run_2d_magnet.py doctor` + `python tests/smoke_test.py` + `python scripts/magnetic_orders.py --regress`。既有产物目录（`<formula>_pipeline/`）建议整目录重建：`magnetic_orders.json` / `J_result.json` 带 schema 版本与 `confidence` 闸门，旧产物不保证被新代码继续接受。
- **卸载**：删目录即卸载。若用过 `pip install .`，用 `pip uninstall 2d-magnet-pipeline` 撤掉 `2d-magnet` 命令与 site-packages 里的 `scripts` 包（**注意**：`scripts` 是个通用包名，若有别的项目也装了同名顶层包，用 `pip show -f 2d-magnet-pipeline` 先确认归属再卸）。
- **不写入的路径**：本流水线不碰注册表、不装服务、不在仓库外留下配置；全部产物都在 CWD 下的 `<formula>_pipeline/`（`--outdir` 可改）。唯一的外部依赖是下载结构时的网络与可选的 `MP_API_KEY` 环境变量。

## CLI 参考（节选，完整见 `python run_2d_magnet.py --help`）

位置参 `<command>`（可选）分阶段运行：`fetch` / `build` / `dft` / `mc` / `solve` / `render` / `report` / `doctor`；省略即完整流水线。见下方「子命令」。

| 组 | 参数 |
|---|---|
| 结构来源 | `--formula` `--mp-id` `--poscar` `--mp-key` `--provider {auto,mp,cod,oqmd,aflow,jarvis}` |
| 命名预设 | `--model {mnps3,cri3,crbr3,fe3gete2,crsbr,gdi2,la2cuo4,k2nif4,crocl}` |
| 材料 | `--element` `--nonmag P S` `--S` `--moment` `--U` |
| J 值 | `--literature <key>` `--J1` `--J2` `--J3` `--MAE` |
| DFT 自旋模式 | `--noncollinear` `--soc` `--saxis 0 0 1`（默认共线 `ISPIN=2`；`--soc` = `LNONCOLLINEAR+LSORBIT+GGA_COMPAT=.FALSE.`，需 SOC 的磁序不开它会硬报错） |
| 弱铁磁倾斜（P2′） | `--cant-theta "8 16"` `--cant-base afm_c` `--cant-azimuth 0.0`（只在 C 心矩形；`--cant-base fm` 与「非 C 心矩形 + `--cant-theta`」硬拒；给了它 `dft` 会**自动**加 `--noncollinear`，但 SOC 仍须显式 `--soc`，否则装配阶段硬拒） |
| 磁性格点 | `--lattice {auto,honeycomb,triangular,cluster,rectangular,rect_centered,square}` `--cluster-cutoff` `--rect-supercell` `--crect-supercell` `--square-supercell` |
| Monte Carlo | `--mc-size` `--mc-tmin` `--mc-tmax` `--mc-nt` `--mc-refine-nt` `--mc-refine-max` `--no-mc-refine` `--mc-seeds` `--fss` `--fss-sizes` |
| 开关 | `--skip-dft` `--skip-mc` `--quiet` `--log` `--outdir` |

### 命名模型预设（`--model`）

`--model <key>` 一键设定 `formula / element / nonmag / lattice / S / moment / U`（材料专属值，取自 `literature_J.json`，与元素表通用兜底不同）：

| key | 体系 | lattice | S / moment / U |
|---|---|---|---|
| `mnps3` | MnPS₃ 蜂窝 Néel 反铁磁 | honeycomb | 2.5 / 5.0 μB / 5.0 eV |
| `cri3` | CrI₃ 蜂窝铁磁 | honeycomb | 1.5 / 3.0 / 3.0 |
| `crbr3` | CrBr₃ 蜂窝铁磁 | honeycomb | 1.5 / 3.0 / 3.0 |
| `fe3gete2` | Fe₃GeTe₂（⚠ 巡游铁磁金属，会被几何判据误判成 cluster） | cluster | 1.0 / 2.0 / 4.0 |
| `crsbr` | CrSBr 矩形晶格铁磁（正交 a≠b，三项各向异性 FM 交换） | rectangular | 1.5 / 3.0 / 4.0 |
| `gdi2` | GdI₂ 三角铁磁（Gd²⁺ 4f7 S=7/2，4f-5d 巡游中介交换，单层 T_C≈241 K 是**理论 MC 估算、非实验值**；文献 J 未核实 ⇒ 留 null） | triangular | 3.5 / 7.0 / 6.0 |
| `la2cuo4` | La₂CuO₄ 方形 G 型反铁磁（Cu²⁺ 3d⁹ S=1/2，**方形定量锚点**，J 已核实于 Coldea PRL 86, 5377 (2001)） | square | 0.5 / 1.0 / 7.0 |
| `k2nif4` | K₂NiF₄（方形晶格类型的**结构原型**，Ni²⁺ 3d⁸ S=1；文献 J 未核实 ⇒ 留 null，只作结构参考） | square | 1.0 / 2.0 / 4.0 |
| `crocl` | CrOCl **C 心矩形**（正交 a≠b 且原胞 2 个 Cr，由 (½,½) C 心矢量相连；J 与 T_N **均未核实** ⇒ 只作结构/序枚举锚点，`mc`/`report` 须显式给 J） | rect_centered | 1.5 / 3.0 / 4.0 |

优先级：显式 CLI > `--model` 预设 > `--literature` 表 > 元素表。`--model` 自带的文献 J 只作**参考**（数量级 sanity-check），落到文献值时会在 stderr 打印 `[ref]` 声明，workflow.md 也会标注「文献参考、非本项目计算」；**真实 J 须跑 VASP + `solve` 反解，或用 `--J1/--J2/--J3` 显式给定**。VASP 侧始终只**导出**输入（`03_dft/` + `runs/<order>/`，含 INCAR/KPOINTS/POTCAR.spec/POSCAR），不替用户提交作业。

### 元素参数表

`--S` / `--moment` / `--U` 未显式给出时，按磁性元素从内置表兜底（**Mn / Cr / Ni / Co / V / Gd**）。Fe 因 Fe²⁺/Fe³⁺/巡游有效自旋三种取值差异太大（S=2/2.5/1）而**故意不列**，须显式给 `--S/--moment`。优先级：显式 CLI > `--model` 预设 > `--literature` 表 > 元素表；落表时会打 stderr 提醒。

### 结构 provider 链

`--provider`（默认 `auto`）选择结构来源：`mp` 走 Materials Project（需 `--mp-key`，会剥单层）；`cod` / `oqmd` / `aflow` / `jarvis` 走公开 REST 接口（**无需 key**，返回体相结构，之后由 `reorient` + `standardize` 处理）。`auto` 先试 MP（有 key 时），失败自动回退 `cod → oqmd → aflow → jarvis`。实现见 `scripts/fetch_structure.py`（CIF/OPTIMADE/AFLUX → POSCAR 三个纯函数转换器 + 可 mock 的 HTTP 层）。说明：这些公开端点可能随时间调整，解析为 best-effort；无网络或全部失败时须 `--poscar` 自带单层 POSCAR。

### Monte Carlo 输出与误差

- 序参量路径（蜂窝铁磁/Néel 反铁磁）的 `05_mc/mc_mt.csv` 多一列 **Binder 累积量 U₄**（`U₄ = 1 − ⟨M⁴⟩/((1+2/3)⟨M²⟩²)`）：有序 → 0.4、高温无序 → 0。跑多个 `--mc-size`，U₄(T) 的交点即有限尺寸标度定出的 T_C（比单一 χ 峰更稳健）。
- `--mc-seeds N`（默认 1）跑 N 个独立种子并平均，峰位报 `± 种子间标准差`（误差棒）。

### 有限尺寸标度（`--fss`）

`--fss` 会在 MC 阶段之后自动跑 `scripts/mc_fss.py`：对 `--fss-sizes`（默认 `12 18 24 30`）逐尺寸算 U₄(T)，取

1. **两两交点**：相邻尺寸的 U₄(T) 在过渡带（U₄∈[0.02,0.38]）内取斜率最陡的符号变化点插值；
2. **交叉验证**：使各尺寸 U₄ 之间方差最小的温度（collapse）。

两者均值即临界温度 `T_C/T_N`，报告 `± 交点间标准差`；结果写 `05_mc/fss_fss.{csv,json,png}`（图与其它渲染路径一样二选一：matplotlib，或 `mc_fss.py` 内置的纯 Pillow `_plot_fss_png`）。独立用法：`python scripts/mc_fss.py <model.json> --sizes 12 18 24`。提示：二维无真正相变（Mermin-Wagner），交点 T_C 本质是有限尺寸伪相变估计，且随种子噪声波动，须 `--mc-seeds > 1` 收敛。

### 子命令（分阶段运行）

`python run_2d_magnet.py <command> [参数]` 只跑某一段，复用已有产物，在阶段边界早退；省略 `command` 走完整流水线。全局开关：`--version`（打印版本号后退出，版本取自 `pyproject.toml`）、`--quiet`（抑制进度输出，但**结论级警告仍打到 stderr**：`[warn] …`、`T_N 未分辨`、`cond(A) 病态`、`过渡带内既无 U4 交点`、`磁态坍缩` 这几类 —— 复审 P0-3 之前它们会被静默吞掉）、`--log <文件>`（把每个阶段子进程的完整 stdout/stderr 追加进去，`--quiet` 下同样有效）。

| 子命令 | 跑什么 | 产出 |
|---|---|---|
| `fetch` | 结构下载（或 `--poscar`）→ `reorient` 加真空 | `01_structure/POSCAR` |
| `build` | 磁性格点识别 → 超胞/磁态枚举（C 心矩形可加 `--cant-theta "8 16"` `--cant-base afm_c` `--cant-azimuth 0.0` 追加 ±θ 弱铁磁倾斜家族；`--cant-base fm` 与「非 C 心矩形 + `--cant-theta`」都被硬拒） | `02_supercell/cases/` |
| `dft` | `build` + DFT 输入组装（自旋模式：默认共线 / `--noncollinear` 3 分量 / `--soc` 加 `LSORBIT+GGA_COMPAT=.FALSE.`；`--kpoints 6 6 1` / `--npar 8` / `--isym 0` 可显式覆盖模板默认，不给就逐字节保持 `3 3 1`+`NPAR=4`；三角走薄适配层 `assemble_dft_tri.py`，转调同一个 `assemble_dft.py`） | `03_dft/` + `runs/<order>/` |
| `ingest` | OUTCAR 回收（`--runs <拷回的目录>` → `cases/`；11 条硬闸门**全过才写**，`--ingest-dry-run` 只校验、`--ingest-lenient` 单条降级） | `cases/<order>/OUTCAR` + `cases/ingest_manifest.json` |
| `kconv` | k 点收敛抽查：`runs/<序>/{INCAR,POSCAR}` → `kconv/runs/<序>@<k网格>/`（**只换 KPOINTS**），OUTCAR 放回同目录后判定「总能量对 k 网格收敛到哪个网格够了」 | `kconv/kconv_plan.json`、`kconv/k_convergence.json` |
| `mc` | `build` + Monte Carlo（J 取文献/CLI，不产 DFT 输入） | `05_mc/mc_mt.csv/.png`、`hamiltonian.md` |
| `solve` | 反解已算 OUTCAR 的 J（需先 `build` 再手动跑 VASP；`cases/magnetic_orders.json` 里有 `cant_orders` 时**自动改派** `solve_J_crect_cant.py --orders cant`，不静默按共线解） | 打印 `J1/J2/J3` (meV) |
| `render` | `build` + 超胞/磁序图（画图二选一：matplotlib，或用纯 Pillow 的等价渲染器 `render_supercell_pil.py` / `render_cluster_orders_pil.py` / `render_tri_orders_pil.py`，driver 按 `have_matplotlib()` 自动探测；`--fss` 的标度图与 FGT 图同理各有一条内置 Pillow 路径） | `figures/supercell_orders.png` |
| `report` | `build` + 汇总卡片（不跑 DFT/MC/render） | `workflow.md` |
| `doctor` | 环境自检（python/依赖/脚本/文献表/**`--model` 预设自检**/MP key/写权限）→ 命中 error catalog 给「症状→修复」 | 控制台诊断 + `KNOWN_ISSUES.md` 参考 |

例：`python run_2d_magnet.py fetch --formula MnPS3 --element Mn --nonmag P S --provider cod`

## J 产物与规格层

八种晶格的 J 提取原先各写各的，壳的定义、系数语义只活在脚本源码里。现在收敛成两层，**新晶格只需写「规格数据」而不必再动公共层**。

### `scripts/magnetic_orders.py` — 公共规格层

把「晶格 + 壳距离 + 磁序」声明成数据，统一算系数矩阵、条件数与残差。九条既有链路重算的系数与它们落盘的 `magnetic_orders.json` **逐位一致**（容差 1e-9，不跑 DFT）：

```powershell
python scripts/magnetic_orders.py --regress                          # 9/9 逐位回归，秒级
python scripts/magnetic_orders.py --emit-spec rect_crsbr --dry-run   # 预览某规格的 order_spec/model_spec
```

一条硬规则：**壳距离必须显式声明，不在运行时自动聚类**。三角走最小镜像、矩形/蜂窝走平移搜索，两种启发式给出的壳本就不同（小超胞下还有别名坍缩），自动聚类无法保证与旧路径同壳。每条键归到「容差内最近的那个显式壳」，没落进任何壳的键会被**点名告警**而不是静默丢弃。

设计矩阵退化（条件数超阈）时**直接报错、不出 J**——这是往后加非共线/SOC 时的安全阀。

### `scripts/j_result.py` — 统一产物 `J_result.json`

```json
{
  "schema": "j_result/1", "material": "CrSBr", "lattice": "rectangular",
  "source": "demo", "confidence": "low", "spin": 1.5, "S2": 2.25, "n_mag": 16,
  "shells": [
    {"name": "J1", "d": 3.5404561399999994, "pair": null, "tol": 0.3,
     "sum_coeff": [-16.0, 16.0, -16.0, 16.0], "z_n": 2.0}
  ],
  "J_meV": {"J1": 2.79, "J2": 2.12, "J3": 3.88},
  "J_tilde_eV": {"J1": 6.2775, "J2": 4.77, "J3": 8.73},
  "cond": 32.0, "residual_eV": 5.551115123125783e-17
}
```

示例里的 `J_meV` 取自 CrSBr 的文献条目（Ke Yang & Hua Wu, PRB 104, 144416 (2021)；原论文命名 `J1=对角 / J2=a / J3=b`，**已按本仓「距离壳升序」重排**）——它只是 `source=demo` 的 schema 样例，**不代表本仓算过 CrSBr**。`confidence` **不可手改**：`source` 不是 `vasp` 只给 `low`；残差 >1e-4 或条件数 >1e3 降 `low`；残差超 1e-2 直接**拒绝写出**（不静默降级）。`validate` 会重算一遍，手写 `"confidence":"high"` 骗不过去。每个壳自带 `tol`，消费端不必自己猜容差。`solve_J_rect.py --jfile-out` 已经会写这个文件。

### 自描述块

`build_rectangular_cases.py` / `build_square_cases.py` 现在把 `order_spec[]`（每个磁序的 `kind`/`coeffs`/`rule`）和 `model_spec`（`lattice`/`policy`/`nshift`/`shells`）直接写进 `02_supercell/cases/magnetic_orders.json`：产物自己说清「我是怎么算的」，规格层一旦与产物漂移立刻可见。其余链路可用 `--emit-spec <KEY> --into <json>` 回填。

### `scripts/poscar_to_mc3.py` 的对齐参数

| 参数 | 作用 |
|---|---|
| `--jfile <J_result.json>` | J 与壳距离都从统一产物读 |
| `--shell-def <json>` | 显式壳 `{shells:[{name,d,tol,pair,J_meV}]}`，替代「按距离排序取最近 3 壳」 |
| `--element X Y` | 多磁性元素（如 FGT 的 Fe 双亚晶格） |
| `--family A=1 B=0,2` | 给 `pair:"AB"` 这类族标签指定 basis 索引 |

旧命令 `--element Cr --S 1.5 --J1 .. --J2 .. --J3 ..` **行为逐位不变**（已与既有 `mc_model.json` 比对）。三个防呆：多个不同 `--S` 直接拒绝（引擎只有一个全局 S²，逐壳 S_eff 留给后续）；`pair:"AB"` 忘配 `--family` 报错并给出两种写法；漏掉中间壳时那条键被点名而不是静默蒸发。

## 目录结构

```
2d_magnet_pipeline/
  run_2d_magnet.py        # 主入口（单文件驱动整个流水线）
  scripts/                # 各阶段 + 分模型解算器（solve_J3 / solve_cluster_J / solve_fgt_J / solve_J_rect / solve_J_square / solve_J_crect / solve_J_crect_cant / solve_J_kagome / ...）
                          #   magnetic_orders.py = 公共规格层（壳/磁序/系数矩阵/条件数）
                          #   j_result.py       = 统一产物 J_result.json 的 schema 与信任度校验
                          #   dmi_matrix.py     = D 矢量 -> 逐壳 3x3 反称交换矩阵（P6；喂 poscar_to_mc3 --shell-def）
                          #   altermagnet_check.py = 交变磁体几何判据（P5：两套亚晶格的连接操作是平移/反演还是旋转）
                          #   solve_J_crect_cant.py = CrOCl 弱铁磁倾斜反解（P2'：闭式 sin²θ 标度 + 倾斜家族 + 自由度审计）
                          #   k_convergence.py  = P4.2：k 点收敛抽查（生成 <序>@<网格> 变体 + 回收后判定"多密才够"）
                          #   readme_to_obsidian.py = 生成本 README 的 Obsidian 阅读版备份（确定性转换 + --check 同步闸门）
                          #   __init__.py 使其可作包安装; literature_J.json 已迁入此处作为包数据
  dsh/2d-magnet/SKILL.md  # agent skill
  dsh/bin/2d-magnet.cmd   # Windows 封装
  tests/                  # smoke_test.py + test_magnetic_orders.py + test_fgt_e2e.py + test_square_cases.py + test_crect_cases.py + test_kagome_cases.py + test_altermagnet.py + test_crect_cant.py + test_ingest.py + check_sdist.py（sdist 内容闸门，不是测试；fixtures/ 内置，无外部路径依赖）
                          #   fixtures/ 里另有 GdI2 / K2NiF4 / La2CuO4 / CrOCl / CrSBr / Nb3I8 / Fe3GeTe2 / Fe3Sn 的 POSCAR 与
                          #   gdi2_nc_soc_120deg.INCAR.norm（非共线 120° 的 INCAR 快照基准）
  NEW_TYPES_PLAN.md       # 新磁性格点类型（方形/kagome/交变磁体…）的扩展计划与工期
  DMI_MTI_MEMO.md         # P6 备忘：DMI/MTI 为何不进主链的判据、能说/不能说什么、可行性估算
  P4_2_KAGOME_RUNBOOK.md  # P4.2 手册：kagome 真实非共线 DFT 的装配/提交/判定/回收/反解/验收（不是结果）
  WORKFLOW.md             # 阶段与脚本来源说明
  KNOWN_ISSUES.md         # 已知问题 / error catalog（doctor 子命令联动的「症状→修复」表）
  README.obsidian.md      # 本 README 的 Obsidian 阅读版备份（**生成物**：正文逐字一致 + properties/内链目录/mermaid 管线图）
  LICENSE                 # MIT
  examples/               # 样例产物快照（见 examples/README.md）：方形 La2CuO4 文献-J 全程 + kagome 合成 dry-run（不需要 VASP）
  templates/              # 扩展骨架（不参与运行）：material_preset.json + new_lattice_checklist.md
  requirements.txt        # 运行期依赖（numpy / numba / matplotlib / requests，另有可选兜底 Pillow）
  pyproject.toml          # 打包基线（`2d-magnet` 控制台入口 + scripts 包数据）
  MANIFEST.in             # sdist 内容清单（文档/手册/skill/样例/夹具/骨架；wheel 只装 run_2d_magnet.py + scripts/）
  .github/workflows/ci.yml # CI：windows-latest + Python 3.12，跑 doctor + regress + smoke + 8 个契约测试，并打包校验 sdist
```

## 测试

```powershell
python tests/smoke_test.py            # 主回归：分类、峰边界、工作流文档、J_result 契约、poscar_to_mc3 新旧路径等价
python tests/test_magnetic_orders.py  # 规格层契约 + 八种晶格逐位回归（8 项，实跑输出 8/8）
python tests/test_square_cases.py     # 方形：构建器 / 规格层 / 手推理想值 三方交叉锚定（28 项）
python tests/test_crect_cases.py      # C 心矩形：basis 退化不可见性 / det 奇异性 / tol=0.2 并壳陷阱 / per_site_bonds
python tests/test_kagome_cases.py     # kagome：q=0 简并 (rank=2) / √3×√3 可解 (rank=4, cond=98.51) / 配位数判据 / mock 闭环（68 项）
python tests/test_fgt_e2e.py          # FGT 双亚晶格反解链路端到端（mock DFT，不跑 VASP）
python tests/test_altermagnet.py      # P5：几何判据分支 / 容差语义 / 真 fixture 结论 / mc_2d 掩码序参量接线（83 项）
python tests/test_crect_cant.py       # P2'：弱铁磁倾斜的闭式 θ 标度 / 倾斜家族自由度审计 / 反解器与产物契约（44 项）
python tests/test_ingest.py           # P4.2：回收通道 ingest_dft.py 的 11 条硬闸门 / 原子性（不过就一件不写）/ --dry-run / --lenient
python tests/check_sdist.py "dist\*.tar.gz"  # 发布闸门：sdist 里该有的 23 项都在、且没混入 __pycache__/_tmp_/_probe_
```

九个测试都**自带直跑入口，不需要 pytest**（本机环境可能没装；`tests/check_sdist.py` 不是测试，是 sdist 内容闸门，需接一个 `.tar.gz` 路径）。测试 fixture 已内置在 `tests/fixtures/`，不依赖本机任何绝对路径；可用环境变量
`MAGNET_FIXTURES` / `MAGNET_FGT_CASES` / `MAGNET_TEST_TMP` 覆盖。

### 脚本级自检（不跑 DFT、秒级、与测试互为交叉核对）

除测试套件外，**每个阶段脚本都自带自检入口**，用来单独核对某一层的契约。两类入口的语义不同：`--selftest` 跑**本脚本自己的**系数/口径/路由不变量，`--demo` 则是**用已知 J 正演再反解**的闭环（应 PASS）。全部自检**只读仓库内 fixture、不联网**（构建器会往 `--outdir` 写一份产物，其余不落盘）。

```powershell
# 数据/口径层（不带位置参数）
python scripts/magnetic_orders.py --regress        # 9 条链路规格层逐位回归（9/9）
python scripts/mag_vectors.py --selftest           # 自旋矢量层：共线/120°/倾斜/螺旋 → 3 分量 MAGMOM（19/19）
python scripts/dft_energy.py --selftest            # OUTCAR 能量口径 + mixed_convention 混用闸门（9/9）
python scripts/dmi_matrix.py --selftest            # P6：D 矢量 → 3×3 反称交换矩阵 + 与 mc_2d._energy 实机对拍（5/5）
python scripts/mc_2d.py --selftest                 # 比热能量口径：_energy 含单离子各向异性 → 与暴力求和一致（3/3；smoke §5b）
python scripts/verify_states.py --selftest         # P0-2：反解前的磁态逐位点核对（符号不一致 / 两态坍缩 / 缺 OUTCAR 三条路径；smoke §5c）
python scripts/mc_fss.py --selftest                # P1-1：FSS 的交点选择（取近全局参照，不取最陡）+ 1/L 外推（4/4；smoke §5f）
python scripts/altermagnet_check.py --selftest     # P5：交变磁体几何判据的分支/容差/候选报告（15/15）

# 构建器（须给 POSCAR + 元素参数；自检在构建完成后才跑，所以要给可写 --outdir）
python scripts/build_square_cases.py tests/fixtures/K2NiF4_POSCAR --element Ni --nonmag K F              --moment 2.0 --spin 1.0 --outdir ../_selfcheck/sq  --selftest
python scripts/build_crect_cases.py  tests/fixtures/CrOCl_POSCAR --element Cr --nonmag Cl O              --moment 3.0 --spin 1.5 --outdir ../_selfcheck/cr  --selftest
python scripts/build_rectangular_cases.py tests/fixtures/CrSBr_POSCAR --element Cr --nonmag S Br         --moment 3.0 --spin 1.5 --outdir ../_selfcheck/rec --selftest
python scripts/build_magnetic_cases_tri.py tests/fixtures/GdI2_POSCAR --element Gd --nonmag I            --moment 7.0 --spin 3.5 --outdir ../_selfcheck/tri --selftest
python scripts/build_cluster_cases.py tests/fixtures/Nb3I8_POSCAR --element Nb --nonmag I                --outdir ../_selfcheck/clu --selftest
python scripts/build_kagome_cases.py tests/fixtures/Fe3Sn_POSCAR --element Fe --nonmag Sn                 --moment 2.0 --spin 1.0 --outdir ../_selfcheck/kag --selftest
python scripts/build_crect_cases.py --selftest      # 免参数快检：用仓库 CrOCl 夹具 + cant 家族样本角 8/16（P2′）

# DFT 输入组装（自检不需要 --root/--element，内部短路）
python scripts/assemble_dft.py --selftest          # 共线/非共线/--soc/--all-orders 四条路径（18/18）
python scripts/assemble_dft_tri.py --selftest      # 三角薄适配层：转调结果与 INCAR 内容（10/10）
python scripts/k_convergence.py --selftest         # P4.2：k 点收敛抽查的网格解析/只换 KPOINTS/闸门/判定（29/29）

# 文档 / 发布闸门
python scripts/readme_to_obsidian.py --check       # README ↔ README.obsidian.md（Obsidian 阅读版备份）是否同步；--selftest 11/11

# 反解器（--demo 需一个 cases 目录；用上面构建器刚生成的目录即可）
python scripts/solve_J_square.py ../_selfcheck/sq --demo --ref la2cuo4
python scripts/solve_J_crect.py  ../_selfcheck/cr --demo --synthetic 5 -0.5 0.2   # CrOCl 文献 J 未核实 ⇒ 必须给 --ref 或 --synthetic
python scripts/solve_J_crect_cant.py ../_selfcheck/cr --demo --synthetic 5 -0.5 0.2        # P2'：倾斜家族 + 闭式 θ 标度（同一目录，反解器不同）
python scripts/solve_J_crect_cant.py ../_selfcheck/cr --orders cant --audit-only            # 只看自由度审计表（不需要 OUTCAR）
python run_2d_magnet.py build --model crocl --poscar tests/fixtures/CrOCl_POSCAR --outdir ../_selfcheck/crocl --cant-theta "8 16"   # P2'：driver 透传倾斜家族（solve 端会自动改派反解器）
python scripts/solve_J_rect.py   ../_selfcheck/rec --demo
python scripts/solve_J_tri.py    ../_selfcheck/tri --demo
python scripts/solve_cluster_J.py ../_selfcheck/clu --demo
python scripts/solve_J3.py       ../_selfcheck/clu --demo                        # 蜂窝：需 fm/neel/zigzag/stripy 四态
python scripts/solve_J_kagome.py ../_selfcheck/kag --demo --synthetic 3 -0.4 0.15 # kagome：Fe3Sn 文献 J 未核实 ⇒ 必须给 --synthetic
```

> `solve` 子命令在反解前**自动**跑 `scripts/verify_states.py`：逐位点比 `magnetization (x)` 的 `tot` 列与 `magmom.txt` 的目标符号，并检测「两个磁态坍缩成同一个态」——后者两态能量相同 ⇒ 残差 ≈ 0 ⇒ `j_result.confidence_of` 会判 **high**，那个高置信度是假的（它只看 residual/cond，看不到磁矩）。不通过就 exit 1、**不写** `J_result.json`（`states_verification.json` 总会留下）；确要看数只能 `--allow-unverified-states`，此时 `J_result.json` 带 `states_verification.verdict = FAIL` + `note`。非共线 / SOC 体系判 `SKIP`（矢量和不能这样逐位点比）。
>
> `--outdir` 只是「一个可写的产物目录」，写哪都行（`../_selfcheck/*` / `%TEMP%` / 任意工作目录）；这里放在仓库上一级是为了让**自检不污染仓库树**。若把仓库本身放在会话工作区里、而工作区只允许写它自己的顶层时，写进仓库内的新目录可能被拒——换个仓库外的 `--outdir` 即可，与流水线本身无关。
>
> 另两个环境细节：① 手动跑脚本时若控制台是 GBK（`chcp 936`），先设 `PYTHONIOENCODING=utf-8`，否则文献注记里的 `⇒` 会触发 `UnicodeEncodeError`（driver 自身已放宽，且它派生的每个子进程都由 `child_env()` 注入 `PYTHONIOENCODING=utf-8`，见局限 13）；② 构建器的 `--selftest` 不是「只自检」，它**先真建一遍超胞**再跑断言，因此必须给 `--outdir`。`doctor` 只做环境探测，读不到上面这些断言结果。

## 已知局限

1. 结构下载走 Materials Project（`--provider mp`，需 key，会剥单层）或公开 provider 链（`cod/oqmd/aflow/jarvis`，无 key，返回体相）。部分材料所有库都没有（如 Fe₃GeTe₂ 的 MP 三元空间为 0 条），需手动建层用 `--poscar`。
2. **磁性格点分类是纯几何判据**（近邻 < `--cluster-cutoff` 就并簇）。连通的多亚晶格铁磁网络（FGT 类）与真分子簇（Nb₃I₈ 类）签名可能相同，会被误判——driver 会打 stderr 警告但**不改路由**，请人工确认。
3. `--J1` 为正向 MC 的必需输入；未核实文献条目（`literature_J.json` 里 `null`）会明确报错而非硬跑。
4. MC 的 T_C/T_N 取比热峰 `argmax(Cv)`；峰落在温度窗口边界时**不再把边界值当结果**，而是显式「未分辨」并自动缩放重扫（`--mc-refine-*`）。
5. 双亚晶格反解：共线、只到 J4、统一 S=1、结构需先弛豫、MAE 未纳入、k 网格 3×3×1 偏粗。
6. DFT 反解要求 meV 级总能差；金属体系建议 `ISMEAR=1/SIGMA=0.2` 且 5 个磁态设置完全一致。
7. 运行目录按 CWD 决定（`./<formula>_pipeline`）；建议在仓库根运行。
8. **非共线 / SOC 只做到「把 VASP 输入写对」**（P3 已完成：driver `--noncollinear` / `--soc` / `--saxis`）：`assemble_dft.py` 现在能写出 `LNONCOLLINEAR = .TRUE.` + `LSORBIT = .TRUE.` + `GGA_COMPAT = .FALSE.` + 3 分量 MAGMOM（角度表在 `scripts/mag_vectors.py`：共线 / 面内 120° / ±θ 倾斜 / 螺旋），三角路径退化成薄适配层、与共线路径共用同一份 INCAR 渲染；能量口径统一到 `scripts/dft_energy.py` 的 `free energy TOTEN`（取最后一次，缺则报错而不是改用 `energy(sigma->0)`；混用会被 `mixed_convention()` 抓住）。**需要 SOC 的磁序（cant/spiral/tilt…）没开 `--soc` 会被硬拒**，不会静默按共线算一个物理上不成立的能量。**但本仓库从未真跑过非共线 DFT**：验收只到「INCAR/MAGMOM 内容 + 快照回归 + mock-OUTCAR 闭环回解 J」（smoke §21），**没有任何真实非共线能量或 T_N 验证**。因此 kagome 的 **q=0/√3×√3 简并**只在**系数层**被解开（P4.1 用含 √3×√3 的序把设计矩阵 rank 从 2 提到 4，见局限 14），CrOCl 的弱铁磁倾斜只到 P2′ 的**系数层 + 反解器**（θ 是输入参数，见局限 17），而**两者的真实 J/T_N 都还没有**（kagome 要 P4.2 的真实非共线 DFT）；P3 只是把前置能力铺好；DMI 的反对称项在 MC 侧已就位、并且已能用**文献 D 值**直接跑格点模型（P6 的 `scripts/dmi_matrix.py`，边界见局限 15），但 D 与 J 的**提取侧**仍需自旋螺旋法。
9. MC 引擎 `mc_2d.py` 的数据结构**已经支持**逐子晶格各向异性与任意 3×3 交换矩阵（含 DMI 的反对称项），瓶颈在 J 提取侧；`poscar_to_mc3` 目前只有一个全局 S²，逐壳 S_eff 未实现（多个不同 `--S` 会被直接拒绝而非静默算错）。
10. **方形晶格（K₂NiF₄/K₂CuF₄/La₂CuO₄ 型）是 3D 氧化物、不可解理**：单层切面丢掉了层间超交换，MC 的 T_C/T_N 只能当量级检查（La₂CuO₄ 实测 T_N = 325 K 是 3D 驱动，2D 单层有限尺寸 χ 峰更低）。方形**必须**用显式平移搜索（`shift_search`）且**必须**含周期 4 的 `p4_a`，否则 J3 会被低估 2 倍 / 与 e0 不可分离。
11. **C 心矩形（CrOCl 类）有三个静默陷阱**：① 两个 Cr 只差 3.03 Å（< 簇判据 3.3 Å）⇒ 旧版分类会当分子簇（现返回 `rect_centered`；方形网并簇返回 `sq_centered` 并**显式报错**）；② J1=3.0347 Å 与 J2=3.1980 Å 只差 0.1632 Å ⇒ `tol=0.2` 会**并壳且不报错**（`build_crect_cases.py` 默认 `tol=0.05`）；③ C 心键连的是胞内两个 basis ⇒ basis 符号对 (n0,n1)-奇偶序不可见，母序与其 basis 翻倍版同选必使 det=0（脚本硬拒）。**两把尺子**：J1 的 3D 键长 3.0347 Å 与 `poscar_to_mc3.py` 的 2D 投影 2.5064 Å 都对但**不可混用** ⇒ driver 用跨口径不变量 `per_site_bonds`(4/2/2) 校验 J 绑定，不一致直接退出。**CrOCl 的基态是弱铁磁倾斜（需 SOC/DMI），本仓库共线 Heisenberg 层不覆盖**；其 J 与 T_N **均未核实**（`literature_J.json` 该条 J 全 null）⇒ `mc`/`report` 必须显式给 `--J1/--J2/--J3`。**非磁元素占位与化学计量有硬校验**：`--nonmag` 漏写胞内元素（如 `O Cl`）或写了胞里没有的元素都直接报错，每胞计量必须为整数并落进 `stoichiometry_per_cell` / `formula_per_cell`（CrOCl = Cr2O2Cl2）——CrOCl 的原胞以 `Cl` 开头，**MAGMOM 必须按元素分组顺序排**，否则会静默错位（已修）。
12. **超胞图的渲染依赖**：首选 matplotlib，但 matplotlib 只 import 就需要 `pyparsing` / `cycler` / `contourpy` / `fontTools` / `kiwisolver`；本机这几个都没装且**离线**（pypi 不通、pip 缓存里也没有对应 wheel）⇒ `import matplotlib` 必失败、任何 PNG 都出不来。所以通用超胞图有第二条等价路径 `scripts/render_supercell_pil.py`（纯 Pillow，同一 CLI / 输入 / 输出，几何与文本规则复用 `render_supercell.py`），driver 用 `have_matplotlib()` 现场二选一。簇磁体、三角两条渲染器与 MC 曲线图（`05_mc/mc_mt.png`）同样各有一条纯 Pillow 等价路径（`render_cluster_orders_pil.py` / `render_tri_orders_pil.py` / `mc_2d.py` 内的 `_plot_mt_png`，三条共用原语 `scripts/pil_fig.py`）⇒ 只要装了 `Pillow`，这四条渲染路径都能出图；**另外两条**（有限尺寸标度图 `mc_fss.py` 内的 `_plot_fss_png`、FGT 双亚晶格磁序图 `build_fgt_cases.py` 内的 `_plot_fgt_png`）也有各自的纯 Pillow 等价实现并同样自动二选一 ⇒ **六条渲染路径全部与 matplotlib 解耦**。所有 Pillow 版共用同一份原语 `scripts/pil_fig.py`（磁序图用点/箭头层，曲线图用 `chart_*` 直角坐标系层）。

13. **控制台编码**：`stdout` 默认 `errors="strict"`，而文献注记/局限说明里有 `⇒` 这类 GBK 编不出的字符 ⇒ 在 GBK 控制台（`chcp 936`）跑 `mc`/`report` 会在打印文献来源时 `UnicodeEncodeError` 直接失败。driver 现已对 `stdout`/`stderr` 调 `reconfigure(errors="replace")`（只放宽错误处理，不改控制台编码），并且**所有子进程都经 `child_env()` 注入 `PYTHONIOENCODING=utf-8`**（子脚本因此不再依赖控制台编码，也不会因打印 `⇒` 崩掉）；想彻底统一，自己在 shell 里设 `PYTHONIOENCODING=utf-8` 亦可。

14. **kagome（Fe₃Sn 类，3 磁位/胞）只做到「系数层」**（P4.1；P4.2 需要真实非共线 DFT）：`detect_lattice` 现在用**枚举周期像后的真实配位数**把 kagome (4,4,4) 与 Fe₃GeTe₂ (4,6,4)、Nb₃I₈ (2,2,2) 分开（只看胞内原子会得到 z=2，与三聚体混为一谈）；磁序枚举含 q=0 一族与 √3×√3 调制。但要注意三点：① **J1 与 J1′ 被并成一条 `J1`**（同距离、不同亚晶格对；3 壳模型取平均，`order_spec` 注释与 `solve_J_kagome.py` 的 `notes` 都写明了），**不许当成「解出了 J1 与 J1′」**；② `--model fe3sn` 的晶格常数/`S=1.0`/`moment=2.0` 是**理想化占位值**，Fe₃Sn 的文献 J **不在** `literature_J.json`（预设刻意写 `literature=None`）⇒ `mc`/`report` 必须显式给 `--J1/--J2/--J3`，J 与 T_N 均无真实验证；③ **q=0 一族 (fm/q0_120/q0_ud/q0_duu) 对 J1/J2 简并**：系数满足 `c1 == c2` ⇒ 那一族的设计矩阵 **rank 只有 2**（`cond≈2.3e34`，J1/J2 不可分离），必须用含 √3×√3 的组合（builder 自带 `selected = fm/q0_120/sq3_updntri/sq3_alt`，rank=4、`cond=98.51`）；换序后**必须重算 rank/cond** ——「`cond` 没超阈值」不等于「三个壳都分得开」。

15. **DMI / 磁性拓扑绝缘体（MTI）只做判别，不进计算主链**（P6；判据、拒绝模板与可行性备忘见 [DMI_MTI_MEMO.md](DMI_MTI_MEMO.md)）：本仓库**不提取 D**（需要 SOC 四态法或自旋螺旋法的第一性计算 + 3 分量 D 矢量）、**不做能带拓扑**（Z₂ / 表面态），MC 侧也没有外场与拓扑荷 ⇒ 报不出 D 分量、`|D|/J`、skyrmion 半径、拓扑荷 `Q`、B–T 相图、层间 `J_c`、Z₂。**唯一可执行的例外**：把**文献 D 值**换算成逐壳 3×3 反称交换矩阵塞进格点模型（`scripts/dmi_matrix.py`：`J_D = -[D]×` 对应 `H_DM = +D·(S_i×S_j)`，矩阵已吸收 `S²`，产物直接喂 `poscar_to_mc3.py --shell-def`；`--selftest` 5/5，含与 `mc_2d._energy` 的实机对拍）——这属于「拿别人的 D 做格点模拟」，**不许说成「本流水线算了 DMI」**，且只能看**零场基态能**。MTI 侧同理：只能给「属候选类 + 本仓库未覆盖其关键物理」，不做层分辨 A 型 AFM 与层间耦合。

16. **交变磁体（altermagnet）只做几何判别 + 标量序参量**（P5）：`scripts/altermagnet_check.py` 只回答一个**可判定的几何问题**——两套磁性亚晶格的连接操作是平移、反演还是旋转：平移连接 ⇒ 常规共线反铁磁（劈裂被对称性禁止）；反演连接（且两套亚晶格共面）⇒ 宇称禁止自旋劈裂，也不是交变磁体；旋转/镜面连接 ⇒ **几何候选**（`--selftest` 15 项，判据容差以 **Å** 计而不是分数坐标元组，抖动 3e-6 的坐标仍能认出对称操作）。**判据到此为止**：是否真有 k 分辨自旋劈裂**要能带计算**，本仓库不做能带 ⇒ 输出里不许出现「这是交变磁体 / 带隙里有劈裂 / 自旋劈裂大小」这类结论。MC 侧新增逐子晶格 Néel 矢量序参量 `mc_2d.py --order sublattice:<mask>`（driver 对应 `--mc-order`，`mask = {0}` 即老的交错序参量；产物仍是 4 列 `mc_mt.csv` ⇒ `summarize_mc`/`mc_fss` 的列契约定不动），但它只给**标量 |L|** 的 Binder U4 交点，**不区分两套亚晶格的 Néel 矢量方向** ⇒ 交变磁体的序参量判据仍属未覆盖。三个已知坑：`stag_signs` 在 `num_basis >= 3` 时把 basis 1..n−1 并成一组（`--order staggered` 静默退化，现已打 warn）；掩码必须留至少一个子晶格在 B 组（否则 `mbx/cb` 除零，脚本硬报错）；`mc_fss.py` 只支持 `auto/staggered/total/sublattice`，没有 `cv`（U4 是它的判据）；它的 `auto` 与 `mc_2d.py` 同源（按 T=0 候选基态能量），遇到条纹基态会**直接报错**并让你显式给 `sublattice`，不会拿错的序参量硬算 U4 交点。真 fixture 的判别结论：**MnPS3 → 几何候选**（镜面连接，残差 0.012 mÅ）、**CrOCl → 宇称禁止**（反演连接，与局限 11 的「基态是弱铁磁倾斜」一致）、**Fe3Sn（kagome）→ 不适用**（3 套亚晶格）。

17. **CrOCl 弱铁磁倾斜（P2′）只做到「系数层 + 反解器」**：`scripts/solve_J_crect_cant.py` 与 `tests/test_crect_cant.py`（44 项）把倾斜家族的系数审计 + 反解跑通，物理内核是一条**闭式**关系——同一母序在倾斜 θ 后的系数等于「共线母序系数 + sin²θ × 一个与 θ 无关的基向量」的合成，所以倾斜构型的系数是**可解析验证**的。但边界同样明确：① **θ 是输入参数，不是自洽解出来的**——闭式给的是「给定 θ 时系数是多少」，不是「基态 θ 是多少」，真实 θ 仍需 SOC 总能量极小（本仓库没做，也没有真实 SOC 能量）；② **一次倾斜 SCF 定不出唯一的 (J, D)**：单壳闭式只有一条冗余的自由度，`--orders cant` 因此取「倾斜家族 + 张成所需的 fm/母序行」而不是「倾斜行自己」（`--audit-only` 会打印每行的系数与自由度）；③ 倾斜序属 `mag_vectors.SOC_ORDERS` ⇒ 装配**必须** `--soc`，否则被硬拒（不会静默按共线算）；④ CrOCl 的文献 J **未核实**（`literature_J.json` 全 null）⇒ 真值路径必须 `--ref` 或 `--synthetic`；⑤ **没有任何真实非共线/SOC 能量**，验收只到「合成能量闭环 + 与独立实现的逐位对拍」；⑥ **覆盖面**：倾斜族只在 4×4×1 超胞上验证（奇数超胞的 cant 路径未跑），θ 与 φ 的取值域由构建器守卫（`|θ| < 90`、`φ ∈ [0, 360)`、`--cant-base fm` 硬拒），`--lsorbit` 接线由 smoke §26 端到端钉住（`assemble_dft` 缺 `--lsorbit` 对 cant 序硬拒；给了则 `runs/<cant 序>/INCAR` 落 `LNONCOLLINEAR = .TRUE.` + `LSORBIT` + `GGA_COMPAT = .FALSE.` 且无 `ISPIN`）；⑦ **driver 已接通，但只是「透传 + 改派」**：`run_2d_magnet.py build --model crocl --cant-theta "8 16"`（可加 `--cant-base` / `--cant-azimuth`）会把倾斜家族交给 `build_crect_cases.py`，`solve` 端一旦在 `cases/magnetic_orders.json` 里看到 `cant_orders` 就**自动改派** `solve_J_crect_cant.py --orders cant`（不静默按共线解）；三条守卫是硬的——非 C 心矩形给 `--cant-theta`、`--cant-base fm`、以及「有了倾斜序却没 `--soc`」（后者在 `dft` 装配阶段硬拒，`build`/`solve` 只警告）；smoke §26 用「驱动 build + 合成 OUTCAR + 驱动 solve」端到端钉住（J 回解误差 < 5e-4 meV）。但**路由不等于物理解**：θ 仍是你给的输入，`solve` 只是把同一组 J 反解出来并顺手做自洽校验（倾斜行不增加秩）；⑧ **易轴是假设不是结论**：倾斜按「面内分量 + z 分量」写 MAGMOM（`SAXIS` 缺省即 `0 0 1`），各向异性大小与易轴方向都**没有**算过，φ 由几何定不出来（`--cant-azimuth` 只是取值，不是解）。

18. **kagome 的物理层（P4.2）只有操作手册 + 装配侧工具，没有任何真实 J / T_N**：[P4_2_KAGOME_RUNBOOK.md](P4_2_KAGOME_RUNBOOK.md) 把「上真实机时」这件事写成可执行步骤——入口条件（真实机时 + **已核实**的文献锚点 + 参数表填满）、装配与已知坑（`ISYM` 默认只写在 `INCAR.relax`、static 不写，但可用 `--isym` 显式覆盖；`NPAR` 默认 4，可用 `--npar` 覆盖；`KPOINTS` 默认 `3 3 1` 且表头已修正（不再写 2×2×1），可用 `--kpoints 6 6 1` 覆盖；`POTCAR` **不生成**，只给 `POTCAR.spec` 的元素顺序；理想化胞未弛豫）、SCF 四条判据与失败协议、回收通道**已自动化**（`ingest_dft.py` 的 11 条硬闸门 + `dft_handoff.py` 的交接包，见手册 §5.0/§5.1；此前「没有从 `03_dft/runs/` 的自动 ingest」的说法已不成立）、反解可信度与「残差不为零先怀疑什么」、§8 的**纯估算**机时（替换成实测前不许当预算）、§9 与原计划的五项偏差（含 P4.2 残余工作量 **≈4–9 人天**的补估）。k 点收敛抽查同样只有**工具**（`python run_2d_magnet.py kconv --outdir <proj> --k-grids 3,4,6,8`，或直调 `scripts/k_convergence.py`）：它把 `runs/<序>/` 复制成 `<序>@<k网格>` 变体（**只换 KPOINTS**，INCAR 逐字节不变）、OUTCAR 放回后判「总能量对 k 网格收敛到哪个网格够了」；它**不代提交作业**，而且总能量收敛只是 J 可信的**必要**条件（J 是能量**差**）——`k_convergence.json` 里因此带着这条免责声明。**手册与工具都不是进度**：路线图上的 P4.2 仍未勾选，`literature_J.json` 里 Fe3Sn **条目已有**（铁磁基态 + T_C = 705 K + 结构/各向异性锚点，见该条 `_note`），但 **J1/J2/J3 按设计仍为 `null`**（原文自承交换常数尚未确立）⇒ `solve_J_kagome.py --ref Fe3Sn` 会被明确拒绝；仓库里没有任何真实 kagome OUTCAR。
19. **反解前的磁态核对只覆盖共线态（P0-2）**：`scripts/verify_states.py` 逐位点比 `magnetization (x)` 的 `tot` 列与 `magmom.txt` 的目标符号（两种 `magmom.txt` 格式都认：collinear 的 `MAGMOM = …` 与 FGT 表格式的 `sigma=`），并检测「两个磁态坍缩成同一个态」；driver 的 `solve` 阶段已把它当**硬闸门**（fail-close，不通过就 exit 1、不写 `J_result.json`）。边界：① 非共线 / SOC 体系判 `SKIP`（矢量和不能逐位点比），即 P3 / P2′ / P4.2 仍只有 `ingest_dft.py` 的净磁化软闸门，**不等于逐位点核对过**；② 它核对的是**磁态身份**，不是 J 的物理正确性（残差、条件数、能量口径各自另有闸门）；③ `--allow-unverified-states` 可放行，但 `J_result.json` 会留下 `states_verification.verdict = FAIL` + `note`，不会静默变成「已验证」。

20. **MC 采样没有诊断量（P1-4，未做）**：`scripts/mc_2d.py` 输出 T / E / Cv / U4 / |m| / χ，但**不给 Metropolis 接受率、不给 τ_int、不做 binning** ⇒ 单条马尔可夫链自身的统计误差**算不出来**；`--therm 2500 --meas 6000` 是经验默认值，没有自洽依据。目前只能靠 `--mc-seeds` 的多初态离散间接估（峰位/交点的 ± 也是这个来源）。要更硬的误差棒得把 `--therm/--meas` 加倍看 T_C 是否移动——这一步**没有自动化**。

21. **能力边界没有机器可读清单（P1-5，未做）**：本文的局限 1–21 与 [dsh/2d-magnet/SKILL.md](dsh/2d-magnet/SKILL.md) 的「本 skill 不做什么」都是自然语言，agent 靠读文档遵守；没有 `capabilities.json` 这类可被程序断言的清单 ⇒ 「报告里不许出现 D 值 / 能带 / Z₂ / 交变磁体结论」这些承诺**无法在 CI 里机器验证**，只能人审。

## 路线图

- [x] 硬编码路径清理 + 测试 fixture 内置 + 打包基线（LICENSE/README/requirements/pyproject）
- [x] `doctor` 环境自检 → 子命令 + error catalog（`KNOWN_ISSUES.md`）
- [x] 元素参数表（Mn/Cr/Ni/Co/V 兜底 S/moment/U）+ `--lattice` 显式化
- [x] `--model` 命名预设（MnPS₃/CrI₃/CrBr₃/Fe₃GeTe₂ 一键设定，文献 J 仅参考）
- [x] Binder U₄ 累积量 + `--mc-seeds` 误差棒
- [x] 自动多尺寸有限尺寸标度（U₄ 交点定位 T_C）
- [x] 结构 provider 链（COD/OQMD/AFLOW/JARVIS）
- [x] 子命令化（`fetch/build/dft/mc/solve/render/report/doctor`）
- [x] `scripts/` 包化，支持非可编辑 `pip install .`
- [x] 矩形晶格分支（CrSBr/CrSCl 类，`detect_lattice` 查 a⊥b 且 a≠b → 4 共线态 + `solve_J_rect`）
- [x] **P0.1** 公共规格层 `magnetic_orders.py`（壳/磁序声明式化，五种晶格系数逐位回归）
- [x] **P0.2** 统一产物 `J_result.json` + 信任度 / 残差 / 条件数闸门
- [x] **P0.3** 产物自描述块 `order_spec[]` / `model_spec`
- [x] **P0.4** `poscar_to_mc3` 扩展 `--jfile` / `--shell-def` / `--element X Y` / `--family`（旧路径逐位不变）
- [x] **P0.5** 规格层与 `poscar_to_mc3` 的契约测试（`smoke_test.py` §16/§17）
- [x] **P1** 方形晶格（K₂NiF₄ / K₂CuF₄ / La₂CuO₄ 类）——`build_square_cases.py` + `solve_J_square.py` + `--lattice square` 路由 + 预设 + §18 回归（J3 必须靠周期 4 的 `p4_a`）
- [x] **P2** C 心矩形 / oxyhalide 家族（CrOCl 类）——`build_crect_cases.py` + `solve_J_crect.py` + `--lattice rect_centered` 路由 + `crocl` 预设 + 规格层 `crect_crocl` 回归 + 跨口径 `per_site_bonds` 校验 + §19
- [x] **P3** 非共线 + SOC 地基——自旋矢量层 `mag_vectors.py` + 统一入口 `assemble_dft.py --noncollinear/--soc/--saxis`（薄适配层 `assemble_dft_tri.py`）+ 能量口径 `dft_energy.py`（六个 `solve_J*` 全部转发）+ §21 快照回归与 mock-OUTCAR 闭环（为 kagome/DMI 与 CrOCl 弱铁磁倾斜铺路；**未跑真实非共线 DFT**，见局限 8）
- [x] **P4.1** 六角 kagome 的**系数层**（纯几何 + 线性代数，不需要 DFT）——`build_kagome_cases.py` + `solve_J_kagome.py` + `detect_lattice` 的**配位数判别**（kagome (4,4,4) vs Fe₃GeTe₂ (4,6,4) vs Nb₃I₈ (2,2,2)）+ **q=0 简并的数值证实**（rank=2、cond≈2.3e34）+ √3×√3 分离（rank=4、cond=98.51）+ `kagome_degenerate` 闸门 + smoke §23 + `tests/test_kagome_cases.py`（68 项）
- [ ] **P4.2** kagome 的**物理层**（Fe₃Sn₂ / Co₃Sn₂S₂ / Mn₃Sn 等）——需要真实非共线 DFT：提交 √3×√3 构型 → 读 `OUTCAR`（`free energy TOTEN`）→ `J_result.json` → MC 求 T_N。**操作手册已就绪**（[P4_2_KAGOME_RUNBOOK.md](P4_2_KAGOME_RUNBOOK.md)：入口条件 / 装配 / SCF 判据 / 回收 / 可信度 / 机时估算 / 五项偏差），**物理层仍未跑**——入口条件是有真实机时 + 已核实的文献锚点（见局限 18）
- [x] **P2′** CrOCl 弱铁磁倾斜（oxyhalide，需 SOC）——**只做系数层 + 反解器**：`build_crect_cases.py` 的倾斜家族 + `scripts/solve_J_crect_cant.py`（闭式 sin²θ 标度、自由度审计 `--audit-only`、`--orders cant` 含 fm/母序张成行）+ 契约测试 `tests/test_crect_cant.py`（44 项）+ smoke §26；**driver 已接通**（`build --cant-theta/--cant-base/--cant-azimuth` 透传 + `solve` 见到 `cant_orders` 自动改派，三条硬守卫）；**θ 是输入参数不是基态解**，且一次倾斜 SCF 定不出唯一 (J, D) ⇒ 不产生「CrOCl 的 J/D」结论（见局限 17）
- [x] **P5** 交变磁体（MnTe / CrSb）——**只做几何判别 + 标量序参量**：`scripts/altermagnet_check.py` 判「两套亚晶格的连接操作是平移 / 反演 / 旋转」，MC 侧新增逐子晶格 Néel 矢量序参量 `mc_2d.py --order sublattice:<mask>`（driver `--mc-order`）+ 契约测试 `tests/test_altermagnet.py`（83 项）；**不做能带 ⇒ 不产生交变磁体结论**（见局限 16）
- [x] **P6** 体相 DMI / 磁性拓扑绝缘体：**只做识别层**（按用户指示跳过 P5 先做 P6）——判据与边界写进 [DMI_MTI_MEMO.md](DMI_MTI_MEMO.md)（一页结论 + G1–G6 逐关复核 + 拒绝模板 + 构型数量与机时估算），`dsh/2d-magnet/SKILL.md` 新增判别节，smoke §24 钉住契约；**不新增 DMI/MTI 计算子链**。另落盘一个受控例外：`scripts/dmi_matrix.py` 把**文献 D 值**换成逐壳 3×3 反称交换矩阵（`--selftest` 5/5），让「含 DMI 的格点模型」今天就能跑 —— 但**没有任何 D 的第一性提取**

> P0–P6 的判定闸门、改动清单与工期见 [NEW_TYPES_PLAN.md](NEW_TYPES_PLAN.md)。

## 发布状态

- **版本**：`0.1.0`（唯一权威声明在 [pyproject.toml](pyproject.toml) 的 `version`；本 README 与 skill 不再各写一份）。`python run_2d_magnet.py --version` 打印 `2d-magnet 0.1.0`，`doctor` 首两行也回显版本与脚本目录——两者都在运行时读取那**一个**声明（已安装时读分发元数据，源码树里读同目录 `pyproject.toml`），代码里没有第二个写死的版本号。
- **许可**：[MIT](LICENSE)（`pyproject.toml` 的 `license` 与之同源，作者署名一致）。
- **P0–P4.1 能力状态**：P0 规格层 / P1 方形 / P2 C 心矩形 / **P3 非共线 + SOC** / **P4.1 kagome 系数层** 均已落地并通过自检与回归（见「路线图」与「脚本级自检」）。**P3 与 P4.1 都是「能力已实现」，不是「物理已验证」**：没有任何真实非共线 DFT 能量，P3 的验收止于 INCAR/MAGMOM 内容 + 快照回归 + mock-OUTCAR 闭环回解 J（smoke §21），P4.1 的验收止于系数矩阵/条件数 + 合成能量回解 + 配位数判据（smoke §23，68 项契约测试），详见局限 8 与局限 14。**P6** 随后落地，但**只是识别层**（判据 + 边界 + 备忘，见 [DMI_MTI_MEMO.md](DMI_MTI_MEMO.md) 与局限 15）——它**不产生 D、不做拓扑判定**，所以不构成任何物理能力；随它落盘的 `scripts/dmi_matrix.py` 是把**文献 D 值**换成格点模型的受控换算器（`--selftest` 5/5），同样不是 DMI 计算能力。**P4.2**（kagome 的真实 J/T_N，需真实非共线 DFT）**只有操作手册**（[P4_2_KAGOME_RUNBOOK.md](P4_2_KAGOME_RUNBOOK.md)），物理层仍未跑；**P2′**（CrOCl 弱铁磁倾斜）已落地为**系数层 + 反解器**（`solve_J_crect_cant.py` + 44 项契约测试）——但 θ 是**输入参数**、一次倾斜 SCF 定不出唯一 (J, D) ⇒ 同样**不是**物理能力（见局限 17）；**P5** 已按用户指示落地为**识别层 + 标量序参量**（几何判据 + `--order sublattice:<mask>`，无能带 ⇒ 不构成物理能力，见局限 16）。
- **可以带发布**：无 CHANGELOG；无 `docs/`（说明集中在 README/WORKFLOW/KNOWN_ISSUES/NEW_TYPES_PLAN 四份，内容自洽）；`templates/` 里放的是扩展骨架（`material_preset.json` / `new_lattice_checklist.md` / `README.md`，**不参与运行**，见该目录自述）；仓库根那份 3 MB 的 `.pptx` 是汇报稿、被 `.gitignore` 的 `*.pptx` 覆盖，且**从未进过版本库**（`git log --all -- "*.pptx"` 为空）⇒ 不存在"要单独开一次提交删掉"这回事（如仍想瘦身可 `git gc --prune=now`）。**已消项**（原列在这里的两条已解决）：CLI 现有 `--version`；FSS 标度图与 FGT 图已有纯 Pillow 兜底（局限 12）。
- **版本控制状态（2026-10）**：P0–P4.1、P2′、P4.2（工具链）、P5 与 P6 的改动**已全部提交**——功能提交 `4839bf8`（`65 files changed, 18063 insertions(+), 438 deletions(-)`，26 个 `M` + 39 个 `A`，无删除），随后是文档状态更新提交 `b371946`，再到发布前总检提交 **`6215c64`**（`发布前总检: examples 样例产物 + CI + 文档/提示修正`，2026-10-16 13:12 +0800，`34 files changed, 1573 insertions(+), 15 deletions(-)`：`examples/` 27 个文件 + `.github/workflows/ci.yml` + 六份文件的文案/提示修正），最后是打包闸门提交 **`84132f6`**（`发布闸门: MANIFEST.in + check_sdist 校验 + CI package 作业`，2026-10-16 13:23 +0800，`5 files changed, 218 insertions(+), 6 deletions(-)`：`MANIFEST.in` + `tests/check_sdist.py` + CI `package` 作业 + smoke §28），以及 Obsidian 阅读版提交 **`d23c243`**（`README 的 Obsidian 阅读版备份 + 同步闸门 (smoke §29)`，2026-10-16 13:36 +0800，`5 files changed, 667 insertions(+), 7 deletions(-)`：`README.obsidian.md` + `scripts/readme_to_obsidian.py`），再到复审第一批修复提交 **`acd127a`**（`vault 友好性 + 复审第一批: 比热能量补单离子项 + 自作用断言 + mc_2d 自检`，2026-10-18 14:01 +0800：`scripts/mc_2d.py` 的比热能量口径 + `--selftest` + smoke §5b，连同 `d23c243` 遗留的 vault 友好性文件共 8 个），以及复审第二批提交 **`add33cd`**（`复审第二批 P0-2: 反解前的磁态硬闸门 (verify_states) + 坍缩检测 + smoke §5c`，2026-10-18 16:00 +0800，`8 files changed`：`scripts/verify_states.py` 新增 + `run_2d_magnet.py` + `tests/smoke_test.py` + `tests/check_sdist.py` + 三份文档）。仓库共 **130 个已跟踪文件**，`add33cd` 提交后工作树是干净的（`git status --porcelain` 无输出；此后只有本轮的 **8 个文件**未提交——`run_2d_magnet.py`、`scripts/mc_2d.py`、`scripts/mc_fss.py`、`tests/smoke_test.py`、`README.md`、`README.obsidian.md`、`KNOWN_ISSUES.md`、`dsh/2d-magnet/SKILL.md`，见下一条）。功能提交含 `templates/`、`DMI_MTI_MEMO.md`、`P4_2_KAGOME_RUNBOOK.md`、`NEW_TYPES_PLAN.md`、P4.1 的 kagome 三件套、P4.2 工具链（`scripts/ingest_dft.py` + `scripts/dft_handoff.py` + `scripts/k_convergence.py`）、P5 的 `scripts/altermagnet_check.py` + `tests/test_altermagnet.py`、P2′ 的 `scripts/solve_J_crect_cant.py` + `tests/test_crect_cant.py`；0.1 期间勘察留下的探针已清理、未进提交。
- **打包闸门（sdist/wheel，已提交 `84132f6`）**：`pyproject.toml` 的 wheel **只装** `run_2d_magnet.py` + `scripts/`（含 `package-data` 的 json），文档、手册、agent skill、样例产物、测试夹具、扩展骨架全靠 [`MANIFEST.in`](MANIFEST.in) 进 sdist——「clone 里有、sdist 里没有」这类发布事故不会被任何运行测试抓到，因此单独立了一道闸门 [`tests/check_sdist.py`](tests/check_sdist.py)（**23 条必须项**用 `endswith` 匹配，并拦下混入 `__pycache__/` / `*.pyc` / `_tmp_` / `_probe_` / `*.egg-info/` 的脏物；smoke §28 用**合成 tar.gz** 反向验证闸门本身会拦缺项与脏物），CI 的 `package` 作业跑 `python -m build` → 在别处只装 wheel 跑一次 smoke（`2d-magnet --version` + 读 `scripts/literature_J.json`）→ `python tests/check_sdist.py "dist\*.tar.gz"`。**诚实声明**：本机没有 `setuptools`（系统 Python 3.12 与 DSH 内置 3.12.14 都没有）且 PyPI 不可达 ⇒ **`python -m build` 没在本地跑过**，这道闸门只有 CI 会实跑；能本地核对的都核了——全部源码在 `ast.parse(feature_version=(3, 9))` 下零语法错（`requires-python = ">=3.9"` 只做到这一层，实际只在 3.12 上跑过）。
- **本轮（8 个文件未提交：复审第三批 P0-3「`--quiet` 把结论级警告吞了」+ 第一批 P1 修正）**：**问题**——子脚本把警告打成 stdout 的普通文本行（全仓 47 处 `[warn]` 里只有 2 处走 stderr），而 `--quiet`（以及封装 `dsh/bin/2d-magnet.cmd` 的默认行为）只转发 stderr、成功时丢弃 stdout ⇒ 「比热 Cv 峰落在温度窗口边界 ⇒ T_N **未分辨**」「过渡带内既无 U4 交点也无 collapse 点」「cond(A) > 1e3 病态」这些**跑完了但结论无效**的告警被静默吞掉，退出码仍是 0。这正是最需要被看见的一类输出。
  - [`run_2d_magnet.py`](run_2d_magnet.py)：新增 `WARN_MARKERS` / `warning_lines()` / `append_log()`；`run()` 在 `--quiet` 下把捕获的 stdout 里的**结论级警告行**重打到 stderr（普通进度行不刷屏）；新增全局 `--log PATH`，把每个阶段子进程的完整 stdout/stderr 追加落盘（`--quiet` 下同样有效，非 `--quiet` 时控制台输出照旧）；`try_mp_fetch()` 与 provider 链的 `fetch_structure.py` 这两个「捕获后就丢」的调用点也补上警告转发与落盘。
  - [`tests/smoke_test.py`](tests/smoke_test.py) 新增 §5d（4 项）：用一个假的阶段脚本验证 `--quiet` 下四类警告都到 stderr、普通行不被误刷、`--log` 落全量、`warning_lines` 判据 5/6。
  - 实测（`--quiet --log`，La2CuO4 8×8、5 个温度点）：stderr 里确实出现 `>>> 磁化率峰处 T ≈ 200.00 K [warn] 落在窗口边界 → 未分辨` 与 `[warn] 峰判据落在温度窗口边界: …该值不能当作 T_C/T_N`；`--log` 文件 133 行、逐阶段分块（`===== reorient.py (exit 0) =====` + `--- stdout ---`）。
  - **P1-2（`--order auto` 不再靠 J1 符号猜）**：[`scripts/mc_2d.py`](scripts/mc_2d.py) 新增 `_candidate_orders()` / `_resolve_auto()` / `resolve_auto_order()`：先按格点几何枚举受支持的候选基态（`fm` / `neel` / `stripe_x` / `stripe_y`；`neel` 只在能二着色时出现，`stripe_*` 只在 `num_basis == 1`），用 `_energy` 比 **T=0 能量**再决定序参量。旧实现只看 `J1`（或各壳 xx 分量之和）的符号 ⇒ `J1>0 + J2<0` 的**条纹反铁磁**会被标成 Néel，于是「交错序参量的峰」根本不是 T_N，而模型照常打印一个温度。现在：条纹基态会打 `>>> [warn] auto 判定的基态是 stripe_x …已自动改用 --order cv (比热峰)` 并真的切到比热路径；候选近简并（差 < 1e-6 meV/格点）也会警告并要求显式给 `--order`。同一个判定被 [`scripts/mc_fss.py`](scripts/mc_fss.py) 复用（它的 FSS 没有比热路径 ⇒ 条纹基态直接报错让人显式给 `sublattice`，而不是拿错的序参量算 U4 交点）。`run_2d_magnet.py` 的结尾摘要（`peak_report` / `summarize_mc`）也改成**按 CSV 表头**判第三列是 `Cv` 还是 `Susceptibility`——否则 `--order auto` 在 `mc_2d` 内部回退到比热后，driver 会把 Cv 峰报成「χ 峰」（同一条静默错标签，代码位置在 `run_2d_magnet.py` 的 `peak_report`）。
  - **P1-1（FSS 的 ± 语义 + 交点选择偏差）**：[`scripts/mc_fss.py`](scripts/mc_fss.py) 的 `±` 改为明确标注「尺寸对交点的**离散度**，不是统计标准误」（这些交点是相关估计，不是独立样本）；一对尺寸出现多个候选交点时不再「取最陡」（那是选择偏差，会系统性偏向噪声陡坡），改取离全局 collapse 估计最近的（`_pick_crossing`）；≥3 个交点时再给一条 `T_c(∞) = T_c + a/L` 的最小二乘外推，并提示 `|a|` 大就是「还没进标度区、外推不可信」；两条不变量都进了 `mc_fss.py --selftest`（4/4）。
  - **P1-3（降级开关写进 agent 规程）**：[`dsh/2d-magnet/SKILL.md`](dsh/2d-magnet/SKILL.md) 新增「降级开关（真实存在，但不许用来绕过闸门）」表：`--allow-unverified-states` / `--ingest-force` / `--ingest-lenient <名>` / `--no-s2` 各松掉了什么、必须在报告里写明用没用什么。原来自我学习协议里三处「回归 → 提交」也改成「回归 → 把改动清单交给用户提交（不要自己 `git commit`）」（P1-6）。
  - [`tests/smoke_test.py`](tests/smoke_test.py) 再新增 §5e（6 项：三种合成模型的 auto 判定 + 条纹基态自动退到比热峰、且警告在 stderr + driver 摘要按 CSV 表头判「Cv 峰 / χ 峰」两种布局）与 §5f（6 项：多候选交点不取最陡、1/L 外推能还原合成直线、报告里不再出现「交点标准差」、`--selftest` 4 项全 PASS、条纹基态在 FSS 的 auto 下直接报错、Pillow 兜底图照常出）。
  - 文档同步：本 README 的全局开关表与用法段、[`dsh/2d-magnet/SKILL.md`](dsh/2d-magnet/SKILL.md) 的「命令封装」段、[`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) 新增 `order_auto_ambiguous`（auto 判不出/近简并/条纹基态时怎么办），另加局限 20（MC 无接受率/τ_int/binning）与局限 21（能力边界没有机器可读清单）两条**未做**项。收尾提交：
  ```powershell
  cd 2d_magnet_pipeline
  git add -A
  git commit -m "复审第三批+第一批 P1: --quiet 不再吞结论级警告 + --log + auto 序参量按 T=0 能量判 + FSS 离散度语义"
  ```
  注意 `git add`/`git commit` 必须在**普通 shell** 里跑：本仓的开发沙箱会拒绝子进程写 `.git/`（`fatal: Unable to create '.git/index.lock': Permission denied`，连 `git add -n` 也一样）。**只要没提交，`git clone` 出来的分发包就缺掉方形/C 心矩形/非共线/kagome 整条链**——这正是「安装与依赖声明」之外最容易被忽略的发布阻塞项。

## License

[MIT](LICENSE)
