# 2d-magnet-pipeline

把单层结构变成交换常数 `J1/J2/J3`，再用蒙特卡洛给出 `T_C` / `T_N` 的估计——每个数字都过一道 fail-closed 闸门。

输入是一个 `POSCAR`（Materials Project、COD / OQMD / AFLOW / JARVIS 公共链，或你自己的文件）。流水线识别磁性格子、建超胞与磁序集、装配 VASP 输入、回读 `OUTCAR`、反解交换常数，然后跑经典 Heisenberg 蒙特卡洛。到 J 有两条路：

- **文献 J 或用户给定的 J → 正向蒙特卡洛**：不需要 VASP；
- **DFT 多磁态能量映射反解 J**：需要 VASP 能量。

支持的晶格共用同一个规格层与同一份产物（`J_result.json`）：honeycomb、三角（120° 非共线）、分子簇、双亚晶格铁磁网络（Fe₃GeTe₂ 类）、矩形（CrSBr 类）、C 心矩形（CrOCl 类）、方形（La₂CuO₄ / K₂NiF₄ 类）、六角 kagome（Fe₃Sn 类）。

做不到的事写在 [`capabilities.json`](capabilities.json) 里（`python scripts/capabilities.py --show` 打印可读版）：

- 不提取 DMI 的 `D` 值（`D` 只能作为输入给出）——见 [`DMI_MTI_MEMO.md`](DMI_MTI_MEMO.md) 与 `scripts/dmi_matrix.py`；
- 不给 k 分辨的能带劈裂，也不做 `Z₂` 分类；
- 不对交变磁体下序参量结论：`scripts/altermagnet_check.py` 只做几何判别（标量序参量）；
- 不代提交 VASP 作业——流水线只导出输入、回读输出。

## 工作流

```text
POSCAR    -> fetch     重新取向、真空层                      01_structure/POSCAR
          -> build     晶格识别、超胞、磁序集                02_supercell/cases/
          -> dft       INCAR / KPOINTS / MAGMOM / POTCAR.spec 03_dft/、runs/<序>/
          [ VASP 在流水线之外运行 ]
          -> ingest    回读 OUTCAR，11 道 fail-closed 闸门    cases/<序>/OUTCAR
          -> solve     反解 J1 J2 J3 (meV)，磁态逐位点核对    02_supercell/cases/J_result.json
          -> mc        Metropolis + Binder U4 -> T_C / T_N    05_mc/mc_mt.csv、mc_diag.json
          -> render    超胞与磁序图                           figures/supercell_orders.png
          -> report    出处卡片                               workflow.md
```

`kconv` 只跑 K 点收敛（`scripts/k_convergence.py`；`--kpoints` 钉住网格）。每个阶段都写下输入、出处与退出状态，不依赖复现不了的现场。

三条规则保证数字可信。不猜参数：`--S/--moment/--U` 只对 Mn/Cr/Ni/Co/V/Gd 有内置表——Fe 被刻意排除，因为它的有效自旋本身有歧义（S = 2 / 2.5 / 1），driver 会停下来问而不是给个值；`literature_J.json` 里 `null` 的条目会直接报错，不会静默跑下去。闸门 fail-closed：磁态核对不过、设计矩阵退化、残差超阈、`OUTCAR` 格式不对，都会中止该阶段；`--quiet` 也不会吞掉说明它的警告。

## 安装

```powershell
py -m venv .venv
.venv\Scripts\activate
py -m pip install -r requirements.txt
py run_2d_magnet.py doctor
```

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python run_2d_magnet.py doctor
```

依赖：Python ≥ 3.9、`numpy`、`numba`、`matplotlib`、`requests`。`matplotlib` 可选——出图会退回 Pillow（`render-fallback` extra），再不行退回文字摘要。`doctor` 报告缺什么，并指出 [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) 里对应的错误键。

代码也可当库用：`import run_2d_magnet`；`pip install .` 之后还有 `2d-magnet` 命令行入口。

## 快速开始

`examples/` 里有两份不用 VASP 就能复现的完整运行（文献 J 进，`T_C` 出）：

```powershell
py run_2d_magnet.py --model la2cuo4 --poscar tests/fixtures/La2CuO4_POSCAR `
   --outdir C:\work\la2 --skip-dft --mc-size 12 12 --mc-seeds 2
py scripts\solve_J_square.py C:\work\la2\02_supercell\cases --demo
```

`examples/la2cuo4_square/` 是这样跑出来的目录树；`examples/kagome_fe3sn_dryrun/` 是一条走完整 DFT 形状的 kagome 空跑（合成能量）。两者各有 `README.md`，列出产出了什么、刻意不产出什么。

## 子命令

| 子命令 | 作用 |
|---|---|
| `fetch` | 从 Materials Project / COD / OQMD / AFLOW / JARVIS 取结构 → `01_structure/POSCAR` |
| `build` | 晶格识别、超胞、磁序集 → `02_supercell/cases/` |
| `dft` | INCAR / KPOINTS / MAGMOM / POTCAR.spec → `03_dft/`、`runs/<序>/` |
| `ingest` | OUTCAR 校验（11 道闸门）→ `cases/<序>/OUTCAR` |
| `kconv` | K 点收敛，只变 KPOINTS |
| `solve` | 多磁态能量映射 → `J_result.json` |
| `mc` | Metropolis / Binder U4 → `T_C` 或 `T_N` |
| `render` | 超胞与磁序图 |
| `report` | 出处卡片 → `workflow.md` |
| `doctor` | 环境自检，对照 `KNOWN_ISSUES.md` |

预设：`la2cuo4`、`crsbr`、`crocl`、`k2nif4`、`fgt`、`mnps3`、`nio`、`gdi2`。driver 选项（`--model`、`--literature`、`--element`、`--S`、`--moment`、`--U`、`--mc-order`、`--mc-seeds`、`--quiet`、`--log`）用 `py run_2d_magnet.py <子命令> --help` 查看。

## 产物

```text
01_structure/POSCAR              取向与真空层处理后的输入
02_supercell/cases/              每个磁序一个目录：POSCAR、magmom.txt、OUTCAR
02_supercell/cases/J_result.json J 与出处：source、confidence、cond(A)、residual、E0
02_supercell/cases/states_verification.json   逐位点磁矩核对（solve 的前置闸门）
03_dft/、runs/<序>/              INCAR、KPOINTS、POTCAR.spec、MAGMOM
05_mc/mc_mt.csv、mc_mt.png       温度、序参量、磁化率、Binder U4
05_mc/mc_diag.json               接受率、tau_int、n_eff、naive 与 binning SEM
figures/supercell_orders.png     带交换键与磁序集的超胞图
workflow.md                      本次运行的出处卡片
```

`J_result.json` 一定记录 J 是怎么来的（`source` = `literature` / `user` / `dft` / `demo`）、设计矩阵条件数、残差，以及磁态有没有逐位点核对过。没有这些出处的数字不算结果。

## 验证（脚本级自检 / 发布闸门）

```powershell
py run_2d_magnet.py doctor
py scripts\capabilities.py --check
py tests\smoke_test.py
py tests\check_sdist.py "dist\*.tar.gz"
```

契约测试是 `tests/` 下八个文件（`test_magnetic_orders.py`、`test_square_cases.py`、`test_crect_cases.py`、`test_kagome_cases.py`、`test_fgt_e2e.py`、`test_altermagnet.py`、`test_crect_cant.py`、`test_ingest.py`），另有 `scripts/magnetic_orders.py --regress` 守位置基线的回归。脚本级自检：

```text
python scripts/mc_2d.py --selftest            能量口径、单离子各向异性、tau_int
python scripts/mc_fss.py --selftest           交点选取、1/L 外推
python scripts/verify_states.py --selftest    逐位点磁矩、坍缩检测
python scripts/dmi_matrix.py --selftest       DMI 矩阵与交换能对拍
python scripts/k_convergence.py --selftest    K 点规划
python scripts/readme_to_obsidian.py --check  README 备份是否同步
```

`solve` 在反解之前自动跑 `scripts/verify_states.py`：两个磁序塌成同一个态、或某个位点的磁矩与目标矛盾，该阶段就以非零退出，不写 `J_result.json`。

## 已知局限

1. 仓库不带真实 DFT 数据。开箱跑出的数字都来自文献或合成能量，并且如实标注——demo 给的 `T_C`/`T_N` 不是预言。
2. 多磁态映射假设共线态映射到 Heisenberg 模型；非共线与 SOC 情形停在系数层（P3），没有经过验证的能量。
3. `S`、`moment`、`U` 是输入。内置表之外的元素必须显式给值：猜一个会静默给出错的物理，所以 driver 宁可拒绝。
4. 各向异性只通过 `--MAE`（默认 0）与 DMI 矩阵进入，两者都是输入，不是反解出来的量。
5. 蒙特卡洛是经典、单自旋更新的 Metropolis，跑在有限团簇上。`T_C`/`T_N` 是该团簇上的估计，不是热力学极限值。
6. 误差棒用 binning 标准误；`mc_diag.json` 给出 `n_eff`、`tau_int` 与接受率，采样不足时能看出来。峰落在温度窗口边界时，`T` 会标成未分辨，而不是照常给一个数。
7. 有限尺寸标度用全局参照交点 + `1/L` 外推，外推值只是参考；没进标度区之前它不可信。
8. `--order auto` 用 T = 0 候选基态（fm / Néel / 条纹）决定用哪个序参量。条纹基态没有 ±1 交错序参量，流水线会自动退回比热峰并给出说明；近简并的模型请显式给 `--order`。
9. kagome 在系数层（P4.1），物理层（P4.2）等真实非共线 DFT——见 [`P4_2_KAGOME_RUNBOOK.md`](P4_2_KAGOME_RUNBOOK.md)。
10. 交变磁体只做几何判别（P5），标量序参量能到子晶格掩码，但不是完整的非共线描述。
11. 倾斜态（P2′）到系数层：`scripts/solve_J_crect_cant.py` 把 θ 当输入，不预测 θ。
12. DMI 以你提供的 `D` 矩阵进入，流水线不从能量里提取 `D`。
13. ingest 校验的是格式、收敛、结构与口径，不是物理。`ingest_report.md` 干净只说明 OUTCAR 可用，不说明计算是对的。
14. 不计算 k 分辨能带劈裂与 `Z₂`；`kconv` 只收敛 K 点网格。
15. `literature_J.json` 里的条目是起点；条目为 `null` 时 driver 拒绝运行，不会编一个值。
16. 磁态核对覆盖共线态；非共线或 SOC 运行会以 `SKIP` 退出（不是通过）。
17. 环境假设：Python ≥ 3.9 只做到语法层检查，实际回归跑在 3.12；`matplotlib` 与 `requests` 可选，缺失时退到有文档的兜底路径。
18. 发布件由 `scripts/make_release.py` 从工作树打包、`tests/check_sdist.py` 校验；PyPI 式的 `sdist`/`wheel` 构建由 CI 的 `package` 作业执行，不在本机。

## 路线图

`✅ 已完成` = 能力已实现、并有自检与回归覆盖，但**不等于物理已验证**；`⬜ 待机时` = 代码与手册已就绪，等真实 DFT 机时。

- ✅ 已完成 · 仓库卫生：路径、夹具、打包基线（LICENSE / README / requirements / pyproject）
- ✅ 已完成 · `doctor` 自检 + `KNOWN_ISSUES.md` 错误目录
- ✅ 已完成 · 元素参数表（Mn/Cr/Ni/Co/V）与显式 `--lattice`
- ✅ 已完成 · **P0** 共享规格层：`--lattice` + `magnetic_orders.json`
- ✅ 已完成 · **P1** 方形晶格：磁序集、ingest、`J1/J2/J3`、蒙特卡洛
- ✅ 已完成 · 多项式 `J1/J2/J3` 最小二乘，带 `cond(A)` 与残差诊断
- ✅ 已完成 · `encut` / `kpoints` / `sigma` 收敛闸门
- ✅ 已完成 · VASP 交接包：INCAR / KPOINTS / POSCAR / POTCAR.spec + 运行索引
- ✅ 已完成 · `doctor` 诊断：缺依赖、缺赝势、k 网格合理性
- ✅ 已完成 · **P2** C 心矩形：CrOCl 类磁序集、非共线 MAGMOM
- ✅ 已完成 · **P2′** 倾斜态：`solve_J_crect_cant.py` + 手册
- ✅ 已完成 · **P3** 非共线与 SOC：`LNONCOLLINEAR` / `LSORBIT` / `SAXIS`、倾斜 MAGMOM
- ✅ 已完成 · Binder 累积量 `U4` 与有限尺寸标度（FSS）+ `1/L` 参考
- ✅ 已完成 · **P4.1** kagome 系数层：`√3×√3`、`q = 0`、`3×3`、子晶格网
- ✅ 已完成 · **P5** 交变磁体（MnTe / CrSb）：几何判别 + 标量子晶格序参量
- ✅ 已完成 · 反解前的逐位点磁态核对与坍缩检测
- ✅ 已完成 · 蒙特卡洛采样诊断：接受率、`tau_int`、`n_eff`、binning SEM
- ✅ 已完成 · 机器可读能力表 `capabilities.json` + 闸门 `scripts/capabilities.py --check`
- ✅ 已完成 · 发布闸门：`check_sdist`、`make_release.py`、CI `package` 作业
- ⬜ 待机时 · **P4.2** kagome 物理层（Fe₃Sn₂ / Co₃Sn₂S₂ / Mn₃Sn）——需要真实非共线 DFT：提交 `√3×√3` 构型、回读 `OUTCAR`（`free energy TOTEN`）、写 `J_result.json`、蒙特卡洛求 `T_N`。手册已就绪：[`P4_2_KAGOME_RUNBOOK.md`](P4_2_KAGOME_RUNBOOK.md)。

## 该读哪份

| 读 | 为了 |
|---|---|
| [`README.md`](README.md) / [`README.zh.md`](README.zh.md) | 本说明书（English / 中文） |
| [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) | 症状 → 修法表，27 个错误键；`doctor` 的输出指向这里 |
| [`WORKFLOW.md`](WORKFLOW.md) | 逐阶段细节：文件、选项、口径 |
| [`dsh/2d-magnet/SKILL.md`](dsh/2d-magnet/SKILL.md) | agent skill：约定、闸门、怎么驱动流水线 |
| [`capabilities.json`](capabilities.json) | 机器可读能力边界——晶格、物理量、支持到哪一层 |
| [`examples/`](examples/) | 两份不用 VASP 就能看的完整产物树 |
| [`P4_2_KAGOME_RUNBOOK.md`](P4_2_KAGOME_RUNBOOK.md) | kagome 物理层手册，等机时就能跑 |
| [`CHANGELOG.md`](CHANGELOG.md) · [`CITATION.cff`](CITATION.cff) | 版本历史 · 引用元数据 |

在 Obsidian 里读：把仓库目录当成 vault 打开，看 `README.zh.obsidian.md`（同一份正文，加了内链目录与管线图）；英文版是 `README.obsidian.md`。两份都是生成物——改完 README 重跑 `python scripts/readme_to_obsidian.py`，`--check` 报告是否漂移。

## 许可

MIT——见 [`LICENSE`](LICENSE)。

## 引用

如果这个流水线对你的发表有帮助，按 [`CITATION.cff`](CITATION.cff) 引用：

```text
Wu, Sihong. 2d-magnet-pipeline: from a 2D magnet structure to exchange constants
and Monte Carlo T_C / T_N. Version 0.1.0, 2026. MIT license.
```
