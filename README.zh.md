# 2d-magnet-pipeline

把单层结构变成交换常数 `J1/J2/J3`，再用蒙特卡洛给出 `T_C` / `T_N` 的估计——每个数字都过一道 fail-closed 闸门。

输入是一个 `POSCAR`（Materials Project、COD / OQMD / AFLOW / JARVIS 公共链，或你自己的文件）。流水线识别磁性格子、建超胞与磁序集、装配 VASP 输入、回读 `OUTCAR`、反解交换常数，然后跑经典 Heisenberg 蒙特卡洛。到 J 有两条路：

- **文献 J 或用户给定的 J → 正向蒙特卡洛**：不需要 VASP；
- **DFT 多磁态能量映射反解 J**：需要 VASP 能量。

## 能力一览

| 能力 | 带来什么 |
|---|---|
| **八类磁性格子，一个规格层** | honeycomb、三角（120° 非共线）、分子簇、双亚晶格铁磁网络（Fe₃GeTe₂ 类）、矩形（CrSBr 类）、C 心矩形（CrOCl 类）、方形（La₂CuO₄ / K₂NiF₄ 类）、六角 kagome（Fe₃Sn 类）——都走同一层规格、都写进同一份 `J_result.json` |
| **J 带着出处，不只是一个数** | `J_result.json` 记录 `source`（`literature` / `user` / `dft` / `demo`）、`cond(A)`、残差、`E0`，以及磁态是否逐位点核对过 |
| **反解之前先核对磁态** | `scripts/verify_states.py` 把每个离子的磁矩与目标比对，并检测两个磁序塌成同一个态；不过就拒绝写 J |
| **ingest 是硬闸门** | 每个 `OUTCAR` 过 11 道检查：格式、收敛、结构、MAGMOM 一致性、能量口径 |
| **估计值自带误差** | Binder 累积量 `U4`、带 `1/L` 参考的有限尺寸标度，以及每个温度的 `tau_int`、`n_eff`、接受率与 binning SEM（`mc_diag.json`）；峰落在窗口边界时标成未分辨，而不是照常给数 |
| **非共线与 SOC** | `LNONCOLLINEAR` / `LSORBIT` / `SAXIS`、倾斜 MAGMOM 磁序集，以及 CrOCl 家族的 `scripts/solve_J_crect_cant.py` |
| **不猜参数** | `--S/--moment/--U` 只对 Mn/Cr/Ni/Co/V/Gd 有内置表；Fe 被刻意排除（有效自旋本身有歧义），driver 会停下来问；`literature_J.json` 里 `null` 的条目直接报错，不静默跑 |
| **能力边界可查询** | [`capabilities.json`](capabilities.json) 写明每类晶格与每个物理量支持到哪一层；`scripts/capabilities.py --check` 保证它与代码同步，`--show` 打印可读版 |
| **同时是 agent skill** | `dsh/2d-magnet/SKILL.md` 把约定、闸门与驱动方式写全，Claude Code / Codex 或任何读指令文件的 agent 都能用 |
| **可离线复现** | `examples/` 有两份完整产物树，不用 VASP 就能检查并重跑；发布包由 `scripts/make_release.py` 确定性打包 |
| **交变磁体判别** | `scripts/altermagnet_check.py` 判断两套亚晶格是否由保自旋的操作联系（交变磁体几何，标量序参量） |

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

`--quiet` 抑制进度输出，但绝不吞掉结论级警告——未分辨的 `T_C`、病态的 `cond(A)`、坍缩的磁态。要看每个阶段的完整输出用 `--log <文件>`。

## 安装

**环境要求**

| 项 | 说明 |
|---|---|
| Python | ≥ 3.9（发布回归跑在 3.12） |
| 依赖包 | `numpy` ≥ 1.24、`numba` ≥ 0.59、`matplotlib` ≥ 3.7、`requests` ≥ 2.31，见 `requirements.txt` |
| 可选包 | `Pillow` ≥ 9.0 —— `matplotlib` 不可用时的出图兜底 |
| VASP | 可选：`dft` / `ingest` / `solve` 需要你自己跑 VASP；文献 J → Monte Carlo 不需要 |
| 网络 | 可选：只有 `fetch` 下载结构需要；`examples/` 与 `tests/fixtures/` 完全离线 |
| 磁盘 | 仓库约 15 MB；每次运行在 `--outdir` 下写几 MB |

**1. 取代码**

```powershell
git clone <repo-url> 2d_magnet_pipeline
cd 2d_magnet_pipeline
```

没有 git：直接下载 Release 里的 `2d-magnet-pipeline-0.1.0.tar.gz` / `.zip` 解压，内容与仓库一致。

第 2–4 步有脚本版。在仓库根目录：

```powershell
powershell -ExecutionPolicy Bypass -File install.ps1
```

```bash
bash install.sh
```

两者都会建 `.venv`、装 `requirements.txt`、跑 `doctor`；换位置用 `-Venv <路径>`（bash 版是 `--venv <路径>`）。Windows 上按上面那样带 `-ExecutionPolicy Bypass` 跑 PowerShell 版——否则激活脚本会被执行策略拦住。

**2. 建并激活虚拟环境**

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
```

```bash
python3 -m venv .venv
. .venv/bin/activate
```

PowerShell 报「禁止运行脚本」时：改用 `cmd` 里的 `.venv\Scripts\activate.bat`，或只对当前会话放开 `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass`。Windows 上建议用 `py` 启动器（它绕开 Microsoft Store 的 `python.exe` 占位程序）；没有 `py` 就用 `python.exe` 的完整路径。

**3. 装依赖**

```powershell
py -m pip install --upgrade pip
py -m pip install -r requirements.txt
```

`matplotlib` 体积最大，并会带入 `pyparsing`、`cycler`、`contourpy`、`fontTools`、`kiwisolver`。装不上时（离线机、精简站点）改装 `Pillow>=9.0`：六条出图路径会自动走纯 Pillow 实现，图照样出（打包 extra：`pip install -e ".[render-fallback]"`）。`requests` 只被 `fetch` 用；没有它，除结构下载外一切照常。

**4. 验证安装**

```powershell
py run_2d_magnet.py doctor
```

`doctor` 每项一行，随后 `---`，再按每个 WARN 给出错误目录里的对应条目。退出码 `0` = 无 WARN；`1` = 至少一项环境问题，不代表代码坏。正常安装会对 Python、`numpy`、`numba`、`scripts/`、`literature_J.json`、`tests/fixtures/` 打 `[OK]`；`MP_API_KEY` 与 `CWD` 两行只是信息。

**5. 可选：装成命令**

```powershell
py -m pip install -e .
2d-magnet --version      # 2d-magnet 0.1.0
2d-magnet doctor
```

**6. 可选：当 agent skill 用**

`dsh/2d-magnet/` 是自包含的 skill（约定 + 闸门 + 驱动方式）。把整个目录拷进你的 harness 技能目录——Claude Code 是 `~/.claude/skills/`，或项目内 `.claude/skills/`；`AGENTS.md` 风格的 harness 直接让 agent 读 `dsh/2d-magnet/SKILL.md`。skill 驱动的是同一个 `run_2d_magnet.py`，不需要另装东西。

**7. 可选：VASP 那一侧**

`dft` 会在 `03_dft/` 与 `runs/<order>/` 写出 `INCAR`、`KPOINTS`、`MAGMOM` 和 `POTCAR.spec`（每元素推荐的 POTCAR 标签）。POTCAR 有许可、**不随仓库分发**：提交前自己从 VASP 发行包里按元素拼到 `runs/<order>/POTCAR`。`ingest` 只读 `OUTCAR`，所以计算本身可以在任意机器上跑。

**常见问题**

| 症状 | 处理 |
|---|---|
| `python` 打开 Microsoft Store / `command not found` | 改用 `py`（Windows 启动器）或 `python.exe` 完整路径 |
| 出图时 `ModuleNotFoundError: No module named 'pyparsing'` | `matplotlib` 装得不完整——装 `Pillow`，自动走兜底渲染 |
| 中文输出乱码、`UnicodeEncodeError` | 运行前设 `$env:PYTHONIOENCODING="utf-8"` |
| `doctor` 报 `CWD 不可写` | 换到可写目录运行，或始终显式给 `--outdir` |
| `fetch` 返回 401/403 或空结构 | 设 `MP_API_KEY` 或传 `--mp-key`，或用 `--poscar` 给自己的文件 |
| `[params] ... 不在内置元素参数表` | 显式给 `--S --moment --U`，或用 `--model` / `--literature` 条目 |

上面每个 WARN 与错误键的完整「症状 → 修法」都在 [`KNOWN_ISSUES.md`](knowledge/KNOWN_ISSUES.md)。

代码也可当库用（`import run_2d_magnet`）；`pip install .` 之后还有 `2d-magnet` 命令行入口。

## 快速开始

`examples/` 里有两份不用 VASP 就能复现的完整运行（文献 J 进，`T_C` 出）：

```powershell
py run_2d_magnet.py --model la2cuo4 --poscar tests/fixtures/La2CuO4_POSCAR `
   --outdir C:\work\la2 --skip-dft --mc-size 12 12 --mc-seeds 2
py scripts\solve_J_square.py C:\work\la2\02_supercell\cases --demo
```

`examples/la2cuo4_square/` 是这样跑出来的目录树；`examples/kagome_fe3sn_dryrun/` 是一条走完整 DFT 形状的 kagome 空跑（合成能量）。两者各有 `README.md`，列出产出了什么。

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
```

`solve` 在反解之前自动跑 `scripts/verify_states.py`：两个磁序塌成同一个态、或某个位点的磁矩与目标矛盾，该阶段就以非零退出，不写 `J_result.json`。

## 该读哪份

| 读 | 为了 |
|---|---|
| [`README.md`](README.md) / [`README.zh.md`](README.zh.md) | 本说明书（English / 中文） |
| [`AGENTS.md`](AGENTS.md) | 给 agent 的五条硬规则（先读这个） |
| [`STRUCTURE.md`](STRUCTURE.md) | 目录约定——新内容该放哪 |
| [`KNOWN_ISSUES.md`](knowledge/KNOWN_ISSUES.md) | 症状 → 修法表，27 个错误键；`doctor` 的输出指向这里 |
| [`WORKFLOW.md`](knowledge/WORKFLOW.md) | 逐阶段细节：文件、选项、口径 |
| [`dsh/2d-magnet/SKILL.md`](dsh/2d-magnet/SKILL.md) | agent skill：约定、闸门、怎么驱动流水线 |
| [`capabilities.json`](capabilities.json) | 机器可读能力边界——晶格、物理量、支持到哪一层 |
| [`DMI_MTI_MEMO.md`](knowledge/DMI_MTI_MEMO.md) | DMI / MTI：怎么判别、能说什么 |
| [`LIMITATIONS.md`](knowledge/LIMITATIONS.md) | 边界说明——每个数能代表什么、不能代表什么（18 条） |
| [`examples/`](examples/) | 两份不用 VASP 就能看的完整产物树 |
| [`P4_2_KAGOME_RUNBOOK.md`](knowledge/P4_2_KAGOME_RUNBOOK.md) | kagome 物理层手册，等机时就能跑 |
| [`CHANGELOG.md`](CHANGELOG.md) · [`CITATION.cff`](CITATION.cff) | 版本历史 · 引用元数据 |


## 许可

MIT——见 [`LICENSE`](LICENSE)。

## 引用

如果这个流水线对你的发表有帮助，按 [`CITATION.cff`](CITATION.cff) 引用：

```text
Wu, Sihong. 2d-magnet-pipeline: from a 2D magnet structure to exchange constants
and Monte Carlo T_C / T_N. Version 0.1.0, 2026. MIT license.
```
