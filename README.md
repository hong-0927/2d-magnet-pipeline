# 2d-magnet-pipeline

给一个二维磁性单层，自动完成：**结构 → 磁交换耦合常数 J → 经典 Heisenberg Monte Carlo 临界温度 T_C / T_N**。

- 支持磁性格点：**蜂窝**（J1/J2/J3）、**三角**（120° 非共线）、**分子簇**、**双亚晶格铁磁网络**（Fe₃GeTe₂ 类，走 DFT 反解）。
- 两种 J 来源：**文献/给定值 → 正向 Monte Carlo**；或 **DFT 多磁态能量映射 → 反解 J**（后者需要 VASP）。
- 附 DeepSeek Harness 的 agent skill（`dsh/2d-magnet/SKILL.md`）。

---

## 快速开始

```powershell
git clone <this repo>
cd 2d_magnet_pipeline
pip install -r requirements.txt     # numpy / numba / matplotlib / requests

# 1) 文献 J 的正向计算（MnPS3，需 Materials Project key 下载结构）
python run_2d_magnet.py --formula MnPS3 --literature MnPS3 --mp-key <MP_API_KEY>

# 2) 已有 POSCAR + 显式 J（跳过下载；J 需自备）
python run_2d_magnet.py --poscar tests/fixtures/MnPS3_POSCAR `
  --element Mn --nonmag P S --S 2.5 --moment 5.0 --U 4.0 --J1 -0.3 --J2 0.1 --J3 0.02

# 3) 双亚晶格反解（FGT 类，需 VASP；详见 dsh/2d-magnet/SKILL.md 或 hamiltonian.md 用法）
python scripts/build_fgt_cases.py <monolayer_POSCAR> --outdir Fe3GeTe2_pipeline
```

> `--formula` 下载结构需要 Materials Project API key（环境变量 `MP_API_KEY` 或 `--mp-key`）。
> 已有 POSCAR 用 `--poscar` 跳过下载。`--J1` 必填（J>0 铁磁 / J<0 反铁磁），除非用 `--literature`。

## 安装

- **Python 3.9+**（`numba` 兼容性），依赖见 [requirements.txt](requirements.txt)。
- 可选：`pip install -e .`（**可编辑**安装）得到 `2d-magnet` 命令。`scripts/` 由入口相对定位，因此只支持可编辑安装；否则直接用 `python run_2d_magnet.py`。
- Windows 封装：`dsh/bin/2d-magnet.cmd`（自动带 `--quiet`；解释器查找顺序 `%MAGNET_PYTHON%` → `python` → `py`）。

## CLI 参考（节选，完整见 `python run_2d_magnet.py --help`）

| 组 | 参数 |
|---|---|
| 结构来源 | `--formula` `--mp-id` `--poscar` `--mp-key` |
| 材料 | `--element` `--nonmag P S` `--S` `--moment` `--U` |
| J 值 | `--literature <key>` `--J1` `--J2` `--J3` `--MAE` |
| 磁性格点 | `--lattice {auto,honeycomb,triangular,cluster}` `--cluster-cutoff` |
| Monte Carlo | `--mc-size` `--mc-tmin` `--mc-tmax` `--mc-nt` `--mc-refine-nt` `--mc-refine-max` `--no-mc-refine` `--mc-seeds` `--fss` `--fss-sizes` |
| 开关 | `--skip-dft` `--skip-mc` `--quiet` `--outdir` `--doctor` |

### 元素参数表

`--S` / `--moment` / `--U` 未显式给出时，按磁性元素从内置表兜底（**Mn / Cr / Ni / Co / V**）。Fe 因 Fe²⁺/Fe³⁺/巡游有效自旋三种取值差异太大（S=2/2.5/1）而**故意不列**，须显式给 `--S/--moment`。优先级：显式 CLI > `--literature` 表 > 元素表；落表时会打 stderr 提醒。

### Monte Carlo 输出与误差

- 序参量路径（蜂窝铁磁/Néel 反铁磁）的 `05_mc/mc_mt.csv` 多一列 **Binder 累积量 U₄**（`U₄ = 1 − ⟨M⁴⟩/((1+2/3)⟨M²⟩²)`）：有序 → 0.4、高温无序 → 0。跑多个 `--mc-size`，U₄(T) 的交点即有限尺寸标度定出的 T_C（比单一 χ 峰更稳健）。
- `--mc-seeds N`（默认 1）跑 N 个独立种子并平均，峰位报 `± 种子间标准差`（误差棒）。

### 有限尺寸标度（`--fss`）

`--fss` 会在 MC 阶段之后自动跑 `scripts/mc_fss.py`：对 `--fss-sizes`（默认 `12 18 24 30`）逐尺寸算 U₄(T)，取

1. **两两交点**：相邻尺寸的 U₄(T) 在过渡带（U₄∈[0.02,0.38]）内取斜率最陡的符号变化点插值；
2. **交叉验证**：使各尺寸 U₄ 之间方差最小的温度（collapse）。

两者均值即临界温度 `T_C/T_N`，报告 `± 交点间标准差`；结果写 `05_mc/fss_fss.{csv,json,png}`。独立用法：`python scripts/mc_fss.py <model.json> --sizes 12 18 24`。提示：二维无真正相变（Mermin-Wagner），交点 T_C 本质是有限尺寸伪相变估计，且随种子噪声波动，须 `--mc-seeds > 1` 收敛。

## 目录结构

```
2d_magnet_pipeline/
  run_2d_magnet.py        # 主入口（单文件驱动整个流水线）
  scripts/                # 各阶段 + 分模型解算器（solve_J3 / solve_cluster_J / solve_fgt_J / ...）
  dsh/2d-magnet/SKILL.md  # agent skill
  dsh/bin/2d-magnet.cmd   # Windows 封装
  literature_J.json       # 已核实文献 J 表（含符号约定说明）
  tests/                  # smoke_test.py + test_fgt_e2e.py（fixtures/ 内置，无外部路径依赖）
  WORKFLOW.md             # 阶段与脚本来源说明
```

## 测试

```powershell
python tests/smoke_test.py        # 单元/回归（分类、峰边界检测、工作流文档等）
python tests/test_fgt_e2e.py      # FGT 双亚晶格反解链路端到端（mock DFT，不跑 VASP）
```

测试 fixture 已内置在 `tests/fixtures/`，不依赖本机任何绝对路径；可用环境变量
`MAGNET_FIXTURES` / `MAGNET_FGT_CASES` / `MAGNET_TEST_TMP` 覆盖。

## 已知局限

1. 结构下载当前只走 Materials Project（无 key 时须自备 POSCAR）；部分材料 MP 没有（如 Fe₃GeTe₂），需手动建层。
2. **磁性格点分类是纯几何判据**（近邻 < `--cluster-cutoff` 就并簇）。连通的多亚晶格铁磁网络（FGT 类）与真分子簇（Nb₃I₈ 类）签名可能相同，会被误判——driver 会打 stderr 警告但**不改路由**，请人工确认。
3. `--J1` 为正向 MC 的必需输入；未核实文献条目（`literature_J.json` 里 `null`）会明确报错而非硬跑。
4. MC 的 T_C/T_N 取比热峰 `argmax(Cv)`；峰落在温度窗口边界时**不再把边界值当结果**，而是显式「未分辨」并自动缩放重扫（`--mc-refine-*`）。
5. 双亚晶格反解：共线、只到 J4、统一 S=1、结构需先弛豫、MAE 未纳入、k 网格 3×3×1 偏粗。
6. DFT 反解要求 meV 级总能差；金属体系建议 `ISMEAR=1/SIGMA=0.2` 且 5 个磁态设置完全一致。
7. 运行目录按 CWD 决定（`./<formula>_pipeline`）；建议在仓库根运行。

## 路线图

- [x] 硬编码路径清理 + 测试 fixture 内置 + 打包基线（LICENSE/README/requirements/pyproject）
- [x] `doctor` 环境自检（`--doctor`）
- [x] 元素参数表（Mn/Cr/Ni/Co/V 兜底 S/moment/U）+ `--lattice` 显式化
- [x] Binder U₄ 累积量 + `--mc-seeds` 误差棒
- [x] 自动多尺寸有限尺寸标度（U₄ 交点定位 T_C）
- [ ] 子命令化（`fetch/build/dft/mc/solve/render/report`）
- [ ] 结构 provider 链（COD/OQMD/AFLOW/JARVIS）
- [ ] `scripts/` 包化，支持非可编辑 `pip install .`

## License

[MIT](LICENSE)
