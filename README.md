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
| Monte Carlo | `--mc-size` `--mc-tmin` `--mc-tmax` `--mc-nt` `--mc-refine-nt` `--mc-refine-max` `--no-mc-refine` |
| 开关 | `--skip-dft` `--skip-mc` `--quiet` `--outdir` |

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
- [ ] 子命令化（`fetch/build/dft/mc/solve/render/report`）+ `doctor` 环境自检
- [ ] `--model` 显式化 + 结构 provider 链（COD/OQMD/AFLOW/JARVIS）+ 元素参数表
- [ ] Binder 累积量 + 有限尺寸标度，给 T_C 误差棒
- [ ] `scripts/` 包化，支持非可编辑 `pip install .`

## License

[MIT](LICENSE)
