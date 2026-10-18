# 变更日志（Changelog）

本文件只记录**对外可见的行为变化**（新增能力、修掉的静默错、口径/命名变更、发布工程）。
格式参考 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，版本号遵循
[语义化版本](https://semver.org/lang/zh-CN/)（`0.x` 阶段：minor 版本之间不保证向后兼容）。

- **能力的权威边界**：[`capabilities.json`](capabilities.json)（机器可读）+ 本 README 的「已知局限」段（人读）
- **症状 → 修复**：[`KNOWN_ISSUES.md`](KNOWN_ISSUES.md)
- **逐脚本细节**：[`WORKFLOW.md`](WORKFLOW.md)、[`dsh/2d-magnet/SKILL.md`](dsh/2d-magnet/SKILL.md)

## [Unreleased]

（暂无）

## [0.1.0] - 2026-10-18

首个对外发布版本。一句话：**给一个二维磁性单层，自动完成 结构 → 磁交换耦合常数 J → 经典 Heisenberg
Monte Carlo 的 T_C / T_N**，八种磁性格点共用一套规格层与统一产物。

### 新增

- **八种磁性格点**：蜂窝（J1/J2/J3）、三角（120° 非共线）、分子簇、双亚晶格铁磁网络（Fe₃GeTe₂ 类，
  走 DFT 反解）、矩形（CrSBr 类）、C 心矩形（CrOCl 类）、方形（La₂CuO₄ / K₂NiF₄ 类）、六角 kagome
  （Fe₃Sn 类）。每类都有 builder + solver + 契约测试。
- **两条 J 来源**：文献/给定值 → 正向 Monte Carlo；DFT 多磁态能量映射 → 反解 J（后者需要 VASP）。
- **统一 driver**：`run_2d_magnet.py` 的 `fetch / build / dft / ingest / kconv / mc / solve / render /
  report / doctor` 子命令；`--quiet` 静默进度、`--log <file>` 落盘每个阶段的完整输出。
- **规格层与可追溯产物**：`scripts/magnetic_orders.py`、`magnetic_orders.json`、`J_result.json`
  （带 `source` / `confidence` / `A` / `shells` / `residual` / `cond` 等 provenance）。
- **回归闸门**：`tests/smoke_test.py`（§1–§29）、8 个契约测试、`scripts/magnetic_orders.py --regress`、
  CI（`.github/workflows/ci.yml`，windows-latest + Python 3.12）。
- **机器可读能力边界**：`capabilities.json` + `scripts/capabilities.py --check/--show`
  （7 类晶格 / 7 个物理量 / 7 条「根本不支持」/ 5 道 fail-close 闸门；`physical_verified` 一律 `false`）。
- **可直接看的复现样例**：[`examples/`](examples/) 两组产物快照（方形 La₂CuO₄ 全程、
  kagome Fe₃Sn 合成空跑），不需要 VASP，也不需要额外依赖。
- **文档**：README、WORKFLOW、SKILL、KNOWN_ISSUES、NEW_TYPES_PLAN、P4.2 runbook、DMI memo，
  以及 Obsidian 阅读版备份 [`README.obsidian.md`](README.obsidian.md)
  （由 `scripts/readme_to_obsidian.py` 生成，`--check` 防漂移，smoke §29 守着）。
- **发布工程**：[`MANIFEST.in`](MANIFEST.in) + [`tests/check_sdist.py`](tests/check_sdist.py)
  （sdist 必须项闸门，smoke §28 用合成 tar.gz 反向验证闸门本身）+ CI `package` 作业
  （`python -m build` → 只装 wheel 跑一次 smoke → 校验 sdist），版本单一来源在 `pyproject.toml`。
- **引用信息**：[`CITATION.cff`](CITATION.cff)；发布记录：[`RELEASE_NOTES_v0.1.0.md`](RELEASE_NOTES_v0.1.0.md)。

### 修复（发布前复审，P0 为「会静默给错数」类）

- **P0-1 比热能量漏算单离子各向异性**：`scripts/mc_2d.py` 的 `_energy` 只算交换项，而 Metropolis 的
  dE 含 `--MAE`，导致「抽样包含各向异性、测量不包含」；现已补 `E_on / N`，并新增
  `--selftest`（与纯 Python 暴力求和逐位点对拍，含「含 A / A=0 一致 / A 确有贡献」三项）。
- **P0-2 反解前不核对逐位点磁矩**：新增 `scripts/verify_states.py`，`solve` 阶段当硬闸门（fail-close）
  自动跑；检出「两个非 fm 磁态坍缩成同一态」（行线性相关 → 残差 0 → 置信度反而虚高），
  证据写进 `states_verification.json` 与 `J_result.json`。确要看数只能 `--allow-unverified-states`，
  此时 provenance 明确保留 `verdict=FAIL` + `note`。
- **P0-3 `--quiet` 吞掉结论级警告**：`run()` 现在把 `[warn] / 未分辨 / 病态 / 坍缩 / 过渡带` 行转发到
  stderr（`--quiet` 下同样可见），失败时完整回显，并新增 `--log` 落盘每个阶段的完整 stdout/stderr。
- **P1-1 有限尺寸标度的 `±` 被当成统计标准误**：改称「有限尺寸离散度，不是标准误」；交点选择从
  「取最陡」改成「取最靠近全局参照」（避开噪声陡坡带来的选择偏差），并加 `1/L` 线性外推参考值。
- **P1-2 `--order auto` 误标序参量**：原判据只看 J1 符号，`J1>0 & J2<0` 的条纹反铁磁会被标成
  「交错（反铁磁）」；现改为按 T=0 候选基态能量（fm / neel / stripe）判定，条纹基态自动退回比热峰并
  警告，`mc_fss` 与 `mc_2d` 共用同一实现（遇条纹基态直接报错，不再给出错误的 T_N）。
- **P1-4 MC 无采样诊断**：每个温度输出 Metropolis 接受率、`τ_int`（Sokal 截断）、
  `n_eff = n_meas/(1+2τ)`、naive 与 binning 两种 SEM（`05_mc/mc_diag.json` +
  `[mc] 采样诊断: …` 一行）；三条阈值（`n_eff < 100` / 接受率 `< 0.05` / naive 低估 `> 3x`）转成警告。
- **P1-6 agent 协议会代替人 `git commit`**：`SKILL.md` 改成「回归 → 把改动清单交给用户提交」，
  并补上真实存在的降级开关表（`--allow-unverified-states` / `--ingest-force` / `--ingest-lenient` / `--no-s2`）
  与「不许用来绕过闸门」的约束。
- **`mc_2d.py` 自作用回绕**：`--size` 小于某个 shell 的 offset 时直接报错（原来会静默生成自环）。
- **driver 把比热峰报成「χ 峰」**：`summarize_mc` 改为按 CSV 表头判标签（`--order auto` 内部回退到 cv 时
  driver 并不知道，靠 `mc_order` 字符串判断会出错）。

### 变更

- **「已知局限」从 21 条收敛到 18 条**：删掉三条已失效的陈述（「T_C 取比热峰」的旧口径、
  「运行目录由 CWD 决定」、「GBK 控制台编码」），并把两条与正文重复的机制叙述压成纯边界；
  全仓「局限 N」交叉引用（README / KNOWN_ISSUES / NEW_TYPES_PLAN / WORKFLOW / P4.2 runbook /
  `scripts/altermagnet_check.py`）同步改号。
- README 新增「该读哪份（文档索引）」「发布状态」「打包闸门（sdist/wheel）」三段；Obsidian 用法改成
  可按步点击的三步（只支持「以文件夹为 vault 打开」，不支持孤立 `.md`）。

### 已知边界（本版本不打算解决，完整清单见 README「已知局限 18 条」）

- 仓库里**没有任何真实 DFT 计算结果**：J / T_C / T_N / D / 能带都只是**系数层或 MC 估算**，
  `capabilities.json` 的 `physical_verified` 一律 `false`。
- DMI / MTI 只到**识别层**（不给 D 值）；交变磁体只到几何判别 + 标量序参量；SOC 非共线链只到
  INCAR / MAGMOM 内容与快照回归。
- P4.2 的 k 点收敛只有手册级 `kconv` 流程，没有「自动进入主链」。
- MC 误差棒只能用 binning SEM（并看 `n_eff`）；只有单链，没有多副本 R̂ / 并行退火。
