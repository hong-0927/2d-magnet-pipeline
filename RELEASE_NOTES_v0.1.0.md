# 2d-magnet-pipeline v0.1.0

首个对外发布版本（2026-10-18）：**给一个二维磁性单层，自动完成 结构 → 磁交换耦合常数 J → 经典
Heisenberg Monte Carlo 的 T_C / T_N**。八种磁性格点共用一套规格层与统一产物。

## 能做什么

- **八种磁性格点**：蜂窝（J1/J2/J3）、三角（120° 非共线）、分子簇、双亚晶格铁磁网络（Fe₃GeTe₂ 类，
  DFT 反解）、矩形（CrSBr 类）、C 心矩形（CrOCl 类）、方形（La₂CuO₄ / K₂NiF₄ 类）、六角 kagome（Fe₃Sn 类）。
- **两条 J 来源**：文献/给定值 → 正向 MC；DFT 多磁态能量映射 → 反解 J（需 VASP）。
- **统一 driver**：`fetch / build / dft / ingest / kconv / mc / solve / render / report / doctor`。
- **可追溯产物**：`magnetic_orders.json`（规格层）+ `J_result.json`（含 source / confidence / A / shells /
  residual / cond 等 provenance）+ `states_verification.json`（反解前逐位点磁态核对）。
- **确定性回归**：`tests/smoke_test.py`（§1–§29）+ 8 个契约测试 + `magnetic_orders.py --regress` + CI（windows-latest / Python 3.12）。

## 明确不做什么

- **没有任何真实 DFT 计算结果**：J / T_C / T_N / D / 能带都只是**系数层或 MC 估算**；
  [`capabilities.json`](capabilities.json) 里 `physical_verified` 一律 `false`，谁改成 `true` 会被 `--check` 与 smoke 拒绝。
- DMI / MTI 只到**识别层**（不给 D 值）；交变磁体只到几何判别 + 标量序参量；SOC 非共线链只到 INCAR / MAGMOM 内容与快照回归。
- P4.2 的 k 点收敛只有手册级 `kconv` 流程，没有自动进入主链。
- MC 误差棒只能用 **binning SEM**（并看 `n_eff`）；只有单链，没有多副本 R̂ / 并行退火。

完整清单见 [README.md](README.md) 的「已知局限 18 条」；逐条症状与修法见 [KNOWN_ISSUES.md](KNOWN_ISSUES.md)。

## 安装与冒烟

```powershell
pip install -r requirements.txt
python run_2d_magnet.py doctor

# 不需要 VASP 的一次完整跑（文献 J → MC），产物在 %TEMP%\la2_example
python run_2d_magnet.py --model la2cuo4 --poscar tests/fixtures/La2CuO4_POSCAR `
  --outdir $env:TEMP\la2_example --skip-dft --mc-size 12 12 --mc-nt 8 --mc-tmax 300 --mc-seeds 2
```

不想先跑、只想看产物长什么样：直接看 [`examples/`](examples/)（两组快照，无需 VASP、无需装依赖）。

## 校验清单（本版本发布前全绿）

```powershell
python tests/smoke_test.py                     # §1–§29 全部通过
python scripts/magnetic_orders.py --regress    # 规格层回归
python scripts/capabilities.py --check         # 能力表与代码/版本/闸门一致
python scripts/readme_to_obsidian.py --check   # Obsidian 备份不漂移
python scripts/mc_2d.py --selftest             # 比热能量含单离子项 (7/7)
python scripts/mc_fss.py --selftest            # FSS 交点选择与 1/L 外推
python scripts/verify_states.py --selftest     # 反解前的磁态闸门
python tests/check_sdist.py "dist\*.tar.gz"    # sdist 必须项闸门
```

## 附件

- `2d-magnet-pipeline-0.1.0.tar.gz`（源码包，136 个文件）/ `.zip`
- `SHA256SUMS`（两个包的 sha256）
- 许可：MIT License © 2026 吴思泓（见 [LICENSE](LICENSE)）

> 注意：本版本**不是** PyPI 上的分发包。源码包由仓库自带的 [`scripts/make_release.py`](scripts/make_release.py)
> 打出（名单取 `git ls-files --cached --others --exclude-standard`，固定前缀/uid/mtime，同一棵树重复构建字节一致）；
> 产物在 `dist/`，不进仓库。PyPI 式 sdist/wheel（`python -m build`）由 CI 的 `package` 作业产出，本机没有
> `setuptools` 且离线，所以没有在本地跑过 —— 第一次 push 后请到 Actions 确认该作业。
