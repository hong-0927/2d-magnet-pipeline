# 边界说明（LIMITATIONS）

这份文件回答一个问题：**这个流水线给出的每个数，边界在哪**。

机器可读的版本在 [`capabilities.json`](capabilities.json)（`python scripts/capabilities.py --show` 打印可读版）；逐条症状与修法在 [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md)。下面是人读版的完整清单，共 18 条。

1. 仓库不带真实 DFT 数据。开箱跑出的数字都来自文献或合成能量，并且如实标注——demo 给的 `T_C`/`T_N` 不是预言。
2. 多磁态映射假设共线态映射到 Heisenberg 模型；非共线与 SOC 情形停在系数层（P3），没有经过验证的能量。
3. `S`、`moment`、`U` 是输入。内置表之外的元素必须显式给值：猜一个会静默给出错的物理，所以 driver 宁可拒绝。
4. 各向异性只通过 `--MAE`（默认 0）与 DMI 矩阵进入，两者都是输入，不是反解出来的量。
5. 蒙特卡洛是经典、单自旋更新的 Metropolis，跑在有限团簇上。`T_C`/`T_N` 是该团簇上的估计，不是热力学极限值。
6. 误差棒用 binning 标准误；`05_mc/mc_diag.json` 给出 `n_eff`、`tau_int` 与接受率，采样不足时能看出来。峰落在温度窗口边界时，`T` 会标成未分辨，而不是照常给一个数。
7. 有限尺寸标度用全局参照交点 + `1/L` 外推，外推值只是参考；没进标度区之前它不可信。
8. `--order auto` 用 T = 0 候选基态（fm / Néel / 条纹）决定用哪个序参量。条纹基态没有 ±1 交错序参量，流水线会自动退回比热峰并给出说明；近简并的模型请显式给 `--order`。
9. kagome 在系数层（P4.1），物理层（P4.2）等真实非共线 DFT——见 [`P4_2_KAGOME_RUNBOOK.md`](P4_2_KAGOME_RUNBOOK.md)。
10. 交变磁体只做几何判别（P5），标量序参量能到子晶格掩码，但不是完整的非共线描述。
11. 倾斜态（P2′）到系数层：`scripts/solve_J_crect_cant.py` 把 θ 当输入，不预测 θ。
12. DMI 以你提供的 `D` 矩阵进入，流水线不从能量里提取 `D`——见 [`DMI_MTI_MEMO.md`](DMI_MTI_MEMO.md)。
13. ingest 校验的是格式、收敛、结构与口径，不是物理。`ingest_report.md` 干净只说明 OUTCAR 可用，不说明计算是对的。
14. 不计算 k 分辨能带劈裂与 `Z₂`；`kconv` 只收敛 K 点网格。
15. `literature_J.json` 里的条目是起点；条目为 `null` 时 driver 拒绝运行，不会编一个值。
16. 磁态核对覆盖共线态；非共线或 SOC 运行会以 `SKIP` 退出（不是通过）。
17. 环境假设：Python ≥ 3.9 只做到语法层检查，实际回归跑在 3.12；`matplotlib` 与 `requests` 可选，缺失时退到有文档的兜底路径。
18. 发布件由 `scripts/make_release.py` 从工作树打包、`tests/check_sdist.py` 校验；PyPI 式的 `sdist`/`wheel` 构建由 CI 的 `package` 作业执行。
