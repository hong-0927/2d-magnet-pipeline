# kagome 非共线 DFT 分支手册（装配 → 交接 → 提交 → 判定 → 回收 → J → T_C）

> **默认路线不需要非共线。** 默认五共线态（`fm/q0_ud/q0_duu/sq3_updntri/sq3_alt`，
> `rank=4, cond≈100`）已经满秩；本手册的**非共线路线**只在两种目的下才需要：
> ① 把 **J1 与 J1′ 分开**（3 壳模型把它们合并成一条，见 §2/§9）；② 做 120° 基态能量对照。
> 若目标只是「拿到 J1/J2/J3 和 T_C」，走默认共线路线即可，不必进本手册。

> **本文件是操作手册, 不是结果。** 本仓库**从未跑过真实 DFT**（`README.md:260`）, 所以下面
> **没有任何实测的 J、T_C 或机时**。它规定的是: 上真实机时之前必须钉死哪些量、怎么装配、
> SCF 算到什么程度算成功、OUTCAR 怎么回收、反解出问题先怀疑什么、以及**什么情况下不许下结论**。
>
> 适用对象: 有 VASP 机时、愿意按本手册逐步留痕的人。**仓库不代提交作业、不写 POTCAR、不提供机时**
> （`README.md:81`、`run_2d_magnet.py:91`）—— 但**交接包与提交脚本可以生成**（`dft_handoff.py`,
> 见 §5.0），**收回通道已自动化**（`ingest_dft.py`, 见 §5.1）。
>
> ## ⚠ 本轮（2026-10 第二次补做）改了三件事, 先读这三条
>
> 1. **锚点名字改了: Fe3Sn 报的是居里温度 T_C, 不是 T_N。** 已核实的来源明确称 Fe3Sn 为
>    *kagome ferromagnet*（`P6_3/mmc`）⇒ 它的磁有序温度是 **T_C**。任务与旧文档里写的
>    「T_N 与文献量级一致」对 Fe3Sn 是**命名错误**。同族不许互相冒充: **FeSn**（`P6/mmm`）才是
>    反铁磁（365 K 以下）, **Fe3Sn₂**（`R-3m`）是**易轴**铁磁（T_C = 612 K, 室温以下有自旋重取向）。
> 2. **文献锚点已核实并落盘**：`scripts/literature_J.json` 现在**有** Fe3Sn 条目（以前没有）。
>    写入的是 **有序温度 / 结构 / 各向异性** 锚点（T_C = 705 K、a = 5.4604 Å、c = 4.3458 Å、
>    2.27 μ_B/Fe、易面 = kagome 面、`ground_state = "fm"`）；**J1/J2/J3 仍然刻意保持 `null`**
>    —— 来源自承 FexSny 的交换常数尚未确立, 它给的是饱和矩与各向异性, **不是**壳分辨的 J。
>    ⇒ `solve_J_kagome.py --ref Fe3Sn` 会被按 null **明确拒绝**（设计行为）。**不要**为了
>    「让 --ref 能用」去填一个来源不明的 J。
> 3. **2D MC 的临界温度与体相 T_C 差多少, 已用本仓库引擎量出来了（不需要任何 DFT）**：
>    见 §7.1。结论是 **M4 的右半在当前模型类里不可通过** —— 不是"暂时没机时", 而是
>    缺层间耦合。**先读 §7.1 再决定要不要开机时。**

---

## 0. 一页速览

| # | 步骤 | 命令（仓库内） | 产物 | 验收判据 | 偏差/失败怎么办 |
|---|------|----------------|------|----------|------------------|
| 1 | 钉死物理量 | —— | §1 的表填满 | 化合物/a/S/μ/U 都有**来源** | 没来源就停, 不许"先用 5.4 Å 试试" |
| 2 | 建构型 | `build_kagome_cases.py … --selftest` | `cases/<order>/{POSCAR,magmom.txt}` + `magnetic_orders.json` | selftest PASS; 配位数 `(4,4,6)` | 配位数不对 = 壳合并, 别硬跑 |
| 3 | 装配输入 | 驱动 `dft` 步 | `03_dft/{INCAR.static,INCAR.relax,KPOINTS,POTCAR.spec}` + `runs/<order>/{INCAR,KPOINTS,POSCAR}` | 非共线 MAGMOM 三分量; `runs/` 下 4 个序齐 | 缺 POTCAR: 自己按 POSCAR 顺序拼 |
| 3.5 | **打交接包** | `dft_handoff.py pack --root <proj> --out <dir>` | 交接包（README/MANIFEST+sha256/CHECKLIST/submit.sh/jobscript.template/POTCAR.spec/runs/） | MANIFEST 记全 sha256; **不含 POTCAR** | 缺 INCAR/KPOINTS/POSCAR 的序被**拒绝打包**（§5.0） |
| 4 | **外部跑 VASP** | 集群上 `bash submit.sh`（自动识别 5 种调度器） | `runs/<order>/OUTCAR` | §4 的四条 SCF 判据**全部**通过 | 不收敛 → 按 §4 的协议, **不许**喂未收敛能量 |
| 5 | **回收** | `run_2d_magnet.py ingest --outdir <proj> [--runs <dir>]`（或 `scripts/ingest_dft.py <runs> <cases>`） | `cases/<order>/OUTCAR` + `cases/ingest_manifest.json` | 10 条硬闸门全过（§5.1） | **不通过就一件都不写**（原子性）; 先 `--dry-run` 看闸门; 每条闸门的「为什么/怎么修」会打印 |
| 6 | 反解 | `solve_J_kagome.py <cases>` / 驱动 `solve` | `J_result.json` | `source=vasp` 且 `confidence=high`（残差<1e-4、cond<1e3） | `confidence=low` 不许拿去做物理结论 |
| 7 | 磁学 | `poscar_to_mc3 --jfile` → `mc_2d.py` / `mc_fss.py` | T_C、U4 交点 | **§7.1**: 2D MC 临界温度**不能**直接对 3D T_C（比值随尺寸 0.74–0.89, 且无收敛值） | §7: 没有锚点 ⇒ **P4.2 不能验收**; §7.1: **有锚点也过不了当前判据**, 必须先加层间 J_c |
| 8 | 留痕 | —— | 实测机时写回 §8; 文献值写进 `literature_J.json` | 有时间戳与出处 | 无出处的文献值 = 不许填 |

---

## 1. 开工前必须钉死的量（每一项都影响可发表性）

| 量 | 代码里的现状 | 必须由人决定/记录 | 不决定会怎样 |
|----|--------------|-------------------|--------------|
| 化合物 | `README.md:295` 提了 Fe₃Sn₂ / Co₃Sn₂S₂ / Mn₃Sn, 但 **fixture 只有 `tests/fixtures/Fe3Sn_POSCAR`**, builder/preset/render 都只认 Fe3Sn | 选**一个**; 换材料要先加 fixture + `--model` 预设 | 拿 Fe3Sn 的几何去讲别的材料的物理 |
| 晶格常数 a | `--model fe3sn` 的 a=5.4 Å 是**理想化占位**（`README.md:268`） | 实验/弛豫来源 + 数值 | 所有 J 都挂在错的几何上 |
| 自旋 S | 预设 `--spin 1.0` | 由价态/磁矩定, 并写进 `J_result.json`（S² 换算用） | J 与 J~ 混用（x1–x4 量级错） |
| 磁矩 μ | 预设 `--moment 2.0` | 同上; 与 MAGMOM 一致 | 构型不是你以为的那个序 |
| U（LDAU） | `assemble_dft.py:143–158`, `--U` **只加在磁性元素**上 | U 值 + 出处（文献/收敛测试）; 是否要 U 本身也要论证 | 3d 化合物不给 U 会得到金属性的假基态 |
| POTCAR/PAW | 仓库**只写 `POTCAR.spec`**（`assemble_dft.py:289–290`）, 从不写 POTCAR | 每种元素的 POTCAR 版本 + 按 POSCAR 顺序拼接 | 顺序错 = 全盘错, 且不报错 |
| k 网格 | `KPOINTS` 默认模板 `3 3 1`（header 已修正, 不再写 "2x2x1 monolayer"）; 可用 `--kpoints 6 6 1` 覆盖 | 至少抽查一次 k 收敛（工具见 §5.3, 记录在 §8） | 采信未收敛的 k 网格 |
| ENCUT | 默认 520（`assemble_dft.py:225`） | 抽查一次（或引用同族 POTCAR 的 ENMAX 规则） | 同上 |
| 超胞/真空 | `--kagome-supercell 3 3 1`（默认）, n0/n1 **必须是 3 的倍数** | 超胞 + 真空层; 记录原子数（本几何 3×3×1 → **75 原子**） | 120° 序放不下 / 层间耦合污染 |
| 并行 | `NPAR` 默认 4（可用 `--npar 8` 覆盖） | 节点/核数 + NPAR/KPAR | 机时浪费或算不动 |

> 一句话: **本手册不替你选参数。** 上面任何一格空着, 就跑不出可解释的结果。

---

## 2. 构型集：为什么是这几个, 以及不许动什么

- builder 产出 **7 个序**; 反解只用 4 个 (`selected_orders`):
  `fm` / `q0_120` / `sq3_updntri` / `sq3_alt`。
  **rank = 4, cond = 98.51, |det| = 2.1e5**（`LATTICE_ROADMAP.md:101`）——设计矩阵状态良好。
- **`q=0` 家族（fm/120/stripy）自身 rank = 2, cond ≈ 2.3e34** ⇒ 只靠它**无法**定 4 个未知量。
  这是几何事实, 不是 bug。
- **退化禁令**: 丢掉一个 selected 序, 或换成一个与已有行线性相关的序, 会触发
  `rank < 4` 或 `cond > 1e4` 的硬闸门（`KNOWN_ISSUES.md:35`）。
  **不要**用"换一个看起来更独立的序"去骗过闸门 —— 闸门报错的意思就是"这批构型没有信息",
  换序只有在**确实**换了方向时才有意义（要重新算 rank/cond, 而不是重跑一遍）。
- **硬要求**: 7 个序全部是**三分量 MAGMOM** ⇒ 装配必须 `--noncollinear`（共线 + 三分量会硬报错,
  `KNOWN_ISSUES.md:30`）。
- **LSORBIT 对现有 7 个序不是必需**: 它们都不在 `mag_vectors.SOC_ORDERS`（cant/spiral/tilt）里,
  所以 `--noncollinear` 单独就够（`SKILL.md:317–323`）。加 `--lsorbit` 会让单点贵 **1.5–3×**;
  但**要分离 J1 与 J1′ 就必须加倾斜/SOC 构型**（`SKILL.md:331` 只写了"属 P4.2+", **仓库没有这些构型**）。
- ⚠ **模型只有 3 个壳**: J1(a/2) / J2(a√3/2) / J3(a), 其中 **J1 是 J1 与 J1′ 的平均**
  （`LATTICE_ROADMAP.md:98`）。所以 P4.2 交付的是**三壳平均** —— **不许**报"分离出了 J1/J1′"，
  也不许把本手册的结果说成"解出了 J1/J1′/Jd 三自由度"（那是 `LATTICE_ROADMAP.md:102` 里被过度
  声明的旧标题, 见 §9）。

---

## 3. 装配（仓库做什么、不做什么）

```bash
# 1) 建构型（示例: Fe3Sn 的理想化胞; 参数按 §1 的表替换）
python scripts/build_kagome_cases.py tests/fixtures/Fe3Sn_POSCAR \
    --element Fe --nonmag Sn --moment 2.0 --spin 1.0 \
    --kagome-supercell 3 3 1 --outdir <outdir>/02_supercell/cases --selftest

# 2) 装配 DFT 输入（驱动 dft 步; 产物落在 <outdir>/03_dft/）
python run_2d_magnet.py --model fe3sn --outdir <outdir> ... dft   # 见驱动 --help 的实际子命令
```

产物:

```
<outdir>/03_dft/INCAR.static      # 静态 SCF 模板
<outdir>/03_dft/INCAR.relax       # 弛豫模板（默认带 ISYM=0）
<outdir>/03_dft/KPOINTS           # 默认 3 3 1（可用 --kpoints 覆盖）
<outdir>/03_dft/POTCAR.spec       # 只是"该用哪些 POTCAR"的清单, 不是 POTCAR
<outdir>/03_dft/runs/<order>/{INCAR,KPOINTS,POSCAR}
```

坑（都是**已知**的, 别把它们当成自己的错）:

1. **`ISYM` 默认只写在 `INCAR.relax`（=0）**; `INCAR.static` 里默认没有 `ISYM` —— 但装配现在支持
   `--isym 0`，直接加上即可（非共线计算必须 `ISYM=0`；至少也要 `ISYM=1` 且确认 VASP 支持该对称性）。
2. **`NPAR` 默认 4**（可用 `--npar 8` 覆盖）: 大胞按节点/核数调整, 别照抄默认。
3. **POSCAR 来自 fixture, 不是弛豫结果**: 要么先自己弛豫（`INCAR.relax`）并把 `CONTCAR` 作为
   后续起点, 要么在记录里**明确声明**"用理想化胞, 不弛豫", 并把它写进局限。
4. **POTCAR 不生成**: 按 POSCAR 的元素顺序自己拼（顺序错不会报错, 只会静默出错）。
5. **仓库不提交作业**: 拷到集群、排队、监控都是你的事。

---

## 4. SCF 成功判据（设置 ≠ 判据）+ 失败协议

`INCAR.static` 固定 `PREC=Normal, EDIFF=1E-6, ALGO=Normal, NELM=200, ISMEAR=0, SIGMA=0.05`
（`assemble_dft.py:47–54`）。**这些只是设置。** 每个 OUTCAR 都要逐条检查:

| 判据 | 怎么看 | 不满足时 |
|------|--------|----------|
| (a) 电子步**真的收敛** | OUTCAR 里到达 `reached required accuracy` / 电子步数 < NELM | 未收敛 ⇒ 见下面的协议 |
| (b) 最后两次能量差 | `\|ΔE\|` 最后两次自由能差 ≲ 1e-4 eV | 同 (a) |
| (c) 磁矩与设定方向一致 | 非共线: 检查各原子磁矩的**分量符号/取向**与 `magmom.txt` 对应 | 序变了 ⇒ 这组能量不能进反解 |
| (d) 没有致命告警 | `ZBRENT`/`DENTET`/`WARNING: Sub-Space-Matrix` 之类的字样 | 逐条判读, 记录在案 |

**失败协议（按顺序试, 每次只改一件事并留痕）**:
1. 提高 `NELM`（200 → 400）; 减小 mixing（`AMIX/BMIX`）;
2. 确认 `ISYM=0`（非共线必须）, `LMAXMIX=4`（有 U/f 时必须）;
3. 换 `ALGO`（Normal → All/Fast）; 用 `ISTART=1/ICHARG=11` 从旧 `WAVECAR/CHGCAR` 续算;
4. 仍不收敛 ⇒ **停下来记录**, 不要"跑一个差不多的数"进反解。

> 反解对能量噪声很敏感（cond≈98 还算温和, 但残差闸门是 1e-4 eV）:
> 喂进去一个没收敛的能量, 得到的不是"差一点的 J", 而是**错的 J 且看不出来**。

---

## 5. 回收（最容易出错的一步）

```bash
# 每个序的 OUTCAR 必须回到这一层（cases/<order>/OUTCAR）
cp <cluster>/<order>/OUTCAR <outdir>/02_supercell/cases/<order>/OUTCAR

# 然后二选一
python scripts/solve_J_kagome.py <outdir>/02_supercell/cases
python run_2d_magnet.py ... solve        # 驱动要求 <outdir>/02_supercell/cases 存在
```

- 驱动的 `solve` 只读 `outdir/02_supercell/cases`, 缺目录直接报错:
  `"solve 需要 02_supercell/cases (先跑 build 或全流水线, 再跑 VASP 填 OUTCAR)"`。
- **收回通道已自动化**（本轮新增）: 见 §5.0（交接包）与 §5.1（ingest 硬闸门）。以前的
  「手工 `cp`」写法已作废 —— 手工 `cp` 是整条链上**最危险**的一步: 未收敛 / 混口径 / 错序 /
  旧文件充数的 OUTCAR 会被**静默**当成好数据, 而 J 看起来完全正常。
- 能量口径: `dft_energy.read_toten` 取 **`free energy TOTEN` 的最后一次出现**;
  混口径（不同 INCAR 的 `sigma→0`/`energy(sigma->0)` 写法）会被 `mixed_convention` 挡下。
- **只用总能量**: 3 壳模型只吃 4 个总能量; 力/应力/磁各向异性都不参与反解。

---

## 5.0 交接包与提交脚本（`dft_handoff.py pack`）

```bash
python scripts/dft_handoff.py pack --root <outdir> --out <handoff_dir> [--scheduler auto|slurm|pbs|sge|lsf|none]
# 把 <handoff_dir> 上传到集群；在集群上 cd 进去, bash submit.sh
```

产出: `README_HANDOFF.md` / `MANIFEST.json`（逐文件 sha256 + 期望回收物）/ `POTCAR.spec` /
`CHECKLIST.md` / `submit.sh` / `jobscript.template` / `runs/<order>/{INCAR,KPOINTS,POSCAR}`。

- `submit.sh` 自动识别 **SLURM / PBS / SGE / LSF / 裸机** 五种调度器, 按元素顺序 `cat` POTCAR,
  逐个序跑 VASP, 跑完自查「OUTCAR 在不在、够不够大」, 再调 `ingest_dft.py --dry-run` 过一遍闸门。
- **三条边界不许越**（`--selftest` 里逐条钉住）:
  ① **不打包 POTCAR**（license 不允许随包分发, 只给元素顺序 + 建议 PAW 的说明）；
  ② **不代提交** —— 本脚本只生成脚本, 提交动作由人在集群上做；
  ③ `pack` 对 `runs/` 做**完整性预检**, 缺 INCAR/KPOINTS/POSCAR 的序直接**拒绝打包**
     （免得人跑到集群才发现少文件）。

---

## 5.1 收回通道：`ingest_dft.py` 的硬闸门

```bash
python run_2d_magnet.py ingest --outdir <outdir> [--runs <从集群拷回的目录>] [--ingest-dry-run]
# 等价直调: python scripts/ingest_dft.py <runs_dir> <cases_dir> [--dry-run] [--force] [--lenient <name>]
```

**先 `--dry-run`**: 只校验不落盘, 先把 10 条闸门过一遍, 再决定要不要真收。

| # | 闸门 | 不过说明什么 |
|---|------|--------------|
| 1 | 四态 OUTCAR 齐全（按 `magnetic_orders.json.selected_orders`） | 缺哪个报哪个（`['sq3_alt']` 这种格式） |
| 2 | 体积 ≥ 下限（`--min-kb`） | 很可能是被 kill 的空壳 |
| 3 | `General timing and accounting` 段在场 | VASP 没跑到正常结束 |
| 4 | SCF 收敛（`aborting loop because EDIFF is reached`） | **未收敛的能量不是能量** |
| 5 | **恰好一个离子步**（静态） | 松弛过的结构与 `cases/` 的 POSCAR/系数不对应 |
| 6 | 能量口径统一（转发 `dft_energy.mixed_convention`） | 混口径的能量不可比 |
| 7 | 结构一致性（格矢 + 原子数跨四态逐位比对, 容差 1e-4） | 有一个序其实算的是**另一个结构** |
| 8 | 四态能量差不退化（`--min-energy-gap`） | 四个序实际上是同一个构型（MAGMOM 没生效） |
| 9 | `magmom.txt` 与 OUTCAR 的磁矩块对得上 | MAGMOM 被 VASP 忽略 / 序写错了 |
| 10 | `--force` 才允许覆盖已有 OUTCAR | 防「旧文件充数」: 你以为 4 个序齐了, 其实其中一个没重算 |

- **原子性**: 任何一条不过 ⇒ **一件都不写**（不会出现"收了 3 个、第 4 个失败"的半成品 cases/）。
- 每条闸门可用 `--lenient <name>` 单独降级为**警告**（会记进 provenance）。**慎用**: 降级的是
  J 的可信度, 不是流程的便利性。
- provenance 落盘: `cases/ingest_manifest.json`（每个 OUTCAR 的 **sha256**、判据命中情况、
  集群侧信息）。这是「这个 J 是哪批 OUTCAR 算出来的」的唯一凭证。

---

## 6. 反解与可信度

```bash
python scripts/solve_J_kagome.py <cases>            # 真值路径: source=vasp
python scripts/solve_J_kagome.py <cases> --demo --synthetic 5 -0.5 0.2   # 无 DFT 自检
```

- `J_result.json` 的 `confidence` 是**算出来的**, 不是填的:
  `source=vasp` 且 `residual < 1e-4 eV` 且 `cond < 1e3` 才是 `high`; 否则 `low`
  （`j_result.py:63–71`, `RESIDUAL_FAIL=1e-2` 时**拒绝写出**）。
- **残差不为零意味着什么**: 三壳 Heisenberg 解释不了这批能量。先怀疑（按概率）
  ① 壳合并（配位数不是 `(4,4,6)`）② 更远的壳 ③ SOC/各向异性 ④ 有的序没收敛。
  **不许**把 D、J1′、更远壳当拟合自由度去"把残差压下去" —— 那是过拟合, 不是物理。
- **噪声放大是负的（已实测, 无需 DFT）**: 用 4 个能量 + 高斯 σ 扰动 refit 400 次,
  σ_E = 1e-6 / 1e-5 / 1e-4 / 1e-3 / 1e-2 eV 对应 max|ΔJ| ≈ 5.5e-5 / 5.9e-4 / 5.5e-3 /
  5.5e-2 / 0.58 meV ⇒ 放大因子 ≈ **0.05–0.06×**（拟合**衰减**能量噪声）。
  ⇒ 限制 J 精度的**不是**线性代数的条件数, 是 DFT 总能量的**系统性误差**。别再担心 cond=98.51。
- 产物往下游走: `poscar_to_mc3 --jfile J_result.json` → `mc_2d.py`（T_C）
  → `mc_fss.py`（Binder U4 交点）。kagome 的 `--shell-def` 匹配到 **6 键/胞** 已经被
  `DMI_MTI_MEMO.md:171–172` 演示过。

---

## 7. T_C 与验收锚点（M4）

- 有序温度来自 **MC**, 不来自 DFT: `mc_2d.py` 与 `mc_fss.py`（U4 交点定临界温度）。
- **锚点名字**: Fe3Sn 是**铁磁体** ⇒ 报 **T_C**。旧文档与任务里的「T_N」对 Fe3Sn 是命名错误
  （见文件头第 1 条）。
- **锚点已落盘**（本轮新增）: `literature_J.json` 现在**有** Fe3Sn 条目 —— T_C = 705 K（单晶实测;
  该文引用的多晶文献值为 743 K, 两个数都记在 `source` 里, **不许混说**）、a = 5.4604 Å、
  c = 4.3458 Å、2.27 μ_B/Fe、易面 = kagome 面、`ground_state = "fm"`、
  **J1/J2/J3 = `null`（刻意）**。
  - 来源: L. Prodan *et al.*, "Large ordered moment with strong easy-plane anisotropy and
    vortex-domain pattern in the kagome ferromagnet Fe3Sn", arXiv:2302.13810 (2023) /
    DOI 10.1063/5.0155295。核实级别写在 `source` 里（哪些数来自摘要页、哪些来自全文表格）——
    **引用前按那句话再核一次**。
  - ⇒ `solve_J_kagome.py --ref Fe3Sn` 仍会被按 null **明确拒绝**。**这是设计行为**:
    「有锚点」不等于「有 J 锚点」。**不要**为了让 `--ref` 能用而填一个来源不明的 J。
- 加锚点的规矩: 逐条核对来源（论文 + 表格号）, 写进 `_note`, 有效位数 ≤ 2–3 位,
  **不许凭记忆填**, 未核实就留 `null`。
- 单材料范围: `README.md:295` 提到 Fe₃Sn₂ / Co₃Sn₂S₂ / Mn₃Sn, 但**没有任何代码路径或 fixture
  支持它们**（见 §9）⇒ 本手册只承诺 **Fe3Sn 一条链**。

---

## 7.1 ⚠ 2D MC 的临界温度**不能**直接对 3D T_C（已实测, 这是当前模型类的硬边界）

**为什么必须先读这一节**: M4 的右半判据是「T_C 与文献量级一致」。而这个判据在**当前模型类里
不可能通过** —— 与「有没有机时」无关。理由与实测如下。

**理由**: `mc_2d.py` 是**二维**格点模型（无层间耦合）。二维**各向同性** Heisenberg 模型由
Mermin–Wagner 定理**没有有限温长程序** ⇒ MC 里看到的峰是**有限尺寸伪临界点**, 不是真相变;
而 Fe3Sn 的 705 K 是**三维体相**有序温度, 正是**层间耦合**撑起来的。两者不是同一个量。

**实测**（本机, 不需要任何 DFT）: 拿**已核实文献 J** 的方形晶格反铁磁 **La₂CuO₄**
（J1 = −111.8 meV, S = 1/2, 实测 **3D** T_N = 325 K）跑有限尺寸 χ 峰扫描
（`mc_2d.py --seeds 2 --cool`, T = 60…330 K, 28 点）:

| L | χ 峰 (K) | T_peak / T_N(325 K) |
|---|----------|---------------------|
| 12 | 290 | 0.89 |
| 18 | 260 | 0.80 |
| 24 | 250 | 0.77 |
| 30 | 240 | **0.74** |

- **峰位随 L 单调下降** ⇒ 没有收敛值; 「换个 `--mc-size` 就换个 T_C」。
- 与仓库既有记录自洽: 仓库既有记录本项目单层 χ 峰 **266.67 K** vs 325 K = **0.82**,
  正落在上表区间内。
- 论文级佐证: 二维各向同性 Heisenberg 的 χ 峰位置按 `T_peak(L) ~ 2πρ_s / ln(cL)` 漂移
  ⇒ 该模型下这个量根本不是热力学量。

**⇒ 任何 P4.2 结论必须同时声明**: ① 报的是**有限尺寸伪临界点**, 不是 T_C; ② 它**随 L 变**;
③ 它**不能**与 705 K 比。**要通过 M4 的右半, 必须先改模型**: 显式加入层间耦合 J_c
（Fe3Sn 是双层堆垛, `c` 方向近邻就是层间）, 再用 U4 交点 + 尺寸外推给收敛值。这是**模型
层面的缺口, 不是机时缺口** —— 先补模型, 再开机时。

---

## 5.2 无 VASP 的 dry-run（把整条链在本地空跑一遍）

**目的**: 上真实机时之前, 用**合成 OUTCAR** 把「装配 → 交接包 → 收回闸门 → 反解 → 凭证」
整条链走通。这样等到真机时, 失败的只可能是物理, 不会是**管道**。

```bash
# ① 建构型 + 装配输入 (不碰 DFT)
python run_2d_magnet.py --model fe3sn --skip-dft --skip-mc --outdir _dry_fe3sn
# ② ③ ④ ⑤ 一条命令跑完: 合成 OUTCAR -> 闸门预演 -> 真收回 -> 反解 -> 往返校验
python scripts/dry_run_kagome.py _dry_fe3sn
python scripts/dry_run_kagome.py _dry_fe3sn --synthetic 5 -0.5 0.2 --mc --keep
```

`dry_run_kagome.py` 会: 按 `coefficients_per_cell` 合成 4 个 OUTCAR（只含闸门会查的标志）
→ `--dry-run` 预演 11 条闸门 → 真收回写 manifest → 反解 → 断言**往返误差 < 1e-6 meV**。
任一步不过就非零退出并指出是哪个序、哪条闸门。

**第二次跑要加 `--force`**（不是 bug）: 合成文件是**确定性**的, 第二次必然与 `cases/` 里已有的
逐位相同 ⇒ 防陈旧闸门 `not_shadowed` 会（正确地）拦下 —— 它分不清"你又跑了一遍 dry-run"和
"你忘了重算 VASP"。这正是它在真数据上要防的事, **别把这条闸门改松**; dry-run 里用 `--force`
显式放行即可。

**产物自曝身份**: `J_result.dryrun.json` 的 `source` 被强制写成 `demo`（⇒ `confidence=low`,
`j_result.validate` 不会报"confidence 自相矛盾"）, 且 `notes` 里带「合成能量, 不是 DFT」的
警告 —— **不许**当 `source=vasp` 的结果引用。`--mc` 会顺带跑一次 `mc_2d.py`, 只为验证管道
（T_C 是合成 J 的伪临界点, **无数值意义**）。

**本地环境的已知限制**: 驱动与 `dry_run_kagome.py` 都会 `os.makedirs(outdir)`;
在**受限沙箱**里 cwd 下建目录会被拒（`PermissionError: [WinError 5]`）。把 `--outdir`
指到**系统临时区**（如 `$env:TEMP\...`）即可 —— 这也是仓库既有测试用 `tempfile` 的原因
（见 `tests/smoke_test.py:81-83` 的注释）。

**仓库里还有这个 dry-run 的契约测试版**（逐项断言, 适合放进 CI）:

```bash
python tests/test_ingest.py          # 闸门/原子性/--lenient/dry-run 的契约
python tests/test_kagome_cases.py    # ingest → solve 的凭证闭环 + J 逐位还原
```

**dry-run 能证明什么 / 不能证明什么**:

| 能证明 | 不能证明 |
|--------|----------|
| 四态装配与设计矩阵自洽（rank=4, cond=98.51） | J1/J2/J3 的**数值**（合成值, 不是 DFT） |
| 交接包完整、提交脚本语法正确、POTCAR 不被误打包 | 集群上 VASP 真能跑（二进制/POTCAR/机时都不在本仓库） |
| 11 条闸门**确实会拦**坏数据（未收敛/混口径/错结构/旧文件） | 真实 OUTCAR 的格式细节（合成文件只含闸门会查的标志） |
| 反解管道 + 凭证链（manifest → J_result.json）通畅 | T_C 的任何数值（见 §7.1: 模型类就不支持） |
| 报错文案能定位到具体序与具体闸门 | 机时估算（§8 仍是估算） |


---

## 5.3 k 点收敛抽查（`k_convergence.py` / `run_2d_magnet.py kconv`）

**为什么单独一步**: `KPOINTS` 默认是**假设**（`3 3 1`）。J 来自**四态能量差**，所以「k 网格够不够密」
不该靠猜、也不该靠「总能量看着收敛了」就宣布 J 收敛 —— 这一步把「多密才够」变成可查的产物。

```bash
# 1) 生成变体：每一序 × 每个网格一个目录（只换 KPOINTS，不提交作业）
python run_2d_magnet.py kconv --outdir <outdir> --k-grids 3,4,6,8 [--k-orders fm,q0_120]
# 等价直调: python scripts/k_convergence.py --root <outdir> --grids 3,4,6,8 [--plan-only|--ingest-only]
```

- 变体落在 `kconv/runs/<序>@<k网格>/{INCAR,KPOINTS,POSCAR}`：**INCAR 逐字节复制**，**只换 KPOINTS**
  （换的是唯一变量；顺带把 `--kpoints` 的默认假设替换成显式网格，见 §3）。
- 在集群上跑完，把 OUTCAR 放回**同一个目录**，再跑同一个命令（不带 `--plan-only`）：
  按 `--min-kb`（默认 200 KB）/ `General timing and accounting` / SCF 收敛三条闸门筛 OUTCAR，
  算每个网格相对**最密**网格的 `ΔE`（meV/atom），判 `converged` / `not_converged`，
  并给出 **`recommended_kpoints` = 满足 `--k-tol`（默认 1.0 meV/atom）的最粗网格**。
- 产物：`kconv/kconv_plan.json`（计划 + 逐目录清单）与 `kconv/k_convergence.json`
  （逐网格能量、ΔE、判定、推荐值、以及下面这条免责声明）。
- ⚠ **边界（引用时一起说）**：① 它判的是**总能量**对 k 网格的收敛，**不是 J 的收敛** —— J 是能量**差**，
  两个序各自收敛**不保证**它们的差收敛（差的收敛通常更快，但那是**经验**，不是定理）；
  ② 它**不代提交作业**，也不替代 §5.1 的闸门（判 VASP 是否成功是 `ingest_dft.py` 的事）；
  ③ 本仓库**没有任何真实 kagome OUTCAR** ⇒ 这一节没有可引用的数字，只有工具与判据。

---

## 8. 机时估算（**纯估算: 假设 × 数量, 未实测**）

**假设块**（改任何一个都要重算）:

| 假设 | 取值 | 说明 |
|------|------|------|
| 原子数 N | **36** | 3×3×1 超胞（`n_atoms=36` = 27 Fe + 9 Sn, 已实测确认。旧版写 75 是错的） |
| ENCUT | 520 eV | 仓库默认 |
| k 网格 | 3×3×1 Γ 中心（可用 `--kpoints` 覆盖） | 仓库默认; 收敛抽查工具见 §5.3, **未实测** |
| 核数 | 32 | 假设值 |
| 单 SCF 时间（非共线, 无 SOC） | **1–3 h** | 经验区间, **不是本机实测** |
| 加 U / f 电子 | ×2–4 | 经验系数 |
| 若加 LSORBIT（SOC） | ×1.5–3 | 经验系数 |
| 规模标定 | N²–N³ | 换超胞时按此放大 |

**同族参照（不是本机实测, 只是量级对照）**: 已核实的 Fe3Sn 文献（Prodan *et al.*
arXiv:2302.13810 / DOI 10.1063/5.0155295）用的是 **ENCUT 800 eV** + Γ 中心 **10×10×10
（弛豫）/ 20×20×20（精确总能量）** + `EDIFF = 1e-8 eV`。仓库默认的 520 eV / 3×3×1
**远低于**该文标准 ⇒ 本节的机时与精度都必须按 §7.1 的结论重新审: 在把 520 eV/3×3×1
的结果当结论之前, **至少要把 ENCUT 与 k 网格各收紧一次并留痕**。

**数量**: `selected_orders` 4 个 ⇒ **4–12 核时×h 量级**（= 4 个 SCF）;
`--all-orders` 7 个 ⇒ 约 1.75×。再加上 §1 里要求的 k/ENCUT 抽查与一次 FM 参考点。

**标定规则（必须做, 否则本节的数字不许引用）**:
1. 先只跑 **`fm` 一个 SCF**（参考点）, 记录: 原子数、ENCUT、k 网格、核数、wall-clock、INCAR 全量;
2. 用实测值**替换**上表的"1–3 h"并写上**日期**, 同时把 `N²–N³` 标定换成实测的两点标定;
3. 在替换之前, **不许把本节当作预算依据**（与 `DMI_MTI_MEMO.md:266–281` 同一条规矩）。

---

## 9. 与原计划的偏差（本手册承认的）

1. `LATTICE_ROADMAP.md:102` 把当年的 item 20 写成"**J1/J1′/Jd 三自由度分离**", 而实际交付的是
   **三个壳**（J1 = J1 与 J1′ 的平均）⇒ 本节 §2 明确**禁止**报"分离出了 J1/J1′"。
2. `LATTICE_ROADMAP.md:154` 只有一个 "P4 7.5–12 人天" 行, 切分后 **P4.2 的残余工作量从未估过**。
   本手册按操作口径给一个**人天级**估计（不含机时等待）: 装配与抽查 ≈ 1–2 人天;
   提交与监控（含不收敛返工）≈ 2–5 人天; 回收 + 反解 + MC + 对文献 ≈ 1–2 人天
   ⇒ **≈ 4–9 人天**, 不确定度 **高**（依赖机时可得性与参数抽查次数）。请把这一行回填到计划表。
3. `README.md:295` 暗示 3 个材料（Fe₃Sn₂ / Co₃Sn₂S₂ / Mn₃Sn）的范围, 代码只支持 Fe3Sn
   ⇒ 本手册把范围**收窄到单材料**; 换材料 = 加 fixture + `--model` 预设 + 重跑 builder。
4. `README.md` 原写"未开工"、`LATTICE_ROADMAP.md:154` 原写"待机时", 两边都没有**入口条件**
   （2026-10 回填: README 现已改成"物理层 = 手册 + 等真实机时", 该行已一并按本手册口径修好）:
   本手册就是入口定义 —— **入口条件 = 有真实机时 + 已核实文献锚点 + §1 的表填满**。
5. 任务措辞常把 P4.2 说成"非共线 **SOC** DFT"。按代码, 现有 7 个序**只需要 `--noncollinear`**;
   LSORBIT 只在你要分离 J1/J1′ 时才需要 —— 不加 SOC 可以把机时门槛降到 1/1.5–1/3（§2）。

---

## 10. 本仓库做不到的（不要写进论文/汇报）

- 不提交作业、不写 POTCAR、不提供机时;
- **J1 与 J1′ 的分离**（需要新构型 + SOC, 尚未实现）;
- **DMI / MTI 的定量结论**（P6 只做识别层, 见 `DMI_MTI_MEMO.md`）;
- **任何"预测 T_N"的说法**: 三壳 Heisenberg + 单一 MC 参数化只能给量级与趋势,
  不许给"预测值 ± 误差";
- 没有实测机时之前, §8 的数字不是预算。

---

## 11. 检查清单（打印出来照着勾）

**上机之前**
- [ ] 化合物 / a / S / μ / U 已定, **每一项都有来源**, 记进 §1
- [ ] 知道 Fe3Sn 报的是 **T_C 不是 T_N**（文件头第 1 条）, 且**没有**把 FeSn/Fe3Sn₂ 的数混进来
- [ ] 读过 §7.1, 并已决定: **要么先加层间 J_c 再上机时**, 要么接受"这次只交 J、不交 T_C"
- [ ] ENCUT 与 k 网格的收敛抽查**已排进计划**（文献对照是 800 eV / 10×10×10–20×20×20, 见 §8）
- [ ] k 网格抽查**已生成变体**（`run_2d_magnet.py kconv --k-grids 3,4,6,8`, 只换 KPOINTS）,
      回收后看过 `recommended_kpoints`（§5.3）, 且知道它**只**证明总能量收敛
- [ ] (无 VASP 时) 手工走一遍 **§5.2 dry-run**, 确认交接包与闸门链路通畅

**装配与交接**
- [ ] POTCAR 已按 POSCAR 顺序拼好（版本也记下来）
- [ ] 弛豫完成, 或明确声明"用理想化胞、不弛豫"
- [ ] `build_kagome_cases.py --selftest` PASS; 配位数 `(4,4,6)`
- [ ] 4 个 `selected_orders` 全部装配; `--noncollinear`（如需分离 J1/J1′ 才加 `--lsorbit`）
- [ ] `dft_handoff.py pack` 出的交接包 SHA256 清单齐全, 且 **不含 POTCAR**

**跑完与回收**
- [ ] 每个 `runs/<order>/OUTCAR` 都满足 §4 的 (a)(b)(c)(d)
- [ ] `run_2d_magnet.py ingest --ingest-dry-run` 先过一遍闸门, 再真收回
- [ ] `cases/ingest_manifest.json` 存在, 且 `checks` 里**没有**被降级的硬闸门
      （如确有降级, 在引结论时明写）

**反解与磁学**
- [ ] `solve_J_kagome.py` → `J_result.json`: `source=vasp`、`confidence=high`、配位数闸门 PASS
- [ ] `J_result.json` 里有 `dft_provenance`（指向 manifest + 逐序 sha256）
- [ ] `poscar_to_mc3 --jfile` → `mc_2d.py` / `mc_fss.py`
- [ ] 报 T_C 时**同时**声明: 有限尺寸伪临界点 / 随 L 变 / 不能直接比 705 K（§7.1）
- [ ] 结果表述符合 §10（不报 J1/J1′ 分离、不报 DMI、不报"预测值"）

**留痕**
- [ ] 实测机时写回 §8（含日期）; 本文档的"估算"字样被替换
- [ ] 任何新核实的文献值写进 `literature_J.json` 并带出处（**J1/J2/J3 仍是 null 就别动**)

---

## 12. 相关文件

| 文件 | 关系 |
|------|------|
| `LATTICE_ROADMAP.md:93–105, 154, 177–178, 183` | P4 的定义、切分理由、M3/M4 判据、入口条件 |
| `README.md` · `LIMITATIONS.md` 第 9 条 | 命令示例与 kagome 物理层的边界声明 |
| `WORKFLOW.md:218, 228` | P4.1 只交付系数层、局限 1e |
| `skill/2d-magnet/SKILL.md:314–334, 360–367` | kagome 分支的操作纪律与局限 |
| `scripts/dft_handoff.py` | §5.0 交接包与提交脚本（`--selftest` 18 项） |
| `scripts/ingest_dft.py` | §5.1 收回通道的 11 条闸门（`--selftest` 18 项） |
| `scripts/k_convergence.py` | §5.3 k 点收敛抽查（生成 `<序>@<k网格>` 变体 + 判定 + 推荐网格; `--selftest` 29 项） |
| `tests/test_ingest.py` | §5.1/§5.2 的契约测试（闸门硬软分工、原子性、--lenient、dry-run） |
| `DMI_MTI_MEMO.md:266–281` | 机时估算的**形式**模板（不是数字） |
| `templates/new_lattice_checklist.md` | 新晶格清单（换材料时用） |
| `tests/test_kagome_cases.py` / `tests/smoke_test.py` §23 | mock 级回归（含 ingest → solve 的凭证闭环） |
