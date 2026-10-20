# 2d-magnet 新结构类型「进主链」落地计划

范围：把技能文档新增的 6 类结构（六角 kagome、体相 DMI、磁性拓扑绝缘体、交变磁体、方形 K2NiF4、矩形/oxyhalide）从「识别 + 路由」推进到「可计算主链」，或明确其只能停留在识别层。
本文只出计划与工期，不含代码。

> **时效声明（默认全共线）**：本文是历史计划，其中「120° 非共线构型」相关段落
> （§P3、§P4 里 `selected = fm/q0_120/...`、cond≈98.51 等）是**当时的快照**。
> 现在三角与 kagome 的默认路线都是**全共线**（`rank=4, cond≈3.7e3 / 101`），
> 120° 只在显式 `--include-noncollinear` / `--noncollinear` 的对照支线里出现；
> 当前口径以 [`../skill/2d-magnet/SKILL.md`](../skill/2d-magnet/SKILL.md) §2.0 与
> [`WORKFLOW.md`](WORKFLOW.md) 为准。

## 0. 勘察结论（计划的前提）

> **时效声明（2026-10 回填）**：本节是立项时的勘察快照，保留原文以便对照。已被后续工作取代的三条：
> ① 第 1 条的 `--order` 当时只有 `auto/staggered/total/cv`，现在还有 `sublattice:<mask>`（P5，见 §P5）；
> ② 第 2 条的 `poscar_to_mc3.py` 现在有 `--shell-def`（逐壳 3×3 `J_mat`）与 `--jfile` 通道，但 `--element` 仍只接受一个磁性元素（P4.1）；
> ③ 第 3 条的 `assemble_dft.py` 现在有 `--noncollinear/--soc/--lsorbit/--saxis`，且对需 SOC 的磁序缺 `--lsorbit` 会**硬报错**（P3，见 §P3）。第 4/5 条与现状一致。

1. **MC 引擎已经足够通用，不是瓶颈。** `2d_magnet_pipeline\scripts\mc_2d.py` 的 `build_arrays()`（mc_2d.py:28-62）直接消费模型 JSON 的 `primitive_vectors` / `basis_frac`（任意多子晶格）/ `onsite`（逐子晶格各向异性）/ `interactions[{b1,b2,offset,J(3×3)}]`，且 `J` 与 `J.T` 分别下发（mc_2d.py:53-54）⇒ **各向异性交换与反对称项（DMI）在数据结构上已可表达**。`--order` 仅 4 档：`auto/staggered/total/cv`（mc_2d.py:264），Binder U4 与 `mc_fss.py` 有限尺寸标度已就绪。
2. **真正的瓶颈在 J 提取侧的「脚本对」，以及 `poscar_to_mc3.py` 的壳→J 映射。** `poscar_to_mc3.py:226-247` 把**距离排序**的第 1/2/3 个壳硬映射为 J1/J2/J3，且一律写成各向同性 `-(J·S²)·I₃`（poscar_to_mc3.py:233-236）；`--element` 只接受**一个**磁性元素（poscar_to_mc3.py:209）。没有 per-shell / per-sublattice-pair / 各向异性 J 的输入通道。
3. **DFT 模板完全没有 SOC。** `assemble_dft.py:43-51` 只有 `ISPIN=2 + LDAU + LMAXMIX=4`；`assemble_dft_tri.py:40-48` 有 `LNONCOLLINEAR=.TRUE.` + 3 分量 MAGMOM，但**没有 `LSORBIT`**。⇒ 任何依赖 DMI / 弱铁磁倾斜 / 拓扑表面态的物理，缺的是 DFT 前置能力，不是 MC 能力。
4. **多子晶格 + 逐壳 J 已有现成模板**：`Fe3GeTe2_pipeline\05_mc\build_fgt_mc_model.py`——`SHELLS` 表按 `(name, pair∈{AA,AB,BB}, d_target, tol, 文献 J, 来源注)` 逐壳给定，用 canonical half-space 去重列键，注释明确 `J_mat = -(J_k·S²)·I₃`、`A = diag(0,0,-MAE)`。这正是应该泛化进主仓的东西。
5. **现有 5 条主链的公共契约（脚本对模式）**：
   - `build_*_cases.py`：POSCAR + 超胞 + 元素/矩 → `cases/<磁序>/{POSCAR,magmom.txt}` + `cases/magnetic_orders.json`（`n_mag / S2 / coefficients_per_cell / coefficients_per_mag / shells_angstrom / selected_orders / noncollinear / inversion_matrix_per_cell`）+ `--selftest`；
   - `solve_*_J.py`：读 `magnetic_orders.json` + 各 `cases/<磁序>/OUTCAR` 的 `free energy TOTEN`，解 A·J̃ = E（`np.linalg.solve`），`J = J̃/S²·1000`（meV），自带 `--demo` 合成能量自检；
   - `render_*`（出图）→ driver `run_2d_magnet.py` 串联 → `poscar_to_mc3` + `mc_2d` → `mc_fss` 精修 → `ship-mc` → `write_hamiltonian`；
   - 回归由 `tests\smoke_test.py` 15 个 section 覆盖（14 = 矩形 CrSBr，15 = 4f GdI2；13 节把「`MODEL_PRESETS` 必须与 `literature_J.json` 材料专属值一致，不得与 `ELEMENT_PARAMS` 通用兜底混淆」钉死），`tests\test_fgt_e2e.py` 覆盖 FGT 双亚晶格 mock。

## 一、判据与六类判定

### 1.1 六道闸门（进主链的判据）

| 闸门 | 内容 | 不通过的后果 |
|---|---|---|
| G1 能量可表示 | 磁能量能写成少量两体交换项之和，且对未知 J **线性**（可用 A·J̃ = E 求解） | 需要双二次/环交换/高阶级数 ⇒ 拟合病态 |
| G2 磁序可枚举 | 有限个共线/非共线磁构型可枚举，且系数矩阵**可解**（条件数可控、残差可校验） | 解不出唯一 J，或解对输入能量噪声过敏 |
| G3 DFT 可交付 | 有可用的 INCAR 模板（含所需 SOC / 非共线 / U 设置），胞尺寸 × k 网格 × 构型数在机时预算内 | 论文级物理但工程上跑不动 |
| G4 MC 序参量可捕获 | 目标 T_C/T_N 的序参量能被现有 `--order`（或小改）定义 | 算出 J 却无法转成主链唯一交付物「临界温度」 |
| G5 有验证锚点 | 有独立文献/实验 T_C 或 J 可对照 | 无基线，eval 无法判定对错 |
| G6 契约可对位 | 参数能落进现有 CLI 契约（element/S/J1..J3/MAE → `cases/` → T_C），产物目录约定不破坏 | 产出物不属于主链契约，只能另立子链 |

判定规则：**G1–G4 全过 + G5 有锚点 + G6 可对位 = 进主链**；缺 G4 或 G3 只需小改 = **半进主链（可算但降级交付）**；G1–G3 任一需要新子链 = **只识别 + 路由**。

### 1.2 六类判定

| 类型 | G1 | G2 | G3 | G4 | G5 | G6 | 判定 | 关键依据（一句话） |
|---|---|---|---|---|---|---|---|---|
| 方形 K2NiF4 | ✅ | ✅ | ✅（共线，模板现成） | ✅（`staggered`） | ✅ | ✅ | **进主链（P1，✅ 已完成）** | 单磁性子晶格 + J1/J2/J3 三个距离壳，a=b 时 x/y 键长简并、合并到同壳**物理上正确** |
| 矩形 / oxyhalide | ✅ | ✅ | ⚠（弱倾斜需 SOC） | ✅ | ⚠ | ✅ | **进主链（P2，✅ 已完成）**；弱 FM 倾斜 → **P2′ 已落地**（系数层 + 反解器；θ 是输入参数、无真实 SOC 能量，见 `LIMITATIONS.md` 11 与 §P2′） | 原判据（「复用已通的 `build_rectangular_cases.py`，增量只是元素/化学计量预设」）**被实测推翻**：CrOCl 是 **(½,½) C 心矩形网、2 个磁原子/胞**，且两 Cr 只差 3.03 Å < 簇判据 ⇒ 实际做法 = 修 `detect_lattice` 簇判据 + 新写 `build_crect_cases.py` / `solve_J_crect.py`（详见 §P2） |
| 六角 kagome | ✅ | ⚠（q=0 与 √3×√3 需非共线） | ⚠（7 个内置序只需 `--noncollinear`；分离 J1/J1′ 才需 SOC） | ⚠（序参量需核对） | ✅ | ✅ | **系数层已进主链（P4.1 ✅）；物理层 = 手册 + 等真实非共线 DFT 机时（P4.2）** | 3 位 basis 在 `reduce_to_primitive`/`find_bonds` 下自然成立，卡点是 **J1 与 J1′ 在同距离上被 3 壳模型合并**（取平均）—— 要分离需要更细的构型 + SOC，本轮**不做**（见 `KAGOME_RUNBOOK.md` §2/§9） |
| 交变磁体 | ✅ | ⚠ | ✅ | ✅（P5：`--order sublattice:<mask>` 逐子晶格 Néel 矢量） | ⚠ | ⚠ | **进识别层（P5 ✅）：几何判别 + 标量序参量，不给 k 分辨自旋劈裂结论** | 需要逐子晶格 Néel 矢量/掩码序参量；且「是否交变磁体」的 k 分辨自旋劈裂判别属研究问题，不属本流水线契约 |
| 体相 DMI | ✅（3×3 J 可编码） | ❌（3 分量 D 矢量拟合成本高） | ⚠（P3 已有 `--lsorbit`；仍无自旋螺旋构型生成器） | ❌（无拓扑荷序参量） | ⚠ | ❌ | **只识别 + 路由，不进主链（P6 ✅ 已做识别层 + 判据备忘）** | 即使算完，主链目标函数是 T_C，DMI 产出物对不上契约；须另立 DMI/skyrmion 子链 |
| 磁性拓扑绝缘体 | ⚠ | ❌（需层分辨 A 型 AFM） | ❌（SOC + slab + 层间，成本高） | ❌（单层 MC 无意义） | ⚠ | ❌ | **只识别 + 路由，不进主链** | MnBi₂Te₄ 类 T_N 由层间耦合控制，与「2D 单层」假设直接冲突，需另立层状子链 |

## 二、按依赖排序的改动清单（脚本 + 文档）

### P0 地基（前置，必做；不改物理，只冻结接口）
1. `scripts\magnetic_orders.py`（**新增公共层**）：磁序规格数据（超胞矩阵、子晶格自旋相位、期望系数）+ 由规格构造系数矩阵 A + 求解 + **条件数/残差校验**。首批以三角、矩形、簇、蜂窝、FGT 五条既有链路为回归对象（数值必须与旧路径一致）。
2. `J_result.json` 统一产物 schema：`{lattice, shells:[{name,d,sum_coeff,z_n}], J:{...}, cond, residual, source:"vasp"|"demo"|"literature", confidence}`；五个 `solve_*_J.py` 适配写出。
3. `cases/magnetic_orders.json` 增补 `order_spec[]`（每个构型的子晶格取向与系数期望），五个 `build_*_cases.py` 适配写出。
4. `scripts\poscar_to_mc3.py` 扩展：`--jfile J_result.json`（逐壳/逐子晶格对/各向异性 J）、`--shell-def`（显式壳定义，绕开距离排序）、多磁性元素 `--element X Y`。（扩展点在 poscar_to_mc3.py:209-212、:226-247、:249-279）
5. `tests\smoke_test.py` 增补「脚本对契约测试」section：每个 build/solve 对的最小回归 + schema 合规检查。

### P1 方形 K2NiF4
> **状态: ✅ 已完成**（`sq_k2nif4` / `sq_la2cuo4` 两条规格；`--regress` 7/7）。实施中的两处偏离计划：
> ① 第 8 项**不新增** `render_square.py`，直接复用 `render_supercell.py`（它已能从 `magnetic_orders.json` 读壳距离与 `coefficients_per_mag`）；
> ② 第 9 项除了方形分支，还必须**先补 `_is_square_net` 判据**——否则 a=b 的方形会被 `detect_lattice` 误路由到 triangular 去解 120° 非共线。
> 另外发现两个**物理陷阱**（计划里没写）：周期 2 的序一律给出相同 c3 ⇒ 不引入周期 4 的 `p4_a` 就无法把 J3 与 e0 分离；J3=2a 在 4 宽超胞里是**半周期**，必须用 `policy="shift_search"`，否则最小镜像把 ±2a 折叠、J3 低估 2 倍。
6. `scripts\build_square_cases.py`（新增）：方形超胞 + G 型(Néel)/stripe/FM 磁序枚举 + 系数矩阵 + `--selftest`。→ **已落盘**（6 共线态 fm/g_type/stripe_a/stripe_b/p4_a/mix_ab，DFT 只提交选定 4 个）。
7. `scripts\solve_J_square.py`（新增）：复用 P0.1 规格层，本体基本只剩规格 + 求解调用。→ **已落盘**（`--demo` / `--ref la2cuo4` 均 PASS，写 `J_result.json`）。
8. `scripts\render_square.py`（新增，可并入 P0 的通用 render）。→ **跳过**：复用 `render_supercell.py`。
9. `run_2d_magnet.py`：方形类从路由分支切到 compute 分支；补 `ELEMENT_PARAMS` / `MODEL_PRESETS` 条目（须与 `literature_J.json` 一致）。→ **已完成**（`la2cuo4` / `k2nif4` 两个预设 + `--lattice square` + `--square-supercell`）。
10. `tests\smoke_test.py` 新增方形 section（同构于现有 14 节矩形）+ ≥20 个磁构型回归用例。→ **已完成**（smoke §18 共 24 项 + `tests\test_square_cases.py` 28 项三方交叉锚定）。
11. 技能文档：`SKILL.md` 路由表、`references\magnet-types.md` 的「可算/不可算」标注、`WORKFLOW.md` 阶段产物。→ **已完成**（另加 `references\magnet-catalog.md` 的 §B₃ 与 README）。

**M2 端到端实测**（无 VASP，走文献 J 正向 MC）：`--model la2cuo4` 判别为 `square` → 4×4×1 / 112 原子 / 16 Cu / cond(A)=40.41 → MC(24×24, 50–700 K, nt=25) 得 χ 峰 **T ≈ 266.67 K**、|dm/dT| 极大 **T ≈ 239.58 K**、Binder U₄ 低温 0.3997 → 高温 0.0005。与 La₂CuO₄ 实测 **T_N = 325 K** 之比 0.82——量级符合，但 T_N 是 3D 有序温度、单层 2D 有限尺寸峰本就该更低（局限已写进 README/WORKFLOW）。

### P2 矩形 / oxyhalide（复用现有 rect 对）
> **状态: ✅ 已完成**（`crect_crocl` 规格；`--regress` **8/8**）。但**计划的前提被实测推翻**，实施有四处在原计划之外：
> ① **oxyhalide 不是「复用现有 rect 对」**：CrOCl 单层（JARVIS JVASP-27906, Pmmn No.59）两个 Cr 占 (0,0,0.534)/(0.5,0.5,0.466)，是 **(½,½) C 心矩形网、2 个磁原子/胞**，不是 CrSBr 那种 1 磁原子/胞的简单正交网。因此**没有**去小改 `build_rectangular_cases.py`（它假设 1 磁原子/胞，`--selftest` 理想值表会报 `系数最大偏差=2.000 FAIL`），而是新写 `build_crect_cases.py` + `solve_J_crect.py`。
> ② **必须先修 `detect_lattice` 的簇判据（真 bug）**：两个 Cr 只差 3.0347 Å < 簇判据 3.3 Å，旧版 `if n_clusters < n_mag: return "cluster"` 会把它当分子簇磁体（Nb₃I₈ 路径）。新增返回 `rect_centered`（正交网并簇）与 `sq_centered`（方形网并簇，暂不支持 ⇒ **显式报错**，不静默套方形路径）。改后 Nb₃I₈ / Fe₃GeTe₂ 的 `('cluster',1,3)` 逐位不变。
> ③ **顺带修两个既有真 bug**：(a) 两个 builder 的 MAGMOM 与 POSCAR **元素分组错位**（MAGMOM 必须按 POSCAR 物种首次出现顺序拼；CrOCl 原胞以 `Cl` 开头才暴露，K₂NiF₄/La₂CuO₄ 恰好是 no-op）；(b) driver **无条件要求 `--J1`**，导致文献 J 为 null 的材料连 `build` 都跑不了 ⇒ 改为只在 `mc`/`report` 阶段强制，否则 stderr 提示「本阶段只用几何/规格层」。理由：null-J 材料的结构与磁序枚举本来就该在没有 J 时可用，否则与「不猜数」政策自相矛盾。
> ④ **发现一个跨脚本口径分歧（原计划没写）**：`poscar_to_mc3.py` 在**面内 2D 投影**里搜键（J1 读到 2.5064 Å），`build_crect_cases.py` 用 **3D 键长**（J1 = 3.0347 Å，含 1.71 Å 层起伏）。两把尺子都对但不可混用；**第一版修法（给 builder 加 `dim2` 假投影）是错的**——它保留 3D 的平移集只丢 z，把每个键数了 3 倍（fm 系数 3D [128,64,64] vs 2D [384,192,192]），已全部回退。最终修法是**跨口径不变量校验**：builder 产物新增 `per_site_bonds`（每格点键数 = 4/2/2），driver 的 `verify_mc_shells()` 在 `poscar_to_mc3.py` 之后比对 MC 模型的 `nn_per_site/nnn/nnnn`，不一致就 **SystemExit**（「J1/J2/J3 可能被绑到了错的壳上——距离自动排序的次序与规格层不同」）。这条校验已在 CrOCl(4/2/2) 与 La₂CuO₄(4/4/4) 两条链路上生效。
>
> 另外两个**陷阱**（与 P1 的陷阱不同源）：C 心键连的是**同一胞内**两个 basis ⇒ basis 符号对一切 (n0,n1)-奇偶序**完全不可见**（`afm_c_a == afm_a` 等逐位相同），母序与其 basis 翻倍版**同选必使 det=0**（脚本硬拒）；J1=3.0347 与 J2=3.1980 只差 0.1632 Å ⇒ `magnetic_shells(tol=0.2)` 会**静默并壳**，`build_crect_cases.py` 默认 `--shell-tol 0.05`。弱 FM 倾斜（需 SOC/DMI）显式留给 P3。
12. `run_2d_magnet.py`：oxyhalide 元素/参数预设（Fe/V/Cr 的 U、S、moment）+ 预设自检。→ **已完成**（`crocl` 预设：Cr / nonmag `O Cl` / lattice `rect_centered` / S=1.5 / moment=3.0 / U=4.0；新增 `--lattice rect_centered` 与 `--crect-supercell`；`sq_centered` 显式报错）。
13. `scripts\build_rectangular_cases.py`：层状 vdW 间隙与非磁元素（O/X）占位检查、化学计量校验（小改）。→ **改为新写 `scripts\build_crect_cases.py`**（同号 13 的意图用新脚本实现：要求原胞**恰好 2 个**磁性原子，非磁 O/Cl 由 `--nonmag` 给；超胞前两维须为偶数且 ≥2；`det<1e-6` 与 `cond>1e4` 硬报错；**非磁元素占位与化学计量校验**（`--nonmag` 漏写胞内元素、或写了胞里没有的元素 ⇒ 硬报错，不静默错位；每胞计量必须为整数并落进 `stoichiometry_per_cell` / `formula_per_cell`，CrOCl 为 Cr2O2Cl2）；`--selftest` 8 项 = basis 简并不变式 / fm 与 afm_c 只差 c1 / 退化组合 det 须为 0 / 7 序互异 / 合成 J 反解 / 4×4×1 解析理想值 / 与规格层逐位交叉校验 / cond 闸门）。`build_rectangular_cases.py` 本身除了 MAGMOM 顺序修正外**未被改动**，保持对 1 磁原子/胞材料的兼容。
14. 技能文档：oxyhalide 归入矩形族的说明 + 与 CrSBr 路径的差异表。→ **已完成**（`WORKFLOW.md` §4 步骤表/§5 目录树/§6 继承表/§7 局限 1d/§8 接手清单 + `README.md` 支持列表与预设表与局限 9 + 本节；差异表 = 「同一个正交网家族，但胞内磁原子数 1 vs 2 ⇒ 判簇判据、壳容差、`afm_c` 序三者全部不同」）。
15. `scripts\solve_J_crect.py`（**计划外新增**）：C 心矩形反解，`--ref`/`--synthetic` 互斥且 `--demo` 必须给其一；由 fm 行导出的每格点配位数必须等于 4/2/2，否则 SystemExit 并提示「壳很可能被 `tol=0.2` 静默合并」；`--demo` PASS = 往返误差 <1e-6 meV 且 z 对 且 cond<1e4。CrOCl 的文献 J 未核实 ⇒ 该脚本**拒绝**用 null 或他材料的值当锚点。

**M2 端到端实测**（无 VASP；J 用 `--synthetic` 合成，**只验证管路**）：`--model crocl` 判别为 `rect_centered` → 4×4×1 / 96 原子 / 32 Cr / cond(A)=64.02 / det=524288 → MC(合成 J1=5, J2=−0.5, J3=0.2 meV) 得 χ 峰 **T ≈ 90.00 K**。**这不是 CrOCl 的物理结果**（J 是编的），合成 J 产生的 `J_result.json`/`05_mc/`/`hamiltonian.md` 已从仓库删除以免误读。要真数需走 DFT+MC 或找到可复核的文献 J。

### P3 非共线基础设施（kagome 与 oxyhalide 倾斜的共同前置）
> **状态: ✅ 已完成**（自旋矢量层 `mag_vectors.py` + 统一入口 `assemble_dft.py` + 能量口径 `dft_energy.py` + smoke §21 快照回归）。实施中的四处偏离/补充：
> ① 第 15 项抽出来的不只是「三角的 MAGMOM 生成器」，而是**整个自旋矢量层** `scripts\mag_vectors.py`：共线 / 面内 120° / `±θ` 倾斜 / 螺旋 + 渲染 / 解析 / 校验一体，三角只是它第一个消费者。`ANGLE_TABLE` 记录**已实现**的三类角度约定（`collinear` / `planar120` / `cant` = ±θ 倾斜），`SOC_ORDERS` 则是 **SOC 门禁的名字表**（还含尚无生成器的 `spiral`/`tilt`/`cant_a|cant_b|cant_diag` —— 光加名字只会让门禁生效，要真用必须自己写生成器）。P4.2 的 √3×√3 与 P2′ 的弱 FM 倾斜都复用 `spin_cant`。
> ② 第 16 项的「统一入口」不是给 `assemble_dft_tri.py` 再抄一份 SOC 代码，而是**反转依赖**：它退化成薄适配层，只拼命令行去调 `assemble_dft.py --noncollinear --all-orders`，两条路径从此**共用同一份 INCAR 渲染**——这正是 §104-1 要求的「从复制代码升级为数据驱动 + 薄适配层」。
> ③ 新增一条**不静默降级**门禁（原计划没写）：磁序落在 `SOC_ORDERS`（cant/spiral/tilt…）却没给 `--lsorbit` ⇒ **硬报错**。理由：缺 SOC 时这些构型在物理上不成立，静默降级会产出「看起来正常但物理错」的 J。同理「3 分量 `magmom.txt` + 共线模式」「`--lsorbit` 缺 `--noncollinear`」也都硬报错。
> ④ 顺带修掉一个**既有真 bug**：旧 `assemble_dft.py` 的 `--species` 必填且**顺序敏感**——与 POSCAR 物种顺序不一致时，`LDAUL/LDAUU/LMAXMIX` 与 MAGMOM 分组会整体错位（P2 已在两个 builder 里修过同类 bug，这是第三个出口）。现在渲染一律按 POSCAR 顺序,`--species` 只做一致性校验（不匹配即报错）。
15. `scripts\assemble_dft_tri.py` 的非共线 MAGMOM 生成器**抽成独立模块** + 子晶格相位/角度表（120°、±θ 倾斜）。→ **已落盘**：新模块 `scripts\mag_vectors.py`（`SOC_ORDERS`/`needs_soc`/`ANGLE_TABLE`/`unit`/`spin_collinear`/`spin_planar`/`spin_120deg`/`spin_cant`/`spin_all`/`magmom_values`/`render_magmom`/`expand_repeats`/`parse_magmom`/`check_magmom`），`--selftest` **19/19 PASS**（含与 `magnetic_orders.planar120` 逐位一致、与真实 `MnI2_pipeline/02_supercell/cases/*/magmom.txt` 自洽：n_mag=25 / n_atoms=75 / 120° 恰为 0·120·240 三色）。
16. `scripts\assemble_dft.py` 增 `--noncollinear` 统一入口 + **`LSORBIT` 开关**（`LNONCOLLINEAR+LSORBIT` + 3 分量 MAGMOM + `SAXIS`，禁止 `ISPIN=2`，需 `GGA_COMPAT=.FALSE.` 类配套）。→ **已落盘**（`assemble_dft.py` 重写：`build_incar(kind, …, magmom_line, noncollinear, lsorbit, saxis, npar)` 单一渲染入口；CLI 增 `--noncollinear/--lsorbit/--saxis/--all-orders/--selftest`；`--selftest` 18/18 PASS。`assemble_dft_tri.py` 变薄适配层（`--print-cmd` 可见它拼出的命令行），`--selftest` 10/10 PASS。驱动新增 `--noncollinear/--soc/--saxis`，日志与 `workflow.md` 都标出「非共线 / 非共线+SOC」）。
17. 非共线能量提取口径统一（`free energy TOTEN`）与收敛判据文档化。→ **已落盘**：新模块 `scripts\dft_energy.py`（`CONVENTION`/`read_energy_detail`/`read_toten`/`mixed_convention`/`check_same_convention`；一律取**最后一次** `free energy TOTEN`，缺了就报错而不是改用 `energy(sigma->0)`；`mixed_convention` 专抓混用口径，`MIX_TOL=1e-3` eV），`--selftest` 9/9 PASS。`solve_J3/solve_J_crect/solve_J_square/solve_J_tri/solve_J_rect/solve_cluster_J` 六个 solver 的 `read_toten` 全部改为转发它，`solve_fgt_J.py`/`fgt_verify_states.py` 的 OUTCAR 分支也改走它（保留 E.txt 兜底）。收敛判据（`EDIFF=1E-6`、`NELM=200`、非共线+SOC 需 `GGA_COMPAT=.FALSE.`、ISPIN 必须删掉）已写进 `WORKFLOW.md` §3 与 `skill\2d-magnet\SKILL.md`。
18. `tests\smoke_test.py` 非共线 INCAR 快照回归（不下发 DFT，仅校验 INCAR/MAGMOM 内容）。→ **已完成**：smoke §21（4 个新脚本 `--selftest` 全通过 + 驱动三角 `--soc` 走通 + `runs/<order>/INCAR` 的 `LNONCOLLINEAR/LSORBIT/GGA_COMPAT` 且**无 `ISPIN`** + MAGMOM 3 分量/面内/角度恰为 0·120·240 + 与快照 fixture `tests\fixtures\gdi2_nc_soc_120deg.INCAR.norm` 一致 + §115 硬报错正反两面 + **M3 mock-OUTCAR 闭环**：用 `coefficients_per_cell` 合成 `free energy TOTEN` 再让 `solve_J_tri.py` 回解，J1/J2/J3 对到 1e-4 meV）。全程不下发 DFT。

### P4 六角 kagome（依赖 P3）— 状态: **P4.1 系数层 ✅；P4.2 物理层：无 VASP 侧能做的已做完，三个 J 与 T_C 仍需真实机时**

> P4.1 / P4.2 的切分是**被迫的**，不是偷工减料：M4 的判据（三个 J 可解 + T_N 与文献量级一致）在本机**永远无法通过** —— 本机没有 VASP、从未跑过真实 DFT。所以把 P4 拆成两段：P4.1 交付**系数层里所有不需要 DFT 的部分**（几何判别、磁序枚举、系数矩阵、rank/cond 闸门、无 DFT 往返闭环、渲染、驱动接线、契约测试），P4.2 等真实机时（提交构型 → 读 OUTCAR → `J_result.json` → MC 求 T_C）。**本轮把 P4.2 侧「不需要 VASP 就能做完并验证」的部分全部做完了**（下面 ⑤–⑩），**物理层仍然没跑**（缺真实机时）。
>
> ① 偏差：第 21 项**没有**新建 `render_kagome.py`。kagome 的磁序图**复用**了 `render_tri_orders_pil.py` / `render_supercell.py`（两者都只吃 `cases/<order>/POSCAR + magmom.txt + magnetic_orders.json`，与具体晶格无关），driver 的 render 分支把 kagome 路由过去即可；新增一个渲染器等于多一份要同步维护的绘图代码，收益为零。**代价已处理**：三角渲染器的默认标题写死了「三角晶格」（对 kagome 是错的），因此给 `render_tri_orders.py` / `render_tri_orders_pil.py` 都加了 `--title` 覆盖开关，driver 的 kagome 分支传自己的标题 —— 复用不等于放任标题张冠李戴。
> ② 偏差：第 20 项的「三自由度分离」必须同时守住「**是哪三个 J**」的口径 —— kagome 的 J1 与 J1′ 是**同距离、不同亚晶格对**，本模型的 3 壳方案把它们并成一条 `J1`（取平均）。这一点写进了 `order_spec` 注释、`solve_J_kagome.py` 的 `notes` 与 `LIMITATIONS.md`，**不许当成「解出了 J1 与 J1′」**。
> ③ 偏差：判据比计划更**硬**。原计划只说「不可解即报错」，实现里加了三道闸门：builder 侧 `--supercell` 前两维必须是 3 的倍数、胞内必须恰好 3 个磁位、配位数必须等于 (4,4,6)；solver 侧 `rank < 4` 或 `cond > 1e4` 直接 `SystemExit`（对应 `ERROR_CATALOG["kagome_degenerate"]`）。
> ④ 新增（原计划没有）：`detect_lattice` 的 **kagome 判别**。Fe3Sn、Fe3GeTe2、Nb3I8 都是「3 个磁性原子/胞 + 最近邻 < cutoff 3.3 A」，旧版一律判 `cluster`，kagome 会被送进**分子簇模型**（物理全错）。判据 = **枚举周期像后的真实配位数**：kagome (4,4,4)（角共享三角形）、Fe3GeTe2 (4,6,4)（混合配位）、Nb3I8 (2,2,2)（三聚体）。只看胞内原子会得到 z=2，把 kagome 与三聚体混为一谈。
> ⑤ **M4 的锚点原措辞对 Fe3Sn 是错的（命名错误，已改）**：任务与旧文档把验收锚点写成「T_N 与文献量级一致」。**Fe3Sn 是铁磁体**（已核实的文献称其为 *kagome ferromagnet*，`P6_3/mmc`），它的磁有序温度是**居里温度 T_C**，不是奈尔温度 T_N ⇒ 锚点必须报 **T_C**。同一族的对照不许互相冒充：**FeSn**（`P6/mmm`）才是反铁磁（365 K 以下，层内仍铁磁耦合），**Fe3Sn₂**（`R-3m`）是**易轴**铁磁（T_C = 612 K，室温以下有自旋重取向 —— 二手来源里那个「Fe3Sn 低温自旋重取向」多半是它的污染）。
> ⑥ **文献锚点已核实并写入 `scripts\literature_J.json`**（原来该文件里**没有** Fe3Sn 条目 ⇒ 旧文档的「入口条件之一」现已从**未满足**变为**部分满足**）：来源 L. Prodan *et al.*, "Large ordered moment with strong easy-plane anisotropy and vortex-domain pattern in the kagome ferromagnet Fe3Sn", arXiv:2302.13810 (2023) / DOI 10.1063/5.0155295。写入的锚点 = 铁磁基态 + `P6_3/mmc` + **T_C = 705 K**（该文单晶实测；其引用的多晶文献值为 743 K，两个数都记进 `source`，不许混说）+ a = 5.4604 Å / c = 4.3458 Å + 饱和矩 2.27 μ_B/Fe + **易面各向异性，且易面就是 kagome 面**。**J1/J2/J3 仍然刻意保持 `null`**：该文自承 FexSny 的交换常数**尚未确立**，其第一性计算给的是饱和矩与各向异性，**不是**壳分辨的 J ⇒ 拿 Fe3Sn 条目当 `solve_J_kagome.py --ref` 会被按 null **明确拒绝**（这是设计行为，不是缺陷）。
> ⑦ **同族度量：2D 模型的临界温度与体相 T_C 差多少，本轮用本仓库自己的引擎量出来了**（不需要任何 DFT）。方法 = 拿**已核实文献 J** 的方形晶格反铁磁 La₂CuO₄（J1 = −111.8 meV, S = 1/2）跑有限尺寸 χ 峰扫描，与其实测 **3D** T_N = 325 K 比。实测（`mc_2d.py`，L = 12/18/24/30，`--seeds 2`，`--cool`）：**χ 峰 T = 290 / 260 / 250 / 240 K** ⇒ 比值 **0.89 / 0.80 / 0.77 / 0.74**，且**随尺寸单调下降**（二维各向同性 Heisenberg 由 Mermin–Wagner 定理无有限温长程序 ⇒ 这个峰是**有限尺寸伪临界点**，不是真相变）。与仓库既有记录自洽（仓库既有记录：本项目单层 χ 峰 266.67 K vs 325 K = 0.82）。
> ⑧ **⑦ 的直接后果（必须写进任何 P4.2 结论）**：kagome 的 2D MC 临界温度**不能**直接与 Fe3Sn 的 705 K 比 —— 本仓模型**没有层间耦合**（`mc_2d` 是二维格点模型），而 Fe3Sn 的 705 K 是**三维体相**有序温度，正是层间耦合撑起来的。同时注意**峰位随 L 单调移动**意味着「换个 `--mc-size` 就换个 T_C」，所以连「2D 模型自己的 T_C」都不是单一数。⇒ M4 的右半在**当前模型类里不可通过**（不是「暂时没机时」），要通过必须先在模型里显式加入层间 J_c（Fe3Sn 是双层堆垛，`c` 方向近邻就是层间），并用 U4 交点 + 尺寸外推给出收敛值。这条比「缺机时」更硬：**即使给了真实 J，现在的 MC 侧也给不出能与 705 K 对的 T_C。**
> ⑨ **回收通道已自动化**（原文档说「**没有**从 `03_dft/runs/` 的自动 ingest」，现已不成立）：新增 `scripts\ingest_dft.py`（`runs/` → `cases/` 的唯一收回通道，`--selftest` 18/18）与驱动 `ingest` 子命令。以前是手工 `cp`，而手工 `cp` 是整条链上**最危险**的一步（未收敛 / 混口径 / 错序 / 旧文件充数的 OUTCAR 会被静默当成好数据）。现在是**硬闸门**：四态齐全、体积下限、`General timing and accounting` 在场、SCF 收敛（`aborting loop because EDIFF is reached`）、**恰好一个离子步**（静态）、能量口径统一（转发 `dft_energy.mixed_convention`）、结构一致性（格矢 + 原子数跨态逐位比对）、四态能量差不退化、`magmom.txt` 与 OUTCAR 磁矩块对得上；**任何一条不过就一件都不写**（原子性），并打印「哪一条、为什么、怎么修」。每条闸门可 `--lenient <name>` 单独降级为警告（会记进 provenance）。产物带 `sha256` 与 provenance → `cases/ingest_manifest.json`。
> ⑩ **外部机时交接已自动化**（原文档只在「明确不做」里写「P4.2 不提交作业」，现已细化）：新增 `scripts\dft_handoff.py`（`pack` / `ingest`，`--selftest` 18/18）：`pack` 打出交接包（`README_HANDOFF.md` + `MANIFEST.json`（逐文件 sha256）+ `POTCAR.spec` + `CHECKLIST.md` + `submit.sh` + `jobscript.template` + `runs/<order>/` 的 INCAR/KPOINTS/POSCAR），`submit.sh` 自动识别 **SLURM / PBS / SGE / LSF / 裸机** 五种调度器并按元素顺序 `cat` POTCAR、跑四个序、跑完自查 OUTCAR 存在且够大、调 `ingest_dft.py --dry-run` 自检。**边界不变**：① **不打包 POTCAR**（license 不允许随包分发，只给元素顺序说明）；② **不代提交** —— 本脚本只生成脚本，提交动作由人在集群上做；③ `pack` 对 `runs/` 做完整性预检，缺文件直接拒（免得人跑到集群才发现少文件）。
> ⑪ **装配侧参数化 + k 点收敛抽查**（原文档把这几项写成「硬编码」/「未做收敛测试」，现已可显式覆盖并有工具）：`scripts\assemble_dft.py` 新增 `--kpoints a b c` / `--npar N` / `--isym N` —— 不给就**逐字节**保持旧默认（`3 3 1` + `NPAR=4` + static 不写 `ISYM`），`KPOINTS` 的过期表头（"for 2x2x1 monolayer supercell"）已改成 `Gamma-centered k-mesh a b c (vacuum along c)`；驱动 `dft` 同步透传，三角薄适配层 `assemble_dft_tri.py` 走同一份渲染（`--selftest` 12/12，`assemble_dft.py --selftest` 29/29）。新增 `scripts\k_convergence.py`（`run_2d_magnet.py kconv --outdir <proj> --k-grids 3,4,6,8`，`--selftest` 29/29）：把 `runs/<序>/` 复制成 `kconv/runs/<序>@<k网格>/`（**INCAR 逐字节不变，只换 KPOINTS** —— 换的是唯一变量），OUTCAR 放回后按体积下限/`General timing and accounting`/SCF 收敛三条闸门筛文件，算每个网格相对**最密**网格的 `ΔE` (meV/atom)、判 `converged` / `not_converged`，并给出 `recommended_kpoints`（`--k-tol` 默认 1.0 meV/atom）。**边界**：它判的是**总能量**对 k 网格的收敛、**不是 J 的收敛**（J 是能量**差**；差的收敛通常更快，但那是经验不是定理），**不代提交作业**，且本仓库没有任何真实 kagome OUTCAR ⇒ 这一节没有可引用的数字（见 `KAGOME_RUNBOOK.md` §5.3）。
19. `scripts\build_kagome_cases.py`（新增）：q=0 上/下三角组合 + √3×√3 120°（3 倍胞）+ 系数矩阵 + `--selftest`。→ **已落盘**：7 个磁序（`fm / q0_120 / q0_ud / q0_duu / sq3_120 / sq3_updntri / sq3_alt`），`selected = fm / q0_120 / sq3_updntri / sq3_alt`（**rank=4, cond=98.51, |det|=2.1e5**）；q=0 一族 (fm/q0_120/q0_ud/q0_duu) 的系数满足 `c1 == c2` ⇒ **rank 只有 2**（J1/J2 在 q=0 子块里不可分离，`cond≈2.3e34`），必须靠 √3×√3 序分离 —— 这条风险从「计划里的猜测」变成了**数值事实**。`--selftest` 覆盖配位数、壳距比 (√3, 2)、3x3 解析理想值、q=0 简并、设计矩阵、合成反解、规格层交叉校验、落盘产物。
20. `scripts\solve_J_kagome.py`（新增）：**3 壳模型 J1/J2/J3 的反解 + 条件数防护**（不可解即报错，不静默出 J）。⚠ **措辞修正**：本条原写「J1/J1′/Jd 三自由度分离」，而实现里 **J1 是 J1 与 J1′ 的同距离平均**、第三个壳是 J3 而不是 Jd ⇒ **不是**三自由度分离，交付物也不许这么写（见 `KAGOME_RUNBOOK.md` §2/§9）。→ **已落盘**：读 `cases/magnetic_orders.json` → 组装 `A` → `rank/cond` 闸门 → 读 `free energy TOTEN`（转发 `dft_energy.read_toten`）→ 写 `J_result.json`（`j_result.build/write`，`lattice="kagome"`）。无 DFT 自检入口是 `--demo --synthetic J1 J2 J3`（往返误差 < 1e-6 meV 才 PASS）—— 按仓库既有约定，builders 用 `--selftest`、solvers 用 `--demo`。
21. `scripts\render_kagome.py`（新增）。→ **改为复用**（见偏差 ①）：driver 把 kagome 的 render 路由到 `render_tri_orders_pil.py`（matplotlib 缺席时走 Pillow 兜底）与 `render_supercell.py`，两条路径都被 smoke §20/§23 钉住。
22. `run_2d_magnet.py`：kagome compute 分支 + 预设。→ **已落盘**：`--model fe3sn` 预设（`lattice="kagome"`、`S=1.0` 占位、`literature=None` —— 不谎称有文献锚点）、`detect_lattice` 判别分支、build/solve/render 三处分派、`ERROR_CATALOG["kagome_degenerate"]`、`--lattice` choices。
23. 技能文档 + smoke section + ≥20 用例。→ **已完成**：smoke §23（判别 + builder `--selftest` + solver `--demo` + 驱动 render 出图 + 落盘工程 + 契约测试）、`tests\test_kagome_cases.py` **68 项全 PASS**（两套独立实现逐位一致：2 个超胞 × 7 序；手推解析理想值；q=0 简并；√3×√3 标签良定义性；跨口径不变量；四类守卫报错；mock-OUTCAR 闭环；落盘复现；驱动集成）、`WORKFLOW.md` / `README.md` / `skill\2d-magnet\SKILL.md` 同批更新。

### P5 交变磁体（可选，半主链）

> **状态: ✅ 已完成（几何判别层 + 标量序参量；**不是**交变磁体的物理能力）**。按用户指示先做 P6 再回头做 P5，三项全部落地。
>
> 24. `scripts\mc_2d.py`：`--order sublattice:<mask>`（逐子晶格 Néel 矢量序参量，Binder U4 同构）。→ **已完成**，且做法比计划**更省**：不改模型 JSON、不加列，只加一个符号数组构造器 `sublattice_signs(num_basis, size, mask)`（`mask = {0}` 与老的 `stag_signs` 在 `num_basis >= 2` 时**逐位一致** ⇒ 蜂窝/FGT 数值不动，是天然回归锚点；`num_basis >= 3` 时能取任意子晶格组合，而老路径只会把 basis 1..n−1 并成一组）。产物仍是 4 列 `T_Kelvin,order_param,Susceptibility,Binder_U4` ⇒ `summarize_mc` / `mc_fss` 的列契约**零破坏**；`mc_fss.py` 同步支持（并显式拒绝 `cv`：U4 是它的判据）；driver 新增 `--mc-order` 透传（默认仍按 lattice 自动选）。掩码空 / 越界 / 覆盖全部子晶格（B 组为空 ⇒ `mbx/cb` 除零）全部硬报错。
> 25. 交变磁体**识别层检查清单**（零净磁化 + 子晶格反平行 + 晶体旋转连接）——只出报告，不下结论。→ **已完成**：新脚本 `scripts\altermagnet_check.py`（`--selftest` 15/15）：把「连接操作」定义成**必须同时是晶体对称操作**（把全部原子按元素映射回原位），再按旋转部分分类：平移 ⇒ 常规共线反铁磁、反演（非共面）⇒ 宇称禁止、旋转/镜面 ⇒ **几何候选**；报告永远带「不可判定：k 分辨自旋劈裂」一行。
> 26. 技能文档：明确「可算 = 子晶格 Néel 的 T_N；不可算 = k 分辨自旋劈裂验证」。→ **已完成**：`skill\2d-magnet\SKILL.md` 新增「分支：交变磁体（P5）——只做几何判别」节 + 局限追加一条；`LIMITATIONS.md` 10 + 脚本/测试清单；`WORKFLOW.md` 目录树 + 继承表 + 局限 `1g.`；smoke §25。
>
> ① 偏差（实施中发现并修掉两个**真问题**，都不在原计划里）：
> (a) **容差语义错了**：初版判据用「分数坐标四舍六入后查元组」，而真实结构文件里同一条对称关系两侧坐标常在第 6~7 位小数上不同 —— 实测 MnPS3 / CrOCl 的 fixture 因此被**静默误报**成「没有任何晶体对称操作连接两套亚晶格」（真实残差只有 0.012 / 0.503 mÅ）。改成按**Å 距离**判定（`nearest_same_element`），并在「无连接」时报告**最接近的候选操作与残差**，让人自己判断是「真的没有」还是「弛豫/写文件精度差一点」；`--tol` 因此以 Å 为单位。
> (b) **反演判定的优先级错了**：初版写成「只有反演连接（`kinds <= {inversion}`）才判宇称禁止」，于是 CrOCl（反演与真旋转/镜面在**同一 t** 下同时连接）被误判成几何候选。改成「**反演只要在场**就判宇称禁止」，并保留**共面 2D 胞的例外**（二维里「反演」与面内 C2 是同一操作，真二维交变磁体天然落旋转支）。
> ② 边界（必须一致地讲）：`altermagnet_check.py` 给的是**几何筛**，不是判据的终点 —— 能带层面的 k 分辨自旋劈裂**要能带计算**，本仓不做 ⇒ 输出只许出现「常规共线反铁磁 / 宇称禁止 / 几何候选（+ 为什么）」，**不许**出现「这是交变磁体」「劈裂多大」。序参量侧只给**标量 |L|** 的 U4 交点，**不分辨两套亚晶格的 Néel 矢量方向**。另一个必须讲清的反直觉点：**镜面连接也可能是常规反铁磁**（实测 MnPS3 = 蜂窝 Néel + 镜面连接 ⇒ 落进「候选」），所以「候选」不表示「是交变磁体」，更不表示「有劈裂」。

### P6 DMI / MTI（只做识别层 + 判据备忘）

> **状态: ✅ 已完成（识别层 + 备忘；按用户指示跳过 P5 先做 P6）**。交付物 = `DMI_MTI_MEMO.md`（22 KB：一页结论 + G1–G6 逐关复核 + 判别判据 + 拒绝模板 + 可行性备忘）+ `skill\2d-magnet\SKILL.md` 的判别节 + smoke §24。**没有**新增 DMI/MTI 计算子链 —— 原计划也没有，这是刻意边界。
>
> ① 偏差（多做了一件，理由充分）：§3.2 的「把文献 D 值写成 3×3 反称交换矩阵塞进 MC」从备忘里的**手推公式**升级为落盘工具 `scripts\dmi_matrix.py`（`--selftest` 5/5，含与 `mc_2d._energy` 的实机对拍）。因为这条等价式的符号极易搞反（`J_D = -[D]_×` 对应 `H_DM = +D·(S_i×S_j)`，反号则 DMI 帮倒忙），而它又决定 MC 结果；顺带修掉 `poscar_to_mc3.py` 读壳定义**不认 UTF-16 BOM** 的坑（PowerShell 的 `python ... > x.json` 会写成 UTF-16LE）。端到端实测：文献 D + 文献 J → `--shell-def` → `poscar_to_mc3` 匹配到 kagome 的 6 键/原胞，输出的 `J` 与 `-[D]_×·S²` 逐位一致、反称性非零。
> ② 边界（必须一致地讲）：P6 **不产生 D、不做拓扑判定**。能拿文献 D 跑「含 DMI 的格点模型」≠「本流水线算了 DMI」；MTI 侧同理（无层间 `J_c`、无能带拓扑）。§3.3 列了不许出现在结论里的输出清单。

27. 技能文档：写明两者为何不进主链、进主链所需的新子链形态。→ **已落盘**：`SKILL.md` 新增「DMI / MTI：怎么判别、能说什么」一节（判别线索、拒绝模板、`dmi_matrix.py` 一行命令、绝不输出的结论清单），并指向 `DMI_MTI_MEMO.md`；README 的「支持/不支持」与已知局限同步。
28. 可行性备忘录（若立项的前置与机时估算）。→ **已落盘**：`DMI_MTI_MEMO.md` —— G1–G6 逐关复核（为何 G1/G2/G4 仍过不了：没有 D 的第一性提取、没有外场/拓扑荷、没有层间分辨）、构型数量**可严格推出**、机时**纯估算**（假设 × 数量，未实测）、风险与替代路径、文献锚点（只用于判别，不代替计算）。

### P2′ 弱铁磁倾斜 + P4.2 手册（本轮补做，2026-10）

> **状态: ✅ 已完成（系数层 + 反解器 + 手册；**不是** CrOCl 的物理能力，也**不是** kagome 的物理层）**。
> 用户指示「把路线图剩余的条目一次做完」，可做的部分（不需要真实机时/DFT 的全部）在本轮落地，剩下两件**本质上需要真实 DFT**的事被写成**入口条件明确的交接物**。

**交付物**

1. **P2′ 系数层 + 反解器**：倾斜家族由 `scripts\build_crect_cases.py` 提供（`magnetic_orders.json` 的 cant 家族，装配走 `--soc`），反解器 `scripts\solve_J_crect_cant.py`（`--orders default|collinear|all|cant|逗号列表`、`--audit-only`、`--demo --ref/--synthetic`、`--json`、`--jfile-out`、`--no-jresult`、`--tol`）。契约测试 `tests\test_crect_cant.py`（**44 项**）+ smoke §26。**driver 已接通**（`build --cant-theta/--cant-base/--cant-azimuth` 透传，`solve` 见到 `cant_orders` 自动改派，三条守卫硬拒），smoke §26 用「驱动 build + 合成 OUTCAR + 驱动 solve」端到端钉住（J 回解误差 < 5e-4 meV）。
2. **P4.2 手册**：[`KAGOME_RUNBOOK.md`](KAGOME_RUNBOOK.md)（12 节：一页速览 / 开工前必须钉死的量 / 构型集与秩条件数 / 装配与五个坑 / SCF 四条判据与失败协议 / 回收 / 反解可信度 / T_N 与验收锚点 / 机时纯估算 / 与原计划的五项偏差 / 做不到的清单 / 检查清单）。

**① 偏差（原计划没有的东西）**

- **P2′ 的物理内核是一条闭式关系**，不是数值拟合：同一母序倾斜 θ 后的系数 = 共线母序系数 + `sin²θ ×` 与 θ 无关的基向量（`canted_coeff`）。这条关系让倾斜构型的系数**可解析验证**（已与测试里的独立实现逐位对拍），也决定了**θ 是输入参数**这件事 —— 原计划把 P2′ 写成「弱 FM 倾斜（可选）」，没说清「倾斜角从哪来」，本轮把它钉死：**基态 θ 要 SOC 总能量极小，本仓库不做**。
- **`--orders cant` 的语义与原计划的直觉相反**：它不是「只取倾斜态」，而是「倾斜家族 **+ 张成所需的 fm/母序行**」。原因：同一母序族的倾斜行彼此线性相关（家族自身 `rank=2`），只给倾斜态会得到一个不可解的方程组 —— 这正是原计划没预见到的自由度冗余（`--audit-only` 会把每行系数与自由度打出来）。
- **P4.2 的交付物是手册而不是代码**：原计划只写了「等真实机时」，没有入口条件、没有装配坑清单、没有 SCF 判据、没有回收路径说明、没有机时估算形式。本轮把它们补成可执行步骤（含**已知会坑人的五处**：`ISYM` 只在 `INCAR.relax`、`NPAR=4` 硬编码、`KPOINTS` 硬编码 `3 3 1` 且 header 过期、`POTCAR` 不生成、理想化胞未弛豫）。
- **补了 P4.2 的残余工作量分账**：原表只有一整条 `P4 7.5–12 人天`，切分后 P4.2 从未单独估过 ⇒ 本轮给 **≈4–9 人天**（不确定度高，不含机时等待），并回填 P6/P5/P2′ 的实测值（见第四节）。
- **driver 接通 P2′ 并顺手修掉一处脆弱的测试断言**：`run_2d_magnet.py` 新增 `--cant-theta/--cant-base/--cant-azimuth` 透传、`has_cant_orders(cases)` 判据、`solve` 的自动改派与三条守卫（`cant_theta_needs_crect` 进 error catalog ⇔ `KNOWN_ISSUES.md` 24 行）；同时 `tests\test_kagome_cases.py` 里「`--lattice` choices 含 kagome」原先靠「源码里第一处 `--lattice` 到下一个 `)` 之间的文本」判断 —— 只要别处（如 error catalog 的修复建议）提到 `--lattice` 就会误报，现改为解析 `run_2d_magnet.py --help` 的 choices 段并逐项核对全部 8 个取值。

**② 边界（必须一致地讲）**

- **P2′ 不产生 CrOCl 的 J/D**：θ 非自洽、一次倾斜 SCF 定不出唯一 (J, D)、SOC 必须开、文献 J 未核实 ⇒ 真值路径必须 `--ref`/`--synthetic`；验收只到「合成能量闭环 + 与独立实现逐位对拍」。
- **P4.2 手册不是进度**：路线图上 P4.2 仍未勾选，`literature_J.json` 里 Fe3Sn **条目已有、J1/J2/J3 仍为 `null`**（见 ⑥）；**入口条件 = 真实机时 + 已核实文献锚点 + 参数表填满**；手册 §8 的机时是纯估算，替换成实测前不许当预算。
- **J1/J1′ 分离、DMI 定量、任何「预测 T_N」的说法**都不在本轮范围内（见第五节的明确不做）。

## 三、工程化建议

1. **脚本对模式应当保留，但必须从「复制代码」升级为「数据驱动 + 薄适配层」。** 现状是 5 个 `build_*_cases.py`（14–31 KB 各一）反复重写 POSCAR 解析、cases 写出、selftest——每加一类都在复制同一批 bug。做法：晶格差异全部集中到 P0.1 的**磁序规格数据**；`build_*` 只声明规格，`solve_*` 只做「读规格 → 组装 A → 求解 → 条件数校验 → 写 J_result.json」。
2. **统一接口冻结为 v2 契约**：`build_*` 入参 `--poscar/--element/--supercell/--outdir/--selftest`，出 `cases/magnetic_orders.json`（含 `order_spec[]`）；`solve_*` 入参 `--cases/--demo/--ref`，出 `J_result.json`；MC 侧 `poscar_to_mc3 --jfile J_result.json`。**每个脚本必须自带无 DFT 的 `--selftest`**（合成能量精确回解，沿用现有 `solve_J3 --demo` 模式）。
3. **单点关键改动优先做**：`poscar_to_mc3` 的 `--jfile/--shell-def/多元素` 是解锁方形、kagome、oxyhalide 三类的共同前提，先做它，三类才可能并行。
4. **验证四级阶梯**：
   - L1 合成能量精确回解（`--selftest`，无 DFT）；
   - L2 双路径交叉：文献 J → MC 的 T_C 与 DFT J → MC 的 T_C 对照（走现有 `literature_J.json` 通道）；
   - L3 端到端 mock：smoke_test 同构 section + 每类一个 mock DFT；
   - L4 网格无关：`mc_fss.py` Binder U4 交点定 T_C/T_N（避免重演 MnI2 的窗宽跨峰缺陷）。
5. **回归门**：契约冻结（P0）后，任何新类型开工的前置是 smoke 全绿 + 现有 5 条链路在 `--jfile` 路径下数值与旧路径一致；每类须有 ≥20 个磁构型回归用例、记录基线分数、上线前跑全流水线回归。
6. **不静默降级**（对应「绝不无声截断」规则）：系数矩阵条件数超阈或残差超阈 ⇒ 停下上报，不得给出 J；需要 SOC 而 INCAR 未开 SOC ⇒ 直接报错，不得按共线近似继续；文献兜底 J 只在数值对不上并且显式声明 `source` 时使用。
7. **产物与命名约定不破坏**：新类型产物继续落 `2d_magnet_pipeline\out_<type>_<date>\`，新增子链不得改动既有目录语义。
8. **文档与代码同步是交付条件**：`SKILL.md` 路由表、`references\magnet-types.md` 标注、`WORKFLOW.md` 阶段产物必须与代码能力同批更新，「可算」标签只在端到端跑通后加。

## 四、工期估算（人天；1 人 = 1 人天）

| 阶段 | 内容 | 人天 | 依赖 | 不确定度 |
|---|---|---|---|---|
| P0 | 地基：规格层 + 契约冻结 + `--jfile` + 契约测试 | 12–17 | — | 中（±3） |
| P1 | 方形 K2NiF4 全链 | 5.5–9 | P0 | 低 | ✅ 已完成 |
| **P2** | 矩形/oxyhalide（不含弱倾斜） | ~~2–3~~ → **实测 ≈8–10** | P0（P2 可与 P1 并行） | ~~低~~ → **中**（计划前提被推翻，见 §P2） | ✅ 已完成 |
| P2′ | oxyhalide 弱 FM 倾斜（可选） | 3–5 | P3 | 高 | **✅ 已落地（实测 ≈1 人天 —— 原估偏悲观 3–5×）**：系数层 + 反解器 `scripts\solve_J_crect_cant.py` + 契约测试 `tests\test_crect_cant.py`（44 项）+ smoke §26。**边界**：θ 是**输入参数**（闭式 `sin²θ` 标度可解析验证）、一次倾斜 SCF **定不出唯一 (J, D)**、装配必须 `--soc` ⇒ 不产生 CrOCl 的 J/D 结论 |
| P3 | 非共线 + SOC 基础设施 | 5–8 | P0 | **最高（±5）** |
| P4 | 六角 kagome 全链 | 7.5–12 | P3 | 高 | **P4.1 ✅（纯几何 + 线性代数）；P4.2 手册 ✅**（[`KAGOME_RUNBOOK.md`](KAGOME_RUNBOOK.md)）——**残余工作量补估 ≈4–9 人天**（不含机时等待，见该手册 §9），物理层待真实机时 |
| P5 | 交变磁体序参量扩展（可选） | 3.5–5.5 | 无（可与 P3/P4 并行） | 中 |
| P6 | DMI/MTI 判据文档 + 可行性备忘 | ~~2–2.5~~ → **实测 ≈0.5** | 无（可全程并行） | 低 | ✅ 已完成（另落盘 `scripts/dmi_matrix.py` + smoke §24） |

- **主链（P0+P1+P2+P3+P4）= 32–49 人天 ≈ 6.5–10 周（1 人）**；2 人合理分工约 **3.5–5.5 周**。
- 含 P2′+P5 的完整版 = 38.5–59.5 人天；P6 文档 +2–2.5 人天。**回填实测（2026-10）**：P6 ≈0.5 人天（原估 2–2.5 ⇒ 偏悲观 4–5×）、P5 ≈1.5–2 人天（原估 3.5–5.5）、**P2′ ≈1 人天**（原估 3–5 ⇒ 偏悲观 3–5×）、**P4.2 手册 ≈0.5 人天 + 物理层残余 ≈4–9 人天**（原 P4 行 7.5–12 从未切分过 ⇒ 这是本轮补的分账）。共同规律：**凡是「只做识别层 / 系数层 / 文档」的部分都严重高估**（2–5×），**凡是要真实 DFT 的都低估**（低估幅度无法估计，因为没有一次真跑过）。

不确定项与风险：
1. **P3 非共线 SOC 无既有工作流**——全项目方差最大项。缓解：先用一个已知体系做最小闭环（单体系、单构型、SOC 开/关能量差量级正确）再推 kagome。
2. **DFT 机时（不是人天，是真实瓶颈）**：非共线 + U + SOC 的 3 倍胞单点能量约为现有共线单点的 3–10 倍；构型数**以 rank/cond 为准**（本仓库 builder 的 `selected_orders` 是 **4** 个；原写的「每类 ≥6 构型」与本仓库实现不一致，已按实测口径修正），每类 ≥2 超胞。建议单列机时预算，并规定「每类最少 3 个构型」的下限，否则 G2 不可解。
3. **kagome 的 J1 与 J1′ 在 q=0 构型下简并**，必须靠 √3×√3 分离；若该构型非共线 DFT 收敛困难，则 J 分解失败。缓解：P0.1 先用几何算出系数矩阵条件数，再决定提交哪些构型（不许先算 DFT 再发现不可解）。**实测补充**：q=0 一族的简并已被数值证实（子块 `rank=2`、`cond≈2.3e34`），含 √3×√3 的 `selected` 组合给出 `rank=4`、`cond=98.51`；而真正的 **J1/J1′ 分离本仓库做不到**（3 壳模型把它们合并成一条 J1）——原文把它与「q=0 简并」混成一件事，实际是两件（见 `KAGOME_RUNBOOK.md` §2）。
4. **文献锚点缺失**：部分 oxyhalide 成员无可靠 T_C，eval 只能对标 J。须在 `J_result.json` 标注 `source` 与置信度，禁止拿 J 对标冒充 T_C 校验。
5. **元素/参数预设污染**：smoke_test 13 节已把「`MODEL_PRESETS` 必须与 `literature_J.json` 材料专属值一致、不得与 `ELEMENT_PARAMS` 通用兜底混淆」钉死，新类型必须同构遵守。
6. **文档能力虚高**：6 类已在文档里是「识别 + 路由」，一旦某类进主链而文档未同批更新，会出现「文档宣称能力 > 实际能力」，这是本技能最容易被误用的故障模式。缓解：文档更新并入该类的完成定义（DoD）。

## 五、实施顺序与可并行部分

硬依赖：`P0.1 → P0.2 → P0.3 →` 之后 P1 / P2 可并行；**P3 是 P4 与 P2′ 的共同前置**；P0.3 是任何新类型在 MC 侧的前置。

里程碑门控：
- **M1 = P0 完成**：契约冻结 + smoke 全绿 + 5 条既有链路在 `--jfile` 路径下数值与旧路径一致 ⇒ 才允许新类型开工。
- **M2 = P1 完成**：方形类端到端（build → solve → MC → T_C）跑通并与文献锚点对齐 ⇒ 地基得到实证（这一步是「模式是否复用成功」的判据）。**已达成**（χ 峰 266.67 K vs La₂CuO₄ T_N 325 K，比值 0.82）。
- **M2b = P2 完成**：C 心矩形（CrOCl）端到端跑通（判别 → 序枚举 → 反解 → MC → 跨口径校验）。**已达成但降级**：CrOCl 的 **J 与 T_N 均未核实**，所以端到端只用 `--synthetic` 合成 J 验证管路，产物不入库、结论不声称是 CrOCl 的物理结果；这一条同时把「计划的前提可能被实测推翻」写进流程（P2 就是先例：原以为只是复用 rect 对，实测发现是 C 心网）。
- **M3 = P3 完成**：一个非共线 SOC 最小闭环成功 ⇒ 解锁 kagome。**已达成但降级**：闭环里的「DFT」环节是 **mock OUTCAR**（用规格层的 `coefficients_per_cell` 合成 `free energy TOTEN`，再让 `solve_J_tri.py` 回解，J 对到 1e-4 meV；见 smoke §21），**没有任何真实非共线 DFT 能量**。所以 P3 解锁的是**能力**（INCAR/MAGMOM 写对、SOC 门禁硬报错、能量口径统一），**不是**「SOC 物理已验证」：P4 的 kagome 构型现在可以设计并提交，但其 J 与 T_N 仍必须等真实机时结果，禁止把 P3 的绿灯当成物理结论。
- **M4 = P4 完成**：kagome 三个 J 可解（条件数可控）+ T_N 与文献量级一致 ⇒ 进主链。**未达成（P4.1 只交付了左半）**：系数层已经可解（`rank=4` / `cond=98.51` / 无 DFT 往返闭环 / q=0 简并已被数值证实），但 **J 与 T_N 都还没有** —— 本机没有 VASP、从未跑过真实非共线 DFT，Fe3Sn 的文献 J 也不在 `literature_J.json` 里（预设刻意写 `literature=None`）。所以 kagome 现在是「**能力就绪、物理未验证**」，与 M3 同款降级：P4.2（提交 √3×√3 构型 → 读 OUTCAR → MC 求 T_N）必须在真实机时上做，禁止把 P4.1 的绿灯当物理结论。**本轮补上了 P4.2 的入口条件与操作手册**（[`KAGOME_RUNBOOK.md`](KAGOME_RUNBOOK.md)）：入口 = 真实机时 + **已核实**的文献锚点（`literature_J.json` 里 Fe3Sn 至今为空 ⇒ 这是**硬条件**，不是提醒）+ 参数表填满；手册 §8 的机时是**纯估算**，在用一个 FM 参考点替换成实测值之前**不许当预算**。
- **M5**：文档路由表与代码能力对齐，「可算」标签只在端到端通过后添加。

并行安排：
- **双人方案**：A 做 P0 + P1（地基 + 方形，串行链在 A 手上）；B 在 M1 后接 P2（oxyhalide 参数与几何）与 P6 文档，随后与 A 一起进 P3。（**单线程实际执行顺序**：P1 ✅ → P2 ✅ → P3。）
- **最大并行机会**：P4.1 的 kagome 磁序枚举与系数矩阵、条件数分析是**纯几何 + 线性代数，不需要 DFT**，可与 P3（SOC 工作流）完全并行——SOC 一旦跑通即可直接提交已设计好的构型。
- **P5（mc_2d 序参量扩展）与 P4 无耦合**，可在 P3 期间由第三人（或 B）并行推进。
- **P6 纯文档**，任何阶段都可插入。

明确不做（写出边界，避免范围蔓延）：DMI 与 MTI 不进主链（只留判据与可行性备忘）；交变磁体不做 k 分辨自旋劈裂验证（研究问题，不属本流水线契约）；**kagome 的 J1/J1′ 分离不做**（3 壳模型合并，需新构型 + SOC）；**弱铁磁倾斜的基态 θ 与 D 不做自洽求解**（θ 是输入参数）；**P4.2 不提交作业、不写 POTCAR、不提供机时**。
