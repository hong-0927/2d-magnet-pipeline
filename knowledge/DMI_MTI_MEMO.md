# DMI / MTI：判据备忘 + 可行性备忘（P6）

> **状态：只做识别层 + 判据备忘。** 本文件是**文档**，不含新的物理计算能力：P6 没有给流水线增加任何
> 「算 D 矢量」或「算拓扑不变量」的代码。`LATTICE_ROADMAP.md` 的 P6 只有两条（第 27 项技能文档、
> 第 28 项可行性备忘），这里把两条都落成可查的判据与数字。
>
> **一句话结论**：体相 DMI 与磁性拓扑绝缘体（MTI）**都不进主链**，理由是**契约对不上**（G6 ❌）而不是
> 「算不动」——主链的终点是 `J` → `T_C/T_N` 的二维海森堡链；DMI 的产物是 3 分量 D 矢量 + 拓扑荷，
> MTI 的产物是层间 `J_c` + 能带拓扑，两者都要**另立子链**。用了它们而套现有 CLI，会得到「看着正常、
> 物理不对」的数（本仓最忌讳的静默降级）。
>
> **数字纪律**：本文件里凡出现机时/成本，一律标注是**估算**并写出假设（§5.3）。**没有任何一个机时数字
> 是实测**——本机没有 VASP、没有任何 OUTCAR，DFT 从未跑过。

---

## 0. 一页结论

| | 体相 DMI（Dzyaloshinskii–Moriya） | 磁性拓扑绝缘体（MTI，MnBi₂Te₄ 类） |
|---|---|---|
| 判定（`LATTICE_ROADMAP.md` §1.2 原文） | **只识别 + 路由，不进主链** | **只识别 + 路由，不进主链** |
| 卡在哪一关 | G2 ❌（3 分量 D 拟合成本高）、G3 ❌（无自旋螺旋提取工作流）、G4 ❌（无拓扑荷序参量）、G6 ❌（产物对不上契约） | G2 ❌（需层分辨 A 型 AFM）、G3 ❌（SOC + slab + 层间）、G4 ❌（单层 MC 无意义）、G6 ❌ |
| P3/P4.1 之后有没有变化 | **G3 部分解除**：`LNONCOLLINEAR`+`LSORBIT` 的 INCAR 层已就位、螺旋/倾斜自旋矢量已就位；仍缺**能量提取与 D 拟合** | **无变化**：层分辨/多 QL 与能带拓扑仍全缺 |
| 今天能做的 | ① 算对称部分 `J`（若磁性格子在主链分支内）② 用 `--noncollinear --soc` 算**倾斜/螺旋构型的能量** ③ **把文献 D 值写成反称 3×3 `J_mat` 塞进 MC**（§3.2；工具已落盘 `scripts/dmi_matrix.py`，`--selftest` 5/5，端到端实测） | ① 把**单层/单个 QL** 当二维层算层内 `J`（但**层内 T_C ≠ 材料的 T_N**）② 用 `mp_fetch_2d.py` 判断层状性与层数 |
| 明确不能做的 | 算 `D`；算 skyrmion 相图/拓扑荷；报「有 DMI」的定量结论 | 算层间 `J_c`；算三维 `T_N`；算 Z₂/能带拓扑/轴子态 |
| 进主链缺什么 | 四态法（或自旋螺旋）提取 + D 矢量拟合 + 拓扑荷序参量 + `D_result.json` | 多 QL slab 构建 + 层分辨磁矩 + `J_c` 反解 + **层状（3D）MC** |

---

## 1. 为什么不在主链：G1–G6 逐关复核（P3 / P4.1 之后）

判定门原文（`LATTICE_ROADMAP.md` §1.1）：**G1–G4 全过 + G5 有锚点 + G6 可对位 = 进主链**；
**G1–G3 任一需要新子链 = 只识别 + 路由**。

### 1.1 体相 DMI

| 关 | 判定 | P3/P4.1 之后的事实（行号可查） |
|---|---|---|
| G1 可编码 | ✅ | 3×3 交换张量在 MC 模型层已支持：`scripts/poscar_to_mc3.py` 的模型 JSON 逐壳 `J_mat`（3×3，文档行 8/30-32，解析 337-348），`scripts/mc_2d.py:44-62` 用 `J_array[(N, maxn, 3, 3)]` 并给反向键存 `J.T` ⇒ **反称部分（= DMI）数学上可承载** |
| G2 参数可拟合 | ❌ | 仓库**没有**任何 D 矢量拟合：`solve_*` 反解的是标量 `J1..J3`；没有四态法/螺旋拟合脚本 |
| G3 自旋构型可算 | ⚠→**部分解除** | INCAR 层已就位（`assemble_dft.py --noncollinear --lsorbit --saxis`）；`scripts/mag_vectors.py` 的 `SOC_ORDERS` 含 `cant*/spiral/tilt`；**仍缺**：四态法所需的「每键 4 个取向」构型生成器、螺旋 `q` 扫描构型生成器、以及把能量差转成 `D` 的提取器 |
| G4 序参量可得 | ❌ | `mc_2d.py` 的序参量只有交错磁化（`stag_signs`）与比热；**没有**拓扑荷 `Q`、skyrmion 半径、螺旋 q；**也没有外磁场**（全仓无 `field`/`h_ext`）⇒ 无法做 B–T 相图 |
| G5 有锚点 | ⚠ | `literature_J.json` 里没有 `D` 字段；文献 D 值只能人工录入 |
| G6 契约可对位 | ❌ | 主链终点是 `J_result.json` → T_C/T_N；D 矢量的天然产物是「D_x,D_y,D_z per bond + 手性 + 拓扑荷」，塞不进 `j_result/1` schema（`scripts/j_result.py`） |

**结论**：契约（G6）先否决，其次是 G2/G4。**不是因为算不动，是因为算了也没地方放、没法验收。**

### 1.2 MTI

| 关 | 判定 | 事实 |
|---|---|---|
| G1 可编码 | ⚠ | 层内 FM 可以编码；**层间 AFM 需要多 QL 胞**，而流水线是「单层 + c 轴真空」假设（`reorient.py` / `assemble_dft.py` 的 `monolayer` 语义） |
| G2 参数可拟合 | ❌ | 需层分辨 A 型 AFM 的逐层磁矩与层间 `J_c`；现有 OUTCAR 解析只有总能量（`dft_energy.py`）与**每离子磁矩**（`fgt_verify_states.py:28/54` 解析最后一个 `magnetization (x)` 块）——**逐层**聚合与 `J_c` 反解不存在 |
| G3 自旋构型可算 | ❌ | SOC + slab + 层间：即使 `--soc` 已能写 INCAR，多 QL 胞（100–250 原子）与高 k 网格的成本属于另一个量级（§5.3） |
| G4 序参量可得 | ❌ | **层状 MC 不存在**：`mc_2d.build_arrays(model, size)` 只吃 `(Lx, Ly)`（`mc_2d.py:33-34`），`poscar_to_mc3.py:501` 写死 `dim=2` ⇒ 无法表达「层内 J + 层间 J_c」的三维格点 |
| G5 有锚点 | ⚠ | MnBi₂Te₄ 的 T_N ≈ 25 K 是文献锚点，但**锚点对照的是层间耦合控制的三维有序温度**，与单层 MC 的输出不是一个量 |
| G6 契约可对位 | ❌ | 拓扑性质（Z₂、表面 Dirac 质量、轴子态）需要能带/波函数工具链，本仓没有任何能带解析（连 DOS 都没有） |

**结论**：MTI 的死穴是 G4+G6——「单层 MC 无意义」这句话在代码层是可验证的事实（`dim=2` + 只有二维尺寸参数）。

### 1.3 与计划原文的差异（要记住的两条）

1. **DMI 的 G3 状态变了**：计划写在 P3 之前（「无 LSORBIT / 无自旋螺旋工作流」）。P3 把 `LNONCOLLINEAR` +
   `LSORBIT` + `SAXIS` + 螺旋/倾斜自旋矢量都做出来了，所以 G3 从 ❌ 变成 ⚠（**输入侧已通，提取侧未做**）。
   **但判定不变**：G2/G4/G6 仍是 ❌ ⇒ 依旧不进主链。
2. **`references\magnet-types.md` 不存在**：计划 §P1 第 11 项与 §三.8 都提到这个文件，仓库里
   **从来没有 `2d_magnet_pipeline/references/`**；「可算 / 不可算」的标注实际落在 `skill/2d-magnet/SKILL.md`
   （agent 入口）与 README/WORKFLOW。本备忘录沿用这个事实，不再指向不存在的路径。

---

## 2. 识别判据（给 agent 用：先判断，再决定要不要拒绝）

### 2.1 体相 DMI 型

**正面判据（满足越多越像；①②是硬条件）**

1. **非中心对称**：晶体/界面缺反演中心。一阶 DMI 要求 SOC **且**无反演对称 —— 中心对称的完美体相结构
   一阶 `D = 0`（仍可能有各向异性交换，但那不是 DMI）。
2. **强 SOC 元素参与磁性或配位**：4d/5d（Ru/Rh/Ir/Pt）、5f、重 p 配体（I、Te、Bi、Se）、
   或 3d 磁性离子 + 重配体（CrI₃ 型、CrTe₂ 型）。
3. **磁序是螺旋/摆线/弱铁磁倾斜**：而不是简单共线 FM/AFM。实验线索：小角中子散射的 `q ≠ 0` 卫星峰、
   拓扑霍尔信号、skyrmion 相、异常大的矫顽力。
4. **磁各向异性与 D 同量级**（MAE 不可忽略）——这决定了 skyrmion 是否稳定。
5. 磁性格子在主链支持的六类之内（否则连 `J` 都算不了，问题不在 DMI）。

**易误判**

- 「有重元素」≠ 有 DMI（中心对称 + 重元素 ⇒ 一阶 D = 0）。
- La₂CuO₄ / K₂NiF₄ 这类中心对称 3d 氧化物：**不要**往 DMI 方向解释，它的倾斜是结构畸变/Dzyaloshinskii 之外的机制。
- CrOCl 的弱铁磁倾斜属 **P2′**（P3 已给 SOC 工具），**不是** DMI。

### 2.2 MTI 型

**正面判据**

1. **层状 vdW**，且**层内 FM、层间反平行（A 型反铁磁）** ⇒ 净磁化 ≈ 0 但每层有矩。
   `scripts/mp_fetch_2d.py` 的 `split_layers` / `layered_metrics`（行 109-206）已能判层状与层数，
   `--list` 里 `layers` 与 `total_magnetization=0` 是**线索**（不充分）。
2. **化学式落在 MnBi₂Te₄ 族**：MnBi₂Te₄ / MnSb₂Te₄ / MnBi₄Te₇ / MnBi₆Te₁₀（Bi/Sb/Te 重 p 元素）。
3. **T_N 很低（~20–30 K）且随层数/堆垛显著变化** ⇒ 由层间 `J_c` 控制（MnBi₂Te₄ 约 24–25 K）。
4. 奇数层剥离后呈弱 FM、偶数层呈 AFM（层数奇偶决定净矩）。
5. 输运/ARPES 线索：Z₂ 非平凡、量子反常霍尔、轴子绝缘体态、表面 Dirac 锥有交换带隙。

**易误判**

- Fe₃GeTe₂：层状 + 重元素 Te，但是**铁磁金属**，不属 MTI（它是本仓已支持的 FGT 双亚晶格分支）。
- CrI₃ / CrSBr：层状 FM 绝缘体/半导体，**Hubbard 型**而非拓扑绝缘体。
- Nb₃I₈：簇磁体，与 MTI 无关。
- **关键**：「层状 + 重元素」**不等于** MTI。判断拓扑非平凡需要能带/波函数工具 —— **本流水线无法判定**，
  这也是为什么 MTI 只能「识别 + 拒绝」，不能「识别 + 折算」。

### 2.3 拒绝模板（结论里必须这么写）

```
该材料属于 <体相 DMI / MTI> 类：<一句判据>。
本流水线对它的输出只到 <可做的部分>，不能给 <D 矢量 / 拓扑荷 / 层间 J_c / 三维 T_N>，
原因：<G6 契约对不上 / 需要层分辨三维模型 / 需要能带工具>。
若要推进，见 DMI_MTI_MEMO.md §4（新子链形态）与 §5（前置与机时估算）。
```

---

## 3. 今天能做什么 / 不能做什么

### 3.1 已就位、可复用的资产（**均已实测**，行号可查）

| 资产 | 位置 | 对这两类材料有什么用 |
|---|---|---|
| SOC 非共线 INCAR 层 | `scripts/assemble_dft.py`（`--noncollinear --lsorbit --saxis`）、`assemble_dft_tri.py`（薄适配层） | 写对 DMI/层间计算所需的 VASP 输入；需要 SOC 的磁序没开 SOC ⇒ 硬报错 |
| 自旋矢量层 | `scripts/mag_vectors.py`（`SOC_ORDERS` 含 `cant*/spiral/tilt`、`spin_cant` 的净矩 = 2 sinθ） | 倾斜/螺旋构型的 3 分量 MAGMOM；DMI 的物理图像（弱倾斜） |
| 能量口径 | `scripts/dft_energy.py`（只认 `free energy TOTEN` 最后一次出现；`mixed_convention` 抓混用） | D 的提取要靠**能量差**，口径统一是前提 |
| J 产物 schema | `scripts/j_result.py`（`j_result/1`、residual/cond 闸门） | 未来 `D_result.json` 的样板（**但不要塞进 `j_result/1`**） |
| 3×3 交换张量（MC 侧） | `scripts/poscar_to_mc3.py`（逐壳 `J_mat`，行 8/30-32/337-348）、`scripts/mc_2d.py:44-62`（`J.T` 反向键）、**`scripts/dmi_matrix.py`（新增：D → `J_mat` 换算器，`--selftest` 5/5）** | **今天就能把文献 D 值写成反称矩阵跑 MC** → §3.2 |
| 单离子各向异性 | `poscar_to_mc3.py --MAE` → `onsite[:,2,2] = -MAE`（行 471-472）；`mc_2d.py:37-39` 读 `onsite` | skyrmion 稳定性需要 MAE 与 D 的竞争 |
| 每离子磁矩解析 | `scripts/fgt_verify_states.py:28/54`（最后一个 `magnetization (x)` 块） | 层分辨磁矩的现成样板（MTI 的第一步） |
| 层状性/层数判断 | `scripts/mp_fetch_2d.py:109-206`（`split_layers` / `layered_metrics`） | MTI 识别线索 + 剥离层数 |
| 非共线快照/闭环回归 | `tests/smoke_test.py` §21；`tests/test_kagome_cases.py`（68 项） | 新子链的测试样板 |

### 3.2 今天就能做的一件事：把**文献 D 值**塞进 MC 跑（已验证的等价形式）

DMI 项 `H_DM = D · (S_i × S_j)` 可以**精确**写成 3×3 交换矩阵形式 `S_iᵀ · J_D · S_j`：

```
J_D = -[D]_×        （[D]_× v = D × v，即 [D]_× = [[0,-Dz,Dy],[Dz,0,-Dx],[-Dy,Dx,0]]）
```

符号由来：`S_iᵀ(-[D]_×)S_j = -(D × S_j)·S_i = -(S_j × S_i)·D = +D·(S_i × S_j)` ✓
（若你的约定是 `H_DM = -D·(S_i×S_j)`，则 `J_D = +[D]_×` —— **符号必须与你的 J 符号约定一起写清楚**，
本仓约定见 `SKILL.md`「J 符号」：`H = -J Σ S·S`，`J>0 = 铁磁`。）

于是：`poscar_to_mc3.py --shell-def` 的每个壳给 `J_mat = J_iso·I₃·(-S²) + J_D`，`mc_2d.py` 的
`J_array` 本来就是 3×3、反向键存 `J.T`，所以**反称部分会被正确消费**（`S_iᵀ J S_j` 对
`J = J_sym + J_anti` 自动分出对称/反对称两部分）。

**但这条路径有三个必须一起声明的限制**：
1. **D 是外来的**（文献值或别人的计算），本仓不产生 D —— 所以 §1.1 的 G2 依旧 ❌。
2. **没有外磁场、没有拓扑荷**：`mc_2d.py` 无 `field`，序参量只有交错磁化与比热 ⇒ 只能看**零场基态**
   （螺旋/摆线/铁磁谁赢）与能量，**做不出 skyrmion 的 B–T 相图**，也报不出 `Q`。
3. 标量 `J` 字段会变成 `NaN`（`poscar_to_mc3.py:481` 的设计如此：含 DMI/各向异性时本来就没有单一标量）。
   下游凡是读标量 `J` 的地方（如 `mc_2d.py:375` 的 `jinfo` 打印）要能容忍。

**工具（已落盘）**：`scripts/dmi_matrix.py` 把上面的换算、S² 吸收和 `--shell-def` JSON 片段一次做完，
`--selftest` 5/5 通过（含与 `mc_2d._energy` 的实机对拍）：

```bat
python scripts\dmi_matrix.py --D 0.31 -0.77 0.55 --J 1.9 --S 1.5 --d 2.70 --tol 0.05 --out shell.json
python scripts\poscar_to_mc3.py <POSCAR> --element Fe --S 1.5 --shell-def shell.json -o model.json
```

实测（Fe3Sn kagome fixture、`D=(0.31,-0.77,0.55)`、`J=1.9`、`S=1.0`）：`poscar_to_mc3` 匹配到
**6 键/原胞**（kagome J1 的 `z=4` ⇒ 每胞 6 键 ✓），输出的 `J = [[-1.9, +0.55, +0.77], ...]` 与
`J_D = -[D]_×·S²` 逐位一致，反称性 `J[0][1]+J[1][0] == J[0][2]+J[2][0] == 0` ✓。

两个实测踩过的坑：
1. **别用 shell 的 `>` 重定向**：PowerShell 里 `python ... > x.json` 会写成 **UTF-16LE + BOM**，
   而 `load_shell_defs` 原本只按 `utf-8-sig` 读 ⇒ `UnicodeDecodeError: 'utf-8' codec can't decode
   byte 0xff in position 0`。现已两边都修：`dmi_matrix.py --out` 自己写 UTF-8；
   `poscar_to_mc3.py` 按前两字节显式分派 `utf-16` / `utf-8-sig`（并在来源行打出 `[utf-16]`）。
2. **S² 只吸收一次**：`dmi_matrix` 已经把 `S²` 乘进 `J_mat`，所以给 `poscar_to_mc3` 的 `--S` 必须与
   建矩阵时一致。`J_mat` 路径**不会**再乘 `s2`（`poscar_to_mc3.py:336-348`），`--no-s2` 只影响标量 J
   的反推与输出里的 `S2` 记录 —— 两边不一致不会报错，只会静默差一个 `S²` 因子。

> 结论：这条路能让「含 DMI 的格点模型」跑起来，**但它不是 DMI 子链**——它是「拿别人的 D 做格点模拟」。
> 写结论时必须这么表述，不许说成「本流水线算了 DMI」。
>
> **契约（改坏了会红）**：`tests/smoke_test.py` §24 钉住了四个不变量 —— `dmi_matrix --selftest` 全绿、
> 换算出的反称矩阵**逐位**落进 mc3 的 6 键/原胞、UTF-16 壳定义可读、`ERROR_CATALOG` 有
> `dmi_needs_shell_def`（点名「`--jfile` 会静默丢掉反称项」这条坑）；另外 README/SKILL.md 必须提到
> 本备忘与 `dmi_matrix.py`，否则 §24 也会红。

### 3.3 明确不支持的输出（不许出现在结论里）

- `D_x, D_y, D_z`（本仓不提取）、`|D|/J` 比值、手性符号的**第一性**判断；
- skyrmion 半径、拓扑荷 `Q`、B–T 相图、临界 D/J；
- MTI 的层间 `J_c`、层分辨磁矩表、三维 `T_N`、Z₂ / 表面态 / 轴子态；
- 任何形式「本材料是 DMI 磁体 / 拓扑绝缘体」的**判定结论**（前者需要 D 的第一性计算，后者需要能带拓扑）。

---

## 4. 若立项：需要的新子链形态

### 4.1 DMI 子链（≈5 件，不改既有产物语义）

1. **`scripts/build_dmi_cases.py`（新增）**：按**四态法**给每个不等价键生成 4 个自旋取向构型
   （两自旋各自沿键垂直平面内的正交方向），外加 FM/基态参考；每个构型一个目录（`POSCAR` + 3 分量
   `magmom.txt`），落 `cases/dmi_configs.json`（记录每构型的**键 → 取向映射**、`SAXIS`、约定）。
   —— 复用 `mag_vectors.py` 的矢量层，不新写角度代码。
2. **`scripts/solve_D_dmi.py`（新增）**：读各构型的 `free energy TOTEN`（转发 `dft_energy.read_toten`）
   → 解出逐键 `D = (D_x,D_y,D_z)`；**必须**带两道闸门：① 能量差的**信噪比**闸门（见 §5.2 的精度要求）
   ② 与螺旋法（若同时做了）的**交叉一致性**闸门；不过就报错，绝不静默出 D。
3. **螺旋法备选路径**：`build_dmi_cases.py --spiral`（`q` 扫描）+ `solve_D_dmi.py --spiral`（取 `±q`
   的反对称部分拟合）。两条路径互为 L2 交叉验证（对应 `LATTICE_ROADMAP.md` §三.4 的验证阶梯）。
4. **`D_result.json`（新 schema，不复用 `j_result/1`）**：`{material, bonds:[{i,j,d, D:[x,y,z], |D|}],
   J_sym:{...}, cond, residual, source:"vasp"|"literature", convention:"H=+D·(Si×Sj)", confidence}`。
   落 `scripts/d_result.py`（镜像 `j_result.py` 的校验风格）。
5. **序参量扩展（MC 侧）**：`mc_2d.py` 加 `--field B`（塞进 `E = ... - gμ_B B·ΣS`，模型格式已留该形式）
   + 拓扑荷序参量 `Q = (1/4π)∫ m·(∂_x m × ∂_y m)`（格点离散版）；`write_hamiltonian.py` 的符号式 H
   也要加 DMI 项（**否则产物自述与模型不符**——它现在只写 `-J Σ S·S + MAE Σ(S^z)²`）。
6. **判据/门禁**：用户必须显式 `--dmi` 才走这条链；`detect_lattice` 保持只管几何（**不要**试图自动判定
   对称性——那需要空间群分析，属另一个工具）。

### 4.2 MTI 子链（≈5 件）

1. **多 QL slab 构建器**：`build_mti_cases.py`——把 N 个 QL 沿 c 叠起来（层内 FM、层间 A 型 AFM），
   层分辨 `magmom.txt`；同时要给 `reorient`/真空处理加「**非单层**」模式（现有链路默认单层 + 真空）。
2. **层分辨磁矩解析**：`dft_energy.py` 加 `read_layer_moments(path, layers)`（样板 = `fgt_verify_states.py:54`）。
3. **`solve_Jc.py`**：层间 `J_c` = (E_FM − E_AFM)/（自旋对数的差）；带残差/一致性闸门。
4. **层状 MC**：模型格式 `dim=2 → dim=3`、`mc_2d.build_arrays` 接受 `(Lx, Ly, Lz)`、
   `poscar_to_mc3` 写层间壳 ⇒ 才能算「由 `J_c` 控制的三维 `T_N`」。这是**中等改动**，不是重写
   （格点图/近邻表机制可复用），但会碰 2D 假设的边界（`mc_fss.py` 的有限尺寸标度也要跟着变）。
5. **拓扑部分**：Z₂ / 表面态 / 轴子态 **建议不进本仓**——交给 Wannier90/WannierTools/pyprocar 专门链；
   本仓最多做「输入准备（slab POSCAR + SOC INCAR）+ 引用外部结果」，并在结论文本里写明来源。

### 4.3 契约影响

- **不改既有语义**：新产物另开 `D_result.json` / `Jc_result.json`；`J_result.json`（`j_result/1`）保持不动。
- **`poscar_to_mc3`**：只需保证 `J_mat` 反称部分被消费（**已成立**）与 `dim=3` 扩展（MTI 才需要）。
- 两类的**门禁**都要走 `ERROR_CATALOG` + `KNOWN_ISSUES.md` 两表一一对应的既有做法。

---

## 5. 可行性备忘（立项前置 + 机时）

### 5.1 前置条件

| 类别 | 需要什么 | 现状 |
|---|---|---|
| DFT 代码 | 支持非共线 + SOC；DMI 还需**约束矩**（VASP 的 `I_CONSTRAINED_M`）或**自旋螺旋**（`LSPIRAL/QSPIRAL` 类功能），或等价的 QE 路径 | **本机无 VASP**，从未跑过任何 DFT |
| 机时 | 见 §5.3（DMI 量级 10¹–10² 小时；MTI 10²–10³ 核时） | 无调度器接入 |
| 人员/知识 | 四态法或螺旋法提取流程；atomistic spin dynamics（skyrmion）| 仓库内无先例 |
| 软件（MC 侧） | DMI 的 MC 已可用（§3.2）；skyrmion 的 B–T 相图需给 `mc_2d` 加外场 + 拓扑荷（§4.1 第 5 件） | 缺 |
| 软件（拓扑侧） | Wannier90 / WannierTools / pyprocar 等 | 缺，且**建议不进本仓** |
| 数值精度 | DMI 的能量差 ~1e-3–1e-2 eV，而总能 ~-1e2 eV ⇒ 需要 `EDIFF=1E-7`（**主链的 1E-6 不够**）、固定 `PREC`/FFT 网格、同一 k 网格 | 现有 INCAR 模板是 `EDIFF=1E-6` |

### 5.2 构型数量（**可严格推出**，不依赖机器）

| 任务 | 需要的 SCF 数 | 说明 |
|---|---|---|
| DMI 四态法 | `4 × n_键` | 每个不等价键 4 个取向组合；kagome 有 4 个不等价键（J1/J1′/J2/J3）⇒ 16 |
| DMI 螺旋法 | `2 × n_q` | 必须取 `±q` 才能分离反对称部分；`n_q = 20–40` ⇒ 40–80 |
| DMI 参考态 | 1–2 | FM/基态能量零点 |
| MTI 层间 `J_c` | 2–4 | 层内 FM × 层间 FM/AFM 两种，附 k 网格/层数收敛 |
| MTI 拓扑表征 | 1 次 dense-k NSCF + Wannier 流程 | 不是 SCF；成本随 k 点数 × 能带数 |

### 5.3 机时估算（**纯估算：假设 × 数量，未实测**）

> **假设块（必须与数字一起引用）**：单次 SOC 非共线 SCF，`N ≈ 60` 原子、`ENCUT = 520 eV`、
> `4×4×1` k 网格、32 核 ⇒ **2–6 h（≈ 64–192 核时）**。缩放律：随电子数/原子数 `N²–N³`、随 k 点数线性；
> SOC 相对共线通常 ×1.5–3；含 f 电子（+U、多组态）再 ×2–4。**这三个系数都是经验区间，不是本机实测。**

| 场景 | 推算 | 量级 |
|---|---|---|
| DMI（四态法，kagome 4 键） | 16 × (2–6 h) × 冗余 1.5–2 | **50–200 h ≈ 10¹–10² 核时级**（数十核时到数百核时） |
| DMI（螺旋法，20–40 q） | 40–80 × (螺旋用原胞，单次更便宜，按 0.5–2 h) | **20–160 h** |
| MTI（3 QL ≈ 100–150 原子） | 单次 12–60 h；`J_c` 2–4 次 + 层数/网格收敛 ×2–4 | **10²–10³ 核时** |
| MTI 拓扑流程 | dense-k NSCF + Wannier + 层数收敛 | **10²–10³ 核时**（且需额外软件授权/学习成本） |

**校准方法（把估算换成实测的唯一正路）**：在你的机器上跑一次 FM 参考 SCF，记录 **wall-clock 与核数**
⇒ 得到 `c₀`（核时/构型）；其余按 §5.3 的缩放律外推；**把实测值写回本表并标注日期**。本表允许被实测替换，
但**在替换之前不许把上表当作预算依据**。

### 5.4 风险与替代路径

| 风险 | 缓解 |
|---|---|
| DMI 能量差淹没在数值噪声里（D ~ 0.1–2 meV vs 总能 ~ -100 eV） | `EDIFF=1E-7`、固定 FFT/PREC、同一 k 网格；四态法与螺旋法**交叉验证**；报 D 时必须给残差 |
| VASP 的自旋螺旋支持随版本/设置变化 | 优先四态法（不依赖特殊功能）；螺旋法作为交叉验证而非唯一路径 |
| MC 侧缺外场/拓扑荷 ⇒ 做不出相图 | 先只报零场基态与能量；要相图必须先做 §4.1 第 5 件 |
| MTI 的多 QL 与「单层 2D」假设冲突 | 明确另立层状子链；不要用单层 MC 的数字冒充 `T_N` |
| 机时/软件都不到位 | **替代路径**：① 用文献 `D` 走 §3.2（今天就能跑，但 D 是外来的）② 只做 §2 的识别 + 明确拒绝（**P6 现状**）③ 交给外部专门工具链并把结果引用进来 |

---

## 6. 文献锚点（用于判别，不用于代替计算）

- DMI 的对称性起源与形式：Dzyaloshinskii (1958)；**Moriya, Phys. Rev. 120, 91 (1960)**（各向异性超交换）。
- **四态法**提取 D：Xiang, Lee, Kuneš 等，Phys. Rev. B **84**, 024404 (2011) 及后续（按构型能量差解 D）。
- 自旋螺旋法提取 J(q)/D(q)：见 VASP 的 `LSPIRAL/QSPIRAL` 相关文档与文献；**本仓未验证**。
- MTI 原型：**Otrokov 等, Nature 576, 416 (2019)**（MnBi₂Te₄ 的预测与观测，A 型反铁磁拓扑绝缘体，`T_N ≈ 25 K`）。
- 具体数字以原文为准；`literature_J.json` 里**没有** D 字段，也**没有** MTI 条目 —— 引用时写明来源与年份。
