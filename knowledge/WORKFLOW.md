# 二维磁性单层通用流水线 — MP → 磁耦合 J → Monte Carlo

给一个二维磁性材料的化学式，自动完成：**Materials Project 下载结构 → 层对齐 → 超胞/四磁态 → 磁耦合常数 J(J1/J2/J3) 能量映射 → 经典 Heisenberg Monte Carlo 临界温度 (T_C / T_N)**。

> 从 `MnPS3_pipeline` 抽象而来（它是已验证实例）；继承了 `CrBr3_pipeline` 的参数化脚本与 `mp_2d` 的下载脚本。之前的 CrI3/CrBr3 流水线只支持 J1/J2 + 铁磁总磁化，本框架扩展了 **J3** 与 **反铁磁交错序参量**，并把 MC 引擎重写成自包含版本（不再依赖外部 `MnSTe_repro` 仓库）。

---

## 1. 一句话运行

```powershell
cd <repo>

# 命名模型预设（最简，一键设定 lattice/S/moment/U；文献 J 仅作参考）：
python run_2d_magnet.py --model mnps3 --poscar tests/fixtures/MnPS3_POSCAR

# 已收录文献 J 的材料（如 MnPS3）：
python run_2d_magnet.py --formula MnPS3 --literature MnPS3 --mp-key <MP_API_KEY>

# 新材料（手动给参数）：
python run_2d_magnet.py --formula VI3 --element V --nonmag I --S 1.0 --moment 2.0 --U 3.0 `
    --J1 -1.0 --J2 -0.1 --J3 0.0 --MAE 0.2 --mp-key <KEY>

# 已有 POSCAR（跳过下载）：
python run_2d_magnet.py --poscar my.vasp --element Cr --nonmag I --S 1.5 --moment 3.0 --U 3.0 `
    --J1 2.4 --J2 0.6

# 无 MP key，改用公开 provider 链（COD/OQMD/AFLOW/JARVIS，返回体相）：
python run_2d_magnet.py --formula MnPS3 --literature MnPS3 --provider cod
```

输出到 `<公式>_pipeline/`（`--outdir` 可改）：`workflow.md`、`hamiltonian.md`（顶层海森堡哈密顿卡片）、`01_structure/`、`02_supercell/`、`03_dft/`、`05_mc/`（含自包含的 `mc_2d.py` + `poscar_to_mc3.py` + `write_hamiltonian.py`）、`figures/`。

---

## 2. 输入参数

| 参数 | 含义 | 例 (MnPS3) |
|---|---|---|
| `--model` | 命名预设，一键设定下方 element/nonmag/lattice/S/moment/U（`mnps3/cri3/crbr3/fe3gete2/crsbr/gdi2/la2cuo4/k2nif4/crocl`） | `mnps3` |
| `--formula` | MP 化学式 | `MnPS3` |
| `--element` | 磁性元素 | `Mn` |
| `--nonmag` | 非磁元素（可多个，空格分隔） | `P S` |
| `--S` | 自旋大小 | `2.5`（Mn²⁺ 5/2） |
| `--moment` | 磁矩初猜 μB | `5.0` |
| `--U` | DFT+U 加磁性元素 d 轨道（f 块元素自动用 f，eV） | `5.0` |
| `--J1/--J2/--J3` | 交换常数 meV，约定 `H = -J S·S`（J>0 铁磁） | `-0.79/-0.04/-0.23` |
| `--MAE` | 单离子各向异性 meV（默认 0） | `0.0` |

**J 符号约定**：`H = -J Σ S·S`，**J>0 = 铁磁，J<0 = 反铁磁**。多数文献用 `H = +J S·S`（J>0=反铁磁），转过来**要变号**。

---

## 3. J 值从哪来（关键）

1. **优先查** `literature_J.json`（**已核实**文献值：MnPS3 / CrSBr / Fe3GeTe2 / La₂CuO₄；CrI₃ / CrBr₃ / GdI₂ / K₂NiF₄ / **CrOCl** 按惯例留 `null` + `_note`）。
2. 否则**联网查文献**（Sivadas / Chittari / Wildes 等），或让用户给，**绝不猜**。
3. 或**真跑 DFT**：`03_dft` + `runs/<order>` 已备好，四磁态静态 + FM 弛豫后
   `scripts/solve_J3.py <cases> --demo` 自检 → 用真实能量反解。

   **自旋模式**（`dft` 阶段可选，P3 起）：
   - 默认 = 共线（`ISPIN=2` + 1 分量 MAGMOM）。
   - `--noncollinear` = 非共线（`LNONCOLLINEAR=.TRUE.` + 3 分量 MAGMOM），用于 120° / 螺旋 / 倾斜。
   - `--soc` = 在非共线上再加 `LSORBIT=.TRUE.` + `GGA_COMPAT=.FALSE.`（可配 `--saxis x y z` 指定自旋量子化轴）；
     驱动日志与 `workflow.md` 都会写明「非共线+SOC」。**需要 SOC 的磁序（cant/spiral/tilt…）没开 `--soc`
     会被硬拒**，不会静默降级成共线去算一个物理上不成立的能量。
   - 三角晶格走 `assemble_dft_tri.py`，它已退化为**薄适配层**（内部转调 `assemble_dft.py --noncollinear
     --all-orders`），两条路径共用同一份 INCAR 渲染，行为不会分叉。

   **能量口径（P3 起统一）**：所有 `solve_*` 一律读 OUTCAR 的 **`free energy  TOTEN`**、取**最后一次**出现；
   缺这一行就报错，**不会**悄悄改用 `energy(sigma->0)`——两种口径混用会让反解出的 J 差一个常数项
   （`dft_energy.py` 的 `mixed_convention()` 专抓这种混用，阈值 1e-3 eV）。E.txt 兜底只在没有 OUTCAR 时启用。

   **收敛判据**（`03_dft/INCAR.static` 默认已给；改了就重新对账）：`EDIFF=1E-6`、`NELM=200`、`ENCUT=520`、
   `ISMEAR=0 + SIGMA=0.05`（金属体系自行改 `ISMEAR=-5`）。非共线 + SOC 必须 `GGA_COMPAT=.FALSE.`
   （否则 SOC 被 GGA 兼容模式抵消）且**删掉 `ISPIN`**（与 `LNONCOLLINEAR` 互斥）。

> `--model` / `--literature` 落回文献 J 时，driver 只在 stderr 打 `[ref]` 声明「文献参考、非本项目计算」，workflow.md 也会标注来源。**文献 J 是数量级 sanity-check，不是本项目结论**：真实 J 必须跑 VASP + `solve` 反解，或用 `--J1/--J2/--J3` 显式给定。VASP 侧只**导出**输入（`03_dft/` + `runs/<order>/`），不代提交——请在集群手动提交。

新材料写进 `literature_J.json`（字段 `formula/element/nonmag/S/moment/U/J1/J2/J3/MAE/ground_state/source`），并在 `run_2d_magnet.py` 的 `MODEL_PRESETS` 里加同名键（含 `lattice`）以便 `--model` 一键调用。**J 未核实的条目 `J1` 留 `null` 并写 `_note`**。

---

## 4. 内部步骤（driver 自动调用）

```
fetch      mp_fetch_2d.py           MP 下载 + 剥单层（取 hull 最小的候选）
           fetch_structure.py       结构 provider 链: COD/OQMD/AFLOW/JARVIS (--provider, 无 key, CIF/OPTIMADE/AFLUX -> POSCAR)
reorient   reorient.py              层法线对齐 z + 加真空（修 MP 单斜/紧胞缺口）
standardize standardize_cell.py     紧凑胞标准化（把胞内两磁性原子变成最近邻；MP 原胞可能非紧凑，如 CrI3 的 8.05 Å）
supercell  build_magnetic_cases3.py 2×2×1 超胞 + 四磁态 + J1/J2/J3 系数
          build_cluster_cases.py   簇磁体：fm+波矢态 + J_intra/J1/J2/J3
          build_magnetic_cases_tri.py 三角：fm + 120°非共线 + stripy/uud + J1/J2/J3
          build_rectangular_cases.py 矩形：fm/afm_a/afm_b/afm_diag 四共线 + J1/J2/J3 (距离壳升序)
           build_square_cases.py    方形：fm/g_type/stripe_a/stripe_b/p4_a/mix_ab 六共线枚举（DFT 只提交选定四态）+ J1/J2/J3 (a / a√2 / 2a)
           build_crect_cases.py     C 心矩形（2 磁原子/胞，如 CrOCl）：7 共线枚举（DFT 只提交选定四态）+ J1/J2/J3 (C 心 3D 键长 / a / b)
                                    + 非磁元素占位与化学计量校验（--nonmag 写错/漏写硬报错，每胞计量落盘）
           build_kagome_cases.py    六角 kagome（3 磁原子/胞、最近邻 a/2 上 4 配位，如 Fe3Sn）：7 磁态
                                    （fm/q0_120/q0_ud/q0_duu/sq3_120/sq3_updntri/sq3_alt，非共线 3 分量）+ J1(a/2)/J2(a√3/2)/J3(a)
                                    + 硬拒 q=0 简并（rank<4 / cond>1e4 / det≈0）；超胞 n0/n1 须为 3 的倍数（P4.1）
dft        assemble_dft.py / assemble_dft_tri.py  生成 INCAR/KPOINTS/POTCAR.spec + runs/<order>/ (三角用非共线)
           自旋模式: 默认共线(ISPIN=2) / --noncollinear(3 分量) / --soc(+LSORBIT,GGA_COMPAT=.FALSE.,--saxis)
           assemble_dft_tri.py 只是适配层, 转调 assemble_dft.py --noncollinear --all-orders (共用同一份 INCAR 渲染)
mc         poscar_to_mc3.py + mc_2d.py  建模型 + Heisenberg MC（自写 fast kernel）
           之后 driver 自己再判一次峰位: 峰落在温度窗口首/末格点 (= 窗口截断或网格跨过窄峰)
           → 自动缩放窗口 + 加密退火重扫, 首轮留 `mc_mt_survey.csv` (避免把边界值当 T_C/T_N)
fss        mc_fss.py (--fss 时)  多尺寸 U₄(T) 有限尺寸标度: 两两交点 + collapse 交叉验证 → T_C/T_N
render     render_supercell.py / render_cluster_orders.py / render_tri_orders.py  磁序俯视图
           (本机缺 matplotlib 时三个分支自动改走等价的 *_pil.py: supercell/cluster/tri)
           kagome 复用 tri 渲染器（非共线 3 分量箭头；render_supercell.py 只认共线, 会画错）
                     标题由 driver 用 --title 覆盖（tri 渲染器默认标题写死「三角晶格」）
ship-mc    (复制) mc_2d.py + poscar_to_mc3.py + write_hamiltonian.py (簇/三角/kagome 再加对应枚举+反解脚本) -> <outdir>/05_mc/
hamiltonian write_hamiltonian.py    生成 <outdir>/hamiltonian.md (Heisenberg 哈密顿卡片, 顶层)
```

脚本均在 `scripts/`，可单独调用。`mc_2d.py` 自包含（numpy+numba+matplotlib）。

磁序图有**两条等价路径**：matplotlib 版（`render_supercell.py` / `render_cluster_orders.py` / `render_tri_orders.py`）与纯 Pillow 版（`render_supercell_pil.py` / `render_cluster_orders_pil.py` / `render_tri_orders_pil.py`，同 CLI / 输入 / 输出 / 配色），共用原语 `scripts/pil_fig.py`。driver 用 `have_matplotlib()` 现场探测，装了 matplotlib 走前者，否则三个分支各自自动改走后者并打 `[warn] matplotlib 不可用 ... 改用 Pillow 兜底`。曲线类图同理各带一条内置 Pillow 路径：MC 曲线 `mc_2d.py` 的 `_plot_mt_png`（两面板版式）、FSS 标度图 `mc_fss.py` 的 `_plot_fss_png`、FGT 磁序图 `build_fgt_cases.py` 的 `_plot_fgt_png`。**六条渲染路径全部与 matplotlib 解耦**：matplotlib 缺席时只要装了 Pillow 就有 PNG。

driver 另支持位置参 `<command>` 分阶段运行：`fetch` / `build` / `dft` / `mc` / `solve` / `render` / `report` / `doctor`。各子命令在阶段边界早退、复用已有产物（如 `mc` 不产 DFT 输入、`report` 不触发 MC/图）；`doctor` 做环境自检并把命中的已知问题映射到 `KNOWN_ISSUES.md` 的「症状→修复」（含 `--model` 预设自检：每个预设声明的 `literature` 锚点必须与 `literature_J.json` 的 formula/element/nonmag/S/moment/U 逐项一致，不一致会以 `preset_mismatch` 命中）；省略 `<command>` 即完整流水线。

---

## 5. 目录结构

```
2d_magnet_pipeline/
├── run_2d_magnet.py     # 单入口 driver
├── WORKFLOW.md          # 本文件
└── scripts/
    ├── __init__.py      # 使 scripts/ 可作包安装 (非可编辑 pip install . 也能定位脚本)
    ├── element_data.py  # f 块元素表 + ldaul_for()：组装脚本按元素判 DFT+U 轨道 (d=2/f=3)
    ├── literature_J.json # 已核实文献 J 表（MnPS3 / CrSBr / Fe3GeTe2 / La₂CuO₄ 已验证；CrI₃ / CrBr₃ / GdI₂ / K₂NiF₄ / CrOCl 未核实留 null + _note）
    ├── mp_fetch_2d.py   / reorient.py / standardize_cell.py
    ├── build_magnetic_cases3.py / solve_J3.py
    ├── poscar_to_mc3.py / mc_2d.py (--order cv 比热判据)
    ├── write_hamiltonian.py # 模型 JSON -> hamiltonian.md (海森堡哈密顿卡片)
    ├── cluster_sites.py # 簇磁体：磁性簇→有效位（Nb3I8 等）
    ├── build_cluster_cases.py # 簇磁体：超胞+磁序枚举+J_intra/J1/J2/J3 映射 (含 DFT runs)
    ├── solve_cluster_J.py # 簇磁体：OUTCAR -> E0/J_intra/J1..Jk
    ├── render_cluster_orders.py # 簇磁体：簇心磁序俯视图
    ├── build_magnetic_cases_tri.py # 三角：超胞 + 磁序枚举(含120°非共线) + J1/J2/J3 映射
    ├── assemble_dft_tri.py # 三角：**薄适配层** —— 转调 assemble_dft.py --noncollinear --all-orders (共用 INCAR 渲染)
    ├── solve_J_tri.py # 三角：OUTCAR -> J1/J2/J3
    ├── render_tri_orders.py # 三角：磁序俯视图 (120° 矢量箭头)
    ├── build_rectangular_cases.py # 矩形(CrSBr/CrSCl)：超胞 + 四共线态(fm/afm_a/afm_b/afm_diag) + J1/J2/J3 映射
    ├── solve_J_rect.py     # 矩形：OUTCAR -> J1/J2/J3（距离壳升序, 对齐 build_rectangular_cases --selftest）
    ├── build_fgt_cases.py # **双亚晶格铁磁网络 (Fe3GeTe2/Fe3GaTe2 类)**：超胞(默认 3x3, 含别名自检)
    │                      #   + 磁序枚举 + 系数矩阵 + 逐原子 MAGMOM 的 INCAR/KPOINTS/runs + 俯视图 + 哈密顿卡片
    ├── solve_fgt_J.py     # **双亚晶格反解**：各磁态总能量 -> J1..J4 + E0（`--demo` 自检, 应 PASS）
    ├── build_square_cases.py # 方形(K₂NiF₄ 型, a=b 且原胞 1 个磁原子)：超胞(默认 4x4x1, 前两维须相等且为 4 的倍数)
    │                      #   + 六共线态 fm/g_type/stripe_a/stripe_b/p4_a/mix_ab (DFT 只提交选定四态) + J1(a)/J2(a√2)/J3(2a) 映射
    │                      #   + 产物含 order_spec/model_spec/design_cond；`--selftest` 与 P0.1 规格层逐位交叉校验
    ├── solve_J_square.py  # 方形反解：OUTCAR -> J1/J2/J3 + E0 → J_result.json（`--demo`/`--ref la2cuo4` 自检, 应 PASS）
    ├── build_kagome_cases.py # **六角 kagome (Fe3Sn 类, 原胞 3 个磁原子: 最近邻 a/2 上 4 配位, 不是 FGT)**
    │                      #   + 超胞(默认 3x3x1, 前两维须为 3 的倍数) + 七磁态 fm/q0_120/q0_ud/q0_duu/sq3_120/sq3_updntri/sq3_alt
    │                      #   + J1(a/2)/J2(a√3/2)/J3(a) 映射; 硬拒 rank<4 (q=0 简并) / cond>1e4 / det≈0
    │                      #   + 产物含 order_spec/model_spec/design_cond/coefficients_per_mag; `--selftest` 与规格层逐位交叉校验
    ├── solve_J_kagome.py  # kagome 反解：OUTCAR -> J1/J2/J3 + E0 → J_result.json（`--demo --synthetic 3 -0.4 0.15` 自检, 应 PASS）
    │                      #   rank<4 / cond>1e4 直接 SystemExit（P4.1 硬门禁, 不静默出 J）; J1 实际是 J1/J1′ 的平均（3 壳模型）
    ├── build_crect_cases.py # **C 心矩形 (CrOCl 类, 正交 a≠b 且原胞 2 个磁原子)**：超胞(默认 4x4x1, 前两维须为偶数)
    │                      #   + 七共线态 fm/afm_c/afm_a/afm_b/afm_diag/afm_c_a/afm_c_diag (DFT 只提交选定四态)
    │                      #   + J1(C 心 3D 键长)/J2(a)/J3(b) 映射, 壳容差默认 0.05（J1/J2 只差 0.1632 Å）
    │                      #   + 硬拒 basis 退化组合(det=0) 与 cond>1e4；产物含 order_spec/**per_site_bonds**/stoichiometry_per_cell
    │                      #   + 非磁元素占位校验: --nonmag 漏写胞内元素或写了胞里没有的元素都硬报错, 每胞计量须为整数
    ├── solve_J_crect.py   # C 心矩形反解：OUTCAR -> J1/J2/J3 + E0 → J_result.json（`--demo --synthetic 5 -0.5 0.2` 自检, 应 PASS）
    │                      #   文献 J 未核实时**拒绝**用 null/他材料的值当锚点（CrOCl 就是这种）
    ├── solve_J_crect_cant.py # **弱铁磁倾斜反解 (P2')**：同一 cases 目录, 但只吃「倾斜家族 + 张成所需的 fm/母序行」
    │                      #   闭式 sin²θ 标度 + 自由度审计(`--audit-only`, 不需要 OUTCAR); θ 是**输入参数**不是基态解
    │                      #   倾斜序属 SOC_ORDERS ⇒ 装配必须 --soc; 一次倾斜 SCF 定不出唯一 (J, D)
    ├── mag_vectors.py   # **自旋矢量层 (P3)**：共线 / 面内 120° / ±θ 倾斜 / 螺旋 -> 3 分量 MAGMOM
    │                    #   + SOC_ORDERS/needs_soc/ANGLE_TABLE + 渲染/解析(expand_repeats)/校验(check_magmom)
    ├── dft_energy.py    # **OUTCAR 能量口径 (P3)**：一律 free energy TOTEN(取最后一次) + mixed_convention 混用检查
    ├── dmi_matrix.py    # **D 矢量 -> 3x3 反称交换矩阵 (P6)**：J_D = -[D]x <=> H_DM = +D·(S_i×S_j)
    │                    #   产物直接喂 `poscar_to_mc3 --shell-def` 的逐壳 J_mat；**本仓不提取 D**（文献/自有 DFT）
    ├── altermagnet_check.py # **交变磁体几何判据 (P5)**：两套亚晶格的连接操作 = 平移(常规AFM) / 反演(宇称禁止) / 旋转(几何候选)
    ├── ingest_dft.py    # **OUTCAR 回收通道 (P4.2)**：runs/ -> cases/，11 条硬闸门**全过才写**(原子性) + sha256 provenance
    ├── dft_handoff.py   # **外部机时交接 (P4.2)**：pack(submit.sh/MANIFEST/POTCAR.spec/CHECKLIST) / ingest；不打包 POTCAR、不代提交
    ├── k_convergence.py # **k 点收敛抽查 (P4.2)**：runs/<序>/ -> kconv/runs/<序>@<k网格>/（只换 KPOINTS，INCAR 逐字节不变）
    │                    #   + OUTCAR 放回后判「总能量对 k 网格收敛到哪个网格够了」（≠ 判 J 收敛，不代提交作业）
    ├── assemble_dft.py  # **VASP 输入统一入口**：默认共线 / --noncollinear / --soc(+--saxis) / --all-orders
    ├── render_supercell.py / render_supercell_pil.py
    ├── render_cluster_orders.py / render_cluster_orders_pil.py   # 簇磁体: 簇心视图 (matplotlib / 纯 Pillow)
    ├── render_tri_orders.py     / render_tri_orders_pil.py       # 三角: 120° 面内自旋箭头 (matplotlib / 纯 Pillow)
    ├── pil_fig.py          # 三条 Pillow 渲染路径共用的原语 (面板变换/虚线/圆/箭头/图例); 新渲染器一律复用它
```

---

## 6. 从哪些项目继承
| 通用框架脚本 | 来源 | 改动 |
|---|---|---|
| `mp_fetch_2d.py` | `mp_2d/` | 原样复用 |
| `build_magnetic_cases3.py` | CrBr3 `build_magnetic_cases(2).py` | J1/J2→J1/J2/J3、多非磁元素 |
| `solve_J3.py` | CrBr3 `solve_J.py` | J1/J2→J1/J2/J3 |
| `render_supercell.py` | CrBr3 `render_supercell.py` | 加多非磁元素 + J1/J2/J3 壳分色/系数标注 |
| `render_supercell_pil.py` | 新写（本仓, 2026-10） | 纯 Pillow 的等价路径：本机 matplotlib 3.11.1 缺 pyparsing/cycler/contourpy/fontTools/kiwisolver 且离线装不上 ⇒ 超胞图改由它出；几何/文本规则复用 `render_supercell.py`（该模块已改成**顶层不 import matplotlib**），CLI 与输出路径完全一致；C 心矩形还会把「MC 面内投影尺子」与「规格层 3D 键长尺子」两个壳距都标进图例 |
| `render_cluster_orders_pil.py` | 新写（本仓, 2026-10） | 簇磁体的纯 Pillow 等价路径（同 CLI/输入/输出）；`render_cluster_orders.py` 已改成**顶层不 import matplotlib**，几何/配色/文本规则全部复用；非磁原子按 matplotlib 版 zorder 画在键之下 |
| `render_tri_orders_pil.py` | 新写（本仓, 2026-10） | 三角的纯 Pillow 等价路径；120° 非共线自旋画成**定长方向箭头**（0.30·\|a\|，只表方向），与 matplotlib 版同步改掉了原 `quiver(scale=0.9)` 在 \|S\|=5 时画出 5.6 Å 长箭头互相纠缠的问题 |
| `pil_fig.py` | 新写（本仓, 2026-10） | 三条 Pillow 渲染路径共用的绘图原语：面板网格/数据→像素变换/手工虚线点线/实心圆（两段画法规避 Pillow 的 `ellipse(width>1)` 十字伪影）/实心箭头（三角头）/图例排布；`render_supercell_pil.py` 早于它（自带一份拷贝），之后新写的 Pillow 渲染器一律复用它 |
| `mag_vectors.py` | 新写（本仓, 2026-10, P3） | 自旋矢量层：共线 / 面内 120° / ±θ 倾斜 / 螺旋 → 3 分量 MAGMOM，含 `SOC_ORDERS`/`needs_soc`/`ANGLE_TABLE` 与渲染/解析/校验（从 `assemble_dft_tri.py` 里抽出，供全部 `assemble_dft*` 与 P4 kagome 复用） |
| `dft_energy.py` | 新写（本仓, 2026-10, P3） | OUTCAR 能量口径唯一出口：一律 `free energy TOTEN` 取最后一次，缺则报错（不静默改用 `energy(sigma->0)`）；`mixed_convention()` 抓两种口径混用；六个 `solve_J*` 的 `read_toten` 全部转发到它 |
| `dmi_matrix.py` | 新写（本仓, 2026-10, P6） | DMI 的等价矩阵形式：`J_D = -[D]×`（⇔ `H_DM = +D·(S_i×S_j)`；反号即帮倒忙），矩阵吸收 `S²`，产物即 `poscar_to_mc3 --shell-def` 的逐壳 `J_mat`（反称部分由引擎的 `J.T` 反向键正确消费）。`--selftest` 5/5，含与 `mc_2d._energy` 的实机对拍（顺带钉死「`_energy` 是每格点能量 ⇒ 键能量 = `_energy × N_sites`」这条坑）。**D 必须外来**（文献或用户自有 DFT）：本仓不提取 D、不做拓扑判定（见 `DMI_MTI_MEMO.md`） |
| `assemble_dft.py` | CrBr3 `assemble_runs.py` 重写 | 参数化 U/LDAUL/物种；**P3 起是唯一的 INCAR 渲染出口**：`--noncollinear / --soc / --saxis / --all-orders`，`ISPIN` 与非共线互斥，3 分量 `magmom.txt` 配共线模式会硬报错 |
| `poscar_to_mc3.py` | MnPS3（源外部 `MnSTe_repro` 仓库） | 加 J3 + 修 bug |
| `mc_2d.py` | 新写 | 自包含 MC + 铁磁/反铁磁序参量 |
| `write_hamiltonian.py` | 新写 | 模型 JSON -> `hamiltonian.md` (海森堡哈密顿卡片) |
| `reorient.py` | 新写 | 修 MP 层法线/真空缺口 |
| `standardize_cell.py` | CrBr3 `standardize_cell.py` | 原样复用：把胞内两磁性原子标准化成最近邻（MP 原胞可能非紧凑，如 CrI3 胞内 8.05 Å = J3 距离，会破坏 M 点枚举） |
| `cluster_sites.py` | 新写 | 簇磁体：磁性簇识别 → 簇质心作为有效磁位（Nb3I8 的 Nb3→1 位） |
| `build_cluster_cases.py` | 新写 | 簇磁体：超胞 + 磁序枚举 + J_intra/J1/J2/J3 映射（含 runs/） |
| `solve_cluster_J.py` | 新写 | 簇磁体：OUTCAR -> E0/J_intra/J1..Jk（对齐 build_cluster_cases --selftest） |
| `render_cluster_orders.py` | 新写 | 簇磁体：簇心磁序俯视图 |
| `build_magnetic_cases_tri.py` | 新写 | 三角：超胞 + 磁序枚举(含 120° 非共线) + J1/J2/J3 映射 |
| `assemble_dft_tri.py` | 新写（P3 起改为**薄适配层**） | 三角：只拼命令行转调 `assemble_dft.py --noncollinear --all-orders`（与共线路径共用同一份 INCAR 渲染）；`--print-cmd` 可看它实际执行的命令 |
| `solve_J_tri.py` | 新写 | 三角：OUTCAR -> J1/J2/J3 |
| `render_tri_orders.py` | 新写 | 三角：磁序俯视图 (120° 矢量箭头) |
| `build_rectangular_cases.py` / `solve_J_rect.py` | 新写（本仓, 2026-10） | 矩形晶格铁磁（CrSBr/CrSCl 类）：4 共线态 fm/afm_a/afm_b/afm_diag 枚举 + 系数矩阵 + OUTCAR 反解 J1/J2/J3（`--selftest`/`--demo` PASS）；J 按距离壳升序（J1=a 轴/J2=b 轴/J3=对角） |
| `mc_2d.py --order cv` | 新写 | 比热 Cv 判据，受挫/非二分晶格（三角）通用定 T_N |
| `mc_2d.py --order sublattice:<mask>` + `altermagnet_check.py` | 新写（本仓, 2026-10, **P5**） | **逐子晶格 Néel 矢量**序参量（`mask={0}` 即老的交错序参量 ⇒ 蜂窝/FGT 逐位等价，是回归锚点；`num_basis>=3` 时可取任意子晶格组合，而 `stag_signs` 只会把 basis 1..n-1 并成一组并打 warn）。产物仍是 4 列 `mc_mt.csv`（`T_Kelvin,order_param,Susceptibility,Binder_U4`）⇒ `summarize_mc`/`mc_fss` 的列契约不变；`mask` 必须留至少一个子晶格在 B 组，否则除零（脚本硬报错）。**只给标量 \|L\|，不给两套亚晶格的 Néel 矢量方向**。配套 `altermagnet_check.py` 判几何（平移/反演/旋转连接，容差以 Å 计），`--selftest` 15/15；判别结论：MnPS3 几何候选、CrOCl 宇称禁止、Fe3Sn 不适用 |
| 峰边界检测 + 自动缩放重扫 | 新写（本仓） | 内置在 `run_2d_magnet.py`：峰落窗口首/末格点 → 缩放加密退火重扫，首轮留 `mc_mt_survey.csv`；`mc_2d.py` 自身也会对边界峰打警告 |
| `build_fgt_cases.py` / `solve_fgt_J.py` | 新写（本仓, 2026-10） | **双亚晶格铁磁网络**（Fe3GeTe2 类）：A(2c)+B(4e) 建模、磁序枚举、系数矩阵(M=[1\|C], cond≈178)、逐原子 MAGMOM、俯视图、哈密顿卡片、总能量反解 J1..J4（`--demo` PASS） |
| `build_square_cases.py` / `solve_J_square.py` | 新写（本仓, 2026-10） | **方形晶格**（K₂NiF₄ / K₂CuF₄ / La₂CuO₄ 型, a=b 且原胞 1 磁原子）：6 共线态 fm/g_type/stripe_a/stripe_b/**p4_a**/mix_ab 枚举（DFT 只提交选定 4 个）、J1(a)/J2(a√2)/J3(2a) 系数、cond(A)≈40、`--selftest` 与 P0.1 规格层逐位交叉校验；**必须用 `policy="shift_search"`**（J3=2a 在 4 宽超胞里是半周期，最小镜像会把 ±2a 折叠 ⇒ J3 低估 2 倍），且**不引入周期 4 的 `p4_a` 就无法把 J3 与 e0 分离**（与 kagome 的 q=0 简并同构，成因不同） |
| `build_crect_cases.py` / `solve_J_crect.py` | 新写（本仓, 2026-10） | **C 心矩形晶格**（CrOCl 类, 正交 a≠b 且原胞 **2** 个磁原子, 两者由 (½,½) C 心矢量相连）：7 共线态 fm/afm_c/afm_a/afm_b/afm_diag/afm_c_a/afm_c_diag 枚举（DFT 只提交选定 4 个）、J1(C 心 3D 键长)/J2(a)/J3(b) 系数、cond(A)=64.02、`--selftest` 与 P0.1 规格层逐位交叉校验；**壳容差必须 0.05**（J1=3.0347 与 J2=3.1980 只差 0.1632 Å，tol=0.2 会静默并壳）；**basis 符号对 (n0,n1)-奇偶序完全不可见**（`afm_c_a == afm_a` 们逐位相同），母序与它的 basis 翻倍版同选必使 det=0，脚本硬拒；规格层新增 `order=afm_c`（胞内两 basis 反平行, C 心 Néel）；`solve_J_crect.py --demo` **必须**配 `--ref` 或 `--synthetic`（CrOCl 文献 J 未核实 ⇒ 拒绝猜锚点） |
| `build_kagome_cases.py` / `solve_J_kagome.py` | 新写（本仓, 2026-10, **P4.1 纯几何层**） | **六角 kagome**（Fe₃Sn 类, 六角 a=b 且原胞 **3** 个磁原子, 最近邻 a/2 上 **4** 配位 —— 配位数判据同时把 FGT（同距离 6 配位）挡在外面）：7 个磁态（q=0 三种 + √3×√3 三种 + fm, 均为 3 分量非共线）、J1(a/2)/J2(a√3/2)/J3(a) 系数、规格层 `kagome_fe3sn`（`policy="shift_search", nshift=1`）、`--selftest` 与规格层 3×3 与 6×6 超胞逐位交叉校验；**q=0 一族（fm/q0_120/q0_ud/q0_duu）天然简并**（c1==c2 且 c3 相同 ⇒ 子块 rank=2, cond≈2.3e34），必须靠含 √3×√3 的组合分离（默认 `selected = fm/q0_120/sq3_updntri/sq3_alt`，rank=4, cond=98.51），builder 与 solver 都对此**硬报错**；超胞 n0/n1 必须是 3 的倍数；**3 壳模型把 J1 与 J1′（同距离、不同亚晶格对）合并** ⇒ `J_result.json` 的 `J1` 是两者平均。**P4.1 只交付系数层**：无任何真实非共线 DFT 能量，J 与 T_N 属 P4.2 |
| `build_crect_cases.py` 的倾斜家族 + `solve_J_crect_cant.py` | 新写（本仓, 2026-10, **P2′ 系数层 + 反解器**） | **CrOCl 弱铁磁倾斜**：同一 cases 目录、第二个反解器（`solve_J_crect.py` 走共线族，本脚本走倾斜族）。物理内核是一条**闭式**关系 —— 同一母序倾斜 θ 后的系数 = 共线母序系数 + `sin²θ ×` 一个与 θ 无关的基向量（`canted_coeff`, 已与独立实现逐位对拍）；`--orders cant` = 倾斜家族 **+ 张成所需的 fm/母序行**（只给倾斜行会退化成 rank 2 的家族，`--audit-only` 打印每行系数与自由度、不需要 OUTCAR）。边界：**θ 是输入参数不是自洽解**（基态 θ 需 SOC 总能量极小，未做）、**一次倾斜 SCF 定不出唯一 (J, D)**、倾斜序在 `SOC_ORDERS` 里 ⇒ 装配必须 `--soc`、CrOCl 文献 J 未核实 ⇒ 真值路径必须 `--ref`/`--synthetic`。契约测试 `tests/test_crect_cant.py`（44 项）+ smoke §26 |
| `P4_2_KAGOME_RUNBOOK.md` + `ingest_dft.py` + `dft_handoff.py` + `k_convergence.py` | 新写（本仓, 2026-10, **P4.2 工具链**） | **kagome 物理层的操作手册与装配侧工具, 不是结果**：入口条件（真实机时 + 已核实文献锚点 + 参数表填满）/ 装配与已知坑（`ISYM` 默认只写在 `INCAR.relax`（static 不写），可用 `--isym` 显式覆盖；`NPAR` 默认 4，可用 `--npar` 覆盖；`KPOINTS` 默认 `3 3 1` 且表头已修正（不再写 2×2×1），可用 `--kpoints 6 6 1` 覆盖；`POTCAR` **不生成**，只给 `POTCAR.spec` 的元素顺序；理想化胞未弛豫）/ SCF 四条判据与失败协议 / 回收通道**已自动化**（`scripts/ingest_dft.py` 的 11 条硬闸门 + `scripts/dft_handoff.py` 的 pack/ingest 交接包，见手册 §5.0/§5.1）/ **k 点收敛抽查**（`scripts/k_convergence.py`，或 `run_2d_magnet.py kconv --k-grids 3,4,6,8`：生成 `<序>@<k网格>` 变体并**只换 KPOINTS**，OUTCAR 放回后判「总能量对 k 网格收敛到哪个网格够了」；**不代提交作业**，且总能量收敛只是 J 可信的**必要**条件）/ 反解可信度与「残差不为零先怀疑什么」/ 机时**纯估算**（替换成实测前不许当预算）/ 与原计划的偏差（含 P4.2 残余工作量 ≈4–9 人天的补估） |

---

## 7. 已知局限（结论里务必声明）

1. **蜂窝/三角/矩形/方形/C 心矩形/簇磁体**均支持磁序枚举 + DFT 输入 + 反解脚本 + 图（三角用非共线 120°；矩形用 4 共线态 fm/afm_a/afm_b/afm_diag，J 按距离壳升序，a≠b；方形用 6 共线态取 4 个 fm/stripe_a/**p4_a**/mix_ab，J1=a/J2=a√2/J3=2a，a=b；C 心矩形用 7 共线态取 4 个 fm/**afm_c**/afm_a/afm_b，J1=C 心 3D 键长/J2=a/J3=b，a≠b 且原胞 2 个磁原子）；**簇磁体**若是**单重态/dimer（如 Nb3X8）**，经典 Heisenberg 的 T_N **无物理意义**（仅机制演示），需两套 J 另建模，见 `dsh/2d-magnet/SKILL.md` 的「## 约定（务必遵守）」节里那条「簇磁体若低温是单重态/dimer（Nb3X8），经典 MC 的 T_N 无物理意义，仅机制演示」（**注意**：本仓从来没有 `references/` 目录，早先指向 `references/cluster-magnets.md` 的写法是悬空引用，见 `DMI_MTI_MEMO.md:67`）。
1c. **方形晶格（K₂NiF₄/K₂CuF₄/La₂CuO₄ 型）是 3D 氧化物、不可解理**：单层切面丢掉了层间超交换 ⇒ MC 的 T_C/T_N 只能当量级检查（La₂CuO₄ 实测 T_N=325 K 是 3D 驱动，2D 单层有限尺寸 χ 峰会低不少）。另外方形**必须**用显式平移搜索且**必须**含周期 4 的 `p4_a`，否则 J3 会低估 2 倍 / 与 e0 不可分离。
1b. **多亚晶格铁磁网络（Fe3GeTe2/Fe3GaTe2/Fe5GeTe2 类）不走通用分支**：几何并簇会把它误判成「分子簇」（Nb3I8 与 FGT 签名相同，都是 cluster/1 簇/3 原子），`build_cluster_cases.py` 会明确报错退出。请用 `build_fgt_cases.py` + `solve_fgt_J.py`。**超胞不能太小**：脚本内建 `r_max < 面内边长/2` 判据（2×2 对 3.99 Å 的键已坍缩）。
1d. **C 心矩形（CrOCl 类）有三个静默陷阱**：① 两个 Cr 只差 3.03 Å（< 簇判据 3.3 Å）⇒ 旧版 `detect_lattice` 会把它当分子簇磁体（现返回 `rect_centered`；方形网并簇返回 `sq_centered` 并**显式报错**，不静默走方形路径）；② J1=3.0347 Å 与 J2=3.1980 Å 只差 0.1632 Å ⇒ `magnetic_shells(tol=0.2)` 会**并成一个壳且不报错**（必须 `tol=0.05`）；③ C 心键连的是胞内两个 basis ⇒ basis 符号对 (n0,n1)-奇偶序不可见，母序与其 basis 翻倍版同选必使 det=0（脚本硬拒）。**两把尺子**：J1 的 3D 键长 3.0347 Å（含 1.71 Å 起伏）与 `poscar_to_mc3.py` 的 2D 投影 2.5064 Å 都对，但不可混用 ⇒ driver 用跨口径不变量 `per_site_bonds`(4/2/2) 校验 J 绑定。另外 **CrOCl 的基态是弱铁磁倾斜（需 SOC/DMI）**，本仓库共线 Heisenberg 层不覆盖，倾斜留给 P3/P2′；CrOCl 的 **J 与 T_N 均未核实**（`literature_J.json` 该条 J 全 null）⇒ `mc`/`report` 必须显式给 `--J1/--J2/--J3`。
1e. **六角 kagome（Fe₃Sn 类）目前只到系数层（P4.1），没有物理结论**：三壳 J1(a/2)/J2(a√3/2)/J3(a) 的系数矩阵、可解性与硬门禁已就位（`build_kagome_cases.py` + `solve_J_kagome.py`），但**没有任何真实非共线 DFT 能量**。四条硬约束：① q=0 一族（`fm/q0_120/q0_ud/q0_duu`）满足 c1==c2 且 c3 相同 ⇒ 子块 rank=2（cond≈2.3e34），J1/J2 在那里**不可分离**，必须用含 √3×√3 的组合（默认 `fm/q0_120/sq3_updntri/sq3_alt`，rank=4，cond=98.51）；② 超胞 n0/n1 必须是 3 的倍数（√3×√3 相位周期 3），配位数必须是 (4,4,6)，否则 builder 硬报错；③ 3 壳模型把 **J1 与 J1′（同距离、不同亚晶格对）合并成一条**，`J_result.json` 的 `J1` 是两者平均 —— 不许说成「解出了 J1/J1′」；④ `--model fe3sn` 的 S=1.0 / moment=2.0 / a=5.4 Å 是理想化占位（Fe₃Sn 文献 J 未核实，预设 `literature=None`）⇒ `mc`/`report` 必须显式给 `--J1/--J2/--J3`。真 J 与 T_N 属 **P4.2**（需真实非共线 DFT）。
1f. **DMI / 磁性拓扑绝缘体（MTI）只做判别，不进计算主链（P6）**：判据、拒绝模板与可行性备忘在 `DMI_MTI_MEMO.md`（同批更新 `dsh/2d-magnet/SKILL.md` 的判别节与 smoke §24）。本仓库**不提取 D**（需 SOC 四态法/自旋螺旋法 + 3 分量 D 矢量）、**不做能带拓扑**（Z₂/表面态），MC 侧无外场、无拓扑荷 ⇒ 报不出 D 分量、`|D|/J`、skyrmion 半径、拓扑荷 `Q`、B–T 相图、层间 `J_c`、Z₂。**唯一受控例外**：把**文献 D 值**换成逐壳 3×3 反称矩阵喂给 MC —— `python scripts\dmi_matrix.py --D … --J … --S … --out shell.json` 再 `poscar_to_mc3.py … --shell-def shell.json`（两个脚本的 `--S` 必须一致；别用 shell 的 `>` 写 JSON，PowerShell 会写成 UTF-16LE）。这是「拿别人的 D 做格点模拟」，**不许说成「本流水线算了 DMI」**，且只能看零场基态能。
1g. **交变磁体（altermagnet）只做几何判别 + 标量序参量（P5），没有能带结论**：`scripts/altermagnet_check.py` 只判两套磁性亚晶格的连接操作 —— 平移 ⇒ 常规共线反铁磁；反演（且两套亚晶格共面）⇒ 宇称禁止自旋劈裂；旋转/镜面 ⇒ **几何候选**。判据到此为止：是否真有 k 分辨自旋劈裂**要能带计算**，本仓库不做能带 ⇒ 不许报「这是交变磁体 / 劈裂多大」。MC 侧 `--order sublattice:<mask>` 给的是**标量 |L|** 的 Binder U4 交点，不分辨两套亚晶格的 Néel 矢量方向。判据容差以 **Å** 计（不是分数坐标元组比较：抖动 3e-6 的坐标仍能认出对称操作）。真 fixture：MnPS3 几何候选（镜面连接，残差 0.012 mÅ）、CrOCl 宇称禁止（反演连接，与 1e/局限 9 的弱铁磁倾斜一致）、Fe3Sn 不适用（3 套亚晶格）。
1h. **CrOCl 弱铁磁倾斜（P2′）只做到系数层 + 反解器，没有 θ 的自洽解**：`solve_J_crect_cant.py` 与 `tests/test_crect_cant.py`（44 项）把倾斜家族的系数与反解跑通，但 ① **θ 是输入参数**（闭式给的是「给定 θ 时系数是多少」，不是「基态 θ 是多少」）；② **一次倾斜 SCF 定不出唯一 (J, D)** —— 单壳闭式只有一条冗余自由度 ⇒ `--orders cant` 取「倾斜家族 + 张成所需的 fm/母序行」，`--audit-only` 打印自由度；③ 倾斜序属 `mag_vectors.SOC_ORDERS` ⇒ 装配**必须** `--soc`（否则硬拒）；④ CrOCl 文献 J 未核实 ⇒ 真值路径必须 `--ref`/`--synthetic`；⑤ **没有任何真实非共线/SOC 能量**，验收只到合成能量闭环 + 与独立实现逐位对拍；⑥ **driver 只做透传 + 改派**（`build --cant-theta/--cant-base/--cant-azimuth` → builder；`solve` 见到 `cant_orders` 自动派 `solve_J_crect_cant.py --orders cant`，三条守卫硬拒：非 C 心矩形给 `--cant-theta`、`--cant-base fm`、有倾斜序却没 `--soc`）——**路由不等于物理解**，θ 仍由用户给定。基态 θ 与 D 仍属未覆盖。
1i. **kagome 的物理层（P4.2）只有手册 + 装配侧工具，没有任何真实 J / T_N**：`P4_2_KAGOME_RUNBOOK.md` 规定入口条件、装配/回收/判定与偏差；回收通道与外部机时交接**已自动化**（`scripts/ingest_dft.py` 的 11 条硬闸门、`scripts/dft_handoff.py` 的 pack/ingest —— 此前「没有从 `03_dft/runs/` 的自动 ingest」的说法已不成立）；k 点收敛抽查也已有工具（`run_2d_magnet.py kconv --outdir <proj> --k-grids 3,4,6,8`，或直调 `scripts/k_convergence.py`），但它**不代提交作业**、也不判「J 是否收敛」——总能量对 k 网格收敛只是 J 可信的**必要**条件（J 是能量**差**），产物 `kconv/k_convergence.json` 里带着这条免责声明。路线图上 P4.2 仍未勾选，`literature_J.json` 里 Fe3Sn **条目已有但 J1/J2/J3 仍为 `null`**（按设计：原文自承交换常数未确立 ⇒ `--ref Fe3Sn` 会被明确拒绝）；手册 §8 的机时是**纯估算**（必须先跑一个 FM 参考点并用实测值替换，否则不许当预算），§9 明确把「J1/J1′ 分离」列为**做不到**。
2. **J 取文献/给定**，未跑 DFT 前非本项目计算值。
3. **MAE 默认 0** → 二维无真正相变（Mermin-Wagner），临界温度为有限尺寸伪相变估计。
4. **`mp_fetch_2d.py` 的 `exfoliate` 有已知缺口**：菱面体层状（层法线沿 [111]）和 MP 已存成单层紧胞会失效；`reorient.py` 修了「层未对齐 + 真空不足」，但多层体相仍需人工剥单层。
5. **k 网格 3×3×1、U 值**为默认/可调假设。
6. **zigzag/stripy 命名**可能和文献互换（M 点两成员标注约定），只影响名字不影响反解。
7. **MC 峰判据必须落在温度窗口内部**才算分辨：峰落在首/末格点 = 窗口截断或网格跨过窄峰，
   该边界值**不是** T_C/T_N。driver 检测到会自动缩放窗口 + 加密退火重扫
   （`--mc-refine-nt` 默认 60 / `--mc-refine-max` 默认 3 / `--no-mc-refine` 关闭），
   首轮原始扫描保留为 `mc_mt_survey.csv/.png`；重扫后仍落边界则结论显式标「未分辨」。
   典型例：MnI2 默认窗口 5–150 K、nt=30（格距 5 K）整个跨过 ~3 K 宽的 Cv 峰，
   曾把 T_N ≈ 7.3 K 报成窗口下限 5.0 K。

---

## 8. Agent 接手检查清单

收到「去 MP 找 XX 结构算磁耦合 + MC」：

1. **先分类磁性格点**：蜂窝（原胞 2 磁原子）/ 三角（1 磁原子）/ 矩形（正交 **a≠b**，如 CrSBr）/ 方形（正交 **a=b**，如 La₂CuO₄、K₂NiF₄）/ **C 心矩形（正交 a≠b 且原胞 2 个磁原子，如 CrOCl）**/ **六角 kagome（六角 a=b 且原胞 3 个磁原子，最近邻 a/2 上 4 配位，如 Fe₃Sn）**/ 簇磁体（磁性原子成簇）；确认价态（S、μB；簇磁体用簇净自旋）。**先数「胞内几个磁性原子」再量 a/b**：1 个磁原子 → a≠b 走 `build_rectangular_cases.py`、a=b 走 `build_square_cases.py`；**2 个磁原子且 a≠b** → 走 `build_crect_cases.py`（不是 `build_rectangular_cases.py`，后者默认 1 磁原子/胞，`--selftest` 理想值表会 FAIL）；**3 个磁原子**则看最近邻配位数：**4** → `build_kagome_cases.py`（kagome），**6** → `build_fgt_cases.py`（双亚晶格铁磁网络，**不是** kagome —— 两者同为「六角 + 3 磁原子/胞」，只靠配位数区分）。`detect_lattice` 已内置这条判据（在并簇之前判），所以 `--model fe3sn` 不会再落到 `cluster`。
2. 查 `literature_J.json`；没有就联网查 J（**转号！**），查不到就问用户或走 DFT。簇磁体需簇间 J。
3. 跑 `run_2d_magnet.py --formula XX --literature ...`（或手动参数）；三角/受挫晶格用 `--order cv`。
4. 读 MC 输出 T_C/T_N（峰判据须落在窗口**内部**；落边界会被自动重扫，仍落边界则标「未分辨」），
   对照实验体相量级（单层偏低）。
5. 在 `workflow.md` 声明：J 来源、MAE、有限尺寸、单层 vs 体相、U/k 假设。
6. **交付 py 文件（ship-mc）**：走 `run_2d_magnet.py` 时由 driver 自动完成（`05_mc/` 内
   `mc_2d.py` + `poscar_to_mc3.py` + `write_hamiltonian.py`，按分支再加枚举+反解脚本，见第 1 节产物清单）。
   **非标准分支（多亚晶格网络，如 Fe3GeTe2）不走 driver ⇒ ship-mc 不会发生，必须手动补齐**
   （首次交付极易漏掉，只剩数据+图而无法脱离主仓重跑），并跑一次
   `python 05_mc/run_rest.py "<05_mc 绝对路径>" 8` 验证自包含：已完成的 run 全 `skip`、退出码 0、
   **目录文件数不变**。参考实例 `Fe3GeTe2_pipeline/05_mc/`（15 个 py）。
