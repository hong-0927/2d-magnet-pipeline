# 新晶格分支 checklist

以 P1（方形 `square`）/ P2（C 心矩形 `rect_centered`）/ P3（非共线 + SOC 地基）/ P4.1（六角 kagome 系数层）四次实战为准。**每一步都有对应的验证闸门**，
不要等最后一起测——规格层没锚定就写构建器，会让误差在下游放大且难以定位。

> 总原则（`NEW_TYPES_PLAN.md` §104–117 + `dsh/2d-magnet/SKILL.md`）：**数据驱动 + 薄适配层**。新的晶格差异应该落在
> 规格层（`SPECS` / `ORDER_RULES`）与一个小构建器里，**不要复制粘贴**一个平行的 builder / solver / renderer。
> 前端（driver）只做路由，不承载物理。

## 0. 先写规格层（`scripts/magnetic_orders.py`）

1. `SPECS["<lattice>_<material>"]`：壳层（按距离升序的 `shells` + 每壳的 J）、磁序集合、`selected_orders`、
   `spin`、以及本材料的 fixture 路径。**J 未核实就写 `null` 并在 `desc` 里说明**，不要填一个"看起来合理"的数。
2. `ORDER_RULES`：新磁序的名字 → 生成规则的声明式描述（如 P3 的 `"120deg": "color=(n0+n1)%3; theta=color*120deg"`）。
   自旋矢量本身放 `scripts/mag_vectors.py`（共线 / 面内 120° / ±θ 倾斜 / 螺旋），**不要在 builder 里手写角度**。
3. `detect_lattice()`：新晶格的几何判据 + **与已有分支的区分性**（P2 的教训：CrOCl 两个 Cr 只差 3.03 Å，
   旧判据会误判成 `cluster` ⇒ 新增 `rect_centered`，并在方形网并簇时返回 `sq_centered` 显式报错而不是静默降级）。
4. 需要平移搜索的晶格必须显式打开（P1 的教训：方形**必须**含周期 4 的 `p4_a`，否则 J3 被低估 2 倍或与 E0 不可分离）。

**闸门**：`python scripts/magnetic_orders.py --regress` 全 PASS（每条链路逐位回归）。

## 1. 构建器（`scripts/build_<lattice>_cases.py`）

照抄最接近的现有构建器（`build_square_cases.py` / `build_crect_cases.py` / `build_kagome_cases.py` / `build_magnetic_cases_tri.py`），但只保留晶格差异：

- 每个磁序一个目录：`POSCAR` + `magmom.txt`（磁原子在前、非磁元素 `0.0` 或 `0.0 0.0 0.0`）+ 胞级声明。
- 写 `magnetic_orders.json`：`selected_orders` / `n_mag` / `n_atoms` / `moment` / `S2` / `coefficients_per_cell` /
  `inversion_matrix_per_cell`（后两者是反解器的输入，**必须来自规格层而不是手抄**）。
- 自检 `--selftest`：用一组**已知的** J 正演出能量再反解回来，容差写死在脚本里（P1/P2 是 0.02）。

**闸门**：`python scripts/build_<lattice>_cases.py tests/fixtures/<FIXTURE>_POSCAR --element X --nonmag ... --outdir ../_selfcheck/<x> --selftest` → `ALL PASS`。

## 2. 反解器（`scripts/solve_J_<lattice>.py`）

- 能量口径**只走** `scripts/dft_energy.py` 的 `read_toten()`（`free energy TOTEN`，取最后一次出现；缺则报错）。
  不要在本脚本里再写一份正则——P3 已把六个 `solve_J*` 全部改成转发。
- 产物**必须**用 `scripts/j_result.py` 的 `build(...)` + `write(...)` 落 `J_result.json`（含 trust / residual / cond 闸门），
  便于 `report` 阶段汇总与 `poscar_to_mc3 --jfile` 消费。
- `--demo`：用已知 J 正演再反解的自检闭环（应 PASS 并打印残差）。

**闸门**：`... --demo` PASS；且 `validate(J_result.json)` 的 `strict=True` 通过。

## 3. driver 路由（`run_2d_magnet.py`）

- `--lattice` 的 choices、`detect_lattice()` 返回值 → 构建器 / 反解器 / 渲染器的**分发表**各加一列。
- 若是三角那样的既有路径，**不要新增分支**：把它退化成薄适配层转调统一入口（P3 的 `assemble_dft_tri.py` 只做 `subprocess.call`）。
- 非共线 / SOC 相关：`--noncollinear` / `--soc` / `--saxis` 的转发与 `child_env()` 已就位，直接用。
- workflow 文本里的阶段行不要写死（用 `_dft_line` / `_dft_line_tri` 这类变量拼），否则文档与行为会漂移。

## 4. 渲染（figures）

- 几何/配色逻辑复用 `render_supercell.py`；**每条渲染路径都要有 matplotlib 与纯 Pillow 两条等价实现**
  （Pillow 版的图表原语在 `scripts/pil_fig.py`），driver 用 `have_matplotlib()` 现场二选一。
- 出图后**肉眼看一眼**（P3 的教训：三角 quiver 在 |S|=5 μB 时画出 5.6 Å 长箭头互相纠缠，数值全对但图不可读）。

## 5. 测试与回归

- `tests/test_<lattice>_cases.py`：**三方交叉锚定**——构建器产物 / 规格层 / 手推理想值，三者逐位对齐
  （P1 方形 28 项、P2 C 心矩形含 basis 退化不可见性、det 奇异性、`tol=0.2` 并壳陷阱）。
- `tests/smoke_test.py` 加一节 `§N`（主回归入口）：至少覆盖路由、产物结构、以及该分支特有的硬报错（§115 不静默降级）。
- fixture 放 `tests/fixtures/`（不要引用本机绝对路径）；非共线快照回归参照
  `tests/fixtures/gdi2_nc_soc_120deg.INCAR.norm` 的「正文逐字 + MAGMOM 归一成描述符」写法。
- 若是**新增错误码**：`run_2d_magnet.py` 的 `ERROR_CATALOG` 加一条（`static` 区分「环境自检可命中」与「运行时才知道」），
  并在 `KNOWN_ISSUES.md` 对应那张表里加同一行（两张表是一一对应的）。

**闸门**：`python tests/smoke_test.py` + `test_magnetic_orders.py` + 新增的 `test_<lattice>_cases.py` 全绿；
`python run_2d_magnet.py doctor` 只允许出现环境类 WARN（缺依赖等），不允许出现逻辑 WARN。

## 6. 文档（与代码同批更新，§117）

- `README.md`：CLI 表、已知局限（新晶格的静默陷阱要写进去）、路线图打勾。
- `WORKFLOW.md`：阶段 → 脚本 → 产物表 + **继承表**（说明新脚本复用了哪一层）+ 目录树。
- `KNOWN_ISSUES.md`：新错误码两表各一行。
- `NEW_TYPES_PLAN.md`：该项状态改 ✅ 并记录**偏离原计划的地方**（这类偏差记录比打勾本身更有价值）。
- `dsh/2d-magnet/SKILL.md`：agent 的入口说明要不要新增一段分支用法。
