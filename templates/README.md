# templates/ —— 扩展本流水线时照抄的骨架

这个目录**不参与运行**（没有任何 `scripts/` 代码读它），它是给「下一个加材料 / 加晶格分支的人（或 agent）」的对照件。
三条扩展路径按改动面从小到大：

| 你要做的事 | 照抄哪个模板 | 主要改动面 | 工时量级 |
| --- | --- | --- | --- |
| 只加一个材料（晶格类型已支持） | `material_preset.json` | `run_2d_magnet.py` 的 `MODEL_PRESETS` + `scripts/literature_J.json` | 分钟级 |
| 加一种新磁序（格子不变） | `material_preset.json` 的 `orders` 段说明 | `scripts/magnetic_orders.py` 的 `ORDER_RULES` / `SPECS` | 小时级 |
| 加一个**新晶格分支** | `new_lattice_checklist.md` | 规格层 + 构建器 + 反解器 + driver 路由 + 渲染 + 测试 + 四份文档 | 天级 |

## `material_preset.json` 怎么用

它对应 `run_2d_magnet.py` 里 `MODEL_PRESETS` 的**一个条目**，字段逐一对齐（`formula` / `element` / `nonmag` /
`lattice` / `S` / `moment` / `U` / `literature` / `desc`）。粘进 `MODEL_PRESETS` 时注意：

- 值必须与 `scripts/literature_J.json` 的同名条目**自洽**（`element` / `nonmag` / `S` / `moment` / `U`），
  否则 `check_presets()` 在 `fetch`/`build` 阶段就会报 `preset_mismatch`（预设是材料专属值，不跟随文献表漂移）。
- `desc` 里**必须**写清 J 是否已核实：写「文献 J 待核实 / 未核实」的材料，跑 `mc`/`report` 时必须显式给
  `--J1/--J2/--J3`，脚本不会替你猜。
- `literature` 只能填 `literature_J.json` 里真实存在的键（或整条省略）。
- JSON 没有注释语法 ⇒ 字段说明写在下面这张表里，模板文件本身保持是**可解析的纯 JSON**。

| 字段 | 含义 | 约束 |
| --- | --- | --- |
| `formula` | 化学式，决定默认输出目录 `<formula>_pipeline` | 与 `literature_J.json` 的键同名时最省事 |
| `element` | 磁性元素（唯一） | 必须在 POSCAR 里存在；MAGMOM 只给这个元素 |
| `nonmag` | 非磁元素列表 | 与 POSCAR 的成分**互斥且完备**：漏写或写多都在构建期硬报错 |
| `lattice` | 晶格分支名 | 必须是 `detect_lattice()` 能返回的名字之一（`honeycomb` / `triangular` / `square` / `rectangular` / `rect_centered` / `cluster`） |
| `S` | 自旋量子数 | 元素表兜底之外的材料专属值 |
| `moment` | 初始磁矩 (μB) | 一般 = 2S（`--moment` 会写进 INCAR 的 MAGMOM 模长） |
| `U` | DFT+U 的 U (eV) | 材料专属（如 MnPS₃ 的 5.0 取自 Sivadas 2015，而非元素表 Mn 的 4.0） |
| `literature` | 文献锚点键 | 可省略；给了就必须存在于 `literature_J.json` |
| `desc` | 一句话说明 + **J 核实状态** | 与文献表不一致会被预设自检抓出 |

## 相关的两份权威文档

- `WORKFLOW.md`：阶段 → 脚本 → 产物的对照表与「继承表」（哪个脚本复用了哪一层的逻辑）。
- `skill/2d-magnet/SKILL.md`：给 agent 的入口说明与自我学习协议（加东西时先看这份的约束）。
