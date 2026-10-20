# 例子 1：方形晶格 —— La₂CuO₄ 单层（文献 J + 本仓 MC）

这个目录是**真实跑出来的快照**（不是手写示意），演示「文献 J → MC → T_C」这条路在方形晶格上的全部产物。
J 取文献值（`literature_J.json` 的 `La2CuO4` 条目，Coldea et al., PRL 86, 5377 (2001)：295 K 自由拟合
J=111.8±4 meV、J'=-11.4±3 meV，已按本仓约定转号成 J1=-111.8 / J2=+11.4 / J3=-2.0 meV）。

## 复现命令（约 1–2 分钟，不需要 VASP）

```powershell
cd 2d_magnet_pipeline
$ex = "$env:TEMP\la2_example"

# 1) 结构 → 系数 → MC → 图（跳过 dft 阶段时加 --skip-dft）
python run_2d_magnet.py --model la2cuo4 --poscar tests/fixtures/La2CuO4_POSCAR `
    --outdir $ex --mc-size 12 12 --mc-nt 8 --mc-tmax 300 --mc-seeds 2

# 2) 反解器产物（内部合成能量, 只为展示 J_result.json 的 schema 与自检闭环）
python scripts/solve_J_square.py $ex\02_supercell\cases --demo
```

## 目录里有什么

| 路径 | 是什么 |
|---|---|
| `01_structure/POSCAR` | 重取向 + 15 Å 真空的单层胞（Cu 1/胞） |
| `02_supercell/cases/magnetic_orders.json` | **规格层**：壳 `J1 3.7992 Å / J2 5.3729 Å / J3 7.5984 Å`、6 个磁序、每序系数、配位数 |
| `02_supercell/cases/J_result.json` | 统一产物：`source=demo`、`confidence=low`、`cond(A)=40.41`、`residual=1.67e-16 eV`（**自检用，不是 DFT 结果**） |
| `02_supercell/cases/fm/{POSCAR,magmom.txt}` | 送 DFT 的一个磁序示例（共线 ⇒ 1 分量 MAGMOM） |
| `03_dft/{INCAR.static,INCAR.relax,KPOINTS,POTCAR.spec}` | 装配出的 VASP 输入模板（`ENCUT=520`、`EDIFF=1E-6`、Cu 的 `LDAU` U=7.0；KPOINTS 默认 3×3×1；**POTCAR 本体不生成**，只给 spec） |
| `runs/fm/{INCAR,KPOINTS,POSCAR}` | 可直接投递的整套（driver 会为每个磁序各写一份） |
| `05_mc/mc_mt.csv`、`mc_mt.png` | MC 结果：T、序参量、Susceptibility、Binder U4（61 行） |
| `figures/supercell_orders.png` | 4×4 超胞的磁序图（fm / stripe_a / p4_a / mix_ab） |
| `hamiltonian.md` | 哈密顿量与符号约定卡片（供人工复核） |

## 结果怎么读（以及不能怎么读）

- MC（12×12、2 seeds、60 个温度点、窗口 131–1200 K）：**χ 峰 ≈ 294 K**，`|dm/dT|` 峰 ≈ 276 K，
  Binder U4 从 0.397（低温有序端）单调降到 ~0.17。峰落在初始窗口上限，driver 自动缩放重扫了一轮
  —— 这是判据（峰必须在窗口内）在工作，不是异常。
- **能说的**：用文献 J 跑出的 T_C 与 La₂CuO₄ 的量级一致（三维 T_N = 325 K，Coldea 2001）。
  这是一次**量级 sanity check**：系数层、MC 引擎、峰判据都通了。
- **不能说的**：这不是「La₂CuO₄ 单层的 T_C 计算值」。理由三条 ——(1) J 是文献值不是本仓算的；
  (2) 12×12 有限尺寸 + 单尺寸、未做 FSS，峰位系统性偏高；(3) `J_result.json` 的
  `source=demo` ⇒ schema 自带规则把它标成 `confidence=low`，不许当结论引用。
- 想要真实数字：走 `run_2d_magnet.py dft` → 上机 → `ingest` 回收 → `solve`（见 `../KAGOME_RUNBOOK.md`
  的同款流程；方形晶格的装配/回收/反解命令与 kagome 完全同形）。
