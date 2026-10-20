# 例子 2：kagome —— Fe₃Sn 的**合成 dry-run**（只证明管道通）

kagome 这条链在本机**没有真实机时**（无 VASP），所以这里给的是「用合成 OUTCAR 空跑整链」的快照：
装配 → 硬闸门 → 回收 → 反解 → 凭证，每一步都真的跑了，但**能量是脚本自己造的**。
它的作用是让你在拿到机时之前，先看清「回收之后会得到哪些文件、长什么样」。

## 复现命令（约 1–2 分钟，不需要 VASP）

```powershell
cd 2d_magnet_pipeline
$ex = "$env:TEMP\kag_example"
New-Item -ItemType Directory -Force -Path $ex\02_supercell\cases | Out-Null

# 1) 构建 √3×√3 构型集（selected 4 个: fm / q0_ud / sq3_updntri / sq3_alt，全共线）
python scripts/build_kagome_cases.py tests/fixtures/Fe3Sn_POSCAR --supercell 3 3 1 `
    --outdir $ex\02_supercell\cases --element Fe --nonmag Sn --moment 2.0 --spin 1.0

# 2) 合成 OUTCAR 空跑整链（--synthetic 是"真值", 只为验证能否反解回来）
python scripts/dry_run_kagome.py $ex --synthetic 5 -0.5 0.2 --mc
```

## 目录里有什么

| 路径 | 是什么 |
|---|---|
| `02_supercell/cases/magnetic_orders.json` | 规格层：3×3×1 超胞（27 个磁位点）、5 个共线磁序（`fm/q0_ud/q0_duu/sq3_updntri/sq3_alt`）、per-site 键数 `(4,4,6)`、`rank(A)=4/4`、`cond(A)≈101` |
| `02_supercell/cases/J_result.dryrun.json` | 反解产物：`source=demo` ⇒ `confidence=low`；`A` 矩阵、`J_meV`、`residual` 都在里面 |
| `02_supercell/cases/ingest_report.md` | **回收报告**：逐序能量、`ionic steps`、SCF 是否收敛、与 fm 之差、sha256 前 12 位、警告清单 |
| `02_supercell/cases/fm/{POSCAR,magmom.txt}` | 送 DFT 的一个磁序（共线标量 MAGMOM，`ISPIN=2`；不需要非共线） |
| `02_supercell/cases/mc_dryrun_mt.csv`、`mc_dryrun_mt.png` | 用**合成** J 跑的一次 MC（曲线形状可看，数值无意义） |
| `figures/supercell_orders.png` | 3×3 超胞的磁序图（共线 ↑/↓ 着色） |

## 结果怎么读（以及不能怎么读）

- **能说的**：管道闭合了。合成能量 → 闸门（体积/收敛/结构/口径）→ `ingest_manifest.json` 凭证
  → 反解 → 往返误差 **0 meV**（合成 J 被原样解回来）→ MC 出曲线。`ingest_report.md` 里
  `LNONCOLLINEAR/LSORBIT/ENCUT/NKPTS/VASP 版本` 都被逐序核对过。
- **不能说的**：任何 J、T_C、T_N 数值。`J_result.dryrun.json` 的 `source=demo` 就是标记，
  MC 那两行也明确写着「伪临界点」。
- 真实机时的下一步：`../P4_2_KAGOME_RUNBOOK.md`。要点是先把 §1 的参数表钉死（化合物 / a / S / U /
  POTCAR / k 网格 / ENCUT / 超胞），再按 §3–§5 装配、提交、回收，最后用 §6–§7 反解并对照 M4 锚点。
- 另有两个与回收有关的工具，这里没跑但值得知道：`scripts/dft_handoff.py`（把 cases + DFT 输入打成
  交接包并附作业脚本，不代提交）与 `scripts/k_convergence.py`（只换 KPOINTS 生成 k 收敛变体，
  回收后判定「多密才够」）。
