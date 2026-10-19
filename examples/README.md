# 样例产物（不需要 VASP）

这两组产物是**快照**，用来回答一个问题：**这条流水线跑完到底会得到什么文件、长什么样。**
它们全部由仓库自带夹具生成，不需要 VASP、不需要联网、不需要额外依赖（缺 matplotlib 时渲染走 Pillow 兜底）。

| 目录 | 怎么生成的 | J 从哪来 | 能不能当物理结论 |
|---|---|---|---|
| `la2cuo4_square/` | `run_2d_magnet.py --model la2cuo4 … --skip-dft`：结构 → 系数 → MC → 图，走完整链路 | **文献条目**（La₂CuO₄：J1=-111.8 / J2=+11.4 / J3=-2.0 meV） | T_C 是本仓 MC 算出来的，但 J 是文献值、MC 尺寸只有 12×12、判据是单次 χ 峰 ⇒ 只当**量级 sanity check** |
| `kagome_fe3sn_dryrun/` | `scripts/dry_run_kagome.py --synthetic 5 -0.5 0.2 --mc`：合成 OUTCAR 空跑整链 | **合成值**（脚本自己造的） | **完全不能**。它证明的是「管道通」（装配 → 硬闸门 → 回收 → 反解 → 凭证），数字没有物理意义 |

每个目录里的 `README.md` 写了复现命令与怎么读那些文件。两组快照的完整工程目录另含
`workflow.md`（逐步阶段报告）与 `02_supercell/cases/<序>/{POSCAR,magmom.txt}`（送 DFT 的输入），
为节省仓库体积没有快照进来。

> 路径字段说明：快照是在临时目录跑出来的，其中的绝对路径（`magnetic_orders.json` 的 `input`、
> `J_result*.json` 的 `source_runs` / `manifest`、`ingest_report.md` 与 `hamiltonian.md` 的来源行）
> 已全部改写成工程内相对路径，`tests/smoke_test.py` §28 有一道闸门禁止再把它带回来。

## 这个流水线会产出什么（有 / 没有 / 为什么）

| 产物 | 谁写 | 落在哪 | 说明 |
|---|---|---|---|
| 结构 `POSCAR`（含重取向） | `reorient.py` / fetch | `01_structure/` | 真空层、取向，注释行写清改了什么 |
| 每个磁序的 `POSCAR` + `magmom.txt` | 各晶格 builder | `02_supercell/cases/<序>/` | **送 DFT 的输入**；共线 1 分量、非共线 3 分量（`mag_vectors.render_magmom`） |
| `magnetic_orders.json`（规格层） | 各 builder | `02_supercell/cases/` | 壳、磁序、每序系数、配位数、可解性（rank/cond） |
| `J_result.json` | 各解算器 | `02_supercell/cases/` | 统一产物：J、cond、residual、`source`/`confidence`（**信任度由代码判定，不能手改**） |
| `INCAR` / `KPOINTS` / `POTCAR.spec` | `assemble_dft.py` | `03_dft/`（模板）+ **`runs/<序>/`（可直接投的整套）** | `--kpoints/--npar/--isym` 可覆盖；**POTCAR 本体不生成**（版权在 VASP 发行方），只给 spec 与拼装顺序 |
| 回收凭证 | `ingest_dft.py` | `cases/ingest_manifest.json`、`cases/ingest_report.md` | 11 条硬闸门全过才写；含 sha256 与 provenance（runs 目录、口径） |
| 交接包 | `dft_handoff.py` | 你指定的 tar 目录 | 打包 `02_supercell/cases` + `03_dft`，附作业脚本；不代提交、不含 POTCAR |
| 超胞磁序图 | `render_*.py` | `figures/supercell_orders.png` | 每个磁序一幅子图；matplotlib 或 Pillow 兜底二选一 |
| 哈密顿量说明 | `write_hamiltonian.py` | `hamiltonian.md` | 键、壳、符号约定、每格点键数核对，供人工复核 |
| MC 结果 | `mc_2d.py` / `mc_fss.py` | `05_mc/mc_mt.csv`、`mc_mt.png`（FSS 另有产物） | T、序参量、χ、Binder U4；峰判据落窗口边界会自动缩放重扫 |
| 阶段报告 | driver | `workflow.md` | 每一步的真实命令与参数（便于复现/排查） |
| k 点收敛抽查 | `k_convergence.py` | `kconv/kconv_plan.json`、`kconv/k_convergence.json` | P4.2：只换 KPOINTS，回收后判定「多密才够」 |

**刻意不产出**的东西：POTCAR 本体、任何 DFT 计算结果（本仓不代跑作业、不含机时）、
k 分辨自旋劈裂 / 能带 / Z₂（超出契约，见 `../LIMITATIONS.md` 14、16）、
DMI 的 D 值（只做识别层与系数层，见 `../DMI_MTI_MEMO.md`）。
