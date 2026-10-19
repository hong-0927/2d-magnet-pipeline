# Known Issues / Error Catalog

运行 `python run_2d_magnet.py doctor` 会做环境自检，并把**命中的已知问题**映射成「症状 → 修复」。
本文件是同一份目录的静态版，覆盖两类条目：**环境可探**（doctor 能当场命中）和**运行时**（跑到某阶段才出现）。

## 环境自检可命中的条目

| key | 症状 | 修复 |
|---|---|---|
| `py39` | 语法/类型报错，或 numpy/numba 版本不兼容 | 升级到 Python ≥ 3.9（如 `py -3.12` 或 conda） |
| `missing_dep` | `ImportError: No module named 'xxx'` | `pip install -r requirements.txt` |
| `packaging` | render 阶段 `ModuleNotFoundError: No module named 'packaging'`（matplotlib 顶层 import 连带需要它） | `pip install packaging` |
| `fig_backend` | `doctor` 报「缺少依赖 matplotlib」或「缺少 Pillow」；render 阶段 `[warn] 画图跳过: No module named 'pyparsing'`，`figures/` 空着没 PNG | 通用超胞图有两条等价路径：matplotlib，或纯 Pillow 的 `scripts/render_supercell_pil.py`（driver 用 `have_matplotlib()` 自动二选一）⇒ 装 `Pillow` 即可出图。要修 matplotlib 需补齐 `pyparsing`/`cycler`/`contourpy`/`fontTools`/`kiwisolver`（matplotlib 在 import 期就要）。簇磁体 `render_cluster_orders.py`、三角 `render_tri_orders.py` 与 MC 曲线图（`mc_2d.py` 内）也都各配了纯 Pillow 等价脚本 `render_cluster_orders_pil.py` / `render_tri_orders_pil.py` / `_plot_mt_png`，driver 与 `mc_2d.py` 同样自动二选一 ⇒ 装了 `Pillow` 这几条渲染路径全部能出图。**另两条也一样**：有限尺寸标度图（`mc_fss.py` 的 `_plot_fss_png`）与 FGT 双亚晶格磁序图（`build_fgt_cases.py` 的 `_plot_fgt_png`）各有纯 Pillow 实现并自动二选一 ⇒ **六条渲染路径全部与 matplotlib 解耦**，共用原语 `scripts/pil_fig.py`（磁序图用点/箭头层，曲线图用 `chart_*` 直角坐标系层） |
| `scripts_missing` | driver 报缺脚本 / 子阶段找不到脚本 | 确认 `scripts/` 目录完整（重新 `git clone` 或解压分发包） |
| `lit_missing` | `--literature <key>` 报找不到条目 / `literature_J.json` 缺失 | 确认 `scripts/literature_J.json` 存在且键名正确 |
| `preset_mismatch` | `doctor` 的 `--model` 预设自检报「预设 X 与文献条目 Y 的 `<字段>` 不一致」，或预设指向不存在的文献条 | 把 `MODEL_PRESETS[X]` 与 `literature_J.json[Y]` 的 `formula/element/nonmag/S/moment/U` 对齐（预设是材料专属值，不跟随文献表漂移）；回归见 `tests/test_crect_cases.py` |
| `mp_key` | fetch（`--formula` 下载）报 401/403 或空结果 | 设 `MP_API_KEY`；或 `--poscar` 用本地结构 / `--provider cod` 走公开库 |
| `cwd_ro` | `PermissionError`：建 `--outdir` / 写文件失败 | `cd` 到可写目录，或 `--outdir` 指到可写路径 |

## 运行时条目（跑到某阶段才出现）

| key | 症状 | 修复 |
|---|---|---|
| `static_incar` | `runs/<order>/INCAR` 无 `IBRION/NSW` → 静默假弛豫，OUTCAR 其实是静态总能 | 用 `03_dft/INCAR.relax` 做弛豫；或 relax_backfill 回填后用位移检查核对 |
| `atom_mismatch` | relax_backfill 报「原子数/物种不符」并拒绝 | 核对 CONTCAR 与 POSCAR.pristine 的物种与原子数一致后再回填 |
| `gbk` | 子进程输出 `UnicodeDecodeError`（GBK 编码） | 设 `PYTHONIOENCODING=utf-8` 再跑（driver 子进程已内置；手动跑脚本需自设） |
| `solve_no_outcar` | `solve_J3` 报「需要 fm/neel/zigzag/stripy 四态 OUTCAR」 | 先跑 DFT 生成各磁态 OUTCAR；或加 `--demo` 自检反解 |
| `fss_no_cross` | FSS 报「未分辨」/ `fss_fss.json` 里 `tc=null` | 加大 `--fss-sizes`（默认 12 18 24 30）或 `--mc-seeds` 压噪 |
| `peak_boundary` | T_C/T_N 标「未分辨」（峰落在温度窗口首/末格点） | 已内置自动缩放重扫；或调 `--mc-tmax` / `--mc-nt` 覆盖真峰 |
| `soc_required` | `assemble_dft.py` 报「磁序 `<order>` 需要 SOC 但没给 `--lsorbit`」/「`magmom.txt` 是 3 分量但当前是共线模式」/「`--lsorbit` 需要 `--noncollinear`」 | driver 侧加 `--noncollinear`（弱铁磁倾斜/DMI/kagome 再加 `--soc`，可配 `--saxis`）；手动跑脚本时给 `--noncollinear --lsorbit`。这是 P3 的硬门禁（§115 不静默降级）：缺 SOC 时倾斜/螺旋构型物理上不成立，不会按共线近似给你一个错的 J |
| `sq_centered` | driver 报「检测到「方形 + 多磁原子/胞」(sq_centered)」并退出 | 该网（a≈b 且胞内多个磁原子被 cutoff 并簇）本仓库未支持；**别硬套 `square`**（会丢掉两个 basis 的反平行自由度 ⇒ 错的 J）。要支持需新写方形版 basis 序 + 显式壳距 + 独立反解脚本 |
| `crect_singular` | `build_crect_cases.py` 报 `det<1e-6` / `cond>1e4` / 「母序与其 basis 翻倍版同选」 | 去掉退化组合（`afm_c_a` 等于 `afm_a`，`afm_c_diag` 等于 `afm_diag`），DFT 只提交 `fm/afm_c/afm_a/afm_b`；壳容差保持 `0.05`（`tol=0.2` 会把 J1=3.0347 与 J2=3.1980 并壳且不报错） |
| `mc_shell_mismatch` | driver 报「`[check]` MC 每格点键数 … != 规格层 `per_site_bonds`」并 SystemExit | J1/J2/J3 很可能被绑到了错的壳上（距离自动排序与规格层不一致）：核对 `build_*_cases.py` 的壳定义与 `poscar_to_mc3` 的口径（CrOCl 应 4/2/2，La₂CuO₄ 应 4/4/4） |
| `energy_convention` | `solve_*` 报「no 'free energy TOTEN' in `<OUTCAR>`」，或 `dft_energy.mixed_convention` 报「两种能量口径混用」 | P3 起所有 `solve_*` 只认 OUTCAR 的 `free energy TOTEN`（取最后一次）：缺它说明计算没跑完或输出被裁剪 ⇒ 补跑/换完整 OUTCAR（完全没有 OUTCAR 时才走 `E.txt` 兜底）。**不要用 `energy(sigma->0)` 凑数**——与 `free energy TOTEN` 混用会让 J 差一个常数项（阈值 1e-3 eV） |
| `kagome_degenerate` | `solve_J_kagome.py` 报「设计矩阵 rank=… < 4 ⇒ 这几个磁序**无法唯一确定** J1/J2/J3」，或「设计矩阵病态: cond=… (>1e4)」，或「`selected_orders=…` 与本脚本期望的 … 不一致」 | kagome 的 q=0 一族（`fm`/`q0_120`/`q0_ud`/`q0_duu`）满足 `c1 == c2` 且 `c3` 相同 ⇒ 子块 rank=2（cond≈2.3e34），J1/J2 在那里不可分离。必须用含 √3×√3 的组合（builder 的 `selected = fm/q0_120/sq3_updntri/sq3_alt`，rank=4，cond=98.51），且超胞 n0/n1 必须是 3 的倍数；换序后要重算 rank/cond——**「cond 没超阈」不等于「三个壳都分得开」**（P4.1 硬门禁，不静默出 J） |
| `dmi_needs_shell_def` | 想跑含 DMI 的格点模型，但用 `poscar_to_mc3 --jfile J_result.json`——**不报错**，反称项却全部丢失（模型退化成纯 Heisenberg） | DMI 只有 3×3 `J_mat` 装得下，而 `J_result.json` 只有标量 `J_meV` ⇒ 必须走 `--shell-def` 显式给 `J_mat`：`python scripts/dmi_matrix.py --D Dx Dy Dz --J … --S … --d … --out shell.json`（`J_D = -[D]×`，对应 `H_DM = +D·(S_i×S_j)`；矩阵已吸收 `S²`），再 `poscar_to_mc3.py … --shell-def shell.json`；两个脚本的 `--S` 必须一致。**D 必须来自文献或你自己的 DFT**（本仓库不提取 D，见 `DMI_MTI_MEMO.md` §3.2/§3.3）；引擎无外场、无拓扑荷 ⇒ 只能看零场基态能，不许报 skyrmion 相图 / 拓扑荷 |
| `sublattice_mask_invalid` | `mc_2d` / `mc_fss` / 驱动 `--mc-order` 报「`--order sublattice:<mask>` 需要掩码」/「掩码有空项」/「子晶格编号越界」/「掩码 `[...]` 覆盖了全部 N 个子晶格 ⇒ B 组为空」 | mask 是**逗号分隔的子晶格编号**（不写就是空掩码）：`sublattice:0` = 老的 `staggered`（双亚晶格体系逐位等价），`sublattice:0,1` 只在 `num_basis >= 3` 时合法（两套亚晶格模型给出它会把 B 组清空 ⇒ `mbx/cb` 除零，脚本硬报错而不静默）。编号是 `i % num_basis` 的 basis 序号（来自 POSCAR 原子顺序），不是超胞格点号；掩码**至少留一个**子晶格在 B 组。序参量只是标量 `|L|`，不分辨两套亚晶格的 Néel 矢量方向；几何判别用 `scripts/altermagnet_check.py`（P5，只判别不下结论，见局限 13） |
| `order_auto_ambiguous` | `mc_2d` 打 `>>> [warn] 候选基态近简并 (…)` 或 `>>> [warn] auto 判定的基态是 stripe_x … 已自动改用 --order cv`；`mc_fss` 直接报 `[error] auto 判定该模型 T=0 基态是 stripe_x (条纹反铁磁)` | `--order auto` 现在的判据是**T=0 候选基态能量比较**（`fm` / `neel` / `stripe_x` / `stripe_y`），不是看 `J1` 符号：条纹（或近简并）时「交错序参量的峰」不是 T_N。显式给 `--order`：铁磁 `total`、Néel `staggered`、逐子晶格 `sublattice:<mask>`、受挫/条纹看 `cv`（**仅 `mc_2d`**，`mc_fss` 用 U4 交点、没有比热路径 ⇒ 条纹基态在 FSS 里必须显式给 `sublattice` 或换 `mc_2d --order cv`） |
| `cant_theta_needs_crect` | 驱动报「给了 `--cant-theta` 但 `lattice=...`，不是 C 心矩形（`rect_centered`）」或「给了 `--cant-theta` 却用 `--cant-base fm`」 | 倾斜家族只在 C 心矩形网（CrOCl 类，2 磁原子/胞）上定义：`--cant-theta "8 16"` 配 `--model crocl`（或 `--lattice rect_centered`）。母序取 `afm_c`（默认）/`afm_a`/`afm_b` —— `fm` 被拒（全同向母序上的 ±θ 只是一次整体自旋旋转，系数行逐位等于 fm 行）。倾斜构型要 `--soc` 装配，缺它 `assemble_dft` 会**硬拒**（不静默按共线算）；`solve` 端按 `magnetic_orders.json` 的 `cant_orders` 自动改派 `solve_J_crect_cant.py --orders cant`。θ 是**输入参数**、不是自洽解出的基态角，φ 由几何定不出来（见 LIMITATIONS.md 11） |
| `states_unverified` | 驱动 `solve` 阶段被 `scripts/verify_states.py` 拦下（exit 1、不写 `J_result.json`）：报磁矩与目标 σ 符号/大小不一致，或「两个磁态坍缩成同一个态」，或缺某态的 OUTCAR | 这三类都是**能量映射的前提不成立**，不是脚本 bug：① 不一致 ⇒ VASP 从初始 MAGMOM 跑到了别的磁态，重跑该态（把初始 MAGMOM 放大量、或检查该序是否被结构对称性强制退化），**不要**加 `--allow-unverified-states` 硬看数；② 坍缩 ⇒ 两个磁序实际是同一个态（常见于反演/平移等价的组合），该方向的 J 不可分辨，换一个真正线性无关的磁序；③ 缺 OUTCAR ⇒ 先 `ingest` 回收。逐态证据在 `02_supercell/cases/states_verification.json`（`verdict` / `per_state` / `collapse_pairs`）。确需出数时用 `--allow-unverified-states`，但 `J_result.json` 会带 `states_verification.verdict = FAIL` + `note`，**不可当已验证结果引用**；非共线 / SOC 体系该闸门判 `SKIP`（见 LIMITATIONS.md 16） |
| `sampling_underresolved` | `mc_2d` 打 `>>> [warn] 峰处 T=… K 的有效独立样本只有 N < 100 (τ_int=… 步)` 或 `峰处 … 的 Metropolis 接受率只有 0.0…`；`05_mc/mc_diag.json` 里 `n_eff` 很小、`sem_binned / sem_naive` 很大 | 这不是 bug，是**这组 `--therm/--meas` 不够**：① `n_eff = n_meas/(1+2τ)`，`--meas` 至少要 `200(1+2τ)` 才谈得上误差棒 —— 按日志给的下界调，或加 `--seeds` 做多种子平均；② 接受率 < 0.05 说明该温度点自旋冻住（`--tmin` 太低，或该尺寸在低温下能垒太高）⇒ 收窄 `--tmin` 或加大 `--size`；③ `naive SEM` 比 `binning SEM` 小 3 倍以上时**只能用 binning**，别拿 naive 当误差棒。峰落窗口边界时这几条告警故意不报（先解决 `peak_boundary`）—— 看 `mc_diag.json` 的 `warnings` 字段比看这条表更准 |
