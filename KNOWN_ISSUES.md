# Known Issues / Error Catalog

运行 `python run_2d_magnet.py doctor` 会做环境自检，并把**命中的已知问题**映射成「症状 → 修复」。
本文件是同一份目录的静态版，覆盖两类条目：**环境可探**（doctor 能当场命中）和**运行时**（跑到某阶段才出现）。

## 环境自检可命中的条目

| key | 症状 | 修复 |
|---|---|---|
| `py39` | 语法/类型报错，或 numpy/numba 版本不兼容 | 升级到 Python ≥ 3.9（如 `py -3.12` 或 conda） |
| `missing_dep` | `ImportError: No module named 'xxx'` | `pip install -r requirements.txt` |
| `packaging` | render 阶段 `ModuleNotFoundError: No module named 'packaging'`（matplotlib 顶层 import 连带需要它） | `pip install packaging` |
| `scripts_missing` | driver 报缺脚本 / 子阶段找不到脚本 | 确认 `scripts/` 目录完整（重新 `git clone` 或解压分发包） |
| `lit_missing` | `--literature <key>` 报找不到条目 / `literature_J.json` 缺失 | 确认 `scripts/literature_J.json` 存在且键名正确 |
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
