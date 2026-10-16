# DFT 回收报告 (ingest_dft.py)

- 生成时间 (UTC): `2026-10-16T04:55:59+00:00`
- 来源 runs: `C:\Users\吴\AppData\Local\Temp\dsh-HxaRfu\_dryrun_runs_10cudyn1`
- 目标 cases: `C:\Users\吴\AppData\Local\Temp\dsh-HxaRfu\kag_example\02_supercell\cases`
- 能量口径: `free energy TOTEN (eV)` (混口径容差 0.001 eV)
- **本报告只证明 OUTCAR 过了格式/收敛/结构/口径闸门, 不证明 J/T_N 的物理正确性。**

## 四个磁态

| 序 | 能量 (eV) | ionic steps | SCF 收敛 | 与 fm 之差 (meV) | sha256 前 12 位 |
|---|---|---|---|---|---|
| fm | -0.2592 | 1 | 是 | 0.000 | `99a40a1399aa` |
| q0_120 | 0.1053 | 1 | 是 | 364.500 | `31eaf5e8aa6d` |
| sq3_alt | -0.0336 | 1 | 是 | 225.600 | `b1dcc2fb1087` |
| sq3_updntri | -0.0216 | 1 | 是 | 237.600 | `572108faa5ec` |

- 能量离散度: **0.3645 eV/胞** (阈值见 manifest `checks`)
- LNONCOLLINEAR: `{'fm': True, 'q0_120': True, 'sq3_updntri': True, 'sq3_alt': True}`   LSORBIT: `{'fm': True, 'q0_120': True, 'sq3_updntri': True, 'sq3_alt': True}`
- ENCUT: `{'fm': 520.0, 'q0_120': 520.0, 'sq3_updntri': 520.0, 'sq3_alt': 520.0}`   NKPTS: `{'fm': 6, 'q0_120': 6, 'sq3_updntri': 6, 'sq3_alt': 6}`   VASP: `{'fm': '6.4.2', 'q0_120': '6.4.2', 'sq3_updntri': '6.4.2', 'sq3_alt': '6.4.2'}`

## 警告 (不拦写入, 但请核对)

- `[accuracy]` (fm) 没打印 `reached required accuracy` (静态 NSW=1 一般不打印, 属正常; 提醒确认一下)
- `[accuracy]` (q0_120) 没打印 `reached required accuracy` (静态 NSW=1 一般不打印, 属正常; 提醒确认一下)
- `[accuracy]` (sq3_alt) 没打印 `reached required accuracy` (静态 NSW=1 一般不打印, 属正常; 提醒确认一下)
- `[accuracy]` (sq3_updntri) 没打印 `reached required accuracy` (静态 NSW=1 一般不打印, 属正常; 提醒确认一下)

## 下一步

```
# 反解 J (kagome)
python run_2d_magnet.py solve --outdir <proj>
# 或直接
python scripts/solve_J_kagome.py <proj>/02_supercell/cases
```
