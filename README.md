# 2d-magnet-pipeline

Turn a monolayer structure into exchange constants `J1/J2/J3` and a Monte Carlo estimate of `T_C` / `T_N` — with every number behind a fail-closed gate.

The input is a `POSCAR`: from the Materials Project, from the public COD / OQMD / AFLOW / JARVIS chain, or your own file. The pipeline detects the magnetic lattice, builds the supercell and the magnetic-order set, assembles VASP inputs, ingests the returned `OUTCAR`s, solves the exchange constants, and runs a classical Heisenberg Monte Carlo. Two routes to J:

- **literature or user-given J → forward Monte Carlo** — no VASP needed;
- **a DFT many-state energy mapping solved back to J** — needs VASP energies.

## Highlights

| Capability | What it gives you |
|---|---|
| **Eight magnetic lattices, one spec layer** | honeycomb, triangular (120° non-collinear), molecular cluster, two-sublattice ferromagnetic network (Fe₃GeTe₂-like), rectangular (CrSBr-like), C-centred rectangular (CrOCl-like), square (La₂CuO₄ / K₂NiF₄-like), hexagonal kagome (Fe₃Sn-like) — all built and solved through the same specification layer, all reporting into the same `J_result.json` |
| **A documented way to J, not just a number** | `J_result.json` records `source` (`literature` / `user` / `dft` / `demo`), `cond(A)`, residual, `E0` and whether the magnetic states were verified per site |
| **Magnetic states verified before inversion** | `scripts/verify_states.py` compares every ion's moment against its target and detects two orders collapsing into one state; `solve` refuses to write J when that check fails |
| **Fail-closed ingestion** | 11 gates on each `OUTCAR`: format, convergence, structure, MAGMOM consistency, energy convention |
| **Estimates with their own error bars** | Binder cumulant `U4`, finite-size scaling with a `1/L` reference, and per-temperature `tau_int`, `n_eff`, acceptance rate and binned SEM in `mc_diag.json`; a peak on the window edge is reported as unresolved instead of being trusted |
| **Non-collinear and SOC support** | `LNONCOLLINEAR` / `LSORBIT` / `SAXIS`, canted MAGMOM sets, and `scripts/solve_J_crect_cant.py` for the canted CrOCl family |
| **No invented parameters** | `--S/--moment/--U` fall back to a built-in table for Mn/Cr/Ni/Co/V/Gd only; Fe is deliberately absent (effective spin is ambiguous), so the driver stops and asks, and `null` literature entries raise an error instead of running |
| **Boundaries you can query** | [`capabilities.json`](capabilities.json) states what each lattice and quantity reaches; `scripts/capabilities.py --check` keeps it in sync with the code, and `--show` prints the readable version |
| **Also an agent skill** | `dsh/2d-magnet/SKILL.md` carries the conventions and gates for Claude Code, Codex or any agent that reads instruction files |
| **Reproducible offline** | `examples/` holds two complete artifact trees you can inspect and regenerate without VASP; release packages are built deterministically by `scripts/make_release.py` |
| **Altermagnet check** | `scripts/altermagnet_check.py` tests whether the two sublattices are related by a spin-preserving operation (交变磁体 geometry, scalar order parameter) |

## Workflow

```text
POSCAR    -> fetch     reorient, vacuum                         01_structure/POSCAR
          -> build     lattice, supercell, magnetic orders      02_supercell/cases/
          -> dft       INCAR / KPOINTS / MAGMOM / POTCAR.spec   03_dft/, runs/<order>/
          [ VASP runs outside the pipeline ]
          -> ingest    OUTCAR back, 11 fail-closed checks       cases/<order>/OUTCAR
          -> solve     J1 J2 J3 (meV), states verified           02_supercell/cases/J_result.json
          -> mc        Metropolis + Binder U4 -> T_C / T_N       05_mc/mc_mt.csv, mc_diag.json
          -> render    supercell and order figures              figures/supercell_orders.png
          -> report    provenance card                          workflow.md
```

`kconv` runs only the K-points convergence step (`scripts/k_convergence.py`; `--kpoints` pins a grid). Every stage writes its inputs, provenance and exit status; nothing is derived from a session you cannot re-run.

`--quiet` suppresses progress output but never the conclusion-level warnings — unresolved `T_C`, ill-conditioned `cond(A)`, collapsed states. Use `--log <file>` to keep the full output of every stage.

## Install

**Requirements**

| Item | Detail |
|---|---|
| Python | ≥ 3.9 (the release regression runs on 3.12) |
| Packages | `numpy` ≥ 1.24, `numba` ≥ 0.59, `matplotlib` ≥ 3.7, `requests` ≥ 2.31 — see `requirements.txt` |
| Optional package | `Pillow` ≥ 9.0 — figure fallback when `matplotlib` cannot be imported |
| VASP | optional: `dft` / `ingest` / `solve` need your own VASP runs; literature-J → Monte Carlo does not |
| Network | optional: only `fetch` downloads structures; `examples/` and `tests/fixtures/` work offline |
| Disk | ≈ 15 MB for the repository; each run writes a few MB under `--outdir` |

**1. Get the code**

```powershell
git clone <repo-url> 2d_magnet_pipeline
cd 2d_magnet_pipeline
```

No git? Download the release archive (`2d-magnet-pipeline-0.1.0.tar.gz` / `.zip` from the Releases page) and unpack it — the archive holds the same tree.

Steps 2–4 are scripted. From the repository root:

```powershell
powershell -ExecutionPolicy Bypass -File install.ps1
```

```bash
bash install.sh
```

Both create `.venv`, install `requirements.txt` and run `doctor`; pass `--venv <path>` (bash: `-Venv <path>` in PowerShell) to use another location. On Windows, run the PowerShell one with `-ExecutionPolicy Bypass` as shown — the venv activation script would otherwise be blocked.

**2. Create and activate a virtual environment**

```powershell
py -m venv .venv
.venv\Scripts\Activate.ps1
```

```bash
python3 -m venv .venv
. .venv/bin/activate
```

If PowerShell refuses to run the activation script (`running scripts is disabled on this system`), either activate from `cmd` with `.venv\Scripts\activate.bat`, or allow scripts for the current session only: `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass`. On Windows use the `py` launcher — it bypasses the Microsoft Store `python.exe` stub; if `py` is not installed, use the full path to `python.exe`.

**3. Install the dependencies**

```powershell
py -m pip install --upgrade pip
py -m pip install -r requirements.txt
```

`matplotlib` is the biggest item and pulls `pyparsing`, `cycler`, `contourpy`, `fontTools` and `kiwisolver`. If it cannot be installed (offline machine, or a minimal site), install `Pillow>=9.0` instead: all six rendering paths then use a pure-Pillow implementation and the figures are still written (packaged extra: `pip install -e ".[render-fallback]"`). `requests` is only used by `fetch`; without it everything except structure download works.

**4. Verify the installation**

```powershell
py run_2d_magnet.py doctor
```

`doctor` prints one line per check, then `---`, then the error-catalog entry for every WARN. Exit code `0` means no WARN; `1` means at least one environment problem — it does not mean the code is broken. A healthy install reports `[OK]` for Python, `numpy`, `numba`, `scripts/`, `literature_J.json` and `tests/fixtures/`; the `MP_API_KEY` and `CWD` lines are informational only.

**5. Optional: install as a command**

```powershell
py -m pip install -e .
2d-magnet --version      # 2d-magnet 0.1.0
2d-magnet doctor
```

**6. Optional: run it as an agent skill**

`dsh/2d-magnet/` is a self-contained skill (conventions, gates, driving instructions). Copy the folder into your harness skill directory — for Claude Code `~/.claude/skills/`, or `.claude/skills/` inside a project; for an `AGENTS.md`-style harness point the agent at `dsh/2d-magnet/SKILL.md`. The skill drives the same `run_2d_magnet.py`; nothing else needs installing.

**7. Optional: the VASP side**

`dft` writes `INCAR`, `KPOINTS`, `MAGMOM` and `POTCAR.spec` (recommended POTCAR label per element) under `03_dft/` and `runs/<order>/`. POTCAR files are licensed and are **not** shipped: assemble them from your VASP distribution into `runs/<order>/POTCAR` before submitting. `ingest` reads only the `OUTCAR`, so the calculation itself can run on any machine.

**Troubleshooting**

| Symptom | Fix |
|---|---|
| `python` opens the Microsoft Store, or `command not found` | use `py` (Windows launcher) or the full path to `python.exe` |
| `ModuleNotFoundError: No module named 'pyparsing'` while rendering | `matplotlib` is incomplete — install `Pillow`, the fallback renderer takes over |
| Garbled Chinese output, `UnicodeEncodeError` | set `$env:PYTHONIOENCODING="utf-8"` before running |
| `doctor` reports `CWD 不可写` | run from a writable directory, or always pass an explicit `--outdir` |
| `fetch` returns 401/403 or an empty structure | set `MP_API_KEY` or pass `--mp-key`, or supply your own file with `--poscar` |
| `[params] ... 不在内置元素参数表` | pass `--S --moment --U` explicitly, or pick a `--model` / `--literature` entry |

Every WARN and error key above has its full symptom → fix text in [`KNOWN_ISSUES.md`](knowledge/KNOWN_ISSUES.md).

The code is also usable as a library (`import run_2d_magnet`) and as the console script `2d-magnet` after `pip install .`.

## Quick start

`examples/` holds two complete runs you can reproduce without VASP — literature J in, `T_C` out:

```powershell
py run_2d_magnet.py --model la2cuo4 --poscar tests/fixtures/La2CuO4_POSCAR `
   --outdir C:\work\la2 --skip-dft --mc-size 12 12 --mc-seeds 2
py scripts\solve_J_square.py C:\work\la2\02_supercell\cases --demo
```

`examples/la2cuo4_square/` is the resulting tree; `examples/kagome_fe3sn_dryrun/` is a kagome run through the full DFT shape with synthetic energies. Each has a `README.md` listing what is produced.

## Commands

| Subcommand | Does |
|---|---|
| `fetch` | structure from Materials Project / COD / OQMD / AFLOW / JARVIS → `01_structure/POSCAR` |
| `build` | lattice detection, supercell, magnetic-order set → `02_supercell/cases/` |
| `dft` | INCAR / KPOINTS / MAGMOM / POTCAR.spec → `03_dft/`, `runs/<order>/` |
| `ingest` | OUTCAR validation (11 gates) → `cases/<order>/OUTCAR` |
| `kconv` | K-points convergence; only KPOINTS is varied |
| `solve` | many-state energy mapping → `J_result.json` |
| `mc` | Metropolis / Binder U4 → `T_C` or `T_N` |
| `render` | supercell and order figures |
| `report` | provenance card → `workflow.md` |
| `doctor` | environment self-check against `KNOWN_ISSUES.md` |

Presets: `la2cuo4`, `crsbr`, `crocl`, `k2nif4`, `fgt`, `mnps3`, `nio`, `gdi2`. Driver options (`--model`, `--literature`, `--element`, `--S`, `--moment`, `--U`, `--mc-order`, `--mc-seeds`, `--quiet`, `--log`) are listed by `py run_2d_magnet.py <subcommand> --help`.

## Outputs

```text
01_structure/POSCAR              oriented, vacuum-trimmed input
02_supercell/cases/              one directory per magnetic order: POSCAR, magmom.txt, OUTCAR
02_supercell/cases/J_result.json J with provenance: source, confidence, cond(A), residual, E0
02_supercell/cases/states_verification.json   per-site moment check that gates solve
03_dft/, runs/<order>/           INCAR, KPOINTS, POTCAR.spec, MAGMOM
05_mc/mc_mt.csv, mc_mt.png       T, order parameter, susceptibility, Binder U4
05_mc/mc_diag.json               acceptance, tau_int, n_eff, naive vs binned SEM
figures/supercell_orders.png     supercell with exchange bonds and the order set
workflow.md                      provenance card for the run
```

## Verification（脚本级自检 / 发布闸门）

```powershell
py run_2d_magnet.py doctor
py scripts\capabilities.py --check
py tests\smoke_test.py
py tests\check_sdist.py "dist\*.tar.gz"
```

The contract suite is eight test files under `tests/` (`test_magnetic_orders.py`, `test_square_cases.py`, `test_crect_cases.py`, `test_kagome_cases.py`, `test_fgt_e2e.py`, `test_altermagnet.py`, `test_crect_cant.py`, `test_ingest.py`), plus `scripts/magnetic_orders.py --regress` for the positional regression baseline. Script-level self-checks:

```text
python scripts/mc_2d.py --selftest            energy convention, single-ion anisotropy, tau_int
python scripts/mc_fss.py --selftest           crossing choice, 1/L extrapolation
python scripts/verify_states.py --selftest    per-site moments, collapse detection
python scripts/dmi_matrix.py --selftest       DMI matrix vs exchange energy
python scripts/k_convergence.py --selftest    K-points planner
```

## Documentation map（该读哪份）

| Read | For |
|---|---|
| [`README.md`](README.md) / [`README.zh.md`](README.zh.md) | this manual (English / 中文) |
| [`AGENTS.md`](AGENTS.md) | the five hard rules for any agent driving this repo |
| [`STRUCTURE.md`](STRUCTURE.md) | the organization convention — where new content goes |
| [`KNOWN_ISSUES.md`](knowledge/KNOWN_ISSUES.md) | symptom → fix table, 27 error keys; `doctor` output points here |
| [`WORKFLOW.md`](knowledge/WORKFLOW.md) | per-stage detail: file by file, option by option |
| [`dsh/2d-magnet/SKILL.md`](dsh/2d-magnet/SKILL.md) | agent skill: conventions, gates, how to drive the pipeline |
| [`capabilities.json`](capabilities.json) | machine-readable capability boundary — lattice, quantity, supported level |
| [`DMI_MTI_MEMO.md`](knowledge/DMI_MTI_MEMO.md) | DMI / MTI: how to tell them apart and what may be claimed |
| [`LIMITATIONS.md`](knowledge/LIMITATIONS.md) | boundary statements — what each number does and does not mean (18 items) |
| [`examples/`](examples/) | two complete artifact trees you can inspect without VASP |
| [`P4_2_KAGOME_RUNBOOK.md`](knowledge/P4_2_KAGOME_RUNBOOK.md) | the kagome physics runbook, ready for machine time |
| [`CHANGELOG.md`](CHANGELOG.md) · [`CITATION.cff`](CITATION.cff) | version history · citation metadata |


## License

MIT — see [`LICENSE`](LICENSE).

## Citation

If this pipeline contributed to work you publish, cite it via [`CITATION.cff`](CITATION.cff):

```text
Wu, Sihong. 2d-magnet-pipeline: from a 2D magnet structure to exchange constants
and Monte Carlo T_C / T_N. Version 0.1.0, 2026. MIT license.
```
