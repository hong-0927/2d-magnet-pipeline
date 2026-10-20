# AGENTS.md

Instructions for any agent (or human) working in this repository. Read this first, then
[`skill/2d-magnet/SKILL.md`](skill/2d-magnet/SKILL.md) for the full driving manual.

## What this repository is

A deterministic pipeline that turns a monolayer structure into exchange constants `J` and a
Monte Carlo `T_C` / `T_N`, for eight lattice types, with every number sitting behind a gate.
It is harness-neutral: anything that can read Markdown and run `python` can drive it.

## Hard rules

1. **Never invent a parameter.** If a value is unknown, the pipeline stops and names the
   missing flag. The element table deliberately omits Fe, Co oxides and other mixed-valence
   species; a literature entry with `J = null` is rejected rather than guessed. Passing
   `--force` / `--lenient` / `--allow-unverified-states` is allowed only when the report says
   so explicitly and says what was weakened.
2. **Gates fail closed.** `solve` runs `scripts/verify_states.py` first and exits non-zero
   without writing `J_result.json` when the per-site moments do not match the target orders.
   Do not route around a red gate; fix the input or state the limitation.
3. **Every number carries provenance.** `source` (`literature` / `demo` / `dft`), `confidence`,
   `residual`, `cond(A)` and `states_verification` live inside `J_result.json`. A number without
   those fields is not reportable.
4. **Capability boundaries are machine-readable.** Before claiming anything about `D`, band
   splitting, `Z2`, or a physical `T_C`, read [`capabilities.json`](capabilities.json)
   (`python scripts/capabilities.py --show`). `physical_verified` is `false` everywhere.
5. **Report what actually ran.** Do not summarise a stage that failed; do not round a
   diagnostic away. `--quiet` keeps conclusion-level warnings on stderr by design.

## Where things live

| Path | What it is |
|---|---|
| `run_2d_magnet.py` | the driver: `fetch build dft ingest kconv mc solve render report doctor` |
| `scripts/` | deterministic tools, all with `--help` and most with `--selftest` |
| `knowledge/` | reference library: workflow, error catalog, limitations, DMI memo, runbook |
| `skill/2d-magnet/SKILL.md` | the agent skill (conventions, gates, per-stage recipes) |
| `tests/` | `smoke_test.py` (29 sections) + 8 contract tests + `check_sdist.py` |
| `examples/` | two complete artifact trees you can inspect without VASP |
| `capabilities.json` | machine-readable capability boundary |

## Before you hand work back

```bash
python tests/smoke_test.py
python tests/test_magnetic_orders.py      # and the other 7 contract tests
python scripts/magnetic_orders.py --regress
python scripts/capabilities.py --check
python tests/check_sdist.py "dist/*.tar.gz"
```

All five must be green. Hand the diff to the user for the commit; do not commit on your own.
