# STRUCTURE.md

The organization convention. Read it before adding content, so new material lands in the
right place and the release gate keeps passing.

## Four content kinds, deliberately distinct

| Kind | Where | Rule |
|---|---|---|
| **Driver + tools** | `run_2d_magnet.py`, `scripts/` | Deterministic. Every script takes flags, prints machine-parseable output, exits non-zero on failure, and has `--help`; most have `--selftest`. No hidden defaults for physical quantities. |
| **Skill / procedure** | `dsh/2d-magnet/SKILL.md` (+ `AGENTS.md`) | How a task is driven: conventions, gates, per-stage recipes, error recovery. Prose, no code. |
| **Reference library** | `knowledge/` | Flat, tool-agnostic science and practice: workflow detail, error catalog, boundary statements, DMI memo, kagome runbook. Read it for ideas; it is not a skill and cannot be executed. |
| **Evidence** | `tests/`, `examples/` | `tests/` holds the gates (smoke sections + contract tests + sdist check). `examples/` holds two complete artifact trees produced by the pipeline, kept so a reader can inspect real file contents without running VASP. |

Faithful reporting is a structural rule, not a stylistic one: numbers live in JSON with
provenance, capability claims live in `capabilities.json`, and long-form boundaries live in
`knowledge/LIMITATIONS.md`.

## Root files

| File | Purpose |
|---|---|
| `README.md` / `README.zh.md` | the manual (English / 中文) |
| `AGENTS.md` | agent guardrails (read first) |
| `STRUCTURE.md` | this convention |
| `LICENSE` | MIT |
| `CITATION.cff` | citation metadata |
| `CHANGELOG.md` / `RELEASE_NOTES_v0.1.0.md` | version history / release body |
| `capabilities.json` | machine-readable capability boundary |
| `pyproject.toml`, `requirements.txt`, `MANIFEST.in` | packaging |
| `install.sh`, `install.ps1` | one-shot environment setup (venv + deps + `doctor`) |

## Adding a lattice type

1. Write the spec layer first (`scripts/magnetic_orders.py`), so orders and the coefficient
   matrix exist before any solver.
2. Add the solver, then a cases test in `tests/`, then a selftest, then the smoke section.
3. Register it in `capabilities.json` with `physical_verified: false` and a `why`.
4. Update the checklist in `templates/new_lattice_checklist.md`, the docs map in both READMEs,
   and the counts in `tests/check_sdist.py` if new required files appear.

The full checklist, including the gates a new lattice has to pass, is
`templates/new_lattice_checklist.md`.
