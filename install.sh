#!/usr/bin/env bash
# One-shot setup: virtual environment, dependencies, environment check.
# Usage:  bash install.sh            (creates .venv in this directory)
#         bash install.sh --venv /path/to/venv
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV="$HERE/.venv"

while [ $# -gt 0 ]; do
  case "$1" in
    --venv) VENV="$2"; shift 2 ;;
    -h|--help) sed -n '2,4p' "$0"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

PY="${PYTHON:-python3}"
command -v "$PY" >/dev/null 2>&1 || { echo "[error] $PY not found; set PYTHON=..." >&2; exit 1; }
"$PY" -c 'import sys; assert sys.version_info >= (3, 9), sys.version' \
  || { echo "[error] Python >= 3.9 required" >&2; exit 1; }

echo "[1/3] venv -> $VENV"
"$PY" -m venv "$VENV"

# shellcheck disable=SC1091
if [ -f "$VENV/bin/activate" ]; then . "$VENV/bin/activate"; else . "$VENV/Scripts/activate"; fi

echo "[2/3] dependencies"
python -m pip install --upgrade pip >/dev/null
python -m pip install -r "$HERE/requirements.txt"

echo "[3/3] doctor"
set +e
python "$HERE/run_2d_magnet.py" doctor
rc=$?
set -e
if [ $rc -ne 0 ]; then
  echo
  echo "[warn] doctor reported WARN lines. They are environment hints, not code errors."
  echo "       Symptom -> fix table: knowledge/KNOWN_ISSUES.md"
else
  echo
  echo "[OK] environment clean. Next: README.md 'Quick start'."
fi
