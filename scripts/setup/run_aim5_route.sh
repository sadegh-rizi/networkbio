#!/usr/bin/env bash
# Run the Aim-5 synthetic tests and structural route audit from WSL.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"
PYTHON="${ROOT}/.venv/bin/python"
[[ -x "${PYTHON}" ]] || { printf 'Missing project venv; no environment changes were made.\n' >&2; exit 1; }
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
"${PYTHON}" -m pytest -q tests/test_metabolic_route.py
"${PYTHON}" scripts/analysis/04b_aim5_metabolic_route.py
printf 'AIM 5 ROUTE FINISHED\n'
