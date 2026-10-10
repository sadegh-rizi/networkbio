#!/usr/bin/env bash
# Run the pre-registered stage-05 robustness controls from WSL.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"
PYTHON="${ROOT}/.venv/bin/python"
[[ -x "${PYTHON}" ]] || { printf 'Missing project venv; no environment changes were made.\n' >&2; exit 1; }
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
"${PYTHON}" -m pytest -q tests
"${PYTHON}" scripts/analysis/06_corneto_robustness.py
printf 'STAGE 05 FINISHED\n'
