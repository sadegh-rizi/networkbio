#!/usr/bin/env bash
# Run the stage-03 synthetic tests and the real PKN build from WSL.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"
PYTHON="${ROOT}/.venv/bin/python"
[[ -x "${PYTHON}" ]] || { printf 'Missing project venv; no environment changes were made.\n' >&2; exit 1; }
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
"${PYTHON}" -m pytest -q tests
"${PYTHON}" scripts/analysis/04_build_pkn.py
printf 'STAGE 03 FINISHED\n'
