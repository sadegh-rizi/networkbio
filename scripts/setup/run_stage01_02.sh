#!/usr/bin/env bash
# WSL: run from any directory; all outputs are project-relative.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"
PYTHON="${ROOT}/.venv/bin/python"
[[ -x "${PYTHON}" ]] || { printf 'Missing project venv; no environment changes were made.\n' >&2; exit 1; }
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
"${PYTHON}" -m pytest -q tests &&
    "${PYTHON}" scripts/analysis/02_preprocess.py &&
    "${PYTHON}" scripts/analysis/03_tf_activities.py &&
    printf 'STAGES 01-02 FINISHED\n'
