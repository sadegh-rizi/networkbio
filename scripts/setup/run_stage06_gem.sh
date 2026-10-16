#!/usr/bin/env bash
# Run the Stage-06 Recon3D FBA/iMAT context analysis from WSL.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"
PYTHON="${ROOT}/.venv/bin/python"
[[ -x "${PYTHON}" ]] || { printf 'Missing project venv; no environment changes were made.\n' >&2; exit 1; }
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
if [[ "${RUN_TESTS:-1}" != "0" ]]; then
  "${PYTHON}" -m pytest -q tests/test_gem.py
fi
"${PYTHON}" scripts/analysis/08_metabolic_context.py "$@"
printf 'STAGE 06 FINISHED\n'
