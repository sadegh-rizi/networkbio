#!/usr/bin/env bash
# Run the stage-05 robustness controls from WSL (see doc/decisions/2026-10-11-stage05-robustness-addendum.md).
# Extra arguments are passed through, e.g. --resume --cosmos.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"
PYTHON="${ROOT}/.venv/bin/python"
[[ -x "${PYTHON}" ]] || { printf 'Missing project venv; no environment changes were made.\n' >&2; exit 1; }
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
"${PYTHON}" -m pytest -q tests
"${PYTHON}" scripts/analysis/06_corneto_robustness.py "$@"
printf 'STAGE 05 FINISHED\n'
