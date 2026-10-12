#!/usr/bin/env bash
# Run stage-04 tests and the preregistered CORNETO inference from WSL.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "${ROOT}"
PYTHON="${ROOT}/.venv/bin/python"
[[ -x "${PYTHON}" ]] || { printf 'Missing project venv; no environment changes were made.\n' >&2; exit 1; }
export MPLBACKEND=Agg
export PYTHONUNBUFFERED=1
"${PYTHON}" -m pytest -q tests
# Extra arguments are passed through, e.g. --resume [--backfill-log <log>].
"${PYTHON}" scripts/analysis/05_corneto_inference.py "$@"
printf 'STAGE 04 FINISHED\n'
