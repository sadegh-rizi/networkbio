#!/usr/bin/env bash
# One command for plan steps 0-2: create the venv, run the tests, the CORNETO
# toy problem and the stage-00 input tables. Run from WSL, repo root:
#   bash scripts/setup/run_setup_and_stage00.sh 2>&1 | tee logs/setup_stage00.log
set -euo pipefail
cd "$(dirname "$0")/../.."
mkdir -p logs
[ -d .venv ] || bash scripts/setup/create_venv.sh
.venv/bin/python -m pytest -q tests
.venv/bin/python scripts/analysis/00_toy_corneto.py
.venv/bin/python scripts/analysis/01_build_inputs.py
echo "ALL STEPS FINISHED"
