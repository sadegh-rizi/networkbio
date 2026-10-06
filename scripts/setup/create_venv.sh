#!/usr/bin/env bash
# Create the project venv (.venv) in the repository root and install the
# pinned packages from requirements.lock.txt. Run from WSL:
#   bash scripts/setup/create_venv.sh
# Refuses to touch an existing .venv. Writes the installed versions to
# logs/pip_freeze_<date>.txt.
set -euo pipefail
cd "$(dirname "$0")/../.."

PY="${PYTHON:-python3}"
"$PY" - <<'PYEOF'
import sys
if sys.version_info < (3, 11):
    sys.exit(f"Python >= 3.11 needed (corneto 1.0.0rc8); found {sys.version.split()[0]}")
print("python", sys.version.split()[0])
PYEOF

if [ -e .venv ]; then
    echo ".venv already exists; not modifying it. Remove it yourself to rebuild." >&2
    exit 1
fi

"$PY" -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.lock.txt

mkdir -p logs
.venv/bin/python -m pip freeze > "logs/pip_freeze_$(date +%F).txt"
.venv/bin/python - <<'PYEOF'
import corneto, decoupler, highspy, cvxpy
print("corneto", corneto.__version__, "| decoupler", decoupler.__version__)
print("cvxpy solvers:", cvxpy.installed_solvers())
PYEOF
echo "venv ready: .venv"
