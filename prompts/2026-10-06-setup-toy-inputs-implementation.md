# Environment, CORNETO toy problem and stage-00 input tables

Date: 2026-10-06
Model: claude-opus-5-5 (Claude, Cowork)
Settings: default
Plan: doc/decisions/2026-10-06-setup-toy-inputs-plan.md
Commit(s): none yet (add after the user commits)
Transcript: not saved

## Prompt actually used

> can you do all this , create the venv and install the packages and build the results

("All this" = the steps proposed in the same session: rewrite the plan with
the user's answers, pin requirements, write and run the CORNETO toy problem,
build the stage-00 input tables.) Follow-up: "can you do this on my wsl now".

## What it produced

- `requirements.txt` (direct pins), `requirements.lock.txt` (uv pip compile,
  Python 3.12), `scripts/setup/create_venv.sh`,
  `scripts/setup/run_setup_and_stage00.sh`
- `scripts/analysis/00_toy_corneto.py`, `scripts/analysis/01_build_inputs.py`,
  `src/inputs.py`, `tests/test_inputs.py`
- Docs: `doc/decisions/2026-10-06-setup-toy-inputs-plan.md` (new), answers
  section in `doc/decisions/2026-09-29-corneto-ionescu-inscs-plan.md`,
  `doc/decisions/README.md`, `doc/analysis-state.md`, `AGENTS.md`
  (environment state), `doc/external-tools.md`, `doc/context/*`,
  `doc/datasets.md`.

## How it was run

The agent could not type into the WSL terminal. It developed and tested the
scripts in its own cloud environment with the same lock file, copied them
into the repo, and the user ran:

    bash scripts/setup/run_setup_and_stage00.sh 2>&1 | tee logs/setup_stage00.log

Result in WSL (Ubuntu 24.04, Python 3.12.3): pytest 5 passed; toy problem
5/5 checks PASS; `ALL STEPS FINISHED`. The five stage-00 summary tables are
byte-identical (md5) between the agent's test run and the WSL run.

## What was changed by hand afterwards

None so far.

## Outputs stamped

`results/corneto_toy/04_inference/toy_checks.json` and
`results/ionescu_corneto/00_inputs/00_inputs.provenance.json` record package
versions and platform; no code commit hash yet (repository not committed).
