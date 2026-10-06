# AGENTS.md

Guidance for AI coding agents working in this repository. Claude Code should also read `CLAUDE.md` and `doc/agent-rules.md`.

**Start with `doc/analysis-state.md`** (current state, open decisions, next steps), then `doc/context/` for detail.

## Project overview

A two-week Network Biology course project. It applies CORNETO (multi-sample, prior-knowledge-based network inference) to public multi-omics data from progressive-MS iNSC lines (Ionescu et al. 2024; Park et al. 2025). Neither the method's outputs nor the datasets are treated as ground truth: see `doc/agent-rules.md`.

Relevant files and directories:

- `data/`: source datasets. Treat as immutable research inputs unless the user explicitly requests a change.
- `doc/datasets.md`: what each dataset is, its design, and what is known vs unknown.
- `doc/agent-rules.md`: canonical shared rules for scientific validity, data handling, statistics, plotting, and verification.
- `doc/analysis-state.md`: single entry point to the current state.
- `doc/decisions/`: confirmed analysis plans, the handoff contract between agents and sessions.
- `doc/external-tools.md`: ledger of third-party software with versions.
- `prompts/`: provenance records for AI-generated code and analysis.
- `requirements.txt`: Python package requirements.

## Instruction priority

1. Follow the user's explicit request.
2. Follow this file and `doc/agent-rules.md`.
3. Preserve the current analysis and data unless a requested change requires otherwise.

When instructions conflict or scientific intent is unclear, stop and state the ambiguity rather than silently choosing a biological or statistical method.

## Planning before implementation

Before writing, editing, or running any code (notebook cells, modules, scripts, or dependency changes), write a short plan: the files you will touch, the approach, and the validation command(s). List every open design, statistical, or biological question implied by the task and ask the user directly; do not silently choose a default for an ambiguous decision. Wait for explicit confirmation of the plan before making any edit.

Confirmed plans go in `doc/decisions/YYYY-MM-DD-<topic>-plan.md`. That document, not the chat, is what the implementing agent reads.

## Environment setup (WSL2)

The user works in Ubuntu under WSL2 with PyCharm. The Claude Science app runs natively on Windows, and its sandbox Python is **not** the project environment. Do not run analysis in the sandbox and report it as a project result. Prepare the commands and let the user run them in WSL.

**Project Python.** A `venv` in the repository root, with packages from `requirements.txt`:

```
"/mnt/d/#University/#MS-SystemsBio/Master_Thesis/repos/networkbio/.venv/bin/python"
```

Always quote that path: the repository lives under directories beginning with `#`, which is the shell comment character. Use the sibling `.venv/bin/pytest` likewise.

**Known state of that environment: created 2026-10-06** (Python 3.12.3, from `requirements.lock.txt` via `scripts/setup/create_venv.sh`; `pip freeze` in `logs/pip_freeze_2026-10-06.txt`). corneto 1.0.0rc8, decoupler 2.2.0, omnipath 1.0.12, cvxpy 1.9.3 with HiGHS and SCIP are installed. Anything not in the lock file is not installed.

**Never `pip install` without asking.** Report what is missing and stop. Do not create, rebuild, upgrade, or change any environment or dependency unless asked.

**CORNETO and its solver.** CORNETO is the method of this project, so it belongs in the project venv (`requirements.txt`), not in an isolated per-tool environment. Which MILP solver to use is an open decision recorded in `doc/decisions/`; test it on a toy example before any real run.

**No R** unless asked. If an R step is ever needed, hand data over on disk with a documented schema; never `reticulate` or `rpy2`.

## Development workflow

- Prefer small, reviewable changes with one stated purpose.
- Before modifying code, inspect the input schema. Do not infer column meanings from names alone.
- Reusable analysis code goes in `src/` as plain functions; notebooks orchestrate, inspect, and visualise.
- Keep paths portable: `pathlib.Path` and project-relative paths, no hard-coded Windows or WSL absolute paths in code.
- Do not add a dependency merely because an agent prefers it. Explain the need and ask before changing `requirements.txt`.

## Context and compute efficiency

- Never read a large table whole to find out what it is. Use `pandas.read_csv(..., nrows=5)` or `usecols=` to probe schema.
- Never `read` a whole `.ipynb`; extract the sources or grep.
- Do not read back a file you just wrote or patched.
- Cache expensive intermediates (activity tables, processed PKN, solver solutions) to `results/` with a provenance sidecar and reload them.
- A MILP solve can be slow and its run time can grow steeply with network size. Estimate the size (nodes, edges, samples) and set a time limit before running.

## Data handling and scientific integrity

- Never delete, rename, overwrite, or normalise in place any file in `data/`.
- Keep source, interim, and derived results distinct.
- The data are derived from patient fibroblasts. Do not attempt to re-identify donors or cross-link them to other sources. Cite the publication and accession for every dataset used in a figure or claim.
- Read `doc/agent-rules.md` before changing preprocessing, activity estimation, the PKN, solver settings, statistics, or figures.

## Verification

Before handing off changes, run the narrowest useful check and report the real result. For scripts, run a relevant import, targeted script, or test. For network inference, the check is a toy problem with a known answer before the real data. Do not claim scientific validity merely because code runs.

## Git workflow

- Start every task on a new branch `task/<short-name>`. Never work on main.
- Verify the tree is clean before starting; stop and report if it is not.
- Commit each generation event separately, with trailers:
  ```
  Prompt: prompts/<file>.md
  Model: <model-version>
  ```
- Never `git add .`; stage named files only. No secrets.
- Never push, force-push, rebase shared branches, `reset --hard`, `clean -fd`, or squash. Merging and pushing are done by the human.
- Never commit anything from `data/`.
