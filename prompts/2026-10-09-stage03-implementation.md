# Stage-03 PKN implementation

Date: 2026-10-09
Model: gpt-5.6-luna
Settings: default
Plan: doc/decisions/2026-10-09-stage03-pkn-plan.md
Commit(s): add after committing
Transcript: not saved

## Prompt actually used

> implement doc/decisions/2026-10-09-stage03-pkn-plan.md

## What it produced

- `src/pkn.py`: cached resource retrieval, P1-P10 filtering, COSMOS grammar
  parsing, HGNC and metabolite mapping, signed BFS, rewiring null and helpers.
- `scripts/analysis/04_build_pkn.py`: stage-03 orchestration and figures.
- `scripts/setup/run_stage03.sh`: tests followed by the stage-03 run.
- `tests/test_pkn.py`: synthetic tests for P2-P10.
- Stage-03 documentation and provenance updates.

## What was changed by hand afterwards

The implementation records the observed mismatch between the NetworkCommons
and OmniPath-hosted COSMOS endpoints and uses NetworkCommons as the registered
primary when both are available.

## Outputs stamped

Stage-03 outputs are under `results/ionescu_corneto/03_pkn/`; the provenance
sidecar records resource checksums, settings, seeds and the code commit.
