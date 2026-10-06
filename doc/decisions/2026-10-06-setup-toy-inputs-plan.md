# Environment, CORNETO toy problem and stage-00 input tables

Date: 2026-10-06
Status: implemented (run in WSL 2026-10-06; not yet committed)

Covers steps 0-2 of `2026-09-29-corneto-ionescu-inscs-plan.md` only. No
preprocessing, activity estimation, PKN or inference decisions are made here.

## Question

Can the project environment run CORNETO with an open-source MILP solver on
a signalling + metabolism graph, and what exactly is in each deposited layer
(samples, lines, replicates, identifiers) before any analysis choice?

## Data available

GSE297192 raw counts and series matrix (Park 2025); Metabolomics Workbench
ST003328, ST003331, ST003332 (Ionescu 2024). See `doc/datasets.md`.
Not used here: ST003330 (extracellular lipidomics), single-cell, WGBS, ATAC.

## Unit of observation

Cell line. Replicates are kept as rows in the stage-00 tables and counted
per line; nothing is averaged here.

## Method

0. Environment: `.venv` in the repo root from `requirements.lock.txt`
   (Python 3.12, corneto 1.0.0rc8 with HiGHS and SCIP, decoupler 2.2.0,
   omnipath 1.0.12). `networkcommons` dropped: it pins graphviz<0.18, which
   conflicts with corneto 1.0.0rc8.
1. Toy problem (`scripts/analysis/00_toy_corneto.py`): a 13-node signed
   graph with a COSMOS-style metabolite node (`Metab__HMDB0000067_c`) and
   decoy edges; `CarnivalFlow(lambda_reg=0.1)`, HiGHS, 60 s limit. Checks:
   planted edges recovered exactly, metabolite edges used, inhibition sign
   correct, input with unknown sign (value 0) inferred, two identical samples
   in one multi-sample problem give the same edges.
2. Stage-00 tables (`scripts/analysis/01_build_inputs.py`, functions in
   `src/inputs.py`): sample sheets with group, treatment, line and replicate;
   feature tables with deposited IDs (KEGG where given) and a 13C-isotopologue
   flag; values as deposited; replicates per line; ID coverage; exclusions.

## Open questions

None for these steps. Questions for later stages are in the main plan.

## Files to change

`requirements.txt`, `requirements.lock.txt` (new), `scripts/setup/create_venv.sh`,
`scripts/setup/run_setup_and_stage00.sh`, `scripts/analysis/00_toy_corneto.py`,
`scripts/analysis/01_build_inputs.py`, `src/inputs.py`, `tests/test_inputs.py`,
and the docs (`analysis-state.md`, `datasets.md`, `external-tools.md`,
`context/*`, `decisions/README.md`, `prompts/`).

## Outputs

`results/corneto_toy/04_inference/{toy_edges.tsv,toy_checks.json}`;
`results/ionescu_corneto/00_inputs/` (see the script docstring).

## Validation

`bash scripts/setup/run_setup_and_stage00.sh` ends with `ALL STEPS FINISHED`:
pytest passes, every toy check prints PASS (exit code 0), and the line counts
printed by `01_build_inputs.py` are 3 Ctrl / 4 PMS for the metabolomics and
lipidomics layers and 3 / 3 for RNA-seq. A toy failure falsifies the
assumption that CORNETO treats metabolite nodes like any other node.
