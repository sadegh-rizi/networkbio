# Analysis state: read this first

Last updated: 2026-10-06. Update this file whenever a decision is made or a
run changes a number. Detail and sources are in `doc/context/`.

## Question

Which subnetwork of a prior-knowledge network does CORNETO select to explain
the PMS-vs-control differences in the Ionescu et al. (2024) iNSC model, and
does it recover the paper's proposed route (HMGCR-driven lipogenic state,
cholesterol-dependent TFs, SASP) without being given it as an input? How
does it compare with the authors' MOFA + GENIE3 analysis?

The effective n is lines: 3 control and 4 PMS (fewer for some layers).
Findings are hypothesis-generating.

## Settled

- **Method focus: CORNETO** (Nature Machine Intelligence 2025), chosen by
  the user over the multilayer, hypergraph and controllability options.
- **Source dataset:** Ionescu et al., *Cell Stem Cell* 2024
  (DOI 10.1016/j.stem.2024.09.014). The Park et al., *Neuron* 2025 data
  (DOI 10.1016/j.neuron.2025.09.022) are the candidate public transcriptomics
  layer.
- **Simvastatin is a separate condition** (user, 2026-10-06). Only the
  lipidomics (ST003328) has it.
- **Environment exists** (2026-10-06): `.venv` in the repo root, Python 3.12.3
  in WSL (Ubuntu 24.04), from `requirements.lock.txt`; corneto 1.0.0rc8,
  HiGHS (highspy 1.15.1) and SCIP available. Installed versions:
  `logs/pip_freeze_2026-10-06.txt` (gitignored).
- **CORNETO handles a metabolite node like any other node**: the toy problem
  passed all five checks in WSL (`results/corneto_toy/04_inference/toy_checks.json`).
  So a COSMOS-style signalling + metabolism PKN can go into `CarnivalFlow`.

## Not yet decided

See the open questions and the 2026-10-06 proposed defaults in
`doc/decisions/2026-09-29-corneto-ionescu-inscs-plan.md`.
The main ones: which layers and conditions are inputs, signalling only or
signalling plus metabolism, how the input nodes are chosen without making the
validation circular, the PKN, the solver, and the course deliverable.

## Known unknowns

- Whether the Park and Ionescu donors are the same 7 lines: not checked
  line by line (control names differ between the papers).
- Deposits checked against the repositories on 2026-09-29; full tables in
  `doc/data-inventory/data_inventory_ionescu_park.xlsx`. In short:
  Metabolomics Workbench ST003328 (42), ST003330 (24), ST003331 (21, 24 h
  13C-glucose only), ST003332 (24) are public, and only the lipidomics
  (ST003328) has an SV arm. MassIVE MSV000095283 (proteomics, secretome) is
  **private**. There is no Ionescu RNA-seq deposit. GSE297192 has 25 samples:
  C1-C3 Ctrl and C4-C6 PMS baseline (split taken from the authors' code, not
  GEO), plus 19 conditioned-medium samples. Park's baseline 3 vs 3 sc/snRNA-seq
  is not in GEO. The processed tier (`01_processed.sh`, 3.85 GB) was
  downloaded on 2026-09-29 and checked; see the download log in
  `doc/datasets.md`.
- Whether to run signalling->metabolism and metabolism->signalling as
  separate problems (as COSMOS does) or as one. The toy problem shows one
  graph with metabolite nodes works; scale is untested.
- The course deadline, deliverable format and grading criteria.

## Stage-00 facts (2026-10-06, `results/ionescu_corneto/00_inputs/`)

- Lines per layer (`summary/replicates_per_line.tsv`): lipidomics 3 Ctrl +
  4 PMS, each untreated and SV, 3 replicates per line and condition;
  intracellular (ST003331) and extracellular (ST003332) metabolomics 3 Ctrl +
  4 PMS, 3 replicates; RNA-seq C1-C6, 3 Ctrl + 3 PMS, one library each.
- Line letters are reused differently in each deposit (ST003328: per group;
  ST003331 and ST003332: subjects A-G across groups; ST003332 controls are
  A, B, D). Lines are **not** matched across deposits.
- Identifiers (`summary/feature_id_coverage.tsv`): ST003331 has KEGG IDs for
  all 158 unlabelled features plus 23 13C-isotopologue rows; ST003332 has
  KEGG IDs for all 136; the lipidomics (642 species) has names and lipid
  classes only, no database IDs. It includes free cholesterol and 30
  cholesteryl esters (`summary/lipid_groups.tsv`).
- Excluded at this stage (`summary/exclusions.tsv`): RNA-seq C7-C24
  (conditioned-medium experiment), GEO sample C18_2 (not in the count
  matrix), the 3 no-cell medium blanks in ST003332.

## Known broken or stale

- `00_toy_corneto.py` uses `lambda_reg=0.1` and a 13-node graph; it says
  nothing about run time on a real PKN.
- `networkcommons` could not be installed (pins graphviz<0.18, conflicts
  with corneto 1.0.0rc8); the COSMOS meta-PKN has to be fetched directly.

## Next steps, in order

1. Confirm the remaining open questions of the main plan (Q4 controls, Q6
   PKN filters, Q7 lambda rule, Q8 limits, Q9 replicates, Q11 deliverable).
2. Stages 01-02 (preprocessing, TF activities): plan written,
   `doc/decisions/2026-10-06-stage01-02-preprocessing-tf-activities-plan.md`
   (status proposed; 6 open questions to confirm before implementation).
3. Stage 03: build and log the PKN (OmniPath, then COSMOS meta-PKN); map
   KEGG IDs to the PKN's metabolite IDs and count what maps.
4. Stage 04: CORNETO runs, with the controls fixed beforehand.
