# Analysis state: read this first

Last updated: 2026-10-09. Update this file whenever a decision is made or a
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
- **Stages 01-02 revised choices (2026-10-09):** plain rows remain registered
  as M+0 because R0 found no deposit statement that they are total pools;
  ST003332 and ST003328 primary values are unscaled log2 peak areas, the
  moderated t is the primary RNA contrast, and CollecTRI complexes are retained
  with a no-complex sensitivity. The implementation is synthetic-tested, and
  the real-data run of stages 01 and 02 completed on 2026-10-09.
- **HGNC ambiguity rule (2026-10-09):** three Ensembl IDs in the HGNC
  2026-10-06 release have multiple distinct approved symbols. They are excluded
  from the one-to-one RNA-seq mapping, not duplicated or assigned arbitrarily,
  and are recorded in the gene map and QC/exclusion tables.
- **Stage 01 real-data run (2026-10-09):** completed successfully after the
  HGNC rule was implemented. The primary output has 14,490 genes passing the
  pre-mapping CPM check, 13,530 retained symbols, and a valid provenance
  sidecar. The prior incomplete output is archived at
  `results/ionescu_corneto/01_preprocessing_failed_2026-10-09_hgnc_ambiguity/`.
- **Stage 02 real-data run (2026-10-09, `logs/stage01_02.log`):** completed;
  30 tests passed in WSL first. Moderated-t prior d0 = 2.85, s0^2 = 0.043;
  CollecTRI 628 TFs after tmin = 5 (AP1 321 targets, NFKB 466), no-complex
  variant 626, DoRothEA A-C 284.
- **Stage 03 initial run (2026-10-09):** the OmniPath result was valid, but the
  COSMOS result was invalidated by the implementation audit because the parser
  treated NetworkCommons reaction indices as Entrez IDs and relabelled numeric
  model metabolites as HMDB IDs. Its outputs are retained in
  `results/ionescu_corneto/03_pkn_failed_2026-10-09_cosmos_parser/`.
- **Stage 03 corrected real-data run (2026-10-09, commit `c7f3719`,
  `results/ionescu_corneto/03_pkn/`):** the wrapper completed with 52 tests
  passed. OmniPath returned 72,272 rows; its primary P1-P6 PKN has 9,870 edges
  and 3,443 nodes. The corrected NetworkCommons COSMOS primary P7/P5/P6 PKN
  has 54,571 edges and 26,510 nodes; its unpruned variant has 81,788 edges
  and 36,363 nodes. Numeric model metabolites remain numeric, pseudo-reaction
  nodes are retained, and the sidecar hashes all 13 stage inputs and the
  implementation code. The OmniPath-hosted COSMOS file was also fetched, but
  the two endpoints differ after cleanup (81,788 versus 65,138 edges) and are
  not treated as interchangeable. Recon3D and 205 Metabolomics Workbench KEGG
  lookups were cached with checksums.
- **Aim-5 structure check (2026-10-09, commit `c7f3719`):** the corrected
  table uses three registered legs. No source-TF-to-TF positive path exists
  within 8 edges in either primary PKN. HMGCR reaches mevalonate in 2-4 edges
  and cholesterol in 21-25 edges in COSMOS under the cap of 80; cholesterol
  reaches several target TFs within 8 edges, while mevalonate and WT1 do not.
  These are structural results about the filtered PKNs, not biological negative
  findings; the complete leg-labelled pair table and rewiring null are in
  `03_pkn/aim5/`.
- **Stage 03 figure correction (2026-10-09, commit `c7f3719`):** path figure
  width is capped at 24 inches so the COSMOS manuscript PNG is bounded; tables
  and path calculations are unchanged.
- **Audit (2026-10-09, `doc/reviews/2026-10-09-stage01-02-implementation-audit.md`):**
  code matches the plans and the limma/decoupler sources; moderated t and the
  ST003332 blank reference were recomputed independently and match. Open:
  log2(CPM+1) has a strong mean-variance trend (limma-trend variant
  recommended); RNA-seq PC1 (69%) splits Ctrl/PMS perfectly, so check SRA run
  metadata for a group-aligned batch; stage 04 needs a pre-stated TF-selection
  rule (ULM padj counts change threefold with the contrast statistic). The
  audit did not look at the six held-out TFs or AP1.

## Not yet decided

See the open questions and the 2026-10-06 proposed defaults in
`doc/decisions/2026-09-29-corneto-ionescu-inscs-plan.md`.
The main ones: which layers and conditions are inputs, signalling only or
signalling plus metabolism, how the input nodes are chosen without making the
validation circular, the solver, and the course deliverable. The PKN filters
and stage-03 resource choice are settled in
`doc/decisions/2026-10-09-stage03-pkn-plan.md`.

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

- `00_toy_corneto.py` uses `lambda_reg=0.1` and a 12-node, 12-edge graph; it
  says nothing about run time on a real PKN. Its figure and Cytoscape tables
  come from `scripts/analysis/00b_plot_toy_network.py`.
- `networkcommons` could not be installed (pins graphviz<0.18, conflicts
  with corneto 1.0.0rc8); the COSMOS meta-PKN has to be fetched directly.

## Next steps, in order

1. Confirm the remaining open questions of the main plan (Q4 controls, Q6
   PKN filters, Q7 lambda rule, Q8 limits, Q11 deliverable). Q9 is resolved
   for stage-01 preprocessing: mean log2 replicates with the 2-of-3 rule.
2. Stages 01-02 (preprocessing, TF activities): done. Plans:
   `doc/decisions/2026-10-06-stage01-02-preprocessing-tf-activities-plan.md`,
   `doc/decisions/2026-10-09-stage01-02-revision-plan.md` and
   `doc/decisions/2026-10-09-hgnc-ambiguous-mapping-plan.md`; code on
   `task/stage01-02`; real run 2026-10-09; HGNC, CollecTRI and DoRothEA
   resources cached under `data/resources/`. Before stage 04, decide the two
   audit points (limma-trend variant or primary; TF-selection rule) and check
   the SRA run metadata for C1-C6.
3. Stage 03: done in commit `c7f3719`. The corrected PKNs, identifier maps,
   leg-labelled Aim-5 structure table, rewiring null and capped-width figures
   are in `results/ionescu_corneto/03_pkn/`.
4. Stage 04: CORNETO runs, with the controls fixed beforehand and the COSMOS
   endpoint mismatch explicitly reported.
