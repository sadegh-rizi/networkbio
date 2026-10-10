# Stage 04 (and pre-registered stage-05 controls): CORNETO inference

Date: 2026-10-10
Status: proposed

Answers main-plan questions Q4 (inputs, controls), Q7 (λ rule), Q8 (solver
limits) and the COSMOS run direction for this project. All choices below
are fixed before any real-data CORNETO run.

## Question

Which signed subnetwork of each prior-knowledge network does CORNETO
(`CarnivalFlow`) select to explain the PMS-vs-Ctrl TF activity changes, and
for COSMOS also the metabolite changes? How stable is it, is it more than
the PKN topology and hub structure would give, and does it contain the
HMGCR → cholesterol → TF route without that route being an input?

## Data available

- TF activity contrasts, stage 02 (`results/ionescu_corneto/02_activities/`):
  CollecTRI (primary), CollecTRI without complexes, DoRothEA A–C; moderated t
  (constant prior), Welch, logFC. Unblinded on 2026-10-09 after E1–E4 were
  confirmed.
- PKNs, stage 03 (commit `c7f3719`): OmniPath primary (9,870 edges / 3,443
  nodes) and variants `nc_rule`, `signor_only`, `unpruned`; COSMOS primary
  (54,571 / 26,510) and `unpruned`; COSMOS `canonical_repair` from
  `2026-10-10-aim5-metabolic-route-plan.md` once that is run.
- Metabolite line tables, stage 01 (3 Ctrl vs 4 PMS lines): ST003331
  intracellular `line_log2.tsv`, ST003332 extracellular
  `line_log2_vs_blank.tsv`, ST003328 lipidomics `line_log2.tsv` (untreated
  columns only). Feature → COSMOS node map: `03_pkn/mapping/`.
- Not available: metabolite contrasts (no stage computes them yet), receptor
  or ligand data, phosphoproteomics, line matching across deposits.

## Unit of observation

Lines. RNA-seq: 3 Ctrl vs 3 PMS (C1–C6). Metabolomics and lipidomics:
3 Ctrl vs 4 PMS, with line letters not matched to the RNA-seq lines. All
measurements enter CORNETO as one group-level contrast "sample" (primary)
or as six per-line RNA samples (secondary multi-sample run). Findings are
hypothesis-generating.

## Change to a confirmed decision

- **E1 deferred (user, 2026-10-10).** Limma-trend is not implemented now.
  Stage 04 uses the constant-prior moderated t (current stage-02 primary).
  Limma-trend becomes a later sensitivity, and the addendum status is
  updated accordingly. E2–E4 apply unchanged.

## Method

### Inputs to CORNETO

1. **Measured TFs (E2).** Universe: CollecTRI TFs kept after tmin = 5 that
   are nodes of the PKN being run. Rank by |ULM score| on the primary
   contrast; take the top k = 50. Compute six leave-one-line-out contrasts
   (moderated t on the remaining 2 vs 3 lines, then ULM with the same
   regulon) and keep a TF only if its sign agrees in all six. Record kept and
   dropped TFs. Values passed to CORNETO: the ULM score (signed, unscaled).
   No TF is added or removed by hand (E3). AP1 and NFKB are not
   measurements in the primary run (E4).
2. **Measured metabolites (COSMOS runs only).** Per study, moderated t
   (constant limma prior, generalised to 3 vs 4) of PMS − Ctrl on line
   means. Primary features: ST003331 `m0_unconfounded`, ST003332
   `m0_unconfounded` (extracellular, mapped to `_e` nodes), and free
   cholesterol from ST003328 untreated. One node per feature: intracellular
   features to the `_c` node if it exists, otherwise the compartment with
   the highest degree; record the choice. Value = moderated t. Sensitivity:
   add `isotopologue_sum_partial` and `m0_confounded` features.
3. **Inputs (unknown perturbation).** All source nodes (in-degree 0) of the
   PKN, entered with value 0 so CORNETO infers their sign and may leave them
   unused. Never inputs: HMGCR, SREBF1, SREBF2, SCAP, INSIG1, E2F1, EGR1,
   WT1, SP1, JUN, JUNB, and the mevalonate and cholesterol nodes (INSIG1 is a
   COSMOS source node and is removed). Counts with the stage-03 PKNs:
   OmniPath 823 source nodes; COSMOS 4,713 (1,017 genes, 3,684 metabolites,
   12 reaction nodes).
4. **Graph pruning.** Only CORNETO's own reachability preprocessing
   (nodes on some input → measurement path). No path-length cap: the
   canonical cholesterol route is 40–42 edges long and a cap would exclude it
   by construction. Record node and edge counts after preprocessing.

### Runs (stage 04, `04_inference/`)

5. **R-omni.** OmniPath primary, TF measurements only.
6. **R-cosmos-tf.** COSMOS primary, the same TF measurements only.
7. **R-cosmos-multi.** COSMOS primary, TFs + metabolites, one joint problem
   (no separate COSMOS forward/backward runs; main-plan Q3).
8. **R-cosmos-repair.** As 7 on the `canonical_repair` variant (tests
   whether the cholesterol route is used once the missing FDFT1 step exists).
9. **R-omni-lines (secondary).** OmniPath primary, multi-sample problem with
   the six per-line TF activity vectors (stage-02 `tf_activity_per_line.tsv`,
   selected TFs only) as six samples with the shared edge penalty.

### Settings

10. **Method.** `CarnivalFlow` from the `corneto` import path used in
    `00_toy_corneto.py`. Before any real run, a smoke step prints
    `help(CarnivalFlow)`, the `solve` keyword names for HiGHS time limit, gap,
    threads and seed, and confirms `corneto.methods.sampler.sample_alternative_solutions`
    exists in 1.0.0rc8. If an API differs from this plan, stop and report.
11. **Solver.** HiGHS, 4 threads, seed 20261010, relative MIP gap 0.01.
    Time limit 1,800 s per OmniPath problem, 3,600 s per COSMOS problem. A
    solution is accepted if optimal, or if the time limit is reached with a
    gap ≤ 0.05; otherwise the run is `unsolved` and is reported, not used.
    SCIP is a cross-check on R-omni at the selected λ only.
12. **λ rule (Q7).** Grid λ ∈ {0, 0.01, 0.03, 0.1, 0.3, 1}. For each run,
    choose the largest λ whose weighted fitting error is ≤ 1.10 × the error
    at λ = 0. Report every grid point (error, number of edges and nodes,
    sign-agreement ratio of measured nodes). The same rule applies to every
    run; λ is not tuned per result.

### Readouts

13. Selected edges and node signs per run; how many measurements are fitted
    with the correct sign; which inputs are used.
14. **Aim 4 (early vs late).** Jaccard of selected protein nodes and edges
    between R-omni, R-cosmos-tf and R-cosmos-multi. Which metabolites are
    explained, and which protein nodes enter only when metabolites are added.
15. **Aim 5.** In R-cosmos-multi and R-cosmos-repair: is HMGCR selected; how
    many canonical reaction nodes are selected; are cholesterol nodes active
    and with which sign; is there a selected route from a cholesterol node to
    E2F1, EGR1, SP1, JUN or JUNB. Report as structural plus data-fit
    statements.

### Robustness and controls (stage 05, `05_robustness/`, all at the λ chosen in 12)

16. **Alternative optima.** `sample_alternative_solutions` on each accepted
    stage-04 run: 30 samples, `percentage=0.03`, `scale=0.03`,
    `rel_opt_tol=0.05` (CORNETO multi-condition tutorial values). Edge and
    node frequencies.
17. **Leave-one-line-out.** For each of the six LOO contrasts, redo TF
    selection (step 1) and R-omni. For COSMOS, R-cosmos-multi only if the
    single R-cosmos-multi run solves in under 1,800 s.
18. **Analysis-choice grid (R-omni).** One factor at a time: statistic
    (Welch, logFC), regulon (DoRothEA A–C, CollecTRI without complexes — this
    is also the AP1/NFKB sensitivity of E4), k (25, 100), PKN variant
    (`nc_rule`, `signor_only`, `unpruned`), λ (neighbouring grid points).
    Metric: edge and node Jaccard with the primary R-omni network. Core
    network = edges present in ≥ 80% of steps 16 and 17 combined.
19. **Label permutation (R-omni).** All 20 assignments of 3 of the six RNA
    lines to "PMS" (the true split, its mirror and 18 others). For each,
    recompute the moderated t, ULM, TF selection and R-omni. Statistics:
    weighted fitting error, network size, Jaccard with the true network.
    Report the rank of the true split; the smallest attainable one-sided p
    is 1/10 after pairing mirrors, which is stated with the result.
20. **Rewired PKN.** Ten degree- and sign-preserving rewirings of OmniPath
    primary (stage-03 function, seed base 20261010) and, if R-cosmos-multi
    solves, five type-preserving rewirings of COSMOS primary (function from
    the Aim-5 route plan). Run R-omni / R-cosmos-multi on each. Compare
    fitting error and whether the selected nodes overlap the real network.
21. **Hub bias.** Logistic regression of node selection (any sample in
    step 16) on log total degree in the real R-omni run, and the same on the
    rewired runs. Report the odds ratio per doubling of degree.

## Open questions

1. Inputs: all source nodes with unknown sign (recommended, step 3), or
   receptors only (OmniPath intercell annotation, expressed genes)?
2. Include free cholesterol from the lipidomics as a COSMOS measurement in
   the primary run (recommended; it is a measurement, not an input), or only
   in a sensitivity run?
3. λ rule: 10% error tolerance (recommended), or CORNETO tutorial's
   overall-agreement metric?
4. Time limits (1,800 / 3,600 s) and the 0.05 acceptance gap: acceptable
   given that the stage-05 grid multiplies them (about 40 OmniPath problems)?
5. Secondary multi-sample run (step 9): in scope now, or after stage 05?

## Files to change

- `src/inference.py` (new): CORNETO data construction, input/measurement
  tables, solve wrapper with acceptance rule, λ grid and selection, solution
  extraction, sampling wrapper.
- `src/activities.py`: generalise `moderated_t_statistics` to unequal group
  sizes (3 vs 3 behaviour unchanged, tested), add leave-one-line-out
  contrast and E2 selection functions.
- `scripts/analysis/05_corneto_inference.py` (stage 04) and
  `scripts/analysis/06_corneto_robustness.py` (stage 05) (new).
- `scripts/setup/run_stage04.sh`, `scripts/setup/run_stage05.sh` (new).
- `tests/test_inference.py` (new): toy PKN with a planted route and decoys
  (adapted from `00_toy_corneto.py`); λ rule on a synthetic error curve;
  acceptance rule; excluded inputs never appear as inputs.
- `tests/test_activities.py`: moderated t 3 vs 4 against an independent
  computation; LOO contrast shapes; E2 selection on a synthetic table.
- `doc/analysis-state.md`, `doc/decisions/README.md`,
  `doc/decisions/2026-10-09-stage02-trend-and-tf-selection-addendum.md`
  (E1 deferral note).

## Outputs

`results/ionescu_corneto/04_inference/`:
- `inputs/{measured_tfs,measured_metabolites,inputs,tf_selection_loo}.tsv`
- `<run>/lambda_<value>/{edges,nodes,fit}.tsv` with `<run>` ∈
  {`omnipath`, `cosmos_tf`, `cosmos_multi`, `cosmos_repair`, `omnipath_lines`}
- `summary/{lambda_grid,selected_lambda,run_status,aim4_overlap,aim5_route}.tsv`
- `figures/` network figures per run at the selected λ.
- `04_inference.provenance.json`

`results/ionescu_corneto/05_robustness/`:
- `sampling/<run>/edge_frequency.tsv`, `loo/<line>/…`, `grid/<factor>/<level>/…`,
  `permutation/<assignment>/…`, `rewired/<pkn>/<replicate>/…`
- `summary/{jaccard,core_network,permutation_rank,rewired,hub_bias}.tsv`,
  `figures/`, `05_robustness.provenance.json`

## Validation

- `pytest` in WSL, then `bash scripts/setup/run_stage04.sh 2>&1 | tee logs/stage04.log`
  and later `run_stage05.sh`.
- Works if: the toy test recovers its planted network; every accepted run
  meets the acceptance rule; no excluded node appears as an input; λ is
  chosen by the rule for every run.
- Falsification: (a) the true label split does not rank above most of the
  19 other assignments → the network reflects the PKN and selection more than
  the PMS contrast; (b) rewired PKNs fit as well as the real PKN → topology
  is not informative; (c) the core network is small relative to single
  solutions → single-solution networks should not be interpreted.
- Expected structural result, stated in advance: without the repair, a
  canonical HMGCR → cholesterol route cannot appear in any COSMOS solution
  because the PKN lacks FPP → presqualene diphosphate. In OmniPath, the
  cholesterol axis has no route to the TFs at all.
