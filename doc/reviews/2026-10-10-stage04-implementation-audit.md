# Stage 04 implementation audit (before the first run)

Date: 2026-10-10
Branch: `task/stage04-inference` (uncommitted files dated 2026-10-10)
Plan: `doc/decisions/2026-10-10-stage04-inference-plan.md` (confirmed)
Scope: `src/inference.py`, the new functions in `src/activities.py`,
`scripts/analysis/05_corneto_inference.py`, `scripts/setup/run_stage04.sh`,
`tests/test_inference.py`. The stage-05 script was not audited.
No stage-04 output exists yet.

## Checked against the installed packages (`.venv`, corneto 1.0.0rc8, cvxpy 1.9.3, highspy 1.15.1)

| Point | Result |
|---|---|
| `CarnivalFlow(lambda_reg).build_from_data(G, data)` | as in the passing toy; defaults `max_flow=1000`, `indirect_rule_penalty=0` |
| Error objective | `|value| × (1 − state × sign(value))` per measurement (`_carnival.py`), so unscaled scores act as weights as planned |
| Objective names | `error_<sample>_<i>` (toy record), so `weighted_fitting_error` sums the right terms |
| Time limit | CORNETO maps `max_seconds` to HiGHS `time_limit` |
| `mip_rel_gap`, `threads`, `random_seed` | passed by cvxpy through `setOptionValue`; valid HiGHS options |
| Time-limited solutions | HiGHS `kTimeLimit` → cvxpy `user_limit`, which is in `SOLUTION_PRESENT`, so values are kept |
| MIP gap | read from `solver_stats.extra_stats` = HiGHS `info.mip_gap` |
| `sample_alternative_solutions` | exists; keyword arguments match |
| `moderated_t_statistics` (3 vs 4) | correct when given a whole feature set |
| `leave_one_line_out_ulm`, `select_stable_top_tfs` | match the plan and an independent recomputation (all top-50 TFs sign-stable for both PKNs) |

## Findings

1. **Blocker: the run stops before solving anything.** `main()` requires
   `03_pkn/cosmos/variants/canonical_repair/pkn.tsv`, which comes from the
   Aim-5 route plan. That plan is not implemented, so a `FileNotFoundError`
   is raised. Either implement the Aim-5 plan first, or skip `cosmos_repair`
   with status `not_run_missing_pkn` and run it later.
2. **Metabolite moderated t uses the wrong prior.** `_metabolite_table`
   estimates the empirical-Bayes prior on the selected features only
   (ST003331: 23, ST003332: 12) and on one feature for cholesterol
   (d0 = ∞). Estimated on all complete features of each study, then subset:
   d0 changes from 4.93 to 2.90 (ST003331) and from 7.37 to 1.94 (ST003332),
   metabolite t-values change by a median of 6% with no sign changes, and
   cholesterol t changes from 5.09 to 6.41 (prior from 642 lipid species).
   The values are measurement weights, so fix this before running.
3. **Node states are reconstructed from edges instead of read from the
   model.** `_fit_table` derives states from selected edges and raises an
   error on any conflict, and `_node_table` uses the first edge seen. The
   error objective uses CORNETO's `vertex_value`. Read `problem.expr.vertex_value`
   (registered by CarnivalFlow) for node states, sign agreement and the node
   table, so the reported fit matches the optimised objective and a conflict
   cannot crash a run partway through the grid.
4. **Aim-5 readout is mislabelled.** `n_canonical_reaction_nodes` counts
   every selected `Gene*` node, not the canonical pathway reactions, and
   `cholesterol_nodes_active` includes measured cholesterol nodes with
   state 0. Count only canonical reactions (the Aim-5 step table) and only
   nodes with state ≠ 0.
5. **Unreachable measurements are not reported.** CORNETO prunes measured
   nodes that no input can reach; they leave the error term but count as
   mismatches in `sign_agreement_ratio`. Report the number of unreachable
   measurements per run and exclude them from the ratio (or state that they
   are included).
6. **Provenance is incomplete.** Files read but not hashed: the
   ST003331/ST003332/ST003328 `features.tsv` and line tables, the CollecTRI
   cache file, and `03_pkn/cosmos/primary/nodes.tsv` (used for degrees).
7. **`inputs.tsv` does not match what is passed to CORNETO** for
   `omnipath_lines` (measured TFs marked as inputs because the lookup uses the
   `PMS_vs_Ctrl` key) and for COSMOS runs (measured metabolites not
   considered). Derive the table from `make_corneto_data`'s returned inputs.
8. **Run time.** Worst case 2 × 6 × 1,800 s + 3 × 6 × 3,600 s ≈ 24 h on one
   machine, and COSMOS memory use is untested (about 4,700 inputs, 54,571
   edges). Add `--runs` and `--lambdas` options for a timing pilot (registered
   settings unchanged), and run OmniPath first.
9. Minor: no prompt record under `prompts/` for this implementation;
   `sample_edges` (stage 05) fails on reshape when no alternative solution is
   accepted.

## Verdict

The CORNETO calls, solver settings, λ rule, acceptance rule and TF selection
are implemented as planned and match the installed library behaviour. Do not
run yet: finding 1 stops the run, and findings 2–4 would produce wrong
weights or misreported readouts.

## Resolution (2026-10-10, fixed by the reviewing agent at the user's request)

| Finding | Change |
|---|---|
| 1 | `cosmos_repair` is skipped and written to `summary/run_status.tsv` as `not_run_missing_pkn` while the repaired PKN does not exist; the other runs proceed. It runs once the Aim-5 route plan has produced the variant. |
| 2 | Metabolite and cholesterol moderated t now use the prior from all complete features of each study (`n_features_prior` recorded in `inputs/measured_metabolites.tsv`). |
| 3 | Node states, sign agreement and `nodes.tsv` come from CORNETO's `vertex_value` (`inference.vertex_states`, `inference.node_table`). Edge/state disagreements are counted (`n_edge_state_conflicts`) instead of raising. |
| 4 | `summary/aim5_route.tsv` counts active reaction nodes whose genes are in the pre-registered canonical list (`CANONICAL_CHOLESTEROL_GENES`, from the Aim-5 plan), and lists cholesterol nodes, HMGCR and targets only with state ≠ 0. |
| 5 | `fit.tsv` reports `n_missing_from_pkn`, `n_unreachable`, `n_reachable`; the sign-agreement ratio uses reachable measurements only; `nodes.tsv` has a `reachability` column. |
| 6 | Provenance hashes the metabolomics features and line tables, the CollecTRI cache file, the COSMOS node table and, when used, the repaired PKN and its node table. |
| 7 | `inputs/inputs.tsv` is built from the inputs actually written to the CORNETO data, per run, with `excluded`, `measured` and `used_as_input` columns. |
| 8 | `--runs` and `--lambdas` (registered grid values only) for a timing pilot; a pilot writes to `results/ionescu_corneto/04_inference_pilot/` and records `pilot: true`. A pilot without λ = 0 selects no λ by design. |
| 9 | `sample_edges` returns an empty table when no alternative is accepted. The prompt record is `prompts/2026-10-10-stage04-audit-fixes.md`. |

Tests: `tests/test_inference.py` gains two tests (a measured source node is
not an input; states from `vertex_value`, an unreachable measurement is
reported and excluded from the ratio). With corneto 1.0.0rc8, cvxpy 1.9.3 and
highspy 1.15.1 (Windows test environment, decoupler not installable there):
`tests/test_inference.py` and `tests/test_pkn.py` pass (19 tests).
`tests/test_activities.py` needs the WSL venv.

Real-scale smoke solve (not a stage-04 result): OmniPath primary, the 50
TFs selected by the E2 rule from the 2026-10-09 provisional contrast table,
λ = 0.1, HiGHS 4 threads. Preprocessed graph 1,626 nodes / 6,850 edges, 823
inputs; optimal in 52 s; all 50 TFs fitted with the correct sign using 65
edges; 0 unreachable measurements, 0 edge/state conflicts.

Observation for interpretation: with 823 unknown-sign inputs the error at
λ = 0.1 is already 0, so the 10% rule reduces to "largest λ with the same
error as λ = 0". A perfect fit is therefore expected and carries little
information on its own; the label-permutation and rewired-PKN controls
(stage 05) are what test whether the network reflects the PMS contrast.
