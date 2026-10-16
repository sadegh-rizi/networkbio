# Stage 06 addendum D9: activity threshold at eps, and in-process FBA readouts

Date: 2026-10-15
Status: confirmed by user (2026-10-15, in chat); implemented (uncommitted)
Amends: `2026-10-13-imat-fba-plan.md` steps 5 and 7, and `2026-10-14-stage06-gpr-and-checks-addendum.md`

## Why

Run 4 of `scripts/analysis/08_metabolic_context.py` finished on 2026-10-15 (22:36). All iMAT
solves were accepted (open: gap 2.6%; defined_chol_plus: 3.5%; defined_chol_minus: 3.4%). The
readouts built from them are affected by a threshold defect.

| Defect | Evidence | Effect |
|---|---|---|
| D9. Activity and direction use different tests at eps | `gem.context_reactions`: `active = |v| >= eps`, `direction = v > eps or v < -eps`. iMAT sets "on" reactions to exactly eps, and the solver returns either `0.01` or `0.0099999999999909`. In `imat/defined_chol_plus/fluxes.tsv` (21,200 rows): 580 rows at 0.0099999… are called inactive, and 2,414 rows at exactly 0.01 are active with direction `inactive`. Only 2,722 rows are clearly above eps | `directional_reaction_nodes` keeps mainly loop reactions at ±1000, so `context_cosmos/pkn.tsv` (25,075 edges) is invalid. It has lost DHCR7 (`Gene8192`), DHCR24 (`Gene4425`) and ER → cytosol `CHSTEROLtrc`, leaving cholesterol fed only by an uptake and a mitochondrial transport step. The R1 differential/enrichment tables and the R2–R4 active-context bounds use the same test and must be recomputed |
| Runtime | The readouts (about 250 LPs) took about 6 h (16:39 → 22:35). Each LP ran in its own process: start-up, Recon3D load and graph build took about 25 s, while the LP itself took 1–2 s | Rerunning the readouts would cost another evening |

## Changes

1. `src/gem.py`: `activity_threshold(eps) = eps − 1e-6` (`ACTIVITY_TOLERANCE`) and
   `is_active(values, eps)`. `context_reactions` uses this threshold for both activity and
   direction, so every active reaction has a direction and every inactive one has none.
2. `scripts/analysis/08_metabolic_context.py`: the step-flux check in `_sterol_flux_ranges` and the
   activity columns in `_signed_vs_flux` use `is_active` instead of `> EPS`.
3. `src/gem_worker.py`: `solve_in_process(task, directory, resume=, root=)` solves FBA tasks in the
   calling process with a cached model (one per model, HGNC and demand-reaction combination).
   - Tasks, signatures, the resume rule and the output files are the same as with
     `solve_in_subprocess`; `result.json` gains `"process": "in_process"`.
   - Checked 2026-10-15: CORNETO's `MultiSampleFBA.build_many` does not modify the graph, so
     successive builds with different bounds give the bound-specific optimum. The edge count is
     re-checked after every solve.
   - A failure in the process (for example HiGHS refusing a second thread-pool size) falls back
     to the subprocess worker.
   - iMAT MILPs always use the subprocess, which keeps their memory isolated.
4. `08_metabolic_context.py`: FBA readouts and checks run in process by default.
   `--fba-subprocess` restores one process per FBA solve. Provenance records
   `activity_threshold` and `fba_in_process`.
5. `tests/test_gem.py` (19 tests, all passing in the sandbox environment, CORNETO 1.0.0rc8):
   - activity and direction at `0.01`, `0.0099999999999909`, `-0.01` and `-0.0099999999999909`;
   - directional COSMOS nodes kept for reactions at eps;
   - in-process FBA with three successive bound sets on one cached model, compared with the
     subprocess result, plus resume;
   - `solve_in_process` refuses iMAT;
   - the earlier toy solves now use the worker's thread count.

## What is reused and what is recomputed

- **Reused with `--resume`, signatures unchanged:**
  - the three iMAT solves;
  - the model checks;
  - medium support.
- **Recomputed, because the active-context bounds change:**
  - flux ranges;
  - cholesterol and simvastatin tasks;
  - `reaction_activity`, `differential_reactions` and `subsystem_enrichment`;
  - `signed_vs_flux`;
  - `context_cosmos/pkn.tsv`.
- **Expected runtime:** minutes rather than hours.

## Rerun

```bash
# keep the D9-affected tables for comparison
mkdir -p results/ionescu_corneto/06_validation/metabolic_context_failed_2026-10-15_d9
cp -r results/ionescu_corneto/06_validation/metabolic_context/{summary,context_cosmos} \
      results/ionescu_corneto/06_validation/metabolic_context_failed_2026-10-15_d9/
RUN_TESTS=1 bash scripts/setup/run_stage06_gem.sh --resume 2>&1 | tee logs/stage06_run5.log
```

The log should show `iMAT …: resumed` for all three media, and `(in process)` for the readout
solves.

## Validation before interpretation

- `context_cosmos/filter_log.tsv`: `n_active_reaction_nodes` is larger than in run 4 (1,840).
- The pruned PKN contains at least one reaction node that produces `Metab__HMDB0000067_c` from
  synthesis (`Gene8192__DHCR7` or `Gene4425__DHCR24`) whenever those reactions are active in
  `reaction_activity.tsv`.
- In `reaction_activity.tsv`, no row is `active = True` with `direction = inactive`.
- Only after these checks: run step 7 with
  `RUN_TESTS=0 bash scripts/setup/run_stage05.sh --resume --cosmos --parts context_cosmos`.

## Open point for interpretation

The direction of a reaction carrying only eps flux is not determined by the data: equally good
iMAT solutions can reverse it. Directional pruning therefore depends on one solution. If the
Golgi → cytosol node (`Gene4449__HDLBP_reverse`) is absent after the fix only because HMR_1916 runs
cytosol → Golgi at eps, report that the shortcut's removal rests on an arbitrary direction. Do not
present it as removed by mass balance.
