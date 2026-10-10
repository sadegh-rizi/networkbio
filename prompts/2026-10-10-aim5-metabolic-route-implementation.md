# Aim-5 metabolic route implementation

Date: 2026-10-10
Model: gpt-5.6-luna
Settings: default
Plan: doc/decisions/2026-10-10-aim5-metabolic-route-plan.md
Commit(s): pending
Transcript: not saved

## Prompt actually used

> After committing Stage 04, create a separate branch and implement the confirmed Aim-5 metabolic-route plan. Use canonical-only repair and the degree-weighted metabolism-only path as the primary readout. Implement the canonical Recon3D step table, KEGG provenance check, COSMOS gap audit, cofactor mapping/removal, repaired PKN, ordered waypoint test, path variants, type-preserving 100-replicate cholesterol-to-TF null, figure, tests, runner, and provenance outputs. Do not modify source data in place or claim a route is recovered if validation fails.

## What it produced

- `src/metabolic_route.py`: canonical step extraction, gap audit, cofactor mapping, repaired PKN, route/path variants, KEGG verification, and typed null.
- `src/pkn.py`: type-preserving signed rewiring helper.
- `scripts/analysis/04b_aim5_metabolic_route.py` and `scripts/setup/run_aim5_route.sh`.
- `tests/test_metabolic_route.py` and type-preserving rewiring coverage in `tests/test_pkn.py`.
- Confirmed Aim-5 plan/state documentation and structural result outputs under `results/ionescu_corneto/03_pkn/aim5_route/`.

## What was changed by hand afterwards

- The full run showed 19 `present_but_unconnected` canonical rows and a failed final restricted-route waypoint; this was recorded as a residual limitation rather than hidden or reclassified as a gap.
- The 100-replicate type-preserving null is run on the primary COSMOS PKN, matching the registered stage-03 null family; repaired-PKN path variants remain reported separately.

## Outputs stamped

`results/ionescu_corneto/03_pkn/aim5_route/aim5_route.provenance.json` records the inputs, settings, code hashes, and output hashes.
