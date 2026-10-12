# Aim-5 metabolic route implementation audit

Date: 2026-10-10
Branch: `task/aim5-metabolic-route` (commit `a4dcbb0` and earlier)
Plan: `doc/decisions/2026-10-10-aim5-metabolic-route-plan.md` (confirmed)
Run audited: `results/ionescu_corneto/03_pkn/aim5_route/` and
`03_pkn/cosmos/variants/canonical_repair/` (2026-10-10, 19:56–20:05).

## What is correct

- The repaired PKN adds exactly `Gene9692__FDFT1` (r0170, ER) and
  `Gene4406__FDFT1` (HMR_1465, cytosol) with four edges in the pathway
  direction; nothing is removed. This is the input to stage-04
  `cosmos_repair` and it does not change after the fixes below.
- The canonical step table is built from an explicit ordered step list
  (`_STEP_SPECS`) matched to Recon3D reactions with canonical genes.
- The primary readout (degree-weighted, metabolism-only, no cofactors) on the
  repaired PKN returns the canonical route: 40 steps to `Metab__HMDB0000067_r`,
  16 pathway enzymes, Bloch branch (desmosterol → DHCR24 → cholesterol).
- The type-preserving null and the KEGG gene check are unaffected.

## Findings

1. **Identifier mapping (main finding).** Substrate and product names were
   taken only from Recon3D HMDB annotations. Recon3D has none for several
   late intermediates (44mctr, 7-dehydrodesmosterol, lathosterol,
   desmosterol, 7-dehydrocholesterol), but COSMOS names them by HMDB ID
   (e.g. `Gene5822__DHCR7` takes `Metab__HMDB0000032_r`). As a result,
   18 step rows had missing names, gap_audit classified 19 rows
   `present_but_unconnected`, the restricted-route test failed at
   Lanosterol → Cholesterol on every variant, and the analysis state reported
   this as a resource limitation. It was a mapping artefact.
2. **Missing values treated as present.** `if row.substrate:` and similar
   tests are true for NaN, so "nan" could enter labels, triples and the
   repair guard.
3. **`HMR_1532` misclassified.** GGPP + PPi ⇌ DMAPP + FPP was matched as a
   DMAPP → FPP step because both compounds occur in it; they are on the same
   side, so it is not a conversion step on the route.
4. **Top-1% hub sensitivity.** Canonical intermediates were not protected
   (only cofactors that are intermediates), so ER cholesterol (degree 180)
   was removed and its path failed; the cut-off was taken over distinct degree
   values rather than metabolite nodes, so only a few hubs were removed.
5. **Misleading metric.** `fraction_ordered_canonical_steps_traversed`
   divided matched steps by all 59 rows (both compartments, both branches),
   so a path through the whole pathway scored 0.17.

## Fixes (made by the reviewing agent at the user's request)

| Finding | Change in `src/metabolic_route.py` |
|---|---|
| 1 | New `infer_names_from_cosmos_reactions`: on each side of each canonical reaction, a name is assigned when exactly one non-cofactor Recon3D metabolite is unresolved and exactly one COSMOS metabolite neighbour is unclaimed; repeated to convergence; conflicts raise. `canonical_step_table(..., edges=...)` fills missing names and records `substrate_source` / `product_source`. The script passes the unpruned COSMOS edges. |
| 2 | `_present()` used wherever substrate/product are tested. |
| 3 | Substrate and product must have opposite stoichiometric signs. |
| 4 | Top 1% of metabolite nodes by degree; canonical intermediates and HMGCR protected. |
| 5 | Replaced by `n_waypoints_in_order` and `fraction_waypoints_in_order`: the share of the 11 ordered intermediates (mevalonate → cholesterol) a path visits in order. HMG-CoA is excluded because paths from the HMGCR gene enter the reaction node directly. |

Tests added to `tests/test_metabolic_route.py`: name inference for an
unannotated intermediate, the same-side guard, and protection of the
cholesterol target in the top-1% variant. `tests/test_metabolic_route.py`,
`tests/test_pkn.py` and `tests/test_inference.py` pass (24 tests, Windows
test environment with corneto 1.0.0rc8, cvxpy 1.9.3, highspy 1.15.1).

## In-memory check on the real files (outputs not overwritten)

| Readout | First run | After fixes |
|---|---|---|
| Step rows / unresolved names | 59 / 18 names missing | 58 / 0 |
| Gap audit | 38 forward, 19 unconnected, 2 absent | 56 forward, 2 absent |
| Waypoints, primary and unpruned | fail at FPP → PSQ-PP and Lanosterol → Cholesterol | fail only at FPP → PSQ-PP |
| Waypoints, repaired | (not reported as complete) | all 12 reachable |
| Restricted route, repaired | none | 40 / 42 steps, 16 enzymes |
| Weighted route, repaired | 40 / 42 steps, 16 enzymes | unchanged; 10 of 11 intermediates in order (GPP skipped by the lumped FDPS reaction) |
| Top-1% variant, repaired, ER target | no path | canonical route, 40 steps |
| Repaired PKN | 4 added edges | identical |

## Rerun check (2026-10-11, 09:44–09:49 run)

The user's rerun reproduces the in-memory check: 58 step rows, all names
resolved (11 inferred from COSMOS neighbours), `HMR_1532` excluded, gap
audit 56 present forward / 2 absent, waypoints 11/12 on primary and unpruned
(first failure FPP → presqualene diphosphate) and 12/12 on the repaired PKN;
degree-weighted and top-1% routes on the repaired PKN: 40 steps to ER
cholesterol, 16 pathway enzymes, 10 of 11 intermediates in order, Bloch
branch. The repaired PKN still adds exactly four edges. All code and output
hashes in the sidecar match, but it records `a4dcbb0` with `git_dirty: true`
(fixes not yet committed).

Type-preserving null (100 replicates, primary PKN): from cytosolic
cholesterol, 33/100 rewired networks have a positive route to E2F1 at least
as short as observed (2 steps), 44/100 to SP1, 55/100 to JUN, 45/100 to JUNB,
66/100 to EGR1. The short cholesterol → TF routes are therefore not more
specific than degree- and type-matched rewiring gives.

One further fix: the restricted-route graph was built by iterating a Python
set, so the choice among equal-length paths depended on string-hash
randomisation (the rerun reported a 15-enzyme route via GGPS1 where the
earlier check gave 16). Edges are now inserted in sorted order; the route is
identical under different `PYTHONHASHSEED` values (40 steps, Kandutsch–Russell
branch). The weighted search was already deterministic. Rerun once more after
committing so the sidecar records a clean commit.

## To do

Rerun `bash scripts/setup/run_aim5_route.sh 2>&1 | tee logs/aim5_route.log`
to regenerate the Aim-5 tables and provenance, and correct the prompt record
`prompts/2026-10-10-aim5-metabolic-route-implementation.md`, whose
"residual limitation" note no longer applies.
