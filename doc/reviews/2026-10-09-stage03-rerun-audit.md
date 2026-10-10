# Stage 03 rerun audit (COSMOS parser correction)

Date: 2026-10-09
Run audited: `results/ionescu_corneto/03_pkn/` (17:31–17:51 local),
`logs/stage03.log`. Previous output archived as
`03_pkn_failed_2026-10-09_cosmos_parser/`.
Follows: `2026-10-09-stage03-implementation-audit.md` and the review section
of `doc/decisions/2026-10-09-stage03-cosmos-correction-plan.md`.

## Verdict

The correction worked. The PKNs and mapping are valid inputs for stage 04.
Two items need action: the run was made from uncommitted code, and the COSMOS
path figure is too wide to use.

## Checks

| Check | Result |
|---|---|
| Tests in run script | 52 passed |
| COSMOS primary | 54,571 edges / 26,510 nodes; edge set identical to the independent in-memory build |
| COSMOS unpruned | 81,788 edges / 36,363 nodes; identical edge set |
| Plain symbols starting with X | present (XBP1, XIAP, XPO1, ...) |
| `Metab__` nodes classified `other` | 0 (`other` = 325 complex or unparsed nodes) |
| P6 removed nodes by type (primary) | gene 1,933; enzyme/reaction 5,599; HMDB metabolite 509; numeric metabolite 1,399; model metabolite 88; other 325; transport/pseudo 0 |
| OmniPath primary | 9,870 edges / 3,443 nodes; identical to the previous run |
| Former false HMDB mappings (HMDB0002142, HMDB0002347, HMDB0004984) | listed only as candidates, `mapped_to_pkn = False`; no mapping to numeric nodes |
| TF coverage (CollecTRI kept TFs) | OmniPath 330, COSMOS 349 |
| Provenance inputs | 13 files incl. three regulon summaries; all hashes match |
| Provenance code hashes | match the current working tree |
| Provenance outputs | all 32 hashes match |
| Provenance commit | `e63552d`, `git_dirty: true` (see action 1) |

Mapped metabolite features (unique features mapped to a COSMOS node):

| Study | Labelling class | Previous run | This run |
|---|---|---:|---:|
| ST003328 | not_labelled | 1 | 1 |
| ST003331 | m0_unconfounded | 21 | 23 |
| ST003331 | isotopologue_sum_partial | 9 | 9 |
| ST003331 | m0_confounded | 82 | 85 |
| ST003332 | m0_unconfounded | 12 | 12 |
| ST003332 | m0_confounded | 50 | 76 |

## Aim-5 results

| Leg | Result | Rewired networks with a path as short (20 replicates) |
|---|---|---|
| HMGCR → mevalonate (c, x / e) | 2 / 4 steps | 0/20 |
| HMGCR → cholesterol (6 compartments) | 21–25 steps | 19–20/20 |
| cholesterol_c → E2F1, SP1, JUN | 2 steps | 2/20, 4/20, 5/20 |
| cholesterol_c → EGR1, JUNB | 3 steps | 5/20 each |
| cholesterol, other compartments → TFs | 4–5 steps | 10–20/20 |
| mevalonate → TFs | none | — |
| HMGCR, SREBF1/2, SCAP, INSIG1 → TFs (both PKNs) | none within 8 steps | — |
| any → WT1 | none (not expressed) | — |

Interpretation:
- HMGCR → mevalonate is a specific, direct reaction edge.
- HMGCR → cholesterol is reachable in essentially every rewired network, so
  reachability carries no information. The observed shortest route is also
  not the biosynthesis pathway (open question 1 in the correction plan).
- cholesterol_c → TFs runs through cholesterol→MAPK1/MAPK3 (and AKT1, CAV1)
  edges. With 20 replicates the empirical p-values are about 0.14–0.29
  ((k+1)/21), so the short route is not more specific than expected from
  degree.
- Caveat: the rewiring preserves degree and sign but not node type, so swaps
  can join metabolites to metabolites or genes to metabolites directly. The
  COSMOS null is therefore approximate.

## Actions

1. **Commit the code used for this run.** The sidecar records `e63552d` with
   `git_dirty: true`. The recorded code hashes equal the current files, so
   committing them unchanged makes the run reproducible. Record the new commit
   hash in `doc/analysis-state.md` next to this run.
2. **Fix the COSMOS path figure.** `_path_figure` sets the width to
   1.4 inches per node; the COSMOS manuscript PNG is 82,443 × 4,037 px and
   cannot be used. Options: cap the width and wrap layers, draw one panel per
   leg, or show only the highlighted path edges. Figure-only change; the
   tables are unaffected.
3. Decide open question 1 (metabolic-leg readout) before stage 04 uses it.
