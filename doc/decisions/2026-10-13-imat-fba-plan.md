# Stage 06: stoichiometric context models (iMAT/FBA) to fix the signed-metabolism shortcuts

Date: 2026-10-13
Status: proposed, awaiting user confirmation (open questions Q1–Q6 have recommended defaults).
Depends on: stages 01–05 (finished). Code from CORNETO 1.0.0rc8 already in `.venv`.

## Why

In stages 04–05 the metabolic layer is COSMOS's signed graph. That causes four issues; this plan
addresses each one:

| Issue (stage 04/05 evidence) | Cause | How this plan addresses it |
|---|---|---|
| I1. Cholesterol rise explained by `HDLBP → Gene4449__HDLBP_reverse → cholesterol` (Recon3D HMR_1916, Golgi ↔ cytosol cholesterol transport, genes HDLBP or NPC1). The simvastatin control on the primary PKN used a 32-edge mitochondrial detour through steroid sulfatase | Signed edges have no mass balance: an active enzyme "produces" its product without consuming a substrate | iMAT/FBA on the stoichiometric Recon3D model, where cholesterol can only appear by synthesis or uptake with balanced fluxes (steps 4–6) |
| I2. FDFT1 gap fixed only by a targeted two-reaction repair | COSMOS lacks 1,616 of 5,938 gene-associated Recon3D reactions (206 of 3,697 enzymatic, 1,410 of 2,241 transport) | Use the complete Recon3D (10,600 reactions, both FDFT1 half-reactions present). No hypothesis-directed repair is needed (step 1) |
| I3. "Is sterol synthesis active in PMS?" was answered only indirectly (HMGCR inactive in a signed fit; HMGCR mRNA t = −0.9) | CarnivalFlow has no notion of flux capacity | Readouts R2–R4: synthesis flux range, de novo capacity with uptake blocked, simvastatin response |
| I4. Metabolic shortcuts and compound hubs in the joint CarnivalFlow fit | Reaction nodes in COSMOS are not constrained by any metabolic context | A context-specific COSMOS: remove reaction nodes whose Recon3D reaction carries no flux in either group's iMAT model, then refit the joint network (step 7) |

What this plan does not fix: the free-input non-identifiability of the signalling layer (one hub,
arbitrary identity) is a property of CarnivalFlow and stays a reported limitation.

## Facts established in a read-only check (2026-10-13)

- `data/resources/pkn/Recon3D.json` is a COBRA JSON model with 10,600 reactions, 5,835 metabolites
  and 2,248 genes. 5,938 reactions have a gene rule. Gene IDs are `<Entrez>_AT<n>`, with the symbol
  in `name` and the Entrez ID in `annotation.ncbigene`.
- All 1,560 exchange reactions (`EX_*`) are open for uptake by default (lower bound < 0), so a medium
  must be defined (step 2). Biomass reactions are `BIOMASS_reaction`, `BIOMASS_maintenance` and
  `BIOMASS_maintenance_noTrTr`. Cholesterol exchange is `EX_chsterol_e`. `HMGCOAR` (reaction 1600)
  has rule `3156_AT1` (HMGCR), and `r0170` and `HMR_1465` (FDFT1) are present.
- cobrapy is not installed, and is not needed. `corneto.io.graph_from_vertex_incidence(S, met_ids,
  rxn_ids)` builds the model from the stoichiometric matrix. Edge attributes `default_lb`,
  `default_ub` and `GPR` are what `MultiSampleFBA`/`MultiSampleIMAT` read, as in
  `corneto.io.cobra_model_to_graph`. GPR rules are evaluated by
  `corneto.methods.metabolism._utils.evaluate_gpr_rules` (AND = min, OR = max).
- `MultiSampleIMAT(lambda_reg, eps, high_expression_threshold, low_expression_threshold,
  use_bigm_constraints, ...)` exposes `build_many(model, gene_expression={cond: {gene: value}},
  reaction_bounds={cond: {rxn: (lb, ub)}}, objectives=...)`. `MultiSampleFBA` has the same
  interface for objectives and bounds.

## Unit of observation

As in stages 01–05, the unit is the cell line. Expression enters as group means (3 Ctrl, 3 PMS
lines). A per-line run (6 conditions) is a sensitivity only (Q5). With no measured fluxes, every
result is a capacity statement: "flux-consistent with the expression pattern", not a measurement.

## Method

1. **Model import** (`src/gem.py`).
   - Build the CORNETO graph from Recon3D.json: S from `reactions[*].metabolites`, bounds from
     `lower_bound`/`upper_bound`, and GPR from `gene_reaction_rule` with gene IDs rewritten to HGNC
     approved symbols. Use the stage-01 HGNC table via `annotation.ncbigene`. Genes without a unique
     approved symbol keep their Entrez ID and are listed.
   - Checks:
     - reaction and metabolite counts equal the JSON;
     - FBA maximising `BIOMASS_maintenance` is feasible and positive with all exchanges open;
     - with `HMGCOAR` and all cholesterol uptake blocked, maximal cholesterol demand is 0, and it is
       positive when either is allowed (a sanity test of the sterol module).
2. **Medium** (two registered variants).
   - (a) **Open**: Recon3D defaults (all uptakes allowed), as an upper-bound reference.
   - (b) **Defined** (primary): uptake allowed only for:
     - compounds detected in the ST003332 no-cell medium blanks (factor "no cell ctrl" in
       `data/Ionescu/MetabolomicsWorkbench/ST003332/`), mapped through KEGG → Recon3D;
     - the registered inorganic set: O2, H2O, H+, Pi, Na+, K+, Cl−, Ca2+, Fe2+, SO4, NH4, CO2;
     - lipoprotein-cholesterol uptake (`EX_chsterol_e`), in two sub-variants:
       - **chol+**: allowed;
       - **chol−**: blocked (Q2).

     Secretion stays open everywhere. If the defined medium cannot support `BIOMASS_maintenance`
     at 10 % of its open-medium maximum, report this. The smallest set of extra exchanges that
     restores it is found by MILP and listed, never added silently.
3. **Expression scores** (Q1).
   - For each group, the mean log2(CPM+1) of the 13,530 filtered genes.
   - Gene score: +1 for the top 25 % of filtered genes in that group, −1 for the bottom 25 % and
     for every Recon3D gene that is not in the filtered set (not expressed), 0 otherwise.
   - Reaction scores come from CORNETO's GPR evaluation. Reactions without a GPR score 0
     (`use_mean_for_missing_reactions=False`).
4. **iMAT context models.**
   - `MultiSampleIMAT(lambda_reg=0, eps=1e-2, use_bigm_constraints=True)` with `build_many` and
     conditions `{Ctrl, PMS}`, each with the same medium bounds. λ = 0 means the conditions are
     independent; structured sparsity across conditions is the sensitivity Q4.
   - HiGHS, 1 % gap, 3,600 s limit, 4 threads, seed 20261010. Accept a final gap ≤ 5 %.
   - One Python process per solve, with resume by task signature, reusing the
     `src/corneto_worker.py` pattern. A small `src/gem_worker.py` is needed because the task type
     differs.
   - Output per condition: reaction activity (flux above eps), flux values, agreement counts (high
     reactions on, low reactions off), and the objective.
5. **Readouts on the context models** (registered before running).
   - **R1. Differential activity:** reactions active in PMS only and in Ctrl only. Test enrichment
     by Recon3D subsystem with a hypergeometric test; the universe is gene-associated reactions;
     use BH FDR.
   - **R2. Sterol synthesis flux range:** within each context model (only the active reaction set
     allowed to carry flux), the minimum and maximum flux through `HMGCOAR`, through each canonical
     step (the 58 rows of `03_pkn/aim5_route/canonical_steps.tsv`) and through the last
     cholesterol-forming steps. Report also whether each of the 34 ordered steps (`step_order`) has
     at least one reaction carrying flux.
   - **R3. Source of cholesterol (metabolic task):** maximise a cholesterol demand reaction in each
     context model:
     - (i) with the medium as defined;
     - (ii) with cholesterol uptake blocked.

     If the maximum is positive in (ii), the cells can make cholesterol de novo; if it is 0 in
     (ii) but positive in (i), they must take it up. Compare Ctrl and PMS.
   - **R4. Simvastatin response:** repeat R3 with the `HMGCOAR` (and `HMGCOARx`) upper bound set
     to 10 % of its R2 maximum. The measured response, a 3.98 log2 fall in free cholesterol in all
     7 lines, predicts that cholesterol production capacity should drop. Report the direction and
     size for each group. This is the stoichiometric version of the stage-05 positive control.
6. **Comparison with the signed analysis.** For the canonical reactions, tabulate:
   - present in COSMOS;
   - used in the stage-04 joint fit;
   - active in the iMAT model (Ctrl, PMS);
   - flux range.

   For HMR_1916 (the "HDLBP" reaction), report whether it is active and in which direction. This
   shows whether the stage-04 shortcut is flux-consistent.
7. **Context-specific COSMOS (fixes I4).**
   - From the primary COSMOS PKN, remove every reaction node `Gene<N>__...` (and `_reverse`) whose
     Recon3D reaction N is inactive in both group models, or active only in the opposite direction.
   - Keep all signalling and metabolite–protein edges, then re-apply stage-03 pruning (drop nodes
     left with no edges).
   - Refit the stage-04 joint problem (cosmos_multi measurements, λ = 0.1, same settings) and the
     stage-05 metabolism-only simvastatin control on this network.
   - Readouts as in stages 04/05:
     - fit error;
     - TFs and metabolites fitted;
     - dominant starting node;
     - whether HMGCR or any canonical reaction is used;
     - how cholesterol is explained;
     - cholesterol → TF routes;
     - edge and node Jaccard with stage 04.
8. **Sensitivities** (only if time allows, in this order):
   - thresholds 10/90 % instead of 25/75 % (Q1);
   - medium variants (a) and chol+/chol− (Q2);
   - λ_reg > 0 across conditions (Q4);
   - per-line conditions (Q5).

## Open questions (recommended default first)

- **Q1. Expression thresholds:** group-specific quartiles (top 25 % = high, bottom 25 % and
  unexpressed = low). Alternatives: 10/90 %, or z-scores across genes.
- **Q2. Cholesterol uptake:** run both chol+ and chol−. The iNSC culture medium and its supplements
  are not stated in our data; check Ionescu et al. methods. B27 supplement contains cholesterol, N2
  does not. If the medium is identified, the matching variant becomes primary.
- **Q3. Context model definition for step 7:** remove reactions inactive in both groups. Stricter
  alternative: inactive in PMS (the condition being explained).
- **Q4. Structured sparsity across conditions:** off (λ_reg = 0) for the primary run.
- **Q5. Per-line iMAT:** sensitivity only. The six-sample CarnivalFlow problem was intractable, and
  per-line expression is noisy.
- **Q6. Isotope labelling** (ST003331 13C-glucose): not used as a flux constraint in this plan.
  Fractional labelling of mevalonate-pathway intermediates is not deposited, so it cannot test
  synthesis directly. Listed as future work.

## Files to change

- `src/gem.py` (new):
  - `load_recon3d(path, hgnc)` → CORNETO graph plus gene map table;
  - `medium_bounds(variant, ...)`;
  - `gene_scores(expression, quantiles)`;
  - `context_reactions(result, eps)`;
  - `prune_cosmos_to_context(pkn, active)`.
- `src/gem_worker.py` (new): one process per FBA/iMAT task, with `task.json`/`result.json`,
  signature and resume, as in `corneto_worker`.
- `scripts/analysis/08_metabolic_context.py` (new): steps 1–6 plus step-7 PKN writing.
- `scripts/analysis/06_corneto_robustness.py`: new part `context_cosmos` for the step-7 refits
  (`--cosmos`).
- `scripts/setup/run_stage06_gem.sh` (new).
- `tests/test_gem.py` (new):
  - toy 6-reaction model: FBA optimum;
  - iMAT selects the high-expression branch;
  - GPR AND/OR evaluation with symbols;
  - blocked uptake makes a demand infeasible;
  - context pruning removes only inactive `Gene<N>` nodes and their `_reverse` copies.
- `doc/decisions/README.md`: add this plan.

## Outputs

`results/ionescu_corneto/06_validation/metabolic_context/`:
- `model_checks.tsv`, `medium/<variant>.tsv`, `gene_map.tsv`;
- `imat/<medium>/<condition>/{task,result}.json` and `fluxes.tsv`;
- `summary/{differential_reactions,subsystem_enrichment,sterol_flux_range,cholesterol_tasks,simvastatin_task,signed_vs_flux}.tsv`;
- `context_cosmos/pkn.tsv` with a filter log;
- provenance sidecar.

Step-7 refits go to `05_robustness/context_cosmos/`.

## Validation before interpretation

- All `tests/test_gem.py` cases pass.
- Model checks in step 1 pass (biomass feasible; sterol sanity test).
- The defined medium supports biomass, or the added exchanges are listed.
- Every iMAT and FBA solve reports status, gap and time, and only accepted solves are interpreted.

## Runtime (estimate, course laptop)

| Step | Estimate |
|---|---|
| FBA/flux-range LPs | seconds to minutes each |
| iMAT on Recon3D (about 10,600 reactions, big-M MILP), 2 conditions × 2 media | roughly 10–60 min per solve (uncertain; the time limit stops it at 1 h) |
| Step-7 joint refit and simvastatin control | about 1–2 h |

Total for the primary analysis: about 3–6 h of solving. The sensitivities in step 8 double this.

## For the report

- **Question 2:** the result enters as one or two sentences: does the flux-consistent PMS model
  synthesise cholesterol, and does the shortcut disappear once the context is applied?
- **Figure 3:** possibly one extra row ("Active in iMAT, PMS").
- **Methods:** one paragraph.
- **Discussion:** "Signs are not chemistry" then describes a tested fix, not a suggestion.
