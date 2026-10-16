# Stage 06 addendum: GPR tokens, unmapped-gene scores and sterol checks

Date: 2026-10-14
Status: confirmed by user (2026-10-14, in chat); implemented (uncommitted)
Amends: `2026-10-13-imat-fba-plan.md`, steps 1 and 3

## Why

The first real-data run of `scripts/analysis/08_metabolic_context.py`
(2026-10-14, output in `results/ionescu_corneto/06_validation/metabolic_context/`)
stopped at `summary/reaction_activity.tsv` (missing directory), and inspection
of the run showed three defects. Its iMAT results must not be interpreted.

| Defect | Evidence | Effect |
|---|---|---|
| D1. Hyphenated symbols break the GPR parser | CORNETO log `Unsupported AST node: BinOp` for complex I/IV rules; `gene_map.tsv` has 12 such symbols (MT-ND1..6, MT-ND4L, MT-CO1..3, MT-CYB, NT5C1B-RDH14) | CORNETO parses GPRs as Python; `MT-ND3` reads as `MT - ND3`. 5 reactions (`NADH2_u10mi`, `CYOR_u10mi`, `CYOOm2i`, `CYOOm3i`, `r0047`) get the default score 0 instead of their (high) expression score |
| D2. Unmapped Entrez IDs are scored as their number | 20 genes with `mapping_status = no_approved_symbol` keep a numeric token (e.g. `2713`); `corneto/methods/metabolism/_utils.py` returns `float(node.value)` for a numeric `ast.Constant` | 75 reactions (52 fatty-acid oxidation) get scores such as 2713, i.e. "highly expressed" |
| D4. Defined medium leaves sinks open (found on the rerun, 2026-10-14) | Rerun `model_checks.tsv`: `cholesterol_demand_defined_chol_minus_hmgcr_blocked = 1000`, `passed = False`. Flux trace: `SK_vitd3_c` (−1000) → `VD3` → `PVD3` (reverse) → 7-dehydrocholesterol → `DHCR72r` → cholesterol. Recon3D has 92 `SK_` and 3 `DM_` reactions with lower bound < 0; 54 supplied metabolites in this check | `medium_bounds` closed only `EX_` uptakes, so the "defined" medium was not defined. All Recon3D `EX_`/`SK_`/`DM_` reactions are `met ->` (coefficient −1), so a negative flux is uptake |
| D3. Sterol sanity check tests the wrong thing | `model_checks.tsv`: `cholesterol_demand_hmgcr_and_uptake_blocked = 1000`, `passed = False`; the run continued | The check used the open medium, where cholesterol sulfate (`EX_chsterols_e`), cholesterol esters (`CHOLESTle_1`) and 7-dehydrocholesterol (`EX_7dhchsterol_e`) still supply cholesterol. The defined media open no other sterol uptake |

## Changes

1. `src/gem.py`: `gpr_token()` replaces non-identifier characters with `_`
   and prefixes identifiers that start with a digit with `ENTREZ_`.
   `load_recon3d` writes these tokens into the GPRs and adds a `gpr_token`
   column to `gene_map`.
2. `scripts/analysis/08_metabolic_context.py`: the expression index is mapped
   with the same `gpr_token` before scoring; a collision raises.
3. **Change to plan step 3 (gene scores):** model genes without a unique
   approved symbol score **0** (unknown), not −1. Their expression cannot be
   looked up, so "not expressed" would be an assumption. Genes with a symbol
   that are absent from the filtered expression set still score −1, as planned.
4. Sterol checks run in the defined media:
   - `defined_chol_minus` + HMGCOAR/HMGCOARx blocked: expected 0;
   - `defined_chol_minus`: expected > 0 (de novo synthesis possible);
   - `defined_chol_plus` + HMGCOAR/HMGCOARx blocked: expected > 0 (uptake).
5. `main()` raises if any model check fails, before any iMAT result is used.
6. `summary/` is created before the summary tables are written.
8. **Change to plan step 2 (medium), D4:** the defined media also close the
   supply direction (lower bound → 0) of every `SK_` and `DM_` reaction.
   Secretion stays open. `minimum_extra_exchanges` may now reopen `SK_`/`DM_`
   as well as `EX_` reactions if biomass needs them, and lists them. The
   medium tables list `SK_`/`DM_` rows too. The open medium is unchanged.
7. `tests/test_gem.py`: tests for `gpr_token`, for scoring of a hyphenated
   and an unmapped gene through CORNETO's `evaluate_gpr`, and for closing
   sink/demand supply in the defined medium.

9. **D5 (found on the third run, 2026-10-14):** with sinks closed, neither
   defined medium supports `BIOMASS_maintenance` (`medium_support.tsv`:
   0 vs 755 open), and `minimum_extra_exchanges` reported `infeasible`. Cause:
   the MILP kept the medium lower bound (0) on every candidate, so opening one
   could not change anything. Fixed: candidate lower bounds are relaxed to the
   model default and switched by the binary. Test added. Whether the selected
   boundary reactions are then added to the defined media (requires new iMAT
   solves) is an open decision for the user.

10. **D6, simvastatin cap on the wrong direction (found 2026-10-15 in the
    completed third run):** Recon3D writes `HMGCOAR` as mevalonate -> HMG-CoA,
    so HMG-CoA reductase activity is negative flux (`HMGCOARx` is written
    forward). The cap limited only the upper bound, so synthesis ran at -857
    and `simvastatin_task.tsv` equalled `cholesterol_tasks.tsv` (125 Ctrl,
    142.9 PMS). Fixed: the direction is read from the mevalonate coefficient
    (`gem.mevalonate_synthesis_sign`) and the cap applies to that direction
    (`gem.cap_synthesis_direction`). `sterol_flux_range.tsv` gains
    `mevalonate_synthesis_capacity`; the HMGCOAR `maximum_flux` is the reverse
    direction and is not synthesis. Bounds used are written to the task tables.
11. **D7, sink leak in the readouts:** `_active_context_bounds` kept medium
    bounds only for `EX_` reactions; active `SK_`/`DM_` reactions were reset to
    Recon3D defaults, reopening their supply in readouts R2-R4. Fixed: all
    boundary reactions keep the medium bounds.
12. **Change to plan step 2, D8 (confirmed by user 2026-10-15):** the defined
    media also allow uptake of L-phenylalanine and L-isoleucine
    (`gem.UNMEASURED_ESSENTIAL_UPTAKES`). Neither is in the ST003332 feature
    list, so they could not be detected in the blanks; all other essential
    amino acids are. Without them `BIOMASS_maintenance` is 0, and the repair
    MILP selected the tripeptide exchange `EX_lyspheile_e` (Lys-Phe-Ile), a
    one-reaction artefact of minimising the count. The iMAT models must be
    re-solved.

## Validation

- `tests/test_gem.py` passes (WSL project venv).
- `model_checks.tsv` all `passed = True`; otherwise the script stops.
- No `Unsupported AST node` warning in the log.

## Rerun

Archive the first-run directory (all of it: its iMAT scores and sterol checks
are invalid) and run Stage 06 again without `--resume`.
