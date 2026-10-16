# Stage-06 GPR token, gene-score and sterol-check fix

Date: 2026-10-14
Model: claude-opus-5-5 (Claude, cloud session linked to the user's computer)
Settings: default
Plan: doc/decisions/2026-10-14-stage06-gpr-and-checks-addendum.md
Commit(s): (add after committing)
Transcript: not saved

## Prompt actually used

> [first Stage-06 run log: CORNETO "Unsupported AST node: BinOp" warnings, then
> OSError "Cannot save file into a non-existent directory: .../summary"]
> I am getting this output ?

Claude diagnosed the crash, the GPR parsing defects and the failed sterol check
from the run's own outputs and proposed four changes; the user replied "yes".

Second event (same session): the rerun stopped on the failed check
`cholesterol_demand_defined_chol_minus_hmgcr_blocked = 1000`; the user pasted
`model_checks.tsv`. Claude traced the flux to `SK_vitd3_c` and closed the
supply direction of `SK_`/`DM_` reactions in the defined media (addendum D4).

Third event (same session, 20:10): the user asked whether the run had
finished. It had stalled at 13:33 UTC; `medium_support.tsv` showed the MILP
bug D5, which Claude fixed (addendum item 9).

Fourth event (2026-10-15): the completed run showed D6-D8 (addendum items
10-12); the user confirmed all three fixes, including adding Phe and Ile.

## What it produced

- `src/gem.py`: `gpr_token()`; tokens used in GPRs; `gpr_token` column in
  `gene_map`; `gene_scores(..., unscored_genes=)` scores unmapped genes 0.
- `scripts/analysis/08_metabolic_context.py`: expression index tokenised;
  sterol checks in defined media; stop on failed checks; `summary/` created.
- `src/gem.py` (second event): `SUPPLY_PREFIXES`/`BOUNDARY_PREFIXES`;
  `medium_bounds` closes `SK_`/`DM_` supply; `minimum_extra_exchanges`
  candidates include them. `08_metabolic_context.py` medium tables list them.
- `src/gem.py` (fourth event): `UNMEASURED_ESSENTIAL_UPTAKES`,
  `mevalonate_synthesis_sign`, `cap_synthesis_direction`; `08_metabolic_context.py`:
  boundary bounds kept in context models, direction-aware simvastatin cap,
  synthesis-capacity column.
- `tests/test_gem.py`: six new tests in total.
- `doc/decisions/2026-10-14-stage06-gpr-and-checks-addendum.md`,
  `doc/decisions/README.md`, `doc/analysis-state.md`.

## What was changed by hand afterwards

None yet.

## Outputs stamped

None. Only an AST syntax check was run (the linked VM has Python 3.10 and
no pytest; the project venv is WSL Python 3.12). Tests not yet run.
