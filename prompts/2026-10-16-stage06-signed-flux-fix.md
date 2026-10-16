# Stage-06 signed flux reporting fix

Date: 2026-10-16
Model: openai/gpt-5.6-luna (OpenCode)
Settings: default
Plan: doc/decisions/2026-10-13-imat-fba-plan.md
Commit(s): add after committing
Transcript: not saved

## Prompt actually used

Commit the current worktree, but fix any clear Stage-06 reporting defect first.
The signed-versus-flux output must include HMR_1916 even though it is not in the
canonical cholesterol-step table, and must determine whether its actual Recon3D
reaction node was selected rather than testing whether any Gene node exists.

## What it produced

- `scripts/analysis/08_metabolic_context.py`: explicitly report HMR_1916 when
  present in the flux table and match its reaction number to selected COSMOS
  reaction nodes.

## What was changed by hand afterwards

- The `present_in_cosmos` field was set from the matched HMR reaction nodes
  rather than unconditionally set to `True`.

## Outputs

The current Stage-06 solver tree remains uncommitted: its provenance points to
the pre-commit dirty state, and the context-specific simvastatin refit failed.
