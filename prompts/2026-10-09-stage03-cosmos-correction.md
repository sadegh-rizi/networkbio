# Stage-03 COSMOS correction

Date: 2026-10-09
Model: gpt-5.6-luna
Settings: default
Plan: doc/decisions/2026-10-09-stage03-cosmos-correction-plan.md
Commit(s): c7f3719
Transcript: not saved

## Prompt actually used

> Continue if you have next steps, or stop and ask for clarification if you are unsure how to proceed.

## Confirmed user decisions

- Keep COSMOS pseudo-reaction nodes without recognized gene tokens.
- Retain multi-gene enzyme nodes when any recognized gene is expressed.
- Keep expression pruning unchanged, leaving WT1 absent from the primary PKN.
- Replace the single eight-edge Aim-5 query with metabolic cap 80 and two
  signalling caps of 8.

## What this generation changes

- Corrects NetworkCommons `Gene<reaction-index>__<suffix>` parsing and rejects
  the incompatible X-prefixed grammar.
- Preserves numeric model metabolite identifiers instead of converting them to
  false HMDB IDs.
- Adds real-resource grammar fixtures and complete stage-input provenance.
- Regenerates the COSMOS PKN, mapping tables, three-leg Aim-5 paths/nulls and
  figures.
