# Stage 01-02 revision plan (plan document, no code)

Date: 2026-10-09
Model: claude-opus-5-5 (Claude Science)
Settings: default
Plan: doc/decisions/2026-10-09-stage01-02-revision-plan.md (this record covers writing it)
Commit(s): none yet (add after the user commits)
Transcript: not saved

## Prompt actually used

> can you review this plan doc/decisions/2026-10-06-stage01-02-preprocessing-tf-activities-plan.md

followed by

> write the revised plan and I will implement it using another agent

## What it produced

- `doc/decisions/2026-10-09-stage01-02-revision-plan.md` (status: proposed;
  decisions D1-D5 to confirm).
- A row in `doc/decisions/README.md`.

Facts in the plan were read from `results/ionescu_corneto/00_inputs/`
(ST003331, ST003332 and ST003328 feature tables; RNA-seq library sizes), the
decoupler 2.2.0 source in `.venv` (`mt/_ulm.py`, `mt/_run.py`,
`op/_collectri.py`), and the Metabolomics Workbench sample-preparation text
recorded earlier in `doc/datasets.md`. Not verified: whether plain metabolite
rows are M+0 only (step R0 of the plan checks the deposit). The implementation
of the parent plan (`src/preprocess.py`, `src/activities.py`) was not audited
line by line.

## What was changed by hand afterwards

_To be filled by the user._

## Outputs stamped

None.
