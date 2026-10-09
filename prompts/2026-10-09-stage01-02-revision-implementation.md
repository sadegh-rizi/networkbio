# Stage-01/02 revision implementation

Date: 2026-10-09
Model: openai/gpt-5.6-luna (OpenCode)
Settings: default
Plan: doc/decisions/2026-10-09-stage01-02-revision-plan.md
Commit(s): uncommitted
Transcript: not saved

## Prompt actually used

> Implement the revised plan here doc/decisions/2026-10-09-stage01-02-revision-plan.md

The request explicitly confirmed D1-D5 after the plan's proposed status. The
implementation agent first performed R0 against the deposited mwTab files.

## What it produced

- `src/labelling.py`: fixed carbon-count table, deposited/missing isotope forms,
  whitelist classes and coverage reporting.
- `src/preprocess.py`: filtered-count handoff, size-factor sensitivity,
  positive peak-area log transform and variant figure stems.
- `src/activities.py`: moderated t, line-subset centring, separate CollecTRI
  complex resource settings, and group-aware diagnostics.
- `scripts/analysis/02_preprocess.py`: R0 record, labelling classes, renamed
  fractional labelling, revised primary normalisations and all variants.
- `scripts/analysis/03_tf_activities.py`: moderated/logFC/Welch contrasts,
  leave-one-line-out activities, size-factor diagnostics and complex variant.
- `tests/test_labelling.py`, `tests/test_preprocess.py`,
  `tests/test_activities.py`: synthetic revision checks and end-to-end fixture.
- Decision, analysis-state, dataset, external-tool and prompt indexes updated.

## R0 finding

The ST003331 and ST003332 mwTab `ST:`, `SP:`, `TR:`, `AN:` and `MS:` sections
document 24-hour [U-13C6]glucose exposure and extraction, but no explicit
statement was found that plain metabolite rows are total pools. The result is
recorded at runtime in `summary/labelling_convention.md`.

## Verification

- Python compile and shell syntax checks passed.
- Labelling/preprocessing tests passed: 15 passed.
- Activity unit tests passed: 7 passed, 1 deselected for the end-to-end test.
- Synthetic stage-01/02 end-to-end test passed, including figure generation,
  provenance/cache reuse, all variants, three resources, contrasts and LOO
  outputs.
- Complete synthetic suite passed: 29 passed in 333.93 seconds.
- Stage-00 feature-table annotation smoke check passed: the nine traced
  metabolites and their deposited/missing forms matched the revised table;
  no transformed stage output was generated.
- The real data pipeline and resource downloads were not run.

## What was changed by hand afterwards

None recorded in this session.

## Outputs stamped

Real outputs are pending. The user must run the WSL resource download and
stage runner after reviewing the synthetic results.

## Update after the real run (2026-10-09)

The statements above that the real pipeline was not run were true when this
record was written. Afterwards the HGNC ambiguity fix was added
(`doc/decisions/2026-10-09-hgnc-ambiguous-mapping-plan.md`) and the user ran
stages 01 and 02 in WSL (`logs/stage01_02.log`, 30 tests passed, stage 02
finished 12:22). Outputs: `results/ionescu_corneto/01_preprocessing/` and
`02_activities/`, each with its provenance sidecar. Audit:
`doc/reviews/2026-10-09-stage01-02-implementation-audit.md`.
