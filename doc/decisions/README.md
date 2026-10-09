# doc/decisions/

Confirmed analysis plans. These are the handoff contract between agents and
between sessions — see `doc/multi-agent-architecture.md`.

Convention: `YYYY-MM-DD-<short-topic>-plan.md`.

## Status field

Every plan carries an explicit status on the second line:

- `proposed`: written, not yet confirmed by the user. Do not implement.
- `confirmed`: user has said go. This is the version to implement.
- `implemented`: code exists; record the commit(s).
- `superseded by <file>`: do not delete, point forward.

## Template

```markdown
# <Short title>

Date: YYYY-MM-DD
Status: proposed | confirmed | implemented (commit <hash>) | superseded by <file>

## Question

What is actually being asked, and why now.

## Data available

Which datasets and files, and what is known versus unknown about them
(cross-reference `doc/datasets.md`). Say what is *not* available.

## Unit of observation

The biological unit, the experimental unit, and the effective n, in lines.
State the dependency structure: paired, nested, repeated.

## Method

Numbered steps, concrete enough to implement without further design
decisions. Include thresholds, seeds, covariates, solver settings, and the
multiple-testing family.

## Open questions

Anything still undecided, listed individually. An empty section on a
`confirmed` plan means every ambiguity was actually resolved.

## Files to change

Explicit list. Anything not on it is scope creep.

## Outputs

Exact paths under `results/<analysis>/<NN_stage>/`.

## Validation

The commands that will be run and what result would count as the step
working. State what would falsify the approach, not only what would confirm it.
```

## Current plans

For what has actually been run and decided, see `doc/analysis-state.md`.

| Plan | Status (as written in the file) | Note |
|---|---|---|
| `2026-09-29-corneto-ionescu-inscs-plan.md` | proposed | Q1, Q12 answered; defaults proposed 2026-10-06 |
| `2026-10-06-setup-toy-inputs-plan.md` | implemented (uncommitted) | venv, toy problem, stage-00 tables |
| `2026-10-06-stage01-02-preprocessing-tf-activities-plan.md` | implemented in part; revised by `2026-10-09-stage01-02-revision-plan.md` | parent implementation retained; revised plan changes primary medium/lipid normalisation and the contrast statistic |
| `2026-10-09-stage01-02-revision-plan.md` | implemented; real-data run completed 2026-10-09 | D1-D5 confirmed by the implementation request; audit in `doc/reviews/2026-10-09-stage01-02-implementation-audit.md` |
| `2026-10-09-hgnc-ambiguous-mapping-plan.md` | implemented; stages 01 and 02 completed 2026-10-09 | ambiguous approved Ensembl IDs are excluded and recorded rather than assigned arbitrarily |
| `2026-10-09-stage03-pkn-plan.md` | implemented (uncommitted) | P1-P10 implemented; real-data run completed 2026-10-09; COSMOS endpoint mismatch recorded |
