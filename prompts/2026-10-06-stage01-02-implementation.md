# Stage-01 preprocessing and stage-02 TF activities

Date: 2026-10-06
Model: openai/gpt-6-astra (OpenCode)
Settings: generation settings not exposed; no project-analysis randomness.
Plan: doc/decisions/2026-10-06-stage01-02-preprocessing-tf-activities-plan.md
Commit(s): implementation uncommitted; baseline is 7b1fa51; branch task/stage01-02.
Transcript: not exported. This record summarises the OpenCode conversation
beginning with the implementation request below; the full tool trace is in
that OpenCode session. A stable exported transcript path is not available.

## Prompt actually used

> implement 2026-10-06-stage01-02-preprocessing-tf-activities-plan.md

The same request was sent twice. The agent found a proposed plan with six
unconfirmed choices and an untracked repository on main. Follow-ups:

> commit the current files and then I confirm working on a new branch

After the baseline commit and branch creation, the agent asked for explicit
confirmation of the six scientific choices (CPM filter, blank-referenced
medium, isotope pooling, replicate averaging, primary/secondary ULM, and
excluding ST003330). The user replied:

> confirm

## What it produced

- `src/preprocess.py`: approved-symbol mapping/count collapse, CPM filtering,
  within-analysis isotopologue pooling and total normalisation, blank
  references, line means with missingness counts, exploratory QC, shared
  provenance/hash-verified output-cache helpers.
- `src/activities.py`: immutable dated regulon caches, gene centring, Welch
  statistics with expression-scale effects and intervals, ULM and target
  counts, explicit raw/BH p-values, exploratory plots and diagnostics.
- `scripts/analysis/02_preprocess.py`, `03_tf_activities.py`: stage orchestration,
  schema/count/missingness checks, outputs and provenance.
- `scripts/download/07_resources.sh`: HGNC download/fallback and provenance.
- `scripts/setup/run_stage01_02.sh`: tests then dependent stage runs.
- `tests/test_preprocess.py`, `tests/test_activities.py`: synthetic numerical
  tests, local end-to-end fixtures, offline cache checks and mocked shell
  download/runner checks.
- Confirmed/implemented plan, analysis state, resource ledger/dataset notes,
  decision index, download instructions and prompt index updated.

## Evidence and implementation details

Read `doc/analysis-state.md`, repository rules, the target plan, dataset
documentation, existing stage-00 parser/script and plotting helpers. Probed
stage-00 TSV schemas with `pandas.read_csv(..., nrows=3)` using project Python;
read stage-00 provenance first. No complete real count matrix was analysed.
Inspected installed decoupler 2.2.0 sources (`_collectri.py`, `_dorothea.py`,
`_ulm.py`, `_run.py`, data extraction and network pruning) before implementing
the wrapper. This established that ULM returns BH-adjusted p-values and that
the default extractor drops empty genes; the implementation records these
behaviours and explicitly retains zero-centred genes. The upstream fixed
seed-0 row shuffle is recorded; synthetic random fixtures use seed 20261006.

Numerical design and the six approved choices are preserved. Additional
reporting files prevent mislabelling adjusted p-values as raw and preserve
the Welch input, uncertainty and requested diagnostics. No held-out TF/path
results were examined to tune settings.

## Verification

- Project environment versions checked with `importlib.metadata.version`:
  pandas 3.0.6, numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2,
  pytest 9.1.1, decoupler 2.2.0, omnipath 1.0.12, statsmodels 0.15.0.
- Initial preprocessing tests: 6 passed, 1 failed because pandas 3 supplies
  a read-only NumPy view to in-place PCA centring. Fixed by requesting a copy.
- First combined test attempt timed out at 120 seconds during first ULM
  startup. Retried with a 600-second limit: **12 passed in 302.45 seconds**,
  including synthetic stage-01/02 execution and real figure exports at
  reduced test DPI. These figures/tables are temporary fixtures, not project
  findings. The production figure DPI settings remain unchanged.
- Shell syntax checks and `git diff --check` passed at that checkpoint.
- Final verification (including added shell failure-path tests):
  `MPLBACKEND=Agg "/mnt/d/#University/#MS-SystemsBio/Master_Thesis/repos/networkbio/.venv/bin/python" -m pytest -q tests --basetemp=/tmp/opencode/stage01-02-tests`
  — **22 passed in 147.44 seconds**. This includes all five existing stage-00
  tests, the synthetic stage-01/02 run, mocked HTTP-404 HGNC fallback,
  immutable/offline cache checks, and runner failure propagation for each
  dependent step. Both new shell scripts pass `bash -n`; tracked changes
  pass `git diff --check`. Subsequent edits only finalised documentation and
  the exact DoRothEA URL recorded in metadata.

The agent did not download resources, run the real preprocessing/activity
pipeline, execute notebooks, install packages or modify stage-00 outputs.
Live resource availability, real mapping coverage, 14,490-gene QC and
biological replicate consistency remain checks for the user's WSL run.

## Outputs stamped

Real outputs are pending. Both stage scripts write a success provenance
sidecar last with inputs/resources, hashes, full settings, package/Python/
platform versions, git HEAD and dirty state, exact source-code hashes and
output hashes. A clean baseline commit is not misrepresented as the code
that generated later uncommitted changes.

## What was changed by hand afterwards

None recorded in this session.
