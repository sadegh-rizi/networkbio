# SRA batch check, stage-02 addendum, stage-03 plan

Date: 2026-10-09
Model: claude-opus-5-5 (Claude Science)
Settings: default
Plans: doc/decisions/2026-10-09-stage02-trend-and-tf-selection-addendum.md,
doc/decisions/2026-10-09-stage03-pkn-plan.md (both proposed)
Commit(s): none yet (add after the user commits)
Transcript: not saved

## Prompt actually used

> done, can you move forward with
> Run the SRA metadata check.
> Write a short addendum for the two audit decisions.
> Draft the stage-03 plan in doc/decisions/ (status: proposed), with the Q6
> filter choices laid out for you to confirm.

## What it produced

- `doc/datasets.md`: section "GSE297192 baseline libraries C1-C6: sequencing
  batch". Sources: ENA filereport (PRJNA1262970), NCBI SRA efetch XML, and
  the first read headers of each `_1.fastq.gz` (HTTP range request on the ENA
  FTP). All six: instrument LH00409, run 296, flow cell 22GWVTLT4, lane 1,
  submission SRA2129539.
- The two plans named above.
- Not yet done: rows for both plans in `doc/decisions/README.md`, the
  `doc/analysis-state.md` next-step update and a row for this record in
  `prompts/README.md`. Those edits failed because the working tree was on
  `main` (baseline `7b1fa51`, before the `task/stage01-02` merge), whose
  copies of these files predate the stage-01/02 work. They are to be added
  after the merge.

Facts used in the plans and how they were checked: `omnipath` 1.0.12 API and
columns (read in `.venv`); NetworkCommons' OmniPath filter and COSMOS
download path (read from its GitHub source); Metabolomics Workbench KEGG ->
HMDB lookups for five IDs (live calls). Not checked: the COSMOS meta-PKN file
itself and its node grammar (host blocked in the sandbox; plan step R1), the
Recon3D JSON, the NetworkCommons data-server URL, and the runtime of the
rewiring null.

## What was changed by hand afterwards

_To be filled by the user._

## Outputs stamped

None.
