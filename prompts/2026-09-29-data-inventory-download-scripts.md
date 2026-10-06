# Data inventory and download scripts

Date: 2026-09-29
Model: claude-opus-5-5 (Claude Science)
Settings: default
Plan: none (direct instruction; no analysis code)
Commit(s): none yet (add after the user commits)
Transcript: not saved

## Prompt actually used

> in D:\#University\#MS-SystemsBio\Master_Thesis\repos\networkbio
> I have the papers in the data folder, Let's say I later want to download the associated data with these papers. Can you give me the experimental design (maybe in the format of a word or powerpoint) and then also give me the list of all the data with their associated GEO numbers (or identifiers) and also the number of samples or metadata for that type of data, and then give me the bash scripts for downloading the data when I want to

Follow-up: "put all the files you generated in the repo".

## What it produced

- `scripts/download/`: `lib.sh`, `01_processed.sh`, `02_atac_fragments.sh`,
  `03_mw_raw.sh`, `04_sra_fastq.sh`, `05_massive_proteomics.sh`,
  `06_reused_public.sh`, `verify.sh`, `download_all.sh`, `manifest.tsv`,
  `sra_runs.tsv`, `README.md`.
- `doc/data-inventory/`: `data_inventory_ionescu_park.xlsx`,
  `samples_ionescu_park.csv`, `files_ionescu_park.csv`,
  `experimental_design_ionescu_park.pptx`, `README.md`.
- This record, and a row in `prompts/README.md`.

Sources queried on 2026-09-29: GEO SOFT and FTP file lists (GSE297192,
GSE297365, GSE297690, GSE251839, GSE251841), NCBI SRA runinfo, ENA
filereport, Metabolomics Workbench REST and mwTab (ST003328, ST003330,
ST003331, ST003332), the MassIVE dataset page (MSV000095283, private),
Zenodo API (13764808, 15581508, 15584700), and the authors' DARG_PMS code
(bulk RNA-seq C1-C3 Ctrl / C4-C6 PMS). Paper text from the PDFs in `data/`.

## Verification

- All shell scripts pass `shellcheck -s bash` (SC1091, SC2155, SC2034 excluded).
- Every URL in `manifest.tsv` answered a HEAD request with 200 and the
  expected size, except the Metabolomics Workbench REST URLs, which reject
  HEAD (406) but answered GET.
- **The scripts have not been executed**: bash does not run in the
  environment that wrote them. Test with `bash scripts/download/01_processed.sh`.
- The slide deck was not rendered in PowerPoint; check the layout.

## What was changed by hand afterwards

_To be filled by the user._

## Outputs stamped

None.
