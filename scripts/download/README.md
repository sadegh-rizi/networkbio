# Download scripts: Ionescu et al. 2024 and Park et al. 2025

Accessions, sample counts and file sizes were checked on 2026-09-29 against
GEO, NCBI SRA, ENA, Metabolomics Workbench, MassIVE and Zenodo. The full
inventory (per dataset, per sample, per file) is in
`doc/data-inventory/data_inventory_ionescu_park.xlsx`.

These scripts were syntax-checked with shellcheck but **not executed**
(no working bash in the environment that wrote them). Try `01_processed.sh`
first; it is small and resumable.

## Layout

Run from the repository root in WSL. Files go under `./data` unless you set
`DATA_ROOT`. Quote paths: the repo path contains `#`.

```
data/
  Park/GEO/<GSE>/                 processed GEO files (+ per_sample/)
  Park/SRA/<GSE>/<GSM>_<title>/   raw FASTQ
  Ionescu/MetabolomicsWorkbench/<ST>/        processed JSON + mwTab
  Ionescu/MetabolomicsWorkbench/<ST>/raw/    raw LC-MS zip
  Ionescu/MassIVE/MSV000095283/   proteomics (when released)
  code/zenodo_<id>/               code archives
```

## Scripts

| Script | Content | Size |
|---|---|---|
| `01_processed.sh` | GEO metadata; bulk RNA count matrices; scRNA .h5; snATAC peak .h5; WGBS CpG tables; MW peak tables (JSON, mwTab, factors); Zenodo code | ~4 GB |
| `02_atac_fragments.sh` | snATAC `fragments.tsv.gz` for GSE297690 | ~6.2 GB |
| `03_mw_raw.sh` | raw LC-MS zips for ST003328/30/31/32 (`--with-duplicate` also gets ST003328_Rawfiles.zip) | ~8.7 GB |
| `04_sra_fastq.sh` | raw reads per GEO series, via ENA (default) or SRA Toolkit | 38 GB to >880 GB |
| `05_massive_proteomics.sh` | MassIVE MSV000095283. The dataset was private on 2026-09-29 | unknown |
| `06_reused_public.sh` | lists (or `--get`) GEO supplements of the datasets Park re-analysed | not checked |
| `verify.sh` | reports missing files and size mismatches against `manifest.tsv` | - |
| `download_all.sh` | 01 (+ 02, 03 with `--medium`), 05, verify | - |

`manifest.tsv` (tier, destination, URL, expected bytes) and `sra_runs.tsv`
(run, GSM, title, condition, sizes, ENA FASTQ URLs) drive the scripts.

## Commands (WSL, from the repo root)

```bash
cd "/mnt/d/#University/#MS-SystemsBio/Master_Thesis/repos/networkbio"
mkdir -p logs
nohup bash scripts/download/01_processed.sh > logs/download_01_processed.log 2>&1 &
bash scripts/download/verify.sh

# optional, larger
nohup bash scripts/download/02_atac_fragments.sh > logs/download_02_atac.log 2>&1 &
nohup bash scripts/download/03_mw_raw.sh > logs/download_03_mw_raw.log 2>&1 &

# raw reads: preview first, then choose series
bash scripts/download/04_sra_fastq.sh --dry-run GSE297192
nohup bash scripts/download/04_sra_fastq.sh GSE297192 > logs/download_04_bulk_fastq.log 2>&1 &
# 10x reprocessing with Cell Ranger needs index reads: use the SRA method
# (needs sra-tools: conda install -c bioconda sra-tools)
nohup bash scripts/download/04_sra_fastq.sh --method sra GSE297365 > logs/download_04_scrna_fastq.log 2>&1 &
```

Requirements: `wget` or `curl`, GNU coreutils (`stat -c`, `df --output`),
bash >= 4.4. `sra-tools` only for `--method sra`.

## Known gaps (details in `doc/data-inventory/data_inventory_ionescu_park.xlsx`, sheet `gaps`)

- Ionescu proteomics / secretome (MSV000095283): private on MassIVE.
- Ionescu bulk RNA-seq (2 Ctrl vs 2 PMS): no accession given in the paper, and none was found.
- Ionescu 13C tracing at 30 min and 6 h, and 13C-glutamine tracing: not in ST003331, which holds only the 24 h glucose time point.
- Park baseline sc/snRNA-seq of 3 Ctrl + 3 PMS lines (26,138 cells) and the multiome: not in GSE297365 or GSE297690, which hold only lines 105 and HD.
- WGBS runs SRR33556538-43 (6 samples added Aug 2026) reported size 0 in SRA and had no ENA FASTQ on 2026-09-29.
- Bulk RNA-seq C1-C24: GEO gives no donor line. C1-C3 = Ctrl and C4-C6 = PMS comes from the authors' code.
