# Datasets

Last verified: 2026-10-06 (Metabolomics Workbench schema); 2026-09-29 (the rest). Facts below were read from the papers' text
(Ionescu et al. from the supplied PDF `mmc7.pdf`, Park et al. from the open
access full text) unless marked **unverified**. Record the download date and
file list here whenever a dataset is added.

## Download log

- **2026-09-29, `scripts/download/01_processed.sh`** (run by the user in WSL;
  log `logs/download_01_processed.log`). 55 files, 3.85 GB, all present with the
  sizes in `scripts/download/manifest.tsv`. Every `.gz` decompresses, the MW
  `.json` files parse, the mwTab files have the MW header, the Zenodo zips pass
  a zip test and the 8 `.h5` files have an HDF5 signature (checked 2026-09-29).
  - `data/Park/GEO/<GSE>/`: SOFT + series matrix for GSE297192, GSE297365,
    GSE297690, GSE251839; `GSE297192_raw_expression.csv.gz`,
    `GSE297192_normalised_expression.csv.gz`, `GSE251839_metadata_WGBS_new_samples.xlsx`.
  - `data/Park/GEO/<GSE>/per_sample/`: 4 scRNA `.h5` (GSE297365), 4 snATAC peak
    `.h5` (GSE297690), 8 WGBS CpG tables + 6 Bismark `.cov` (GSE251839). The
    ATAC fragment files and their indexes come with `02_atac_fragments.sh`.
  - `data/Ionescu/MetabolomicsWorkbench/<ST>/`: `_data.json`, `_factors.json`,
    `_summary.json` and one mwTab per analysis for ST003328, ST003330,
    ST003331, ST003332.
  - `data/code/`: the three Zenodo code archives.
- Not yet downloaded: ATAC fragments (`02`), raw LC-MS (`03`), raw reads (`04`),
  MassIVE proteomics (`05`, private).

### Park single-cell and WGBS files: schema (inspected 2026-09-29)

- **scRNA-seq, GSE297365** (4 files): 10x HDF5, CellRanger 7.1.0, reference
  `refdata-cellranger-GRC38-3.1.0`, 60,623 features. These are the **raw
  (unfiltered)** matrices, with 1.76-1.95 M barcodes each, so cells have to be
  called first. Barcodes with >1,000 UMIs: 105 ABT-263 9,252; 105 DMSO 3,847;
  HD ABT-263 11,513; HD DMSO 9,531. Design: 1 PMS line x 1 Ctrl line x
  DMSO/ABT-263, no replicate libraries.
- **snATAC-seq, GSE297690** (4 files): 10x HDF5, cellranger-atac 2.1.0, GRCh38,
  raw peak x barcode matrices (483-509 k barcodes). **Each library has its own
  peak set** (278,809 / 276,049 / 243,265 / 241,277 peaks), so the matrices
  cannot be joined as they are; a shared peak set needs the fragment files
  (`02_atac_fragments.sh`). Barcodes with >1,000 fragments in peaks: 7,528-11,612.
- **WGBS, GSE251839**, two formats:
  - Batch 1 (8 files, `GSM79888xx_WGBS_*.txt.gz`: lines 103, 113, DR, 111,
    fibroblast + iNSC): `chr:start,percent` per 100-bp tile (14.5 M tiles in
    the 103 iNSC file, 24 chromosomes). **No read counts**, so coverage-aware
    tests (methylKit, DSS) cannot be run from these files.
  - Batch 2 (6 files, `*.bismark.cov.gz`: lines 105, HD, fibroblast / iNSC /
    iPS-NSC): Bismark coverage, `chr start end percent n_meth n_unmeth` per
    cytosine, strands separate (53.6 M rows in 105 iNSC, 453 contigs including
    alt/unplaced). Median depth per row 9; 49% of rows >= 10 reads (from every
    50th row).
  - The two batches have different formats and resolution. Comparing them needs
    per-CpG data for batch 1 too, i.e. re-alignment of the raw reads.

### Metabolomics Workbench files: schema (inspected 2026-10-06)

Parsed by `src/inputs.py`; tidy tables in `results/ionescu_corneto/00_inputs/metabolomics/`.

- Each study: `_data.json` (one record per analysis x metabolite: name, RefMet
  name, units "peak area", values per sample), `_factors.json` (sample,
  factors), and one mwTab per analysis (ion mode). The mwTab METABOLITES block
  holds `KEGG ID` (ST003331, ST003332) or `LipidGroup` (ST003328; no database
  IDs). Subject IDs are only in the mwTab SUBJECT_SAMPLE_FACTORS rows.
- **ST003328 lipidomics**: 642 species (217 negative, 425 positive mode); 42
  samples = 3 Ctrl ("AMC A-C") + 4 PMS lines ("PMA A", "PMA B", "PMS C",
  "PMS D"; the factors field says PMS for all four) x untreated/SV x P1-P3.
  Includes free cholesterol and 30 cholesteryl esters.
- **ST003331 intracellular metabolomics**: 181 rows (112 + 69), of which 23 are
  13C isotopologues ("Citrate 13C2" ...); 21 run IDs DN47-xx, subjects A-C
  Ctrl and D-G PMS, 3 each. Contains mevalonate and phosphomevalonate.
- **ST003332 extracellular metabolomics**: 136 rows (79 + 57); 24 samples =
  Ctrl "AMC A, B, D" + PMS "A-D" x S7-S9, plus 3 no-cell medium blanks
  ("NC 7-9"). Subject IDs A-G run across groups.
- Line letters do not identify the same line across studies (for example
  ST003328 PMS "A" is subject A, ST003332 PMS "A" is subject D). Treat line
  identity as within-study only until the paper's Table S1 says otherwise.
- The deposited documentation describes "technical triplicates per line".

The deposits themselves (sample counts, groups, files, sizes, gaps) were
checked against the repositories on 2026-09-29. The per-sample and per-file
tables are in `doc/data-inventory/` (`data_inventory_ionescu_park.xlsx`;
design slides in `experimental_design_ionescu_park.pptx`). Download scripts
are in `scripts/download/`. Where the inventory and the "unverified" notes
below disagree, the inventory is the newer source.

## Cohort (shared model)

Directly reprogrammed iNSCs from skin fibroblasts of people with PMS and
healthy controls. Ionescu: 3 controls and 4 PMS (3 secondary progressive,
1 primary progressive), aged 25 to 63, both sexes (5 male, 2 female), with
fibroblasts from the New York Stem Cell Foundation repository and from
the University of Wuerzburg (two control lines). Park describes 3 healthy
controls aged 25 to 63 and cites Ionescu for how the lines were made and
quality-controlled.

**Donor identity across the two papers is unverified.** Control lines are
named A and B in Ionescu and C1 and C2 in Park. Compare the two Table S1
files before stating that the lines are the same.

Limitations stated by the authors: no isogenic controls, differing age and
genetic background, small donor number.

## Ionescu et al. 2024, Cell Stem Cell 31:1574-1590

Article: DOI 10.1016/j.stem.2024.09.014. The file `mmc7.pdf` in the
repository root is the article plus supplementary figures S1-S5. The
supplementary tables (S1-S5) are separate files, not in this PDF.

| Layer | Design (from the paper) | Deposit |
|---|---|---|
| Intracellular proteomics | Ctrl vs PMS | MassIVE MSV000095283 |
| Secretomics (conditioned medium proteomics) | Ctrl vs PMS, with and without simvastatin (SV) | MassIVE MSV000095283 |
| Intracellular metabolomics | steady state; 13C-glucose (30 min, 6 h, 24 h) and 13C-glutamine (24 h) tracing | Metabolomics Workbench ST003331 |
| Extracellular metabolomics | steady state | Metabolomics Workbench ST003332 |
| Intracellular lipidomics | Ctrl vs PMS, with and without SV | Metabolomics Workbench ST003328 |
| Extracellular lipidomics | Ctrl vs PMS | Metabolomics Workbench ST003330 |
| Bulk RNA-seq | 2 Ctrl vs 2 PMS lines, at least 2 replicates each | not listed in the paper's deposited data |
| Functional assays | Seahorse, lipid droplets, senescence, neuron toxicity | not deposited |

Code: Zenodo 10.5281/zenodo.13764808 (MOFA and covariation network).

Design notes:

- Most experiments: n = 3 Ctrl and n = 4 PMS lines, each in at least 3
  replicates. The authors averaged replicates per line before MOFA and
  GENIE3.
- The authors' MOFA used the intracellular lipidome (642 features) and the
  secretome (241 features) from untreated and SV-treated lines, 3 factors.
- Proteomics: minimum-value imputation then median normalisation; only
  proteins seen in more than 20% of lines were kept. Metabolomics and
  lipidomics: per-total normalisation and standardisation.
- The secretome has many missing values from detection limits.
- Cell-level proteomics searched against the human UniProt database;
  conditioned-medium proteomics against a mouse database (as written in the
  Methods). Check this before joining identifiers. **Unverified** which
  identifiers the deposit uses.

**Unverified:** which deposited files exist and in what format; whether the
deposits include SV arms for metabolomics and proteomics.

## Park et al. 2025, Neuron

DOI 10.1016/j.neuron.2025.09.022 (PMC12834227).

| Layer | Accession |
|---|---|
| Bulk RNA-seq | GEO GSE297192 |
| Whole-genome bisulfite sequencing (fibroblasts and iNSCs) | GEO GSE251839 |
| scRNA-seq | GEO GSE297365 |
| snATAC-seq (some multiome) | GEO GSE297690 |

Code: Zenodo 10.5281/zenodo.15581507 and 10.5281/zenodo.15584699. An
interactive app is at `mohorianulab.org/shiny/pluchino/DARG_PMS/`.

**Unverified:** number of lines and replicates per accession, whether the
matrices are raw counts or processed, and the genome annotation used.

## Traps to check before using any file

- Row/column orientation and identifier type (gene symbol, UniProt,
  Ensembl, metabolite name, HMDB or KEGG ID).
- Whether values are raw, normalised or log-transformed, and whether
  missing values are NA or zero.
- Metabolite and lipid annotation level (identified vs putative).
- Sample naming: map every column to a line and condition explicitly.
