# Stage 01 preprocessing and stage 02 TF activities

Date: 2026-10-06
Status: proposed

Do not implement until the status says `confirmed`. Every default below was
proposed by the agent, not chosen by the user; the open questions list the
ones the user must confirm. Parent plan: `2026-09-29-corneto-ionescu-inscs-plan.md`
(steps 2-3). Previous stage: `2026-10-06-setup-toy-inputs-plan.md`.

## Question

Turn the stage-00 tables into (a) one normalised value per feature per line
(and per condition) for every layer, and (b) per-line transcription-factor
(TF) activities from the RNA-seq, so that stage 03 (PKN) and stage 04
(CORNETO) have measured-node values. Nothing here is inference; no
PMS-vs-control claim is made at this stage.

## Data available

All inputs are the stage-00 outputs in `results/ionescu_corneto/00_inputs/`
(read `00_inputs.provenance.json` first). Do not re-read `data/` except
through those tables.

| Layer | Stage-00 files | Lines | Replicates | IDs |
|---|---|---|---|---|
| RNA-seq GSE297192 | `rnaseq/counts_raw_baseline.tsv.gz`, `rnaseq/samples.tsv` | 3 Ctrl (C1-C3), 3 PMS (C4-C6) | 1 library per line | Ensembl gene ID, no version suffix; 78,932 rows |
| Intracellular metabolomics ST003331 | `metabolomics/ST003331/{samples,features,values}` | 3 Ctrl, 4 PMS | 3 | KEGG for all 158 unlabelled rows; 23 13C-isotopologue rows (9 KEGG IDs) |
| Extracellular metabolomics ST003332 | `metabolomics/ST003332/...` | 3 Ctrl, 4 PMS + 3 no-cell blanks | 3 | KEGG for all 136 rows |
| Lipidomics ST003328 | `metabolomics/ST003328/...` | 3 Ctrl, 4 PMS, each untreated and SV | 3 | names + lipid class only |

Facts checked by the planning agent on the stage-00 tables (2026-10-06):

- No KEGG ID appears in more than one analysis (ion mode) and no KEGG ID has
  more than one unlabelled row, in ST003331 or ST003332. Every isotopologue
  row has an unlabelled row with the same KEGG ID in the same analysis. Match
  isotopologues to their pool **by (analysis_id, kegg_id)**, never by name
  (names differ, e.g. "phosphomevalonate M-2H" vs "p-mevalonate 13C4 M-2H").
- Missing values (NaN; there are no zeros): ST003331 16 features (at most 6
  of 21 samples), ST003332 10 features (at most 8 of 24), ST003328 3 features
  (1 sample each).
- In ST003332, 87 of 136 metabolites have a median over cell-conditioned
  samples below 1.5x the median of the no-cell blanks. Conditioned medium
  is mostly fresh-medium components, and consumption shows up as values
  **below** the blank. So the blanks are a reference, not a filter.
- RNA-seq: 14,490 genes have CPM >= 1 in at least 3 of the 6 libraries
  (16,323 at CPM >= 0.5). Use this as the expected count for validation.

Not available: paired samples across layers; RNA-seq replicates within a
line; any SV arm outside the lipidomics; database IDs for lipids.

## Unit of observation

The cell line. Replicates are technical/culture replicates of a line and are
averaged to one value per line (per condition) at the end of stage 01.
Ctrl vs PMS lines are independent; untreated vs SV within a lipidomics line
is paired. Lines are matched **within one study only** (`line_key` =
`study:group:letter`); never join line keys across studies.

## Method

Use `pandas`, `numpy`, `decoupler` 2.2.0 and `omnipath` 1.0.12 from `.venv`
(already installed; do not add packages). Reusable functions go in
`src/preprocess.py` and `src/activities.py`; scripts only orchestrate.
All randomness: none expected; if any appears, seed 20261006.

### Resources (new download script, run by the user in WSL)

R1. `scripts/download/07_resources.sh` downloads into `data/resources/`
    (new folder; never edit existing `data/` files):
    - HGNC complete set, quarterly release 2026-10-01:
      `https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2026-10-01.txt`
      → `data/resources/hgnc/hgnc_complete_set_2026-10-01.txt`. If that file
      does not exist, use 2026-07-01 and record which one was used. (URL
      pattern taken from the pysec2pri source; not reachable from the
      planning agent's sandbox, so it is unverified.)
R2. CollecTRI and DoRothEA are fetched in Python by
    `scripts/analysis/03_tf_activities.py` the first time it runs
    (`decoupler.op.collectri(organism="human")`,
    `decoupler.op.dorothea(organism="human", levels=["A","B","C"])`) and cached
    as `data/resources/regulons/{collectri,dorothea_ABC}_human_<YYYY-MM-DD>.tsv.gz`.
    Later runs read the cache. Record download date, row count, sha256.

### Stage 01: preprocessing (`scripts/analysis/02_preprocess.py`)

RNA-seq
1. Gene symbols: map `ensembl_gene_id` to `symbol` with the HGNC table
   (columns `ensembl_gene_id`, `symbol`; status "Approved" only). Unmapped
   genes are dropped and counted. If several Ensembl IDs map to one symbol,
   sum their counts and log how many symbols this affects.
2. Expression filter: keep a gene if CPM >= 1 in >= 3 of the 6 libraries.
   CPM = count / library size x 1e6, library size from the unfiltered
   matrix.
3. Normalise: log2(CPM + 1) on the filtered genes, library sizes recomputed
   on the filtered matrix. No batch correction (one batch as far as known).
4. Output gene x line matrix (6 columns, C1-C6, with group in a sample sheet).

Intracellular metabolomics (ST003331)
5. Pool isotopologues: for each (analysis_id, kegg_id), sum the unlabelled
   row and all its 13C rows per sample = total pool. NaN + value = value; a
   pool is NaN only if all its rows are NaN. Also write the fractional
   labelling (labelled sum / pool) as a separate table; it is not used
   downstream in this plan.
6. Per-total normalisation within each analysis_id: divide each sample by the
   sum of its non-NaN pooled values in that analysis, multiply by the
   median of those sums across samples. Then log2.

Extracellular metabolomics (ST003332)
7. Per-total normalisation within each analysis_id over cell-conditioned
   samples and blanks together (same rule as 6), then log2.
8. Exchange value: log2 value minus the mean log2 value of the 3 blanks for
   that metabolite (negative = net uptake, positive = net release). Keep the
   blanks out of all line-level tables. Do not filter metabolites on the
   blanks.

Lipidomics (ST003328)
9. Per-total normalisation within each analysis_id across all 42 samples
   (untreated and SV together), then log2. Keep the lipid class column.

All metabolomics and lipidomics layers
10. No imputation. A replicate value that is NaN stays NaN.
11. Line level: mean of the log2 replicate values per line (per condition for
    ST003328), over non-NaN replicates. If fewer than 2 of 3 replicates are
    non-NaN, the line value is NaN and the case goes into the exclusion table.
12. Write replicate-level and line-level tables. Line-level column names:
    `<group>_<letter>` (plus `_SV` / `_untreated` for ST003328), and a sample
    sheet mapping each column to line_key, group, treatment, n_replicates.

QC (exploratory figures only, via `src/plot_config.save_figure`)
13. Per layer: PCA of the replicate-level log2 matrix (features with no NaN,
    each feature centred), points coloured by group with a second cue
    (marker shape), labelled by line letter; replicate-correlation heatmap
    (Pearson on log2). Captions state n in lines. No statistical tests.

### Stage 02: TF activities (`scripts/analysis/03_tf_activities.py`)

14. Input: the stage-01 log2(CPM+1) gene x line matrix (6 lines).
15. Primary: per-line activities. Centre each gene across the 6 lines
    (subtract its mean over lines), then run `decoupler.mt.ulm` with
    CollecTRI, `tmin=5`, on a DataFrame of lines x genes. Score = ULM t-value.
    Sign convention: positive = the TF's targets are, weight-signed, higher in
    that line than the 6-line mean.
16. Secondary: one contrast. Per gene, Welch t-statistic PMS (C4-C6) vs Ctrl
    (C1-C3) on log2(CPM+1); then ULM on that single vector with CollecTRI,
    `tmin=5`. Report it, but state that with 3 vs 3 lines it has little power.
17. Sensitivity (proposal Aim 1): repeat 15 and 16 with DoRothEA levels A-C.
    Same filenames, variant in the path (see Outputs).
18. Do not threshold, select or rank TFs for inference here; selection of
    measured nodes is a stage-04 decision with its own plan.
19. Do not look at the paper's six TFs, SREBF1/2 or HMGCR to tune any
    setting in this plan. They are the held-out test (main plan, Aim 5).
20. Exploratory figures: heatmap of per-line activities for the 30 TFs with
    the highest variance across lines (rows clustered, columns ordered by
    group); scatter of CollecTRI vs DoRothEA per-line scores for shared TFs.

### Provenance

21. Each stage writes `<stage>.provenance.json`: input files with sha256,
    resource files with download date and sha256, every parameter above,
    package versions (`importlib.metadata.version`), Python, platform, and
    `git rev-parse HEAD` if available (else "uncommitted").

## Open questions

The user must confirm or change each before the status becomes `confirmed`.

1. Gene filter: CPM >= 1 in >= 3 of 6 libraries (gives 14,490 genes before
   symbol mapping).
2. Medium (ST003332): per-total normalisation then log2 ratio to the blank
   mean, as in steps 7-8. The authors used per-total normalisation without
   a blank reference; confirm the blank-referenced version is wanted.
3. Isotopologues: total pool (step 5) rather than unlabelled-only.
4. Replicates (main plan Q9): average log2 values per line (step 11), with the
   2-of-3 rule for missing values.
5. TF activities: per-line ULM on gene-centred data as primary (step 15), the
   3-vs-3 contrast as secondary (step 16).
6. Whether ST003330 (extracellular lipidomics) should be added. This plan
   leaves it out, as the proposal does.

## Files to change

New:
- `scripts/download/07_resources.sh` (HGNC; add its row to `scripts/download/README.md`)
- `src/preprocess.py`, `src/activities.py`
- `scripts/analysis/02_preprocess.py` (stage 01), `scripts/analysis/03_tf_activities.py` (stage 02)
- `tests/test_preprocess.py`, `tests/test_activities.py`
- `scripts/setup/run_stage01_02.sh` (runs pytest, then the two scripts)
- `prompts/2026-MM-DD-stage01-02-implementation.md`

Edit:
- `doc/analysis-state.md`, `doc/external-tools.md` (HGNC, CollecTRI, DoRothEA
  rows with dates), `doc/datasets.md` (resources), `doc/decisions/README.md`
  (status), `prompts/README.md`.

Anything else is out of scope. Do not change `requirements*.txt`, `data/`
files that exist, stage-00 code or outputs.

## Outputs

```
results/ionescu_corneto/01_preprocessing/
  rnaseq/gene_map.tsv                 ensembl_gene_id, symbol, status, kept
  rnaseq/log2cpm.tsv                  symbol x line (C1-C6)
  rnaseq/samples.tsv
  metabolomics/ST003331/{features.tsv, replicate_log2.tsv, line_log2.tsv, line_samples.tsv, fractional_labelling.tsv}
  metabolomics/ST003332/{features.tsv, replicate_log2_vs_blank.tsv, line_log2_vs_blank.tsv, line_samples.tsv}
  metabolomics/ST003328/{features.tsv, replicate_log2.tsv, line_log2.tsv, line_samples.tsv}
  summary/exclusions.tsv              layer, feature or line, reason
  summary/counts.tsv                  features in/out at every step, per layer
  figures/pca_<layer>.{pdf,_presentation.png,_manuscript.png}
  figures/replicate_correlation_<layer>.{...}
  01_preprocessing.provenance.json
results/ionescu_corneto/02_activities/
  collectri/tf_activity_per_line.tsv  TF x line, ULM score
  collectri/tf_pvalue_per_line.tsv
  collectri/tf_activity_contrast.tsv  TF, score, pvalue (PMS vs Ctrl)
  collectri/regulon_summary.tsv       n TFs, n edges before/after tmin and gene filter
  dorothea_ABC/<same four filenames>
  figures/tf_activity_heatmap_collectri.{...}
  figures/collectri_vs_dorothea.{...}
  02_activities.provenance.json
```

## Validation

Commands (the user runs them in WSL; the implementing agent does not run
pipelines on the user's machine):

```bash
cd "/mnt/d/#University/#MS-SystemsBio/Master_Thesis/repos/networkbio"
bash scripts/download/07_resources.sh
bash scripts/setup/run_stage01_02.sh 2>&1 | tee logs/stage01_02.log
```

Tests that must pass (synthetic data, no network):
- Isotopologue pooling: two rows (unlabelled + 13C) sum correctly; a NaN row
  does not make the pool NaN; pooling keys on (analysis_id, kegg_id).
- Per-total normalisation: equal sample totals after normalisation within an
  analysis; NaN ignored in the total.
- Line averaging: the 2-of-3 rule gives NaN for 1 non-NaN replicate.
- ULM sign: a toy net where TF_A activates g1-g5 and represses g6-g10; a line
  with g1-g5 up and g6-g10 down gets a positive score, the opposite line a
  negative one.

Checks on the real run (print them and write them into `summary/counts.tsv`):
- 6 RNA-seq lines; 7 lines in ST003331 and ST003332; 14 line-conditions in
  ST003328.
- Genes after the expression filter = 14,490 before symbol mapping.
- Stop and report if fewer than 90% of filtered genes map to an approved
  symbol, if any ST003331 pool cannot be formed, or if any line has fewer
  than 2 replicates for more than 10% of its features.
- Every TF in `tf_activity_per_line.tsv` has at least 5 targets among the
  filtered genes (`regulon_summary.tsv`).

What would show the approach is wrong: per-line TF activities that track
library size or the order of C1-C6 rather than anything else (check the
correlation of each line's mean absolute score with its library size and
report it); PCA in which replicates of one line do not cluster together;
CollecTRI and DoRothEA per-line scores uncorrelated for shared TFs
(report the median Spearman correlation across lines).
