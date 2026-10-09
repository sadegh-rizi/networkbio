# Stage 01-02 revision: isotope labelling, normalisation variants, contrast statistic

Date: 2026-10-09
Status: implemented; real-data run completed 2026-10-09 (`logs/stage01_02.log`); audited in `doc/reviews/2026-10-09-stage01-02-implementation-audit.md`

Revises `2026-10-06-stage01-02-preprocessing-tf-activities-plan.md` (the
"parent" below). Everything in the parent stays in force unless this file
changes it. Where the two disagree, this file wins. The parent's confirmed
choice 2 (ST003332 per-total normalisation as primary) is changed by D2 and
needs re-confirmation.

Origin: review of the parent plan on 2026-10-09, checked against the stage-00
tables, the decoupler 2.2.0 source in `.venv`, and the Metabolomics Workbench
sample-preparation text. No real stage-01/02 outputs exist yet
(`results/ionescu_corneto/` holds only `00_inputs/`), so nothing has to be
archived.

## Question

Same as the parent: produce line-level values for each layer and per-line TF
activities. This revision fixes five problems found in review:

1. ST003331/ST003332 come from a 24 h [U-13C6]glucose experiment. Plain
   metabolite rows are very likely the unlabelled (M+0) isotopologue only,
   so for glucose-derived metabolites they measure the unlabelled fraction,
   not the abundance. The parent treats them as abundances.
2. Isotopologue sums are incomplete (several isotopologues not deposited),
   but the parent calls them "total pool".
3. Per-total normalisation removes global differences in samples that were
   already extracted at fixed cell number (cells) or fixed volume (medium),
   and in the medium it mixes blanks into the totals.
4. The Welch t contrast with 3 vs 3 lines (about 2-4 df) is dominated by genes
   with small sample variance.
5. CollecTRI in decoupler 2.2.0 keeps the AP1 and NFKB complexes as sources
   by default (`remove_complexes=False`), which affects how JUN/JUNB (two of
   the held-out TFs) are read in Aim 5. This must be decided before any
   result is seen.

## Data available

As in the parent. Facts added by this review (checked 2026-10-09 on
`results/ionescu_corneto/00_inputs/`):

- ST003331 has 23 isotopologue rows for 9 KEGG IDs. Deposited labelled
  forms, and the forms that are missing given the carbon count:

  | Metabolite | KEGG | C | Deposited 13C forms | Missing 13C forms |
  |---|---|---|---|---|
  | 2-Oxoglutarate | C00026 | 5 | 2, 3, 4, 5 | 1 |
  | cis-Aconitate | C00417 | 6 | 2, 3, 4 | 1, 5, 6 |
  | Citrate | C00158 | 6 | 2, 3, 4, 5 | 1, 6 |
  | Fumarate | C00122 | 4 | 2, 3, 4 | 1 |
  | Glucose | C00031 | 6 | 6 | 1, 2, 3, 4, 5 |
  | Itaconate | C00490 | 5 | 2 | 1, 3, 4, 5 |
  | Malate | C00149 | 4 | 2, 3, 4 | 1 |
  | Phosphomevalonate | C01107 | 6 | 4 | 1, 2, 3, 5, 6 |
  | Succinate | C00042 | 4 | 2, 3, 4 | 1 |

- Every other ST003331 metabolite has only a plain row. These include
  pyruvate, 3-phosphoglycerate, fructose 1,6-bisphosphate, glucose
  6-phosphate, alanine, aspartate, glutamate, serine, UDP-glucose, ribose
  1-phosphate, nucleotides, acetyl-carnitine, fatty acids and (R)-mevalonate.
  Lactate is not in the cell data.
- ST003332 (medium) has no isotopologue rows at all, although the medium
  contained 20.4 mM [U-13C6]glucose. Plain rows include D-Glucose, Lactate,
  Pyruvate, Citrate, L-alanine, L-glutamate and L-glutamine.
- Sample preparation (Metabolomics Workbench `SP:SAMPLEPREP_SUMMARY`): cells
  extracted at 2e6 cells/mL (ST003331, ST003328); medium 20 uL into fixed
  solvent volume (ST003332).
- Lipidomics ST003328: 217 + 425 rows; no metabolite name occurs in both ion
  modes, so no cross-mode deduplication is needed.
- RNA-seq counts per library: Ctrl 26.9-29.9 M, PMS 27.1-32.5 M
  (`00_inputs/rnaseq/library_sizes.tsv`). Library size is partly aligned with
  group.

Not known: whether plain rows are M+0 only. The naming ("D-Glucose" next to
"Glucose 13C6", "Citrate" next to "Citrate 13C2-5") strongly suggests it.
This plan treats plain rows as M+0 unless step R0 finds otherwise.

## Unit of observation

Unchanged: the cell line; replicates averaged on the log2 scale per line
(per condition for ST003328); lines matched within one study only.

## Decisions to confirm (recommended default first)

- **D1. Plain rows are M+0.** Treat plain rows in ST003331/ST003332 as the
  unlabelled isotopologue; flag every metabolite that glucose carbon can
  reach within 24 h as `m0_confounded` (step M2). Alternative: treat them as
  abundances, as the parent does (not recommended).
- **D2. Medium normalisation.** Primary: no per-total normalisation (fixed
  volume); exchange value = log2 peak area minus the mean log2 of the three
  blanks. Sensitivity: the parent's per-total (blanks included) version.
- **D3. Cell and lipid normalisation.** ST003331 primary: per-total (as the
  authors did), sensitivity: none. ST003328 primary: none (fixed cell
  number, so a global simvastatin effect is kept), sensitivity: per-total
  across all 42 samples (the parent's version).
- **D4. Contrast statistic.** Primary: limma-style moderated t (step T2).
  Sensitivities: mean log2 difference (logFC), and the parent's Welch t.
- **D5. CollecTRI complexes.** Keep `remove_complexes=False` (decoupler
  default) as primary, and register now, before any result: in Aim 5 an
  "AP1" activity counts as evidence for JUN/JUNB, and "NFKB" is reported but
  is not one of the six TFs. Sensitivity: `remove_complexes=True`.
  DoRothEA A-C has no complexes and is unchanged.

RNA-seq normalisation is not a decision: log2(CPM + 1) stays primary and a
size-factor version is added as a sensitivity (step N1).

## Method

Only the changes and additions to the parent are listed. Parent step
numbers are cited as P<n>.

### R0. Check the isotopologue convention (before coding)

R0. Search the deposited mwTab files
    (`data/Ionescu/MetabolomicsWorkbench/ST003331/*.mwtab.txt`, and ST003332)
    for any statement on how unlabelled rows were reported (sections `MS:`,
    `AN:`, `ST:`, comments). Record the finding, or "no statement found", in
    `01_preprocessing/summary/labelling_convention.md`. If the deposit says
    plain rows are total pools, stop and report; D1 must then be revisited.

### Metabolomics labelling (replaces P5 for naming; adds new tables)

M1. Isotopologue sum (P5, renamed). Keep the P5 arithmetic (sum the plain
    row and its 13C rows per sample by (analysis_id, kegg_id); NaN + value =
    value; NaN only if all rows are NaN). Call the result
    `isotopologue_sum_deposited`, never "total pool". For each of the 9
    traced metabolites write the deposited and missing 13C forms, using the
    carbon counts in the table above (hard-coded in `src/labelling.py` and
    checked by a test against the stage-00 feature table).
M2. Labelling class. Every ST003331 and ST003332 feature gets one class:
    - `isotopologue_sum_partial`: the 9 traced metabolites (ST003331 only).
    - `m0_unconfounded`: KEGG ID in the whitelist below. Glucose carbon cannot
      reach these in human cells within 24 h, so the plain row is close to
      the abundance.
    - `m0_confounded`: everything else, including any metabolite whose
      synthesis can take up glucose carbon. This is deliberately conservative.

    Whitelist (`src/labelling.py`, one reason string per entry):
    essential amino acids L-histidine C00135, L-leucine/isoleucine C00123,
    L-lysine C00047, L-methionine C00073, L-phenylalanine C00079, L-threonine
    C00188, L-tryptophan C00078, L-valine C00183; L-tyrosine C00082 (from
    Phe); L-methionine S-oxide C02989; L-selenomethionine C05335; L-carnitine
    C00318 (from Lys and Met); ascorbate C00072 and dehydroascorbate C05422
    (not synthesised in humans); pantothenol C00864 and D-4'-phosphopantothenate
    C03492 (vitamin B5); linoleic acid C01595 (essential fatty acid); and
    tryptophan or tyrosine catabolites without added carbon: N-formylkynurenine
    C02700, anthranilate C00108, picolinic acid C10164, quinolinic acid
    C03722, 2-aminomuconate C00322, indole C00463, indoxyl C05658,
    indolepyruvate C00331, indole-3-acetaldehyde C00637,
    5-hydroxyindoleacetate C05635, 3-methyleneoxindole C02796,
    L-adrenaline C00788.
    The implementing agent must not add to this list. Proposed additions go
    into the implementation record for the user to decide. Write which
    whitelist IDs were found in each study.
M3. Medium glucose. In ST003332, D-Glucose (C00031) gets class
    `m0_confounded` with the note "13C6-glucose medium; uptake not measured
    by the M+0 row". Lactate, pyruvate, alanine and citrate fall into
    `m0_confounded` by M2.
M4. Fractional labelling (P5). Rename the output to
    `fractional_labelling_uncorrected.tsv` and state in its header comment
    that there is no natural-abundance correction and that missing
    isotopologues bias it downwards.
M5. Interpretation (replaces the P8 sentence). The sign rule "negative = net
    uptake, positive = net release" applies only to `m0_unconfounded`
    features. For `m0_confounded` features the exchange value is the change
    in the unlabelled fraction. Say this in the output header and in
    `summary/counts.tsv`.
M6. All line-level and replicate-level tables keep every feature. The class
    goes into each layer's `features.tsv` (new column `labelling_class`, plus
    `deposited_13c` and `missing_13c` for the traced nine). Measured-node
    selection stays a stage-04 decision, but the stage-04 plan must not use
    `m0_confounded` features as abundances.

### Normalisation variants (changes P6, P7, P9; adds RNA-seq sensitivity)

N1. RNA-seq sensitivity: DESeq2 median-of-ratios size factors on the
    filtered symbol-level counts (Anders and Huber 2010). Geometric means over
    genes with all six counts > 0; size factor = median over those genes of
    count / geometric mean. Value = log2(count / size factor + 1). Write
    `rnaseq/variants/sizefactor/log2norm.tsv` and the six size factors. The
    primary `rnaseq/log2cpm.tsv` is unchanged.
N2. ST003331: primary = P6 (per-total within analysis_id, then log2);
    variant `none` = log2 of the pooled values.
N3. ST003332: primary = no per-total; log2 peak area, then exchange value =
    log2 minus the mean of the available blank log2 values (parent rules
    for missing blanks unchanged). Variant `pertotal_authors` = P7 + P8
    exactly as in the parent.
N4. ST003328: primary = log2 peak area, no per-total; variant
    `pertotal_all42` = P9 exactly as in the parent. Untreated vs SV is
    compared only within a line (paired); this stage computes no contrast.
N5. Layout: primary tables keep the parent's paths. Each variant goes in
    `<layer>/variants/<variant>/` with the same filenames. QC PCA (P13) is
    drawn for primary and variants.

### TF activities (changes P15-P17)

T1. Per-line activities (P15) unchanged, with one addition: the function
    takes the set of lines to use and centres on those lines only. Write
    leave-one-line-out activities, `<resource>/loo/without_<line>.tsv`
    (6 files per resource), each recomputed with centring on the remaining
    five lines, never sliced from the 6-line result.
T2. Contrast, primary = moderated t (Smyth 2004, limma `eBayes`). In plain
    numpy/scipy, no new package:
    - Per gene on the stage-01 log2(CPM + 1): group means, mean difference
      d = PMS - Ctrl, pooled variance s2 with df_g = 6 - 2 = 4.
    - Prior: z = log(s2), e = z - digamma(df_g/2) + log(df_g/2). Let
      e_bar = mean(e) and v = var(e) - trigamma(df_g/2). If v > 0, solve
      trigamma(d0/2) = v for d0 (`scipy.optimize.brentq`, trigamma =
      `scipy.special.polygamma(1, x)`), and set
      s0^2 = exp(e_bar + digamma(d0/2) - log(d0/2)). If v <= 0, set d0 = inf
      and s0^2 = exp(e_bar) (limma `fitFDist`); the posterior variance is
      then s0^2 for every gene.
    - Posterior variance s_post^2 = (d0 s0^2 + df_g s2) / (d0 + df_g).
      Moderated t = d / sqrt(s_post^2 (1/3 + 1/3)), df = d0 + df_g.
    - Exclude genes with s2 = 0 from the prior fit (log undefined) but keep
      them in the output with the moderated t computed from s_post^2.
    - Record d0, s0^2, and the number of genes used in `summary/counts.tsv`.
    Run ULM (CollecTRI, tmin=5) on the moderated-t vector:
    `<resource>/tf_activity_contrast.tsv` (primary, same filename as the
    parent).
T3. Contrast sensitivities: ULM on the mean log2 difference
    (`tf_activity_contrast_logfc.tsv`) and on the parent's Welch t
    (`tf_activity_contrast_welch.tsv`). `gene_contrast.tsv` gains columns
    for the pooled s2, s_post^2 and moderated t.
T4. CollecTRI complexes: primary `remove_complexes=False`; variant
    `collectri_nocomplex/` with `remove_complexes=True`, same five files. In
    `regulon_summary.tsv` mark the complex sources (`is_complex`). Counting
    targets for any TF is allowed; looking at any of the six TFs' activities
    to choose between variants is not (parent step 19).
T5. DoRothEA A-C (P17): unchanged, plus the T1 leave-one-out and the T2/T3
    contrast files.

### Diagnostics (changes P "what would show the approach is wrong")

G1. The library-size check stays, but report group means of library size
    and state in `diagnostics.json` that library size is partly aligned
    with group, so at n = 6 the correlation cannot separate an artefact from
    biology. Add the same correlation for the size-factor variant.
G2. Report, for each metabolomics layer, the share of features per labelling
    class, and the PCA for primary vs variant normalisation side by side.

## Open questions

None for D1-D5. Stage 03 still decides mapping KEGG/HMDB IDs to the PKN's
metabolite identifiers and which `labelling_class` values may become measured
nodes.

## Files to change

New:
- `src/labelling.py` (carbon-count table, whitelist with reasons, class
  assignment)
- `tests/test_labelling.py`
- `prompts/2026-10-09-stage01-02-revision-implementation.md`

Edit:
- `src/preprocess.py`, `src/activities.py`
- `scripts/analysis/02_preprocess.py`, `scripts/analysis/03_tf_activities.py`
- `tests/test_preprocess.py`, `tests/test_activities.py`
- `doc/decisions/2026-10-06-stage01-02-preprocessing-tf-activities-plan.md`:
  change only its status line to
  `Status: implemented in part; revised by 2026-10-09-stage01-02-revision-plan.md`
- `doc/decisions/README.md` (rows for both plans), `doc/analysis-state.md`,
  `prompts/README.md`

Unchanged: `requirements*.txt`, existing `data/` files, stage-00 code and
outputs, `scripts/download/07_resources.sh`, `scripts/setup/run_stage01_02.sh`
(unless a new test file must be added to its pytest call).

## Outputs

Parent paths stay. Additions:

```
results/ionescu_corneto/01_preprocessing/
  summary/labelling_convention.md
  metabolomics/ST003331/features.tsv          + labelling_class, deposited_13c, missing_13c
  metabolomics/ST003331/fractional_labelling_uncorrected.tsv   (renamed)
  metabolomics/ST003331/variants/none/{replicate_log2,line_log2,line_samples}.tsv
  metabolomics/ST003332/features.tsv          + labelling_class
  metabolomics/ST003332/variants/pertotal_authors/{replicate_log2_vs_blank,line_log2_vs_blank,line_samples}.tsv
  metabolomics/ST003328/variants/pertotal_all42/{replicate_log2,line_log2,line_samples}.tsv
  rnaseq/variants/sizefactor/{log2norm.tsv,size_factors.tsv}
results/ionescu_corneto/02_activities/
  collectri/tf_activity_contrast.tsv          now moderated t (primary)
  collectri/tf_activity_contrast_logfc.tsv
  collectri/tf_activity_contrast_welch.tsv
  collectri/loo/without_<line>.tsv            6 files
  collectri_nocomplex/<same files as collectri/>
  dorothea_ABC/<same files as collectri/>
```

Note the changed meaning of the primary files: ST003332 and ST003328
primary tables are no longer per-total normalised, and
`tf_activity_contrast.tsv` is no longer the Welch version.

## Validation

Tests (synthetic, no network), in addition to the parent's:
- Labelling: the carbon table reproduces the "Missing 13C forms" column
  above from the stage-00 ST003331 feature table; every ST003331/ST003332
  row gets exactly one class; whitelist IDs absent from both studies are
  listed, not silently ignored.
- Isotopologue sum: unchanged parent tests pass under the new name.
- Normalisation: variant `none` equals log2 of the input; ST003332 primary
  exchange value equals log2(x) - mean log2(blanks); the `pertotal_authors`
  variant reproduces the parent's output on the same fixture.
- Size factors: on a fixture where library B = 2 x library A for every gene,
  the size-factor ratio is 2 and the normalised values are equal.
- Moderated t: (a) simulate 10,000 genes with true variance
  sigma^2 = s0^2 d0 / chi2(d0) (d0 = 6, s0^2 = 1) and observed
  s2 = sigma^2 chi2(4) / 4 (seed 20261006); the estimated d0 is within 25%
  of 6 and s0^2 within 15% of 1. (b) When all genes have the same s2, v <= 0,
  d0 = inf, and the moderated t equals d / sqrt(s0^2 (2/3)) with
  s0^2 = exp(log(s2) - digamma(2) + log(2)). (c) For a gene with s2 far below
  s0^2 the moderated |t| is smaller than the ordinary pooled |t|.
- Leave-one-out: activities for 5 lines equal a fresh run on those 5 lines
  and differ from slicing the 6-line output on a fixture with a non-zero
  left-out line.
- Complexes: on a toy CollecTRI table containing "AP1", the primary keeps it
  and the variant drops it.

Real run (user, in WSL):

```bash
cd "/mnt/d/#University/#MS-SystemsBio/Master_Thesis/repos/networkbio"
bash scripts/setup/run_stage01_02.sh 2>&1 | tee logs/stage01_02.log
```

Checks printed and written to `summary/counts.tsv`: all parent checks;
23 isotopologue rows in 9 traced metabolites; class counts per layer;
d0 and s0^2 of the moderated t; number of AP1/NFKB targets.

What would show the revision is wrong or unnecessary:
- R0 finds that plain rows are total pools (then D1 and M2-M5 are
  withdrawn).
- Primary and variant normalisations give essentially the same line-level
  PCA and Ctrl-vs-PMS ordering for a layer (then the variant was not needed
  for that layer; report it, do not delete it).
- The moderated-t and logFC contrasts give TF rankings almost identical to
  Welch (then point 4 of the Question did not matter for these data).
