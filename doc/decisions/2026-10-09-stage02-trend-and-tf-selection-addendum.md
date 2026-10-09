# Stage 02 addendum: limma-trend contrast and the stage-04 TF-selection rule

Date: 2026-10-09
Status: proposed

Adds to `2026-10-09-stage01-02-revision-plan.md`. Origin: audit findings F1
and F3 (`doc/reviews/2026-10-09-stage01-02-implementation-audit.md`). Both
decisions have to be fixed **before** anyone looks at the contrast scores of
E2F1, EGR1, WT1, SP1, JUN, JUNB or AP1.

## Blinding record (fill in before confirming)

What has been exposed so far:
- `logs/stage01_02.log` and `02_activities/summary/diagnostics.json` print
  the row order of the exploratory heatmap (top 30 CollecTRI TFs by variance
  across the six lines). JUN, JUNB, SP1, SP3 and AP1 are in that list. The
  auditor saw the list, but no per-line or contrast values.
- `02_activities/figures/tf_activity_heatmap_collectri.*` shows per-line
  scores for those TFs.
- User: has the heatmap, `tf_activity_per_line.tsv` or any
  `tf_activity_contrast*.tsv` been opened? _To be filled by the user._

Per-line variance does not show the PMS-vs-Ctrl direction, so the contrast is
still blind if only the heatmap was seen. Record it either way.

## Decisions to confirm (recommended default first)

- **E1. Contrast statistic: limma-trend becomes the primary.** The pooled
  variance of log2(CPM + 1) falls about 4.5-fold from the lowest to the
  highest expression quintile (median 0.112 to 0.025; Spearman with mean
  expression -0.48). limma's guidance for log-CPM with library sizes within
  a factor of 3 (here 26.9-32.5 M, ratio 1.2) is `eBayes(trend = TRUE)`. The
  current constant-prior moderated t becomes a sensitivity. Alternative:
  keep the constant prior as primary and add trend as a sensitivity.
- **E2. Measured TFs for stage 04: fixed top-k by |score|, filtered for
  leave-one-line-out sign stability.** Universe: TFs of the primary resource
  (CollecTRI with complexes) kept after tmin = 5 that map to a node of the
  stage-04 PKN. Rank by |ULM score| on the primary contrast; take the top
  k = 50 (sensitivities k = 25 and 100). Keep a selected TF only if its sign
  is the same in all six leave-one-line-out contrasts; report how many were
  dropped. Adjusted p-values are reported but never used for selection.
  Alternatives: top-k without the stability filter; a |score| threshold.
- **E3. The six held-out TFs and AP1 get no special treatment.** They enter
  stage 04 only if E2 selects them; none is added or removed by hand.
- **E4. Complex sources (AP1, NFKB) are not stage-04 measurements in the
  primary run**, because they are not PKN nodes. Their activities stay in the
  stage-02 tables, and the D5 rule (AP1 counts as JUN/JUNB evidence in
  Aim 5) still applies to the activity layer. Whether a sensitivity run
  attaches AP1 to its member genes is a stage-04 decision, taken before
  E2's selection is looked at.

## Method

T6. limma-trend prior (limma `fitFDist` with a covariate, as used by
    `eBayes(trend = TRUE)`):
    - Covariate A_g = mean log2(CPM + 1) of gene g over the six lines.
    - e_g = log(s2_g) - digamma(df_g/2) + log(df_g/2), df_g = 4, for genes with
      s2_g > 0.
    - Fit e_g on a natural cubic spline basis in A_g with 4 degrees of freedom
      including the intercept: knots at the minimum, the 1/3 and 2/3
      quantiles and the maximum of A (same column space as R
      `splines::ns(A, df = 4, intercept = TRUE)`). Use the truncated-power
      basis N1 = 1, N2 = A, N_{k+2} = d_k - d_3 for k = 1, 2, with
      d_k(A) = [(A - xi_k)_+^3 - (A - xi_4)_+^3] / (xi_4 - xi_k) (Hastie,
      Tibshirani and Friedman, ESL eq. 5.4-5.5). Ordinary least squares in numpy.
    - Residual variance r = sum(residual^2) / (n_genes - 4);
      v = r - trigamma(df_g/2). If v > 0: solve trigamma(d0/2) = v (brentq, as
      in T2), s0_g^2 = exp(fitted_g + digamma(d0/2) - log(d0/2)). If v <= 0:
      d0 = inf, s0_g^2 = exp(fitted_g).
    - Posterior and moderated t as in T2, with gene-specific s0_g^2.
    - Record d0, the four spline coefficients and s0^2 at the A quartiles in
      `summary/counts.tsv`.
T7. Files: `<resource>/tf_activity_contrast.tsv` is now limma-trend;
    the T2 version moves to `<resource>/tf_activity_contrast_constant_prior.tsv`.
    `gene_contrast.tsv` gains `trend_prior_s0_squared`, `trend_posterior_variance`
    and `trend_moderated_t`.
T8. Leave-one-line-out contrasts (for E2): for each line, recompute the
    limma-trend moderated t on the remaining five lines (2 vs 3, df_g = 3,
    prior refitted on those five) and run ULM.
    `<resource>/loo_contrast/without_<line>.tsv`, 6 files per resource. The
    selection itself is applied in stage 04, after the PKN mapping.
T9. Provenance: the stage-02 code list hashes the revision plan, the HGNC
    plan and this addendum (audit F5 for stage 02).

The stage-02 cache will refuse to reuse `02_activities/` because the code
changes. Move it aside first, for example to
`results/ionescu_corneto/02_activities_archived_2026-10-09_constant_prior/`.
Stage 01 is unchanged and stays cached.

## Files to change

- `src/activities.py`: natural-spline basis, `moderated_t_statistics(..., trend=...)`,
  leave-one-line-out contrast helper.
- `scripts/analysis/03_tf_activities.py`: T7-T9.
- `tests/test_activities.py`: tests below.
- `doc/analysis-state.md`, `doc/decisions/README.md`, `prompts/` record.

## Validation

Synthetic tests:
- Constant prior (the revision plan's missing test (a)): sigma^2 =
  s0^2 d0 / chi2(d0), d0 = 6, s0^2 = 1, s2 = sigma^2 chi2(4) / 4, 10,000 genes,
  seed 20261006; estimated d0 within 25% of 6 and s0^2 within 15% of 1.
- Trend: the same simulation with log s0^2(A) = -1 - 0.5 A and A uniform on
  [0, 8]; the fitted s0^2 at A = 2, 4 and 6 is within 20% of the truth, and d0
  is within 25% of 6.
- No trend: on the constant-prior simulation, trend and constant d0 agree
  within 10%.
- The spline basis has rank 4, and its fitted values equal those of any
  other basis with the same knots (check against a B-spline basis built with
  `scipy.interpolate.BSpline` plus the natural boundary constraints, or
  against a hard-coded R `ns()` fixture).
- Leave-one-out contrast: the 5-line result equals a fresh run on those five
  lines and differs from any slice of the 6-line output.

Real run checks: d0 and s0^2 at the quartiles printed; Spearman between
trend and constant-prior TF scores reported (no TF names printed).
