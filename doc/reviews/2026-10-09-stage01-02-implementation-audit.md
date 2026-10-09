# Audit: stage 01-02 implementation and real run

Date: 2026-10-09
Auditor: claude-opus-5-5 (Claude Science); implementation by another agent
Scope: `doc/decisions/2026-10-09-stage01-02-revision-plan.md` (and its parent,
and `2026-10-09-hgnc-ambiguous-mapping-plan.md`) against `src/labelling.py`,
`src/preprocess.py`, `src/activities.py`, `scripts/analysis/02_preprocess.py`,
`scripts/analysis/03_tf_activities.py`, the tests, `logs/stage01_02.log` and the
outputs in `results/ionescu_corneto/01_preprocessing/` and `02_activities/`.
Blinding: the activities of the six held-out TFs (E2F1, EGR1, WT1, SP1, JUN,
JUNB) and of AP1 were not looked at. All checks below are about method and
data, not about which TFs come out.

## Verdict

The implementation matches the plan and the method sources. Two independent
recomputations reproduce the outputs to floating-point precision. Nothing has
to be rerun. There are two method-level points to decide before stage 03/04
(F1, F2), one interpretation rule to fix before stage 04 (F3), and some small
gaps in tests, provenance and documentation (F4-F8).

## What was checked and how

| Step | Source checked against | Result |
|---|---|---|
| R0 isotopologue convention | mwTab text (implementation record) | No total-pool statement found; recorded in `summary/labelling_convention.md` |
| M1 pooling | plan M1 | 181 ST003331 rows - 23 isotopologue rows = 158 features; 9 traced pools |
| M2/M3 classes | plan whitelist | Whitelist copied exactly, one reason per ID; 29/29 IDs found in ST003331, 21 in ST003332 (absent IDs listed); D-Glucose note set |
| M4 fractional labelling | plan M4 | Renamed; header comment states no natural-abundance correction |
| N1 size factors | Anders & Huber 2010 | Geometric means over genes with all counts > 0, median ratio; 13,244 genes used; factors 0.94-1.13 |
| N2-N4 normalisation | plan D2/D3 | Primary/variant assignment correct for all three studies; pooling precedes per-total, so traced metabolites count once in the totals |
| N3 medium | independent recomputation | `replicate_log2_vs_blank` = log2(x) - mean log2(blanks), max abs difference 3.6e-15 |
| HGNC mapping | HGNC plan | 3 Ensembl IDs with >1 approved symbol excluded and listed; 14,490 genes pass CPM before mapping, 13,530 symbols after (93.4% mapped) |
| T2 moderated t | limma `fitFDist`/`eBayes` (Smyth 2004) | Formula matches line by line (bias-corrected log variances, var with n-1, trigamma inversion, s0^2, v <= 0 branch, posterior variance). Independent recomputation: d0 = 2.853, s0^2 = 0.0434, max abs difference in t 7e-15 |
| ULM | decoupler 2.2.0 `mt/_ulm.py`, `mt/_run.py` | Score = slope t (`tval=True`) over all supplied genes (adjacency built over all features); raw p reconstructed with df = n_genes - 2 is exact; BH per observation |
| T1 leave-one-out | plan T1 | Recomputed with centring on 5 lines (differs from slicing the 6-line output: max abs difference 1.41, r = 0.993) |
| T4 complexes | decoupler `op.collectri` | Primary keeps AP1 (321 targets) and NFKB (466); variant has 1,183 TFs = 1,185 - 2 |
| T5 DoRothEA | decoupler `op.dorothea` | Weights = mode of regulation / confidence (A 1, B 2, C 3), the decoupler default, recorded in settings |
| Cache | `reuse_stage` | Stage 01 written 11:27 after the last code edit (11:05); reuse at 12:21 requires identical inputs, parameters, code hashes and versions |

## Findings

**F1. The moderated-t prior ignores the mean-variance trend (method; the plan's
gap, not the implementer's).** The pooled variance of log2(CPM + 1) falls with
expression: median s2 by expression quintile 0.112, 0.072, 0.041, 0.030, 0.025
(Spearman with mean expression -0.48). One prior s0^2 = 0.043 therefore shrinks
low-expressed genes towards too small a variance (inflating their |t|) and
high-expressed genes towards too large a one. limma's answer for log-CPM is
`eBayes(trend = TRUE)` (limma-trend): s0^2 is a smooth function of the gene's
average log expression. Recommendation: add `tf_activity_contrast_trend.tsv`
as a sensitivity (fit e against mean log2 expression with a lowess, then the
same d0 estimate on the residual variance), and decide now, before looking at
any held-out TF, whether it replaces the primary. The decision rests on a data
property, not on TF results, so changing the primary is allowed.

**F2. RNA-seq PC1 is 69% of the variance and splits Ctrl from PMS perfectly.**
`figures/pca_GSE297192.*`. Library size is also partly aligned with group
(group means 28.1 M vs 29.5 M). At n = 3 vs 3, disease cannot be separated from
a processing batch that coincides with group. Check the SRA run metadata for
the six libraries (instrument, run date, flow cell, `LoadDate`) and record the
result in `doc/datasets.md`. If library prep or sequencing runs are aligned
with group, every TF contrast inherits that and the stage-04 interpretation
has to say so.

**F3. ULM p-values on a contrast vector are not line-level evidence.** Number
of CollecTRI TFs with padj < 0.05: moderated t 68, logFC 134, Welch 32. The
ULM p-value treats the ~13,500 genes as independent observations and tests
regulon enrichment in one gene-level vector; it says nothing about
replication across lines, and it changes threefold with the contrast
statistic. Stage 04 must not select measured TFs by a padj threshold unless
the rule is written into the stage-04 plan first (for example a fixed top-k by
|score|, or a leave-one-line-out stability criterion).

**F4. Tests the plan asked for that are missing or weaker.**
- Moderated t (a), the simulation that recovers d0 and s0^2, is missing. Only
  the v <= 0 case and the shrinkage direction are tested. Risk is low because
  the formula was verified against limma and recomputed above.
- Leave-one-out: no unit test compares a 5-line run with slicing. Verified on
  the real outputs instead (above).
- The carbon-table test uses a synthetic fixture; the check against the real
  stage-00 table was a one-off smoke check outside the suite.

**F5. Provenance hashes the wrong plans.** Stage 01 hashes only the parent plan
(`PLAN` in `02_preprocess.py`), not the revision or HGNC plans, so editing
those would not invalidate the stage-01 cache. Stage 02 hashes only the
revision plan. Low risk; fix when next touching the scripts.

**F6. M5 header not written.** The uptake/release rule (valid only for
`m0_unconfounded`) appears in `summary/counts.tsv` and
`labelling_convention.md`, but not in the header of the ST003332 line table.
Acceptable if the stage-04 plan cites it; otherwise add the header.

**F7. G2 not side by side.** Primary and variant PCAs are separate files
rather than one panel pair. Cosmetic.

**F8. Status lines are stale.** The revision plan and `doc/analysis-state.md`
still say "real-data run pending"; the HGNC plan says "stage 02 pending"; the
implementation record says the real pipeline was not run. All were overtaken
by the 12:22 run.

## Did the revision matter?

The plan listed what would show it was unnecessary:
- Contrast statistic: CollecTRI TF scores, Spearman moderated vs Welch 0.963,
  moderated vs logFC 0.910, logFC vs Welch 0.813; top-50 overlap 41/50
  (moderated vs Welch) and 39/50 (vs logFC). The choice changes the edges of
  the ranking, not its core.
- Size factors vs CPM: size factors 0.94-1.13; mean |activity| still
  correlates with library size (r = 0.57 CPM, 0.68 size factor), which is
  uninterpretable here because library size is partly aligned with group.

## Not checked

- Stage-01 PCA and correlation figures other than `pca_GSE297192` and
  `pca_ST003332`; these two look as expected (replicates cluster by line,
  blanks separate on PC1 of the medium data).
- The parent plan's steps that this revision did not change, beyond what the
  code read-through covered.
- The download script `07_resources.sh`, except that its outputs carry
  checksums that the stage scripts verify.
