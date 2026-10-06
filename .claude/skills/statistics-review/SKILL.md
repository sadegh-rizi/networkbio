---
name: statistics-review
description: Use before proposing, implementing, interpreting, or reviewing statistical tests, differential abundance, activity estimation, permutation nulls, model selection, or enrichment in the networkbio project.
---

# Statistics review

Read `doc/agent-rules.md` first. This skill is the checklist to run through before a test or a model is written and again before a result is interpreted.

## Before running

1. State the question, outcome, contrast, covariates, unit of analysis, and test in one sentence each.
2. The unit is the cell line. Say the effective n per layer (at most 3 control and 4 PMS). Replicates are averaged or modelled, never counted as n.
3. Say whether observations are independent (control vs PMS lines) or paired (simvastatin vs vehicle within a line), and use the matching test.
4. Say which factors are confounded with disease (age, sex, genetic background) and that they cannot be adjusted for with 7 lines.
5. Define the multiple-testing family and the correction (FDR for high-dimensional features) before testing. Report raw and adjusted values and the universe tested.

## Small-sample points that apply here

- With 3 vs 4 lines there are 35 distinct label permutations, so the smallest possible permutation p-value is 1/35 (about 0.029) for a statistic taken in absolute value. An exact two-sided rank-sum (Wilcoxon / Mann-Whitney) test with 3 vs 4 lines and no ties cannot go below 2/35 (about 0.057), and with 2 vs 2 lines not below 1/3. Do not present an FDR-adjusted table from such tests as if it had the resolution of a large study.
- Moderated tests (for example empirical-Bayes variance shrinkage) help with few replicates but do not create independent samples.
- Report effect sizes with intervals. At this n, an interval is more informative than a p-value.
- Mass-spectrometry data: missing values follow detection limits and are informative. Say how they were handled (filtering by detection frequency, imputation) and check that conclusions are not driven by imputed values.
- The authors' own filters (features seen in more than 20% of lines) are a documented choice, not a default to copy blindly.

## Activity estimation

- Activities (TF, pathway) are model outputs. Record the resource and version, statistic, minimum targets, and identifier mapping.
- A contrast built from 2 vs 2 lines cannot support a p-value claim. Prefer per-line activities.

## Networks and nulls

- An inferred network is not tested by its own fit. Validate against controls fixed in advance, use randomised inputs or randomised prior networks as a null, and check stability across seeds and parameters.
- A permutation of line labels is limited to 35 arrangements here; combine it with other nulls.

## Before interpreting

- Is the result exploratory or confirmatory? At this n, say exploratory.
- Would a different reasonable normalisation, imputation, or parameter give a different answer? Run the sensitivity check or state that it was not run.
- Do not claim causality, mechanism, or "drivers" from an association or an inferred network.
