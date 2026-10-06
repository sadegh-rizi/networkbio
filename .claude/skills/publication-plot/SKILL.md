---
name: publication-plot
description: Use when creating, revising, or reviewing exploratory or publication-oriented plots and tables for the networkbio project.
---

# Publication plot

Read `AGENTS.md` and `doc/agent-rules.md` first.

## Plotting workflow

1. Establish the analysis purpose: QC, exploration, model diagnostic, effect estimate, or final communication. Do not use a final-looking plot to conceal an unvalidated exploratory result.
2. Check the plotted unit, group definitions, sample count, transformation, and missing-data handling. State n as cell lines (donors), not only replicates or cells.
3. Use the Matplotlib object-oriented API for reusable figure code. Prefer a colorblind-safe palette and encode group differences with labels/markers as well as color when practical.
4. Show distributions and observations where sample size permits. Do not default to bars of mean ± SEM.
5. Label axes with units and transformations, name groups unambiguously, make legends readable, and avoid deceptive truncated axes.
6. Save or propose a clear output path outside `data/`; use PDF/SVG for final vector figures. Do not overwrite existing figures without permission.
7. Save every figure via `plot_config.save_figure`: a vector PDF plus a presentation-resolution PNG and a manuscript-resolution PNG (`_presentation.png` / `_manuscript.png` suffixes). Exact DPI/sizing per venue are placeholders (`PRESENTATION_DPI`, `MANUSCRIPT_DPI` in `src/plot_config.py`) meant to be tuned later, not fixed now — see `doc/agent-rules.md`'s "Figures and tables" section.

## Project conventions

- Shared style, palettes, and category ordering live in `src/plot_config.py`. Use them rather than hard-coding colors per figure. Cell-type and phenotype-state colors must be consistent across every figure in the project — a reader comparing two panels should not have to re-read two legends.
- Figures belong in `results/<analysis>/<stage>/figures/`, never in `data/`.
- Rasterize the point layer of large scatter plots (`rasterized=True`) and keep text/axes vector. A 300,000-cell UMAP as pure vector produces an unopenable PDF.

## Completion criteria

State the data source, transformation, grouping, number of observations (at the correct level), output path, and whether the figure was actually generated. Do not claim statistical significance from visual separation alone.

## Pitfalls

- **A UMAP is not evidence.** Do not annotate one with a p-value, and do not let cluster separation on a UMAP stand in for a tested effect. Inter-cluster distance on a UMAP is not interpretable.
- A violin or box plot over tens of thousands of cells implies an n it does not have. If the comparison is between conditions, plot the per-donor pseudobulk points on top, or plot the donor-level summary instead.
- PCA score plots require explained-variance labels and clear handling of grouping/metadata.
- Loading and regulon plots require feature/annotation provenance and are not causal evidence.
- Heatmaps need documented scaling, clustering, distance/linkage choices, and feature selection rules. State whether the genes shown were selected using the same data being displayed.
- Spatial plots need the scale, the FOV/slide boundaries, and an honest statement of what segmentation produced the cells.

<!--
Deliberately omitted: a general matplotlib/seaborn tutorial.
Modern coding models already know those APIs; loading a tutorial costs a
meaningful fraction of a session's context every time this skill triggers, with
no benefit. Keep this file limited to project-specific rules and conventions.
If a genuinely non-obvious plotting technique is needed repeatedly, add that
specific snippet — not a library overview.
-->

## Network figures

- State the PKN, sparsity parameter, solver, and number of sampled solutions in the caption or a linked table.
- Say whether each edge is a single-solution edge or a frequency across solutions, and encode the frequency (for example by width) if so.
- Encode edge sign with more than colour (line style or arrowhead), and say what node colour means (input, intermediate, measured, direction of change).
- A layout is not evidence. Do not read distance in a force-directed layout as biological distance.
