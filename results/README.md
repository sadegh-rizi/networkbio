# results/

Derived outputs. Nothing here is a source of truth: everything is
regenerable from `data/` plus the code and the plan that produced it.

## Layout

Never write into a flat directory. Group each analysis into numbered stages:

```
results/<analysis>/
  00_inputs/          # one-time conversions of raw tables + provenance sidecars
  01_preprocessing/
  02_activities/
      figures/
  03_pkn/
  04_inference/
      figures/
  05_robustness/
  06_validation/
  07_downstream/
```

Tables sit at the stage level; that stage's figures go in its `figures/`
subfolder.

## Variants go in the path, not the filename

When outputs vary along an experimental axis (PKN version, sparsity
parameter, solver, condition set, layer subset), put that axis in the
directory path and keep filenames identical across variants.

## Figures

Each figure is three files in the same folder: `<name>.pdf`,
`<name>_presentation.png`, `<name>_manuscript.png`, all written by
`plot_config.save_figure`.

## Large objects

Pickles, parquet files and solver outputs are gitignored. Each one needs a
committed `<name>.provenance.json` sidecar recording the source accession,
the code commit, the parameters, the seed, the solver and its version, and the
package versions. A file with no sidecar is an orphan.
