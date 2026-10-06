---
name: python-analysis
description: Use when adding or refactoring reusable Python analysis logic, parsing source data files, or moving repeated notebook code into modules for the networkbio project.
---

# Python analysis

Read `AGENTS.md` and `doc/agent-rules.md` first. The project runs in WSL2 with a venv that may not yet contain the packages you want (see `AGENTS.md`). Never assume an import works.

## Workflow

1. Inspect the input schema and existing path handling before editing. Do not infer column semantics from names alone; `doc/datasets.md` lists the traps.
2. Probe schema cheaply: `pandas.read_csv(..., nrows=5)` or `usecols=`. Never load a large table to find out what its columns are.
3. Put reusable logic in `src/` as plain functions with explicit inputs and outputs. Notebooks orchestrate and plot.
4. Map every sample column to a line and condition explicitly, in code, and fail loudly on an unmapped column.
5. Keep track of feature counts through each filter and log them; do not drop features silently.
6. Set and record seeds. Write results to `results/<analysis>/<NN_stage>/` with a `.provenance.json` sidecar for any non-trivial intermediate.
7. Test with small synthetic fixtures in `tests/`, not the real data.

## Do not

- Add a dependency without asking. Do not `pip install`.
- Hard-code absolute paths.
- Overwrite anything in `data/`.
