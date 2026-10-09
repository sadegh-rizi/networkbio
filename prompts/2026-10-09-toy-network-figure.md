# Toy CORNETO network figure and Cytoscape tables

Date: 2026-10-09
Model: claude-opus-5-5 (Claude Science)
Settings: default
Written in a temporary worktree (`.worktrees/toy-network-viz`, branch
`viz/toy-network` from `7b1fa51`) while another agent was editing
`task/stage01-02`. The worktree's link file could not be used from WSL, so the
files were copied into the main tree and committed on `task/stage01-02`.
Commit(s): none yet (add after the user commits)
Transcript: not saved

## Prompt actually used

> can you do it using the toy_edges? (be careful that another agent is
> implementing the plan in the directory, so you can add another working tree
> in the repos and work there)

## What it produced

- `src/network_plot.py`: reusable drawing of a CORNETO/CARNIVAL solution on
  its signed PKN (layered layout, node states from edge values with a
  consistency check, edge-frequency table, legend handles).
- `scripts/analysis/00b_plot_toy_network.py`: reads
  `results/corneto_toy/04_inference/toy_edges.tsv` and `toy_checks.json`,
  takes the PKN/measurements from the constants in `00_toy_corneto.py` (parsed
  with `ast`, CORNETO not imported), writes `figures/toy_network.*` and the
  Cytoscape tables `toy_network_edges.tsv` (36 rows = 12 PKN edges x 3 runs)
  and `toy_network_nodes.tsv` (48 rows).
- `tests/test_network_plot.py`: 9 synthetic tests.

## Verification

Run on Windows with the sandbox Python (matplotlib 3.11.2, pandas 2.3.3), not
the project `.venv`: script exit 0, 0 measured-sign mismatches in all three
panels, `pytest tests/test_network_plot.py tests/test_inputs.py` 14 passed.
The figure was inspected visually. The committed outputs should be regenerated
in WSL with `.venv/bin/python scripts/analysis/00b_plot_toy_network.py`.

## Notes for later audits

- `toy_checks.json` reports `regularization_edge_has_signal = 8` while 7 PKN
  edges are selected; the extra one is presumably the boundary (input) edge
  CORNETO adds for R1. Check how boundary edges enter the λ penalty before
  comparing edge counts across real runs.
- The toy has 12 nodes, not 13 as `doc/analysis-state.md` says.

## What was changed by hand afterwards

_To be filled by the user._

## Outputs stamped

None.
