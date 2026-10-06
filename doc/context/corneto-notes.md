# CORNETO notes

Last verified: 2026-10-06 (API section); 2026-09-29 (the rest). Everything here was read in search results and
abstracts during planning, not in the package documentation or code. Treat
the API-level statements as **unverified** until checked against the pinned
version.

## Read in the paper and preprint abstracts

- CORNETO is a framework that reformulates a range of network inference
  methods as mixed-integer optimisation problems using network flows and
  structured sparsity, so that several samples can be inferred jointly.
  It supports undirected, directed and signed (hyper)graphs.
- It reimplements and extends CARNIVAL to the multi-sample case. CARNIVAL
  infers signalling networks from TF activities estimated from
  transcriptomics and known TF targets, and reconstructs networks from known
  inputs such as drug targets or receptors.
- Reported benefit: with multi-sample data, sparser networks at the same
  fitting error, capturing shared interactions; less variability across
  solver seeds than single-sample runs.
- Metabolic use cases (iMAT-style and flux-balance-type problems) are
  demonstrated separately from signalling.

Sources: Nature Machine Intelligence 2025,
https://www.nature.com/articles/s42256-025-01069-9 ; bioRxiv 2024,
https://www.biorxiv.org/content/10.1101/2024.10.26.620390v1 ; repository
https://github.com/saezlab/corneto ; docs https://corneto.org/.

## Seen in a documentation tutorial snippet (not run)

The multi-condition CARNIVAL tutorial builds a data object where each
sample maps node names to `value`, `role` (input or output) and `mapping`;
inputs with unknown direction get value 0; a `lambda_reg` parameter sets the
sparsity, and the tutorial chooses it by sampling multiple solutions and
scoring agreement with the data. Source:
https://corneto.org/v1.0.0-beta.3/tutorials/contrib/multi_condition_tutorial/multi_condition_carnival.html

## Verified in the installed code (corneto 1.0.0rc8, 2026-10-06)

- `corneto.methods.CarnivalFlow(lambda_reg=..., max_flow=1000, ...)` is the
  multi-sample CARNIVAL; `build_from_data(graph, data)` (plain `build` is
  deprecated) then `P.solve(solver="highs", max_seconds=...)`.
- Graph from signed SIF tuples: `corneto.io.load_graph_from_sif_tuples`.
- Data: `cn.Data.from_cdict({sample: {node: {"value", "role": "input"|"output", "mapping": "vertex"}}})`.
  An input with value 0 has its sign chosen by the solver (checked in the toy).
- Preprocessing prunes the graph to vertices on input→output paths per sample
  and adds boundary edges; selected edges are `P.expr.edge_value` on
  `method.processed_graph` (|value| > 0.5).
- A metabolite node (`Metab__HMDB…_c`) is handled like any other node
  (`results/corneto_toy/04_inference/toy_checks.json`).

## Unverified

- Run time and memory on the COSMOS meta-PKN with real inputs.
- Which solver is practical at real scale (HiGHS and SCIP both installed;
  16 GB RAM, WSL limit not measured for this project).
