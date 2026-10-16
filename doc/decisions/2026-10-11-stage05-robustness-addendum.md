# Stage 05 robustness: addendum to the stage-04 plan

Date: 2026-10-11
Status: confirmed by user on 2026-10-12. Controls A-D remain post hoc by design.
Parent plan: `doc/decisions/2026-10-10-stage04-inference-plan.md` (confirmed), steps 16-21.
Code: `scripts/analysis/06_corneto_robustness.py`, `src/corneto_worker.py`; tests in `tests/test_robustness.py`.

The stage-04 plan was confirmed before stage 04 ran. This addendum records (i) where the stage-05
implementation differs from the plan text and why, and (ii) four controls added after the stage-04
results were inspected. Because the additions were designed after the results were seen, the report
must describe them as post hoc.

## Execution

- Every CORNETO solve runs in its own Python process (`src/corneto_worker.py`). Stage 04 was killed
  twice on 2026-10-11 when several COSMOS problems were solved in one process. Alternative-solution
  sampling needs the live CORNETO problem object, so a sampling task solves and samples in the same
  worker.
- Finished solves are reused with `--resume` when the task signature matches. The signature covers
  the PKN checksum, measurements, inputs, lambda, time limit, MIP gap and solver settings.
- The default run covers the OmniPath controls plus the simvastatin (SV) and receptor controls,
  about 3-4 h in total. COSMOS controls take 45-95 min per solve and run only with `--cosmos`:

```bash
python scripts/analysis/06_corneto_robustness.py                 # OmniPath, SV, receptors
python scripts/analysis/06_corneto_robustness.py --resume --cosmos
```

## Differences from the plan text

1. **Alternative solutions (step 16).** The plan asks for 30 alternatives for every accepted stage-04
   run. Each alternative is a full re-solve, so 30 COSMOS alternatives would take 20-45 h per run.
   Proposed: 30 alternatives for `omnipath`, and 10 for `cosmos_multi` with `--cosmos`
   (`--cosmos-samples` changes this). `cosmos_tf` and `cosmos_repair` are not sampled.
   `omnipath_lines` is skipped. Stage 04 selected lambda = 0 for it because no lambda > 0 found
   an incumbent within 1,800 s. Without the sparsity penalty every sign-consistent edge set is
   optimal, so alternatives would say nothing about stability.
   Edge frequencies are computed over the re-solved optimum plus every alternative within the 5 %
   objective tolerance. The re-solved optimum is compared with the stage-04 network
   (`resolve_edge_jaccard_with_stage04`), because stage 04 kept no problem objects.
2. **Leave one line out (step 17).** Applying the LOO stability filter again inside a LOO contrast
   would need 1-vs-3 contrasts, which the moderated t cannot estimate. The stage-04 six-LOO sign
   filter is applied to each LOO contrast instead. For `cosmos_multi` the metabolite measurements
   come from different cell lines and are kept unchanged.
3. **Permutations (step 19).** Inputs are free-sign sources, so most label assignments fit
   perfectly (error 0) and the fitting error cannot rank them. The table also reports two
   statistics computed before CORNETO: the number of LOO-stable TFs and the summed |score| of the
   selected TFs. It gives the rank of the true split for each statistic, with the number of ties.
   The mirror assignment (C4-C6 as controls) is labelled because its network should match the true one.
4. **Rewired PKNs (step 20).** Overlap is measured against the real stage-04 *selected* network
   (node and edge Jaccard). The draft script compared against the whole PKN, which is near zero by
   construction. COSMOS rewiring uses `pkn.rewire_type_preserving_signed_graph` (duplicate removed).
5. **Hub bias (step 21).** Measured TFs enter every network as outputs, so they are removed from the
   logistic regression of node selection on log2 degree.
6. **Core network.** Every edge's frequency now has the same denominator: the number of ensemble
   networks plus the number of accepted LOO networks. In the draft, an edge seen only in LOO
   networks got frequency 1/1.
7. **Choice grid (step 18).** For the Welch and logFC levels, the LOO sign filter still uses the
   moderated-t LOO contrasts, because stage 02 stores no per-statistic LOO. Unsolved grid points are
   reported as missing (NaN), not as Jaccard 0.

## Added controls

A. **Tighter-gap re-solve** (`--cosmos`, part `tight_gap`). `cosmos_multi` is re-solved at
   lambda 0, 0.01, 0.03, 0.1 and 0.3 with `mip_rel_gap = 1e-4` (3,600 s each). Reason: at the 1 %
   gap the stage-04 errors are not monotone in lambda (0.01: 10.57 > 0.03: 8.66), and the 10 %
   rule selected lambda = 0.1. The table reports both error sets and the lambda the rule picks
   from the tight solves. Lambda = 1 is not re-solved: its stage-04 error (28.8) is far outside the
   10 % band.

B. **Simvastatin positive control** (part `sv_control`). Known input HMGCR = -1. The single
   measurement is free cholesterol in ST003328, SV vs untreated, paired within each of the 7 lines.
   It is a one-sample moderated t, with the prior fitted to all 642 complete lipid species
   (`activities.one_sample_moderated_t`): mean log2 difference -3.98, t = -11.2. Lambda = 0.01.
   The control runs on four PKNs: full primary COSMOS, full repaired COSMOS, and the
   metabolism-only cofactor-free graph of each. The metabolism-only graph keeps metabolite-reaction
   edges and enzyme-own-reaction edges (`metabolic_route.metabolism_graph`) and drops the
   non-protected cofactors registered in the Aim-5 plan.
   Trial solves on 2026-10-11 (600 s limit, not part of the registered output):
   - Full primary and full repaired PKN: no incumbent within 600 s. After pruning, about
     14,000 nodes and 32,000 edges reach from HMGCR to cholesterol. The registered run uses 3,600 s.
   - Metabolism-only repaired PKN: accepted (gap 4.4 %), error 0, 52 edges. The HMGCR to cholesterol
     route passes through mevalonate and uses 15 of the 20 canonical genes (HMGCR, MVK, PMVK, MVD,
     FDPS, FDFT1, SQLE, LSS, CYP51A1, TM7SF2, NSDHL, DHCR24, EBP, SC5D, DHCR7), in canonical order
     (Kandutsch-Russell-like). So CORNETO recovers the known route once the signalling and
     transport shortcuts are removed and the FDFT1 gap is repaired.

C. **Receptor-only inputs** (part `receptors`). The OmniPath run at the stage-04 lambda, with inputs
   restricted to OmniPath intercell `receptor` proteins in the PKN (expressed by construction),
   excluding measured and excluded nodes. The resource is cached at
   `data/resources/pkn/omnipath_intercell_receptor.tsv` with a checksum sidecar
   (URL `https://omnipathdb.org/intercell?categories=receptor&format=tsv`). Not yet run, so the
   query and column filter (`category == parent == receptor`, protein entities) are untested
   against the live service.

D. **No compound-layer inputs** (`--cosmos`, part `no_drug_inputs`). In both multi-omic stage-04
   networks, one unmeasured free input carries about half the fit:
   - In `cosmos_multi` (lambda = 0.1), the tamoxifen node `Metab__5376_c` (PubChem CID 5376)
     explains 39 of 79 reachable measured nodes. Every route into SP1, JUN, JUNB and EGR1 starts
     there, through MAPK1/MAPK14.
   - In `cosmos_repair` (lambda = 0.1), cinnamaldehyde (`Metab__HMDB0003441_c`) takes this role
     and explains 44 of 79.

   Both are source metabolites whose out-edges go only to proteins. They belong to the COSMOS
   compound-target layer and take part in no PKN reaction. Of the 4,711 free inputs offered to
   stage-04 `cosmos_multi`, 3,683 are metabolites and 3,246 of those are compound-layer nodes.
   The control re-solves `cosmos_multi` twice, with inputs filtered in two ways; the removed
   nodes stay in the PKN:
   - (i) `compound_layer`: without the compound-layer metabolite inputs.
   - (ii) `all_metabolites`: without any metabolite input.

   Fit and network are compared with stage 04. The hub swap between `cosmos_multi` and
   `cosmos_repair` already shows that the upstream input is not identifiable from these data.

## Outputs

`results/ionescu_corneto/05_robustness/summary/`: `sampling_status.tsv`, `edge_frequencies.tsv`,
`loo_status.tsv`, `choice_grid.tsv`, `permutations.tsv`, `permutation_rank.tsv`, `rewired.tsv`,
`receptor_inputs.tsv`, `sv_cholesterol_contrast.tsv`, `sv_control.tsv`, `tight_gap.tsv`,
`no_drug_inputs.tsv`, `hub_bias.tsv`, `core_network.tsv`. Each solve's directory holds `task.json`,
`result.json`, `edges.tsv`, `nodes.tsv`, `fit.tsv` and, when sampled, `alternatives.tsv`.

## Changes to shared code

- `src/inference.py`: `make_corneto_data` / `solve_carnival` accept `inputs={node: value}`. These
  replace the default input set: 0 means unknown sign, ±1 a known perturbation. Explicit inputs may
  be internal nodes and bypass the exclusion list. Stage-04 behaviour is unchanged when `inputs` is
  omitted. `sample_edges` now reads every sample column, so multi-sample problems such as
  `omnipath_lines` are handled; the draft read only sample 0.
- `src/activities.py`: new `one_sample_moderated_t` (paired design).
- `src/corneto_worker.py`: new.

## Results of item D (run finished 2026-10-13, `summary/no_drug_inputs.tsv`)

| Variant | Inputs removed | Error (stage 04: 8.93) | TFs fitted | Metabolites fitted | Dominant input | Inferred-node Jaccard with stage 04 |
|---|---|---|---|---|---|---|
| `compound_layer` | 3,246 | 9.43 | 50/50 | 19/36 | LAMA5 (33 of 79) | 0.42 |
| `all_metabolites` | 3,683 | 11.66 | 50/50 | 16/36 | EFNB3 (25 of 77) | 0.36 |

- **Single-hub structure:** removing the compound starting nodes keeps it; the hub moves to a protein.
- **Kinase core:** MAPK1, MAPK14, PRKACA, FYN and AKT1 stay in both variants.
- **Cholesterol:** explained by the HDLBP transport reaction in both, and HMGCR is never active.
- **Cholesterol to TF routes:** these appear only once the compound hubs are gone:
  - `compound_layer`: cholesterol → MAPK3 → CEBPB → STAT3 → JUNB;
  - `all_metabolites`: cholesterol → MAPK3 → NCOA1 → FOS → JUN → EGR1.

## Item E: late versus early integration (added 2026-10-13, `--cosmos`, part `late_integration`)

- **Late network:** one metabolite-only CORNETO solve on COSMOS (36 metabolite measurements, default
  free inputs) at the joint run's lambda (0.1), merged with the stage-04 `cosmos_tf` network at the
  same lambda. This merge is the late-integration network.
- **Early network:** the stage-04 `cosmos_multi` network at the same lambda.
- **Readouts** (`summary/late_integration.tsv`), for the TF-only, metabolite-only, merged and joint
  networks:
  - edges, nodes and weakly connected components;
  - edges crossing between the metabolic (metabolite or reaction) and protein layers;
  - measured-metabolite-to-measured-TF directed paths;
  - for the merge: node-state conflicts between the two networks and edge Jaccard with the joint
    network;
  - for the joint network: edges found in neither single-layer network.
- **Network properties** without new solves: `scripts/analysis/07_network_properties.py` writes to
  `results/ionescu_corneto/07_downstream/network_properties/`:
  - degree bias of node selection, for real and randomised networks;
  - out-degree percentiles of the starting nodes used;
  - centrality of the kinase core and the hypothesis nodes;
  - metabolite-to-TF distances.
