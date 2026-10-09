# Stage 03: prior-knowledge networks, identifier mapping, Aim-5 structure check

Date: 2026-10-09
Status: proposed

Implements step 4 of `2026-09-29-corneto-ionescu-inscs-plan.md` and settles
its open question Q6 (PKN), split below into P1-P10. No CORNETO run, no λ and
no TF activity values are used in this stage.

## Question

1. Build two signed, directed PKNs and log every filter: a protein-only
   signalling PKN (OmniPath) and a signalling + metabolism PKN (COSMOS
   meta-PKN).
2. Map the measured features onto PKN nodes: TFs from the stage-02 regulons
   (names only), metabolites by KEGG ID, cholesterol from the lipidomics.
3. Aim 5, first half: does a signed path from the cholesterol-synthesis axis
   to the six cholesterol-dependent TFs exist in each PKN, how short is it,
   does it pass through a metabolite node, and how unusual is that compared
   with degree-preserving rewired PKNs?

## Data available

- OmniPath interactions through the `omnipath` 1.0.12 client in `.venv`:
  `omnipath.interactions.OmniPath.get(genesymbols=True)` (dataset
  `omnipath`; columns include `consensus_direction`,
  `consensus_stimulation`, `consensus_inhibition`, `curation_effort`).
- COSMOS meta-PKN (OmniPath signalling + STITCH + Recon3D; 117,065 signed
  interactions in the cosmosR documentation). Two copies of the same 2020
  build: NetworkCommons downloads `prior_knowledge/meta_network.sif` from its
  data server (`networkcommons/data/network/_moon.py`; server
  `https://commons.omnipathdb.org/`, URL not opened from the sandbox), and
  OmniPath hosts `https://metapkn.omnipathdb.org/metapkn__20200122.txt`
  (blocked in the sandbox, not opened). Metabolite nodes look like
  `Metab__HMDB0000067_c` (HMDB ID plus compartment); the grammar of enzyme
  and gene nodes is **not yet verified** (step R1).
- Recon3D model with metabolite cross-references (BiGG,
  `http://bigg.ucsd.edu/static/models/Recon3D.json`; not opened yet). This is
  the namespace the COSMOS metabolic layer was built from.
- Metabolomics Workbench REST `compound/kegg_id/<KEGG>/all` (checked
  2026-10-09): cholesterol C00187 -> HMDB0000067, mevalonate C00418 ->
  HMDB0000227, 2-oxoglutarate C00026 -> HMDB0000208, but D-glucose C00031 ->
  HMDB0304632 (not the HMDB0000122 used in most models). One ID per KEGG
  compound, so it is a fallback only.
- Stage 01: 13,530 expressed symbols (`rnaseq/log2cpm.tsv` index);
  metabolomics `features.tsv` with `kegg_id` and `labelling_class` (ST003331
  158 features, ST003332 136); lipidomics names only (free cholesterol and 30
  cholesteryl esters).
- Stage 02: `regulon_summary.tsv` per resource (TF names and target counts;
  no activity values are read in this stage).
- HGNC complete set (`data/resources/hgnc/`) with `prev_symbol` and
  `alias_symbol`.

Not available: OmnipathR's newer `cosmos_pkn()` build (R only, 30-40 min
build). Out of scope for the course.

## Unit of observation

Not a statistical stage. The structural null (P10) treats one rewired PKN as
one draw.

## Decisions to confirm (recommended default first)

**P1. Signalling PKN source.** OmniPath `omnipath` dataset (literature-curated
causal interactions from many resources). Variant: SIGNOR-only subset of the
same table (rows whose `sources` include SIGNOR). The main plan's "SIGNOR
first" is read as "SIGNOR as the sensitivity".

**P2. Direction and sign.** Keep rows with `consensus_direction` true and
exactly one of `consensus_stimulation` / `consensus_inhibition` true; sign
+1 or -1 accordingly. Rows with both signs are dropped and counted.
Variant `nc_rule`: the NetworkCommons rule (either sign true; sign +1 if
stimulation is true, so two-sign rows become +1).

**P3. Evidence.** `curation_effort >= 2` (the NetworkCommons default).
No variant (one sensitivity per decision is enough for two weeks).

**P4. Complexes and identifiers.** Gene symbols; drop complex nodes
(`COMPLEX:` UniProt IDs, or symbols joined with `_`), and count the dropped
edges. Duplicate source-target pairs that disagree in sign are dropped and
counted; identical duplicates are merged. Self-loops are dropped.

**P5. Symbol harmonisation.** PKN symbols that are HGNC approved symbols are
kept. Otherwise map through `prev_symbol`, then `alias_symbol`, only when the
mapping is to exactly one approved symbol; otherwise leave the node as is
(it is then removed by P6 if not expressed). Log each count.

**P6. Expression pruning.** Remove gene/protein nodes whose approved symbol
is not among the 13,530 expressed symbols (the stage-01 filter, CPM >= 1 in at
least 3 libraries), as COSMOS does. Metabolite nodes are not pruned. Variant
`unpruned`. For COSMOS enzyme or reaction nodes, the gene symbol is parsed
from the node name (R1) and the node is removed if that gene is not
expressed.

**P7. COSMOS meta-PKN version and clean-up.** The 2020 build as distributed by
NetworkCommons (`meta_network.sif`), cleaned as NetworkCommons
`meta_network_cleanup` does: drop self-loops, average the sign of duplicate
source-target pairs, keep only edges whose result is exactly +1 or -1. Then
P5 and P6 on its gene nodes. If the NetworkCommons server is unreachable from
WSL, use the OmniPath-hosted `metapkn__20200122.txt` and record which was
used. The two must give the same edge set after clean-up; this is checked if
both can be fetched.

**P8. Metabolite mapping.** KEGG -> HMDB from the Recon3D metabolite
annotations first (all HMDB IDs listed for that KEGG ID), then Metabolomics
Workbench as a fallback; each mapping records its source. HMDB IDs are
normalised to the 7-digit form (`HMDB0000067`). Compartments: cell extracts
(ST003331) map to `_c` and `_m` nodes (the NetworkCommons default `["c", "m"]`);
medium (ST003332) maps to `_e`. Lipidomics: only free cholesterol is mapped
(name match to `Cholesterol`, HMDB0000067); cholesteryl esters and other
species are listed as unmapped. Every mapping row carries the feature's
`labelling_class`; stage 04 may only use `m0_unconfounded` and
`isotopologue_sum_partial` features as measurements (revision plan M6).

**P9. Aim-5 structure check.**
- Sources: HMGCR, SREBF1, SREBF2, SCAP, INSIG1. Targets: E2F1, EGR1, WT1, SP1,
  JUN, JUNB. In COSMOS also report the route through the mevalonate/cholesterol
  metabolite nodes (HMDB0000227, HMDB0000067, any compartment).
- For each source-target pair: shortest directed path length, the signs it can
  have (breadth-first search on (node, sign) states, so a path's sign is the
  product of its edge signs), the number of shortest paths of each sign, and
  whether any shortest path passes a metabolite node. Path length capped at
  8 edges.
- Registered now: this plan reads Ionescu et al. 2024 as predicting a net
  positive sign from HMGCR to the six TFs (statin lowers cholesterol, which
  lowers their activity). The user confirms this reading against the paper
  before confirming the plan. "Route present" means at least one positive-sign path of
  length <= 8. Report every pair either way; do not tune P1-P8 after seeing
  this table.
- No activity values in this stage: the figure colours nodes by type
  (source, target, metabolite, other), not by stage-02 scores.

**P10. Structural null.** Degree-preserving, sign-preserving rewiring: repeated
directed double-edge swaps (A->B, C->D becomes A->D, C->B; each edge keeps its
own sign; swaps creating self-loops or duplicate edges are rejected),
3 x |E| accepted swaps per replicate. 100 replicates for the OmniPath PKN, 20
for COSMOS (larger; runtime not yet measured). Seeds 20261009 + replicate
index. For each pair in P9: fraction of replicates whose shortest
positive-sign path is no longer than the observed one. Rewired PKNs are not
stored; stage 04 regenerates them from the same seeds and code.

## Method

R1. Inspect before coding the parser: download the COSMOS meta-PKN and
    record in `03_pkn/summary/cosmos_node_grammar.md` the node-name patterns
    (metabolites, genes, enzymes, reaction or direction suffixes) with counts
    and five examples each. If gene symbols cannot be parsed unambiguously,
    stop and report.
S1. Fetch and cache resources with checksum sidecars, the pattern of
    `src/activities.py::load_regulon`: OmniPath interactions, COSMOS
    meta-PKN, Recon3D JSON, and the Metabolomics Workbench lookups for the
    KEGG IDs present in ST003331/ST003332 (one request per ID, cached as a
    table). Location `data/resources/pkn/` and `data/resources/metabolite_ids/`.
S2. OmniPath PKN: P2 -> P3 -> P4 -> P5 -> P6, logging nodes and edges after each
    step (`filter_log.tsv`). Variants: `signor_only`, `nc_rule`, `unpruned`
    (each changes one step).
S3. COSMOS PKN: P7 clean-up -> P5 -> P6, same logging. Variant `unpruned`.
S4. Node tables: node, type (gene, metabolite, enzyme or reaction, other),
    approved symbol, expressed, in/out degree.
S5. Mapping tables (P8, plus TF coverage): which stage-02 TFs (names from
    `regulon_summary.tsv`, kept after tmin) exist as nodes in each PKN, per
    resource; for AP1 and NFKB, which member genes exist.
S6. Aim-5 path table and null (P9, P10).
S7. Figure: the union of the shortest positive and negative paths from P9, per
    PKN, drawn with `src/network_plot.py` (path edges highlighted, other edges
    among those nodes dashed). Caption states PKN, filters, path cap, and
    that no data values are shown.

## Open questions

None for this stage beyond confirming P1-P10. Left to the stage-04 plan: λ
grid and selection rule (Q7), solver time limit and gap (Q8), measurement
value encoding (score or sign), the AP1/NFKB sensitivity (addendum E4),
COSMOS run direction (one problem vs signalling->metabolism and
metabolism->signalling), and the label-permutation control.

## Files to change

New:
- `src/pkn.py`: fetch/cache, filters with logging, symbol harmonisation,
  pruning, COSMOS name parsing, metabolite mapping, signed BFS, rewiring.
- `scripts/analysis/04_build_pkn.py` (stage 03).
- `scripts/setup/run_stage03.sh` (pytest, then the stage, like
  `run_stage01_02.sh`).
- `tests/test_pkn.py`.
- `prompts/<date>-stage03-implementation.md`.

Edit: `doc/analysis-state.md`, `doc/datasets.md` (resource versions),
`doc/decisions/README.md`, `prompts/README.md`, this plan's status line.

Unchanged: stage 00-02 code and outputs, `data/` except the new
`data/resources/` subfolders.

## Outputs

```
results/ionescu_corneto/03_pkn/
  03_pkn.provenance.json
  omnipath/primary/{pkn.tsv,nodes.tsv,filter_log.tsv}
  omnipath/variants/{signor_only,nc_rule,unpruned}/{pkn.tsv,nodes.tsv,filter_log.tsv}
  cosmos/primary/{pkn.tsv,nodes.tsv,filter_log.tsv}
  cosmos/variants/unpruned/{pkn.tsv,nodes.tsv,filter_log.tsv}
  mapping/{metabolite_mapping.tsv,lipid_mapping.tsv,tf_node_coverage.tsv}
  aim5/{paths.tsv,path_edges.tsv,null_summary.tsv}
  summary/{counts.tsv,cosmos_node_grammar.md}
  figures/aim5_paths_{omnipath,cosmos}.{pdf,_presentation.png,_manuscript.png}
```

`pkn.tsv` columns: `source`, `sign`, `target` (loadable with
`corneto.io.load_graph_from_sif_tuples`).

## Validation

Synthetic tests (no network):
- P2: a row with both consensus signs is dropped under the primary rule and
  becomes +1 under `nc_rule`; a row without consensus direction is dropped.
- P4: complex nodes, self-loops and sign-conflicting duplicates are removed
  and counted; identical duplicates are merged.
- P5: a previous symbol mapping to one approved symbol is renamed; an alias
  shared by two approved symbols is not.
- P6: an unexpressed gene node is removed together with its edges; a
  metabolite node with the same degree is kept.
- P7: duplicate pairs with signs +1 and -1 average to 0 and are dropped.
- P8: one KEGG ID with two HMDB IDs maps to both, in each listed compartment;
  the mapping source is recorded; the labelling class is carried over.
- P9: on a toy graph with an inhibition edge, the BFS returns the right
  shortest length and sign, finds paths of both signs when both exist, and
  respects the cap.
- P10: rewiring keeps every node's in- and out-degree, the number of +1 and
  -1 edges, and creates no self-loops or duplicates; the same seed gives the
  same graph.

Real run (user, in WSL):

```bash
bash scripts/setup/run_stage03.sh 2>&1 | tee logs/stage03.log
```

Checks printed and written to `summary/counts.tsv`: node and edge counts after
every filter for every PKN; how many of HMGCR, SREBF1/2, SCAP, INSIG1 and the
six TFs survive pruning; metabolite mapping coverage per study and labelling
class; TF coverage per resource.

What would show the approach is wrong or uninformative:
- One of the five source genes or six TFs is not expressed or not in a PKN
  (then Aim 5 cannot be tested in that PKN; report, do not add the node).
- Signed paths exist for every pair in nearly all rewired PKNs at the same or
  shorter length. Then path existence carries no evidence in a PKN this
  dense, and Aim 5 has to rest on what stage 04 selects. This is a reportable
  result, not a failure.
- Fewer than about 10 `m0_unconfounded` or `isotopologue_sum_partial`
  metabolites map to COSMOS nodes. Then the metabolic layer adds little as
  measurements, and that limits Aim 4.
