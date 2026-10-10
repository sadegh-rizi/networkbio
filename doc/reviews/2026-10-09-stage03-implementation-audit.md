# Audit: stage 03 implementation (PKNs, mapping, Aim-5 structure check)

Date: 2026-10-09
Auditor: claude-opus-5-5 (Claude Science); implementation by another agent
(gpt-5.6-luna), commit `75597e3` on `task/stage03-pkn`.
Scope: `doc/decisions/2026-10-09-stage03-pkn-plan.md` against `src/pkn.py`,
`scripts/analysis/04_build_pkn.py`, `tests/test_pkn.py`, the cached resources
in `data/resources/pkn/` and `data/resources/metabolite_ids/`, and the outputs
in `results/ionescu_corneto/03_pkn/` (written 2026-10-09 15:14-15:16; a second
run by the user was still in progress at 15:38, `logs/stage03.log`).
Blinding: stage 03 reads no activity values; none were looked at here.

## Verdict

The OmniPath half works and is correct. The COSMOS half does **not** work:
the COSMOS node parser assumes a node grammar that the selected file does not
use, so expression pruning keeps or removes enzyme nodes according to
unrelated genes, and numeric metabolite IDs are turned into wrong HMDB IDs.
The COSMOS PKN, its Aim-5 rows, its null and its figure must be regenerated
after the fix (B1-B5). Separately, the plan's 8-edge cap cannot reach
cholesterol through the metabolic layer (P-a), which is the planner's error,
not the implementer's.

## What was checked

| Item | How | Result |
|---|---|---|
| P2 direction/sign (OmniPath) | code | primary keeps consensus direction AND exactly one consensus sign; `nc_rule` keeps either and gives +1 to two-sign rows. Correct |
| P3-P6 (OmniPath) | code, `omnipath/primary/filter_log.tsv` | 72,272 rows -> P2 71,116 (134 two-sign, 1,022 without consensus direction) -> P3 15,278 -> complexes 14,640 -> duplicates 14,448 -> expression 9,870 edges / 3,443 nodes. Order and logging as planned |
| P5 symbols | code | approved, then unique previous symbol, then unique alias; re-deduplicated afterwards. Correct |
| Signed BFS | code | breadth-first search on (node, sign) states; distances and shortest-walk counts are correct; "passes a metabolite" is propagated along shortest walks. See caveat C1 |
| Rewiring (P10) | code | directed double-edge swaps; each edge keeps its sign; in/out degree and sign counts preserved; self-loops and duplicates rejected; seeded. Correct |
| P8 mapping | code, `mapping/metabolite_mapping.tsv` | Recon3D first, Workbench fallback, compartments c/m (cells, lipids) and e (medium), labelling class carried from the stage-01 feature tables. D-glucose maps to HMDB0000122 through Recon3D; free cholesterol to HMDB0000067. Correct, apart from B2 |
| COSMOS source check (P7) | provenance | both files fetched; after clean-up they differ (81,788 vs 65,138 edges). They are different builds in different namespaces (below). The plan's "must give the same edge set" did not hold; the implementation recorded this and used the NetworkCommons file |
| COSMOS grammar (R1) | raw files | see B1 |
| Blinding | code | only `regulon_summary.tsv` (names, target counts) is read from stage 02. Correct |

The two COSMOS downloads:

| File | Genes | Enzyme/reaction nodes | Metabolites |
|---|---|---|---|
| NetworkCommons `meta_network.sif` (selected; 83,557 rows) | gene symbols (`HMGCR`) | `Gene<reaction no.>__<SYMBOL or pseudo-name>[_reverse]`, e.g. `Gene1600__HMGCR`, `Gene10001__SLC7A6_TRANSPORTER1` | `Metab__HMDB0000067_c` (3,491 nodes) and numeric `Metab__2676_c` (3,970 nodes) |
| OmniPath `metapkn__20200122.txt` (67,042 rows) | `X<Entrez>` | `XGene<reaction no.>__<Entrez or pseudo-name>` | `XMetab__227___c____` (bare HMDB number) |

Neither matches the 117,065 interactions in the cosmosR documentation, so
the NetworkCommons file's build date is unknown; record it as "NetworkCommons
`meta_network.sif`, downloaded 2026-10-09, sha256 77468042...", not as the
2020 build.

## Bugs (implementation)

**B1. COSMOS enzyme nodes are assigned to the wrong genes.** `parse_cosmos_node`
matches `(?:X)?Gene(\d+)__(.+)` and looks up the number as an Entrez ID. In the
selected file the number is a reaction index and the gene is the suffix
(`Gene1600__HMGCR`; Entrez 1600 is a different gene). Examples from
`cosmos/primary/nodes.tsv`: `Gene10001__SLC7A6_TRANSPORTER1 -> MED6`,
`Gene10007__... -> GNPDA1`, `Gene10013__... -> HDAC6`. P6 then keeps or drops
each of the 19,530 enzyme nodes according to that unrelated gene (6,772 kept).
The raw file does connect the cholesterol axis correctly
(`HMGCR -> Gene1600__HMGCR(_reverse) -> Metab__HMDB0000227_c`); the bug
breaks it. Fix: take the gene from the suffix; strip `_reverse`; split on `_`
and keep the tokens that are HGNC approved symbols after P5; nodes with no
gene token (transport, exchange and demand pseudo-reactions such as `EX_...`,
`DM_...`, `...EXCHANGE1`) get type `transport_or_pseudo` (decision D-1).
Also make the parser raise on `X`-prefixed nodes, so the other grammar cannot
be read silently.

**B2. Numeric metabolite IDs are converted into HMDB IDs.** `Metab__2676_c`
becomes `Metab__HMDB0002676_c`. These 3,970 numeric IDs are not HMDB numbers
(most likely PubChem CIDs from STITCH; unverified), so the conversion creates
false identities. Three measured features mapped only because of it:
`Metab__HMDB0002142_c`, `Metab__HMDB0002347_c`, `Metab__HMDB0004984_c` do not
occur in the raw file. Fix: keep numeric IDs unchanged, typed
`metabolite_numeric_id`, and never treat them as HMDB. Optional: map measured
features to them through the PubChem CID the Workbench lookup already returns.

**B3. Provenance cannot tie outputs to code.** The sidecar records
`code_commit 6f3ff82` (the merge), but the code was uncommitted at run time and
`src/pkn.py` was edited again at 15:25, after the outputs (15:16). There are no
code hashes (stages 01-02 have them), and the inputs omit the stage-01 feature
tables and the stage-02 regulon summaries that were read. Fix: reuse
`preprocess.stage_provenance` / `reuse_stage` as stages 01-02 do.

**B4. Tests encode the assumed grammar, not the real one.** The synthetic tests
passed because they used the same wrong grammar. Add a fixture of real lines
from the NetworkCommons file (`HMGCR 1 Gene1600__HMGCR`,
`Gene1600__HMGCR_reverse 1 Metab__HMDB0000227_c`, a numeric metabolite, a
transporter pseudo-node) and assert: the enzyme gene is HMGCR, the
HMGCR -> mevalonate path has length 2, the numeric node stays numeric.

**B5. The figure is empty.** With no paths, both figures show only "No path
within the eight-edge cap". Expected until B1 and P-a are fixed.

## Plan problems (planner's errors)

**P-a. An 8-edge cap cannot reach cholesterol in COSMOS.** Each metabolic
reaction is two edges (metabolite -> enzyme node -> metabolite). Mevalonate to
cholesterol takes many reactions, so HMGCR -> cholesterol is far longer than 8
edges, while cholesterol has direct signalling edges (for example to AKT1,
MAPK1, CAV1 in the raw file). Replace the single capped query by three legs:
(i) HMGCR enzyme -> mevalonate and -> cholesterol (metabolic leg, cap 80,
report length); (ii) cholesterol and mevalonate nodes -> the six TFs (cap 8);
(iii) SREBF1/2, SCAP, INSIG1 -> the six TFs (cap 8, as now). The length of leg
(i) is itself an Aim-5 result: CORNETO's sparsity penalty charges every edge,
so a route this long is strongly disfavoured, whatever the data say.

**P-b. Expression pruning removes WT1.** WT1 is not among the 13,530 expressed
symbols, so it is pruned from both PKNs, although CollecTRI scores its
activity from its targets. WT1 can then never be a measured node in stage 04.
Decision D-2.

**C1 (caveat, not a bug).** A shortest walk of a given sign in the (node, sign)
state graph can revisit a node with the other sign, so "positive shortest
length" may belong to a non-simple walk. Report whether each listed path is
simple (check in `path_edges.tsv` once paths exist).

## Decisions needed before the fix

- **D-1. Enzyme nodes without a gene token** (transport, exchange, demand
  pseudo-reactions): keep them like metabolites (recommended: they carry
  metabolite flow between compartments and are not genes), or drop them.
- **D-1b. Multi-gene enzyme nodes:** keep if any token is expressed
  (recommended; isoenzymes) or only if all are.
- **D-2. WT1 and other unexpressed TFs:** keep pruning as registered and report
  WT1 as structurally unreachable, with the `unpruned` variant as the check
  (recommended; changing the rule after seeing which TF it removes is the
  kind of post-hoc edit the plan forbids), or exempt regulon TFs from pruning
  for all TFs (a new rule, registered before rerunning).
- **P-a** as above (three legs).

## Results that stand

- **OmniPath: no signed path from the cholesterol axis to any of the six TFs.**
  Path search was run on the primary PKN; the out-edges were also checked in
  the `unpruned` variant (same picture), and `signor_only` is a subset of the
  primary. `nc_rule` (134 extra two-sign rows) was not checked. HMGCR, SREBF1
  and SREBF2 have no outgoing edges; SCAP points only to SREBF1/2; INSIG1 only to SCAP and SERBP1.
  The axis is a sink in a protein-signalling PKN. This is a real Aim-5
  result: without a metabolite layer, CARNIVAL/CORNETO cannot reach the
  proposed route from HMGCR, because the PKN has no such edges.
- **TF coverage:** of the 628 CollecTRI TFs kept in stage 02, 330 are nodes
  of the OmniPath primary PKN and 349 of the (buggy) COSMOS primary. E2F1,
  EGR1, SP1, JUN and JUNB are present in both; WT1 in neither (P-b).
- **Metabolite mapping** (to the buggy COSMOS node set; will change slightly):
  distinct features mapped ST003331: 21 `m0_unconfounded`, 9
  `isotopologue_sum_partial`, 82 `m0_confounded`; ST003332: 12
  `m0_unconfounded`, 50 `m0_confounded`; lipidomics: free cholesterol. The
  plan's "fewer than about 10 usable metabolites" failure condition is not
  met (30 usable in cells, 12 in medium).

## Not checked

The tests other than for the grammar point; the figures beyond their empty
content; the rewiring runtime; the OmniPath variants beyond their filter logs.
