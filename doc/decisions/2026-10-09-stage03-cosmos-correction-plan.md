# Stage 03 COSMOS correction and Aim-5 scope plan

Date: 2026-10-09
Status: implemented (commit c7f3719)

## Question

Can the stage-03 COSMOS PKN and Aim-5 structure analysis be corrected so that
the NetworkCommons node grammar is interpreted from the actual resource, the
derived PKN is reproducible from all files it reads, and the proposed
cholesterol route is tested at a biologically meaningful path scale?

## Findings motivating the correction

- The selected NetworkCommons resource uses reaction-index nodes such as
  `Gene1600__HMGCR`; `1600` is not an Entrez identifier. The gene symbol is in
  the suffix.
- Numeric nodes such as `Metab__2676_c` are model identifiers, not HMDB IDs,
  and must not be rewritten as `HMDB0002676`.
- The original eight-edge query cannot traverse the full metabolic route from
  an HMGCR reaction to cholesterol.
- The stage-03 provenance sidecar must include the stage-01 feature tables and
  stage-02 regulon summaries read by the script, in addition to resources and
  the expression table.

## Confirmed decisions

- Keep COSMOS transport, exchange, demand, and other pseudo-reaction nodes
  without a recognized gene token; label them `transport_or_pseudo` and do not
  expression-prune them.
- For multi-gene enzyme nodes, retain the node when any recognized gene token
  is expressed.
- Keep the registered expression-pruning rule. WT1 remains absent from the
  primary PKN when it is not expressed; report this and retain the unpruned
  variant as the sensitivity check.
- Replace the single Aim-5 query with three registered legs:
  1. source enzyme to mevalonate and cholesterol, cap 80;
  2. mevalonate and cholesterol to the six target TFs, cap 8;
  3. the five source genes to the six target TFs, cap 8.
- The current WT1 and path-cap decisions are structural analysis choices, not
  claims that the route is biologically present.

## Implementation

Files to change:

- `src/pkn.py`: correct NetworkCommons parsing, numeric metabolite typing,
  pseudo-reaction handling, multi-gene expression checks, and grammar report.
- `scripts/analysis/04_build_pkn.py`: implement the three Aim-5 legs, update
  path figures and settings, and use complete stage input provenance.
- `tests/test_pkn.py`: add real NetworkCommons grammar fixtures and assertions
  for HMGCR, numeric metabolites, pseudo-reactions, and multi-gene pruning.
- `doc/decisions/2026-10-09-stage03-pkn-plan.md`: register the corrected
  grammar and three-leg Aim-5 definition.
- `doc/analysis-state.md`: record that the initial COSMOS result was invalid,
  the correction, and the regenerated counts.
- `doc/reviews/2026-10-09-stage03-implementation-audit.md`: record resolution
  status and the final verification result if this existing audit is retained.
- `prompts/2026-10-09-stage03-cosmos-correction.md`: record AI provenance.

No dependency, source `data/` file, stage-00/01/02 code, or existing stage-00
through stage-02 output will be modified.

## Outputs

Regenerate all files under `results/ionescu_corneto/03_pkn/`, including both
COSMOS variants, mapping tables, the three-leg Aim-5 tables and figures,
filter summaries, and `03_pkn.provenance.json`. Existing outputs will be
overwritten only after the corrected tests pass; they are derived results,
not source data.

## Validation

1. Run the targeted PKN tests.
2. Run the full suite with `pytest -q tests`.
3. Run the complete stage command:

   ```bash
   bash scripts/setup/run_stage03.sh 2>&1 | tee logs/stage03.log
   ```

4. Verify the regenerated sidecar identifies the implementation commit and
   hashes the stage-01 feature tables, stage-01 RNA expression table, and all
   stage-02 regulon summaries actually read by the script.
5. Inspect the corrected COSMOS counts, node grammar, three-leg path table,
   null summary, and figures before stage 04.

## Implementation review (2026-10-09, before rerun)

A read-only check of the implementation against the real NetworkCommons file,
the HGNC table and the stage-01 expressed genes found three defects. All three
are fixed in `src/pkn.py`, with regression tests in `tests/test_pkn.py`:

1. **X-prefix guard.** The parser raised on every node starting with "X",
   which includes 20 plain HGNC symbols in the resource (XBP1, XIAP, XPO1,
   XRCC1, ...), so the COSMOS step would have stopped. Only the metaPKN
   grammar (`X<digits>`, `XGene<digits>`, `XMetab__`) is now rejected.
2. **Model metabolites with underscores.** IDs such as `Metab__gd1b2_hs_g`
   or `Metab__2hibup_S_r` (797 nodes, 4,713 unpruned edges) were classified
   `other` and all removed by P6, contrary to the decision to keep
   metabolites. Any `Metab__<id>_<compartment>` node is now a metabolite;
   only a whole-ID HMDB match (spaces allowed, e.g. `Metab__HMDB10384 _c`)
   is normalised to an HMDB ID.
3. **Log counts.** `P5_hgnc` type counts were summed over edge endpoints,
   not unique nodes. They are now unique-node counts for every type, and
   `P6_expression` reports the removed nodes per type (`other` is now the 325
   protein-complex or unparsed nodes, removed as for OmniPath).

Expected after the rerun (computed in memory with the fixed module):
COSMOS primary 54,571 edges / 26,510 nodes; unpruned 81,788 / 36,363;
no `Metab__` node classified `other`. Aim-5 shortest lengths in primary:
HMGCR→mevalonate_c 2; HMGCR→cholesterol_c 21; cholesterol_c→E2F1, SP1, JUN 2
and →EGR1, JUNB 3 (via MAPK1/MAPK3); mevalonate→TFs none; WT1 absent (not
expressed); SREBF1/2, SCAP, INSIG1→TFs none.

## Open questions

1. **Aim-5 metabolic leg readout.** The shortest HMGCR→cholesterol_c path is
   not the biosynthesis route. It runs HMG-CoA→mitochondrial ketogenesis
   (HMGCL, OXCT1)→succinate dehydrogenase→fatty-acid steps→palmitate→ABCA1
   protein→`Gene5715__ABCA1`→cholesterol, through cofactor and side
   metabolites. All 18 canonical enzymes (MVK to DHCR24) are present and
   expressed. `positive_route_present` therefore does not show that the
   mevalonate pathway is in the PKN. Options: (a) remove currency metabolites
   before the path search; (b) test the ordered canonical route explicitly;
   (c) report the current readout as reachability only. Undecided; this also
   limits how the leg's rewiring null can be interpreted.

Stage 04 CORNETO controls, lambda path, solver, time limit, and inference
validation remain separate decisions.
