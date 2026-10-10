# Aim 5: HMGCR → cholesterol route in the COSMOS PKN (cofactor removal and ordered-route test)

Date: 2026-10-10
Status: proposed

Resolves open question 1 of `2026-10-09-stage03-cosmos-correction-plan.md`
(metabolic-leg readout). Stage-03 outputs are not changed; this adds a new
output folder.

## Question

Does the COSMOS meta-PKN contain the cholesterol biosynthesis route from
HMGCR to cholesterol, and can a path search recover that route instead of a
shortcut? The stage-03 metabolic leg (`hmgcr_to_metabolite`) found a 21-step
positive route to cytosolic cholesterol, but it runs through ketone, TCA and
fatty-acid metabolism and an ABCA1 protein node, and 19–20 of 20 rewired
networks have a route as short. Reachability therefore says nothing about the
pathway.

## Data available

- COSMOS primary and unpruned PKNs, `results/ionescu_corneto/03_pkn/cosmos/`
  (commit `c7f3719`): 54,571 / 81,788 edges.
- Raw NetworkCommons file `data/resources/pkn/cosmos_meta_network.sif`.
- Recon3D (`data/resources/pkn/Recon3D.json`, cached in stage 03): 10,600
  reactions, metabolite annotations (BiGG, HMDB, KEGG, ChEBI) and gene rules
  (`<Entrez>_AT<n>`).
- HGNC table (`data/resources/hgnc/`) for symbol ↔ Entrez.
- Not available: atom mappings, so a reaction graph links compounds through
  any shared co-substrate.

### Facts established in a read-only check (2026-10-10)

These guide the design. The implementation must recompute them.

1. COSMOS `Gene<N>__<symbols>` nodes are Recon3D reaction number N (1-based
   order in `Recon3D.json`), e.g. Gene1600 = HMGCOAR, Gene4992 = HMR_1467,
   Gene9762 = r0575. `_reverse` follows Recon3D's written direction. Recon3D
   writes HMGCOAR as mevalonate → HMG-CoA, so the step's direction must be
   taken from stoichiometry, not from the suffix.
2. Of 5,938 Recon3D reactions with a gene rule, 1,616 (27%) have no COSMOS
   node. A stoichiometric coefficient ≥2 does not explain this (9% of absent
   and 10% of present reactions have one).
3. Squalene synthase (FDFT1) is two half-reactions in Recon3D:
   FPP → presqualene diphosphate (`r0170` ER, `HMR_1465` cytosol) and
   presqualene diphosphate → squalene (`r0575` ER, `HMR_1467` cytosol).
   COSMOS contains only the second. Nothing produces presqualene diphosphate
   (`Metab__HC01118_*`) from FPP.
4. Restricted to the canonical enzymes, their metabolites and transports
   between them, HMGCR reaches farnesyl diphosphate but not cholesterol;
   cholesterol is reached backwards only up to squalene. Adding the single
   ER half-reaction `r0170` gives a 40-step route to `Metab__HMDB0000067_r`
   (42 to `_c`) through MVK, PMVK, MVD, FDPS, FDFT1, SQLE, LSS, CYP51A1,
   TM7SF2, NSDHL, HSD17B7, DHCR24, EBP, SC5D and DHCR7.
5. Most classic currency metabolites are already absent from the COSMOS file
   (no nodes for H⁺, ADP, NADH, NADP(H), CoA, Pi, CO₂, acetyl-CoA). The
   remaining hubs include FAD/FADH₂ (degree 171), H₂O₂ (146), ATP (66),
   GTP, NH₄⁺, ubiquinone and S-adenosylhomocysteine. The highest-degree
   metabolite is cholesterol in the ER (180), so a degree cut-off without
   protection would delete the target.
6. Unweighted shortest paths stay non-canonical in every variant tried
   (21–26 steps): removing cofactors moves the detour to
   fumarate → arginine → AKT1 → PRKACA → ABCA1; a metabolism-only graph moves
   it to succinate → 2-oxoglutarate → sulfate → steroid sulfatase (STS).
   A degree-weighted search on the metabolism-only graph without cofactors
   returns the canonical route (16 canonical enzymes) only after the `r0170`
   repair.

## Unit of observation

PKN structure; no samples are used. Results are structural statements about
the resource, not biological findings.

## Method

1. **Canonical gene set.** Pre-registered list (textbook pathway; KEGG
   hsa00900 terpenoid backbone and hsa00100 steroid biosynthesis): HMGCR,
   MVK, PMVK, MVD, IDI1, IDI2, FDPS, FDFT1, SQLE, LSS, CYP51A1, TM7SF2, LBR,
   MSMO1, NSDHL, HSD17B7, EBP, SC5D, DHCR7, DHCR24. At run time, fetch the
   two KEGG gene lists (`https://rest.kegg.jp/link/hsa/hsa00900`, `.../hsa00100`),
   cache them with provenance, and fail if a listed gene is in neither.
2. **Canonical step table.** From Recon3D, take every reaction whose gene
   rule contains a canonical gene (Entrez via HGNC). Orient each reaction in
   the pathway direction: the substrate is the pathway intermediate that
   appears in the preceding step (start: HMG-CoA → mevalonate). Ignore
   co-substrates in the cofactor list (step 4). Keep both sterol branches
   (Bloch via desmosterol/DHCR24 last; Kandutsch–Russell via
   7-dehydrocholesterol/DHCR7 last) and all compartments. Output one row per
   reaction: step order, Recon3D ID and number, substrate, product,
   compartment, gene(s), branch.
3. **Gap audit.** For each step row, classify against COSMOS primary and
   unpruned: `present_forward` (substrate → Gene<N> → product edges exist in
   the pathway direction), `present_reverse_only`, `pruned_unexpressed`
   (in unpruned only), `absent_from_cosmos` (no Gene<N> node). Expected:
   `r0170` and `HMR_1465` absent. Also write the genome-wide count of absent
   gene-associated reactions (fact 2) by Recon3D subsystem.
4. **Cofactor list.** Pre-registered BiGG base IDs, mapped to COSMOS nodes
   through Recon3D HMDB annotations and model IDs, in every compartment:
   h, h2o, atp, adp, amp, gtp, gdp, utp, udp, ctp, cdp, nad, nadh, nadp,
   nadph, fad, fadh2, coa, pi, ppi, o2, co2, hco3, nh4, h2o2, na1, k, cl,
   q10, q10h2, amet, ahcys. Canonical intermediates are protected even if
   listed. Write the mapped node list with degrees, and the unmapped IDs.
   Sensitivity: additionally remove metabolites with metabolic degree in the
   top 1%, except protected nodes.
5. **Repaired PKN variant.** Add only the canonical steps classified
   `absent_from_cosmos`, as `Gene<N>__<symbols>` nodes with Recon3D's number
   and substrate → reaction → product edges in the pathway direction
   (expected: `Gene9692__FDFT1` for r0170, `Gene4406__FDFT1` for HMR_1465).
   No other reaction is added. Written as a separate PKN variant
   (`cosmos/variants/canonical_repair/`), never replacing the primary.
6. **Ordered-route test.** On primary, unpruned and repaired PKNs:
   (a) route-restricted reachability from HMGCR to each cholesterol node using
   only canonical enzyme nodes, their substrates/products and transport nodes
   whose input and output are canonical intermediates; (b) a waypoint test that
   reports, for each ordered intermediate (HMG-CoA, mevalonate,
   phosphomevalonate, diphosphomevalonate, IPP/DMAPP, GPP, FPP, presqualene
   diphosphate, squalene, 2,3-oxidosqualene, lanosterol, cholesterol), whether
   it is reachable from the previous one inside the restricted graph. The
   first failing waypoint is the gap.
7. **Path-search variants.** From HMGCR to each cholesterol node, on
   primary and repaired PKNs: (a) unweighted full graph (stage-03 readout),
   (b) metabolism-only graph (reaction/transport ↔ metabolite edges, and a
   gene only to its own reaction nodes), (c) (b) minus cofactors,
   (d) (c) with degree weights (cost of entering a metabolite = its degree in
   the metabolism-only graph; Croes et al. 2006, J Mol Biol), (e) (d) with the
   top-1% degree sensitivity. Readouts per path: steps, cost, number of
   canonical enzymes on it, fraction of ordered canonical steps traversed,
   and the first non-canonical node.
8. **Type-preserving null for `metabolite_to_tf`.** The stage-03 rewiring
   swaps edges regardless of node type. Add a swap that only exchanges
   targets between edges of the same (source type, target type) class, keep
   degree and sign, and rerun the null for the cholesterol → TF pairs with
   100 replicates (stage 03 used 20; smallest attainable p was 1/21).
   Seed base 20261010.

## Open questions

1. Repair scope: only canonical steps absent from COSMOS (recommended), or
   also report a variant adding all 1,616 absent reactions? Recommended: only
   canonical steps; the genome-wide gap is reported as a count.
2. Whether the degree-weighted search (7d) is reported as the primary
   readout of the metabolic leg, with 7a–c as context (recommended), or only
   the ordered-route test (6).

## Files to change

- `src/metabolic_route.py` (new): canonical step table, gap audit, cofactor
  mapping, repaired variant, restricted and weighted path search.
- `src/pkn.py`: type-preserving rewiring function (new function; existing
  `rewire_signed_graph` unchanged).
- `scripts/analysis/04b_aim5_metabolic_route.py` (new).
- `scripts/setup/run_aim5_route.sh` (new).
- `tests/test_metabolic_route.py` (new): synthetic two-branch pathway with a
  missing step, a cofactor hub and a shortcut; checks gap detection, repair,
  waypoint failure, that weighting prefers the pathway, and that type-
  preserving swaps keep class counts.
- `doc/analysis-state.md`, `doc/decisions/README.md`.

## Outputs

`results/ionescu_corneto/03_pkn/aim5_route/`:
`canonical_steps.tsv`, `gap_audit.tsv`, `recon3d_absent_by_subsystem.tsv`,
`cofactor_nodes.tsv`, `waypoints.tsv`, `path_variants.tsv`,
`path_variant_edges.tsv`, `null_metabolite_to_tf_typed.tsv`,
`figures/canonical_route.{pdf,png}` (pathway diagram with gap marked),
`aim5_route.provenance.json`. Repaired PKN:
`results/ionescu_corneto/03_pkn/cosmos/variants/canonical_repair/{pkn,nodes,filter_log}.tsv`.

## Validation

- `bash scripts/setup/run_aim5_route.sh 2>&1 | tee logs/aim5_route.log`,
  after `pytest tests/test_metabolic_route.py`.
- Works if: facts 1–4 and 6 are reproduced; the gap audit lists exactly the
  FPP → presqualene diphosphate step as absent among canonical steps (or
  lists more, which is then the result); the repaired variant differs from
  primary by only the added nodes and edges.
- Falsified (design is wrong) if the restricted graph reaches cholesterol
  without repair, or if the weighted search returns the canonical route on
  the unrepaired PKN; then the "missing step" explanation is wrong and the
  step table must be rechecked.
