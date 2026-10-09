# HGNC ambiguous Ensembl mappings

Date: 2026-10-09
Status: implemented; stages 01 and 02 completed 2026-10-09

## Question

The real HGNC 2026-10-06 complete set contains multiple distinct approved
symbols for some Ensembl IDs. How should stage 01 map those IDs without
duplicating their counts or choosing an arbitrary symbol?

## Confirmed decision

Exclude an Ensembl ID from the one-to-one mapping when it has more than one
distinct `Approved` HGNC symbol. Do not duplicate its count or select a symbol.
Keep the count row in `rnaseq/gene_map.tsv` with reason
`ambiguous_approved_mapping`, list it in `summary/exclusions.tsv`, and record
the HGNC ambiguity and excluded-row counts in `summary/counts.tsv`.

The 2026-10-06 release contains three such IDs; all three occur in the RNA-seq
count matrix. Two have nonzero counts and are therefore excluded from the
symbol-level expression matrix rather than being assigned an arbitrary symbol.

## Files and method

- `src/preprocess.py`: identify IDs with more than one distinct approved
  symbol, remove only those IDs from the mapping table, and add explicit
  status, reason, and QC metrics.
- `scripts/analysis/02_preprocess.py`: record the rule in stage provenance
  parameters; existing exclusion and count-table orchestration is retained.
- `tests/test_preprocess.py`: test that ambiguous rows are excluded and
  reported without changing the existing CPM denominator behavior.
- `doc/datasets.md` and this decision record: document the handling.

## Validation

- Synthetic mapping regression test.
- `".venv/bin/pytest" -q tests/test_preprocess.py`.
- `".venv/bin/pytest" -q -k 'not stages_end_to_end_on_synthetic_tables'`.
- `".venv/bin/python" -m compileall -q src scripts/analysis tests`.
- Small chunked inspection of the real HGNC and RNA-seq tables.
- Real stage-01 run completed after the failed output was moved aside.

## Open questions

None for this mapping decision.
