# doc/data-inventory/

Inventory of the public deposits for Ionescu et al. 2024 (*Cell Stem Cell*) and
Park et al. 2025 (*Neuron*), checked against the repositories on 2026-09-29
(GEO, NCBI SRA, ENA, Metabolomics Workbench, MassIVE, Zenodo). AI-generated;
see `prompts/2026-09-29-data-inventory-download-scripts.md`.

| File | Content |
|---|---|
| `data_inventory_ionescu_park.xlsx` | Sheets: README, datasets (one row per deposit), samples (one row per GEO GSM / MW sample), files (URL + size), gaps (described but not deposited), reused_public (post-mortem datasets Park re-analysed) |
| `samples_ionescu_park.csv` | Same as the `samples` sheet (162 rows: 51 GEO samples, 111 Metabolomics Workbench samples) |
| `files_ionescu_park.csv` | Same as the `files` sheet |
| `experimental_design_ionescu_park.pptx` | 9 slides: shared iNSC model, design of each paper, deposited samples, availability, gaps, download plan |

Download scripts that use these accessions: `scripts/download/`.

Group labels come from the repository metadata, except where the
`group_source` column says otherwise (bulk RNA-seq C1-C6 split from the
authors' code; WGBS groups for lines 105/HD inferred from GSE297690).
Line IDs are not comparable across deposits; see the `gaps` sheet.
