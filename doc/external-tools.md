# External tools

Every third-party package or resource that produces a result in this
analysis, with the exact version. A tool that is not in this table is not
reproducible. This is the same job `prompts/` does for AI-generated code:
bind an output to the exact thing that produced it.

## Rules

- **Pin exact versions**, and for resources (regulons, prior-knowledge networks) record the download date and source URL as well. Databases change without a version number.
- Packages that are the method of this project (CORNETO and its solver stack) go in the project venv and `requirements.txt`. Packages that carry conflicting dependencies or are run as a CLI can be isolated instead; ask first.
- Document the handoff schema for every tool boundary: what goes in (orientation, units, sign convention) and what comes out.
- Record why the tool was chosen.

## Ledger

| Tool or resource | Source | Version / commit / date | Used for | Handoff in → out | Added | Why |
|---|---|---|---|---|---|---|
| CORNETO | PyPI `corneto[os]` | 1.0.0rc8 (pre-release, 2026-10-02) | multi-sample network inference (`CarnivalFlow`) | `cn.Data` (per sample: node → value, role input/output) + signed graph → edge values per sample | 2026-10-06 | method chosen by the user; toy problem passed |
| HiGHS (highspy) | PyPI | 1.15.1 | MILP solver for CORNETO via cvxpy 1.9.3 | | 2026-10-06 | open source; used in the CORNETO tutorial |
| SCIP (pyscipopt) | PyPI | 6.2.1 | alternative MILP solver | | 2026-10-06 | comes with `corneto[os]`; for cross-checks |
| decoupler | PyPI | 2.2.0 | TF activity from expression (ULM) | lines x gene-centred log2(CPM+1), or moderated-t/logFC/Welch vectors → TF x line/contrast ULM scores, raw and BH-adjusted p-values | 2026-10-06 (synthetic ULM verified; real run pending) | TF inputs for signalling inference |
| omnipath (client) | PyPI | 1.0.12 | PKN and regulon download (CollecTRI, OmniPath) | | 2026-10-06 (installed, not yet used) | resource access |
| networkcommons | PyPI | not installed | would provide `get_cosmos_pkn` | | rejected 2026-10-06 | pins graphviz<0.18, incompatible with corneto 1.0.0rc8 |
| Prior-knowledge network | _to decide_ | _download date_ | network for inference | | planned | |
| HGNC complete set | [quarterly archive](https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2026-10-06.txt) | release 2026-10-06; downloaded 2026-10-09; fallback 2026-07-07 only if absent | approved gene-symbol mapping; ambiguous IDs excluded and recorded | Ensembl count rows → summed counts per approved symbol | 2026-10-09 | dated, reproducible identifier mapping |
| CollecTRI | [Zenodo 8192729](https://zenodo.org/records/8192729/files/CollecTRI_regulons.csv?download=1), via `decoupler.op.collectri` | human, complexes retained (`remove_complexes=False`); download pending | primary TF regulon | signed source/target/weight edges → ULM, >=5 observed targets | 2026-10-06 | primary resource in confirmed stage-02 plan |
| CollecTRI no-complex | same Zenodo source, via `decoupler.op.collectri` | human, AP1/NFKB removed (`remove_complexes=True`); download pending | complex-retention sensitivity | same signed ULM handoff as CollecTRI | 2026-10-09 | registered before inspecting held-out TF activities |
| DoRothEA | [OmniPath interactions](https://omnipathdb.org/interactions), via `decoupler.op.dorothea` | human, A-C; default signed confidence weights (divisors 1/2/3); download pending | regulon sensitivity | same gene universe and ULM settings as CollecTRI | 2026-10-06 | resource sensitivity in proposal Aim 1 |

Download scripts/functions write resource sidecars with date, row count and
sha256. Each stage copies those records into its provenance alongside exact
package versions. In decoupler 2.2.0, ULM returns BH-adjusted p-values;
`src/activities.py` explicitly reconstructs the raw t-distribution p-values.
Correction families are all retained TFs within each resource and each line
or contrast, with no significance-based selection. RNA-seq size-factor
normalisation and the moderated-t contrast are implemented in numpy/scipy;
they do not add a package dependency.

## Checking a tool before trusting it

1. It runs on a small synthetic input where you know the right answer.
2. Its version is the one the paper or tutorial used, if you are reproducing one.
3. Its assumptions about your input match reality: sign convention, scale, orientation, identifier type.
4. It is deterministic, or its seed and solver settings are recorded.

A tool that produces plausible output on real data has not been validated.
