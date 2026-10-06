# networkbio

Course project for Network Biology (2 weeks): knowledge-driven network
inference with **CORNETO** (Rodriguez-Mier et al., *Nature Machine
Intelligence* 2025) on public multi-omics data from progressive multiple
sclerosis (PMS) induced neural stem cells (iNSCs).

The question is which signalling / regulatory subnetwork, drawn from a
prior-knowledge network, explains the PMS-vs-control differences that
Ionescu et al. (*Cell Stem Cell* 2024) attribute to an HMGCR-driven
cholesterol-synthesis state, and how that network compares with the
authors' own MOFA + GENIE3 analysis. See `doc/analysis-state.md` for the
current state and `doc/decisions/` for the plan.

## Repository contents

- `data/`: source datasets, one directory per study. Immutable, **not**
  committed. See `doc/datasets.md`.
- `notebooks/`: exploratory analyses.
- `src/`: reusable functions (I/O, preprocessing, activity estimation,
  network inference wrappers, plotting).
- `tests/`: automated tests on synthetic fixtures.
- `results/`: derived outputs in numbered pipeline stages.
- `doc/`: dataset notes, agent rules, analysis state, plans and decisions.
- `prompts/`: provenance records for AI-generated code and analysis.
- `skills-to-install/`: agent skills to copy into `.claude/skills/` by hand
  (agents do not write to `.claude/`).

## Data

Two papers from the same lab on the same iNSC model. Accessions were read
from the papers' data-availability sections; download status and file
contents are recorded in `doc/datasets.md`.

| Study | Modalities | Accessions |
|---|---|---|
| Ionescu et al. 2024, *Cell Stem Cell* | intracellular and secreted proteomics; intra- and extracellular metabolomics and lipidomics; 13C tracing | MassIVE MSV000095283; Metabolomics Workbench ST003331, ST003332, ST003328, ST003330 |
| Park et al. 2025, *Neuron* | bulk RNA-seq, scRNA-seq, snATAC-seq, WGBS | GEO GSE297192, GSE297365, GSE297690, GSE251839 |

The modalities are matched at the level of cell line, not of sample:
different experiments, passages and replicates. Effective n is 3 control
and 4 PMS lines (fewer for some layers). Results are hypothesis-generating.

## Setup

Commands below are for **WSL / Ubuntu**, run from the repository root.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

The environment is not created yet, and `requirements.txt` lists the
network-inference packages as commented-out lines to install deliberately.
See `AGENTS.md` (environment section) before installing anything.

## Usage

Run tests from the repository root with the environment active:

```bash
pytest
```

## Working with AI agents

Read `AGENTS.md` and `doc/agent-rules.md` first. Confirmed plans live in
`doc/decisions/`; provenance for AI-generated code lives in `prompts/`.
