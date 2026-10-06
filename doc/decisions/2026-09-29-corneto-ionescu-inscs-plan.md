# CORNETO network inference on PMS iNSC multi-omics

Date: 2026-09-29
Status: proposed

Nothing in this document is confirmed. The method steps are a candidate
design written from the discussion so far; the open questions below decide
what is actually implemented.

## Question

Which subnetwork of a prior-knowledge network (PKN) does CORNETO select to
explain the differences between PMS and control iNSC lines? Does it recover
the route proposed by Ionescu et al. 2024 (an HMGCR-driven lipogenic state,
cholesterol-dependent TFs, and the SASP) when that route is not supplied as
an input? How does the result compare with the authors' MOFA + GENIE3
analysis?

Why now: two-week Network Biology course project. The user chose CORNETO over
the multilayer-spectral, hypergraph and controllability options.

## Data available

See `doc/datasets.md`. Candidate layers:

- Transcriptomics for TF activities: Park et al. bulk RNA-seq (GSE297192) or
  pseudobulk of scRNA-seq (GSE297365). The Ionescu bulk RNA-seq (2 vs 2 lines)
  is not listed as deposited.
- Metabolomics (intra- and extracellular), lipidomics, and proteomics
  including secretomics: Ionescu deposits (MassIVE and Metabolomics
  Workbench).
- Simvastatin arm: confirmed only for lipidomics and secretomics.

Not available: paired measurements on the same sample across layers; an
isogenic control; phospho-proteomics (the deposited layers are not
phospho-specific as far as read).

## Unit of observation

The cell line. 3 control and 4 PMS lines at most, fewer for some layers.
Replicates within a line are averaged or modelled, not counted. Control vs
PMS is between-line (independent); SV vs vehicle is within-line (paired).
Age, sex and genetic background are confounded with disease status.

## Method (candidate)

1. Toy problem: run a small CORNETO example with a known answer and confirm
   the solver works in WSL. Record versions in `doc/external-tools.md`.
2. Build per-line input tables from the chosen layers (identifiers mapped,
   sign convention documented).
3. Estimate TF activities per line with a stated resource and statistic.
4. Build the PKN (source, date, filters recorded; node and edge counts
   logged).
5. Define input nodes, measured nodes, and the positive and negative controls
   **before** any inference (see open question 4).
6. Run single-condition inference, then CORNETO's multi-sample inference across
   the chosen conditions.
7. Select the regularisation parameter with a pre-stated criterion; sample
   multiple solutions; report edge frequencies.
8. Compare with the authors' MOFA/GENIE3 results and with the controls.
9. Optional: driver-node (maximum matching) analysis on the inferred network.

## Open questions

Each one has to be answered by the user (or supervisor) before this plan can
become `confirmed`.

1. **Conditions.** Control vs PMS only, or also PMS + SV where the layer has
   it?
2. **Layers.** Which of transcriptomics, metabolomics, lipidomics, and
   proteomics or secretomics enter as measured nodes?
3. **Signalling only, or signalling plus metabolism?** Whether CORNETO can
   handle one combined problem, or only separate ones, is not yet verified
   against its documentation.
4. **Input nodes and circularity.** Which nodes are inputs? If HMGCR or the
   SREBP axis is given as an input, recovering the paper's route is not
   evidence. Options: unbiased inputs (for example receptors or all
   candidate upstream nodes) with the known route held out as the test.
5. **Transcriptomics source.** Park bulk (line count unverified) or
   pseudobulk of the scRNA-seq; and whether Park and Ionescu are the same
   donors (Table S1 comparison).
6. **PKN.** Source, version and filters (signed and directed edges only?).
7. **Regularisation parameter.** Selection criterion and grid.
8. **Solver.** Which MILP solver, and the time limit and gap.
9. **Replicates.** Average per line before inference, or keep replicates as
   separate samples in the multi-sample problem?
10. **Optional controllability step.** In or out of scope for two weeks?
11. **Course deliverable.** Format, deadline and grading criteria are unknown.
12. **Where CORNETO lives.** Project venv (needs approval) as proposed in
    `AGENTS.md`, or an isolated environment.

## Answers and proposed defaults (2026-10-06)

Status stays `proposed`. Answers from the user are marked (user); the rest
are proposed defaults to confirm. Layer and tool choices follow the course
proposal (`doc/proposal/proposal.tex`, 2026-10-02).

- Q1 Conditions (user): simvastatin (SV) is a separate condition. Only the
  lipidomics (ST003328) has an SV arm (3 Ctrl + 4 PMS lines x untreated/SV),
  so SV samples carry lipid measurements only. In SV samples HMGCR = -1 is a
  known input; in the baseline contrast HMGCR must not be an input.
- Q2 Layers (proposal): GSE297192 bulk RNA-seq (TF activities), ST003331
  intracellular and ST003332 extracellular metabolomics, ST003328 lipidomics.
  Proteomics is private on MassIVE and is out.
- Q3 Signalling + metabolism (checked): CORNETO has no COSMOS-specific
  method, but the COSMOS meta-PKN encodes metabolism as signed, directed
  edges, so `CarnivalFlow` on it is the joint problem. The toy problem
  (`scripts/analysis/00_toy_corneto.py`) tests this with a metabolite node.
  COSMOS runs signalling->metabolism and metabolism->signalling separately;
  whether to do the same is still open.
- Q4 Inputs and controls (proposed): baseline inputs are either all source
  nodes with unknown sign (value 0) or receptors inferred from RNA-seq; never
  HMGCR, SREBF1/2 or the paper's TFs. Positive control: signed path from
  cholesterol synthesis to the paper's TFs exists in the PKN, and how often
  those nodes appear across sampled solutions. Negative controls:
  degree-preserving rewired PKN, the 35 label permutations, a planted route;
  an "expected absent" pathway still to be chosen from the paper.
- Q5 Transcriptomics (proposal): Park bulk C1-C6 (3 vs 3; split from the
  authors' code). Donor identity across papers still unverified.
- Q6 PKN (proposed): OmniPath (SIGNOR first), human, consensus direction,
  signed (stimulation XOR inhibition), then the COSMOS meta-PKN for the joint
  problem; prune to expressed genes and to nodes on input->measurement paths;
  log node/edge counts after each filter.
- Q7 lambda (open): grid and selection rule to be fixed before inference.
- Q8 Solver (proposal): HiGHS (open source); time limit and gap open.
- Q9 Replicates (open).
- Q10 Controllability: out (not in the proposal).
- Q11 Deliverable (open).
- Q12 CORNETO location (user, 2026-10-06): project venv; see
  `2026-10-06-setup-toy-inputs-plan.md`.

## Files to change

_To be filled once the questions are answered._ Expected: `src/` modules for
input building, activity estimation and inference wrappers; `notebooks/`;
`requirements.txt` and `doc/external-tools.md`; `doc/datasets.md` after
download.

## Outputs

_To be filled._ Under `results/ionescu_corneto/` in the numbered stages of
`doc/agent-rules.md`.

## Validation

- The toy problem returns its known network.
- Input tables: every column maps to a line and condition; counts per line
  per layer are logged.
- Positive control: the expected route is recovered without being given.
  Negative control: an expectation that must not appear (to be defined in
  question 4). What would falsify the approach: the route is only found when
  its nodes are inputs, or the result changes qualitatively across
  reasonable parameter and seed choices.
