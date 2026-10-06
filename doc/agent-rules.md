# Shared agent rules: networkbio

These rules apply to every coding agent. They are safeguards for reproducible exploratory research, not a substitute for the supervisor's or the user's statistical and biological decisions.

## Scope and data integrity

- `data/` contains source research artifacts from public repositories. Never modify them in place, and never commit them.
- Inspect table orientation, identifiers, sample columns, and missing/zero encodings before writing preprocessing code. Do not assume rows are features or columns are samples without checking the file.
- Maintain a traceable link from each derived table or figure to its source dataset, accession, sample subset, filtering rules, transformations, and code version.
- Use project-relative paths.
- Record the accession, publication, and download date for every dataset in `doc/datasets.md`. A processed table from a repository has had undocumented decisions applied to it: say which are known and which are not.
- The data derive from patient fibroblasts. Do not attempt re-identification.

## Study-design safeguards

**The unit of observation is the cell line (donor).** The Ionescu cohort is 3 control and 4 PMS lines. Each line was measured in several replicates (mostly at least 3, fewer for bulk RNA-seq), and replicates are technical or culture replicates of one line, not independent samples. Pool or model them so that the effective n is stated as lines. The scRNA-seq and snATAC-seq data are cells nested in lines: thousands of cells from one line are one observation of that line.

**The layers are matched at the level of line, not of sample.** Proteomics, metabolomics, lipidomics and RNA-seq come from different experiments, passages and replicates, and (for the Park data) a different paper. Do not describe any joint analysis as if all layers were measured on the same material. Whether Park and Ionescu used the same donors is inferred from the cohort description and has to be checked against both Table S1 files before it is stated as fact.

**Disease status is confounded with age, sex and genetic background.** The authors state this limitation themselves. With 7 lines it cannot be adjusted for. Report it wherever a PMS-vs-control effect is discussed.

**Treatment coverage differs by layer.** Simvastatin (SV) data are confirmed for lipidomics and secretomics only. Do not assume an SV arm exists for any other layer without checking the deposit.

Before any inferential analysis, establish and state:

1. the biological unit of observation and the experimental unit, which differ here;
2. the sample-to-line mapping, condition, treatment, replicate, and batch fields for each layer;
3. whether observations are independent, paired, nested, or repeated (control vs PMS lines are independent; SV vs vehicle within a line is paired);
4. the intended contrast, covariates, exclusion rules, and missing-data handling.

## Preprocessing and quality control

For each layer record what is known versus unknown about the deposited table: raw or processed, normalisation, filtering, imputation, identifier system, and software version. For mass-spectrometry layers in particular, missing values are informative (detection limits). The secretome table has many missing values; the authors filtered proteins expressed in fewer than 20% of lines and imputed the rest. Do not repeat or replace an imputation silently.

Propose filtering, normalisation, transformation, and imputation separately. For each, give the rationale, the affected features, and how it will be recorded. Keep an exclusion table with counts and reasons.

## Activity estimation

- Transcription-factor (and, if used, kinase or pathway) activities are model outputs, not measurements. Record the resource (for example a TF-target collection), its version and date, the statistic, the minimum targets per TF, and how identifiers were mapped.
- With 2 control and 2 PMS lines in the Ionescu bulk RNA-seq, an activity contrast has very little power. A per-line activity table is preferable to a single contrast statistic if the network method can use it. State which is used and why.
- Signs matter to sign-consistent network methods. Record the sign convention and check it on a known case.

## Network inference (prior knowledge + optimisation)

- A prior-knowledge network (PKN) is an input, not a result. The inferred network can only contain edges that the PKN contains. Record the PKN source, version or download date, the filters applied (organism, signed, directed, confidence), and the node and edge counts before and after each filter.
- The inferred subnetwork depends on the chosen input nodes, the measured nodes and their values, the sparsity or regularisation parameters, and the solver settings. Record all of them with each result. A network without its settings is not a result.
- Solutions are generally not unique. Report the solver, its version, the seed, the time limit and optimality gap, and whether the run finished at optimality. Sample multiple solutions, or vary the seed, and report edge frequencies rather than one solution as if it were the answer.
- Select regularisation parameters with a criterion stated **before** looking at which edges appear, and do not tune them until a preferred edge appears.
- **Validation must not be circular.** If the known mechanism (for example HMGCR to SASP-associated TFs) is supplied as an input, then recovering it is not evidence. Fix the positive and negative controls in the plan before running: which nodes are inputs, which are measured, which paths are expected, and which are not.
- Compare against a null where feasible. With 3 control and 4 PMS lines there are only 35 distinct ways to relabel lines into groups of 3 and 4, so a label-permutation p-value cannot be smaller than 1/35, and adjacent permutations are strongly correlated. State this limit; use node/edge-level randomisation of the input values or PKN as a complementary null.
- Do not infer causality from an inferred network. It is a hypothesis about which prior-knowledge edges are sufficient to explain the input pattern.
- Compare with the authors' MOFA + GENIE3 results as a baseline, and say what each can and cannot claim.
- If a controllability (driver-node) analysis is added, it inherits every bias of the network it is run on, and its minimum driver sets are usually not unique. Report it as one representative solution or as frequencies.

## Statistical analysis

- Write the question, outcome, predictor/contrast, covariates, unit of analysis, and model/test before running it.
- Distinguish exploratory work from confirmatory claims. At n = 3 vs 4 lines, nearly everything is exploratory.
- Report effect size and uncertainty alongside p-values. Define the multiple-testing family before testing and use an explicit correction (normally FDR for high-dimensional features). Report raw and adjusted p-values and the correction universe.
- Do not try multiple tests, transformations, subgroups, or cutoffs until a significant result appears.
- For enrichment, record the identifier mapping, universe/background, method, database version and date, and how unmapped identifiers were handled. The background is what was measured in that layer, not the whole genome.
- For MOFA, GENIE3, WNN or any factor or covariation analysis: document tuning and selection, assess stability where feasible, and do not infer causality from loadings or correlations.

## Figures and tables

- Use a colourblind-safe palette and do not encode groups by colour alone when a second cue is practical.
- Label axes with units and transformations, state sample size as lines (not only replicates or cells), and include a legend.
- Show individual observations when sample sizes permit; do not default to bar charts of mean ± SEM.
- Networks: state the settings that produced them (PKN, λ, solver, seed count), whether edges are single-solution or frequency-weighted, and what node colours and edge signs encode.
- Save every figure as vector PDF plus a presentation-resolution PNG and a manuscript-resolution PNG through `plot_config.save_figure` (`src/plot_config.py`).
- Separate exploratory figures from final inferential figures. Titles and captions must not make a stronger biological claim than the data support.

## Results directory organisation

- Never write results into one flat directory. Group `results/<analysis>/` into numbered stage subfolders: `00_inputs/`, `01_preprocessing/`, `02_activities/`, `03_pkn/`, `04_inference/`, `05_robustness/`, `06_validation/`, `07_downstream/`. Tables sit at the stage level; figures go in that stage's `figures/` subfolder.
- When outputs vary along an experimental axis (PKN version, λ, solver, condition set, layer subset), put that axis in the **directory path and keep filenames identical across variants**.
- Each of the three `save_figure` outputs belongs in the same folder.
- Large intermediates live in `results/` but are gitignored. Their provenance is recorded in a committed sidecar (`<name>.provenance.json`) with accession, code commit, parameters, seed, and package versions.

## Python and notebook quality

- Use `pathlib.Path`, explicit function inputs/outputs, meaningful names, and type hints on reusable public functions.
- Use deterministic seeds for stochastic steps and record them. Solver behaviour can still differ across versions and thread counts; pin versions where a result matters.
- Keep notebook cells focused on orchestration and inspection; move logic into `src/` once it is reused.
- Preserve notebook outputs unless asked to clear them.
- Default to the simplest correct implementation: plain functions over classes, few parameters over speculative flexibility. See `.claude/skills/research-code-style/SKILL.md`.

## Completion standard

Before handoff, report modified files, real verification commands and results, whether data or notebooks were executed, the analytical assumptions, and residual risks. A passing script is not evidence that a biological or statistical conclusion is valid.
