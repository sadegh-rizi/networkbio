---
name: careful-planning
description: Use before starting any implementation in the networkbio project (notebook cells, modules, scripts, or dependency changes) to resolve open design questions and get an explicit plan confirmed before code is written.
---

# Careful planning

This is the first gate for any non-trivial implementation task in this repository. Read `AGENTS.md` and `doc/agent-rules.md` first; this skill governs how you approach them, not what they say.

## Before writing any code

1. Restate the task and the intended outcome in your own words.
2. Identify every design, scientific, or technical decision the task implies that the user has not already specified. Check in particular:
   - which layers and conditions (control, PMS, simvastatin) and which accession, and whether a layer exists for the condition you need (`doc/datasets.md`);
   - the biological unit (cell line) versus replicates and cells, and the dependency structure (independent between lines, paired within a line);
   - for network inference: the PKN and its filters, input and measured nodes, sign convention, sparsity parameter and how it is chosen, solver, time limit, seeds, and how the result will be validated without circularity;
   - which method, test, or model, and why versus the alternatives;
   - normalisation, imputation, missing-value handling, covariates, contrasts, exclusions, and multiple-testing correction;
   - whether the change belongs in a notebook cell or a module;
   - output location and naming for any derived file, table, or figure;
   - whether any new dependency is required and whether it is installed;
   - compute cost: expected peak RAM and runtime (MILP problems can be slow), and whether an intermediate needs caching;
   - scope: what files will be touched, and what is explicitly out of scope.
3. Ask the user every open question directly. Do not silently pick a default for an ambiguous biological, statistical, or scope decision.
4. Once questions are answered, write the plan to `doc/decisions/YYYY-MM-DD-<topic>-plan.md`: the question, the data actually available, the method, the files to change, the outputs, the validation command(s), and an explicit list of anything still open.
5. Wait for explicit confirmation of the plan before making any edit or running any code. A plan you posted is not confirmation; the user has to say go.

## After confirmation

Proceed to implementation and invoke the relevant domain skill (`python-analysis`, `network-inference`, `statistics-review`, `publication-plot`).

## Pitfalls

- A plausible-sounding default is not a resolved question. Tutorial defaults for a sparsity parameter or a solver time limit are exactly this trap.
- Do not bundle open questions into one broad "let me know if this looks wrong"; ask a specific list.
- Do not skip this gate for "quick" changes that touch inference settings, preprocessing, or data handling.
- Do not plan an analysis whose runtime or memory you have not estimated.

## Completion criteria

Before implementing, you should be able to point to the open questions you surfaced, the user's answers, and the confirmed plan document in `doc/decisions/`. If any is missing, the gate has not been satisfied.
