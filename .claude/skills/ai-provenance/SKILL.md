---
name: ai-provenance
description: Use whenever AI-generated code, analysis, or figures are committed in this repository — records the prompt, model, and generation event behind a change so it can be traced later, since re-running the same prompt does not reproduce the same code.
---

# AI provenance

LLM code generation is non-deterministic even at fixed settings — running
the same prompt again does not reliably reproduce the same code. So
provenance here is not a recipe to replay; it is a record of an event: the
prompt, the exact code it produced, and the result, bound together by a
stable identifier. The git commit is that identifier.

## Before committing AI-generated code

1. Save the prompt to `prompts/`, not just chat history. Plain text or a
   short Markdown/YAML file, one per substantive generation (a new module, a
   non-trivial function, a generated analysis or plan). The point is that
   the exact text used is on disk and diffable, not that every prompt needs
   a bespoke format.
2. Reference the prompt file from the commit message: its path (and a
   version/hash if the prompt file has since changed), plus the model
   name/version actually used (e.g. `claude-sonnet-5`, `gpt-5.6`) and, if
   relevant, temperature/seed. `git blame` on any line should be able to
   answer "what prompt produced this" via the commit it belongs to.
3. For agent sessions (Claude Code, Claude Science, Codex) that did
   more than a trivial edit, export and keep the session transcript
   alongside the prompt file, or at minimum note where it can be recovered.
   For an agent, the provenance is the whole trace of what ran and what fed
   into it, not just the opening prompt — a wrong result can originate in
   planning or tool use, not only the final generation step.
4. Stamp derived outputs (figures, tables, exported results) with the
   commit hash that produced them — in the filename, a sidecar
   `.json`/`.txt` metadata file, or embedded plot metadata — plus the
   input-data identity (accession + file, or a hash), so a figure can be
   walked backward: figure -> commit -> prompt file -> exact prompt and model.
5. For notebook cells that came from AI generation, add a short markdown
   cell immediately above stating what was asked (or pointing at the
   `prompts/` file) and which model produced it — literate co-location is
   the tightest binding available for exploratory notebook work, per
   `doc/agent-rules.md`'s traceability requirement.

## Large intermediates need a sidecar, not git

Large results (pickles, parquet, solver outputs) in `results/` are gitignored, so git
cannot record their provenance. Write a `<name>.provenance.json` next to each
one — and commit that — recording: source dataset and accession, the code
commit that produced it, the parameters used, the seed, and the package
versions that matter (corneto, cvxpy, the solver, decoupler). A result file
with no sidecar is an orphan: nobody can tell later which QC thresholds are
baked into it.

## Minimal version, if adopting all of this at once is too much

A `prompts/` directory, one commit per AI-generated change whose message
references the prompt file, and outputs stamped with their commit hash.
That alone answers "which prompt produced which code and which result"
using git and nothing new.

## Pitfalls

- Do not treat the prompt text alone as sufficient provenance — without the
  actual generated code and the model/settings pinned alongside it, the
  prompt cannot reconstruct anything, since the same prompt on a different
  model version (or even the same model, rerun) can produce different code.
- Do not let prompts live only in chat/session history that isn't exported
  — if the session is gone, the prompt is gone.
- Do not skip this for "quick" AI-generated changes that touch QC
  thresholds, normalization, integration, statistics, or figures — those are
  exactly the changes `doc/agent-rules.md` already requires be traceable to
  their source and transformations; this skill is how that traceability
  actually gets recorded for the AI-generation step specifically.

## Completion criteria

For any non-trivial AI-generated change, you should be able to point to:
the prompt file in `prompts/`, the commit that applied the generated code
(with model/settings in its message), and — if the change produced a figure,
table, or large result file — the commit hash stamped on that output or its sidecar.
