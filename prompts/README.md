# prompts/

Saved prompts behind non-trivial AI-generated code, analysis, or plan
documents in this repository. See `skills-to-install/ai-provenance/SKILL.md`
for the scheme.

Convention: `YYYY-MM-DD-short-topic-implementation.md`, containing at
minimum:

- The task/prompt text actually used.
- Model name/version and any non-default settings.
- The commit hash(es) of the code it produced (add after committing).

This directory records events, not a script to replay; the same prompt run
again is not guaranteed to reproduce the same code.

## Write the record in the same session that generates the code

A session that produced code but no `prompts/` file is not finished. The
`Prompt:` and `Model:` commit trailers are unwritable without the file.

## File types

| Pattern | What it is |
|---|---|
| `YYYY-MM-DD-<topic>-implementation.md` | The provenance record. Required. |
| `session-<id>.md` | Raw agent transcript. Optional, large. |
| `UNRECORDED-<topic>.md` | A known, acknowledged gap. Rename once recovered. |

## Template

```markdown
# <Topic>

Date: YYYY-MM-DD
Model: <model-name-and-version>
Settings: <non-default settings, if any>
Plan: doc/decisions/<file>.md
Commit(s): <hash> (add after committing)
Transcript: prompts/session-<id>.md (or: not saved)

## Prompt actually used

<verbatim text>

## What it produced

<files created or changed>

## What was changed by hand afterwards

<edits the human made on top of the generated code>

## Outputs stamped

<figures / tables and where their commit hash is recorded>
```

## Current state

| Date | Topic | Record |
|---|---|---|
| 2026-09-29 | repository scaffold | `2026-09-29-repository-scaffold.md` |
| 2026-09-29 | data inventory and download scripts | `2026-09-29-data-inventory-download-scripts.md` |
| 2026-10-02 | course proposal and bibliography | `2026-10-02-course-proposal.md` |
| 2026-10-06 | venv, CORNETO toy problem, stage-00 inputs | `2026-10-06-setup-toy-inputs-implementation.md` |
| 2026-10-06 | stage-01 preprocessing and stage-02 TF activities | `2026-10-06-stage01-02-implementation.md` |
| 2026-10-09 | stage-01/02 revision plan (review fixes) | `2026-10-09-stage01-02-revision-plan.md` |
| 2026-10-09 | stage-01/02 revision implementation (review fixes) | `2026-10-09-stage01-02-revision-implementation.md` |
| 2026-10-09 | toy CORNETO network figure (branch `viz/toy-network`) | `2026-10-09-toy-network-figure.md` |
