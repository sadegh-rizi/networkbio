# How agents are used in this project

Carried over in shortened form from the melanoma thesis repository, where
the conventions were developed.

## The load-bearing idea: plans are the interface

`doc/decisions/*.md` is the handoff protocol. A plan written by one agent,
confirmed by the user, and implemented by another is a durable, reviewable
contract. It survives context compaction and switching between tools, and it
makes "did you do what we agreed" a checkable question. The matching
provenance convention (`prompts/*.md` plus commit trailers) closes the loop:
plan → prompt → commit → result.

## Roles

| Role | Tool | Notes |
|---|---|---|
| Design, literature and statistics review | Claude Science (Windows app) | Reads the repo, writes plans and code into granted folders. Its sandbox cannot run the WSL project environment. |
| Implementation and runs | Claude Code or another agent in WSL, and the user | The user runs every pipeline command. |
| Audit | a different model family from the implementer | Checks the implementation against the plan file, not the chat. |

## Deliberately not done

- No orchestrator agent that spawns other agents; the user is the orchestrator.
- No agent-to-agent messaging; plans and results go through files.
- No same-model self-review; an agent reviewing its own work reproduces its own blind spots.

## Folder access in Claude Science on Windows

The repository path contains `#`. Grants that worked in the sibling
repositories: a read-write grant on a repository with no `.claude` folder, or
a read-only root plus read-write subfolders when `.claude` exists. Agents do
not write into `.claude/`; skills are staged in `skills-to-install/` and
copied by the user. `git` does not run in the Claude Science sandbox, so the
user runs git in WSL.
