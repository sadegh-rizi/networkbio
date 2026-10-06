# doc/context/

Established knowledge about this project that is not obvious from the code:
measured facts, where the labels come from, why decisions were made, and
what went wrong. It exists so the project's memory lives in git rather than
in one tool's chat history.

Start with `doc/analysis-state.md`. These files are the detail behind it.

| File | Topic |
|---|---|
| `corneto-notes.md` | what was read about CORNETO, with sources, and what is still unverified |
| `compute-environment.md` | WSL, the Windows sandbox, and solver notes |
| `corrections-log.md` | retracted claims and fixed bugs, so they are not reintroduced |

Conventions: each file has a "Last verified" date and names the result
files or sources its numbers came from. When a number changes, update it here
in the same commit as the rerun. Anything computed interactively and not
reproducible from `src/` is flagged as such.
