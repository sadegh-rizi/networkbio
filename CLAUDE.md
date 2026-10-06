# Claude Code instructions: networkbio (CORNETO on PMS iNSC multi-omics)

Read these files before proposing or making changes:

1. **`doc/analysis-state.md`** first. It is the current state in one file: what is settled, what is open, what to do next.
2. `AGENTS.md`
3. `doc/agent-rules.md`
4. `doc/context/` for the topic you are touching (index in `doc/context/README.md`), including **`corrections-log.md`** (claims already retracted; do not reintroduce them)
5. `doc/datasets.md` when the task touches any dataset
6. The relevant `doc/decisions/*.md` plan when the task implements one

## Working mode

- The analysis is a Python script and notebook workflow in WSL2 (PyCharm with a WSL interpreter). The network-inference method is CORNETO; it is a first-party dependency of this project, not an external CLI tool.
- **Sadegh runs every pipeline command himself.** Write or edit the scripts, then give the exact commands, labelled as WSL or Windows, with `nohup … > logs/<name>.log 2>&1 &` for long runs and dependent steps chained by `&&`. Do not execute pipeline runs.
- Re-read a file from disk immediately before editing it.
- Every number stated as a finding must name the result file it came from. Say "reasoning" when it is not from a source that was actually read.
- When a run or decision changes the state, update `doc/analysis-state.md` (and the relevant `doc/context/` file) in the same change.
- Do not alter dependencies or environments without approval. Do not `pip install`.
- Make only the requested change. Do not modify raw input data, commit, push, or clean unrelated files unless explicitly asked.
- Do not write to `.claude/`. Skills live in `skills-to-install/` and are copied in by hand.
- For a non-trivial change, first summarize the proposed files, assumptions, and validation command(s). Ask when the analysis question is genuinely ambiguous.

## Required finish report

State:

1. files changed;
2. commands actually run and their real results;
3. whether notebooks were executed;
4. assumptions and remaining data, biological, or statistical risks.

Project skills (copy from `skills-to-install/` into `.claude/skills/`):

- `/careful-planning` before starting any implementation;
- `/network-inference` before building or interpreting a CORNETO / prior-knowledge network;
- `/statistics-review` before statistical inference, permutation nulls, or model selection;
- `/python-analysis` for reusable Python logic;
- `/publication-plot` for any figure;
- `/research-code-style` whenever writing or reviewing code;
- `/ai-provenance` before committing AI-generated code, analysis, or figures.
