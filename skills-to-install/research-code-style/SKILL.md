---
name: research-code-style
description: Use whenever writing or reviewing Python code in this project — keeps generated code at the complexity level of early-stage scientific/exploratory research rather than defaulting to production-software patterns, so a domain scientist can read, run, and modify it without a software-engineering background.
---

# Research code style

LLMs default to production-software instincts — classes, config objects,
broad abstraction, defensive error handling, CLI wrappers — because that's
the majority of their training distribution. Early-stage research code has
a different goal: a person who knows the biology/statistics but not
necessarily software architecture needs to read it top to bottom, see the
actual computation, and change one thing without breaking five others.
Default to the simplest thing that is still correct and readable; add
structure only when the code demonstrably needs it, not because it's good
practice in general.

## Prefer, by default

- Plain functions over classes. A class needs an actual reason (real
  internal state that changes across multiple method calls, or several
  methods sharing the same construction). "One function wrapped in a class
  with an `__init__` and a single `run()` method" is not a reason.
- A script/notebook cell that reads in the order a person would explain the
  analysis, over a "framework" with config objects, factories, or
  plugin-style dispatch built for hypothetical future variants that don't
  exist yet.
- Explicit, inlined logic over a deep chain of tiny wrapper functions that
  each just call the next — if the actual computation can't be seen without
  following several files, it's too indirect for this stage. This applies
  hard to analysis pipelines: a function that just calls a library routine
  with the same arguments every time is not a useful abstraction.
- Failing loudly (let the exception propagate) over broad `try/except`
  that swallows or logs and continues — in exploratory research code, a
  visible crash on bad input is more useful than a silently degraded
  result. Matches `doc/agent-rules.md`'s "avoid... bare except."
- Few, meaningful parameters over many optional kwargs "for flexibility" —
  add a parameter when a second real use case needs it, not speculatively.
- Docstrings/comments that explain *why* a domain choice was made (why this
  QC threshold, why pseudobulk, why this resolution) over exhaustively
  documenting *what* the code already obviously shows.
- One module sized to a single readable analysis step, over splitting into
  many small files before there is a real reuse reason.

## Still required, non-negotiable

Standard Python hygiene doesn't go away: `pathlib.Path` over string paths,
type hints on public function signatures, meaningful names, no mutable
default arguments, no silent type coercion, deterministic seeds for
anything stochastic. Simple and sloppy are not the same thing — the goal is
code a domain scientist can trust and modify, which still has to be correct
and unambiguous.

Two additions specific to single-cell work:

- Do not mutate an `AnnData` in place in a way that makes a cell
  non-idempotent. Write to a named layer or a new object; re-running a cell
  should give the same answer as running it once.
- Anything with a random seed — solver runs, permutation tests, sampling,
  permutation tests — takes the seed as an explicit argument. Never rely on
  a library's global default.

## When to add structure

Add a class, a config object, an abstraction layer, or split a module only
once one of these has actually happened, not preemptively:

- the same logic is now called from 3+ places with variations a parameter
  can't cleanly express;
- the same piece of state genuinely needs to persist and change across
  several calls;
- a second person or tool is about to reuse this code in a context you can
  name today, not "might."

## Pitfalls

- Do not add logging frameworks, CLI argument parsers, or settings files
  for a script with one caller (you, running a notebook cell) — a plain
  function argument or a notebook variable is the CLI.
- Do not build a "pipeline class" that wraps a linear sequence of library
  calls. The linear sequence *is* the readable form; wrapping it hides the
  analysis and makes changing one step harder, not easier.
- Do not vectorize/parallelize/optimize before it's needed — but note the
  inverse also applies here: a Python loop over 300,000 cells is not
  "simple," it's broken. Use the array operation when the data size demands
  it, and say so in a comment.
- Do not build a generalized/reusable utility for a hypothetical future
  analysis — build the specific thing this analysis needs; generalize
  later if it's actually reused.
- Do not let structure hide the scientific logic — a reviewer with domain
  knowledge but not a software background should be able to read the
  module top-to-bottom and see the same steps as the relevant plan doc, in
  the same order.

## Completion criteria

Before treating analysis code as done, check: could someone who knows the
biology/statistics but not software architecture read this in one pass and
see what it computes? Could they change one parameter (a threshold, a test,
a resolution, a column name) without touching unrelated code? If either
answer is no, simplify rather than document around it.
