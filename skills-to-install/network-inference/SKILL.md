---
name: network-inference
description: Use before building, running, or interpreting a prior-knowledge network inference (CORNETO / CARNIVAL-style) in the networkbio project: choosing inputs and measured nodes, the PKN, sign conventions, sparsity, solver settings, validation, and reporting.
---

# Network inference with a prior-knowledge network

Read `doc/agent-rules.md` ("Network inference" section) and `doc/context/corneto-notes.md` first. Most API statements in the notes are unverified; check the pinned version before relying on them.

## What the method does, in one paragraph

Given a prior-knowledge network (signed, directed edges), a set of input nodes, and measured values on other nodes (for example TF activities), the method selects a sparse subnetwork whose sign-consistent paths explain the measured pattern. In the multi-sample case, the samples are solved jointly with a penalty that encourages sharing edges across samples. The output is a subset of PKN edges and node states per sample. It contains no edge that the PKN does not contain.

## Checklist before a run

1. **Toy problem first.** A small network with a known answer, run end to end on this machine, with the solver named.
2. **Inputs and measurements written down.** Which nodes are inputs, which are measured, what value and sign each gets, and where each value came from (file, column, statistic).
3. **Circularity check.** Is the result you plan to call a validation supplied, directly or indirectly, as an input? If yes, it is not a validation. Fix the positive and negative controls in the plan.
4. **PKN provenance.** Source, download date, filters, node and edge counts, and how many measured nodes are actually in the PKN. Measured nodes missing from the PKN are unexplained by construction; report how many.
5. **Sparsity parameter.** Chosen by a stated criterion before looking at preferred edges. Report the whole path, not one value.
6. **Solver settings.** Solver and version, time limit, optimality gap, threads, seed. Note whether each run reached optimality.
7. **Multiple solutions.** Sample solutions or seeds; report edge frequencies.

## Interpreting

- The network is a hypothesis about sufficient prior-knowledge edges, not evidence of causation.
- Compare with the baseline (the authors' MOFA + GENIE3) and say what each can claim.
- Inspect what drives the result: which measured nodes contribute most, and whether removing one changes the network.
- A hub or a long path is often a property of the PKN's connectivity. Compare with a randomised PKN of the same degree distribution.
- Report failures too: runs that hit the time limit, measured nodes that could not be reached.

## Optional: driver nodes

Maximum-matching driver-node analysis on an inferred network inherits the network's biases, and minimum driver sets are usually not unique. Report one representative set together with how often each node appears across sampled matchings, and say that direction and sign errors in the input network propagate into it.
