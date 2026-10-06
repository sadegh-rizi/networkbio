# tests/

Automated tests for the reusable functions in `src/`.

Tests use **synthetic fixtures**, not the real data. Build a small table with
known values, assert the function does the specific thing it claims, and keep
the whole suite fast enough to run on every change.

What is worth testing here, in rough priority order:

1. Sample-to-line mapping: every column maps to a line and a condition, and an
   unmapped column raises.
2. Aggregation of replicates to lines: that it groups by the right key and that
   the effective n comes out as lines, not replicates.
3. Filtering and exclusion accounting: that the reported number of removed
   features matches what was removed.
4. Sign conventions and identifier mapping in the input tables for the network
   method.
5. A small network problem with a known answer (needs the solver; mark it so it
   can be skipped when the solver is absent).
6. Anything with a seed: the same seed gives the same answer.

Run from the repository root:

```bash
pytest
```
