# Commit Aim-5 and Stage-05 changes

Date: 2026-10-12
Status: confirmed by user

## Question

How should the current mixed worktree be validated and committed without
including failed or OOM execution artifacts?

## Confirmed scope

- Keep the Stage-05 robustness addendum as written, including its stated
  post-hoc status, process-isolated workers, sampling limits, and controls A-D.
- Commit the Aim-5 identifier-mapping correction, tests, audit documentation,
  regenerated route outputs, and provenance.
- Commit the successful Stage-04 outputs and the Stage-04/05 worker, runner,
  activity, inference, report, test, and documentation changes.
- Do not stage or delete
  `results/ionescu_corneto/04_inference_failed_2026-10-11_gc_import/` or
  `results/ionescu_corneto/04_inference_killed_2026-10-11_oom/`.
- Record the untracked journal-club presentation separately from the code
  changes if it is included; it is not treated as generated analysis output.

## Validation

- Run `git diff --check`; the repaired rows in the generated Aim-5 waypoint
  table intentionally retain a trailing tab for the required empty fifth
  field, so those warnings are expected and the five-column schema is kept.
- Run the focused inference, metabolic-route, and robustness tests using the
  repository `.venv`.
- Run Python compilation checks for changed Python modules and
  `git diff --check`.
- Verify that successful output sidecars contain input, code, package, and
  execution metadata, and that no failed/OOM directory is staged.
- Do not run the multi-hour Stage-05 real-data robustness analysis in this
  commit-validation step.

## Open questions

None. The user approved the Stage-05 scope and exclusion of failed/OOM
directories on 2026-10-12.
