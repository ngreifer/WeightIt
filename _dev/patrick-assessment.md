# patrick in the Test Suite

How to assess, convert, and verify tests with *patrick* is in the `r-patrick` skill
(`~/.config/agents/skills/r-patrick/`). This file keeps only what is specific to
*WeightIt*.

## Models and Helpers

`test-method_ebal.R` is the model for a `weight.mat` loop, `test-method_energy.R` for a
loop with per-case warning state, and `test-vcov_arg.R` for a table of model fitters
whose calls `update()` re-evaluates. In `helpers.R`, `skip_if_method_unavailable(method)`
skips a case whose method's packages are missing, and
`methods_for(treat.type, installed = FALSE)` builds a method column that keeps those
methods so they skip visibly. A `skip_if_not_installed("rootSolve")` also decides the
default solver for ebal and cbps, so keep it where it is rather than moving it into a case.

## Files Left As They Are

`test-method_glm_re.R` (heterogeneous blocks, no loops), `test-method_bart_re.R` (three
treatment-type tests with different checks), `test-weightitMSM_re.R` (its loops walk
outputs or generate data), and `test-trim_ESS_full_rank.R`, `test-sbps.R`, and
`test-ps_arg.R`, where each test checks a different behavior.
