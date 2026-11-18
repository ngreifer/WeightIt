# Tasks

## In Progress

## To Do

- [ ] When *bartisan* reaches CRAN, replace the dbarts-based `method = "bart"` with the bartisan-based one following `_dev/bartisan-migration.md`.

- [ ] Review the notes in `_dev/longitudinal-treatments-vignette-notes.md` and decide on the open questions about `vignette("longitudinal-treatments")`.

- [ ] Review the notes in `_dev/weighting-methods-vignette-notes.md` and decide on the open questions about `vignette("weighting-methods")`.

## Done

- [x] Fix GitHub issue #97: 16 test failures with `NOT_CRAN=true` on x86_64 Linux, where `capabilities("long.double")` is `TRUE` and so `eps` is 1e-5 (on Apple silicon it is 1e-3, which hid them). Running the whole suite here with `eps` forced to 1e-5 and *scs*, *clarabel*, *highs*, and *lpSolve* hidden reproduced every one. Test-only fixes: the `glm()` references in `test-glm_weightit.R`, `test-weights_arg.R`, `test-multinom_weightit.R`, and `test-ordinal_weightit.R` (and the cloglog propensity score model in `test-by_mest.R`) are fit with `epsilon = 1e-12`, because `sandwich()` and `vcov()` on a `glm()` fit use the working weights from the previous coefficients; the `tols` boundary checks in `test-method_ebal.R` and `test-method_energy.R` use a `tols.eps` of 1e-3 as `test-method_optweight.R` already did; the `br = TRUE` variance comparisons with `sandwich()` on a `brglmFit` object use fixed tolerances (.05 and .1), since those differ by design; and the `norm = "linf"` checks in `test-method_optweight.R` skip without *highs*, since *osqp* does not reach exact balance on that problem. Also silenced a *MASS* warning in `test-ordinal_weightit.R`. With the fixes, the whole suite passes at `eps = 1e-5` with those solvers hidden. Whether to add *scs* and *highs* to Suggests is still open; the issue has not been answered on GitHub (2026-11-18).

- [x] Stop `method = "optweight"` from printing the *clarabel* solver's log with `verbose = FALSE`. The three fitting functions in `R/weightit2optweight.R` passed `verbose = TRUE` to *optweight* and relied on `capture.output()`, which cannot catch *clarabel*'s output (it writes to stdout outside R's connections); they now pass the user's `verbose`. Present in 2.1.0, so a NEWS item was added. Test added in `test-method_optweight.R` (mocks *optweight*'s fitting functions), and it fails on the old code (2026-10-27).

- [x] Fix GitHub issue #96: the `norm = "entropy"` test in `test-censoring.R` errored when neither *scs* nor *clarabel* was installed; it now skips in that case. *scs* was not added to Suggests pending a decision; the issue has not been answered on GitHub (2026-10-27).

- [x] `sbps()` gives an informative error when `obj` did not come from `weightit()` (e.g., from `as.weightit()`) and `obj2` or `formula` is missing; previously it failed with "no terms component nor attribute". With both supplied, such objects already worked and give the same result as the `weightit()` objects holding the same weights. Documented on `?sbps`, test added (2026-10-15).

- [x] Fix the two bugs reported from *mediation*: the five model fitting functions failed for `as.weightit()`/`as.weightitMSM()` objects without `s.weights` (`.process_fit()`), and `predict()` after `ordinal_weightit()` failed for one row of `newdata` (`get_pp()`), which also dropped the row name with `level` for `multinom_weightit()`. Probing for other manifestations found a third: `ordinal_weightit()` and ordered treatments in `weightit(method = "glm")` failed when the top category held most units (starting values split above the top category). Tests added in four files; each fails on the old code (2026-10-14).

- [x] Remove the effects of individual treatments from `vignette("longitudinal-treatments")`, so its estimands are the comparisons between regimes only (2026-10-09).

- [x] Break the saturated MSM equation in `vignette("longitudinal-treatments")` over two lines so it fits the pkgdown margins, and replace the worked code in the longitudinal section of `vignette("estimating-effects")` with a pointer to the new vignette (2026-10-09).

- [x] Write a bartisan-based `method = "bart"` (`_dev/weightit2bart_bartisan.R`) with its roxygen page, a test script, timing, and a migration guide (2026-10-08).

- [x] Add cross-references to `vignette("weighting-methods")` and `vignette("longitudinal-treatments")` in the README, the help pages (`weightit()`, `weightitMSM()`, `.cens()`, `msmdata`, `.weightit_methods`, `get_w_from_ps()`, `ESS()`, `trim()`, `calibrate()`, and the eleven method pages), and the other vignettes (2026-10-08).

- [x] Write `vignette("longitudinal-treatments")`, a worked analysis of a longitudinal treatment with `weightitMSM()`, with a second analysis adding censoring weights (2026-10-08).

- [x] Write `vignette("weighting-methods")`, a prose guide to the weighting methods modeled on MatchIt's `vignette("matching-methods")`, with a summary table and a selection flowchart (2026-10-08).
