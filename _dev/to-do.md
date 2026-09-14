## To Do
* Create function for estimating treatment effects to remove marginaleffects as dependency (possibly in a new package)
* Implement RieszBoost for GBM weighting
* Drop `"cobalt.treat"` from `.treat_classes` in `R/treat.R`, once *cobalt* registers its `[` method on `treat` rather than on `cobalt.treat`. It is on `cobalt.treat` only because *WeightIt* 2.0.0 registered a competing `[.treat`; *cobalt* cannot move it while that version is the one on CRAN. See `cobalt/_dev/cens-transition.md`.

## Larger future directions

### Architecture

* **Unify the "weight one sample to fixed targets" problem.** IPCW (`.cens()`), the ATT, survey calibration, and transportability to an external population are all the same shape: weight one group so its covariate moments match a supplied target vector. `weightit2optweight.cens()` already goes through `optweight::optweight.svy.fit()`, and `weightit2ebal.cens()` is the ATE branch restricted to one group. Factoring this into a single internal solver interface (targets + tolerances + a divergence) would let every optimization-based method support all four cases from one implementation, instead of each gaining a bespoke `.cens` branch.
* **Make `Mparts` a real class rather than a list of closures.** The M-estimation parts are passed around as a bare list and manipulated by helpers that must stay mutually consistent (`.invert_num_Mpart()`, `.expand_Mparts_by()`, `.att_out_to_cens()`, the stacking in `.compute_vcov()`). `.expand_Mparts_by()` silently ignores the `Xtreat`/`A` it is handed in favor of ones it closed over, which is invisible from the outside and easy to trip over. A class with `stack()`/`invert()`/`expand()` methods and validity checks would make the contract explicit and the failure loud.
* **Replace the name-mangling method dispatch with a registry.** Methods are found by pasting a suffix onto a function name (`weightit2glm` + `.cens`/`.multi`/`.cont`), with `.weightit_methods` as a parallel table of capabilities that has to be updated by hand in lockstep. A proper registry (or S3 class per method) would let other packages contribute methods, and would remove the need to keep the table and the function names in sync.
* **Split `weightitMSM2cbps()`.** It now handles five treatment types inline in one function, with per-type behavior spread across a covariate-prep `switch()`, a `coef_ind` `switch()`, and five parallel closure lists. It is also less featureful than the point-treatment version (no `link` support, empirical rather than closed-form `Sigma()`). Factor the closure construction into per-type constructors and close the feature gap.
* **Scale the distance-based methods past ~10k units.** `method = "energy"` and `method = "cfd"` build dense n x n matrices and hand them to a QP. Nystrom approximation or random Fourier features would make them usable on large data, where they are currently infeasible.

### New weighting methods

* **Kernel balancing (`method = "kbal"`)** as a first-class method rather than the user-defined-function example in `?method_user`.
* **Transportability / generalizability weights**, targeting an external population's covariate distribution rather than the study sample. Falls out almost free from the targets refactor above.
* **Cross-fitted / debiased weight estimation**, where nuisance models are fit out-of-fold. Relevant to the machine-learning methods (`gbm`, `bart`, `super`), whose plug-in weights have no valid asymptotic variance at present.
* **Weight uncertainty from a posterior** rather than from M-estimation, e.g. propagating BART posterior draws through the weights so effect estimates reflect the full posterior. Would complement the existing M-estimation and bootstrap machinery.
* **Survival-outcome weighting beyond IPCW**, e.g. time-varying censoring weights on a long-format dataset rather than the current one-row-per-unit wide format.

### Diagnostics and testing

* ~~**Target-based balance assessment.**~~ Done in *cobalt* 5.0.0, which compares the weighted uncensored units against the at-risk sample and handles a censoring model sitting among longitudinal treatments. `DESCRIPTION` and `test-censoring.R` now name that version rather than a placeholder.
* **A fast test tier that runs on CRAN.** Under CRAN's own conditions the suite reports `[ FAIL 0 | WARN 0 | SKIP 145 | PASS 46 ]` and finishes in about 6 seconds: 145 blocks skip, and the 46 expectations that survive come from a handful of blocks. A small fast subset that always runs would catch platform-specific breakage that currently only surfaces locally. (Count the skips from the check's `tests/testthat.Rout` rather than by grepping for `skip_on_cran()` -- many `test_that()` calls sit inside `for` loops over methods, so the static count and the runtime count are different numbers.)
* **Coverage in CI**, with the `NOT_CRAN=true` and non-parallel caveats baked in (see below) so the number is meaningful.

## Running the checks on this machine

Two environment traps, both of which have cost whole check runs.

### `R CMD check` hangs at "checking use of S3 registration"

The symptom is a check that sits on that line forever at **0% CPU** — not slow, stopped.
That step preloads the base and recommended packages, `tcltk` among them, and
`loadNamespace("tcltk")` blocks indefinitely here: `/opt/X11` exists but
`XQuartz.app` does not, and `DISPLAY` is set to an XQuartz launchd socket that nothing
is serving, so Tk waits on a connection that will never be accepted.

The workaround is one line:

```sh
unset DISPLAY   # then R CMD build / R CMD check as usual
```

`tcltk` then loads in about a second with a harmless `no display name and no $DISPLAY
environment variable` warning. The real fix is to reinstall XQuartz or drop the stale
`DISPLAY` from the environment.

Worth knowing because the same root cause wears three different faces:

- this hang, in WeightIt's own check;
- `checking Rd files ... WARNING` / `checking for unstated dependencies in examples ...
  WARNING` in *cobalt*, both reading `bal.tab.cem.match.Rd:31: couldn't connect to
  display ""` — `cem` loads `tcltk`. These show up in *cobalt*'s revdep results against
  any version of WeightIt and are not real problems; they disappear with `DISPLAY`
  unset;
- `cem` hanging *cobalt*'s build outright, which `cobalt/_dev/refactor-notes.md`
  § "`R CMD check`" already records as "`cem` must not be installed". That note predates
  knowing why, and could say so now.

Diagnosing it is quick if you know the shape: loading every package in `Suggests` takes
about a second each, so a probe that loads them one at a time under a timeout finds
nothing. The culprit is in the *preload* list, not in our dependencies.

### `R CMD build` cannot run while `revdepcheck` is going

`R CMD build` copies the entire package directory to a temp dir **before** applying
`.Rbuildignore`, so a concurrent `revdepcheck` rewriting `revdep/library.noindex/`
breaks the copy:

```
cp: WeightIt/revdep/library.noindex/CIMTx: No such file or directory
 ERROR
copying to build directory failed
```

`^revdep$` being in `.Rbuildignore` does not help — the exclusion is applied to the copy,
not to the copying. Build from a clean export instead:

```sh
rsync -a --exclude revdep --exclude .git /path/to/WeightIt/ /tmp/src/WeightIt/
R CMD build /tmp/src/WeightIt
```

### Timing NOTEs under load

`checking examples` reports elapsed time, so running a check while `revdepcheck` has
four workers going inflates it. `calibrate`'s example crossed the 5s threshold that way
at 5.29s elapsed against 3.04s CPU; run on its own it is ~2.6s user / ~2.9s elapsed.
Compare CPU against elapsed before believing a timing NOTE.

### The test suite

`devtools::test()` needs `NOT_CRAN=true` to run more than a handful of blocks: 277 of
286 `test_that()` blocks are behind `skip_on_cran()`, which is also why the `checking
tests` step of `R CMD check` finishes in seconds and proves little. Run the suite
separately from the check.
