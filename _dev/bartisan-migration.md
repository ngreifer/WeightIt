# Migrating `method = "bart"` from dbarts to bartisan

A guide for replacing the dbarts-based BART weighting method with the bartisan-based one once *bartisan* is on CRAN. Written 2026-10-08 against WeightIt 2.1.0.9000 and bartisan 0.1.0 (development). Read it start to finish before changing anything; the steps are ordered so the package stays buildable after each one.

## What exists already

- `_dev/weightit2bart_bartisan.R`: the replacement for `R/weightit2bart.R`, with its roxygen page (`?method_bart`) and the four method functions (`weightit2bart()`, `.cens`, `.multi`, `.cont`) plus the helpers they need (`.make_bartisan_args()`, `.scale_s_weights()`, `.make_bartisan_data()`, `.fit_bartisan()`, `.bartisan_family_takes_weights()`, `.bartisan_grid_density()`, `.bartisan_density()`).
- `_dev/bartisan-inject.R`: loads the development package and swaps the new functions into its namespace, patching `.weightit_methods[["bart"]]`, so the new method can be exercised as `method = "bart"` without touching `R/`.
- `_dev/test-weightit2bart_bartisan.R`: testthat expectations covering every treatment type, the options, the argument translations, censoring, random intercepts, and `weightitMSM()`, ending with a comparison against the dbarts implementation. All of them pass as of the date above.
- `_dev/bartisan-timing.R`: timings by sample size.

## How the new method differs from the old one

These are the behavior changes a user would notice, and each needs a NEWS item or a line in the documentation; the roxygen in the new file already describes the new behavior.

| Aspect | dbarts-based (old) | bartisan-based (new) |
|---|---|---|
| Fitting function | `dbarts::bart2()`; `stan4bart::stan4bart()` with random effects | `bartisan::bartisan()` for everything |
| Decision rules | hard (standard BART) | soft by default (`gate = "smoothstep"`); `gate = "hard"` recovers standard BART |
| Binary link | probit (fixed) | `link` argument: `"logit"` (default), `"probit"`, `"cloglog"` |
| Multi-category treatments | one binary BART per level; probabilities not normalized | one `multinomial()` fit; probabilities sum to 1; `link`, `reference`, `replicates` arguments |
| Continuous treatments | BART for the conditional mean, parametric or kernel `density` on the residuals | conditional density estimated by the model: `family = bartisan::dpm()` (default), `gaussian()`, `gaussian_ls()`, or any bartisan family; `density`, `bw`, `adjust`, `kernel`, `n` ignored with a warning |
| Sampling weights | not supported (`s.weights_ok = FALSE`) | supported, scaled to mean 1 before fitting; with `dpm()` they enter only the marginal-density numerator (see Decisions) |
| Random effects | random intercepts and slopes, several grouping factors, via stan4bart | random intercepts only, several grouping factors and nesting, in bartisan itself; a random slope is an error from bartisan |
| Reproducibility | `seed` argument, or `set.seed()` with `n.threads = 1` | `set.seed()`; `seed` ignored with a warning; chains in parallel still reproducible |
| Parallelism | `n.threads` | a *future* plan with `chains > 1`; `n.threads` ignored |
| Control arguments | `bart2()` and `dbartsControl()` names | `bartisan_control()` names; `n.trees`, `n.samples`, `n.burn`, `n.chains` translated |
| `include.obj` | a `bart` fit, or a list of them for multi-category | one `bartisan_fit`, including for multi-category |
| `use.offset` | probit GLM (binary), linear (continuous) | GLM with the chosen `link` (binary), linear (continuous) |
| Missing covariates | `missing = "ind"` | `missing = "ind"` (bartisan's native handling of missing predictors is not used; see Follow-ups) |
| M-estimation | no | no |

## Steps

### 1. Confirm the environment

- Install the CRAN release of *bartisan* and record its version; `packages_needed` and `package_versions_needed` below take it.
- Run `Rscript _dev/test-weightit2bart_bartisan.R` from the package root against the CRAN version. Every block should pass. If *bartisan*'s API changed since 0.1.0, the places most likely to break are `predict(type = "prob")` (multi-category probabilities), `predict(type = "density", newdata, offset, iterations)` (continuous), `fitted()` on a binomial fit, and the names of `bartisan_control()`'s arguments, which `.make_bartisan_args()` reads with `rlang::fn_fmls_names()`.

### 2. Replace the method file

- `git mv _dev/weightit2bart_bartisan.R R/weightit2bart.R` (overwriting the old file). The roxygen page keeps the name `method_bart`, so every existing link to it still resolves.
- Delete `.make_stan4bart_args()` from `R/functions_for_processing.R`; the new file is its only caller. Grep to confirm: `rg -n make_stan4bart_args R/`.
- In `.make_re_data_formula()` (same file), the `bart` argument wrapped the fixed effects in `bart(...)` for stan4bart. The new method calls it with `bart = FALSE`. Remove the argument and the branch unless another method uses it (`rg -n "bart = TRUE" R/` should find nothing after the move).
- Delete `_dev/bartisan-inject.R`, `_dev/bartisan-timing.R`, and, after step 6, `_dev/test-weightit2bart_bartisan.R`.

### 3. DESCRIPTION

- In `Suggests`, replace `dbarts (>= 0.9-29)` and `stan4bart` with `bartisan (>= <CRAN version>)`.
- Nothing else in DESCRIPTION mentions BART.

### 4. The method registry

In `R/weightit_methods.R`, the `"bart"` entry becomes

```r
packages_needed = "bartisan",
package_versions_needed = c(bartisan = "<CRAN version>"),
s.weights_ok = TRUE,
density_ok = FALSE,
```

with `re_ok = TRUE`, `ps = TRUE`, `msm_valid = TRUE`, `subclass_ok = TRUE`, and `stabilize_ok = TRUE` unchanged. Note that `density_ok` is documentation only (nothing in `R/` reads it), so the warning for a supplied `density` lives in `weightit2bart.cont()` itself.

### 5. Documentation outside the method page

- `vignettes/installing-packages.Rmd`, section "Propensity score weighting using BART": replace the *dbarts* and *stan4bart* instructions with *bartisan* (`pak::pak("bartisan")`; development version from `ngreifer/bartisan`; it requires compilation). Drop the paragraph on random effects needing a second package.
- `vignettes/weighting-methods.Rmd`:
  - The BART section (two paragraphs under "Propensity Score Weighting Using BART"): the implementation now uses *bartisan*; decision rules are soft by default (cite Linero and Yang 2018, key `lineroBayesianRegressionTree2018` or whatever Better BibTeX gives; the entry is not yet in `references.bib`); for continuous treatments the conditional density is estimated by a Dirichlet process mixture (cite George et al. 2019 if the entry is added); random intercepts are fit by bartisan directly; sampling weights are now supported, so delete "Sampling weights are not supported"; `set.seed()` replaces `seed`.
  - The machine learning intro says M-estimation is unavailable for these methods, which remains true.
  - The summary table row for `"bart"`: Required package *bartisan*; drop "No sampling weights" from Limitations; consider adding "soft trees" to Strengths.
  - The flowchart SVG does not mention packages and needs no change.
- `vignettes/longitudinal-treatments.Rmd` and `vignettes/WeightIt.Rmd` do not mention BART.
- `README.Rmd`: the methods table row for `"bart"` is unchanged (all four columns stay checked). Regenerate `README.md` with `rmarkdown::render("README.Rmd")` only if the text changes.
- `R/weightit.R`: the methods table and the `formula` parameter ("random effects terms ... for methods that accept them") need no change. Check `?weightit.fit`'s `ps` item, which lists the methods that return a propensity score; `"bart"` stays on the list.
- `R/weightit2super.R`'s See Also ("stacking predictions from several machine learning methods, including BART") is still accurate.
- `NEWS.md`: one item under the development version, in the house style (state the new behavior; one short clause for the old). Cover, in this order: the backend change and soft trees; multi-category treatments fit as one multinomial model; continuous treatments' density estimated by `dpm()` with `family` replacing `density`; sampling weights now accepted; random intercepts fit natively and random slopes no longer supported; `set.seed()` replacing `seed`; the translated `n.*` arguments; `include.obj` returning a single fit for multi-category treatments; the default link changing from probit to logit.
- Run `devtools::document()` and check `man/method_bart.Rd` renders (`?method_bart` after `load_all()`), in particular the `\pkgfun{bartisan}{...}` links, which require *bartisan* installed.

### 6. Tests

- `tests/testthat/test-method_bart.R`: replace with the expectations in `_dev/test-weightit2bart_bartisan.R`, converted to the package's conventions: `skip_if_not_installed("bartisan")` and `skip_on_cran()` at the top of each block, the fixture read with `test_path("fixtures", "test_data.rds")`, and the fast settings (`num_draws = 150L, num_burn = 100L`) kept, since the defaults (800 and 200) make the file slow. The old file's `patrick` grids over estimands and `density` values can stay as grids over estimands; the `density` grid no longer applies. Drop the test of `bart2()` argument exceptions (`samplerOnly`, `keepCall`, and so on), which no longer mean anything. Keep the "balance improved" checks rather than exact-balance checks.
- `tests/testthat/test-method_bart_re.R`: rewrite for bartisan random intercepts: the fit is a `bartisan_fit` with a `ranef` component; two grouping factors work; a random slope errors with a message from bartisan. Remove the stan4bart argument-translation checks.
- `tests/testthat/test-censoring.R`, around line 410: `skip_if_not_installed("dbarts")` becomes `"bartisan"`; the `n.trees`, `n.samples`, `n.burn`, `n.chains`, `n.threads` arguments still work through translation but are better rewritten as `num_trees`, `num_draws`, `num_burn`.
- `tests/testthat/test-weightitMSM.R`: line 98 passes `seed = 123` and `n.threads` for bart, which now produce a warning (`seed`) and are ignored (`n.threads`); replace with `set.seed()` before the fit and the `num_*` arguments. Line 270's test that `s.weights` errors for bart must be inverted or removed, since sampling weights are now supported.
- Run `devtools::test(filter = "bart")`, then `devtools::test(filter = "censoring|weightitMSM")`, then the full suite before release.

### 7. Checks

- `devtools::check()`; the examples on `?method_bart` are wrapped in `\donttest{}` and gated on `rlang::is_installed("bartisan")`.
- Build the pkgdown site and read the `method_bart` reference page and the three vignettes that mention BART.
- Add `_dev/formula.test.R` and `_dev/multilevel-data.Rmd` to the list of things to clean up or leave; both mention dbarts but are not part of the package.

## Decisions made in the new implementation, to confirm

1. **Default link is logit**, matching `method = "glm"` and bartisan's binomial default, where the old method's probit was dbarts's only option. Users comparing old and new results will see small differences from this alone.
2. **Sampling weights with `dpm()`.** `bartisan::dpm()` refuses prior weights (a weight would have to be a multiplicity in the Dirichlet process). Rather than error, the method fits the conditional density without them and uses them only in the weighted average that gives the marginal density; the help page says so and points to `gaussian()` and `gaussian_ls()` for fits that use them. The alternatives were to error, which would make the default family unusable with sampling weights, or to switch families silently, which seemed worse. If you would rather the default family be `gaussian()` whenever `s.weights` vary, that is a two-line change in `weightit2bart.cont()`.
3. **The marginal density numerator** is computed as in `.get_w_from_gps_internal_cont()`: the conditional density at 50 grid values spanning the treatment, averaged over units with the sampling weights, and interpolated with `spline()` to the observed values. Because `predict(type = "density")` on new rows has to evaluate every tree in every draw for every row, the draws are thinned to 25 (`max_draws`) and at most 2000 units are averaged over (`max_units`, a random subset), so each grid value averages up to 50,000 density evaluations and the cost stops growing with the sample size; both are internal defaults of `.bartisan_grid_density()` and could be exposed if wanted. The denominator is computed on the fitted data with every draw, which is cheap because the fit stores its additive predictor.
4. **Translated arguments.** `n.trees`, `n.samples`, `n.burn`, and `n.chains` are translated so existing code and the tests keep working; `n.threads` is dropped silently and `seed` with a warning. A later release could deprecate the translations.
5. **`.make_covs_full_rank()`** is still applied for binary and multi-category treatments, as before, though BART does not need it.
6. **The data frame is embedded in the call** that `bartisan()` stores, as the old method embedded the matrix in `bart2()`'s call, so a fit kept with `include.obj = TRUE` carries a copy of the covariates. The weights and offset are referred to by column name (`.s.weights`, `.offset`) so they are not copied into the call separately.
7. **bartisan's own warnings pass through**, in particular the one about the leaf scale when the treatment is nearly separable, which recommends `update_sigma_mu = FALSE`; that argument can be passed through `weightit()`.
8. **Missing covariates** still use `missing = "ind"`. bartisan handles missing predictors natively (missingness incorporated in attributes), which could be offered as a second `missing` option later; it would need `.weightit_methods[["bart"]]$missing` extended and the indicator step skipped.

## Results of the comparison (lalonde, defaults, one chain)

| Case | Max imbalance, unweighted | dbarts | bartisan | ESS dbarts | ESS bartisan | Seconds dbarts | Seconds bartisan |
|---|---|---|---|---|---|---|---|
| Binary ATE | 2.27 | .74 | .67 | 351 | 286 | .2 | .8 |
| Binary ATT | 2.27 | .42 | .36 | 240 | 220 | .2 | .8 |
| Multi-category ATE | .62 | .19 | .29 | 352 | 340 | .5 | 3.0 |
| Continuous (re75) | .55 | .48 | .52 | 38 | 592 | .1 | 2.0 |

Imbalance is the largest absolute standardized mean difference, or treatment-covariate correlation for the continuous case. The continuous row shows the main practical gain: the DPM density gives weights with an effective sample size of 593 of 614 units where the normal density on BART residuals gave 38, at similar balance (re75 is a poor continuous treatment, being zero-inflated, so neither balances it well; on the test fixture's `Ac` the correlation drops from .35 to .14 with ESS 436 of 750). Multi-category balance was somewhat worse on lalonde in this run; the multinomial fit is one model for three categories where the old method fit three separate ones, and one run on one dataset is not a verdict.

## Timing by sample size (test fixture, defaults, one chain)

| n | Binary | Multi-category | Continuous | of which the `dpm()` fit |
|---|---|---|---|---|
| 750 | 1.3 s | 4.1 s | 3.5 s | 1.8 s |
| 3000 | 4.0 s | 14.5 s | 20.5 s | 4.2 s |

For comparison, the dbarts implementation fits lalonde (614 units) in .2, .5, and .1 seconds. The continuous case is dominated by the marginal-density numerator: `predict()` on new rows re-evaluates the soft trees for every draw, whereas the denominator reuses the additive predictor the fit stores. With the thinning described under Decisions, that part is about 15 seconds and no longer grows with the sample size beyond 2000 units. The multi-category fit is the slowest of the three families per observation, since it carries one forest per level. All three are bartisan at its defaults (800 draws, 200 warmup, 50 trees, soft rules); `gate = "hard"` or fewer draws will cut the times, as will `chains > 1` with a future plan on a machine with spare cores.

## Follow-ups worth considering after the migration

- Expose `max_draws` and `max_units` of the marginal-density computation, or compute the marginal density without re-evaluating the trees: the grid rows have the same predictors as the fitted data, so their additive predictor is already in `fit$eta`, and only the error density at `a - eta` is needed. For `dpm()` that is one shifted error density per draw, which `error_density()` reports only as a posterior average; a bartisan function returning the per-draw mixture density at given values, or `predict(type = "density")` accepting a response vector for the fitted rows, would make the numerator nearly free. The same applies to `gaussian()` (`fit$aux` has the sigma draws) but not to `gaussian_ls()`, whose scale varies with the predictors.
- Offer bartisan's native missing-value handling as `missing = "mia"` (or similar).
- Decide whether `gate = "hard"` should be the default for speed, with soft rules opt-in; the current default follows bartisan.
- Decide whether to add the two references the roxygen page cites (Linero 2025; Linero and Yang 2018; George et al. 2019) to `vignettes/references.bib` for the weighting-methods vignette.
