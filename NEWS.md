WeightIt News and Updates
======

# `WeightIt` 1.3.0

* Added `anova()` methods for `glm_weightit`, `multinom_weightit`, `ordinal_weightit`, and `coxph_weightit` objects to perform Wald tests for comparing nested models. The models do not have to be symbolically nested.

* Added the new user-facing object `.weightit_methods`, which contains information on each method and the options allowed with it. This is used within `WeightIt` for checking arguments but can also be used by other package developers who call functions in `WeightIt`. See `help(".weightit_methods")` for details.

* `plot.weightit()` can be used with `method = "optweight"` to display the dual variables.

* `missing` no longer allows partial matching.

* `moments` can now be set to 0 when `quantile` is supplied to ensure balance on the quantiles without the moments for the methods that accepts `quantiles`. Thanks to @BERENZ for the suggestion.

* For `ordinal_weightit` objects, `summary()` now has the option to omit thresholds from the output.

* Fixed a bug in `ordinal_weightit()` where the Hessian (and therefore the HC0 robust variance) were calculated incorrectly when come coefficients were aliased (i.e., due to linearly dependent predictors).

* Fixed a bug in `print.summary.glm_weightit()` when confidence intervals were requested. A new printing function is used that produces slightly nicer tables.

* Fixes to vignettes and tests to satisfy CRAN checks.

* Minor bug, performance, and readability fixes.

# `WeightIt` 1.2.0

* Added two new functions, `multinom_weightit()` and `ordinal_weightit()` for multinomial logistic regression and ordinal regression with capabilities to estimate a covariance matrix that accounts for estimation of the weights using M-estimation. Previously, multinomial logistic regression could be requested using `glm_weightit()` with `family = "multinomial"`; this has been deprecated.

* M-estimation can now be used for weighting with ordinal regression for weights with multi-category ordered treatments with `method = "glm"`.

* M-estimation can now be used with bias-reduced regression as implemented in `brglm2` for the propensity score (`method = "glm"` with `link = "br.{.}"`) and for the outcome model (`glm_weightit()` with `br = TRUE`). Thanks to Ioannis Kosmidis for supplying some starter code to implement this.

* For any weighting methods with continuous treatments that support a `density` argument to specify the numerator and denominator densities of the weights, `density` can now be specified as `"kernel"` to request kernel density estimation. Previously, this was requested by setting `use.kernel = TRUE`, which is now deprecated.

* Standard errors are now correctly computed when an offset is included in `glm_weightit()`. Thanks to @zeynepbaskurt. (#63)

* Improved robustness of `get_w_from_ps()` to propensity scores of 0 and 1.

* Updates to `weightit()` with `method = "gbm"`:

    * `use.offset` is now tunable.
    * The same random seed is used across specifications as requested by @mitchellcameron123. (#64)
    * For binary and multi-category treatments with cross-validation used as the criterion, `class.stratify.cv` is now set to `TRUE` by default to stratify on treatment.
    * For continuous treatments, the default density now corresponds to the distribution requested.
    * `plot()` can be used on the output of `weightit()` to display the results of the tuning process; see `help("plot.weightit")` for details.
    * Fixed a bug where `distribution` was not included in the output when tuned.
    * Fixed a bug when propensity scores were estimated to be 0 or 1. Thanks to @mitchellcameron123. Propensity scores are now shifted slightly away from 0 or 1. (#64)

* When using `weightit()` with `method = "super"` for binary and multi-category treatments, cross-validation now stratifies on treatment, as recommended by [Phillips et al. (2023)](https://doi.org/10.1093/ije/dyad023).

* Fixed a bug and clarified some error messages when using ordered treatments with `method = "glm"`. Thanks to Steve Worthington for pointing them out.

* Updated the help page of `get_w_from_ps()` to include formulas for the weights.

# `WeightIt` 1.1.0

* Added a new function, `coxph_weightit()`, for fitting Cox proportional hazards models in the weighted sample, with the option of accounting for estimation of the weights in computing standard errors via bootstrapping. This function uses the `summary()` and `print()` methods for `glm_weightit` objects, which are different from those for `coxph` objects.

* `glm_weightit()` gets a new `print()` method that omits some invalid statistics displayed by the `print()` method for `glm` objects and displays the type of standard error estimated.

* `summary.glm_weightit()` (which is also used for `coxph_weightit` objects) gets a new argument, `transform`, which can be used to transform the displayed coefficients and confidence interval bounds (if requested), e.g., by exponentiating them.

* M-estimation is now supported for `method = "glm"` with continuous treatments.

* A new estimator is now used for `method = "cbps"` with longitudinal treatments (i.e., using `weightitMSM()`). Previously, the weights from CBPS applied to each time point were multiplied together. Now, balance at all time points is optimized using a single set of weights. This implementation is close to that described by [Huffman and van Gameren (2018)](https://doi.org/10.1515/jci-2017-0002), not that of [Imai and Ratkovic (2015)](https://doi.org/10.1080/01621459.2014.956872).

* A new estimator is now used for `method = "cbps"` with continuous treatments. The unconditional mean and variance are now included as parameters to be estimated. For the just-identified CBPS, this will typically improve balance, but results will depart from those found using `CBPS::CBPS()`.

* For point treatments (i.e., using `weightit()`), the `stabilize` argument has some new behavior. It can now be be specified as a formula, and the stabilization factor is estimated separately and included in the M-estimation if allowed. It can now only be used when `estimand = "ATE"` (weights for other estimands should not be stabilized).

* For binary treatments with `method = "glm"`, `link` can now be specified as `"flic"` or `"flac"` to use Firth corrected logistic regression as implemented in the `logistf` package.

* With `method = "gbm"`, an error is now thrown if `criterion` (formerly known as `stop.method`) is supplied as anything other than a string.

* For binary and continuous treatments with `method = "gbm"`, a new argument, `use.offset`, can be supplied, which, if `TRUE`, uses the linear predictor from a generalized linear model as an offset in the boosting model, which can improve performance.

* Added a section on conducting moderation analysis to the estimating effect vignette (`vignette("estimating-effects")`).

* Fixed a bug when using M-estimation for sequential treatments with `weightitMSM()` and `stabilize = TRUE`. Standard errors incorrectly accounted for estimation of the stabilization factor; they are now correct.

* Fixed a bug when using `method = "ipt"` for the ATE.

* Fixed a bug when some coefficients were aliased for `glm_weightit()`. Thanks to @kkwi5241.

* Updated kernel balancing example in `method_user`.

* Improved warnings and errors for bad models throughout the package.

# `WeightIt` 1.0.0

* Added a new function, `glm_weightit()` (along with wrapper `lm_weightit()`) and associated methods for fitting generalized linear models in the weighted sample, with the option of accounting for estimation of the weights in computing standard errors via M-estimation or two forms of bootstrapping. `glm_weightit()` also supports multinomial logistic regression in addition to all models supported by `glm()`. Cluster-robust standard errors are supported, and output is compatible with any functions that accept `glm()` objects. Not all weighting methods support M-estimation, but for those that do, a new component is added to the `weightit` output object. Currently, GLM propensity scores, entropy balancing, just-identified CBPS, and inverse probability tilting (described below) support M-estimation-based standard errors with `glm_weightit()`.

* Added inverse probability tilting (IPT) as described by Graham, Pinto, and Egel (2012), which can be requested by setting `method = "ipt"`. Thus is similar to entropy balancing and CBPS in that it enforces exact balance and yields a propensity score, but has some theoretical advantages to both methods. IPT does not rely on any other packages and runs very quickly.

* Estimating covariate balancing propensity score weights (i.e., `method = "cbps"`) no longer depends on the `CBPS` package. The default is now the just-identified versions of the method; the over-identified version can be requested by setting `over = TRUE`. The ATT for multi-category treatments is now supported, as are arbitrary numbers of treatment groups (`CBPS` only natively support up to 4 groups and only the ATE for multi-category treatments). For binary treatments, generalized linear models other than logistic regression are now supported (e.g., probit or Poisson regression).

* New function `calibrate()` to apply Platt scaling to calibrate propensity scores as recommended by [Gutman et al. (2024)](https://doi.org/10.1097/EDE.0000000000001733).

* A new argument `quantile` can be supplied to `weightit()` with all the methods that accept `moments` and `int` (`"ebal"`, `"cbps"`, `"ipt"`, `"npcbps"`, `"optweight"`, and `"energy"`). This allows one to request balance on the quantiles of the covariates, which can add some robustness as demonstrated by [Beręsewicz (2023)](https://arxiv.org/abs/2310.11969).

* `as.weightit()` now has a method for `weightit.fit` objects, which now have additional components included in the output.

* `trim()` now has a `drop` argument; setting to `TRUE` sets the weights of all trimmed units to 0 (effectively dropping them).

* When using `weightit()` with a continuous treatment and a `method` that estimates the generalized propensity score (e.g., `"glm"`, `"gbm"`, `"super"`), sampling weights are now be incorporated into the density when `use.kernel = FALSE` (the default) when supplied to `s.weights`. Previously they were ignored in calculating the density, but have always been and remain used in the modeling the treatment (when allowed).

* Fixed a bug when `criterion` was not specified when using `method = "gbm"`.

* Fixed a bug when `ps` was supplied for continuous treatments. Thanks to @taylordunn. (#53)

* Warning messages now display immediately rather than at the end of evaluation.

* The vignettes have been changed to use a slightly different estimator for weighted g-computation. The estimated weights are no longer to be included in the call to `avg_comparisons()`, etc.; that is, they are only used to fit the outcome model. This makes the estimators more consistent with other software, including `teffects ipwra` in Stata, and most of the literature on weighted g-computation. Note this will not effect any estimates for the ATT or ATC and will only yield at most minor changes for the ATE. For other estimands (e.g., ATO), the weights are still to be included.

* The word "multinomial" to describe treatments with more than two categories has been replaced with "multi-category" in all documentation and messages.

* Transferred all help files to Roxygen and reorganized package scripts.

* Reorganization of some functions.

# `WeightIt` 0.14.2

* Fixed a bug when using `estimand = "ATC"` with multi-category treatments. (#47)

* Fixed a bug in the Estimating Effects vignette. (#46)

# `WeightIt` 0.14.1

* `cobalt` version 4.5.1 or greater is now required.

* Fixed a bug when using balance Super Learner with `cobalt` 4.5.1.

* Added a section to the Estimating Effects vignette (`vignette("estimating-effects")`) on estimating the effect of a continuous treatment after weighting.

# `WeightIt` 0.14.0

* Added energy balancing for continuous treatments, requested using `method = "energy"`, as described in [Huling et al. (2023)](https://doi.org/10.1080/01621459.2023.2213485). These weights minimize the distance covariance between the treatment and covariates while maintaining representativeness. This method supports exact balance constraints, distributional balance constraints, and sampling weights. The implementation is similar to that in the `independenceWeights` package. See `?method_energy` for details.

* Added a new vignette on estimating effects after weighting, accessible using `vignette("estimating-effects", package = "WeightIt")`. The new workflow relies on the `marginaleffects` package. The main vignette (`vignette("WeightIt")`) has been modernized as well.

* Added a new dataset, `msmdata`, to demonstrate capabilities for longitudinal treatments. `twang` is no longer a dependency.

* Methods that use a balance criterion to select a tuning parameter, in particular GBM and balance Super Learner, now rely on `cobalt`'s `bal.init()` and `bal.compute()` functionality, which adds new balance criteria. The `stop.method` argument for these functions has been renamed to `criterion` and `help("stop.method")` has been removed; the same page is now available at `help("bal.compute", package = "cobalt")`, which describes the additional statistics available. This also fixes some bugs that were present in some balance criteria.

* Renamed `method = "ps"` to `method = "glm"`. `"ps"` continues to work as it always had for back compatibility. `"glm"` is a more descriptive name since many methods use propensity scores; what distinguishes this method is that it uses generalized linear models.

* Using `method = "ebcw"` for empirical balancing calibration weighting is no longer available because the `ATE` package has been removed. Use `method = "ebal"` for entropy balancing instead, which is essentially identical.

* Updated the `trim()` documentation to clarify the form of trimming that is implemented (i.e., winsorizing). Suggested by David Novgorodsky.

* Fixed bugs when some `s.weights` are equal to zero with `method = "ebal"`, "`cbps"`, and `"energy"`. Suggested by @statzhero. (#41)

* Improved performance of `method = "energy"` for the ATT.

* Fixed a bug when using `method = "energy"` with `by`.

* With `method = "energy"`, setting `int = TRUE` automatically sets `moments = 1` if unspecified.

* Errors and warnings have been updated to use `chk`.

* The missingness indicator approach now imputes the variable median rather than 0 for missing values. This will not change the performance of most methods, but change others, and doesn't affect balance assessment.

# `WeightIt` 0.13.1

* For ordinal multi-category treatments, setting `link = "br.logit"` now uses `brglm2::bracl()` to fit a bias-reduced ordinal regression model.

* Added the vignette "Installing Supporting Packages" to explain how to install the various packages that might be needed for `WeightIt` to use certain methods, including when the package is not on CRAN. See the vignette at `vignette("installing-packages")`.

* Fixed a bug that would occur when a factor or character predictor with a single level was passed to `weightit()`.

* Improved the code for entropy balancing, fixing a bug when using `s.weights` with a continuous treatment and improving messages when the optimization fails to converge. (#33)

* Improved robustness of documentation to missing packages.

* Updated the logo, thanks to [Ben Stillerman](https://stillben.com).

# `WeightIt` 0.13.0

* Fixed a bug that would occur when the `formula.tools` package was loaded, which would occur most commonly when `logistf` was loaded. It would cause the error `The treatment and covariates must have the same number of units.` (#25)

* Fixed a bug where the `info` component would not be included in the output of `weightit()` when using `method = "super"`.

* Added the ability to specify `num.formula` as a list of formulas in `weightitMSM()`. This is primarily to get around the fact that when `stabilize = TRUE`, a fully saturated model with all treatments is used to compute the stabilization factor, which, for many time points, is time-consuming and may be impossible (especially if not all treatment combinations are observed). Thanks to @maellecoursonnais for bringing up this issue (#27).

* `ps.cont()` has been retired since the same functionality is available using `weightit()` with `method = "gbm"` and in the `twangContinuous` package.

* With `method = "energy"`, a new argument, `lambda`, can be supplied, which puts a penalty on the square of the weights to control the effective sample size. Typically this is not needed but can help when the balancing is too aggressive.

* With `method = "energy"`, `min.w` can now be negative, allowing for negative weights.

* With `method = "energy"`, `dist.mat` can now be supplied as the name of a method to compute the distance matrix: `"scaled_euclidean"`, `"mahalanobis"`, or `"euclidean"`.

* Support for negative weights added to `summary()`. Negative weights are possible (though not by default) when using `method = "energy"` or `method = "optweight"`.

* Fixed a bug where `glm()` would fail to converge with `method = "ps"` for binary treatments due to bad starting values. (#31)

* `miss = "saem"` can once again be used with `method = "ps"` when missing values are present in the covariates.

* Fixed bugs with processing input formulas.

* An error is now thrown if an incorrect `link` is supplied with `method = "ps"`.

# `WeightIt` 0.12.0

* The use of `method = "twang"` has been retired and will now give an error message. Use `method = "gbm"` for nearly identical functionality with more options, as detailed at `?method_gbm`.

* With multinomial treatments with `link = "logit"` (the default), if the `mclogit` package is installed, it can be requested for estimating the propensity score by setting the option `use.mclogit = TRUE`, which uses `mclogit::mblogit()`. It should give the same results as the default, which uses `mlogit`, but can be faster and so is recommended.

* Added a `plot()` method for `summary.weightitMSM` objects that functions just like `plot.summary.weightit()` for each time point.

* Fixed a bug in `summary.weightit()` where the labels of the top weights were incorrect. Thanks to Adam Lilly.

* Fixed a bug in `sbps()` when using a stochastic search (i.e., `full.search = FALSE` or more than 8 moderator levels). (#17)

* Fixed a bug that would occur when all weights in a treatment group were `NA`. Bad weights (i.e., all the same) now produce a warning rather than an error so the weights can be diagnosed manually. (#18)

* Fixed a bug when using `method = "energy"` with `estimand = "ATE"` and `improved = TRUE` (the default). The between-treatment energy distance contribution was half of what it should have been; this has now been corrected.

* Added L1 median measure as a balance criterion. See `?stop.method` for details.

* Fixed a bug where logical treatments would yield an error. (#21)

* Fixed a bug where `Warning: Deprecated` would appear sometimes when `purrr` (part of the `tidyverse`) was loaded. (#22) Thanks to MrFlick on StackOverflow for the [solution](https://stackoverflow.com/a/66897921/6348551).

# `WeightIt` 0.11.0

* Added support for estimating propensity scores using Bayesian additive regression trees (BART) with `method = "bart"`. This method fits a BART model for the treatment using functions in the `dbarts` package to estimate propensity scores that are used in weights. Binary, multinomial, and continuous treatments are supported. BART uses Bayesian priors for its hyperparameters, so no hyperparameter tuning is necessary to get well-performing predictions.

* Fixed a bug when using `method = "gbm"` with `stop.method = "cv{#}"`.

* Fixed a bug when setting `estimand = "ATC"` for methods that produce a propensity score. In the past, the output propensity score was the probability of being in the control group; now, it is the probability of being in the treated group, as it is for all other estimands. This does not affect the weights.

* Setting `method = "twang"` is now deprecated. Use `method = "gbm"` for improved performance and increased functionality. `method = "twang"` relies on the `twang` package; `method = "gbm"` calls `gbm` directly.

* Using `method = "ebal"` no longer requires the `ebal` package. Instead, `optim()` is used, as it has been with continuous treatments. Balance is a little better, but some options have been removed. 

* When using `method = "ebal"` with continuous treatments, a new argument, `d.moments`, can now be specified. This controls the number of moments of the covariate and treatment distributions that are constrained to be the same in the weighted sample as they are in the original sample. Vegetabile et al. (2020) recommend setting `d.moments` to at least 3 to ensure generalizability and reduce bias due to effect modification.

* Made some minor changes to `summary.weightit()` and `plot.summary.weightit()`. Fixed how negative entropy was computed.

* The option `use.mnlogit` in `weightit()` with multi-category treatments and `method = "ps"` has been removed because `mnlogit` appears uncooperative.

* Fixed a bug (#16) when using `method = "cbps"` with factor variables, thanks to @danielebottigliengo.

* Fixed a bug when using binary factor treatments, thanks to Darren Stewart.

* Cleaned up the documentation.

# `WeightIt` 0.10.2

* Fixed a bug where treatment values were accidentally switched for some methods.

# `WeightIt` 0.10.1

* With `method = "gbm"`, added the ability to tune hyperparameters like `interaction.depth` and `distribution` using the same criteria as is used to select the optimal tree. A summary of the tuning results is included in `info` in the `weightit` output object.

* Fixed a bug where `moments` and `int` were ignored unless both were specified.

* Effective sample sizes now print only up to two digits (believe me, you don't need three) and print more cleanly with whole numbers.

* Fixed a bug when using `by`, thanks to @frankpopham. (#11)

* Fixed a bug when using `weightitMSM` with methods that process `int` and `moments` (though you probably shouldn't use them anyway). Thanks to Sven Reiger.

* Fixed a bug when using `method = "npcbps"` where weights could be excessively small and mistaken for all being the same. The weights now sum to the number of units.

# `WeightIt` 0.10.0

* Added support for energy balancing with `method = "energy"`. This method minimizes the energy distance between samples, which is a multivariate distance measure. This method uses code written specifically for `WeightIt` (i.e., it does not call a package specifically designed for energy balancing) using the `osqp` package for the optimization (same as `optweight`). See Huling & Mak (2020) for details on this method. Also included is an option to require exact balance on moments of the covariates while minimizing the energy distance. The method works for binary and multinomial treatments with the ATE, ATT, or ATC. Sampling weights are supported. Because the method requires the calculation and manipulation of a distance matrix for all units, it can be slow and/or memory intensive for large datasets.

* Improvements to `method = "gbm"` and to `method = "super"` with `SL.method = "method.balance"`. A new suite of `stop.method`s are allowed. For binary treatments, these include the energy distance, sample Mahalanobis distance, and pseudo-R2 of the weighted treatment model, among others. See `?stop.method` for allowable options. In addition, performance for both is quite a bit faster.

* With multinomial treatments with `link = "logit"` (the default), if the `mnlogit` package is installed, it can be requested for estimating the propensity score by setting the option `use.mnlogit = TRUE`. It should give the same results as the default, which uses `mlogit`, but can be faster for large datasets.

* Added option `estimand = "ATOS"` for the "optimal subset" treatment effect as described by Crump et al. (2009). This estimand finds the subset of units who, with ATE weights applied, yields a treatment effect with the lowest variance, assuming homoscedasticity (and other assumptions). It is only available for binary treatments with `method = "ps"`. In general it makes more sense to use `estimand = "ATO"` if you want a low-variance estimate and don't care about the target population, but I added this here for completeness. It is available in `get_w_from_ps()` as well.

* `make_full_rank()` is now faster.

* Cleaning up of some error messages.

* Fixed a bug when using `link = "log"` for `method = "ps"` with binary treatments.

* Fixed a bug when using `method = "cbps"` with continuous treatments and sampling weights. Previously the returned weights included the sampling weights multiplied in; now they are separated, as they are in all other scenarios and for all other methods.

* Improved processing of non-0/1 binary treatments, including for `method = "gbm"`. A guess will be made as to which treatment is considered "treated"; this only affects produced propensity scores but not weights.

* Changed default value of `at` in `trim()` from .99 to 0.

* Added output for the number of weights equal to zero in `summary.weightit`. This can be especially helpful when using `"optweight"` or `"energy"` methods or when using `estimand = "ATOS"`.

# `WeightIt` 0.9.0

* Added support for entropy balancing (`method = "ebal"`) for continuous treatments as described by Tübbicke (2020). Relies on hand-written code contributed by Stefan Tübbicke rather than another R package. Sampling weights and base weights are both supported as they are with binary and multi-category treatments.

* Added support for Balance SuperLearner as described by Pirracchio and Carone (2018) with `method = "super"`. Rather than using NNLS to choose the optimal combination of predictions, you can now optimize balance. To do so, set `SL.method = "method.balance"`. You will need to set an argument to `stop.method`, which works identically to how it does for `method = "gbm"`. For example, for `stop.method = "es.max"`, the predicted values given will be the combination of predicted values that minimizes the largest absolute standardized mean difference of the covariates in the sample weighted using the predicted values as propensity scores.

* Changed some of the statistics displayed when using `summary()`: the weight ratio is gone (because weights can be 0, which is not problematic but would explode the ratio), and the mean absolute deviation and entropy of the weights are now present.

* Added `crayon` for prettier printing of `summary()` output.

# `WeightIt` 0.8.0

* Formula interfaces now accept `poly(x, .)` and other matrix-generating functions of variables, including the `rms`-class-generating functions from the `rms` package (e.g., `pol()`, `rcs()`, etc.) (the `rms` package must be loaded to use these latter ones) and the `basis`-class-generating functions from the `splines` package (i.e., `bs()` and `ns()`). A bug in an early version of this was found by @ahinton-mmc.

* Added support for marginal mean weighting through stratification (MMWS) as described by Hong (2010, 2012) for `weightit()` and `get_w_from_ps()` through the `subclass` argument (see References at `?get_w_from_ps`). With this method, subclasses are formed based on the propensity score and weights are computed based on the number of units in each subclass. MMWS can be used with any method that produces a propensity score. The implementation here ensures all subclasses have a least one member by filling in empty subclasses with neighboring units.

* Added `stabilize` option to `get_w_from_ps()`. 

* A new `missing` argument has been added to `weightit()` to choose how missing data in the covariates is handled. For most methods, only `"ind"` (i.e., missing indicators with single-value imputation) is allowed, but for `"ps"`, `"gbm"`, and `"twang"`, other methods are possible. For `method = "ps"`, a stochastic approximation of the EM algorithm (SAEM) can be used through the `misaem` package by setting `missing = "saem"`. 

* For continuous treatments with the `"ps"`, `"gbm"`, and `"super"` methods (i.e., where the conditional density of the treatment needs to be estimated), the user can now supply their own density as a string or function rather than using the normal density or kernel density estimation. For example, to use the density of the t-distribution with 3 degrees of freedom, one can set `density = "dt_3"`. T-distributions often work better than normal distributions for extreme values of the treatment.

* Some methods now have an `info` component in the output object. This contains information that might be useful in diagnosing or reporting the method. For example, when `method = "gbm"`, `info` contains the tree that was used to compute the weights and the balance resulting from all the trees, which can be plotted using `plot()`. When `method = "super"`, `info` contains the coefficients in the stacking model and the cross-validation risk of each of the component methods.

* For `method = "gbm"`, the best tree can be chosen using cross validation rather than balance by setting `stop.method = "cv5"`, e.g., to do 5-fold cross-validation.

* For `method = "gbm"`, a new optional argument `start.tree` can be set to select the tree at which balance begins to be computed. This can speed things up when you know that the best tree is not within the first 100 trees, for example.

* When using `method = "gbm"` with multi-category treatments and estimands other than the `ATE`, `ATT`, or `ATC` are used with standardized mean differences as the stopping rule, the mean differences will be between the weighted overall sample and each treatment group. Otherwise, some efficiency improvements.

* When using `method = "ps"` with multi-category treatments, the use of `use.mlogit = FALSE` to request multiple binary regressions instead of multinomial regression is now documented and an associated bug is now fixed, thanks to @ahinton-mmc.

* When use `method = "super"`, one can now set `discrete = TRUE` to use discrete SuperLearner instead of stacked SuperLearner, but you probably shouldn't.

* `moments` and `int` can now be used with `method = "npcbps"`.

* Performance enhancements.

# `WeightIt` 0.7.1

* Fixed bug when using `weightit()` inside another function that passed a `by` argument explicitly. Also changed the syntax for `by`; it must now either be a string (which was always possible) or a one-sided formula with the stratifying variable on the right-hand side. To use a variable that is not in `data`, you must use the formula interface. 

* Fixed bug when trying to use `ps` with `by` in `weightit()`.

# `WeightIt` 0.7.0

* Added new `sbps()` function for estimating subgroup balancing propensity score weights, including both the standard method and a new smooth version.

* Setting `method = "gbm"` and `method = "twang"` will now do two different things. `method = "gbm"` uses `gbm` and `cobalt` functions to estimate the weights and is much faster, while `method = "twang"` uses `twang` functions to estimate the weights. The results are similar between the two methods. Prior to this version, `method = "gbm"` and `method = "twang"` both did what `method = "twang"` does now. 

* Bug fixes when `stabilize = TRUE`, thanks to @ulriksartipy and Sven Rieger.

* Fixes for using `base.weight` argument with `method = "ebal"`. Now the supplied vector should have a length equal to the number of units in the dataset (in contrast to its use in `ebalance`, which requires a length equal to the number of control units).

* Restored dependency on `cobalt` for examples and vignette.

* When `method = "ps"` and the treatment is ordered (i.e., ordinal), `MASS::polr()` is used to fit an ordinal regression. Make the treatment un-ordered to to use multinomial regression instead.

* Added support for using bias-reduced fitting functions when `method = "ps"` as provided by the `brglm2` package. These can be accessed by changing the `link` to, for example, `"br.logit"` or `"br.probit"`. For multinomial treatments, setting `link = "br.logit"` fits a bias-reduced multinomial regression model using `brglm2::brmultinom()`. This can be helpful when regular maximum likelihood models fail to converge, though this may also be a sign of lack of overlap.

# `WeightIt` 0.6.0

* Bug fixes. Functions now work better when used inside other functions (e.g., `lapply`).

* Behavior of `weightit()` in the presence of non-`NULL` `focal` has changed. When `focal` is specified, `estimand` is assumed to be `ATT`. Previously, `focal` would be ignored unless `estimand = "ATT"`.

* Processing of `estimand` and `focal` is improved. Functions are smarter about guessing which group is the focal group when one isn't specified, especially with non-numeric treatments. `focal` can now be used with `estimand = "ATC"` to indicate which group is the control group, so `"ATC"` and `"ATT"` now function more similarly. 

* Added function `get_w_from_ps()` to transform propensity scores into weights (instead of having to go through `weightit()`).

* Added functions `as.weightit()` and `as.weightitMSM()` to convert weights and treatments and other components into `weightit` objects so that `summary.weightit()` can be used on them.

* Updated documentation to describe how missing data in the covariates is handled. Some bugs related to missing data have been fixed as well, thanks to Yong Hao Pua.

* `ps.cont()` had the "z-transformed correlation" options removed to simplify output. This function and its supporting functions will be deprecated as soon as the new version of `twang` is released.

* When using `method = "ps"` or `method = "super"` with continuous treatments, setting `use.kernel = TRUE` and `plot = TRUE`, the plot is now made with `ggplot2` rather than the base R plots.

* Added `plot.summary.weightit()` to plot the distribution of weights (a feature also in `optweight`).

* Removed dependency on `cobalt` temporarily, which means the examples and vignette won't run. 

* Added `ggplot2` to Imports.

# `WeightIt` 0.5.1

* Fixed a bug when using the `ps` argument in `weightit()`.

* Fixed a bug when setting `include.obj = TRUE` in `weightitMSM()`.

* Added warnings for using certain methods with longitudinal treatments as they are not validated and may lead to incorrect inferences.

# `WeightIt` 0.5.0

* Added `super` method to estimate propensity scores using the `SuperLearner` package.

* Added `optweight` method to estimate weights using optimization (but you should probably just use the `optweight` package).

* `weightit()` now uses the correct formula to estimate weights for the ATO with multinomial treatments as described by Li & Li (2018).

* Added `include.obj` option in `weightit()` and `weightitMSM()` to include the fitted object in the output object for inspection. For example, with `method = "ps"`, the `glm` object containing the propensity score model will be included in the output.

* Rearranged the help pages. Each method now has its own documentation page, linked from the `weightit` help page.

* Propensity scores are now included in the output for binary treatments with `gbm` and `cbps` methods. Thanks to @Blanch-Font for the suggestion.

* Other bug fixes and minor changes.

# `WeightIt` 0.4.0

* Added `trim()` function to trim weights.

* Added `ps.cont()` function, which estimates generalized propensity score weights for continuous treatments using generalized boosted modeling, as in `twang`. This function uses the same syntax as `ps()` in `twang`, and can also be accessed using `weightit()` with `method = "gbm"`. Support functions were added to make it compatible with `twang` functions for assessing balance (e.g., `summary`, `bal.table`, `plot`). Thanks to Donna Coffman for enlightening me about this method and providing the code to implement it.

* The input formula is now much more forgiving, allowing objects in the environment to be included. The `data` argument to `weightit()` is now optional. To simplify things, the output object no longer contains a `data` field.

* Under-the-hood changes to facilitate adding new features and debugging. Some aspects of the output objects have been slightly changed, but it shouldn't affect use for most users.

* Fixed a bug where variables would be thrown out when `method = "ebal"`.

# `WeightIt` 0.3.2

* Added new `moments` and `int` options for some `weightit()` methods to easily specify moments and interactions of covariates.

* Fixed bug when using objects not in the data set in `weightit()`. Behavior has changed to include transformed covariates entered in formula in `weightit()` output.

* Fixed bug resulting from potential collinearity when using `ebal` or `ebcw`.

* Added a vignette.

# `WeightIt` 0.3.1

* Edits to code and help files to protect against missing `CBPS` package.

* Corrected sampling weights functionality so they work correctly. Also expanded sampling weights to be able to be used with all methods, including those that do not natively allow for sampling weights (e.g., `ATE`).

* Minor bug fixes and spelling corrections.

# `WeightIt` 0.3.0

* Added `weightitMSM()` function (and supporting `print()` and `summary()` functions) to estimate weights for marginal structural models with time-varying treatments and covariates.

* Fixed some bugs, including when using CBPS with continuous treatments, and when using `focal` incorrectly.

# `WeightIt` 0.2.0

* Added `method = "sbw"` for stable balancing weights (now removed and replaced with `method = "optweight"`)

* Allowed for estimation of multinomial propensity scores using multiple binary regressions if `mlogit` is not installed

* Allowed for estimation of multinomial CBPS using multiple binary CBPS for more than 4 groups

* Added README and NEWS

# `WeightIt` 0.1.0

* First version!
