WeightIt News and Updates
======

# WeightIt 0.10.2

* Fixed a bug where treatment values were accidentally switched for some methods.

# WeightIt 0.10.1

* With `method = "gbm"`, added the ability to tune hyperparameters like `interaction.depth` and `distribution` using the same critera as is used to select the optimal tree. A summary of the tuning results is included in `info` in the `weightit` output object.

* Fixed a bug where `moments` and `int` were ignored unless both were specified.

* Effective sample sizes now print only up to two digits (believe me, you don't need three) and print more cleanly with whole numbers.

* Fixed a bug when using `by`, thanks to @frankpopham. (#11)

* Fixed a bug when using `weightitMSM` with methods that process `int` and `moments` (though you probably shouldn't use them anyway). Thanks to Sven Reiger.

* Fixed a bug when using `method = "npcbps"` where weights could be excessively small and mistaken for all being them. The weights now sum to the number of units.

# WeightIt 0.10.0

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

# WeightIt 0.9.0

* Added support for entropy balancing (`method = "ebal"`) for continuous treatments as described by Tübbicke (2020). Relies on hand-written code contributed by Stefan Tübbicke rather than another R package. Sampling weights and base weights are both supported as they are with binary and multi-category treatments.

* Added support for Balance SuperLearner as described by Pirracchio and Carone (2018) with `method = "super"`. Rather than using NNLS to choose the optimal combination of predictions, you can now optimize balance. To do so, set `SL.method = "method.balance"`. You will need to set an argument to `stop.method`, which works identically to how it does for `method = "gbm"`. For example, for `stop.method = "es.max"`, the predicted values given will be the combination of predicted values that minimizes the largest absolute standardized mean difference of the covariates in the sample weighted using the predicted values as propensity scores.

* Changed some of the statistics displayed when using `summary()`: the weight ratio is gone (because weights can be 0, which is not problematic but would explode the ratio), and the mean absolute deviation and entropy of the weights are now present.

* Added `crayon` for prettier printing of `summary()` output.

# WeightIt 0.8.0

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

# WeightIt 0.7.1

* Fixed bug when using `weightit()` inside another function that passed a `by` argument explicitly. Also changed the syntax for `by`; it must now either be a string (which was always possible) or a one-sided formula with the stratifying variable on the right-hand side. To use a variable that is not in `data`, you must use the formula interface. 

* Fixed bug when trying to use `ps` with `by` in `weightit()`.

# WeightIt 0.7.0

* Added new `sbps()` function for estimating subgroup balancing propensity score weights, including both the standard method and a new smooth version.

* Setting `method = "gbm"` and `method = "twang"` will now do two different things. `method = "gbm"` uses `gbm` and `cobalt` functions to estimate the weights and is much faster, while `method = "twang"` uses `twang` functions to estimate the weights. The results are similar between the two methods. Prior to this version, `method = "gbm"` and `method = "twang"` both did what `method = "twang"` does now. 

* Bug fixes when `stabilize = TRUE`, thanks to @ulriksartipy and Sven Rieger.

* Fixes for using `base.weight` argument with `method = "ebal"`. Now the supplied vector should have a length equal to the number of units in the dataset (in contrast to its use in `ebalance`, which requires a length equal to the number of control units).

* Restored dependency on `cobalt` for examples and vignette.

* When `method = "ps"` and the treatment is ordered (i.e., ordinal), `MASS::polr()` is used to fit an ordinal regression. Make the treatment un-ordered to to use multinomial regression instead.

* Added support for using bias-reduced fitting functions when `method = "ps"` as provided by the `brglm2` package. These can be accessed by changing the `link` to, for example, `"br.logit"` or `"br.probit"`. For multinomial treatments, setting `link = "br.logit"` fits a bias-reduced multinomial regression model using `brglm2::brmultinom()`. This can be helpful when regular maximum likelihood models fail to converge, though this may also be a sign of lack of overlap.

# WeightIt 0.6.0

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

# WeightIt 0.5.1

* Fixed a bug when using the `ps` argument in `weightit()`.

* Fixed a bug when setting `include.obj = TRUE` in `weightitMSM()`.

* Added warnings for using certain methods with longitudinal treatments as they are not validated and may lead to incorrect inferences.

# WeightIt 0.5.0

* Added `super` method to estimate propensity scores using the `SuperLearner` package.

* Added `optweight` method to estimate weights using optimization (but you should probably just use the `optweight` package).

* `weightit()` now uses the correct formula to estimate weights for the ATO with multinomial treatments as described by Li & Li (2018).

* Added `include.obj` option in `weightit()` and `weightitMSM()` to include the fitted object in the output object for inspection. For example, with `method = "ps"`, the `glm` object containing the propensity score model will be included in the output.

* Rearranged the help pages. Each method now has its own documentation page, linked from the `weightit` help page.

* Propensity scores are now included in the output for binary treatments with `gbm` and `cbps` methods. Thanks to @Blanch-Font for the suggestion.

* Other bug fixes and minor changes.

# WeightIt 0.4.0

* Added `trim()` function to trim weights.

* Added `ps.cont()` function, which estimates generalized propensity score weights for continuous treatments using generalized boosted modeling, as in `twang`. This function uses the same syntax as `ps()` in `twang`, and can also be accessed using `weightit()` with `method = "gbm"`. Support functions were added to make it compatible with `twang` functions for assessing balance (e.g., `summary`, `bal.table`, `plot`). Thanks to Donna Coffman for enlightening me about this method and providing the code to implement it.

* The input formula is now much more forgiving, allowing objects in the environment to be included. The `data` argument to `weightit()` is now optional. To simplify things, the output object no longer contains a `data` field.

* Under-the-hood changes to facilitate adding new features and debugging. Some aspects of the output objects have been slightly changed, but it shouldn't affect use for most users.

* Fixed a bug where variables would be thrown out when `method = "ebal"`.

# WeightIt 0.3.2

* Added new `moments` and `int` options for some `weightit()` methods to easily specify moments and interactions of covariates.

* Fixed bug when using objects not in the data set in `weightit()`. Behavior has changed to include transformed covariates entered in formula in `weightit()` output.

* Fixed bug resulting from potentially colinearity when using `ebal` or `ebcw`.

* Added a vignette.

# WeightIt 0.3.1

* Edits to code and help files to protect against missing `CBPS` package.

* Corrected sampling weights functionality so they work correctly. Also expanded sampling weights to be able to be used with all methods, including those that do not natively allow for sampling weights (e.g., `ATE`).

* Minor bug fixes and spelling corrections.

# WeightIt 0.3.0

* Added `weightitMSM()` function (and supporting `print()` and `summary()` functions) to estimate weights for marginal structural models with time-varying treatments and covariates.

* Fixed some bugs, including when using CBPS with continuous treatments, and when using `focal` incorrectly.

# WeightIt 0.2.0

* Added `method = "sbw"` for stable balancing weights (now removed and replaced with `method = "optweight"`)

* Allowed for estimation of multinomial propensity scores using multiple binary regressions if `mlogit` is not installed

* Allowed for estimation of multinomial CBPS using multiple binary CBPS for more than 4 groups

* Added README and NEWS

# WeightIt 0.1.0

* First version!
