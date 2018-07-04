WeightIt News and Updates
======

Version 0.4.0

* Added `trim()` function to trim weights.

* Added `ps.cont()` function, which estimates generalized propensity score weights for continuous treatments using generalized boosted modeling, as in `twang`. This function uses the same syntax as `ps()` in `twang`, and can also be accessed using `weightit()` with `method = "gbm"`. Support functions were added to make it compatible with `twang` functions for assessing balance (e.g., `summary`, `bal.table`, `plot`). Thanks to Donna Coffman for enlightening me about this method and providing the code to implement it.

* The input formula is now much more forgiving, allowing objects in the environment to be included. The `data` argument to `weightit()` is now optional. To simplify things, the output object no longer contains a `data` field.

* Under-the-hood changes to facilitate adding new features and debugging. Some aspects of the output objects have been slightly changed, but it shouldn't affect use for most users.

* Fixed a bug where variables would be thrown out when `method = "ebal"`.

* Added support for sampling weights with stable balancing weighting and empirical balancing calibration weighting.

Version 0.3.2

* Added new `moments` and `int` options for some `weightit()` methods to easily specify moments and interactions of covariates.

* Fixed bug when using objects not in the data set in `weightit()`. Behavior has changed to include transformed covariates entered in formula in `weightit()` output.

* Fixed bug resulting from potentially colinearity when using `ebal` or `ebcw`.

* Added a vignette.

Version 0.3.1

* Edits to code and help files to protect against missing `CBPS` package.

* Corrected sampling weights functionality so they work correctly. Also expanded sampling weights to be able to be used with all methods, including those that do not natively allow for sampling weights (e.g., `sbw` and `ATE`)

* Minor bug fixes and spelling corrections.

Version 0.3.0

* Added `weightitMSM()` function (and supporting `print()` and `summary()` functions) to estimate weights for marginal structural models with time-varying treatments and covariates.

* Fixed some bugs, including when using CBPS with continuous treatments, and when using `focal` incorrectly.

Version 0.2.0

* Added `method = "sbw"` for stable balancing weights

* Allowed for estimation of multinomial propensity scores using multiple binary regressions if `mlogit` is not installed

* Allowed for estimation of multinomial CBPS using multiple binary CBPS for more than 4 groups

* Added README and NEWS

Version 0.1.0

* First version!
