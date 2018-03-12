WeightIt News and Updates
======

Version 0.3.2

* Added new `moments` and `int` options for some `weightit()` methods to easily specify moments and interactions of covariates.

* Fixed bug when using objects not in the data set in `weightit()`. Behavior has changed to include transformed covariates entered in formula in `weightit()` output.

* Fixed bug resulting from potentially colinearity when using `ebal` or `ebcw`.

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
