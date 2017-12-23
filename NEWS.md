WeightIt News and Updates
======
Version 0.3.0

* Added `weightitMSM()` function (and supporting `print()` and `summary()` functions) to estimate weights for marginal structural models with time-varying treatments and covariates.

* Bug fixes and minor improvements

Version 0.2.1

* Fixed some bugs, including when using CBPS with continuous treatments, and when using `focal` incorrectly.

Version 0.2.0

* Added `method = "sbw"` for stabilized balancing weights

* Allowed for estimation of multinomial propensity scores using multiple binary regressions if `mlogit` is not installed

* Allowed for estimation of multinomial CBPS using multiple binary CBPS for more than 4 groups

* Added README and NEWS

Version 0.1.0

* First version!
