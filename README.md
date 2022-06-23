
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WeightIt: Weighting for Covariate Balance in Observational Studies <img src="man/figures/logo.png" align="right" width="150"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/WeightIt?color=00622B)](https://CRAN.R-project.org/package=WeightIt)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/WeightIt?color=00622B)](https://cran.r-project.org/package=WeightIt)
<!-- badges: end -->

------------------------------------------------------------------------

### Overview

`WeightIt` is a one-stop package to generate balancing weights for point
and longitudinal treatments in observational studies. Contained within
`WeightIt` are methods that call on other R packages to estimate
weights. The value of `WeightIt` is in its unified and familiar syntax
used to generate the weights, as each of these other packages have their
own, often challenging to navigate, syntax. `WeightIt` extends the
capabilities of these packages to generate weights used to estimate the
ATE, ATT, ATC, and other estimands for binary or multinomial treatments,
and treatment effects for continuous treatments when available. In these
ways, `WeightIt` does for weighting what `MatchIt` has done for
matching, and `MatchIt` users will find the syntax familiar.

For a complete vignette, see the
[website](https://ngreifer.github.io/WeightIt/articles/WeightIt.html)
for `WeightIt` or `vignette("WeightIt")`.

To install and load `WeightIt`, use the code below:

``` r
#CRAN version
install.packages("WeightIt")

#Development version
remotes::install_github("ngreifer/WeightIt")

library("WeightIt")
```

The workhorse function of `WeightIt` is `weightit()`, which generates
weights from a given formula and data input according to methods and
other parameters specified by the user. Below is an example of the use
of `weightit()` to generate propensity score weights for estimating the
ATE:

``` r
data("lalonde", package = "cobalt")

W <- weightit(treat ~ age + educ + nodegree + 
                married + race + re74 + re75, 
              data = lalonde, method = "ps", 
              estimand = "ATE")
W
```

    A weightit object
     - method: "ps" (propensity score weighting)
     - number of obs.: 614
     - sampling weights: none
     - treatment: 2-category
     - estimand: ATE
     - covariates: age, educ, nodegree, married, race, re74, re75

Evaluating weights has two components: evaluating the covariate balance
produced by the weights, and evaluating whether the weights will allow
for sufficient precision in the eventual effect estimate. For the first
goal, functions in the `cobalt` package, which are fully compatible with
`WeightIt`, can be used, as demonstrated below:

``` r
library("cobalt")

bal.tab(W, un = TRUE)
```

    Call
     weightit(formula = treat ~ age + educ + nodegree + married + 
        race + re74 + re75, data = lalonde, method = "ps", estimand = "ATE")

    Balance Measures
                    Type Diff.Un Diff.Adj
    prop.score  Distance  1.7569   0.1360
    age          Contin. -0.2419  -0.1676
    educ         Contin.  0.0448   0.1296
    nodegree      Binary  0.1114  -0.0547
    married       Binary -0.3236  -0.0944
    race_black    Binary  0.6404   0.0499
    race_hispan   Binary -0.0827   0.0047
    race_white    Binary -0.5577  -0.0546
    re74         Contin. -0.5958  -0.2740
    re75         Contin. -0.2870  -0.1579

    Effective sample sizes
               Control Treated
    Unadjusted  429.    185.  
    Adjusted    329.01   58.33

For the second goal, qualities of the distributions of weights can be
assessed using `summary()`, as demonstrated below.

``` r
summary(W)
```

                     Summary of weights

    - Weight ranges:

               Min                                   Max
    treated 1.1721 |---------------------------| 40.0773
    control 1.0092 |-|                            4.7432

    - Units with 5 most extreme weights by group:
                                                    
                  68     116      10     137     124
     treated 13.5451 15.9884 23.2967 23.3891 40.0773
                 597     573     381     411     303
     control  4.0301  4.0592  4.2397  4.5231  4.7432

    - Weight statistics:

            Coef of Var   MAD Entropy # Zeros
    treated       1.478 0.807   0.534       0
    control       0.552 0.391   0.118       0

    - Effective Sample Sizes:

               Control Treated
    Unweighted  429.    185.  
    Weighted    329.01   58.33

Desirable qualities include small coefficients of variation close to 0
and large effective sample sizes.

The table below contains the available methods in `WeightIt` for
estimating weights for binary, multinomial, and continuous treatments
using various methods and functions from various packages. See
`vignette("installing-packages")` for information on how to install
these packages.

| Treatment type  | Method (`method =`)                                 | Package        |
|-----------------|-----------------------------------------------------|----------------|
| **Binary**      | Binary regression PS (`"ps"`)                       | various        |
| \-              | Generalized boosted modeling PS (`"gbm"`)           | `gbm`          |
| \-              | Covariate Balancing PS (`"cbps"`)                   | `CBPS`         |
| \-              | Non-Parametric Covariate Balancing PS (`"npcbps"`)  | `CBPS`         |
| \-              | Entropy Balancing (`"ebal"`)                        | \-             |
| \-              | Empirical Balancing Calibration Weights (`"ebcw"`)  | `ATE`          |
| \-              | Optimization-Based Weights (`"optweight"`)          | `optweight`    |
| \-              | SuperLearner PS (`"super"`)                         | `SuperLearner` |
| \-              | Bayesian additive regression trees PS (`"bart"`)    | `dbarts`       |
| \-              | Energy Balancing (`"energy"`)                       | \-             |
| **Multinomial** | Multinomial regression PS (`"ps"`)                  | various        |
| \-              | Generalized boosted modeling PS (`"gbm"`)           | `gbm`          |
| \-              | Covariate Balancing PS (`"cbps"`)                   | `CBPS`         |
| \-              | Non-Parametric Covariate Balancing PS (`"npcbps"`)  | `CBPS`         |
| \-              | Entropy Balancing (`"ebal"`)                        | \-             |
| \-              | Empirical Balancing Calibration Weights (`"ebcw"`)  | `ATE`          |
| \-              | Optimization-Based Weights (`"optweight"`)          | `optweight`    |
| \-              | SuperLearner PS (`"super"`)                         | `SuperLearner` |
| \-              | Bayesian additive regression trees PS (`"bart"`)    | `dbarts`       |
| \-              | Energy Balancing (`"energy"`)                       | \-             |
| **Continuous**  | Generalized linear model GPS (`"ps"`)               | \-             |
| \-              | Generalized boosted modeling GPS (`"gbm"`)          | `gbm`          |
| \-              | Covariate Balancing GPS (`"cbps"`)                  | `CBPS`         |
| \-              | Non-Parametric Covariate Balancing GPS (`"npcbps"`) | `CBPS`         |
| \-              | Entropy Balancing (`"ebal"`)                        | \-             |
| \-              | Optimization-Based Weights (`"optweight"`)          | `optweight`    |
| \-              | SuperLearner GPS (`"super"`)                        | `SuperLearner` |
| \-              | Bayesian additive regression trees GPS (`"bart"`)   | `dbarts`       |

In addition, `WeightIt` implements the subgroup balancing propensity
score using the function `sbps()`. Several other tools and utilities are
available.

Please submit bug reports or other issues to
<https://github.com/ngreifer/WeightIt/issues>. If you would like to see
your package or method integrated into `WeightIt`, or for any other
questions or comments about `WeightIt`, please contact the author. Fan
mail is greatly appreciated.
