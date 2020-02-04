
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WeightIt

[![CRAN\_Status\_Badge](http://r-pkg.org/badges/version-last-release/WeightIt?color=0047ab)](https://cran.r-project.org/package=WeightIt)
[![CRAN\_Downloads\_Badge](http://cranlogs.r-pkg.org/badges/WeightIt?color=0047ab)](https://cran.r-project.org/package=WeightIt)

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

For a complete vignette, see the [CRAN
page](https://CRAN.R-project.org/package=WeightIt) for `WeightIt`.

To install and load `WeightIt`, use the code below:

``` r
install.packages("WeightIt")  #CRAN version
devtools::install_github("ngreifer/WeightIt")  #Development version
library("WeightIt")
```

The workhorse function of `WeightIt` is `weightit()`, which generates
weights from a given formula and data input according to methods and
other parameters specified by the user. Below is an example of the use
of `weightit()` to generate propensity score weights for estimating the
ATE:

``` r
data("lalonde", package = "cobalt")
W <- weightit(treat ~ age + educ + nodegree + married + race + re74 + re75, data = lalonde, 
    method = "ps", estimand = "ATE")
print(W)
```

    A weightit object
     - method: "ps" (propensity score weighting)
     - number of obs.: 614
     - sampling weights: none
     - treatment: 2-category
     - estimand: ATE
     - covariates: age, educ, nodegree, married, race, re74, re75

Evaluating weights has two components: evaluating the covariate balance
produces by the weights, and evaluating whether the weights will allow
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
    Unadjusted 429.000 185.000
    Adjusted   329.008  58.327

For the second goal, qualities of the distributions of weights can be
assessed using `summary()`, as demonstrated below.

``` r
summary(W)
```

    Summary of weights:
    
    - Weight ranges:
               Min                                   Max
    treated 1.1721 |---------------------------| 40.0773
    control 1.0092 |-|                            4.7432
    
    - Units with 5 greatest weights by group:
                                                    
                 137     124     116      68      10
     treated 13.5451 15.9884 23.2967 23.3891 40.0773
                 597     573     411     381     303
     control  4.0301  4.0592  4.2397  4.5231  4.7432
    
              Ratio Coef of Var
    treated 34.1921      1.4777
    control  4.7002      0.5519
    overall 39.7134      1.3709
    
    - Effective Sample Sizes:
               Control Treated
    Unweighted 429.000 185.000
    Weighted   329.008  58.327

Desirable qualities include ratios close to 1, coefficients of variation
close to 0, and large effective sample sizes.

The table below contains the available methods in `WeightIt` for
estimating weights for binary, multinomial, and continuous treatments
using various methods and functions from various
packages.

| Treatment type  | Method (`method =`)                                                | Function                | Package          |
| --------------- | ------------------------------------------------------------------ | ----------------------- | ---------------- |
| **Binary**      | Binary regression PS (`"ps"`)                                      | `glm()`                 | `base`           |
| \-              | Generalized boosted modeling PS (`"gbm"`/`"twang"`)                | `gbm.fit()`/`ps()`      | `gbm`/`twang`    |
| \-              | Covariate Balancing PS (`"cbps"`)                                  | `CBPS()`                | `CBPS`           |
| \-              | Non-Parametric Covariate Balancing PS (`"npcbps"`)                 | `npCBPS()`              | `CBPS`           |
| \-              | Entropy Balancing (`"ebal"`)                                       | `ebalance()`            | `ebal`           |
| \-              | Empirical Balancing Calibration Weights (`"ebcw"`)                 | `ATE()`                 | `ATE`            |
| \-              | Optimization-Based Weights (`"optweight"`)                         | `optweight()`           | `optweight`      |
| \-              | SuperLearner PS (`"super"`)                                        | `SuperLearner()`        | `SuperLearner`   |
| **Multinomial** | Multiple binary regression PS (`"ps"`)                             | `glm()`                 | `base`           |
| \-              | Multinomial regression PS (`"ps"`)                                 | `mlogit()`              | `mlogit`         |
| \-              | Bayesian multinomial regression PS (`"ps", link = "bayes.probit"`) | `MNP()`                 | `MNP`            |
| \-              | Generalized boosted modeling PS (`"gbm"`/`"twang"`)                | `gbm.fit()`/`mnps()`    | `gbm`/`twang`    |
| \-              | Covariate Balancing PS (`"cbps"`)                                  | `CBPS()`                | `CBPS`           |
| \-              | Non-Parametric Covariate Balancing PS (`"npcbps"`)                 | `npCBPS()`              | `CBPS`           |
| \-              | Entropy Balancing (`"ebal"`)                                       | `ebalance()`            | `ebal`           |
| \-              | Empirical Balancing Calibration Weights (`"ebcw"`)                 | `ATE()`                 | `ATE`            |
| \-              | Optimization-Based Weights (`"optweight"`)                         | `optweight()`           | `optweight`      |
| \-              | SuperLearner PS (`"super"`)                                        | `SuperLearner()`        | `SuperLearner`   |
| **Continuous**  | Generalized linear model PS (`"ps"`)                               | `glm()`                 | `base`           |
| \-              | Generalized boosted modeling PS (`"gbm"`/`"twang"`)                | `gbm.fit()`/`ps.cont()` | `gbm`/`WeightIt` |
| \-              | Covariate Balancing PS (`"cbps"`)                                  | `CBPS()`                | `CBPS`           |
| \-              | Non-Parametric Covariate Balancing PS (`"npcbps"`)                 | `npCBPS()`              | `CBPS`           |
| \-              | Entropy Balancing (`"ebal"`)                                       | `optim()`               | `base`           |
| \-              | Optimization-Based Weights (`"optweight"`)                         | `optweight()`           | `optweight`      |
| \-              | SuperLearner PS (`"super"`)                                        | `SuperLearner()`        | `SuperLearner`   |

In addition, `WeightIt` implements the subgroup balancing propensity
score using the function `sbps()`. Several other tools and utilities are
available.

Please submit bug reports or other issues to
<https://github.com/ngreifer/WeightIt/issues>. If you would like to see
your package or method integrated into `WeightIt`, or for any other
questions or comments about `WeightIt`, please contact Noah Greifer at
<noah.greifer@gmail.com>. Fan mail is greatly appreciated.
