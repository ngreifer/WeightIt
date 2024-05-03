
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
and longitudinal treatments in observational studies. Support is
included for binary, multi-category, and continuous treatments, a
variety of estimands including the ATE, ATT, ATC, ATO, and others, and
support for a wide variety of weighting methods, including those that
rely on parametric modeling, machine learning, or optimization.
`WeightIt` also provides functionality for fitting regression models in
weighted samples that account for estimation of the weights in
quantifying uncertainty. `WeightIt` uses a familiar formula interface
and is meant to complement `MatchIt` as a package that provides a
unified interface to basic and advanced weighting methods.

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
ATT:

``` r
data("lalonde", package = "cobalt")

W <- weightit(treat ~ age + educ + nodegree + 
                married + race + re74 + re75, 
              data = lalonde, method = "glm", 
              estimand = "ATT")
W
```

    #> A weightit object
    #>  - method: "glm" (propensity score weighting with GLM)
    #>  - number of obs.: 614
    #>  - sampling weights: none
    #>  - treatment: 2-category
    #>  - estimand: ATT (focal: 1)
    #>  - covariates: age, educ, nodegree, married, race, re74, re75

Evaluating weights has two components: evaluating the covariate balance
produced by the weights, and evaluating whether the weights will allow
for sufficient precision in the eventual effect estimate. For the first
goal, functions in the `cobalt` package, which are fully compatible with
`WeightIt`, can be used, as demonstrated below:

``` r
library("cobalt")

bal.tab(W, un = TRUE)
```

    #> Balance Measures
    #>                 Type Diff.Un Diff.Adj
    #> prop.score  Distance  1.7941  -0.0205
    #> age          Contin. -0.3094   0.1188
    #> educ         Contin.  0.0550  -0.0284
    #> nodegree      Binary  0.1114   0.0184
    #> married       Binary -0.3236   0.0186
    #> race_black    Binary  0.6404  -0.0022
    #> race_hispan   Binary -0.0827   0.0002
    #> race_white    Binary -0.5577   0.0021
    #> re74         Contin. -0.7211  -0.0021
    #> re75         Contin. -0.2903   0.0110
    #> 
    #> Effective sample sizes
    #>            Control Treated
    #> Unadjusted  429.       185
    #> Adjusted     99.82     185

For the second goal, qualities of the distributions of weights can be
assessed using `summary()`, as demonstrated below.

``` r
summary(W)
```

    #>                   Summary of weights
    #> 
    #> - Weight ranges:
    #> 
    #>            Min                                  Max
    #> treated 1.0000         ||                    1.0000
    #> control 0.0092 |---------------------------| 3.7432
    #> 
    #> - Units with the 5 most extreme weights by group:
    #>                                            
    #>               5      4      3      2      1
    #>  treated      1      1      1      1      1
    #>             597    573    381    411    303
    #>  control 3.0301 3.0592 3.2397 3.5231 3.7432
    #> 
    #> - Weight statistics:
    #> 
    #>         Coef of Var   MAD Entropy # Zeros
    #> treated       0.000 0.000   0.000       0
    #> control       1.818 1.289   1.098       0
    #> 
    #> - Effective Sample Sizes:
    #> 
    #>            Control Treated
    #> Unweighted  429.       185
    #> Weighted     99.82     185

Desirable qualities include small coefficients of variation close to 0
and large effective sample sizes.

Finally, we can estimate the effect of the treatment using a weighted
outcome model, accounting for estimation of the weights in the standard
error of the effect estimate:

``` r
fit <- lm_weightit(re78 ~ treat, data = lalonde,
                   weightit = W)

summary(fit, ci = TRUE)
```

    #> 
    #> Call:
    #> lm_weightit(formula = re78 ~ treat, data = lalonde, weightit = W)
    #> 
    #> Coefficients:
    #>             Estimate Std. Error z value  Pr(>|z|)  2.5 % 97.5 %
    #> (Intercept)     5135      583.8   8.797 1.411e-18 3990.9   6279
    #> treat           1214      798.2   1.521 1.282e-01 -350.3   2778
    #> Standard error: HC0 robust (adjusted for estimation of weights)

The table below contains the available methods in `WeightIt` for
estimating weights for binary, multinomial, and continuous treatments
using various methods and functions from various packages. Many of these
methods do not require any other package to use (i.e., those with “-” in
the Package column). See `vignette("installing-packages")` for
information on how to install packages that are used.

| Treatment type     | Method (`method =`)                                 | Package        |
|--------------------|-----------------------------------------------------|----------------|
| **Binary**         | Binary regression PS (`"glm"`)                      | various        |
| \-                 | Generalized boosted modeling PS (`"gbm"`)           | `gbm`          |
| \-                 | Covariate balancing PS (`"cbps"`)                   | \-             |
| \-                 | Non-parametric covariate balancing PS (`"npcbps"`)  | `CBPS`         |
| \-                 | Entropy Balancing (`"ebal"`)                        | \-             |
| \-                 | Inverse probability tilting (`"ipt"`)               | \-             |
| \-                 | Optimization-based Weights (`"optweight"`)          | `optweight`    |
| \-                 | SuperLearner PS (`"super"`)                         | `SuperLearner` |
| \-                 | Bayesian additive regression trees PS (`"bart"`)    | `dbarts`       |
| \-                 | Energy balancing (`"energy"`)                       | \-             |
| **Multi-category** | Multinomial regression PS (`"glm"`)                 | various        |
| \-                 | Generalized boosted modeling PS (`"gbm"`)           | `gbm`          |
| \-                 | Covariate balancing PS (`"cbps"`)                   | \-             |
| \-                 | Non-Parametric covariate balancing PS (`"npcbps"`)  | `CBPS`         |
| \-                 | Entropy balancing (`"ebal"`)                        | \-             |
| \-                 | Inverse probability tilting (`"ipt"`)               | \-             |
| \-                 | Optimization-based weights (`"optweight"`)          | `optweight`    |
| \-                 | SuperLearner PS (`"super"`)                         | `SuperLearner` |
| \-                 | Bayesian additive regression trees PS (`"bart"`)    | `dbarts`       |
| \-                 | Energy balancing (`"energy"`)                       | \-             |
| **Continuous**     | Generalized linear model GPS (`"glm"`)              | \-             |
| \-                 | Generalized boosted modeling GPS (`"gbm"`)          | `gbm`          |
| \-                 | Covariate balancing GPS (`"cbps"`)                  | \-             |
| \-                 | Non-Parametric covariate balancing GPS (`"npcbps"`) | `CBPS`         |
| \-                 | Entropy balancing (`"ebal"`)                        | \-             |
| \-                 | Optimization-based weights (`"optweight"`)          | `optweight`    |
| \-                 | SuperLearner GPS (`"super"`)                        | `SuperLearner` |
| \-                 | Bayesian additive regression trees GPS (`"bart"`)   | `dbarts`       |
| \-                 | Distance covariance optimal weighting (`"energy"`)  | \-             |

In addition, `WeightIt` implements the subgroup balancing propensity
score using the function `sbps()`. Several other tools and utilities are
available, including `trim()` to trim or truncate weights.

Please submit bug reports, questions, comments, or other issues to
<https://github.com/ngreifer/WeightIt/issues>. If you would like to see
your package or method integrated into `WeightIt`, please contact the
author. Fan mail is greatly appreciated.
