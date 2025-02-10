
<!-- README.md is generated from README.Rmd. Please edit that file -->

# WeightIt: Weighting for Covariate Balance in Observational Studies <img src="man/figures/logo.png" align="right" width="150"/>

[![CRAN
status](https://www.r-pkg.org/badges/version/WeightIt?color=00622B)](https://CRAN.R-project.org/package=WeightIt)
[![CRAN_Downloads_Badge](https://cranlogs.r-pkg.org/badges/WeightIt?color=00622B)](https://cran.r-project.org/package=WeightIt)

### Overview

*WeightIt* is a one-stop package to generate balancing weights for point
and longitudinal treatments in observational studies. Support is
included for binary, multi-category, and continuous treatments, a
variety of estimands including the ATE, ATT, ATC, ATO, and others, and
support for a wide variety of weighting methods, including those that
rely on parametric modeling, machine learning, or optimization.
*WeightIt* also provides functionality for fitting regression models in
weighted samples that account for estimation of the weights in
quantifying uncertainty. *WeightIt* uses a familiar formula interface
and is meant to complement `MatchIt` as a package that provides a
unified interface to basic and advanced weighting methods.

For a complete vignette, see the
[website](https://ngreifer.github.io/WeightIt/articles/WeightIt.html)
for *WeightIt* or `vignette("WeightIt")`.

To install and load *WeightIt* , use the code below:

``` r
#CRAN version
pak::pkg_install("WeightIt")

#Development version
pak::pkg_install("ngreifer/WeightIt")

library("WeightIt")
```

The workhorse function of *WeightIt* is `weightit()`, which generates
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
goal, functions in the *cobalt* package, which are fully compatible with
*WeightIt*, can be used, as demonstrated below:

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

Desirable qualities include large effective sample sizes, which imply
low variability in the weights (and therefore increased precision in
estimating the treatment effect).

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
    #>             Estimate Std. Error z value Pr(>|z|)  2.5 % 97.5 %    
    #> (Intercept)   5135.1      583.8   8.797   <1e-06 3990.9 6279.2 ***
    #> treat         1214.1      798.2   1.521    0.128 -350.3 2778.4    
    #> Standard error: HC0 robust (adjusted for estimation of weights)

The tables below contains the available methods in *WeightIt* for
estimating weights for binary, multi-category, and continuous
treatments. Many of these methods do not require any other package to
use; see `vignette("installing-packages")` for information on how to
install packages that are used.

#### Binary Treatments

| Method | `method` |
|----|----|
| Binary regression PS | [`"glm"`](https://ngreifer.github.io/WeightIt/reference/method_glm.html) |
| Generalized boosted modeling PS | [`"gbm"`](https://ngreifer.github.io/WeightIt/reference/method_gbm.html) |
| Covariate balancing PS | [`"cbps"`](https://ngreifer.github.io/WeightIt/reference/method_cbps.html) |
| Non-Parametric covariate balancing PS | [`"npcbps"`](https://ngreifer.github.io/WeightIt/reference/method_npcbps.html) |
| Entropy balancing | [`"ebal"`](https://ngreifer.github.io/WeightIt/reference/method_ebal.html) |
| Inverse probability tilting | [`"ipt"`](https://ngreifer.github.io/WeightIt/reference/method_ipt.html) |
| Stable balancing weights | [`"optweight"`](https://ngreifer.github.io/WeightIt/reference/method_optweight.html) |
| SuperLearner PS | [`"super"`](https://ngreifer.github.io/WeightIt/reference/method_super.html) |
| Bayesian additive regression trees PS | [`"bart"`](https://ngreifer.github.io/WeightIt/reference/method_bart.html) |
| Energy balancing | [`"energy"`](https://ngreifer.github.io/WeightIt/reference/method_energy.html) |

#### Multi-Category Treatments

| Method | `method` |
|----|----|
| Multinomial regression PS | [`"glm"`](https://ngreifer.github.io/WeightIt/reference/method_glm.html) |
| Generalized boosted modeling PS | [`"gbm"`](https://ngreifer.github.io/WeightIt/reference/method_gbm.html) |
| Covariate balancing PS | [`"cbps"`](https://ngreifer.github.io/WeightIt/reference/method_cbps.html) |
| Non-parametric covariate balancing PS | [`"npcbps"`](https://ngreifer.github.io/WeightIt/reference/method_npcbps.html) |
| Entropy balancing | [`"ebal"`](https://ngreifer.github.io/WeightIt/reference/method_ebal.html) |
| Inverse probability tilting | [`"ipt"`](https://ngreifer.github.io/WeightIt/reference/method_ipt.html) |
| Stable balancing weights | [`"optweight"`](https://ngreifer.github.io/WeightIt/reference/method_optweight.html) |
| SuperLearner PS | [`"super"`](https://ngreifer.github.io/WeightIt/reference/method_super.html) |
| Bayesian additive regression trees PS | [`"bart"`](https://ngreifer.github.io/WeightIt/reference/method_bart.html) |
| Energy balancing | [`"energy"`](https://ngreifer.github.io/WeightIt/reference/method_energy.html) |

#### Continuous Treatments

| Method | `method` |
|----|----|
| Generalized linear model GPS | [`"glm"`](https://ngreifer.github.io/WeightIt/reference/method_glm.html) |
| Generalized boosted modeling GPS | [`"gbm"`](https://ngreifer.github.io/WeightIt/reference/method_gbm.html) |
| Covariate balancing GPS | [`"cbps"`](https://ngreifer.github.io/WeightIt/reference/method_cbps.html) |
| Non-parametric covariate balancing GPS | [`"npcbps"`](https://ngreifer.github.io/WeightIt/reference/method_npcbps.html) |
| Entropy balancing | [`"ebal"`](https://ngreifer.github.io/WeightIt/reference/method_ebal.html) |
| Stable balancing weights | [`"optweight"`](https://ngreifer.github.io/WeightIt/reference/method_optweight.html) |
| SuperLearner GPS | [`"super"`](https://ngreifer.github.io/WeightIt/reference/method_super.html) |
| Bayesian additive regression trees GPS | [`"bart"`](https://ngreifer.github.io/WeightIt/reference/method_bart.html) |
| Distance covariance optimal weighting | [`"energy"`](https://ngreifer.github.io/WeightIt/reference/method_energy.html) |

In addition, *WeightIt* implements the subgroup balancing propensity
score using the function `sbps()`. Several other tools and utilities are
available, including `trim()` to trim or truncate weights, `calibrate()`
to calibrate propensity scores, `get_w_from_ps()` to compute weights
from propensity scores.

*WeightIt* provides functions to fit weighted models that account for
the uncertainty in estimating the weights. These include
`glm_weightit()` for fitting generalized linear models,
`ordinal_weightit()` for ordinal regression models,
`multinom_weightit()` for multinomial regression models, and
`coxph_weightit()` for Cox proportional hazards models. Several methods
are available for computing the parameter variances, including
asymptotically correct M-estimation-based variances, robust variances
that treat the weights as fixed, and traditional and fractional weighted
bootstrap variances. Clustered variances are supported. See
`vignette("estimating-effects")` for information on how to use these
after weighting to estimate treatment effects.

Please submit bug reports, questions, comments, or other issues to
<https://github.com/ngreifer/WeightIt/issues>. If you would like to see
your package or method integrated into *WeightIt*, please contact the
author. Fan mail is greatly appreciated.
