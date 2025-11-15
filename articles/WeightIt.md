# Using WeightIt to Estimate Balancing Weights

## Introduction

*WeightIt* contains several functions for estimating and assessing
balancing weights for observational studies. These weights can be used
to estimate the causal parameters of marginal structural models. I will
not go into the basics of causal inference methods here. For good
introductory articles, see Austin
([2011](#ref-austinIntroductionPropensityScore2011)), Austin and Stuart
([2015](#ref-austinMovingBestPractice2015)), Robins, Hernán, and
Brumback ([2000](#ref-robinsMarginalStructuralModels2000)), or Thoemmes
and Ong ([2016](#ref-thoemmesPrimerInverseProbability2016)).

Typically, the analysis of an observational study might proceed as
follows: identify the covariates for which balance is required; assess
the quality of the data available, including missingness and measurement
error; estimate weights that balance the covariates adequately; and
estimate a treatment effect and corresponding standard error or
confidence interval. This guide will go through these steps for two
observational studies: estimating the causal effect of a point treatment
on an outcome, and estimating the causal parameters of a marginal
structural model with multiple treatment periods. This is not meant to
be a definitive guide, but rather an introduction to the relevant
issues.

## Balancing Weights for a Point Treatment

First we will use the Lalonde dataset to estimate the effect of a point
treatment. We’ll use the version of the data set that comes with the
*cobalt* package, which we will use later on as well. Here, we are
interested in the average treatment effect on the treated (ATT).

``` r
library("cobalt")
data("lalonde", package = "cobalt")
head(lalonde)
```

    ##   treat age educ   race married nodegree re74 re75    re78
    ## 1     1  37   11  black       1        1    0    0  9930.0
    ## 2     1  22    9 hispan       0        1    0    0  3595.9
    ## 3     1  30   12  black       0        0    0    0 24909.5
    ## 4     1  27   11  black       0        1    0    0  7506.1
    ## 5     1  33    8  black       0        1    0    0   289.8
    ## 6     1  22    9  black       0        1    0    0  4056.5

We have our outcome (`re78`), our treatment (`treat`), and the
covariates for which balance is desired (`age`, `educ`, `race`,
`married`, `nodegree`, `re74`, and `re75`). Using *cobalt*, we can
examine the initial imbalance on the covariates:

``` r
bal.tab(treat ~ age + educ + race + married + nodegree +
          re74 + re75,
        data = lalonde,
        estimand = "ATT",
        thresholds = c(m = .05))
```

    ## Balance Measures
    ##                Type Diff.Un      M.Threshold.Un
    ## age         Contin.  -0.309 Not Balanced, >0.05
    ## educ        Contin.   0.055 Not Balanced, >0.05
    ## race_black   Binary   0.640 Not Balanced, >0.05
    ## race_hispan  Binary  -0.083 Not Balanced, >0.05
    ## race_white   Binary  -0.558 Not Balanced, >0.05
    ## married      Binary  -0.324 Not Balanced, >0.05
    ## nodegree     Binary   0.111 Not Balanced, >0.05
    ## re74        Contin.  -0.721 Not Balanced, >0.05
    ## re75        Contin.  -0.290 Not Balanced, >0.05
    ## 
    ## Balance tally for mean differences
    ##                     count
    ## Balanced, <0.05         0
    ## Not Balanced, >0.05     9
    ## 
    ## Variable with the greatest mean difference
    ##  Variable Diff.Un      M.Threshold.Un
    ##      re74  -0.721 Not Balanced, >0.05
    ## 
    ## Sample sizes
    ##     Control Treated
    ## All     429     185

Based on this output, we can see that all variables are imbalanced in
the sense that the standardized mean differences (for continuous
variables) and differences in proportion (for binary variables) are
greater than .05 for all variables. In particular, `re74` and `re75` are
quite imbalanced, which is troubling given that they are likely strong
predictors of the outcome. We will estimate weights using
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
to try to attain balance on these covariates.

First, we’ll start simple, and use inverse probability weights from
propensity scores generated through logistic regression. We need to
supply
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
with the formula for the model, the data set, the estimand (ATT), and
the method of estimation (`"glm"` for generalized linear model
propensity score weights).

``` r
library("WeightIt")
W.out <- weightit(treat ~ age + educ + race + married + nodegree +
                    re74 + re75,
                  data = lalonde,
                  estimand = "ATT",
                  method = "glm")
W.out #print the output
```

    ## A weightit object
    ##  - method: "glm" (propensity score weighting with GLM)
    ##  - number of obs.: 614
    ##  - sampling weights: none
    ##  - treatment: 2-category
    ##  - estimand: ATT (focal: 1)
    ##  - covariates: age, educ, race, married, nodegree, re74, re75

Printing the output of
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
displays a summary of how the weights were estimated. Let’s examine the
quality of the weights using
[`summary()`](https://rdrr.io/r/base/summary.html). Weights with low
variability are desirable because they improve the precision of the
estimator. This variability is presented in several ways, but the most
important is the effective sample size (ESS) computed from the weights,
which we hope is as close to the original sample size as possible. What
constitutes a “large enough” ESS is mostly relative, though, and must be
considered with respect other constraints, including covariate balance.

``` r
summary(W.out)
```

    ##                   Summary of weights
    ## 
    ## - Weight ranges:
    ## 
    ##           Min                                 Max
    ## treated 1.            ||                    1.   
    ## control 0.009 |---------------------------| 3.743
    ## 
    ## - Units with the 5 most extreme weights by group:
    ##                                     
    ##             1     2    3     4     5
    ##  treated    1     1    1     1     1
    ##           412   388  226   196   118
    ##  control 3.03 3.059 3.24 3.523 3.743
    ## 
    ## - Weight statistics:
    ## 
    ##         Coef of Var   MAD Entropy # Zeros
    ## treated       0.000 0.000   0.000       0
    ## control       1.818 1.289   1.098       0
    ## 
    ## - Effective Sample Sizes:
    ## 
    ##            Control Treated
    ## Unweighted  429.       185
    ## Weighted     99.82     185

These weights have quite high variability, and yield an ESS of close to
100 in the control group. Let’s see if these weights managed to yield
balance on our covariates.

``` r
bal.tab(W.out, stats = c("m", "v"),
        thresholds = c(m = .05))
```

    ## Balance Measures
    ##                 Type Diff.Adj         M.Threshold V.Ratio.Adj
    ## prop.score  Distance   -0.021     Balanced, <0.05       1.032
    ## age          Contin.    0.119 Not Balanced, >0.05       0.458
    ## educ         Contin.   -0.028     Balanced, <0.05       0.664
    ## race_black    Binary   -0.002     Balanced, <0.05           .
    ## race_hispan   Binary    0.000     Balanced, <0.05           .
    ## race_white    Binary    0.002     Balanced, <0.05           .
    ## married       Binary    0.019     Balanced, <0.05           .
    ## nodegree      Binary    0.018     Balanced, <0.05           .
    ## re74         Contin.   -0.002     Balanced, <0.05       1.321
    ## re75         Contin.    0.011     Balanced, <0.05       1.394
    ## 
    ## Balance tally for mean differences
    ##                     count
    ## Balanced, <0.05         9
    ## Not Balanced, >0.05     1
    ## 
    ## Variable with the greatest mean difference
    ##  Variable Diff.Adj         M.Threshold
    ##       age    0.119 Not Balanced, >0.05
    ## 
    ## Effective sample sizes
    ##            Control Treated
    ## Unadjusted  429.       185
    ## Adjusted     99.82     185

For nearly all the covariates, these weights yielded very good balance.
Only `age` remained imbalanced, with a standardized mean difference
greater than .05 and a variance ratio greater than 2. Let’s see if we
can do better. We’ll choose a different method: entropy balancing
([Hainmueller 2012](#ref-hainmuellerEntropyBalancingCausal2012)), which
guarantees perfect balance on specified moments of the covariates while
minimizing the negative entropy (a measure of dispersion) of the
weights.

``` r
W.out <- weightit(treat ~ age + educ + race + married + nodegree +
                    re74 + re75,
                  data = lalonde,
                  estimand = "ATT",
                  method = "ebal")
summary(W.out)
```

    ##                   Summary of weights
    ## 
    ## - Weight ranges:
    ## 
    ##           Min                                Max
    ## treated 1.       ||                         1.  
    ## control 0.019 |---------------------------| 9.42
    ## 
    ## - Units with the 5 most extreme weights by group:
    ##                                    
    ##              1     2   3     4    5
    ##  treated     1     1   1     1    1
    ##            423   412 226   196  118
    ##  control 7.127 7.501   8 9.036 9.42
    ## 
    ## - Weight statistics:
    ## 
    ##         Coef of Var   MAD Entropy # Zeros
    ## treated       0.000 0.000   0.000       0
    ## control       1.834 1.287   1.101       0
    ## 
    ## - Effective Sample Sizes:
    ## 
    ##            Control Treated
    ## Unweighted  429.       185
    ## Weighted     98.46     185

The variability of the weights has not changed much, but let’s see if
there are any gains in terms of balance:

``` r
bal.tab(W.out, stats = c("m", "v"),
        thresholds = c(m = .05))
```

    ## Balance Measures
    ##                Type Diff.Adj     M.Threshold V.Ratio.Adj
    ## age         Contin.        0 Balanced, <0.05       0.410
    ## educ        Contin.        0 Balanced, <0.05       0.664
    ## race_black   Binary        0 Balanced, <0.05           .
    ## race_hispan  Binary       -0 Balanced, <0.05           .
    ## race_white   Binary        0 Balanced, <0.05           .
    ## married      Binary        0 Balanced, <0.05           .
    ## nodegree     Binary       -0 Balanced, <0.05           .
    ## re74        Contin.        0 Balanced, <0.05       1.326
    ## re75        Contin.       -0 Balanced, <0.05       1.335
    ## 
    ## Balance tally for mean differences
    ##                     count
    ## Balanced, <0.05         9
    ## Not Balanced, >0.05     0
    ## 
    ## Variable with the greatest mean difference
    ##  Variable Diff.Adj     M.Threshold
    ##      re75       -0 Balanced, <0.05
    ## 
    ## Effective sample sizes
    ##            Control Treated
    ## Unadjusted  429.       185
    ## Adjusted     98.46     185

Indeed, we have achieved perfect balance on the means of the covariates.
However, the variance ratio of `age` is still quite high. We could
continue to try to adjust for this imbalance, but if there is reason to
believe it is unlikely to affect the outcome, it may be best to leave it
as is. (You can try adding `I(age^2)` to the formula and see what
changes this causes.)

Now that we have our weights stored in `W.out`, let’s estimate our
treatment effect in the weighted sample. The functions
[`lm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
and friends make it easy to fit (generalized) linear models that account
for estimation of of the weights in their standard errors. We can then
use functions in *marginaleffects* to perform g-computation to extract a
treatment effect estimation from the outcome model.

``` r
# Fit outcome model
fit <- lm_weightit(re78 ~ treat * (age + educ + race + married +
                                     nodegree + re74 + re75),
                   data = lalonde, weightit = W.out)
```

``` r
# G-computation for the treatment effect
library("marginaleffects")

avg_comparisons(fit, variables = "treat",
                newdata = subset(treat == 1))
```

    ## 
    ##  Estimate Std. Error    z Pr(>|z|)   S 2.5 % 97.5 %
    ##      1273        770 1.65   0.0983 3.3  -236   2783
    ## 
    ## Term: treat
    ## Type: probs
    ## Comparison: 1 - 0

Our confidence interval for `treat` contains 0, so there isn’t evidence
that `treat` has an effect on `re78`. Several types of standard errors
are available in *WeightIt*, including analytical standard errors that
account for estimation of the weights using M-estimation, robust
standard errors that treat the weights as fixed, and bootstrapping.
These are described in detail at
[`vignette("estimating-effects")`](https://ngreifer.github.io/WeightIt/articles/estimating-effects.md).

## Balancing Weights for a Longitudinal Treatment

*WeightIt* can estimate weights marginal structural models with
longitudinal treatment as well. This time, we’ll use the sample data set
`msmdata` to estimate our weights. Data must be in “wide” format, with
one row per unit.

``` r
data("msmdata")

head(msmdata)
```

    ##   X1_0 X2_0 A_1 X1_1 X2_1 A_2 X1_2 X2_2 A_3 Y_B
    ## 1    2    0   1    5    1   0    4    1   0   0
    ## 2    4    0   1    9    0   1   10    0   1   1
    ## 3    4    1   0    5    0   1    4    0   0   1
    ## 4    4    1   0    4    0   0    6    1   0   1
    ## 5    6    1   1    5    0   1    6    0   0   1
    ## 6    5    1   0    4    0   1    4    0   1   0

We have a binary outcome variable (`Y_B`), pre-treatment time-varying
variables (`X1_0` and `X2_0`, measured before the first treatment,
`X1_1` and `X2_1` measured between the first and second treatments, and
`X1_2` and `X2_2` measured between the second and third treatments), and
three time-varying binary treatment variables (`A_1`, `A_2`, and `A_3`).
We are interested in the joint, unique, causal effects of each treatment
period on the outcome. At each treatment time point, we need to achieve
balance on all variables measured prior to that treatment, including
previous treatments.

Using *cobalt*, we can examine the initial imbalance at each time point
and overall:

``` r
library("cobalt") #if not already attached

bal.tab(list(A_1 ~ X1_0 + X2_0,
             A_2 ~ X1_1 + X2_1 +
               A_1 + X1_0 + X2_0,
             A_3 ~ X1_2 + X2_2 +
               A_2 + X1_1 + X2_1 +
               A_1 + X1_0 + X2_0),
        data = msmdata, stats = c("m", "ks"),
        which.time = .all)
```

    ## Balance by Time Point
    ## 
    ##  - - - Time: 1 - - - 
    ## Balance Measures
    ##         Type Diff.Un KS.Un
    ## X1_0 Contin.   0.690 0.276
    ## X2_0  Binary  -0.325 0.325
    ## 
    ## Sample sizes
    ##     Control Treated
    ## All    3306    4194
    ## 
    ##  - - - Time: 2 - - - 
    ## Balance Measures
    ##         Type Diff.Un KS.Un
    ## X1_1 Contin.   0.874 0.340
    ## X2_1  Binary  -0.299 0.299
    ## A_1   Binary   0.127 0.127
    ## X1_0 Contin.   0.528 0.201
    ## X2_0  Binary  -0.060 0.060
    ## 
    ## Sample sizes
    ##     Control Treated
    ## All    3701    3799
    ## 
    ##  - - - Time: 3 - - - 
    ## Balance Measures
    ##         Type Diff.Un KS.Un
    ## X1_2 Contin.   0.475 0.212
    ## X2_2  Binary  -0.594 0.594
    ## A_2   Binary   0.162 0.162
    ## X1_1 Contin.   0.573 0.237
    ## X2_1  Binary  -0.040 0.040
    ## A_1   Binary   0.100 0.100
    ## X1_0 Contin.   0.361 0.148
    ## X2_0  Binary  -0.040 0.040
    ## 
    ## Sample sizes
    ##     Control Treated
    ## All    4886    2614
    ##  - - - - - - - - - - -

[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
indicates significant imbalance on most covariates at most time points,
so we need to do some work to eliminate that imbalance in our weighted
data set. We’ll use the
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
function to specify our weight models. The syntax is similar both to
that of
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
for point treatments and to that of
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html)
for longitudinal treatments. We’ll use `method = "glm"` and
`stabilize = TRUE` for stabilized propensity score weights estimated
using logistic regression.

``` r
Wmsm.out <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                             A_2 ~ X1_1 + X2_1 +
                               A_1 + X1_0 + X2_0,
                             A_3 ~ X1_2 + X2_2 +
                               A_2 + X1_1 + X2_1 +
                               A_1 + X1_0 + X2_0),
                        data = msmdata, method = "glm",
                        stabilize = TRUE)

Wmsm.out
```

    ## A weightitMSM object
    ##  - method: "glm" (propensity score weighting with GLM)
    ##  - number of obs.: 7500
    ##  - sampling weights: none
    ##  - number of time points: 3 (A_1, A_2, A_3)
    ##  - treatment:
    ##     + time 1: 2-category
    ##     + time 2: 2-category
    ##     + time 3: 2-category
    ##  - covariates:
    ##     + baseline: X1_0, X2_0
    ##     + after time 1: X1_1, X2_1, A_1, X1_0, X2_0
    ##     + after time 2: X1_2, X2_2, A_2, X1_1, X2_1, A_1, X1_0, X2_0
    ##  - stabilized; stabilization factors:
    ##     + baseline: (none)
    ##     + after time 1: A_1
    ##     + after time 2: A_1, A_2, A_1:A_2

[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
estimates separate weights for each time period and then takes the
product of the weights for each individual to arrive at the final
estimated weights. Printing the output of
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
provides some details about the function call and the output. We can
take a look at the quality of the weights with
[`summary()`](https://rdrr.io/r/base/summary.html), just as we could for
point treatments.

``` r
summary(Wmsm.out)
```

    ##                         Time 1                        
    ##                   Summary of weights
    ## 
    ## - Weight ranges:
    ## 
    ##           Min                                 Max
    ## treated 0.153 |---------------------------| 57.08
    ## control 0.109 |--------|                    20.46
    ## 
    ## - Units with the 5 most extreme weights by group:
    ##                                            
    ##            3172   2446   2115   2025   1938
    ##  treated 22.101 24.128   25.7 27.786 57.079
    ##            2943   2778   2726   1121    832
    ##  control 12.894  13.09 14.523 14.705 20.465
    ## 
    ## - Weight statistics:
    ## 
    ##         Coef of Var   MAD Entropy # Zeros
    ## treated       1.779 0.775   0.573       0
    ## control       1.331 0.752   0.486       0
    ## 
    ## - Mean of Weights = 0.981
    ## 
    ## - Effective Sample Sizes:
    ## 
    ##            Control Treated
    ## Unweighted    3306    4194
    ## Weighted      1193    1007
    ## 
    ##                         Time 2                        
    ##                   Summary of weights
    ## 
    ## - Weight ranges:
    ## 
    ##           Min                                 Max
    ## treated 0.109 |---------------------------| 57.08
    ## control 0.15  |--------|                    20.49
    ## 
    ## - Units with the 5 most extreme weights by group:
    ##                                            
    ##            2902   2266   1954   1869   1779
    ##  treated 22.101 24.128   25.7 27.786 57.079
    ##            3389   3023   3021    911    620
    ##  control 14.523 14.705 14.808 16.231 20.486
    ## 
    ## - Weight statistics:
    ## 
    ##         Coef of Var   MAD Entropy # Zeros
    ## treated       1.797 0.779   0.580       0
    ## control       1.359 0.750   0.488       0
    ## 
    ## - Mean of Weights = 0.991
    ## 
    ## - Effective Sample Sizes:
    ## 
    ##            Control Treated
    ## Unweighted    3701  3799. 
    ## Weighted      1300   898.2
    ## 
    ##                         Time 3                        
    ##                   Summary of weights
    ## 
    ## - Weight ranges:
    ## 
    ##           Min                                 Max
    ## treated 0.109 |---------------------------| 57.08
    ## control 0.208 |-----------|                 25.7 
    ## 
    ## - Units with the 5 most extreme weights by group:
    ##                                            
    ##            1991   1558   1254   1244   1195
    ##  treated 20.583 22.101 24.128 27.786 57.079
    ##            4479   4025   4021   2455    112
    ##  control 14.705 14.808  16.97 20.486   25.7
    ## 
    ## - Weight statistics:
    ## 
    ##         Coef of Var   MAD Entropy # Zeros
    ## treated       2.008 0.931   0.753       0
    ## control       1.269 0.672   0.407       0
    ## 
    ## - Mean of Weights = 1.040.97
    ## 
    ## - Effective Sample Sizes:
    ## 
    ##            Control Treated
    ## Unweighted    4886  2614. 
    ## Weighted      1871   519.8

Displayed are summaries of how the weights perform at each time point
with respect to variability. Next, we’ll examine how well they perform
with respect to covariate balance.

``` r
bal.tab(Wmsm.out, stats = c("m", "ks"),
        which.time = .none)
```

    ## Balance summary across all time points
    ##        Times    Type Max.Diff.Adj Max.KS.Adj
    ## X1_0 1, 2, 3 Contin.        0.033      0.018
    ## X2_0 1, 2, 3  Binary        0.018      0.018
    ## X1_1    2, 3 Contin.        0.087      0.039
    ## X2_1    2, 3  Binary        0.031      0.031
    ## A_1     2, 3  Binary        0.130      0.130
    ## X1_2       3 Contin.        0.104      0.054
    ## X2_2       3  Binary        0.007      0.007
    ## A_2        3  Binary        0.154      0.154
    ## 
    ## Effective sample sizes
    ##  - Time 1
    ##            Control Treated
    ## Unadjusted    3306    4194
    ## Adjusted      1193    1007
    ##  - Time 2
    ##            Control Treated
    ## Unadjusted    3701  3799. 
    ## Adjusted      1300   898.2
    ##  - Time 3
    ##            Control Treated
    ## Unadjusted    4886  2614. 
    ## Adjusted      1871   519.8

By setting `which.time = .none` in
[`bal.tab()`](https://ngreifer.github.io/cobalt/reference/bal.tab.html),
we can focus on the overall balance assessment, which displays the
greatest imbalance for each covariate across time points. We can see
that our estimated weights balance all covariates all time points with
respect to means and KS statistics. Now we can estimate our treatment
effects.

First, we fit a marginal structural model for the outcome using
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
with the `weightit` object supplied:

``` r
# Fit outcome model
fit <- glm_weightit(Y_B ~ A_1 * A_2 * A_3 * (X1_0 + X2_0),
                    data = msmdata,
                    weightit = Wmsm.out,
                    family = binomial)
```

Then, we compute the average expected potential outcomes under each
treatment regime using
[`marginaleffects::avg_predictions()`](https://marginaleffects.com/man/r/predictions.html):

``` r
library("marginaleffects")

p <- avg_predictions(fit, variables = c("A_1", "A_2", "A_3"))

p
```

    ## 
    ##  A_1 A_2 A_3 Estimate Std. Error    z Pr(>|z|)     S 2.5 % 97.5 %
    ##    0   0   0    0.687     0.0166 41.4   <0.001   Inf 0.654  0.719
    ##    0   0   1    0.521     0.0379 13.7   <0.001 140.3 0.447  0.595
    ##    0   1   0    0.491     0.0213 23.1   <0.001 389.1 0.449  0.532
    ##    0   1   1    0.438     0.0295 14.8   <0.001 163.2 0.380  0.496
    ##    1   0   0    0.602     0.0211 28.5   <0.001 590.8 0.561  0.644
    ##    1   0   1    0.544     0.0314 17.3   <0.001 221.0 0.482  0.605
    ##    1   1   0    0.378     0.0163 23.2   <0.001 393.1 0.346  0.410
    ##    1   1   1    0.422     0.0261 16.1   <0.001 192.3 0.371  0.473
    ## 
    ## Type: probs

We can compare the expected potential outcomes under each regime using
[`marginaleffects::hypotheses()`](https://marginaleffects.com/man/r/hypotheses.html).
To get all pairwise comparisons, supply the
[`avg_predictions()`](https://marginaleffects.com/man/r/predictions.html)
output to `hypotheses(., ~ pairwise)`. To compare individual regimes, we
can use
[`hypotheses()`](https://marginaleffects.com/man/r/hypotheses.html),
identifying the rows of the
[`avg_predictions()`](https://marginaleffects.com/man/r/predictions.html)
output. For example, to compare the regimes with no treatment for all
three time points vs. the regime with treatment for all three time
points, we would run

``` r
hypotheses(p, "b8 - b1 = 0")
```

    ## 
    ##  Hypothesis Estimate Std. Error     z Pr(>|z|)    S  2.5 % 97.5 %
    ##     b8-b1=0   -0.265     0.0308 -8.61   <0.001 56.9 -0.325 -0.204

These results indicate that receiving treatment at all time points
reduces the risk of the outcome relative to not receiving treatment at
all.

## References

Austin, Peter C. 2011. “An Introduction to Propensity Score Methods for
Reducing the Effects of Confounding in Observational Studies.”
*Multivariate Behavioral Research* 46 (3): 399–424.
<https://doi.org/10.1080/00273171.2011.568786>.

Austin, Peter C., and Elizabeth A. Stuart. 2015. “Moving Towards Best
Practice When Using Inverse Probability of Treatment Weighting (IPTW)
Using the Propensity Score to Estimate Causal Treatment Effects in
Observational Studies.” *Statistics in Medicine* 34 (28): 3661–79.
<https://doi.org/10.1002/sim.6607>.

Hainmueller, J. 2012. “Entropy Balancing for Causal Effects: A
Multivariate Reweighting Method to Produce Balanced Samples in
Observational Studies.” *Political Analysis* 20 (1): 25–46.
<https://doi.org/10.1093/pan/mpr025>.

Robins, James M., Miguel Ángel Hernán, and Babette Brumback. 2000.
“Marginal Structural Models and Causal Inference in Epidemiology.”
*Epidemiology* 11 (5): 550–60.
<https://doi.org/10.1097/00001648-200009000-00011>.

Thoemmes, Felix J., and Anthony D. Ong. 2016. “A Primer on Inverse
Probability of Treatment Weighting and Marginal Structural Models.”
*Emerging Adulthood* 4 (1): 40–59.
<https://doi.org/10.1177/2167696815621645>.
