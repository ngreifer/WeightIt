# Propensity Score Weighting Using SuperLearner

This page explains the details of estimating weights from
SuperLearner-based propensity scores by setting `method = "super"` in
the call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
This method can be used with binary, multi-category, and continuous
treatments.

In general, this method relies on estimating propensity scores using the
SuperLearner algorithm for stacking predictions and then converting
those propensity scores into weights using a formula that depends on the
desired estimand. For binary and multi-category treatments, one or more
binary classification algorithms are used to estimate the propensity
scores as the predicted probability of being in each treatment given the
covariates. For continuous treatments, regression algorithms are used to
estimate generalized propensity scores as the conditional density of
treatment given the covariates. This method relies on
[`SuperLearner::SuperLearner()`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
from the [SuperLearner](https://CRAN.R-project.org/package=SuperLearner)
package.

### Binary Treatments

For binary treatments, this method estimates the propensity scores using
[`SuperLearner::SuperLearner()`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
. The following estimands are allowed: ATE, ATT, ATC, ATO, ATM, and
ATOS. Weights can also be computed using marginal mean weighting through
stratification for the ATE, ATT, and ATC. See
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
for details.

### Multi-Category Treatments

For multi-category treatments, the propensity scores are estimated using
several calls to
[`SuperLearner::SuperLearner()`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
, one for each treatment group; the treatment probabilities are not
normalized to sum to 1. The following estimands are allowed: ATE, ATT,
ATC, ATO, and ATM. The weights for each estimand are computed using the
standard formulas or those mentioned above. Weights can also be computed
using marginal mean weighting through stratification for the ATE, ATT,
and ATC. See
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
for details.

### Continuous Treatments

For continuous treatments, the generalized propensity score is estimated
using
[`SuperLearner::SuperLearner()`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
. In addition, kernel density estimation can be used instead of assuming
a normal density for the numerator and denominator of the generalized
propensity score by setting `density = "kernel"`. Other arguments to
[`density()`](https://rdrr.io/r/stats/density.html) can be specified to
refine the density estimation parameters. `plot = TRUE` can be specified
to plot the density for the numerator and denominator, which can be
helpful in diagnosing extreme weights.

### Longitudinal Treatments

For longitudinal treatments, the weights are the product of the weights
estimated at each time point.

### Sampling Weights

Sampling weights are supported through `s.weights` in all scenarios.

### Missing Data

In the presence of missing data, the following value(s) for `missing`
are allowed:

- `"ind"` (default):

  First, for each variable with missingness, a new missingness indicator
  variable is created which takes the value 1 if the original covariate
  is `NA` and 0 otherwise. The missingness indicators are added to the
  model formula as main effects. The missing values in the covariates
  are then replaced with the covariate medians (this value is arbitrary
  and does not affect estimation). The weight estimation then proceeds
  with this new formula and set of covariates. The covariates output in
  the resulting `weightit` object will be the original covariates with
  the `NA`s.

### M-estimation

M-estimation is not supported.

## Details

SuperLearner works by fitting several machine learning models to the
treatment and covariates and then taking a weighted combination of the
generated predicted values to use as the propensity scores, which are
then used to construct weights. The machine learning models used are
supplied using the `SL.library` argument; the more models are supplied,
the higher the chance of correctly modeling the propensity score. It is
a good idea to include parametric models, flexible and tree-based
models, and regularized models among the models selected. The predicted
values are combined using the method supplied in the `SL.method`
argument (which is nonnegative least squares by default). A benefit of
SuperLearner is that, asymptotically, it is guaranteed to perform as
well as or better than the best-performing method included in the
library. Using Balance SuperLearner by setting
`SL.method = "method.balance"` works by selecting the combination of
predicted values that minimizes an imbalance measure.

## Note

Some methods formerly available in SuperLearner are now in
SuperLearnerExtra, which can be found on GitHub at
<https://github.com/ecpolley/SuperLearnerExtra>.

The `criterion` argument used to be called `stop.method`, which is its
name in twang. `stop.method` still works for backward compatibility.
Additionally, the criteria formerly named as `es.mean`, `es.max`, and
`es.rms` have been renamed to `smd.mean`, `smd.max`, and `smd.rms`. The
former are used in twang and will still work with
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
for backward compatibility.

As of version 1.2.0, the default behavior for binary and multi-category
treatments is to stratify on the treatment when performing
cross-validation to ensure all treatment groups are represented in
cross-validation. To recover previous behavior, set
`cvControl = list(stratifyCV = FALSE)`.

## Additional Arguments

- `discrete`:

  if `TRUE`, uses discrete SuperLearner, which simply selects the best
  performing method. Default `FALSE`, which finds the optimal
  combination of predictions for the libraries using `SL.method`.

An argument to `SL.library` **must** be supplied. To see a list of
available entries, use
[`SuperLearner::listWrappers()`](https://rdrr.io/pkg/SuperLearner/man/listWrappers.html)
.

All arguments to
[`SuperLearner::SuperLearner()`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
can be passed through
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md),
with the following exceptions:

- `obsWeights` is ignored because sampling weights are passed using
  `s.weights`.

- `method` in `SuperLearner()` is replaced with the argument `SL.method`
  in
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md).

For binary and multi-category treatments, the following arguments may be
supplied:

- `subclass`:

  `integer`; the number of subclasses to use for computing weights using
  marginal mean weighting through stratification (MMWS). If `NULL`,
  standard inverse probability weights (and their extensions) will be
  computed; if a number greater than 1, subclasses will be formed and
  weights will be computed based on subclass membership. See
  [`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
  for details and references.

For continuous treatments, the following arguments may be supplied:

- `density`:

  A function corresponding to the conditional density of the treatment.
  The standardized residuals of the treatment model will be fed through
  this function to produce the numerator and denominator of the
  generalized propensity score weights. If blank,
  [`dnorm()`](https://rdrr.io/r/stats/Normal.html) is used as
  recommended by Robins et al. (2000). This can also be supplied as a
  string containing the name of the function to be called. If the string
  contains underscores, the call will be split by the underscores and
  the latter splits will be supplied as arguments to the second argument
  and beyond. For example, if `density = "dt_2"` is specified, the
  density used will be that of a t-distribution with 2 degrees of
  freedom. Using a t-distribution can be useful when extreme outcome
  values are observed (Naimi et al., 2014).

  Can also be `"kernel"` to use kernel density estimation, which calls
  [`density()`](https://rdrr.io/r/stats/density.html) to estimate the
  numerator and denominator densities for the weights. (This used to be
  requested by setting `use.kernel = TRUE`, which is now deprecated.)

- `bw`, `adjust`, `kernel`, `n`:

  If `density = "kernel"`, the arguments to
  [`density()`](https://rdrr.io/r/stats/density.html). The defaults are
  the same as those in
  [`density()`](https://rdrr.io/r/stats/density.html) except that `n` is
  10 times the number of units in the sample.

- `plot`:

  If `density = "kernel"`, whether to plot the estimated densities.

### Balance SuperLearner

In addition to the methods allowed by `SuperLearner()`, one can specify
`SL.method = "method.balance"` to use "Balance SuperLearner" as
described by Pirracchio and Carone (2018), wherein covariate balance is
used to choose the optimal combination of the predictions from the
methods specified with `SL.library`. Coefficients are chosen (one for
each prediction method) so that the weights generated from the weighted
combination of the predictions optimize a balance criterion, which must
be set with the `criterion` argument, described below.

- `criterion`:

  a string describing the balance criterion used to select the best
  weights. See
  [`cobalt::bal.compute()`](https://ngreifer.github.io/cobalt/reference/bal.compute.html)
  for allowable options for each treatment type. For binary and
  multi-category treatments, the default is `"smd.mean"`, which
  minimizes the average absolute standard mean difference among the
  covariates between treatment groups. For continuous treatments, the
  default is `"p.mean"`, which minimizes the average absolute Pearson
  correlation between the treatment and covariates.

Note that this implementation differs from that of Pirracchio and Carone
(2018) in that here, balance is measured only on the terms included in
the model formula (i.e., and not their interactions unless specifically
included), and balance results from a sample weighted using the
estimated predicted values as propensity scores, not a sample matched
using propensity score matching on the predicted values. Binary and
continuous treatments are supported, but currently multi-category
treatments are not.

## Additional Outputs

- `info`:

  For binary and continuous treatments, a list with two entries, `coef`
  and `cvRisk`. For multi-category treatments, a list of lists with
  these two entries, one for each treatment level.

  `coef`

  :   The coefficients in the linear combination of the predictions from
      each method in `SL.library`. Higher values indicate that the
      corresponding method plays a larger role in determining the
      resulting predicted value, and values close to zero indicate that
      the method plays little role in determining the predicted value.
      When `discrete = TRUE`, these correspond to the coefficients that
      would have been estimated had `discrete` been `FALSE`.

  `cvRisk`

  :   The cross-validation risk for each method in `SL.library`. Higher
      values indicate that the method has worse cross-validation
      accuracy. When `SL.method = "method.balance"`, the sample weighted
      balance statistic requested with `criterion`. Higher values
      indicate worse balance.

- `obj`:

  When `include.obj = TRUE`, the SuperLearner fit(s) used to generate
  the predicted values. For binary and continuous treatments, the output
  of the call to
  [`SuperLearner::SuperLearner()`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html)
  . For multi-category treatments, a list of outputs to calls to
  [`SuperLearner::SuperLearner()`](https://rdrr.io/pkg/SuperLearner/man/SuperLearner.html).

## References

### Binary treatments

Pirracchio, R., Petersen, M. L., & van der Laan, M. (2015). Improving
Propensity Score Estimators’ Robustness to Model Misspecification Using
Super Learner. *American Journal of Epidemiology*, 181(2), 108–119.
[doi:10.1093/aje/kwu253](https://doi.org/10.1093/aje/kwu253)

### Continuous treatments

Kreif, N., Grieve, R., Díaz, I., & Harrison, D. (2015). Evaluation of
the Effect of a Continuous Treatment: A Machine Learning Approach with
an Application to Treatment for Traumatic Brain Injury. *Health
Economics*, 24(9), 1213–1228.
[doi:10.1002/hec.3189](https://doi.org/10.1002/hec.3189)

### Balance SuperLearner (`SL.method = "method.balance"`)

Pirracchio, R., & Carone, M. (2018). The Balance Super Learner: A robust
adaptation of the Super Learner to improve estimation of the average
treatment effect in the treated based on propensity score matching.
*Statistical Methods in Medical Research*, 27(8), 2504–2518.
[doi:10.1177/0962280216682055](https://doi.org/10.1177/0962280216682055)

See
[`method_glm`](https://ngreifer.github.io/WeightIt/reference/method_glm.md)
for additional references.

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md),
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)

## Examples

``` r
data("lalonde", package = "cobalt")

#Note: for time, all exmaples use a small set of
#      learners. Many more should be added if
#      possible, including a variety of model
#      types (e.g., parametric, flexible, tree-
#      based, regularized, etc.)

# Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", estimand = "ATT",
                SL.library = c("SL.glm", "SL.stepAIC",
                               "SL.glm.interaction")))
#> Loading required package: nnls
#> A weightit object
#>  - method: "super" (propensity score weighting with SuperLearner)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74

summary(W1)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>          Min                                 Max
#> treated 1.          ||                     1.   
#> control 0.01 |---------------------------| 3.967
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>              1     2     3     4     5
#>  treated     1     1     1     1     1
#>            404   226   224   111    84
#>  control 2.052 2.135 2.277 2.691 3.967
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.000 0.000   0.000       0
#> control       0.962 0.718   0.386       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    223.01     185

cobalt::bal.tab(W1)
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.1346
#> age         Contin.  -0.0574
#> educ        Contin.   0.0241
#> married      Binary  -0.0050
#> nodegree     Binary   0.0219
#> re74        Contin.  -0.0295
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    223.01     185

# Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", estimand = "ATE",
                SL.library = c("SL.glm", "SL.stepAIC",
                               "SL.glm.interaction")))
#> A weightit object
#>  - method: "super" (propensity score weighting with SuperLearner)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 3-category (black, hispan, white)
#>  - estimand: ATE
#>  - covariates: age, educ, married, nodegree, re74

summary(W2)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>          Min                                  Max
#> black  1.399 |-----------------|           12.866
#> hispan 1.98   |--------------------------| 18.542
#> white  1.066 |----|                         5.   
#> 
#> - Units with the 5 most extreme weights by group:
#>                                           
#>            203    166    155    153    152
#>   black  7.234  7.849  9.856  11.77 12.866
#>             59     43     37     36     18
#>  hispan 14.316 14.355 14.397 14.509 18.542
#>            285    205    177     94      3
#>   white  4.303  4.308  4.381  4.942      5
#> 
#> - Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.620 0.385   0.131       0
#> hispan       0.417 0.340   0.089       0
#> white        0.394 0.324   0.072       0
#> 
#> - Effective Sample Sizes:
#> 
#>             black hispan  white
#> Unweighted 243.    72.   299.  
#> Weighted   175.81  61.46 258.91

cobalt::bal.tab(W2)
#> Balance summary across all treatment pairs
#>             Type Max.Diff.Adj
#> age      Contin.       0.1814
#> educ     Contin.       0.0754
#> married   Binary       0.0309
#> nodegree  Binary       0.0132
#> re74     Contin.       0.0332
#> 
#> Effective sample sizes
#>             black hispan  white
#> Unadjusted 243.    72.   299.  
#> Adjusted   175.81  61.46 258.91

# Balancing covariates with respect to re75 (continuous)
# assuming t(8) conditional density for treatment
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", density = "dt_8",
                SL.library = c("SL.glm", "SL.ridge",
                               "SL.glm.interaction")))
#> A weightit object
#>  - method: "super" (propensity score weighting with SuperLearner)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: continuous
#>  - covariates: age, educ, married, nodegree, re74

summary(W3)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>       Min                                  Max
#> all 0.044 |---------------------------| 22.616
#> 
#> - Units with the 5 most extreme weights:
#>                                      
#>        485   484    483    431    354
#>  all 9.587 16.21 18.466 20.328 22.616
#> 
#> - Weight statistics:
#> 
#>     Coef of Var  MAD Entropy # Zeros
#> all       1.375 0.51   0.352       0
#> 
#> - Effective Sample Sizes:
#> 
#>             Total
#> Unweighted 614.  
#> Weighted   212.62

cobalt::bal.tab(W3)
#> Balance Measures
#>             Type Corr.Adj
#> age      Contin.   0.0342
#> educ     Contin.   0.0372
#> married   Binary   0.0605
#> nodegree  Binary  -0.0627
#> re74     Contin.   0.0405
#> 
#> Effective sample sizes
#>             Total
#> Unadjusted 614.  
#> Adjusted   212.62

# Balancing covariates between treatment groups (binary)
# using balance SuperLearner to minimize the maximum
# KS statistic
(W4 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", estimand = "ATT",
                SL.library = c("SL.glm", "SL.stepAIC",
                               "SL.lda"),
                SL.method = "method.balance",
                criterion = "ks.max"))
#> A weightit object
#>  - method: "super" (propensity score weighting with SuperLearner)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74

summary(W4)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                 Max
#> treated 1.                  ||              1.   
#> control 0.022 |---------------------------| 2.044
#> 
#> - Units with the 5 most extreme weights by group:
#>                                    
#>             1     2   3     4     5
#>  treated    1     1   1     1     1
#>           410   226 224   111    84
#>  control 1.33 1.436 1.5 1.637 2.044
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.000 0.000    0.00       0
#> control       0.823 0.701    0.33       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    255.99     185

cobalt::bal.tab(W4, stats = c("m", "ks"))
#> Balance Measures
#>                Type Diff.Adj KS.Adj
#> prop.score Distance   0.0199 0.0944
#> age         Contin.   0.0459 0.2764
#> educ        Contin.  -0.0360 0.0601
#> married      Binary   0.0044 0.0044
#> nodegree     Binary   0.0080 0.0080
#> re74        Contin.  -0.0275 0.2839
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    255.99     185
```
