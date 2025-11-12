# Generate Balancing Weights for Longitudinal Treatments

`weightitMSM()` allows for the easy generation of balancing weights for
marginal structural models for time-varying treatments using a variety
of available methods for binary, continuous, and multi-category
treatments. Many of these methods exist in other packages, which
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
calls; these packages must be installed to use the desired method.

## Usage

``` r
weightitMSM(
  formula.list,
  data = NULL,
  method = "glm",
  stabilize = FALSE,
  by = NULL,
  s.weights = NULL,
  num.formula = NULL,
  missing = NULL,
  verbose = FALSE,
  include.obj = FALSE,
  keep.mparts = TRUE,
  is.MSM.method,
  weightit.force = FALSE,
  ...
)
```

## Arguments

- formula.list:

  a list of formulas corresponding to each time point with the
  time-specific treatment variable on the left hand side and
  pre-treatment covariates to be balanced on the right hand side. The
  formulas must be in temporal order, and must contain all covariates to
  be balanced at that time point (i.e., treatments and covariates
  featured in early formulas should appear in later ones). Interactions
  and functions of covariates are allowed.

- data:

  an optional data set in the form of a data frame that contains the
  variables in the formulas in `formula.list`. This must be a wide data
  set with exactly one row per unit.

- method:

  a string of length 1 containing the name of the method that will be
  used to estimate weights. See
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  for allowable options. The default is `"glm"`, which estimates the
  weights using generalized linear models.

- stabilize:

  `logical`; whether or not to stabilize the weights. Stabilizing the
  weights involves fitting a model predicting treatment at each time
  point from treatment status at prior time points. If `TRUE`, a fully
  saturated model will be fit (i.e., all interactions between all
  treatments up to each time point), essentially using the observed
  treatment probabilities in the numerator (for binary and
  multi-category treatments). This may yield an error if some
  combinations are not observed. Default is `FALSE`. To manually specify
  stabilization model formulas, e.g., to specify non-saturated models,
  use `num.formula`. With many time points, saturated models may be
  time-consuming or impossible to fit.

- by:

  a string containing the name of the variable in `data` for which
  weighting is to be done within categories or a one-sided formula with
  the stratifying variable on the right-hand side. For example, if
  `by = "gender"` or `by = ~gender`, a separate propensity score model
  or optimization will occur within each level of the variable
  `"gender"`. Only one `by` variable is allowed; to stratify by multiply
  variables simultaneously, create a new variable that is a full cross
  of those variables using
  [`interaction()`](https://rdrr.io/r/base/interaction.html).

- s.weights:

  a vector of sampling weights or the name of a variable in `data` that
  contains sampling weights. These can also be matching weights if
  weighting is to be used on matched data. See the individual pages for
  each method for information on whether sampling weights can be
  supplied.

- num.formula:

  optional; a one-sided formula with the stabilization factors (other
  than the previous treatments) on the right hand side, which adds, for
  each time point, the stabilization factors to a model saturated with
  previous treatments. See Cole & Hernán (2008) for a discussion of how
  to specify this model; including stabilization factors can change the
  estimand without proper adjustment, and should be done with caution.
  Can also be a list of one-sided formulas, one for each time point.
  Unless you know what you are doing, we recommend setting
  `stabilize = TRUE` and ignoring `num.formula`.

- missing:

  `character`; how missing data should be handled. The options and
  defaults depend on the `method` used. Ignored if no missing data is
  present. It should be noted that multiple imputation outperforms all
  available missingness methods available in
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  and should probably be used instead. Consider the
  [MatchThem](https://CRAN.R-project.org/package=MatchThem) package for
  the use of
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  with multiply imputed data.

- verbose:

  `logical`; whether to print additional information output by the
  fitting function.

- include.obj:

  `logical`; whether to include in the output a list of the fit objects
  created in the process of estimating the weights at each time point.
  For example, with `method = "glm"`, a list of the `glm` objects
  containing the propensity score models at each time point will be
  included. See the help pages for each method for information on what
  object will be included if `TRUE`.

- keep.mparts:

  `logical`; whether to include in the output components necessary to
  estimate standard errors that account for estimation of the weights in
  [`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md).
  Default is `TRUE` if such parts are present. See the individual pages
  for each method for whether these components are produced. Set to
  `FALSE` to keep the output object smaller, e.g., if standard errors
  will not be computed using
  [`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md).

- is.MSM.method:

  `logical`; whether the method estimates weights for multiple time
  points all at once (`TRUE`) or by estimating weights at each time
  point and then multiplying them together (`FALSE`). This is only
  relevant for user-specified functions.

- weightit.force:

  `logical`; several methods are not valid for estimating weights with
  longitudinal treatments, and will produce an error message if
  attempted. Set to `TRUE` to bypass this error message.

- ...:

  other arguments for functions called by
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  that control aspects of fitting that are not covered by the above
  arguments. See Details at
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md).

## Value

A `weightitMSM` object with the following elements:

- weights:

  The estimated weights, one for each unit.

- treat.list:

  A list of the values of the time-varying treatment variables.

- covs.list:

  A list of the covariates used in the fitting at each time point. Only
  includes the raw covariates, which may have been altered in the
  fitting process.

- data:

  The data.frame originally entered to `weightitMSM()`.

- estimand:

  "ATE", currently the only estimand for MSMs with binary or
  multi-category treatments.

- method:

  The weight estimation method specified.

- ps.list:

  A list of the estimated propensity scores (if any) at each time point.

- s.weights:

  The provided sampling weights.

- by:

  A data.frame containing the `by` variable when specified.

- stabilization:

  The stabilization factors, if any.

When `keep.mparts` is `TRUE` (the default) and the chosen method is
compatible with M-estimation, the components related to M-estimation for
use in
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
are stored in the `"Mparts.list"` attribute. When `by` is specified,
`keep.mparts` is set to `FALSE`.

## Details

Currently only "wide" data sets, where each row corresponds to a unit's
entire variable history, are supported. You can use
[`reshape()`](https://rdrr.io/r/stats/reshape.html) or other functions
to transform your data into this format; see example below.

In general, `weightitMSM()` works by separating the estimation of
weights into separate procedures for each time period based on the
formulas provided. For each formula, `weightitMSM()` simply calls
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
to that formula, collects the weights for each time period, and
multiplies them together to arrive at longitudinal balancing weights.

Each formula should contain all the covariates to be balanced on. For
example, the formula corresponding to the second time period should
contain all the baseline covariates, the treatment variable at the first
time period, and the time-varying covariates that took on values after
the first treatment and before the second. Currently, only wide data
sets are supported, where each unit is represented by exactly one row
that contains the covariate and treatment history encoded in separate
variables.

The `"cbps"` method, which calls `CBPS()` in CBPS, will yield different
results from `CBMSM()` in CBPS because `CBMSM()` takes a different
approach to generating weights than simply estimating several
time-specific models.

## References

Cole, S. R., & Hernán, M. A. (2008). Constructing Inverse Probability
Weights for Marginal Structural Models. *American Journal of
Epidemiology*, 168(6), 656–664.
[doi:10.1093/aje/kwn164](https://doi.org/10.1093/aje/kwn164)

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
for information on the allowable methods

[`summary.weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/summary.weightit.md)
for summarizing the weights

## Examples

``` r
data("msmdata")
(W1 <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                        A_2 ~ X1_1 + X2_1 +
                          A_1 + X1_0 + X2_0,
                        A_3 ~ X1_2 + X2_2 +
                          A_2 + X1_1 + X2_1 +
                          A_1 + X1_0 + X2_0),
                   data = msmdata,
                   method = "glm"))
#> A weightitMSM object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 7500
#>  - sampling weights: none
#>  - number of time points: 3 (A_1, A_2, A_3)
#>  - treatment:
#>     + time 1: 2-category
#>     + time 2: 2-category
#>     + time 3: 2-category
#>  - covariates:
#>     + baseline: X1_0, X2_0
#>     + after time 1: X1_1, X2_1, A_1, X1_0, X2_0
#>     + after time 2: X1_2, X2_2, A_2, X1_1, X2_1, A_1, X1_0, X2_0
summary(W1)
#>                         Time 1                        
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                   Max
#> treated 1.079 |---------------------------| 403.483
#> control 1.276 |-------------------|         284.764
#> 
#> - Units with the 5 most extreme weights by group:
#>                                                 
#>             3172    3065    2025    1938     731
#>  treated 166.992 170.555 196.414 213.193 403.483
#>             2301    1275    1145    1121     832
#>  control 155.625 168.964  172.42 245.882 284.764
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       1.914 0.816   0.649       0
#> control       1.706 0.862   0.670       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted 3306.    4194. 
#> Weighted    845.79   899.4
#> 
#>                         Time 2                        
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                   Max
#> treated 1.079 |---------------------------| 403.483
#> control 1.276 |----------------|            245.882
#> 
#> - Units with the 5 most extreme weights by group:
#>                                                 
#>             2902    1869    1779    1509    1313
#>  treated 168.964 170.555 196.414 284.764 403.483
#>             2684    2549    1250     911     620
#>  control 155.625 166.992  172.42 213.193 245.882
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       1.892 0.819   0.652       0
#> control       1.748 0.869   0.686       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted 3701.   3799.  
#> Weighted    912.87  829.87
#> 
#>                         Time 3                        
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                   Max
#> treated 1.079 |---------------------------| 403.483
#> control 1.276 |---------|                   148.155
#> 
#> - Units with the 5 most extreme weights by group:
#>                                                 
#>             1991    1254     893     668     468
#>  treated 196.414 213.193 245.882 284.764 403.483
#>             4479    4021    2455    2427     112
#>  control  88.072  97.827 104.623 121.845 148.155
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       1.832 0.975   0.785       0
#> control       1.254 0.683   0.412       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted 4886.   2614.  
#> Weighted   1900.26  600.12
#> 
cobalt::bal.tab(W1)
#> Balance summary across all time points
#>        Times    Type Max.Diff.Adj
#> X1_0 1, 2, 3 Contin.       0.0342
#> X2_0 1, 2, 3  Binary       0.0299
#> X1_1    2, 3 Contin.       0.0657
#> X2_1    2, 3  Binary       0.0299
#> A_1     2, 3  Binary       0.0262
#> X1_2       3 Contin.       0.0643
#> X2_2       3  Binary       0.0096
#> A_2        3  Binary       0.0054
#> 
#> Effective sample sizes
#>  - Time 1
#>            Control Treated
#> Unadjusted 3306.    4194. 
#> Adjusted    845.79   899.4
#>  - Time 2
#>            Control Treated
#> Unadjusted 3701.   3799.  
#> Adjusted    912.87  829.87
#>  - Time 3
#>            Control Treated
#> Unadjusted 4886.   2614.  
#> Adjusted   1900.26  600.12

# Using stabilization factors
W2 <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                        A_2 ~ X1_1 + X2_1 +
                          A_1 + X1_0 + X2_0,
                        A_3 ~ X1_2 + X2_2 +
                          A_2 + X1_1 + X2_1 +
                          A_1 + X1_0 + X2_0),
                   data = msmdata,
                   method = "glm",
                   stabilize = TRUE,
                   num.formula = list(~ 1,
                                      ~ A_1,
                                      ~ A_1 + A_2))

# Same as above but with fully saturated stabilization factors
# (i.e., making the last entry in 'num.formula' A_1*A_2)
W3 <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                        A_2 ~ X1_1 + X2_1 +
                          A_1 + X1_0 + X2_0,
                        A_3 ~ X1_2 + X2_2 +
                          A_2 + X1_1 + X2_1 +
                          A_1 + X1_0 + X2_0),
                   data = msmdata,
                   method = "glm",
                   stabilize = TRUE)
```
