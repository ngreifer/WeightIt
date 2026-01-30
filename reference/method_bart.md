# Propensity Score Weighting Using BART

This page explains the details of estimating weights from Bayesian
additive regression trees (BART)-based propensity scores by setting
`method = "bart"` in the call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
This method can be used with binary, multi-category, and continuous
treatments.

In general, this method relies on estimating propensity scores using
BART and then converting those propensity scores into weights using a
formula that depends on the desired estimand. This method relies on
[`dbarts::bart2()`](https://rdrr.io/pkg/dbarts/man/bart.html) from the
[dbarts](https://CRAN.R-project.org/package=dbarts) package.

### Binary Treatments

For binary treatments, this method estimates the propensity scores using
[`dbarts::bart2()`](https://rdrr.io/pkg/dbarts/man/bart.html) . The
following estimands are allowed: ATE, ATT, ATC, ATO, ATM, and ATOS.
Weights can also be computed using marginal mean weighting through
stratification for the ATE, ATT, and ATC. See
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
for details.

### Multi-Category Treatments

For multi-category treatments, the propensity scores are estimated using
several calls to
[`dbarts::bart2()`](https://rdrr.io/pkg/dbarts/man/bart.html) , one for
each treatment group; the treatment probabilities are not normalized to
sum to 1. The following estimands are allowed: ATE, ATT, ATC, ATO, and
ATM. The weights for each estimand are computed using the standard
formulas or those mentioned above. Weights can also be computed using
marginal mean weighting through stratification for the ATE, ATT, and
ATC. See
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
for details.

### Continuous Treatments

For continuous treatments, weights are estimated as \\w_i = f_A(a_i) /
f\_{A\|X}(a_i)\\, where \\f_A(a_i)\\ (known as the stabilization factor)
is the unconditional density of treatment evaluated the observed
treatment value and \\f\_{A\|X}(a_i)\\ (known as the generalized
propensity score) is the conditional density of treatment given the
covariates evaluated at the observed value of treatment. The shape of
\\f_A(.)\\ and \\f\_{A\|X}(.)\\ is controlled by the `density` argument
described below (normal distributions by default), and the predicted
values used for the mean of the conditional density are estimated using
BART as implemented in
[`dbarts::bart2()`](https://rdrr.io/pkg/dbarts/man/bart.html) . Kernel
density estimation can be used instead of assuming a specific density
for the numerator and denominator by setting `density = "kernel"`. Other
arguments to [`density()`](https://rdrr.io/r/stats/density.html) can be
specified to refine the density estimation parameters.

### Longitudinal Treatments

For longitudinal treatments, the weights are the product of the weights
estimated at each time point.

### Sampling Weights

Sampling weights are not supported.

### Missing Data

In the presence of missing data, the following value(s) for `missing`
are allowed:

- `"ind"` (default):

  First, for each variable with missingness, a new missingness indicator
  variable is created which takes the value 1 if the original covariate
  is `NA` and 0 otherwise. The missingness indicators are added to the
  model formula as main effects. The missing values in the covariates
  are then replaced with the covariate medians. The weight estimation
  then proceeds with this new formula and set of covariates. The
  covariates output in the resulting `weightit` object will be the
  original covariates with the `NA`s.

### M-estimation

M-estimation is not supported.

## Details

BART works by fitting a sum-of-trees model for the treatment or
probability of treatment. The number of trees is determined by the
`n.trees` argument. Bayesian priors are used for the hyperparameters, so
the result is a posterior distribution of predicted values for each
unit. The mean of these for each unit is taken for use in computing the
(generalized) propensity score. Although the hyperparameters governing
the priors can be modified by supplying arguments to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
that are passed to the BART fitting function, the default values tend to
work well and require little modification (though the defaults differ
for continuous and categorical treatments; see the
[`dbarts::bart2()`](https://rdrr.io/pkg/dbarts/man/bart.html)
documentation for details). Unlike many other machine learning methods,
no loss function is optimized and the hyperparameters do not need to be
tuned (e.g., using cross-validation), though performance can benefit
from tuning. BART tends to balance sparseness with flexibility by using
very weak learners as the trees, which makes it suitable for capturing
complex functions without specifying a particular functional form and
without overfitting.

### Reproducibility

BART has a random component, so some work must be done to ensure
reproducibility across runs. See the *Reproducibility* section at
[`dbarts::bart2()`](https://rdrr.io/pkg/dbarts/man/bart.html) for more
details. To ensure reproducibility, one can do one of two things:

1.  supply an argument to `seed`, which is passed to
    [`dbarts::bart2()`](https://rdrr.io/pkg/dbarts/man/bart.html) and
    sets the seed for single- and multi-threaded uses, or

2.  call [`set.seed()`](https://rdrr.io/r/base/Random.html) and set
    `n.threads = 1` to use single-threading. Note that to ensure
    reproducibility on any machine, regardless of the number of cores
    available, one should use single-threading by setting
    `n.threads = 1` and either supply `seed` or call
    [`set.seed()`](https://rdrr.io/r/base/Random.html).

## Additional Arguments

All arguments to
[`dbarts::bart2()`](https://rdrr.io/pkg/dbarts/man/bart.html) can be
passed through
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md),
with the following exceptions:

- `test`, `weights`,`subset`, `offset.test` are ignored

- `combine.chains` is always set to `TRUE`

- `sampleronly` is always set to `FALSE`

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

## Additional Outputs

- `obj`:

  When `include.obj = TRUE`, the `bart2` fit(s) used to generate the
  predicted values. With multi-category treatments, this will be a list
  of the fits; otherwise, it will be a single fit. The predicted
  probabilities used to compute the propensity scores can be extracted
  using [`fitted()`](https://rdrr.io/r/stats/fitted.values.html).

## References

Hill, J., Weiss, C., & Zhai, F. (2011). Challenges With Propensity Score
Strategies in a High-Dimensional Setting and a Potential Alternative.
*Multivariate Behavioral Research*, 46(3), 477–513.
[doi:10.1080/00273171.2011.570161](https://doi.org/10.1080/00273171.2011.570161)

Chipman, H. A., George, E. I., & McCulloch, R. E. (2010). BART: Bayesian
additive regression trees. *The Annals of Applied Statistics*, 4(1),
266–298. [doi:10.1214/09-AOAS285](https://doi.org/10.1214/09-AOAS285)

Note that many references that deal with BART for causal inference focus
on estimating potential outcomes with BART, not the propensity scores,
and so are not directly relevant when using BART to estimate propensity
scores for weights.

See
[`method_glm`](https://ngreifer.github.io/WeightIt/reference/method_glm.md)
for additional references on propensity score weighting more generally.

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md),
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)

[`method_super`](https://ngreifer.github.io/WeightIt/reference/method_super.md)
for stacking predictions from several machine learning methods,
including BART.

## Examples

``` r
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "bart", estimand = "ATT"))
#> A weightit object
#>  - method: "bart" (propensity score weighting with BART)
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
#>           Min                                 Max
#> treated 1.       ||                         1.   
#> control 0.003 |---------------------------| 9.501
#> 
#> - Units with the 5 most extreme weights by group:
#>                                      
#>              1     2     3    4     5
#>  treated     1     1     1    1     1
#>            423   407   384  224   189
#>  control 2.334 2.819 2.907 3.32 9.501
#> 
#> - Weight statistics:
#> 
#>         Coef of Var  MAD Entropy # Zeros
#> treated       0.000 0.00   0.000       0
#> control       1.794 0.92   0.716       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    101.86     185

cobalt::bal.tab(W1)
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.4942
#> age         Contin.   0.0583
#> educ        Contin.  -0.0228
#> married      Binary  -0.0320
#> nodegree     Binary   0.0358
#> re74        Contin.  -0.0553
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    101.86     185

#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "bart", estimand = "ATE"))
#> A weightit object
#>  - method: "bart" (propensity score weighting with BART)
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
#> black  1.228 |----------------|             9.034
#> hispan 2.592    |------------------------| 13.343
#> white  1.068 |-------------|                7.829
#> 
#> - Units with the 5 most extreme weights by group:
#>                                           
#>            192    171    166    164    152
#>   black  6.926  7.379  7.863  8.257  9.034
#>             69     67     59     50     37
#>  hispan 11.872 12.715 12.871 12.892 13.343
#>             15      7      6      5      3
#>   white  4.429  4.939  5.435  7.594  7.829
#> 
#> - Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.579 0.376   0.127       0
#> hispan       0.388 0.317   0.075       0
#> white        0.447 0.316   0.080       0
#> 
#> - Effective Sample Sizes:
#> 
#>             black hispan  white
#> Unweighted 243.    72.   299.  
#> Weighted   182.12  62.71 249.39

cobalt::bal.tab(W2)
#> Balance summary across all treatment pairs
#>             Type Max.Diff.Adj
#> age      Contin.       0.1919
#> educ     Contin.       0.1580
#> married   Binary       0.0530
#> nodegree  Binary       0.0255
#> re74     Contin.       0.1139
#> 
#> Effective sample sizes
#>             black hispan  white
#> Unadjusted 243.    72.   299.  
#> Adjusted   182.12  62.71 249.39

#Balancing covariates with respect to re75 (continuous)
#assuming t(3) conditional density for treatment
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "bart", density = "dt_3"))
#> A weightit object
#>  - method: "bart" (propensity score weighting with BART)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: continuous
#>  - covariates: age, educ, married, nodegree, re74

summary(W3)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>      Min                                  Max
#> all 0.08 |---------------------------| 23.063
#> 
#> - Units with the 5 most extreme weights:
#>                                     
#>        487   486   484    469    310
#>  all 7.974 7.992 9.919 18.119 23.063
#> 
#> - Weight statistics:
#> 
#>     Coef of Var   MAD Entropy # Zeros
#> all       1.181 0.475   0.283       0
#> 
#> - Effective Sample Sizes:
#> 
#>             Total
#> Unweighted 614.  
#> Weighted   256.64

cobalt::bal.tab(W3)
#> Balance Measures
#>             Type Corr.Adj
#> age      Contin.   0.0292
#> educ     Contin.   0.0495
#> married   Binary   0.0771
#> nodegree  Binary  -0.0781
#> re74     Contin.   0.1150
#> 
#> Effective sample sizes
#>             Total
#> Unadjusted 614.  
#> Adjusted   256.64
```
