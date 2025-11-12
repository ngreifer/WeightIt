# Trim (Winsorize) Large Weights

Trims (i.e., winsorizes) large weights by setting all weights higher
than that at a given quantile to the weight at the quantile or to 0.
This can be useful in controlling extreme weights, which can reduce
effective sample size by enlarging the variability of the weights. Note
that by default, no observations are fully discarded when using
`trim()`, which may differ from the some uses of the word "trim" (see
the `drop` argument below).

## Usage

``` r
trim(x, ...)

# S3 method for class 'weightit'
trim(x, at = 0, lower = FALSE, drop = FALSE, ...)

# Default S3 method
trim(x, at = 0, lower = FALSE, treat = NULL, drop = FALSE, ...)
```

## Arguments

- x:

  a `weightit` object or a vector of weights.

- ...:

  not used.

- at:

  `numeric`; either the quantile of the weights above which weights are
  to be trimmed. A single number between .5 and 1, or the number of
  weights to be trimmed (e.g., `at = 3` for the top 3 weights to be set
  to the 4th largest weight).

- lower:

  `logical`; whether also to trim at the lower quantile (e.g., for
  `at = .9`, trimming at both .1 and .9, or for `at = 3`, trimming the
  top and bottom 3 weights). Default is `FALSE` to only trim the higher
  weights.

- drop:

  `logical`; whether to set the weights of the trimmed units to 0 or
  not. Default is `FALSE` to retain all trimmed units. Setting to `TRUE`
  may change the original targeted estimand when not the ATT or ATC.

- treat:

  a vector of treatment status for each unit. This should always be
  included when `x` is numeric, but you can get away with leaving it out
  if the treatment is continuous or the estimand is the ATE for binary
  or multi-category treatments.

## Value

If the input is a `weightit` object, the output will be a `weightit`
object with the weights replaced by the trimmed weights (or 0) and will
have an additional attribute, `"trim"`, equal to the quantile of
trimming.

If the input is a numeric vector of weights, the output will be a
numeric vector of the trimmed weights, again with the aforementioned
attribute.

## Details

`trim()` takes in a `weightit` object (the output of a call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md))
or a numeric vector of weights and trims (winsorizes) them to the
specified quantile. All weights above that quantile are set to the
weight at that quantile unless `drop = TRUE`, in which case they are set
to 0. If `lower = TRUE`, all weights below 1 minus the quantile are
trimmed. In general, trimming weights can increase imbalance but also
decreases the variability of the weights, improving precision at the
potential expense of unbiasedness (Cole & Hernán, 2008). See Lee,
Lessler, and Stuart (2011) and Thoemmes and Ong (2015) for discussions
and simulation results of trimming weights at various quantiles. Note
that trimming weights can also change the target population and
therefore the estimand.

When using `trim()` on a numeric vector of weights, it is helpful to
include the treatment vector as well. The helps determine the type of
treatment and estimand, which are used to specify how trimming is
performed. In particular, if the estimand is determined to be the ATT or
ATC, the weights of the target (i.e., focal) group are ignored, since
they should all be equal to 1. Otherwise, if the estimand is the ATE or
the treatment is continuous, all weights are considered for trimming. In
general, weights for any group for which all the weights are the same
will not be considered in the trimming.

## References

Cole, S. R., & Hernán, M. Á. (2008). Constructing Inverse Probability
Weights for Marginal Structural Models. *American Journal of
Epidemiology*, 168(6), 656–664.
[doi:10.1093/aje/kwn164](https://doi.org/10.1093/aje/kwn164)

Lee, B. K., Lessler, J., & Stuart, E. A. (2011). Weight Trimming and
Propensity Score Weighting. *PLoS ONE*, 6(3), e18174.
[doi:10.1371/journal.pone.0018174](https://doi.org/10.1371/journal.pone.0018174)

Thoemmes, F., & Ong, A. D. (2016). A Primer on Inverse Probability of
Treatment Weighting and Marginal Structural Models. *Emerging
Adulthood*, 4(1), 40–59.
[doi:10.1177/2167696815621645](https://doi.org/10.1177/2167696815621645)

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)

## Examples

``` r
library("cobalt")
data("lalonde", package = "cobalt")

(W <- weightit(treat ~ age + educ + married +
                 nodegree + re74, data = lalonde,
               method = "glm", estimand = "ATT"))
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74

summary(W)
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
#>  control 1.33 1.437 1.5 1.637 2.044
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

#Trimming the top and bottom 5 weights
trim(W, at = 5, lower = TRUE)
#> Trimming the top and bottom 5 weights where treat is not 1.
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74
#>  - weights trimmed at the top and bottom 5

#Trimming at 90th percentile
(W.trim <- trim(W, at = .9))
#> Trimming weights where treat is not 1 to 90%.
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74
#>  - weights trimmed at 90%

summary(W.trim)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                   Max
#> treated 1.                                 || 1.   
#> control 0.022   |-------------------------|   0.941
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>              1     2     3     4     5
#>  treated     1     1     1     1     1
#>             79    84   100   111   118
#>  control 0.941 0.941 0.941 0.941 0.941
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.000 0.000   0.000       0
#> control       0.766 0.682   0.303       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    270.58     185
#Note that only the control weights were trimmed

#Trimming a numeric vector of weights
all.equal(trim(W$weights, at = .9, treat = lalonde$treat),
          W.trim$weights)
#> Trimming weights where treat is not 1 to 90%.
#> [1] TRUE

#Dropping trimmed units
(W.trim <- trim(W, at = .9, drop = TRUE))
#> Setting weights beyond 90% where treat is not 1 to 0.
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74
#>  - weights trimmed at 90% and units dropped

summary(W.trim)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>         Min                                   Max
#> treated   1                              || 1.   
#> control   0   |-------------------------|   0.941
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>              1     2     3     4     5
#>  treated     1     1     1     1     1
#>            171   184   188   281   282
#>  control 0.941 0.941 0.941 0.941 0.941
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.000 0.000   0.000       0
#> control       0.881 0.757   0.303      40
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    241.72     185
#Note that we now have zeros in the control group

#Using made up data and as.weightit()
treat <- rbinom(500, 1, .3)
weights <- rchisq(500, df = 2)
W <- as.weightit(weights, treat = treat,
                 estimand = "ATE")

summary(W)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                  Max
#> treated 0.005 |----------------------|      11.484
#> control 0.004 |---------------------------| 13.9  
#> 
#> - Units with the 5 most extreme weights by group:
#>                                           
#>              75     46    39     33      9
#>  treated  6.131  6.444 6.591  9.601 11.484
#>             344    329   309    255      7
#>  control 10.478 11.104 11.27 13.524   13.9
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.973 0.684   0.398       0
#> control       1.039 0.772   0.454       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  361.    139.  
#> Weighted    173.91   71.65

summary(trim(W, at = .95))
#> Trimming weights to 95%.
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                 Max
#> treated 0.005 |---------------------------| 6.404
#> control 0.004 |---------------------------| 6.404
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>             75     9    33    39    46
#>  treated 6.131 6.404 6.404 6.404 6.404
#>              2     3     7    16    37
#>  control 6.404 6.404 6.404 6.404 6.404
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.875 0.668   0.361       0
#> control       0.918 0.745   0.402       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  361.    139.  
#> Weighted    196.13   78.98
```
