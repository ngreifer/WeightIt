# Calibrate Propensity Score Weights

`calibrate()` calibrates propensity scores used in weights. This
involves fitting a new propensity score model using logistic or isotonic
regression with the previously estimated propensity score as the sole
predictor. Weights are computed using this new propensity score.

## Usage

``` r
calibrate(x, ...)

# Default S3 method
calibrate(x, treat, s.weights = NULL, data = NULL, method = "platt", ...)

# S3 method for class 'weightit'
calibrate(x, method = "platt", ...)
```

## Arguments

- x:

  a `weightit` object or a vector of propensity scores. Only binary
  treatments are supported.

- ...:

  not used.

- treat:

  a vector of treatment status for each unit. Only binary treatments are
  supported.

- s.weights:

  a vector of sampling weights or the name of a variable in `data` that
  contains sampling weights.

- data:

  an optional data frame containing the variable named in `s.weights`
  when supplied as a string.

- method:

  `character`; the method of calibration used. Allowable options include
  `"platt"` (default) for Platt scaling as described by Gutman et
  al. (2024) and `"isoreg"` for isotonic regression as described by van
  der Laan et al. (2024).

## Value

If the input is a `weightit` object, the output will be a `weightit`
object with the propensity scores replaced with the calibrated
propensity scores and the weights replaced by weights computed from the
calibrated propensity scores.

If the input is a numeric vector of weights, the output will be a
numeric vector of the calibrated propensity scores.

## References

Gutman, R., Karavani, E., & Shimoni, Y. (2024). Improving Inverse
Probability Weighting by Post-calibrating Its Propensity Scores.
*Epidemiology*, 35(4).
[doi:10.1097/EDE.0000000000001733](https://doi.org/10.1097/EDE.0000000000001733)

van der Laan, L., Lin, Z., Carone, M., & Luedtke, A. (2024). Stabilized
Inverse Probability Weighting via Isotonic Calibration. arXiv.
<https://arxiv.org/abs/2411.06342>

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)

## Examples

``` r
library("cobalt")
data("lalonde", package = "cobalt")

#Using GBM to estimate weights
(W <- weightit(treat ~ age + educ + married +
                 nodegree + re74, data = lalonde,
               method = "gbm", estimand = "ATT",
               criterion = "smd.max"))
#> A weightit object
#>  - method: "gbm" (propensity score weighting with GBM)
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
#> treated 1.       ||                          1.  
#> control 0.004 |---------------------------| 10.61
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>              1     2     3     4     5
#>  treated     1     1     1     1     1
#>            423   407   384   400   189
#>  control 4.352 5.009 5.858 5.858 10.61
#> 
#> - Weight statistics:
#> 
#>         Coef of Var  MAD Entropy # Zeros
#> treated       0.000 0.00   0.000       0
#> control       2.422 1.06   1.016       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted   429.      185
#> Weighted      62.6     185

#Calibrating the GBM propensity scores
Wc <- calibrate(W)

#Calibrating propensity scores directly
PSc <- calibrate(W$ps, treat = lalonde$treat)
```
