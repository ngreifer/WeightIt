# Compute effective sample size of weighted sample

Computes the effective sample size (ESS) of a weighted sample, which
represents the size of an unweighted sample with approximately the same
amount of precision as the weighted sample under consideration.

## Usage

``` r
ESS(w)
```

## Arguments

- w:

  a vector of weights.

## Details

The ESS is calculated as \\(\sum w)^2/\sum w^2\\. It is invariant to
multiplicative scaling of the weights (i.e., multiplying all weights by
a nonzero scalar).

## References

McCaffrey, D. F., Ridgeway, G., & Morral, A. R. (2004). Propensity Score
Estimation With Boosted Regression for Evaluating Causal Effects in
Observational Studies. *Psychological Methods*, 9(4), 403–425.
[doi:10.1037/1082-989X.9.4.403](https://doi.org/10.1037/1082-989X.9.4.403)

Shook‐Sa, B. E., & Hudgens, M. G. (2020). Power and sample size for
observational studies of point exposure effects. *Biometrics*,
biom.13405. [doi:10.1111/biom.13405](https://doi.org/10.1111/biom.13405)

## See also

[`summary.weightit()`](https://ngreifer.github.io/WeightIt/reference/summary.weightit.md)

## Examples

``` r
library("cobalt")
#>  cobalt (Version 4.6.2, Build Date: 2026-01-29)
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "glm", estimand = "ATE"))
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATE
#>  - covariates: age, educ, married, nodegree, re74

summary(W1)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                  Max
#> treated 1.556  |--------------------------| 73.332
#> control 1.022 ||                             3.044
#> 
#> - Units with the 5 most extreme weights by group:
#>                                            
#>             184    182    181    172    124
#>  treated 11.228 11.344 12.085 26.178 73.332
#>             410    226    224    111     84
#>  control   2.33  2.437    2.5  2.637  3.044
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       1.609 0.555   0.403       0
#> control       0.247 0.211   0.029       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.    185.  
#> Weighted    404.35   51.73

ESS(W1$weights[W1$treat == 0])
#> [1] 404.3484
ESS(W1$weights[W1$treat == 1])
#> [1] 51.73462
```
