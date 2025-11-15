# Nonparametric Covariate Balancing Propensity Score Weighting

This page explains the details of estimating weights from nonparametric
covariate balancing propensity scores by setting `method = "npcbps"` in
the call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
This method can be used with binary, multi-category, and continuous
treatments.

In general, this method relies on estimating weights by maximizing the
empirical likelihood of the data subject to balance constraints. This
method relies on
[`CBPS::npCBPS()`](https://rdrr.io/pkg/CBPS/man/npCBPS.html) from the
[CBPS](https://CRAN.R-project.org/package=CBPS) package.

### Binary Treatments

For binary treatments, this method estimates the weights using
[`CBPS::npCBPS()`](https://rdrr.io/pkg/CBPS/man/npCBPS.html). The ATE is
the only estimand allowed. The weights are taken from the output of the
`npCBPS` fit object.

### Multi-Category Treatments

For multi-category treatments, this method estimates the weights using
[`CBPS::npCBPS()`](https://rdrr.io/pkg/CBPS/man/npCBPS.html). The ATE is
the only estimand allowed. The weights are taken from the output of the
`npCBPS` fit object.

### Continuous Treatments

For continuous treatments, this method estimates the weights using
[`CBPS::npCBPS()`](https://rdrr.io/pkg/CBPS/man/npCBPS.html). The
weights are taken from the output of the `npCBPS` fit object.

### Longitudinal Treatments

For longitudinal treatments, the weights are the product of the weights
estimated at each time point. **NOTE: the use of npCBPS with
longitudinal treatments has not been validated!**

### Sampling Weights

Sampling weights are **not** supported with `method = "npcbps"`.

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

Nonparametric CBPS involves the specification of a constrained
optimization problem over the weights. The constraints correspond to
covariate balance, and the loss function is the empirical likelihood of
the data given the weights. npCBPS is similar to [entropy
balancing](https://ngreifer.github.io/WeightIt/reference/method_ebal.md)
and will generally produce similar results. Because the optimization
problem of npCBPS is not convex it can be slow to converge or not
converge at all, so approximate balance is allowed instead using the
`cor.prior` argument, which controls the average deviation from zero
correlation between the treatment and covariates allowed.

## Additional Arguments

- `moments`:

  `integer`; the highest power of each covariate to be balanced. For
  example, if `moments = 3`, each covariate, its square, and its cube
  will be balanced. Can also be a named vector with a value for each
  covariate (e.g., `moments = c(x1 = 2, x2 = 4)`). Values greater than 1
  for categorical covariates are ignored. Default is 1 to balance
  covariate means.

- `int`:

  `logical`; whether first-order interactions of the covariates are to
  be balanced. Default is `FALSE`.

- `quantile`:

  a named list of quantiles (values between 0 and 1) for each continuous
  covariate, which are used to create additional variables that when
  balanced ensure balance on the corresponding quantile of the variable.
  For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures
  the 25th, 50th, and 75th percentiles of `x1` in each treatment group
  will be balanced in the weighted sample. Can also be a single number
  (e.g., `.5`) or a vector (e.g., `c(.25, .5, .75)`) to request the same
  quantile(s) for all continuous covariates. Only allowed with binary
  and multi-category treatments.

All arguments to `npCBPS()` can be passed through
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).

All arguments take on the defaults of those in `npCBPS()`.

## Additional Outputs

- `obj`:

  When `include.obj = TRUE`, the nonparametric CB(G)PS model fit. The
  output of the call to
  [`CBPS::npCBPS()`](https://rdrr.io/pkg/CBPS/man/npCBPS.html).

## References

Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity
score for a continuous treatment: Application to the efficacy of
political advertisements. *The Annals of Applied Statistics*, 12(1),
156–177. [doi:10.1214/17-AOAS1101](https://doi.org/10.1214/17-AOAS1101)

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md),
[`method_cbps`](https://ngreifer.github.io/WeightIt/reference/method_cbps.md)

[`method_optweight`](https://ngreifer.github.io/WeightIt/reference/method_optweight.md),
which can also be used to perform npCBPS by setting `norm = "log"`. In
generally, this `"optweight"` implementation is more stable and
flexible.

[`CBPS::npCBPS()`](https://rdrr.io/pkg/CBPS/man/npCBPS.html) for the
fitting function

## Examples

``` r
# Examples take a long time to run
data("lalonde", package = "cobalt")
# \donttest{
  #Balancing covariates between treatment groups (binary)
  (W1 <- weightit(treat ~ age + educ + married +
                    nodegree + re74, data = lalonde,
                  method = "npcbps", estimand = "ATE"))
#> A weightit object
#>  - method: "npcbps" (non-parametric covariate balancing propensity score weighting)
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
#>           Min                                 Max
#> treated 0.559 |---------------------------| 9.886
#> control 0.559 |---|                         2.129
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>            182   181   172    69    58
#>  treated 3.363 4.199 8.369  8.44 9.886
#>            410   226   224   111    84
#>  control 1.645 1.663 1.741 1.825 2.129
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       1.143 0.512   0.302       0
#> control       0.269 0.230   0.035       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.    185.  
#> Weighted    400.19   80.48

  cobalt::bal.tab(W1)
#> Balance Measures
#>             Type Diff.Adj
#> age      Contin.   0.0295
#> educ     Contin.  -0.0135
#> married   Binary   0.0407
#> nodegree  Binary  -0.0121
#> re74     Contin.   0.0704
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    400.19   80.48

  #Balancing covariates with respect to race (multi-category)
  (W2 <- weightit(race ~ age + educ + married +
                    nodegree + re74, data = lalonde,
                  method = "npcbps", estimand = "ATE"))
#> A weightit object
#>  - method: "npcbps" (non-parametric covariate balancing propensity score weighting)
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
#>          Min                                 Max
#> black  0.643  |--------------------------| 9.266
#> hispan 0.286 |---------------|             5.499
#> white  0.469  |----|                       2.51 
#> 
#> - Units with the 5 most extreme weights by group:
#>                                      
#>           241   166   163   153   152
#>   black 2.531 2.722 3.066 4.479 9.266
#>            67    43    39    36    28
#>  hispan 1.998 2.245 3.077 4.286 5.499
#>           291   285   258   205     6
#>   white 1.949 1.968 2.015  2.13  2.51
#> 
#> - Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.747 0.412   0.159       0
#> hispan       0.818 0.461   0.219       0
#> white        0.413 0.331   0.079       0
#> 
#> - Effective Sample Sizes:
#> 
#>             black hispan  white
#> Unweighted 243.    72.   299.  
#> Weighted   156.22  43.38 255.47

  cobalt::bal.tab(W2)
#> Balance summary across all treatment pairs
#>             Type Max.Diff.Adj
#> age      Contin.       0.0314
#> educ     Contin.       0.0415
#> married   Binary       0.0215
#> nodegree  Binary       0.0153
#> re74     Contin.       0.0397
#> 
#> Effective sample sizes
#>             black hispan  white
#> Unadjusted 243.    72.   299.  
#> Adjusted   156.22  43.38 255.47
# }
```
