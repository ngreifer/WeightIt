# Inverse Probability Tilting

This page explains the details of estimating weights using inverse
probability tilting by setting `method = "ipt"` in the call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
This method can be used with binary and multi-category treatments.

In general, this method relies on estimating propensity scores using a
modification of the usual generalized linear model score equations to
enforce balance and then converting those propensity scores into weights
using a formula that depends on the desired estimand. This method relies
on code written for WeightIt using
[`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html)
.

### Binary Treatments

For binary treatments, this method estimates the weights using formulas
described by Graham, Pinto, and Egel (2012). The following estimands are
allowed: ATE, ATT, and ATC. When the ATE is requested, the optimization
is run twice, once for each treatment group.

### Multi-Category Treatments

For multi-category treatments, this method estimates the weights using
modifications of the formulas described by Graham, Pinto, and Egel
(2012). The following estimands are allowed: ATE and ATT. When the ATE
is requested, estimation is performed once for each treatment group.
When the ATT is requested, estimation is performed once for each
non-focal (i.e., control) group.

### Continuous Treatments

Inverse probability tilting is not compatible with continuous
treatments.

### Longitudinal Treatments

For longitudinal treatments, the weights are the product of the weights
estimated at each time point. This method is not guaranteed to yield
exact balance at each time point. **NOTE: the use of inverse probability
tilting with longitudinal treatments has not been validated!**

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

M-estimation is supported for all scenarios. See
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
and
[`vignette("estimating-effects")`](https://ngreifer.github.io/WeightIt/articles/estimating-effects.md)
for details.

## Details

Inverse probability tilting (IPT) involves specifying estimating
equations that fit the parameters of two or more generalized linear
models with a modification that ensures exact balance on the covariate
means. These estimating equations are solved, and the estimated
parameters are used in the (generalized) propensity score, which is used
to compute the weights. Conceptually and mathematically, IPT is very
similar to entropy balancing and just-identified CBPS. For the ATT and
ATC, entropy balancing, just-identified CBPS, and IPT will yield
identical results. For the ATE or when `link` is specified as something
other than `"logit"`, the three methods differ.

Treatment effect estimates for binary treatments are consistent if the
true propensity score is a logistic regression or the outcome model is
linear in the covariates and their interaction with treatments. For
entropy balancing, this is only true for the ATT, and for
just-identified CBPS, this is only true if there is no effect
modification by covariates. In this way, IPT provides additional
theoretical guarantees over the other two methods, though potentially
with some cost in precision.

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
  quantile(s) for all continuous covariates.

- `link`:

  the link used to determine the inverse link for computing the
  (generalized) propensity scores. Default is `"logit"`, which is used
  in the original description of the method by Graham, Pinto, and Egel
  (2012), but `"probit"`, `"cauchit"`, `"cloglog"`, `"loglog"`, `"log"`,
  and `"clog"` are also allowed. Note that negative weights are possible
  with these last two and they should be used with caution. An object of
  class `"link-glm"` can also be supplied. The argument is passed to
  [`quasibinomial()`](https://rdrr.io/r/stats/family.html).

The `stabilize` argument is ignored.

## Additional Outputs

- `obj`:

  When `include.obj = TRUE`, the output of the call to
  [`optim()`](https://rdrr.io/r/stats/optim.html), which contains the
  coefficient estimates and convergence information. For ATE fits or
  with multi-category treatments, a list of
  [`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html)
  outputs, one for each weighted group.

## References

### `estimand = "ATE"`

Graham, B. S., De Xavier Pinto, C. C., & Egel, D. (2012). Inverse
Probability Tilting for Moment Condition Models with Missing Data. *The
Review of Economic Studies*, 79(3), 1053–1079.
[doi:10.1093/restud/rdr047](https://doi.org/10.1093/restud/rdr047)

### `estimand = "ATT"`

Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust
difference-in-differences estimators. *Journal of Econometrics*, 219(1),
101–122.
[doi:10.1016/j.jeconom.2020.06.003](https://doi.org/10.1016/j.jeconom.2020.06.003)

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)

[method_ebal](https://ngreifer.github.io/WeightIt/reference/method_ebal.md)
and
[method_cbps](https://ngreifer.github.io/WeightIt/reference/method_cbps.md)
for entropy balancing and CBPS, which work similarly.

## Examples

``` r
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ipt", estimand = "ATT"))
#> A weightit object
#>  - method: "ipt" (inverse probability tilting)
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
#> treated 1.                 ||               1.   
#> control 0.017 |---------------------------| 2.263
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>              1     2     3     4     5
#>  treated     1     1     1     1     1
#>            410   404   224   111    84
#>  control 1.464 1.485 1.576 1.743 2.263
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.000 0.000   0.000       0
#> control       0.839 0.707   0.341       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    252.12     185

cobalt::bal.tab(W1)
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.0164
#> age         Contin.  -0.0000
#> educ        Contin.   0.0000
#> married      Binary  -0.0000
#> nodegree     Binary  -0.0000
#> re74        Contin.  -0.0000
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    252.12     185

#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ipt", estimand = "ATE"))
#> A weightit object
#>  - method: "ipt" (inverse probability tilting)
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
#> black  1.57   |------------|               15.787
#> hispan 1.711  |--------------------------| 29.07 
#> white  1.102 |--|                           4.693
#> 
#> - Units with the 5 most extreme weights by group:
#>                                          
#>            203    166    163   153    152
#>   black  6.567   6.77  7.096 9.976 15.787
#>             67     43     39    36     28
#>  hispan 17.436 21.673 23.033 24.23  29.07
#>            291    285    258   205      6
#>   white  3.841  3.912  3.934 4.177  4.693
#> 
#> - Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.618 0.398   0.133       0
#> hispan       0.618 0.442   0.164       0
#> white        0.389 0.316   0.070       0
#> 
#> - Effective Sample Sizes:
#> 
#>             black hispan  white
#> Unweighted 243.     72.  299.  
#> Weighted   176.11   52.3 259.76

cobalt::bal.tab(W2)
#> Balance summary across all treatment pairs
#>             Type Max.Diff.Adj
#> age      Contin.            0
#> educ     Contin.            0
#> married   Binary            0
#> nodegree  Binary            0
#> re74     Contin.            0
#> 
#> Effective sample sizes
#>             black hispan  white
#> Unadjusted 243.     72.  299.  
#> Adjusted   176.11   52.3 259.76
```
