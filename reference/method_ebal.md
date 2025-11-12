# Entropy Balancing

This page explains the details of estimating weights using entropy
balancing by setting `method = "ebal"` in the call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
This method can be used with binary, multi-category, and continuous
treatments.

In general, this method relies on estimating weights by minimizing the
negative entropy of the weights subject to exact moment balancing
constraints. This method relies on code written for WeightIt using
[`optim()`](https://rdrr.io/r/stats/optim.html).

### Binary Treatments

For binary treatments, this method estimates the weights using
[`optim()`](https://rdrr.io/r/stats/optim.html) using formulas described
by Hainmueller (2012). The following estimands are allowed: ATE, ATT,
and ATC. When the ATE is requested, the optimization is run twice, once
for each treatment group.

### Multi-Category Treatments

For multi-category treatments, this method estimates the weights using
[`optim()`](https://rdrr.io/r/stats/optim.html). The following estimands
are allowed: ATE and ATT. When the ATE is requested,
[`optim()`](https://rdrr.io/r/stats/optim.html) is run once for each
treatment group. When the ATT is requested,
[`optim()`](https://rdrr.io/r/stats/optim.html) is run once for each
non-focal (i.e., control) group.

### Continuous Treatments

For continuous treatments, this method estimates the weights using
[`optim()`](https://rdrr.io/r/stats/optim.html) using formulas described
by Tübbicke (2022) and Vegetabile et al. (2021).

### Longitudinal Treatments

For longitudinal treatments, the weights are the product of the weights
estimated at each time point. This method is not guaranteed to yield
exact balance at each time point. **NOTE: the use of entropy balancing
with longitudinal treatments has not been validated and should not be
done!**

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

Entropy balancing involves the specification of an optimization problem,
the solution to which is then used to compute the weights. The
constraints of the primal optimization problem correspond to covariate
balance on the means (for binary and multi-category treatments) or
treatment-covariate covariances (for continuous treatments), positivity
of the weights, and that the weights sum to a certain value. It turns
out that the dual optimization problem is much easier to solve because
it is over only as many variables as there are balance constraints
rather than over the weights for each unit, and it is unconstrained.

Zhao and Percival (2017) found that entropy balancing for the ATT of a
binary treatment actually involves the estimation of the coefficients of
a logistic regression propensity score model but using a specialized
loss function different from that optimized with maximum likelihood.
Entropy balancing is doubly robust (for the ATT) in the sense that it is
consistent either when the true propensity score model is a logistic
regression of the treatment on the covariates or when the true outcome
model for the control units is a linear regression of the outcome on the
covariates, and it attains a semi-parametric efficiency bound when both
are true. Entropy balancing will always yield exact mean balance on the
included terms.

## Additional Arguments

- `base.weights`:

  a vector of base weights, one for each unit. These correspond to the
  base weights \$q\$ in Hainmueller (2012). The estimated weights
  minimize the Kullback entropy divergence from the base weights,
  defined as \\\sum w \log(w/q)\\, subject to exact balance constraints.
  These can be used to supply previously estimated weights so that the
  newly estimated weights retain the some of the properties of the
  original weights while ensuring the balance constraints are met.
  Sampling weights should not be passed to `base.weights` but can be
  included in a
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  call that includes `s.weights`.

- `reltol`:

  the relative tolerance for convergence of the optimization. Passed to
  the `control` argument of
  [`optim()`](https://rdrr.io/r/stats/optim.html). Default is `1e-10`.

- `maxit`:

  the maximum number of iterations for convergence of the optimization.
  Passed to the `control` argument of
  [`optim()`](https://rdrr.io/r/stats/optim.html). Default is 1000 for
  binary and multi-category treatments and 10000 for continuous and
  longitudinal treatments.

- `solver`:

  the solver to use to estimate the parameters. Allowable options
  include `"multiroot"` to use
  [`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html)
  and `"optim"` to use
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html). `"multiroot"`
  is the default when rootSolve is installed, as it tends to be much
  faster and more accurate; otherwise, `"optim"` is the default and
  requires no dependencies. Regardless of `solver`, the output of
  [`optim()`](https://rdrr.io/r/stats/optim.html) is returned when
  `include.obj = TRUE` (see below).

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

- `d.moments`:

  `integer`; with continuous treatments, the number of moments of the
  treatment and covariate distributions that are constrained to be the
  same in the weighted sample as in the original sample. For example,
  setting `d.moments = 3` ensures that the mean, variance, and skew of
  the treatment and covariates are the same in the weighted sample as in
  the unweighted sample. `d.moments` should be greater than or equal to
  `moments` and will be automatically set accordingly if not (or if not
  specified). Vegetabile et al. (2021) recommend setting
  `d.moments = 3`, even if `moments` is less than 3. This argument
  corresponds to the tuning parameters \\r\\ and \\s\\ in Vegetabile et
  al. (2021) (which here are set to be equal). Ignored for binary and
  multi-category treatments.

The `stabilize` argument is ignored; in the past it would reduce the
variability of the weights through an iterative process. If you want to
minimize the variance of the weights subject to balance constraints, use
`method = "optweight"`.

## Additional Outputs

- `obj`:

  When `include.obj = TRUE`, the output of the call to
  [`optim()`](https://rdrr.io/r/stats/optim.html), which contains the
  dual variables and convergence information. For ATE fits or with
  multi-category treatments, a list of
  [`optim()`](https://rdrr.io/r/stats/optim.html) outputs, one for each
  weighted group.

## References

### Binary Treatments

#### `estimand = "ATT"`

Hainmueller, J. (2012). Entropy Balancing for Causal Effects: A
Multivariate Reweighting Method to Produce Balanced Samples in
Observational Studies. *Political Analysis*, 20(1), 25–46.
[doi:10.1093/pan/mpr025](https://doi.org/10.1093/pan/mpr025)

Zhao, Q., & Percival, D. (2017). Entropy balancing is doubly robust.
*Journal of Causal Inference*, 5(1).
[doi:10.1515/jci-2016-0010](https://doi.org/10.1515/jci-2016-0010)

#### `estimand = "ATE"`

Källberg, D., & Waernbaum, I. (2023). Large Sample Properties of Entropy
Balancing Estimators of Average Causal Effects. *Econometrics and
Statistics*.
[doi:10.1016/j.ecosta.2023.11.004](https://doi.org/10.1016/j.ecosta.2023.11.004)

### Continuous Treatments

Tübbicke, S. (2022). Entropy Balancing for Continuous Treatments.
*Journal of Econometric Methods*, 11(1), 71–89.
[doi:10.1515/jem-2021-0002](https://doi.org/10.1515/jem-2021-0002)

Vegetabile, B. G., Griffin, B. A., Coffman, D. L., Cefalu, M., Robbins,
M. W., & McCaffrey, D. F. (2021). Nonparametric estimation of population
average dose-response curves using entropy balancing weights for
continuous exposures. *Health Services and Outcomes Research
Methodology*, 21(1), 69–110.
[doi:10.1007/s10742-020-00236-2](https://doi.org/10.1007/s10742-020-00236-2)

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)

[method_ipt](https://ngreifer.github.io/WeightIt/reference/method_ipt.md)
and
[method_cbps](https://ngreifer.github.io/WeightIt/reference/method_cbps.md)
for inverse probability tilting and CBPS, which work similarly.

## Examples

``` r
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATT"))
#> A weightit object
#>  - method: "ebal" (entropy balancing)
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
#> treated 1.         ||                      1.   
#> control 0.04 |---------------------------| 5.247
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>              1     2     3     4     5
#>  treated     1     1     1     1     1
#>            410   404   224   111    84
#>  control 3.396 3.443 3.655 4.043 5.247
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
#>             Type Diff.Adj
#> age      Contin.        0
#> educ     Contin.        0
#> married   Binary       -0
#> nodegree  Binary       -0
#> re74     Contin.       -0
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    252.12     185

#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATE"))
#> A weightit object
#>  - method: "ebal" (entropy balancing)
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
#> black  0.553   |-------------------------| 5.35 
#> hispan 0.141 |----------------|            3.332
#> white  0.398  |-------|                    1.923
#> 
#> - Units with the 5 most extreme weights by group:
#>                                      
#>           203   166   163   153   152
#>   black 2.521 2.549 2.806 3.555  5.35
#>            67    43    39    36    28
#>  hispan 2.046  2.53 2.632 2.705 3.332
#>           291   285   258   205     6
#>   white 1.711 1.723 1.743 1.774 1.923
#> 
#> - Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.590 0.413   0.131       0
#> hispan       0.609 0.440   0.163       0
#> white        0.371 0.306   0.068       0
#> 
#> - Effective Sample Sizes:
#> 
#>             black hispan  white
#> Unweighted 243.    72.   299.  
#> Weighted   180.47  52.71 262.93

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
#> Unadjusted 243.    72.   299.  
#> Adjusted   180.47  52.71 262.93

#Balancing covariates and squares with respect to
#re75 (continuous), maintaining 3 moments of the
#covariate and treatment distributions
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", moments = 2,
                d.moments = 3))
#> A weightit object
#>  - method: "ebal" (entropy balancing)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: continuous
#>  - covariates: age, educ, married, nodegree, re74

summary(W3)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>     Min                                 Max
#> all   0 |---------------------------| 17.65
#> 
#> - Units with the 5 most extreme weights:
#>                                    
#>        484   200   180    171   166
#>  all 6.745 7.616 8.523 10.266 17.65
#> 
#> - Weight statistics:
#> 
#>     Coef of Var   MAD Entropy # Zeros
#> all       1.252 0.614   0.423       0
#> 
#> - Effective Sample Sizes:
#> 
#>             Total
#> Unweighted 614.  
#> Weighted   239.32

cobalt::bal.tab(W3, poly = 2,
                stats = c("c", "m"))
#> Balance Measures
#>             Type Corr.Adj Diff.Target.Adj
#> age      Contin.  -0.0000              -0
#> educ     Contin.  -0.0001              -0
#> married   Binary  -0.0001               0
#> nodegree  Binary   0.0005               0
#> re74     Contin.  -0.0001               0
#> age²     Contin.  -0.0000              -0
#> educ²    Contin.  -0.0001              -0
#> re74²    Contin.  -0.0000               0
#> 
#> Effective sample sizes
#>             Total
#> Unadjusted 614.  
#> Adjusted   239.32
```
