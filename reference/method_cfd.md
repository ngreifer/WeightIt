# Characteristic Function Distance Balancing

This page explains the details of estimating weights using
characteristic function distance (CFD) balancing by setting
`method = "cfd"` in the call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
This method can be used with binary, multi-category, and continuous
treatments.

In general, this method relies on estimating weights by minimizing a
scalar measure of covariate balance, the CFD. The CFD is related to the
maximum mean discrepancy and captures covariate balance for the joint
covariate distribution as determined by a specific choice of kernel.
This method relies on code written for WeightIt using
[`osqp::osqp()`](https://rdrr.io/pkg/osqp/man/osqp.html) from the
[osqp](https://CRAN.R-project.org/package=osqp) package to perform the
optimization. This method may be slow or memory-intensive for large
datasets.

### Binary Treatments

For binary treatments, this method estimates the weights using `osqp()`
using formulas described by Santra, Chen, and Park (2026). The following
estimands are allowed: ATE, ATT, and ATC.

### Multi-Category Treatments

For multi-category treatments, this method estimates the weights using
`osqp()` using formulas described by Santra, Chen, and Park (2026). The
following estimands are allowed: ATE and ATT.

### Continuous Treatments

CFD balancing is not compatible with continuous treatments.

### Longitudinal Treatments

For longitudinal treatments, the weights are the product of the weights
estimated at each time point. This method is not guaranteed to yield
optimal balance at each time point. **NOTE: the use of CFD balancing
with longitudinal treatments has not been validated!**

### Sampling Weights

Sampling weights are supported through `s.weights` in all scenarios. In
some cases, sampling weights will cause the optimization to fail due to
lack of convexity or infeasible constraints.

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

CFD balancing is a method of estimating weights using optimization
without a propensity score. The weights are the solution to a
constrained quadratic optimization problem where the objective function
concerns covariate balance as measured by the CFD between groups.

CFD balancing for binary and multi-category treatments involves
minimizing the CFD between the treatment groups and between each
treatment group and a target group (e.g., the full sample for the ATE).
The CFD is a scalar measure of the difference between two multivariate
distributions. The performance of CFD balance depends on the choice of
kernel, controlled by the `kernel` argument. Each kernel corresponds to
different assumptions about the form of the true outcome model. See
Santra et al. (2026) for a comparison of these different kernels.
Setting `kernel = "energy"` is equivalent to entropy balancing, which
can also be requested by using
[`method = "energy"`](https://ngreifer.github.io/WeightIt/reference/method_energy.md).

The primary benefit of CFD balancing is that all features of the
covariate distribution are balanced, not just means, as with other
optimization-based methods like entropy balancing. Still, it is possible
to add additional balance constraints to require balance on individual
terms using the `moments` argument, just like with entropy balancing.
CFD balancing can sometimes yield weights with high variability; the
`lambda` argument can be supplied to penalize highly variable weights to
increase the effective sample size at the expense of balance.

### Reproducibility

Although there are no stochastic components to the optimization, a
feature turned off by default is to update the optimization based on how
long the optimization has been running, which will vary across runs even
when a seed is set and no parameters have been changed. See the
discussion [here](https://github.com/osqp/osqp-r/issues/19) for more
details. To ensure reproducibility by default, `adaptive_rho_interval`
is set to 10. See
[`osqp::osqpSettings()`](https://rdrr.io/pkg/osqp/man/osqpSettings.html)
for details.

## Note

Sometimes the optimization can fail to converge because the problem is
not convex. A warning will be displayed if so. In these cases, try
simply re-fitting the weights without changing anything (but see the
*Reproducibility* section above). If the method repeatedly fails, you
should try another method or change the supplied parameters (though this
is uncommon). Increasing `max_iter` or changing `adaptive_rho_interval`
might help.

If it seems like the weights are balancing the covariates but you still
get a failure to converge, this usually indicates that more iterations
are needs to find the optimal solutions. This can occur when `moments`
or `int` are specified. `max_iter` should be increased, and setting
`verbose = TRUE` allows you to monitor the process and examine if the
optimization is approaching convergence.

If `min.w` is positive and you still get a warning about the presence of
negative weights, try setting `eps` to a smaller number (e.g., to
`1e-12`).

As of version 1.5.0, `polish` is now set to `TRUE` by default. This
should yield slightly improved solutions but may be a little slower.

## Additional Arguments

The following following additional arguments can be specified:

- `kernel`:

  the name of the kernel used to characterize the CFD. Allowable optiosn
  include `"gaussian"` for the multivariate Gaussian kernel (the
  default), `"matern"` for the multivariate Matern kernel, `"energy"`
  for the energy distance kernel, `"laplace"` for the univariate
  Laplacian kernel, and `"t"` for the univariate t-dsitribution kernel.

- `nu`:

  for `kernel = "matern"`, the \\\nu\\ parameter used to control
  smoothness. The default value is 3/2. For any values other than 1/2,
  3/2, and 5/2, the GPBayes package is required to compute the Matern
  kernel. For `kernel = "t"`, the degrees of freedom for the univariate
  t-distributions used in the kernel. The default value is 5. Ignored
  for other kernels.

- `nsim`:

  for `kernel = "t"`, the number of simulations to use to compute the
  t-distribution kernel. Default is 5000. Greater is better but takes
  longer and uses more memory.

- `lambda`:

  a positive numeric scalar used to penalize the square of the weights.
  This value divided by the square of the total sample size is added to
  the diagonal of the quadratic part of the loss function. Higher values
  favor weights with less variability. Default is .0001, which is
  essentially 0.

- `moments`:

  `integer`; the highest power of each covariate to be balanced. For
  example, if `moments = 3`, each covariate, its square, and its cube
  will be balanced. Can also be a named vector with a value for each
  covariate (e.g., `moments = c(x1 = 2, x2 = 4)`). Values greater than 1
  for categorical covariates are ignored. Default is 0 to impose no
  constraint on balance.

- `int`:

  `logical`; whether first-order interactions of the covariates are to
  be balanced. Default is `FALSE`.

- `tols`:

  when `moments` is positive, a number corresponding to the maximum
  allowed standardized mean difference (for binary and multi-category
  treatments) or treatment-covariate correlation (for continuous
  treatments) allowed. Default is 0. Ignored when `moments = 0`.

- `min.w`:

  the minimum allowable weight. Negative values (including `-Inf`) are
  allowed. Default is `1e-8`.

For binary and multi-category treatments, the following additional
arguments can be specified:

- `improved`:

  `logical`; whether to include an additional term in the CFD objective
  function to minimize the distance between pairs of groups when
  `estimand = "ATE"`. Default is `TRUE`.

- `quantile`:

  a named list of quantiles (values between 0 and 1) for each continuous
  covariate, which are used to create additional variables that when
  balanced ensure balance on the corresponding quantile of the variable.
  For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures
  the 25th, 50th, and 75th percentiles of `x1` in each treatment group
  will be balanced in the weighted sample. Can also be a single number
  (e.g., `.5`) or a vector (e.g., `c(.25, .5, .75)`) to request the same
  quantile(s) for all continuous covariates.

The `moments` argument functions differently for `method = "cfd"` from
how it does with some other methods. When unspecified or set to zero,
CFD balancing weights are estimated as described by Santra et al. (2026)
for binary and multi-category treatments. When `moments` is set to an
integer larger than 0, additional balance constraints on the requested
moments of the covariates are also included, guaranteeing exact moment
balance on these covariates while minimizing the CFD of the weighted
sample. This involves exact balance on the means of the entered
covariates. The constraint on exact balance can be relaxed using the
`tols` argument.

Any other arguments will be passed to
[`osqp::osqpSettings()`](https://rdrr.io/pkg/osqp/man/osqpSettings.html)
. Some defaults differ from those in `osqpSettings()`; see
*Reproducibility* section.

## Additional Outputs

- `obj`:

  When `include.obj = TRUE`, the output of the call to
  [`osqp::solve_osqp()`](https://rdrr.io/pkg/osqp/man/solve_osqp.html) ,
  which contains the dual variables and convergence information.

## References

Santra, D., Chen, G., & Park, C. (2026). Distributional Balancing for
Causal Inference: A Unified Framework via Characteristic Function
Distance (arXiv:2601.15449). arXiv.
[doi:10.48550/arXiv.2601.15449](https://doi.org/10.48550/arXiv.2601.15449)

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)

## Examples

``` r
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cfd", estimand = "ATE"))
#> A weightit object
#>  - method: "cfd" (characteristic function distance balancing)
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
#>         Min                                  Max
#> treated   0 |------------------|          11.685
#> control   0 |---------------------------| 16.628
#> 
#> - Units with the 5 most extreme weights by group:
#>                                            
#>             184    176    172    131     82
#>  treated  8.282  8.608  9.163 10.923 11.685
#>             423    372    136    118     37
#>  control 10.134 10.217 10.518 14.278 16.628
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       2.169 1.453   1.534       0
#> control       1.972 1.274   1.262       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted   429.   185.  
#> Weighted      87.9   32.57

cobalt::bal.tab(W1)
#> Balance Measures
#>             Type Diff.Adj
#> age      Contin.  -0.0308
#> educ     Contin.   0.0019
#> married   Binary  -0.0080
#> nodegree  Binary   0.0037
#> re74     Contin.  -0.0222
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted   429.   185.  
#> Adjusted      87.9   32.57

#Using a different kernel:
(W1b <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cfd", estimand = "ATE",
                kernel = "matern", nu = 5/2))
#> Warning: Some weights are negative; these cannot be used in most model fitting
#> functions.
#> A weightit object
#>  - method: "cfd" (characteristic function distance balancing)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATE
#>  - covariates: age, educ, married, nodegree, re74

summary(W1b)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>         Min                                  Max
#> treated  -0 |---------------------------| 11.213
#> control  -0 |----------------------|       9.265
#> 
#> - Units with the 5 most extreme weights by group:
#>                                        
#>            185   184   181    82     44
#>  treated 6.945 7.517 7.836 8.774 11.213
#>            423   326   136   118     37
#>  control 6.947 8.895  9.04 9.212  9.265
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       2.068 1.445   1.339       0
#> control       1.449 0.988   0.778       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.    185.  
#> Weighted    138.59   35.21

cobalt::bal.tab(W1b)
#> Balance Measures
#>             Type Diff.Adj
#> age      Contin.  -0.0088
#> educ     Contin.   0.0007
#> married   Binary  -0.0027
#> nodegree  Binary   0.0013
#> re74     Contin.  -0.0061
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    138.59   35.21

#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cfd", estimand = "ATT",
                focal = "black"))
#> Warning: Some weights are negative; these cannot be used in most model fitting
#> functions.
#> A weightit object
#>  - method: "cfd" (characteristic function distance balancing)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 3-category (black, hispan, white)
#>  - estimand: ATT (focal: black)
#>  - covariates: age, educ, married, nodegree, re74

summary(W2)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>        Min                                  Max
#> black    1   ||                           1.   
#> hispan   0 |-----------------|           10.955
#> white   -0 |---------------------------| 17.069
#> 
#> - Units with the 5 most extreme weights by group:
#>                                           
#>              1      2      3      4      5
#>   black      1      1      1      1      1
#>             69     59     48     28      2
#>  hispan  5.033  5.133  5.336   6.45 10.955
#>             17     15      7      4      3
#>   white 11.269 12.203 13.019 16.934 17.069
#> 
#> - Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.000 0.000   0.000       0
#> hispan       1.963 1.356   1.337       0
#> white        2.485 1.446   1.268       1
#> 
#> - Effective Sample Sizes:
#> 
#>            black hispan white
#> Unweighted   243  72.   299. 
#> Weighted     243  15.01  41.8

cobalt::bal.tab(W2)
#> Balance summary across all treatment pairs
#>             Type Max.Diff.Adj
#> age      Contin.       0.0842
#> educ     Contin.       0.0427
#> married   Binary       0.0042
#> nodegree  Binary       0.0098
#> re74     Contin.       0.0081
#> 
#> Effective sample sizes
#>            hispan white black
#> Unadjusted  72.   299.    243
#> Adjusted    15.01  41.8   243
```
