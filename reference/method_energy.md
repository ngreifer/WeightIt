# Energy Balancing

This page explains the details of estimating weights using energy
balancing by setting `method = "energy"` in the call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
This method can be used with binary, multi-category, and continuous
treatments.

In general, this method relies on estimating weights by minimizing an
energy statistic related to covariate balance. For binary and
multi-category treatments, this is the energy distance, which is a
multivariate distance between distributions, between the treatment
groups. For continuous treatments, this is the sum of the distance
covariance between the treatment variable and the covariates and the
energy distances between the treatment and covariates in the weighted
sample and their distributions in the original sample. This method
relies on code written for WeightIt using
[`osqp::osqp()`](https://rdrr.io/pkg/osqp/man/osqp.html) from the
[osqp](https://CRAN.R-project.org/package=osqp) package to perform the
optimization. This method may be slow or memory-intensive for large
datasets.

### Binary Treatments

For binary treatments, this method estimates the weights using `osqp()`
using formulas described by Huling and Mak (2024). The following
estimands are allowed: ATE, ATT, and ATC.

### Multi-Category Treatments

For multi-category treatments, this method estimates the weights using
`osqp()` using formulas described by Huling and Mak (2024). The
following estimands are allowed: ATE and ATT.

### Continuous Treatments

For continuous treatments, this method estimates the weights using
`osqp()` using formulas described by Huling, Greifer, and Chen (2023).

### Longitudinal Treatments

For longitudinal treatments, the weights are the product of the weights
estimated at each time point. This method is not guaranteed to yield
optimal balance at each time point. **NOTE: the use of energy balancing
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

Energy balancing is a method of estimating weights using optimization
without a propensity score. The weights are the solution to a constrain
quadratic optimization problem where the objective function concerns
covariate balance as measured by the energy distance and (for continuous
treatments) the distance covariance.

Energy balancing for binary and multi-category treatments involves
minimizing the energy distance between the treatment groups and between
each treatment group and a target group (e.g., the full sample for the
ATE). The energy distance is a scalar measure of the difference between
two multivariate distributions and is equal to 0 when the two
distributions are identical.

Energy balancing for continuous treatments involves minimizing the
distance covariance between the treatment and the covariates; the
distance covariance is a scalar measure of the association between two
(possibly multivariate) distributions that is equal to 0 when the two
distributions are independent. In addition, the energy distances between
the treatment and covariate distributions in the weighted sample and the
treatment and covariate distributions in the original sample are
minimized.

The primary benefit of energy balancing is that all features of the
covariate distribution are balanced, not just means, as with other
optimization-based methods like entropy balancing. Still, it is possible
to add additional balance constraints to require balance on individual
terms using the `moments` argument, just like with entropy balancing.
Energy balancing can sometimes yield weights with high variability; the
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

- `dist.mat`:

  the name of the method used to compute the distance matrix of the
  covariates or the numeric distance matrix itself. Allowable options
  include `"scaled_euclidean"` for the Euclidean (L2) distance on the
  scaled covariates (the default), `"mahalanobis"` for the Mahalanobis
  distance, and `"euclidean"` for the raw Euclidean distance.
  Abbreviations allowed. Note that some user-supplied distance matrices
  can cause the R session to abort due to a bug within osqp, so this
  argument should be used with caution. A distance matrix must be a
  square, symmetric, numeric matrix with zeros along the diagonal and a
  row and column for each unit. Can also be supplied as the output of a
  call to [`dist()`](https://rdrr.io/r/stats/dist.html).

- `lambda`:

  a positive numeric scalar used to penalize the square of the weights.
  This value divided by the square of the total sample size is added to
  the diagonal of the quadratic part of the loss function. Higher values
  favor weights with less variability. Note this is distinct from the
  lambda value described in Huling and Mak (2024), which penalizes the
  complexity of individual treatment rules rather than the weights, but
  does correspond to lambda from Huling et al. (2023). Default is .0001,
  which is essentially 0.

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

For binary and multi-category treatments, the following additional
arguments can be specified:

- `improved`:

  `logical`; whether to use the improved energy balancing weights as
  described by Huling and Mak (2024) when `estimand = "ATE"`. This
  involves optimizing balance not only between each treatment group and
  the overall sample, but also between each pair of treatment groups.
  Huling and Mak (2024) found that the improved energy balancing weights
  generally outperformed standard energy balancing. Default is `TRUE`;
  set to `FALSE` to use the standard energy balancing weights instead
  (not recommended).

- `quantile`:

  a named list of quantiles (values between 0 and 1) for each continuous
  covariate, which are used to create additional variables that when
  balanced ensure balance on the corresponding quantile of the variable.
  For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures
  the 25th, 50th, and 75th percentiles of `x1` in each treatment group
  will be balanced in the weighted sample. Can also be a single number
  (e.g., `.5`) or a vector (e.g., `c(.25, .5, .75)`) to request the same
  quantile(s) for all continuous covariates.

For continuous treatments, the following additional arguments can be
specified:

- `d.moments`:

  The number of moments of the treatment and covariate distributions
  that are constrained to be the same in the weighted sample as in the
  original sample. For example, setting `d.moments = 3` ensures that the
  mean, variance, and skew of the treatment and covariates are the same
  in the weighted sample as in the unweighted sample. `d.moments` should
  be greater than or equal to `moments` and will be automatically set
  accordingly if not (or if not specified).

- `dimension.adj`:

  `logical`; whether to include the dimensionality adjustment described
  by Huling et al. (2023). If `TRUE`, the default, the energy distance
  for the covariates is weighted \\\sqrt{p}\\ times as much as the
  energy distance for the treatment, where \\p\\ is the number of
  covariates. If `FALSE`, the two energy distances are given equal
  weights. Default is `TRUE`.

The `moments` argument functions differently for `method = "energy"`
from how it does with other methods. When unspecified or set to zero,
energy balancing weights are estimated as described by Huling and Mak
(2024) for binary and multi-category treatments or by Huling et al.
(2023) for continuous treatments. When `moments` is set to an integer
larger than 0, additional balance constraints on the requested moments
of the covariates are also included, guaranteeing exact moment balance
on these covariates while minimizing the energy distance of the weighted
sample. For binary and multi-category treatments, this involves exact
balance on the means of the entered covariates; for continuous
treatments, this involves exact balance on the treatment-covariate
correlations of the entered covariates. The constraint on exact balance
can be relaxed using the `tols` argument.

Any other arguments will be passed to
[`osqp::osqpSettings()`](https://rdrr.io/pkg/osqp/man/osqpSettings.html).
Some defaults differ from those in `osqpSettings()`; see
*Reproducibility* section.

## Additional Outputs

- `obj`:

  When `include.obj = TRUE`, the output of the call to
  [`osqp::solve_osqp()`](https://rdrr.io/pkg/osqp/man/solve_osqp.html),
  which contains the dual variables and convergence information.

## References

### Binary and multi-category treatments

Huling, J. D., & Mak, S. (2024). Energy balancing of covariate
distributions. *Journal of Causal Inference*, 12(1).
[doi:10.1515/jci-2022-0029](https://doi.org/10.1515/jci-2022-0029)

### Continuous treatments

Huling, J. D., Greifer, N., & Chen, G. (2023). Independence weights for
causal inference with continuous treatments. *Journal of the American
Statistical Association*, 0(ja), 1–25.
[doi:10.1080/01621459.2023.2213485](https://doi.org/10.1080/01621459.2023.2213485)

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)

## Author

Noah Greifer, using code from Jared Huling's
[independenceWeights](https://CRAN.R-project.org/package=independenceWeights)
package for continuous treatments.

## Examples

``` r
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "energy", estimand = "ATE"))
#> A weightit object
#>  - method: "energy" (energy balancing)
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
#> treated   0 |---------------|              7.527
#> control   0 |---------------------------| 12.926
#> 
#> - Units with the 5 most extreme weights by group:
#>                                        
#>            184   183   181   172    124
#>  treated 3.797 3.946 4.725 5.738  7.527
#>            423   382   189   118      7
#>  control 5.096  5.16 5.163 6.539 12.926
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       1.148 0.852   0.551       0
#> control       1.162 0.765   0.467       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.    185.  
#> Weighted    182.79   80.09

cobalt::bal.tab(W1)
#> Balance Measures
#>             Type Diff.Adj
#> age      Contin.  -0.0235
#> educ     Contin.   0.0052
#> married   Binary  -0.0055
#> nodegree  Binary   0.0006
#> re74     Contin.  -0.0061
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.    185.  
#> Adjusted    182.79   80.09

#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "energy", estimand = "ATT",
                focal = "black"))
#> A weightit object
#>  - method: "energy" (energy balancing)
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
#> hispan   0 |---------|                    4.686
#> white    0 |---------------------------| 12.674
#> 
#> - Units with the 5 most extreme weights by group:
#>                                        
#>             1     2     3      4      5
#>   black     1     1     1      1      1
#>            69    48    40     28      2
#>  hispan 3.637 3.932  4.11  4.545  4.686
#>           184    15    10      7      4
#>   white 7.642 7.659 9.391 11.166 12.674
#> 
#> - Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.000 0.000   0.000       0
#> hispan       1.165 0.879   0.618       0
#> white        1.733 1.073   0.967       0
#> 
#> - Effective Sample Sizes:
#> 
#>            black hispan  white
#> Unweighted   243   72.  299.  
#> Weighted     243   30.8  74.86

cobalt::bal.tab(W2)
#> Balance summary across all treatment pairs
#>             Type Max.Diff.Adj
#> age      Contin.       0.0389
#> educ     Contin.       0.0455
#> married   Binary       0.0094
#> nodegree  Binary       0.0027
#> re74     Contin.       0.0259
#> 
#> Effective sample sizes
#>            hispan  white black
#> Unadjusted   72.  299.     243
#> Adjusted     30.8  74.86   243

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                    nodegree + re74, data = lalonde,
                  method = "energy", moments = 1,
                  tols = .01))
#> A weightit object
#>  - method: "energy" (energy balancing)
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
#> all   0 |---------------------------| 9.578
#> 
#> - Units with the 5 most extreme weights:
#>                                  
#>        501   486  180   171   166
#>  all 6.738 6.739 7.46 9.256 9.578
#> 
#> - Weight statistics:
#> 
#>     Coef of Var   MAD Entropy # Zeros
#> all       1.194 0.809   0.553       0
#> 
#> - Effective Sample Sizes:
#> 
#>             Total
#> Unweighted 614.  
#> Weighted   253.28

cobalt::bal.tab(W3, poly = 2)
#> Balance Measures
#>             Type Corr.Adj
#> age      Contin.   0.0037
#> educ     Contin.  -0.0025
#> married   Binary   0.0100
#> nodegree  Binary  -0.0050
#> re74     Contin.   0.0100
#> age²     Contin.  -0.0072
#> educ²    Contin.   0.0096
#> re74²    Contin.  -0.0348
#> 
#> Effective sample sizes
#>             Total
#> Unadjusted 614.  
#> Adjusted   253.28
```
