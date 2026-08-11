# Weighting from Supplied Propensity Scores

This page explains the details of computing weights from propensity
scores that have already been estimated, requested by supplying them to
the `ps` argument of
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
instead of setting `method`. This can be used with binary,
multi-category, and continuous treatments, as well as with censoring.

No model is fit; the supplied scores are converted to weights exactly as
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
would, and the result is a full `weightit` object rather than a bare
vector of weights. `formula` must still contain the treatment variable,
but the covariates on the right hand side play no role in the
computation. Any `method` supplied alongside `ps` is ignored unless it
is a user-supplied function.

Because nothing is estimated, no M-estimation components are produced,
and
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
will treat the weights as fixed.

### Binary Treatments

`ps` may be a numeric vector of the probability of being in the treated
group, a one-column matrix or data frame of the same, or a two-column
matrix or data frame with one column per treatment level. When two
columns are supplied and their names match the treatment levels, those
names determine which column is which; otherwise the columns are assumed
to be in the order of the treatment levels. All estimands allowed by
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
are available: `"ATE"`, `"ATT"`, `"ATC"`, `"ATO"`, `"ATM"`, and
`"ATOS"`.

### Multi-Category Treatments

`ps` may be a matrix or data frame with one column per treatment level,
or a numeric vector (or one-column matrix) giving each unit's
probability of being in the treatment group it is actually in. The
estimands `"ATE"`, `"ATT"`, `"ATC"`, `"ATO"`, and `"ATM"` are available.

### Continuous Treatments

`ps` is interpreted as the conditional mean of the treatment (i.e., the
fitted values of a treatment model) rather than as a probability. The
generalized propensity score is formed by feeding the standardized
residuals through the requested `density`, exactly as in
[`method_glm`](https://ngreifer.github.io/WeightIt/reference/method_glm.md).

### Longitudinal Treatments

For longitudinal treatments, `ps` is supplied as a list with one entry
per time point, and the weights are the product of the time-specific
weights.

### Censoring Weights

Wrapping the censoring indicator in
[`.cens()`](https://ngreifer.github.io/WeightIt/reference/dot-cens.md)
requests inverse probability of censoring weights. `ps` is then the
probability of *being censored*; the weights are `1/P(C = 0 | X)` for
the units still under observation and exactly 0 for the censored units.

### Sampling Weights

Sampling weights are supported and are applied to the resulting weights,
but they play no part in computing them, since no model is fit.

## Additional Arguments

For binary and multi-category treatments, the following additional
argument can be specified:

- `subclass`:

  `integer`; the number of subclasses to use for computing weights using
  marginal mean weighting through stratification (MMWS). If `NULL`,
  standard inverse probability weights (and their extensions) will be
  computed; if a number greater than 1, subclasses will be formed and
  weights will be computed based on subclass membership. See
  [`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
  for details and references.

For continuous treatments, the following additional arguments may be
supplied:

- `density`:

  A function corresponding to the conditional density of the treatment.
  The standardized residuals of the treatment will be fed through this
  function to produce the denominator of the generalized propensity
  score weights. If blank,
  [`dnorm()`](https://rdrr.io/r/stats/Normal.html) is used as
  recommended by Robins et al. (2000). This can also be supplied as a
  string containing the name of the function to be called. If the string
  contains underscores, the call will be split by the underscores and
  the latter splits will be supplied as arguments to the second argument
  and beyond. For example, if `density = "dt_2"` is specified, the
  density used will be that of a t-distribution with 2 degrees of
  freedom. Using a t-distribution can be useful when extreme treatment
  values are observed (Naimi et al., 2014).

  Can also be `"kernel"` to use kernel density estimation, which calls
  [`density()`](https://rdrr.io/r/stats/density.html) to estimate the
  denominator density for the weights. (This used to be requested by
  setting `use.kernel = TRUE`, which is now deprecated.)

- `bw`, `adjust`, `kernel`, `n`:

  If `density = "kernel"`, the arguments to
  [`density()`](https://rdrr.io/r/stats/density.html). The defaults are
  the same as those in
  [`density()`](https://rdrr.io/r/stats/density.html).

## Additional Outputs

- `obj`:

  No fit object is produced, since no model is fit; `include.obj` has no
  effect.

## References

See
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
for references on each estimand and on marginal mean weighting through
stratification.

Naimi, A. I., Moodie, E. E. M., Auger, N., & Kaufman, J. S. (2014).
Constructing Inverse Probability Weights for Continuous Exposures: A
Comparison of Methods. *Epidemiology*, 25(2), 292–299.
[doi:10.1097/EDE.0000000000000053](https://doi.org/10.1097/EDE.0000000000000053)

Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural
Models and Causal Inference in Epidemiology. *Epidemiology*, 11(5),
550–560.

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)

[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md),
which performs the same computation but returns only the weights

[`method_glm`](https://ngreifer.github.io/WeightIt/reference/method_glm.md),
for estimating the propensity scores within
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
instead of supplying them

## Examples

``` r
data("lalonde", package = "cobalt")

# Estimate the propensity score separately, then supply it
fit <- glm(treat ~ age + educ + married + nodegree + re74,
           data = lalonde, family = binomial)

(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                ps = fitted(fit), estimand = "ATT"))
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74

summary(W1)
#>                   Summary of weights
#> 
#> ─ Weight ranges:
#> 
#>           Min                                 Max
#> Treated 1.                    │             1.   
#> Control 0.022 ╞═══════════════════════════╡ 2.044
#> 
#> ─ Units with the 5 most extreme weights by group:
#>                                    
#>             5     4   3     2     1
#>  Treated    1     1   1     1     1
#>           411   595 269   409   296
#>  Control 1.33 1.437 1.5 1.637 2.044
#> 
#> ─ Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> Treated       0.    0.       0.         0
#> Control       0.823 0.701    0.33       0
#> 
#> ─ Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    255.99     185

# Marginal mean weighting through stratification
(W2 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                ps = fitted(fit), estimand = "ATT",
                subclass = 20))
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74

summary(W2)
#>                   Summary of weights
#> 
#> ─ Weight ranges:
#> 
#>           Min                               Max
#> Treated 1.          │                         1
#> Control 0.083 ╞═══════════════════════════╡   8
#> 
#> ─ Units with the 5 most extreme weights by group:
#>                             
#>            5   4   3   2   1
#>  Treated   1   1   1   1   1
#>          296 269 612 589 264
#>  Control 2.4 2.4   4   4   8
#> 
#> ─ Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> Treated       0.    0.      0.          0
#> Control       1.367 0.747   0.499       0
#> 
#> ─ Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    149.75     185
```
