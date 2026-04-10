# Fitting (Weighted) Multinomial Regression Models

`multinom_weightit()` fits a multinomial logistic regression model with
a covariance matrix that accounts for estimation of weights, if
supplied. By default, this function uses M-estimation to construct a
robust covariance matrix using the estimating equations for the
weighting model and the outcome model when available.

## Usage

``` r
multinom_weightit(
  formula,
  data,
  link = "logit",
  weightit = NULL,
  vcov = NULL,
  cluster,
  R = 500L,
  offset,
  start = NULL,
  control = list(...),
  x = FALSE,
  y = TRUE,
  contrasts = NULL,
  fwb.args = list(),
  ...
)
```

## Arguments

- formula:

  an object of class [`formula`](https://rdrr.io/r/stats/formula.html)
  (or one that can be coerced to that class): a symbolic description of
  the model to be fitted.

- data:

  a data frame containing the variables in the model. If not found in
  data, the variables are taken from `environment(formula)`, typically
  the environment from which the function is called.

- link:

  a string corresponding to the desired link function. Currently, only
  `"logit"` is allowed.

- weightit:

  a `weightit` or `weightitMSM` object; the output of a call to
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  or
  [`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
  If not supplied, an unweighted model will be fit.

- vcov:

  string; the method used to compute the variance of the estimated
  parameters. Allowable options include `"asympt"`, which uses the
  asymptotically correct M-estimation-based method that accounts for
  estimation of the weights when available; `"const"`, which uses the
  usual maximum likelihood estimates (only available when `weightit` is
  not supplied); `"HC0"`, which computes the robust sandwich variance
  treating weights (if supplied) as fixed; `"BS"`, which uses the
  traditional bootstrap (including re-estimation of the weights, if
  supplied); `"FWB"`, which uses the fractional weighted bootstrap as
  implemented in
  [`fwb::fwb()`](https://ngreifer.github.io/fwb/reference/fwb.html)
  (including re-estimation of the weights, if supplied); and `"none"` to
  omit calculation of a variance matrix. If `NULL` (the default), will
  use `"asympt"` if `weightit` is supplied and M-estimation is available
  and `"HC0"` otherwise. See the `vcov_type` component of the outcome
  object to see which was used.

- cluster:

  optional; for computing a cluster-robust variance matrix, a variable
  indicating the clustering of observations, a list (or data frame)
  thereof, or a one-sided formula specifying which variable(s) from the
  fitted model should be used. Note the cluster-robust variance matrix
  uses a correction for small samples, as is done in
  [`sandwich::vcovCL()`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)
  by default. Cluster-robust variance calculations are available only
  when `vcov` is `"asympt"`, `"HC0"`, `"BS"`, or `"FWB"`.

- R:

  the number of bootstrap replications when `vcov` is `"BS"` or `"FWB"`.
  Default is 500. Ignored otherwise.

- offset:

  optional; a numeric vector containing the model offset. See
  [`offset()`](https://rdrr.io/r/stats/offset.html). An offset can also
  be present in the model formula.

- start:

  optional starting values for the coefficients.

- control:

  a list of parameters for controlling the fitting process.

- x, y:

  logical values indicating whether the response vector and model matrix
  used in the fitting process should be returned as components of the
  returned value.

- contrasts:

  an optional list defining contrasts for factor variables. See
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).

- fwb.args:

  an optional list of further arguments to supply to
  [`fwb::fwb()`](https://ngreifer.github.io/fwb/reference/fwb.html) when
  `vcov = "FWB"`.

- ...:

  arguments to be used to form the default control argument if it is not
  supplied directly.

## Value

A `multinom_weightit` object.

Unless `vcov = "none"`, the `vcov` component contains the covariance
matrix adjusted for the estimation of the weights if requested and a
compatible `weightit` object was supplied. The `vcov_type` component
contains the type of variance matrix requested. If `cluster` is
supplied, it will be stored in the `"cluster"` attribute of the output
object, even if not used.

The `model` component of the output object (also the
[`model.frame()`](https://rdrr.io/r/stats/model.frame.html) output) will
include two extra columns when `weightit` is supplied: `(weights)`
containing the weights used in the model (the product of the estimated
weights and the sampling weights, if any) and `(s.weights)` containing
the sampling weights, which will be 1 if `s.weights` is not supplied in
the original
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
call.

## Details

`multinom_weightit()` implements multinomial logistic regression using a
custom function in WeightIt that optionally computes a coefficient
variance matrix that can be adjusted to account for estimation of the
weights if a `weightit` or `weightitMSM` object is supplied to the
`weightit` argument. This implementation is less robust to failures than
other multinomial logistic regression solvers and should be used with
caution. Estimation of coefficients should align with that from
[`mlogit::mlogit()`](https://rdrr.io/pkg/mlogit/man/mlogit.html) and
[`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html)
but might differ from
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html) due to
the relaxed convergence thresholds of the latter.

When no argument is supplied to `weightit` or there is no `"Mparts"`
attribute in the supplied object, the default variance matrix returned
will be the "HC0" sandwich variance matrix, which is robust to
misspecification of the outcome family (including heteroscedasticity).
Otherwise, the default variance matrix uses M-estimation to additionally
adjust for estimation of the weights. When possible, this often yields
smaller (and more accurate) standard errors. See the individual methods
pages to see whether and when an `"Mparts"` attribute is included in the
supplied object. To request that a variance matrix be computed that
doesn't account for estimation of the weights even when a compatible
`weightit` object is supplied, set `vcov = "HC0"`, which treats the
weights as fixed.

Bootstrapping can also be used to compute the coefficient variance
matrix; when `vcov = "BS"` or `vcov = "FWB"`, which implement the
traditional resampling-based and fractional weighted bootstrap,
respectively, the entire process of estimating the weights and fitting
the outcome model is repeated in bootstrap samples (if a `weightit`
object is supplied). This accounts for estimation of the weights and can
be used with any weighting method. It is important to set a seed using
[`set.seed()`](https://rdrr.io/r/base/Random.html) to ensure
replicability of the results. The fractional weighted bootstrap is more
reliable but requires the weighting method to accept sampling weights
(which most do, and you'll get an error if it doesn't). Setting
`vcov = "FWB"` and supplying `fwb.args = list(wtype = "multinom")` also
performs the resampling-based bootstrap but with the additional features
fwb provides (e.g., a progress bar and parallelization).

## See also

- [`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
  for fitting generalized linear models that adjust for estimation of
  the weights.

- [`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/ordinal_weightit.md)
  for fitting ordinal regression models that adjust for estimation of
  the weights.

- [`coxph_weightit()`](https://ngreifer.github.io/WeightIt/reference/coxph_weightit.md)
  for fitting Cox proportional hazards models that adjust for estimation
  of the weights.

- [`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html)
  for fitting multinomial regression models that do not account for
  estimation of the weights.

## Examples

``` r
data("lalonde", package = "cobalt")

# Logistic regression ATT weights
w.out <- weightit(treat ~ age + educ + married + re74,
                  data = lalonde, method = "glm",
                  estimand = "ATT")

# Multinomial logistic regression outcome model
# that adjusts for estimation of weights
lalonde$re78_3 <- factor(findInterval(lalonde$re78,
                                      c(0, 5e3, 1e4)))

fit <- multinom_weightit(re78_3 ~ treat,
                         data = lalonde,
                         weightit = w.out)

summary(fit)
#> 
#> Call:
#> multinom_weightit(formula = re78_3 ~ treat, data = lalonde, weightit = w.out)
#> 
#> Coefficients:
#>               Estimate Std. Error z value Pr(>|z|)    
#> 2~(Intercept) -0.90398    0.14906  -6.064   <1e-06 ***
#> 2~treat        0.05006    0.23633   0.212    0.832    
#> 3~(Intercept) -1.02170    0.15699  -6.508   <1e-06 ***
#> 3~treat        0.12015    0.23835   0.504    0.614    
#> Standard error: HC0 robust (adjusted for estimation of weights)
#> 
```
