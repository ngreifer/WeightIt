# Fitting Weighted Generalized Linear Models

`glm_weightit()` is used to fit generalized linear models with a
covariance matrix that accounts for estimation of weights, if supplied.
`lm_weightit()` is a wrapper for `glm_weightit()` with the Gaussian
family and identity link (i.e., a linear model). `ordinal_weightit()`
fits proportional odds ordinal regression models. `multinom_weightit()`
fits multinomial logistic regression models. `coxph_weightit()` fits a
Cox proportional hazards model and is a wrapper for
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html). By
default, these functions use M-estimation to construct a robust
covariance matrix using the estimating equations for the weighting model
and the outcome model when available.

## Usage

``` r
glm_weightit(
  formula,
  data,
  family = gaussian,
  weightit = NULL,
  vcov = NULL,
  cluster,
  R = 500L,
  offset,
  start = NULL,
  etastart,
  mustart,
  control = list(...),
  x = FALSE,
  y = TRUE,
  contrasts = NULL,
  fwb.args = list(),
  br = FALSE,
  ...
)

lm_weightit(
  formula,
  data,
  weightit = NULL,
  vcov = NULL,
  cluster,
  R = 500L,
  offset,
  x = FALSE,
  y = TRUE,
  contrasts = NULL,
  fwb.args = list(),
  ...
)

ordinal_weightit(
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

coxph_weightit(
  formula,
  data,
  weightit = NULL,
  vcov = NULL,
  cluster,
  R = 500L,
  x = FALSE,
  y = TRUE,
  fwb.args = list(),
  ...
)
```

## Arguments

- formula:

  an object of class "formula" (or one that can be coerced to that
  class): a symbolic description of the model to be fitted. For
  `coxph_weightit()`, see
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) for
  how this should be specified.

- data:

  a data frame containing the variables in the model. If not found in
  data, the variables are taken from `environment(formula)`, typically
  the environment from which the function is called.

- family:

  a description of the error distribution and link function to be used
  in the model. This can be a character string naming a family function,
  a family function or the result of a call to a family function. See
  [family](https://rdrr.io/r/stats/family.html) for details of family
  functions.

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
  be preset in the model formula.

- start:

  optional starting values for the coefficients.

- etastart, mustart:

  optional starting values for the linear predictor and vector of means.
  Passed to [`glm()`](https://rdrr.io/r/stats/glm.html).

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

- br:

  `logical`; whether to use bias-reduced regression as implemented by
  [`brglm2::brglmFit()`](https://rdrr.io/pkg/brglm2/man/brglmFit.html) .
  If `TRUE`, arguments passed to `control` or ... will be passed to
  [`brglm2::brglmControl()`](https://rdrr.io/pkg/brglm2/man/brglmControl.html)
  .

- ...:

  arguments to be used to form the default control argument if it is not
  supplied directly.

- link:

  for `ordinal_weightit()` and `multinom_weightit()`, a string
  corresponding to the desired link function. For `ordinal_weightit()`,
  any allowed by [`binomial()`](https://rdrr.io/r/stats/family.html) are
  accepted; for `multinom_weightit()`, only `"logit"` is allowed.
  Default is `"logit"` for ordinal or multinomial logistic regression,
  respectively.

## Value

For `lm_weightit()` and `glm_weightit()`, a `glm_weightit` object, which
inherits from `glm`. For `ordinal_weightit()` and `multinom_weightit()`,
an `ordinal_weightit` or `multinom_weightit` object, respectively. For
`coxph_weightit()`, a `coxph_weightit` object, which inherits from
`coxph`. See
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) for
details.

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

`glm_weightit()` is essentially a wrapper for
[`glm()`](https://rdrr.io/r/stats/glm.html) that optionally computes a
coefficient variance matrix that can be adjusted to account for
estimation of the weights if a `weightit` or `weightitMSM` object is
supplied to the `weightit` argument. When no argument is supplied to
`weightit` or there is no `"Mparts"` attribute in the supplied object,
the default variance matrix returned will be the "HC0" sandwich variance
matrix, which is robust to misspecification of the outcome family
(including heteroscedasticity). Otherwise, the default variance matrix
uses M-estimation to additionally adjust for estimation of the weights.
When possible, this often yields smaller (and more accurate) standard
errors. See the individual methods pages to see whether and when an
`"Mparts"` attribute is included in the supplied object. To request that
a variance matrix be computed that doesn't account for estimation of the
weights even when a compatible `weightit` object is supplied, set
`vcov = "HC0"`, which treats the weights as fixed.

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
fwb provides (e.g., a progress bar and parallelization) at the expense
of needing to have fwb installed.

`multinom_weightit()` implements multinomial logistic regression using a
custom function in WeightIt. This implementation is less robust to
failures than other multinomial logistic regression solvers and should
be used with caution. Estimation of coefficients should align with that
from [`mlogit::mlogit()`](https://rdrr.io/pkg/mlogit/man/mlogit.html)
and
[`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html).

`ordinal_weightit()` implements proportional odds ordinal regression
using a custom function in WeightIt. Estimation of coefficients should
align with that from
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html).

`coxph_weightit()` is a wrapper for
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) to
fit weighted survival models. It differs from `coxph()` in that the
`cluster` argument (if used) should be specified as a one-sided formula
(which can include multiple clustering variables) and uses a small
sample correction for cluster variance estimates when specified.
Currently, M-estimation is not supported, so bootstrapping (i.e.,
`vcov = "BS"` or `"FWB"`) is the only way to correctly adjust for
estimation of the weights. By default, the robust variance is estimated
treating weights as fixed, which is the same variance returned when
`robust = TRUE` in `coxph()`.

Functions in the sandwich package can be to compute standard errors
after fitting, regardless of how `vcov` was specified, though these will
ignore estimation of the weights, if any. When no adjustment is done for
estimation of the weights (i.e., because no `weightit` argument was
supplied or there was no `"Mparts"` component in the supplied object),
the default variance matrix produced by `glm_weightit()` should align
with that from `sandwich::vcovHC(. type = "HC0")` or
`sandwich::vcovCL(., type = "HC0", cluster = cluster)` when `cluster` is
supplied. Not all types are available for all models.

## See also

[`lm()`](https://rdrr.io/r/stats/lm.html) and
[`glm()`](https://rdrr.io/r/stats/glm.html) for fitting generalized
linear models without adjusting standard errors for estimation of the
weights.
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) for
fitting Cox proportional hazards models without adjusting standard
errors for estimation of the weights.

## Examples

``` r
data("lalonde", package = "cobalt")

# Logistic regression ATT weights
w.out <- weightit(treat ~ age + educ + married + re74,
                  data = lalonde, method = "glm",
                  estimand = "ATT")

# Linear regression outcome model that adjusts
# for estimation of weights
fit1 <- lm_weightit(re78 ~ treat, data = lalonde,
                    weightit = w.out)

summary(fit1)
#> 
#> Call:
#> lm_weightit(formula = re78 ~ treat, data = lalonde, weightit = w.out)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   5515.9      376.5  14.649   <1e-06 ***
#> treat          833.2      669.1   1.245    0.213    
#> Standard error: HC0 robust (adjusted for estimation of weights)
#> 

# Linear regression outcome model that treats weights
# as fixed
fit2 <- lm_weightit(re78 ~ treat, data = lalonde,
                    weightit = w.out,
                    vcov = "HC0")

summary(fit2)
#> 
#> Call:
#> lm_weightit(formula = re78 ~ treat, data = lalonde, weightit = w.out, 
#>     vcov = "HC0")
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   5515.9      353.5  15.604   <1e-06 ***
#> treat          833.2      676.6   1.232    0.218    
#> Standard error: HC0 robust
#> 

# Can also just call summary() with `vcov` option
summary(fit1, vcov = "HC0")
#> 
#> Call:
#> lm_weightit(formula = re78 ~ treat, data = lalonde, weightit = w.out)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   5515.9      353.5  15.604   <1e-06 ***
#> treat          833.2      676.6   1.232    0.218    
#> Standard error: HC0 robust
#> 
# Linear regression outcome model that bootstraps
# estimation of weights and outcome model fitting
# using fractional weighted bootstrap with "Mammen"
# weights
set.seed(123)
fit3 <- lm_weightit(re78 ~ treat, data = lalonde,
                    weightit = w.out,
                    vcov = "FWB",
                    R = 50, #should use way more
                    fwb.args = list(wtype = "mammen"))

summary(fit3)
#> 
#> Call:
#> lm_weightit(formula = re78 ~ treat, data = lalonde, weightit = w.out, 
#>     vcov = "FWB", R = 50, fwb.args = list(wtype = "mammen"))
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   5515.9      361.5  15.257   <1e-06 ***
#> treat          833.2      644.4   1.293    0.196    
#> Standard error: fractional weighted bootstrap
#> 
# Multinomial logistic regression outcome model
# that adjusts for estimation of weights
lalonde$re78_3 <- factor(findInterval(lalonde$re78,
                                      c(0, 5e3, 1e4)))

fit4 <- multinom_weightit(re78_3 ~ treat,
                          data = lalonde,
                          weightit = w.out)

summary(fit4)
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

# Ordinal probit regression that adjusts for estimation
# of weights
fit5 <- ordinal_weightit(ordered(re78_3) ~ treat,
                         data = lalonde,
                         link = "probit",
                         weightit = w.out)

summary(fit5)
#> 
#> Call:
#> ordinal_weightit(formula = ordered(re78_3) ~ treat, data = lalonde, 
#>     link = "probit", weightit = w.out)
#> 
#> Coefficients:
#>       Estimate Std. Error z value Pr(>|z|)
#> treat   0.0554     0.1114   0.498    0.619
#> Standard error: HC0 robust (adjusted for estimation of weights)
#> 
#> Thresholds:
#>     Estimate Std. Error z value Pr(>|z|)    
#> 1|2  0.16926    0.07479   2.263   0.0236 *  
#> 2|3  0.82476    0.08193  10.066   <1e-06 ***
```
