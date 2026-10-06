# Fitting (Weighted) Cox Proportional Hazards Models

`coxph_weightit()` fits a Cox proportional hazards model with a
covariance matrix that accounts for estimation of weights, if supplied,
and is a wrapper for functions in the survival package. By default, this
function uses M-estimation to construct a robust covariance matrix using
the estimating equations for the weighting model and the outcome model
when available.

## Usage

``` r
coxph_weightit(
  formula,
  data,
  weightit = NULL,
  vcov = NULL,
  cluster,
  R = 500L,
  control = list(...),
  x = FALSE,
  y = TRUE,
  fwb.args = list(),
  br = FALSE,
  weights,
  ...
)
```

## Arguments

- formula:

  an object of class [`formula`](https://rdrr.io/r/stats/formula.html)
  (or one that can be coerced to that class): a symbolic description of
  the model to be fitted. Should include a
  [`Surv()`](https://rdrr.io/pkg/survival/man/Surv.html) term as the
  response. See
  [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) for
  how this should be specified.

- data:

  a data frame containing the variables in the model. If not found in
  data, the variables are taken from `environment(formula)`, typically
  the environment from which the function is called.

- weightit:

  a `weightit` or `weightitMSM` object; the output of a call to
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  or
  [`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
  If neither `weightit` nor `weights` is supplied, an unweighted model
  will be fit.

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
  [`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html)
  by default. Cluster-robust variance calculations are available only
  when `vcov` is `"asympt"`, `"HC0"`, `"BS"`, or `"FWB"`.

- R:

  the number of bootstrap replications when `vcov` is `"BS"` or `"FWB"`.
  Default is 500. Ignored otherwise.

- control:

  a list of parameters for controlling the fitting process, passed to
  [`survival::coxph.control()`](https://rdrr.io/pkg/survival/man/coxph.control.html)
  .

- x, y:

  logical values indicating whether the response vector and model matrix
  used in the fitting process should be returned as components of the
  returned value.

- fwb.args:

  an optional list of further arguments to supply to
  [`fwb::fwb()`](https://ngreifer.github.io/fwb/reference/fwb.html) when
  `vcov = "FWB"`.

- br:

  `logical`; whether to use bias reduction, i.e., to maximize the
  partial likelihood penalized by the Jeffreys invariant prior as
  described by Heinze and Schemper (2001) rather than the partial
  likelihood. This yields estimates with smaller asymptotic bias that
  are always finite, even when the maximum partial likelihood estimates
  are not (e.g., under a monotone likelihood). Default is `FALSE`. See
  Details.

- weights:

  an optional vector of weights to be used in the fitting process, which
  are treated as fixed. Can be supplied as for
  [`glm()`](https://rdrr.io/r/stats/glm.html) (e.g., as the unquoted
  name of a variable in `data`), as a numeric vector, or as a string
  containing the name of a variable in `data`. Only one of `weights` and
  `weightit` can be supplied; a `weightit` or `weightitMSM` object
  supplied to `weights` is treated as though it had been supplied to
  `weightit`. When `weights` is supplied, the default `vcov` is `"HC0"`,
  and bootstrapping holds the weights fixed rather than re-estimating
  them.

- ...:

  other arguments passed to
  [`survival::coxph.control()`](https://rdrr.io/pkg/survival/man/coxph.control.html)
  .

## Value

A `coxph_weightit` object, which inherits from `coxph`. See
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
call. When `weights` is supplied instead, only `(weights)` is included,
containing those weights.

## Details

`coxph_weightit()` is essentially a simplified version of
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) to
fit weighted survival models that optionally computes a coefficient
variance matrix that can be adjusted to account for estimation of the
weights if a `weightit` or `weightitMSM` object is supplied to the
`weightit` argument. It differs from `coxph()` in a few ways:

- the `cluster` argument (if used) should be specified as a one-sided
  formula (which can include multiple clustering variables) and uses a
  small sample correction for cluster variance estimates when specified

- Special formula components, such as `strata()`, `cluster()`,
  `pspline()`, `frailty()`, `ridge()`, and `tt()` are not allowed

- Only right censoring is allowed, and only two-state models are allowed
  (i.e., the `Surv()` component of `formula` must be of the form
  `Surv(time, event)`)

- Time-varying predictors are not allowed and there must be one
  observation per unit (and the `id` and `istate` arguments to `coxph()`
  are ignored)

- Weights of 0 are allowed (`coxph()` rejects them). Such units are
  omitted from the model fit, which is exact rather than an
  approximation: a unit with a weight of 0 contributes neither its own
  term nor anything to any risk set denominator, both of which are
  weighted by the same weights. Missing values in the model variables
  are tolerated for these units, which makes it possible to fit an
  outcome model after censoring, where the event time is unascertained
  for censored units (see
  [`.cens()`](https://ngreifer.github.io/WeightIt/reference/dot-cens.md)).
  Missing values in units with a nonzero weight produce an error.

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

### Bias reduction

When `br = TRUE`, the coefficients maximize the partial likelihood
penalized by the Jeffreys invariant prior, which is equivalent to
solving the bias-reducing adjusted score equations of Firth (1993) and
removes the first-order term in the asymptotic bias of the estimates.
The estimates are finite even when the likelihood is monotone, i.e.,
when some maximum partial likelihood estimates are infinite, which
happens most often in small samples with heavy censoring and strongly
predictive covariates (Heinze & Schemper, 2001). Without weights,
estimation should align with that from
[`coxphf::coxphf()`](https://rdrr.io/pkg/coxphf/man/coxphf.html) . The
penalized likelihood is maximized by Newton-Raphson iterations using the
information matrix, which converge more slowly than those for the
unpenalized likelihood, so the maximum number of iterations is 100
rather than 20 unless `iter.max` is supplied.

As in
[`multinom_weightit()`](https://ngreifer.github.io/WeightIt/reference/multinom_weightit.md)
and
[`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/ordinal_weightit.md),
the weights are scaled to have a mean of 1 among the units with a
nonzero weight before fitting (Mukhopadhyay, 2020), which makes the
estimates invariant to multiplying the weights by a constant, and the
reported variance matrix uses the information matrix at the estimates,
i.e., the adjustment is treated as fixed. M-estimation and bootstrapping
can be used with `br = TRUE` just as they can without it. Note that
under a monotone likelihood, the distribution of the estimates can be
far from normal, so Wald confidence intervals can perform poorly (Heinze
& Schemper, 2001); bootstrapping may be preferable in that case.

## References

Firth, D. (1993). Bias reduction of maximum likelihood estimates.
*Biometrika*, 80(1), 27–38.
[doi:10.1093/biomet/80.1.27](https://doi.org/10.1093/biomet/80.1.27)

Heinze, G., & Schemper, M. (2001). A solution to the problem of monotone
likelihood in Cox regression. *Biometrics*, 57(1), 114–119.
[doi:10.1111/j.0006-341X.2001.00114.x](https://doi.org/10.1111/j.0006-341X.2001.00114.x)

Mukhopadhyay, P. K. (2020). Firth's penalized likelihood for
proportional hazards regressions for complex surveys. *Survey
Methodology*, 46(2), 215–241.

## See also

- [`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html) for
  fitting Cox proportional hazards models without adjusting standard
  errors for estimation of the weights.

- [`coxphf::coxphf()`](https://rdrr.io/pkg/coxphf/man/coxphf.html) for
  fitting bias-reduced Cox proportional hazards models that do not
  account for estimation of the weights.

- [`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
  for fitting generalized linear models that adjust for estimation of
  the weights.

- [`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/ordinal_weightit.md)
  and
  [`multinom_weightit()`](https://ngreifer.github.io/WeightIt/reference/multinom_weightit.md)
  for fitting ordinal and multinomial regression models that adjust for
  estimation of the weights.

## Examples

``` r
# See `vignette("estimating-effects")` for an example
```
