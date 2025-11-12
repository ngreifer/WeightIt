# Methods for `glm_weightit()` objects

This page documents methods for objects returned by
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`lm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`multinom_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
and
[`coxph_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md).
[`predict()`](https://rdrr.io/r/stats/predict.html) methods are
described at
[`predict.glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/predict.glm_weightit.md)
and [`anova()`](https://rdrr.io/r/stats/anova.html) methods are
described at
[`anova.glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/anova.glm_weightit.md).

## Usage

``` r
# S3 method for class 'glm_weightit'
summary(object, ci = FALSE, level = 0.95, transform = NULL, vcov = NULL, ...)

# S3 method for class 'multinom_weightit'
summary(object, ci = FALSE, level = 0.95, transform = NULL, vcov = NULL, ...)

# S3 method for class 'ordinal_weightit'
summary(
  object,
  ci = FALSE,
  level = 0.95,
  transform = NULL,
  thresholds = TRUE,
  vcov = NULL,
  ...
)

# S3 method for class 'coxph_weightit'
summary(object, ci = FALSE, level = 0.95, transform = NULL, vcov = NULL, ...)

# S3 method for class 'glm_weightit'
print(x, digits = max(3L, getOption("digits") - 3L), ...)

# S3 method for class 'glm_weightit'
vcov(object, complete = TRUE, vcov = NULL, ...)

# S3 method for class 'glm_weightit'
estfun(x, asympt = TRUE, ...)

# S3 method for class 'glm_weightit'
update(object, formula. = NULL, ..., evaluate = TRUE)
```

## Arguments

- object, x:

  an output from one of the above modeling functions.

- ci:

  `logical`; whether to display Wald confidence intervals for estimated
  coefficients. Default is `FALSE`. (Note: this argument can also be
  supplied as `conf.int`.)

- level:

  when `ci = TRUE`, the desired confidence level.

- transform:

  the function used to transform the coefficients, e.g., `exp` (which
  can also be supplied as a string, e.g., `"exp"`); passed to
  [`match.fun()`](https://rdrr.io/r/base/match.fun.html) before being
  used on the coefficients. When `ci = TRUE`, this is also applied to
  the confidence interval bounds. If specified, the standard error will
  be omitted from the output. Default is no transformation.

- vcov:

  either a string indicating the method used to compute the variance of
  the estimated parameters for `object`, a function used to extract the
  variance, or the variance matrix itself. Default is to use the
  variance matrix already present in `object`. If a string or function,
  arguments passed to `...` are supplied to the method or function.
  (Note: for [`vcov()`](https://rdrr.io/r/stats/vcov.html), can also be
  supplied as `type`.)

- ...:

  for [`vcov()`](https://rdrr.io/r/stats/vcov.html) or
  [`summary()`](https://rdrr.io/r/base/summary.html) or
  [`confint()`](https://rdrr.io/r/stats/confint.html) with `vcov`
  supplied, other arguments used to compute the variance matrix
  depending on the method supplied to `vcov`, e.g., `cluster`, `R`, or
  `fwb.args`. For [`update()`](https://rdrr.io/r/stats/update.html),
  additional arguments to the call or arguments with changed values. See
  [`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
  for details.

- thresholds:

  `logical`; whether to include thresholds in the
  [`summary()`](https://rdrr.io/r/base/summary.html) output for
  `ordinal_weightit` objects. Default is `TRUE`.

- digits:

  the number of *significant* digits to be passed to
  [`format`](https://rdrr.io/r/base/format.html)`(`[`coef`](https://rdrr.io/r/stats/coef.html)`(x), .)`
  when [`print()`](https://rdrr.io/r/base/print.html)ing.

- complete:

  `logical`; whether the full variance-covariance matrix should be
  returned also in case of an over-determined system where some
  coefficients are undefined and `coef(.)` contains `NA`s
  correspondingly. When `complete = TRUE`,
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) is compatible with
  [`coef()`](https://rdrr.io/r/stats/coef.html) also in this singular
  case.

- asympt:

  `logical`; for `estfun()`, whether to use the asymptotic empirical
  estimating functions that account for estimation of the weights (when
  `Mparts` is available). Default is `TRUE`. Set to `FALSE` to ignore
  estimation of the weights. Ignored when `Mparts` is not available or
  no argument was supplied to `weightit` in the fitting function.

- formula.:

  changes to the model formula, passed to the `new` argument of
  [`update.formula()`](https://rdrr.io/r/stats/update.formula.html).

- evaluate:

  `logical`; whether to evaluate the call (`TRUE`, the default) or just
  return it.

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns a
`summary.glm_weightit()` object, which has its own
[`print()`](https://rdrr.io/r/base/print.html) method. For
[`coxph_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
objects, the [`print()`](https://rdrr.io/r/base/print.html) and
[`summary()`](https://rdrr.io/r/base/summary.html) methods are more like
those for `glm` objects than for `coxph` objects.

Otherwise, all methods return the same type of object as their generics.

## Details

[`vcov()`](https://rdrr.io/r/stats/vcov.html) by default extracts the
parameter covariance matrix already computed by the fitting function,
and [`summary()`](https://rdrr.io/r/base/summary.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) uses this covariance
matrix to compute standard errors and Wald confidence intervals
(internally calling
[`confint.lm()`](https://rdrr.io/r/stats/confint.html)), respectively.
Supplying arguments to `vcov` or `...` will compute a new covariance
matrix. If `cluster` was supplied to the original fitting function, it
will be incorporated into any newly computed covariance matrix unless
`cluster = NULL` is specified in
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`summary()`](https://rdrr.io/r/base/summary.html), or
[`confint()`](https://rdrr.io/r/stats/confint.html). For other arguments
(e.g., `R` and `fwb.args`), the defaults are those used by
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md).
Note that for `vcov = "BS"` and `vcov = "FWB"` (and `vcov = "const"` for
`multinom_weightit` or `ordinal_weightit` objects), the environment for
the fitting function is used, so any changes to that environment may
affect calculation. It is always safer to simply recompute the fitted
object with a new covariance matrix than to modify it with the `vcov`
argument, but it can be quicker to just request a new covariance matrix
when refitting the model is slow.

[`update()`](https://rdrr.io/r/stats/update.html) updates a fitted model
object with new arguments, e.g., a new model formula, dataset, or
variance matrix. When only arguments that control the computation of the
variance are supplied, only the variance will be recalculated (i.e., the
parameters will not be re-estimated). When `data` is supplied,
`weightit` is not supplied, and a `weightit` object was originally
passed to the model fitting function, the `weightit` object will be
re-fit with the new dataset before the model is refit using the new
weights and new data. That is, calling `update(obj, data = d)` is
equivalent to calling
`update(obj, data = d, weightit = update(obj$weightit, data = d))` when
a `weightit` object was supplied to the model fitting function.
Similarly, supplying `s.weights` or `weights` passes the argument
through to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
to be refit. When `s.weights` or `weights` are supplied and no
`weightit` object is present, a fake one containing just the supplied
weights will be created.

`estfun()` extracts the empirical estimating functions for the fitted
model, optionally accounting for the estimation of the weights (if
available). This, along with `bread()`, is used by
[`sandwich::sandwich()`](https://sandwich.R-Forge.R-project.org/reference/sandwich.html)
to compute the robust covariance matrix of the estimated coefficients.
See
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
and [`vcov()`](https://rdrr.io/r/stats/vcov.html) above for more
details.

## See also

[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
for the page documenting
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`lm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`multinom_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
and
[`coxph_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md).
[`summary.glm()`](https://rdrr.io/r/stats/summary.glm.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`confint()`](https://rdrr.io/r/stats/confint.html) for the relevant
methods pages.
[`predict.glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/predict.glm_weightit.md)
for computing predictions from the models.
[`anova.glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/anova.glm_weightit.md)
for comparing models using a Wald test.

[`sandwich::estfun()`](https://sandwich.R-Forge.R-project.org/reference/estfun.html)
and
[`sandwich::bread()`](https://sandwich.R-Forge.R-project.org/reference/bread.html)
for the `estfun()` and `bread()` generics.

## Examples

``` r
## See examples at ?glm_weightit
```
