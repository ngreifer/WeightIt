# Predictions for `glm_weightit` objects

[`predict()`](https://rdrr.io/r/stats/predict.html) generates
predictions for models fit using
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`multinom_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
or
[`coxph_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md).
This page only details the
[`predict()`](https://rdrr.io/r/stats/predict.html) methods after using
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
or
[`multinom_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md).
See
[`survival::predict.coxph()`](https://rdrr.io/pkg/survival/man/predict.coxph.html)
for predictions when fitting Cox proportional hazards models using
[`coxph_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md).

## Usage

``` r
# S3 method for class 'glm_weightit'
predict(object, newdata = NULL, type = "response", na.action = na.pass, ...)

# S3 method for class 'ordinal_weightit'
predict(
  object,
  newdata = NULL,
  type = "response",
  na.action = na.pass,
  values = NULL,
  level = NULL,
  ...
)

# S3 method for class 'multinom_weightit'
predict(
  object,
  newdata = NULL,
  type = "response",
  na.action = na.pass,
  values = NULL,
  level = NULL,
  ...
)
```

## Arguments

- object:

  a `glm_weightit` object.

- newdata:

  optionally, a data frame in which to look for variables with which to
  predict. If omitted, the fitted values applied to the original dataset
  are used.

- type:

  the type of prediction desired. Allowable options include
  `"response"`, for predictions on the scale of the original response
  variable (also `"probs"`); `"link"`, for predictions on the scale of
  the linear predictor (also `"lp"`); `"class"`, for the modal predicted
  category for ordinal and multinomial models; `"mean"`, for the
  expected value of the outcome for ordinal and multinomial models; and
  `"stdlv"`, for the standardized latent variable values for ordinal
  models. See Details for more information. The default is `"response"`
  for all models, which differs from
  [`stats::predict.glm()`](https://rdrr.io/r/stats/predict.glm.html).

- na.action:

  function determining what should be done with missing values in
  `newdata`. The default is to predict `NA`.

- ...:

  further arguments passed to or from other methods.

- values:

  when `type = "mean"`, the numeric values each level corresponds to.
  Should be supplied as a named vector with outcome levels as the names.
  If `NULL` and the outcome levels can be converted to numeric, those
  will be used. See Details.

- level:

  when `type = "response"` for ordinal and multinomial models, an
  optional string or number corresponding to the outcome level for which
  the predictions are to be produced. If `NULL` (the default), a matrix
  of predictions for all levels will be produced.

## Value

A numeric vector containing the desired predictions, except for the
following circumstances when an ordinal or multinomial model was fit:

- when `type = "response"` and `levels = NULL`, a numeric matrix with a
  row for each unit and a column for each level of the outcome with the
  predicted probability of the corresponding outcome in the cells

- when `type = "class"`, a factor with the modal predicted class for
  each unit; for ordinal models, this will be an ordered factor.

## Details

For generalized linear models other than ordinal and multinomial models,
see [`stats::predict.glm()`](https://rdrr.io/r/stats/predict.glm.html)
for more information on how predictions are computed and which arguments
can be specified. Note that standard errors cannot be computed for the
predictions using `predict.glm_weightit()`.

For ordinal and multinomial models, setting `type = "mean"` computes the
expected value of the outcome for each unit; this corresponds to the sum
of the values supplied in `values` weighted by the predicted probability
of those values. If `values` is omitted,
[`predict()`](https://rdrr.io/r/stats/predict.html) will attempt to
convert the outcome levels to numeric values, and if this cannot be
done, an error will be thrown. `values` should be specified as a named
vector, e.g., `values = c(one = 1, two = 2, three = 3)`, where `"one"`,
`"two"`, and `"three"` are the original outcome levels, and 1, 2, and 3
are the numeric values they correspond to. This method only makes sense
to use if the outcome levels meaningfully correspond to numeric values.

For ordinal models, setting `type = "link"` (also `"lp"`) computes the
linear predictor without including the thresholds. This can be
interpreted as the prediction of the latent variable underlying the
ordinal response. This cannot be used with multinomial models. Setting
`type = "stdlv"` standardizes these predictions by the implied standard
deviation of the ordinal responses, which is a function of the link
function, the original covariates, and the coefficient estimates.

## See also

[`stats::predict.glm()`](https://rdrr.io/r/stats/predict.glm.html) for
predictions from generalized linear models.
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
for the fitting function.
[`survival::predict.coxph()`](https://rdrr.io/pkg/survival/man/predict.coxph.html)
for predictions from Cox proportional hazards models.

## Examples

``` r
data("lalonde", package = "cobalt")

# Logistic regression model
fit1 <- glm_weightit(
  re78 > 0 ~ treat * (age + educ + race + married +
                        re74 + re75),
  data = lalonde, family = binomial, vcov = "none")

summary(predict(fit1))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.2747  0.6986  0.7868  0.7671  0.8429  1.0000 

# G-computation using predicted probabilities
p0 <- predict(fit1, type = "response",
              newdata = transform(lalonde,
                                  treat = 0))

p1 <- predict(fit1, type = "response",
              newdata = transform(lalonde,
                                  treat = 1))

mean(p1) - mean(p0)
#> [1] 0.0921469

# Multinomial logistic regression model
lalonde$re78_3 <- factor(findInterval(lalonde$re78,
                                      c(0, 5e3, 1e4)),
                         labels = c("low", "med", "high"))

fit2 <- multinom_weightit(
  re78_3 ~ treat * (age + educ + race + married +
                      re74 + re75),
  data = lalonde, vcov = "none")

# Predicted probabilities
head(predict(fit2))
#>         low       med      high
#> 1 0.4897563 0.1814665 0.3287772
#> 2 0.5027746 0.3133358 0.1838896
#> 3 0.5453691 0.1974158 0.2572151
#> 4 0.5809239 0.2219530 0.1971231
#> 5 0.5964825 0.2995950 0.1039225
#> 6 0.6217877 0.2666334 0.1115788

# Predicted probabilities for a single level
head(predict(fit2, level = "low"))
#>         1         2         3         4         5         6 
#> 0.4897563 0.5027746 0.5453691 0.5809239 0.5964825 0.6217877 

# Class assignment accuracy
mean(predict(fit2, type = "class") == lalonde$re78_3)
#> [1] 0.5635179

# G-computation using expected value of the outcome
values <- c("low" = 2500,
            "med" = 7500,
            "high" = 12500)

p0 <- predict(fit2, type = "mean", values = values,
              newdata = transform(lalonde,
                                  treat = 0))

p1 <- predict(fit2, type = "mean", values = values,
              newdata = transform(lalonde,
                                  treat = 1))

mean(p1) - mean(p0)
#> [1] 677.4256
# \donttest{
# Ordinal logistic regression
fit3 <- ordinal_weightit(
  re78 ~ treat * (age + educ + race + married +
                    re74 + re75),
  data = lalonde, vcov = "none")

# G-computation using expected value of the outcome;
# using original outcome values
p0 <- predict(fit3, type = "mean",
              newdata = transform(lalonde,
                                  treat = 0))

p1 <- predict(fit3, type = "mean",
              newdata = transform(lalonde,
                                  treat = 1))

mean(p1) - mean(p0)
#> [1] 943.9958
# }
```
