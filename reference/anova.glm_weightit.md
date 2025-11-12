# Methods for `glm_weightit()` objects

[`anova()`](https://rdrr.io/r/stats/anova.html) is used to compare
nested models fit with
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
`mutinom_weightit()`,
[`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
or
[`coxph_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
using a Wald test that incorporates uncertainty in estimating the
weights (if any).

## Usage

``` r
# S3 method for class 'glm_weightit'
anova(
  object,
  object2,
  test = "Chisq",
  method = "Wald",
  tolerance = 1e-07,
  vcov = NULL,
  ...
)
```

## Arguments

- object, object2:

  an output from one of the above modeling functions. `object2` is
  required.

- test:

  the type of test statistic used to compare models. Currently only
  `"Chisq"` (the chi-square statistic) is allowed.

- method:

  the kind of test used to compare models. Currently only `"Wald"` is
  allowed.

- tolerance:

  for the Wald test, the tolerance used to determine if models are
  symbolically nested.

- vcov:

  either a string indicating the method used to compute the variance of
  the estimated parameters for `object`, a function used to extract the
  variance, or the variance matrix itself. Default is to use the
  variance matrix already present in `object`. If a string or function,
  arguments passed to `...` are supplied to the method or function.
  (Note: for [`vcov()`](https://rdrr.io/r/stats/vcov.html), can also be
  supplied as `type`.)

- ...:

  other arguments passed to the function used for computing the
  parameter variance matrix, if supplied as a string or function, e.g.,
  `cluster`, `R`, or `fwb.args`.

## Value

An object of class `"anova"` inheriting from class `"data.frame"`.

## Details

[`anova()`](https://rdrr.io/r/stats/anova.html) performs a Wald test to
compare two fitted models. The models must be nested, but they don't
have to be nested symbolically (i.e., the names of the coefficients of
the smaller model do not have to be a subset of the names of the
coefficients of the larger model). The larger model must be supplied to
`object` and the smaller to `object2`. Both models must contain the same
units, weights (if any), and outcomes. The variance-covariance matrix of
the coefficients of the smaller model is not used.

## See also

[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
for the page documenting
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`lm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
[`multinom_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md),
and
[`coxph_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md).
[`anova.glm()`](https://rdrr.io/r/stats/anova.glm.html) for model
comparison of `glm` objects.

## Examples

``` r
data("lalonde", package = "cobalt")

# Model comparison for any relationship between `treat`
# and `re78` (not the same as testing for the ATE)
fit1 <- glm_weightit(
  re78 ~ treat * (age + educ + race + married + nodegree +
                    re74 + re75), data = lalonde
)

fit2 <- glm_weightit(
  re78 ~ age + educ + race + married + nodegree +
    re74 + re75, data = lalonde
)

anova(fit1, fit2)
#> 
#> Wald test
#> Variance: HC0 robust
#> 
#> Model 1: re78 ~ treat * (age + educ + race + married + nodegree + re74 + re75)
#> Model 2: re78 ~ age + educ + race + married + nodegree + re74 + re75
#> 
#>   Res.Df Df  Chisq Pr(>Chisq)  
#> 1    596                       
#> 2    605  9 17.563    0.04059 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Using the usual maximum likelihood variance matrix
anova(fit1, fit2, vcov = "const")
#> 
#> Wald test
#> Variance: maximum likelihood
#> 
#> Model 1: re78 ~ treat * (age + educ + race + married + nodegree + re74 + re75)
#> Model 2: re78 ~ age + educ + race + married + nodegree + re74 + re75
#> 
#>   Res.Df Df  Chisq Pr(>Chisq)  
#> 1    596                       
#> 2    605  9 19.761    0.01944 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Using a bootstrapped variance matrix
anova(fit1, fit2, vcov = "BS", R = 100)
#> 
#> Wald test
#> Variance: traditional bootstrap
#> 
#> Model 1: re78 ~ treat * (age + educ + race + married + nodegree + re74 + re75)
#> Model 2: re78 ~ age + educ + race + married + nodegree + re74 + re75
#> 
#>   Res.Df Df  Chisq Pr(>Chisq)  
#> 1    596                       
#> 2    605  9 17.524    0.04111 *
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Model comparison between spline model and linear
# model; note they are nested but not symbolically
# nested
fit_s <- glm_weightit(re78 ~ splines::ns(age, df =4),
                      data = lalonde)

fit_l <- glm_weightit(re78 ~ age,
                      data = lalonde)

anova(fit_s, fit_l)
#> 
#> Wald test
#> Variance: HC0 robust
#> 
#> Model 1: re78 ~ splines::ns(age, df = 4)
#> Model 2: re78 ~ age
#> 
#>   Res.Df Df  Chisq Pr(>Chisq)   
#> 1    609                        
#> 2    612  3 14.166   0.002688 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
