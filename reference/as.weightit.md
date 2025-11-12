# Create a `weightit` object manually

This function allows users to get the benefits of a `weightit` object
when using weights not estimated with
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
These benefits include diagnostics, plots, and direct compatibility with
cobalt for assessing balance.

## Usage

``` r
as.weightit(x, ...)

# S3 method for class 'weightit.fit'
as.weightit(x, covs = NULL, ...)

# Default S3 method
as.weightit(
  x,
  treat,
  covs = NULL,
  estimand = NULL,
  s.weights = NULL,
  ps = NULL,
  ...
)

as.weightitMSM(x, ...)

# Default S3 method
as.weightitMSM(
  x,
  treat.list,
  covs.list = NULL,
  estimand = NULL,
  s.weights = NULL,
  ps.list = NULL,
  ...
)
```

## Arguments

- x:

  required; a `numeric` vector of weights, one for each unit, or a
  `weightit.fit` object from
  [`weightit.fit()`](https://ngreifer.github.io/WeightIt/reference/weightit.fit.md).

- ...:

  additional arguments. These must be named. They will be included in
  the output object.

- covs:

  an optional `data.frame` of covariates. For using WeightIt functions,
  this is not necessary, but for use with cobalt it is. Note that when
  using with a `weightit.fit` object, this should not be the matrix
  supplied to the `covs` argument of
  [`weightit.fit()`](https://ngreifer.github.io/WeightIt/reference/weightit.fit.md)
  unless there are no factor/character variables in it. Ideally this is
  the original, unprocessed covariate data frame with factor variables
  included.

- treat:

  a vector of treatment statuses, one for each unit. Required when `x`
  is a vector of weights.

- estimand:

  an optional `character` of length 1 giving the estimand. The text is
  not checked.

- s.weights:

  an optional `numeric` vector of sampling weights, one for each unit.

- ps:

  an optional `numeric` vector of propensity scores, one for each unit.

- treat.list:

  a list of treatment statuses at each time point.

- covs.list:

  an optional list of `data.frame`s of covariates of covariates at each
  time point. For using WeightIt functions, this is not necessary, but
  for use with cobalt it is.

- ps.list:

  an optional list of `numeric` vectors of propensity scores at each
  time point.

## Value

An object of class `weightit` (for `as.weightit()`) or `weightitMSM`
(for `as.weightitMSM()`).

## Examples

``` r
treat <- rbinom(500, 1, .3)
weights <- rchisq(500, df = 2)

W <- as.weightit(weights, treat = treat, estimand = "ATE")
summary(W)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                  Max
#> treated 0.004 |-------------------------|   11.167
#> control 0.006 |---------------------------| 12.028
#> 
#> - Units with the 5 most extreme weights by group:
#>                                        
#>            147   126    58    21     19
#>  treated 6.888 7.372 7.614  8.83 11.167
#>            294   280   109    59     44
#>  control 8.195 8.826 8.843 9.904 12.028
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.980 0.736   0.410       0
#> control       0.977 0.739   0.422       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  343.    157.  
#> Weighted    175.77   80.33

# See ?weightit.fit for using as.weightit() with a
# weightit.fit object.
```
