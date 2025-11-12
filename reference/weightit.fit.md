# Generate Balancing Weights with Minimal Input Processing

`weightit.fit()` dispatches one of the weight estimation methods
determined by `method`. It is an internal function called by
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
and should probably not be used except in special cases. Unlike
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
`weightit.fit()` does not accept a formula and data frame interface and
instead requires the covariates and treatment to be supplied as a
numeric matrix and atomic vector, respectively. In this way,
`weightit.fit()` is to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
what [`lm.fit()`](https://rdrr.io/r/stats/lmfit.html) is to
[`lm()`](https://rdrr.io/r/stats/lm.html): a thinner, slightly faster
interface that performs minimal argument checking.

## Usage

``` r
weightit.fit(
  covs,
  treat,
  method = "glm",
  s.weights = NULL,
  by.factor = NULL,
  estimand = "ATE",
  focal = NULL,
  stabilize = FALSE,
  ps = NULL,
  missing = NULL,
  verbose = FALSE,
  include.obj = FALSE,
  ...
)
```

## Arguments

- covs:

  a numeric matrix of covariates.

- treat:

  a vector of treatment statuses.

- method:

  a string containing the name of the method that will be used to
  estimate weights. See
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  for allowable options. The default is `"glm"` for propensity score
  weighting using a generalized linear model to estimate the propensity
  score.

- s.weights:

  a numeric vector of sampling weights. See the individual pages for
  each method for information on whether sampling weights can be
  supplied.

- by.factor:

  a factor variable for which weighting is to be done within levels.
  Corresponds to the `by` argument in
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md).

- estimand:

  the desired estimand. For binary and multi-category treatments, can be
  `"ATE"`, `"ATT"`, `"ATC"`, and, for some methods, `"ATO"`, `"ATM"`, or
  `"ATOS"`. The default for both is `"ATE"`. This argument is ignored
  for continuous treatments. See the individual pages for each method
  for more information on which estimands are allowed with each method
  and what literature to read to interpret these estimands.

- focal:

  when `estimand` is set to `"ATT"` or `"ATC"`, which group to consider
  the "treated" or "control" group. This group will not be weighted, and
  the other groups will be weighted to resemble the focal group. If
  specified, `estimand` will automatically be set to `"ATT"` (with a
  warning if `estimand` is not `"ATT"` or `"ATC"`). See section
  *`estimand` and `focal`* in Details at
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md).

- stabilize:

  `logical`; whether or not to stabilize the weights. For the methods
  that involve estimating propensity scores, this involves multiplying
  each unit's weight by the proportion of units in their treatment
  group. Default is `FALSE`. Note this differs from its use with
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md).

- ps:

  a vector of propensity scores. If specified, `method` will be ignored
  and set to `"glm"`.

- missing:

  `character`; how missing data should be handled. The options depend on
  the `method` used. If `NULL`, `covs` will be checked for `NA` values,
  and if present, `missing` will be set to `"ind"`. If `""`, `covs` will
  not be checked for `NA` values; this can be faster when it is known
  there are none.

- verbose:

  `logical`; whether to print additional information output by the
  fitting function.

- include.obj:

  `logical`; whether to include in the output any fit objects created in
  the process of estimating the weights. For example, with
  `method = "glm"`, the `glm` objects containing the propensity score
  model will be included. See the individual pages for each method for
  information on what object will be included if `TRUE`.

- ...:

  other arguments for functions called by `weightit.fit()` that control
  aspects of fitting that are not covered by the above arguments.

## Value

A `weightit.fit` object with the following elements:

- weights:

  The estimated weights, one for each unit.

- treat:

  The values of the treatment variable.

- estimand:

  The estimand requested.

- method:

  The weight estimation method specified.

- ps:

  The estimated or provided propensity scores. Estimated propensity
  scores are returned for binary treatments and only when `method` is
  `"glm"`, `"gbm"`, `"cbps"`, `"ipt"`, `"super"`, or `"bart"`. The
  propensity score corresponds to the predicted probability of being
  treated; see section *`estimand` and `focal`* in Details at
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
  for how the treated group is determined.

- s.weights:

  The provided sampling weights.

- focal:

  The focal treatment level if the ATT or ATC was requested.

- fit.obj:

  When `include.obj = TRUE`, the fit object.

- info:

  Additional information about the fitting. See the individual methods
  pages for what is included.

The `weightit.fit` object does not have specialized
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), or
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods. It is
simply a list containing the above components. Use
[`as.weightit()`](https://ngreifer.github.io/WeightIt/reference/as.weightit.md)
to convert it to a `weightit` object, which does have these methods. See
Examples.

## Details

`weightit.fit()` is called by
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
after the arguments to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
have been checked and processed. `weightit.fit()` dispatches the
function used to actually estimate the weights, passing on the supplied
arguments directly. `weightit.fit()` is not meant to be used by anyone
other than experienced users who have a specific use case in mind. The
returned object contains limited information about the supplied
arguments or details of the estimation method; all that is processed by
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md).

Less argument checking or processing occurs in `weightit.fit()` than
does in
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
which means supplying incorrect arguments can result in errors, crashes,
and invalid weights, and error and warning messages may not be helpful
in diagnosing the problem. `weightit.fit()` does check to make sure
weights were actually estimated, though.

`weightit.fit()` may be most useful in speeding up simulation simulation
studies that use
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
because the covariates can be supplied as a numeric matrix, which is
often how they are generated in simulations, without having to go
through the potentially slow process of extracting the covariates and
treatment from a formula and data frame. If the user is certain the
arguments are valid (e.g., by ensuring the estimated weights are
consistent with those estimated from
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
with the same arguments), less time needs to be spent on processing the
arguments. Also, the returned object is much smaller than a `weightit`
object because the covariates are not returned alongside the weights.

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
which you should use for estimating weights unless you know better.

[`as.weightit()`](https://ngreifer.github.io/WeightIt/reference/as.weightit.md)
for converting a `weightit.fit` object to a `weightit` object.

## Examples

``` r
library("cobalt")
data("lalonde", package = "cobalt")

# Balancing covariates between treatment groups (binary)
covs <- lalonde[c("age", "educ", "race", "married",
                  "nodegree", "re74", "re75")]
## Create covs matrix, splitting any factors using
## cobalt::splitfactor()
covs_mat <- as.matrix(splitfactor(covs))

WF1 <- weightit.fit(covs_mat, treat = lalonde$treat,
                    method = "glm", estimand = "ATT")
str(WF1)
#> List of 10
#>  $ weights  : num [1:614] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ treat    : 'treat' int [1:614] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "treat.type")= chr "binary"
#>   ..- attr(*, "treated")= int 1
#>  $ estimand : chr "ATT"
#>  $ method   : chr "glm"
#>  $ ps       : num [1:614] 0.639 0.225 0.678 0.776 0.702 ...
#>  $ s.weights: num [1:614] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ focal    : int 1
#>  $ missing  : chr ""
#>  $ fit.obj  : NULL
#>  $ info     : Named list()
#>  - attr(*, "Mparts")=List of 6
#>   ..$ psi_treat :function (Btreat, Xtreat, A, SW)  
#>   ..$ wfun      :function (Btreat, Xtreat, A)  
#>   ..$ dw_dBtreat:function (Btreat, Xtreat, A, SW)  
#>   ..$ Xtreat    : num [1:614, 1:9] 1 1 1 1 1 1 1 1 1 1 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:614] "1" "2" "3" "4" ...
#>   .. .. ..$ : chr [1:9] "(Intercept)" "age" "educ" "race_hispan" ...
#>   ..$ A         : int [1:614] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..$ btreat    : Named num [1:9] 0.214 0.156 0.424 -2.082 -3.065 ...
#>   .. ..- attr(*, "names")= chr [1:9] "(Intercept)" "age" "educ" "race_hispan" ...
#>  - attr(*, "class")= chr "weightit.fit"

# Converting to a weightit object for use with
# summary() and bal.tab()
W1 <- as.weightit(WF1, covs = covs)
W1
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, race, married, nodegree, re74, re75
summary(W1)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                 Max
#> treated 1.            ||                    1.   
#> control 0.009 |---------------------------| 3.743
#> 
#> - Units with the 5 most extreme weights by group:
#>                                     
#>             1     2    3     4     5
#>  treated    1     1    1     1     1
#>           412   388  226   196   118
#>  control 3.03 3.059 3.24 3.523 3.743
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.000 0.000   0.000       0
#> control       1.818 1.289   1.098       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted     99.82     185
bal.tab(W1)
#> Balance Measures
#>                 Type Diff.Adj
#> prop.score  Distance  -0.0205
#> age          Contin.   0.1188
#> educ         Contin.  -0.0284
#> race_black    Binary  -0.0022
#> race_hispan   Binary   0.0002
#> race_white    Binary   0.0021
#> married       Binary   0.0186
#> nodegree      Binary   0.0184
#> re74         Contin.  -0.0021
#> re75         Contin.   0.0110
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted     99.82     185
```
