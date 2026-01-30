# Estimate Balancing Weights

`weightit()` allows for the easy generation of balancing weights using a
variety of available methods for binary, continuous, and multi-category
treatments. Many of these methods exist in other packages, which
`weightit()` calls; these packages must be installed to use the desired
method.

## Usage

``` r
weightit(
  formula,
  data = NULL,
  method = "glm",
  estimand = "ATE",
  stabilize = FALSE,
  focal = NULL,
  by = NULL,
  s.weights = NULL,
  ps = NULL,
  missing = NULL,
  verbose = FALSE,
  include.obj = FALSE,
  keep.mparts = TRUE,
  ...
)
```

## Arguments

- formula:

  a formula with a treatment variable on the left hand side and the
  covariates to be balanced on the right hand side. See
  [`glm()`](https://rdrr.io/r/stats/glm.html) for more details.
  Interactions and functions of covariates are allowed.

- data:

  an optional data set in the form of a data frame that contains the
  variables in `formula`.

- method:

  a string of length 1 containing the name of the method that will be
  used to estimate weights. See Details below for allowable options. The
  default is `"glm"` for propensity score weighting using a generalized
  linear model to estimate the propensity score.

- estimand:

  the desired estimand. For binary and multi-category treatments, can be
  `"ATE"`, `"ATT"`, `"ATC"`, and, for some methods, `"ATO"`, `"ATM"`, or
  `"ATOS"`. The default for both is `"ATE"`. This argument is ignored
  for continuous treatments. See the individual pages for each method
  for more information on which estimands are allowed with each method
  and what literature to read to interpret these estimands.

- stabilize:

  whether or not and how to stabilize the weights. If `TRUE`, each
  unit's weight will be multiplied by a standardization factor, which is
  the the unconditional probability (or density) of each unit's observed
  treatment value. If a formula, a generalized linear model will be fit
  with the included predictors, and the inverse of the corresponding
  weight will be used as the standardization factor. Can only be used
  with continuous treatments or when `estimand = "ATE"`. Default is
  `FALSE` for no standardization. See also the `num.formula` argument at
  [`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
  For continuous treatments, weights are already stabilized, so setting
  `stabilize = TRUE` will be ignored with a warning (supplying a formula
  still works).

- focal:

  when `estimand` is set to `"ATT"` or `"ATC"`, which group to consider
  the "treated" or "control" group. This group will not be weighted, and
  the other groups will be weighted to resemble the focal group. If
  specified, `estimand` will automatically be set to `"ATT"` (with a
  warning if `estimand` is not `"ATT"` or `"ATC"`). See section
  *`estimand` and `focal`* in Details below.

- by:

  a string containing the name of the variable in `data` for which
  weighting is to be done within categories or a one-sided formula with
  the stratifying variable on the right-hand side. For example, if
  `by = "gender"` or `by = ~gender`, a separate propensity score model
  or optimization will occur within each level of the variable
  `"gender"`. Only one `by` variable is allowed; to stratify by multiply
  variables simultaneously, create a new variable that is a full cross
  of those variables using
  [`interaction()`](https://rdrr.io/r/base/interaction.html).

- s.weights:

  a vector of sampling weights or the name of a variable in `data` that
  contains sampling weights. These can also be matching weights if
  weighting is to be used on matched data. See the individual pages for
  each method for information on whether sampling weights can be
  supplied.

- ps:

  a vector of propensity scores or the name of a variable in `data`
  containing propensity scores. If not `NULL`, `method` is ignored
  unless it is a user-supplied function, and the propensity scores will
  be used to create weights. `formula` must include the treatment
  variable in `data`, but the listed covariates will play no role in the
  weight estimation. Using `ps` is similar to calling
  [`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
  directly, but produces a full `weightit` object rather than just
  producing weights.

- missing:

  `character`; how missing data should be handled. The options and
  defaults depend on the `method` used. Ignored if no missing data is
  present. It should be noted that multiple imputation outperforms all
  available missingness methods available in `weightit()` and should
  probably be used instead. Consider the
  [MatchThem](https://CRAN.R-project.org/package=MatchThem) package for
  the use of `weightit()` with multiply imputed data.

- verbose:

  `logical`; whether to print additional information output by the
  fitting function.

- include.obj:

  `logical`; whether to include in the output any fit objects created in
  the process of estimating the weights. For example, with
  `method = "glm"`, the `glm` objects containing the propensity score
  model will be included. See the individual pages for each method for
  information on what object will be included if `TRUE`.

- keep.mparts:

  `logical`; whether to include in the output components necessary to
  estimate standard errors that account for estimation of the weights in
  [`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md).
  Default is `TRUE` if such parts are present. See the individual pages
  for each method for whether these components are produced. Set to
  `FALSE` to keep the output object smaller, e.g., if standard errors
  will not be computed using
  [`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md).

- ...:

  other arguments for functions called by `weightit()` that control
  aspects of fitting that are not covered by the above arguments. See
  Details.

## Value

A `weightit` object with the following elements:

- weights:

  The estimated weights, one for each unit.

- treat:

  The values of the treatment variable.

- covs:

  The covariates used in the fitting. Only includes the raw covariates,
  which may have been altered in the fitting process.

- estimand:

  The estimand requested.

- method:

  The weight estimation method specified.

- ps:

  The estimated or provided propensity scores. Estimated propensity
  scores are returned for binary treatments and only when `method` is
  `"glm"`, `"gbm"`, `"cbps"`, `"ipt"`, `"super"`, or `"bart"`. The
  propensity score corresponds to the predicted probability of being
  treated; see section *`estimand` and `focal`* in Details for how the
  treated group is determined.

- s.weights:

  The provided sampling weights.

- focal:

  The focal treatment level if the ATT or ATC was requested.

- by:

  A data.frame containing the `by` variable when specified.

- obj:

  When `include.obj = TRUE`, the fit object.

- info:

  Additional information about the fitting. See the individual methods
  pages for what is included.

When `keep.mparts` is `TRUE` (the default) and the chosen method is
compatible with M-estimation, the components related to M-estimation for
use in
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
are stored in the `"Mparts"` attribute. When `by` is specified,
`keep.mparts` is set to `FALSE`.

## Details

The primary purpose of `weightit()` is as a dispatcher to functions that
perform the estimation of balancing weights using the requested
`method`. Below are the methods allowed and links to pages containing
more information about them, including additional arguments and outputs
(e.g., when `include.obj = TRUE`), how missing values are treated, which
estimands are allowed, and whether sampling weights are allowed.

|                                                                                    |                                                                            |
|------------------------------------------------------------------------------------|----------------------------------------------------------------------------|
| [`"glm"`](https://ngreifer.github.io/WeightIt/reference/method_glm.md)             | Propensity score weighting using generalized linear models                 |
| [`"gbm"`](https://ngreifer.github.io/WeightIt/reference/method_gbm.md)             | Propensity score weighting using generalized boosted modeling              |
| [`"cbps"`](https://ngreifer.github.io/WeightIt/reference/method_cbps.md)           | Covariate Balancing Propensity Score weighting                             |
| [`"npcbps"`](https://ngreifer.github.io/WeightIt/reference/method_npcbps.md)       | Non-parametric Covariate Balancing Propensity Score weighting              |
| [`"ebal"`](https://ngreifer.github.io/WeightIt/reference/method_ebal.md)           | Entropy balancing                                                          |
| [`"ipt"`](https://ngreifer.github.io/WeightIt/reference/method_ipt.md)             | Inverse probability tilting                                                |
| [`"optweight"`](https://ngreifer.github.io/WeightIt/reference/method_optweight.md) | Stable balancing weights                                                   |
| [`"super"`](https://ngreifer.github.io/WeightIt/reference/method_super.md)         | Propensity score weighting using SuperLearner                              |
| [`"bart"`](https://ngreifer.github.io/WeightIt/reference/method_bart.md)           | Propensity score weighting using Bayesian additive regression trees (BART) |
| [`"energy"`](https://ngreifer.github.io/WeightIt/reference/method_energy.md)       | Energy balancing                                                           |
| [`"cfd"`](https://ngreifer.github.io/WeightIt/reference/method_cfd.md)             | Characteristic function distance balancing                                 |

`method` can also be supplied as a user-defined function; see
[`method_user`](https://ngreifer.github.io/WeightIt/reference/method_user.md)
for instructions and examples. Setting `method = NULL` computes unit
weights.

### `estimand` and `focal`

For binary and multi-category treatments, the argument to `estimand`
determines what distribution the weighted sample should resemble. When
set to `"ATE"`, this requests that each group resemble the full sample.
When set to `"ATO"`, `"ATM"`, or `"ATOS"` (for the methods that allow
them), this requests that each group resemble an "overlap" sample. When
set to `"ATT"` or `"ATC"`, this requests that each group resemble the
treated or control group, respectively (termed the "focal" group).
Weights are set to 1 for the focal group.

How does `weightit()` decide which group is the treated and which group
is the control? For binary treatments, several heuristics are used. The
first is by checking whether a valid argument to `focal` was supplied
containing the name of the focal group, which is the treated group when
`estimand = "ATT"` and the control group when `estimand = "ATC"`. If
`focal` is not supplied, guesses are made using the following criteria,
evaluated in order:

- If the treatment variable is `logical`, `TRUE` is considered treated
  and `FALSE` control.

- If the treatment is numeric (or a string or factor with values that
  can be coerced to numeric values), if 0 is one of the values, it is
  considered the control, and otherwise, the lower value is considered
  the control (with the other considered treated).

- If exactly one of the treatment values is `"t"`, `"tr"`, `"treat"`,
  `"treated"`, or `"exposed"`, it is considered the treated (and the
  other control).

- If exactly one of the treatment values is `"c"`, `"co"`, `"ctrl"`,
  `"control"`, or `"unexposed"`, it is considered the control (and the
  other treated).

- If the treatment variable is a factor, the first level is considered
  control and the second treated.

- The lowest value after sorting with
  [`sort()`](https://rdrr.io/r/base/sort.html) is considered control and
  the other treated.

To be safe, it is best to code your binary treatment variable as `0` for
control and `1` for treated. Otherwise, `focal` should be supplied when
requesting the ATT or ATC. For multi-category treatments, `focal` is
required when requesting the ATT or ATC; none of the heuristics above
are used.

### Citing WeightIt

When using `weightit()`, please cite both the WeightIt package (using
`citation("WeightIt")`) and the paper(s) in the references section of
the method used.

## See also

[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
for estimating weights with sequential (i.e., longitudinal) treatments
for use in estimating marginal structural models (MSMs).

[`weightit.fit()`](https://ngreifer.github.io/WeightIt/reference/weightit.fit.md),
which is a lower-level dispatcher function that accepts a matrix of
covariates and a vector of treatment statuses rather than a formula and
data frame and performs minimal argument checking and processing. It may
be useful for speeding up simulation studies for which the correct
arguments are known. In general, `weightit()` should be used.

[`summary.weightit()`](https://ngreifer.github.io/WeightIt/reference/summary.weightit.md)
for summarizing the weights

## Examples

``` r
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "glm", estimand = "ATT"))
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
#> - Weight ranges:
#> 
#>           Min                                 Max
#> treated 1.                  ||              1.   
#> control 0.022 |---------------------------| 2.044
#> 
#> - Units with the 5 most extreme weights by group:
#>                                    
#>             1     2   3     4     5
#>  treated    1     1   1     1     1
#>           410   226 224   111    84
#>  control 1.33 1.437 1.5 1.637 2.044
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.000 0.000    0.00       0
#> control       0.823 0.701    0.33       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    255.99     185
bal.tab(W1)
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.0199
#> age         Contin.   0.0459
#> educ        Contin.  -0.0360
#> married      Binary   0.0044
#> nodegree     Binary   0.0080
#> re74        Contin.  -0.0275
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    255.99     185

#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATE"))
#> A weightit object
#>  - method: "ebal" (entropy balancing)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 3-category (black, hispan, white)
#>  - estimand: ATE
#>  - covariates: age, educ, married, nodegree, re74
summary(W2)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>          Min                                 Max
#> black  0.553   |-------------------------| 5.35 
#> hispan 0.141 |----------------|            3.332
#> white  0.398  |-------|                    1.923
#> 
#> - Units with the 5 most extreme weights by group:
#>                                      
#>           203   166   163   153   152
#>   black 2.521 2.549 2.806 3.555  5.35
#>            67    43    39    36    28
#>  hispan 2.046  2.53 2.632 2.705 3.332
#>           291   285   258   205     6
#>   white 1.711 1.723 1.743 1.774 1.923
#> 
#> - Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.590 0.413   0.131       0
#> hispan       0.609 0.440   0.163       0
#> white        0.371 0.306   0.068       0
#> 
#> - Effective Sample Sizes:
#> 
#>             black hispan  white
#> Unweighted 243.    72.   299.  
#> Weighted   180.47  52.71 262.93
bal.tab(W2)
#> Balance summary across all treatment pairs
#>             Type Max.Diff.Adj
#> age      Contin.            0
#> educ     Contin.            0
#> married   Binary            0
#> nodegree  Binary            0
#> re74     Contin.            0
#> 
#> Effective sample sizes
#>             black hispan  white
#> Unadjusted 243.    72.   299.  
#> Adjusted   180.47  52.71 262.93

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cbps"))
#> A weightit object
#>  - method: "cbps" (covariate balancing propensity score weighting)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: continuous
#>  - covariates: age, educ, married, nodegree, re74
summary(W3)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>      Min                                  Max
#> all 0.01 |---------------------------| 20.946
#> 
#> - Units with the 5 most extreme weights:
#>                                        
#>         485    484    483    482    481
#>  all 10.209 13.112 13.974 17.816 20.946
#> 
#> - Weight statistics:
#> 
#>     Coef of Var   MAD Entropy # Zeros
#> all       1.454 0.535   0.396       0
#> 
#> - Effective Sample Sizes:
#> 
#>             Total
#> Unweighted 614.  
#> Weighted   197.36
bal.tab(W3)
#> Balance Measures
#>             Type Corr.Adj
#> age      Contin.       -0
#> educ     Contin.        0
#> married   Binary       -0
#> nodegree  Binary       -0
#> re74     Contin.        0
#> 
#> Effective sample sizes
#>             Total
#> Unadjusted 614.  
#> Adjusted   197.36
```
