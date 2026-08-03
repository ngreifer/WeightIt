# Estimate Balancing Weights

`weightit()` allows for the easy generation of balancing weights using a
variety of available methods for binary, continuous, and multi-category
treatments. Some of these methods require functions in other packages,
which `weightit()` calls; these packages must be installed to use the
desired method.

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
  unit's weight will be multiplied by a stabilization factor, which is
  the the unconditional probability (or density) of each unit's observed
  treatment value. If a formula, a generalized linear model will be fit
  with the included predictors, and the inverse of the corresponding
  weight will be used as the standardization factor. Can only be used
  when `estimand = "ATE"` or with continuous treatments. Default is
  `FALSE` for no stabilization. See also the `num.formula` argument at
  [`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
  For continuous treatments, weights are already stabilized, so setting
  `stabilize = TRUE` will be ignored with a warning (supplying a formula
  still works).

- focal:

  when `estimand` is set to `"ATT"` or `"ATC"`, which group to consider
  the "treated" or "control" group, respectively. This group will not be
  weighted, and the other groups will be weighted to resemble the focal
  group. If specified, `estimand` will automatically be set to `"ATT"`
  (with a warning if `estimand` is not `"ATT"` or `"ATC"`). See section
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

  an optional vector of sampling weights or the name of a variable in
  `data` that contains sampling weights. See the individual pages for
  each method for information on whether sampling weights can be
  supplied.

- ps:

  an optional vector of propensity scores or the name of a variable in
  `data` containing propensity scores. If supplied, `method` is ignored
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
  fitting function. Default is `FALSE` to suppress output.

- include.obj:

  `logical`; whether to include in the output any fit objects created in
  the process of estimating the weights. For example, with
  `method = "glm"`, the `glm` objects containing the propensity score
  model will be included. See the individual pages for each method for
  information on what object will be included if `TRUE`. Default is
  `FALSE` to keep the returned object small.

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
  one that estimates propensity scores. The propensity score corresponds
  to the predicted probability of being treated; see section *`estimand`
  and `focal`* in Details for how the treated group is determined.

- s.weights:

  The provided sampling weights, or a vector of 1s of none are provided.

- focal:

  The focal treatment level if the ATT or ATC was requested.

- by:

  A data frame containing the `by` variable when specified.

- obj:

  When `include.obj = TRUE`, the fit object.

- info:

  Additional information about the fitting. See the individual methods
  pages for what is included.

When `keep.mparts` is `TRUE` (the default) and the chosen method is
compatible with M-estimation, the components related to M-estimation for
use in
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
are stored in the `"Mparts"` attribute. When `by` is specified, the
per-stratum M-estimation components are instead combined and stored in
the `"Mparts.list"` attribute; the resulting standard errors produced by
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
are asymptotically equivalent to those from estimating the weights from
a single model in which the `by` variable is fully interacted with all
the covariates. The same is true for
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md),
where the equivalent model interacts the `by` variable with all the
covariates at every time point.

## Details

The primary purpose of `weightit()` is as a dispatcher to functions that
perform the estimation of balancing weights using the requested
`method`. Below are the methods allowed and links to pages containing
more information about them, including additional arguments and outputs
(e.g., when `include.obj = TRUE`), how missing values are treated, which
estimands are allowed, and whether sampling weights are allowed.

|  |  |
|----|----|
| `method` | Name |
| [`"glm"`](https://ngreifer.github.io/WeightIt/reference/method_glm.md) | Propensity score weighting using generalized linear models |
| [`"gbm"`](https://ngreifer.github.io/WeightIt/reference/method_gbm.md) | Propensity score weighting using generalized boosted modeling |
| [`"cbps"`](https://ngreifer.github.io/WeightIt/reference/method_cbps.md) | Covariate Balancing Propensity Score weighting |
| [`"npcbps"`](https://ngreifer.github.io/WeightIt/reference/method_npcbps.md) | Non-parametric Covariate Balancing Propensity Score weighting |
| [`"ebal"`](https://ngreifer.github.io/WeightIt/reference/method_ebal.md) | Entropy balancing |
| [`"ipt"`](https://ngreifer.github.io/WeightIt/reference/method_ipt.md) | Inverse probability tilting |
| [`"optweight"`](https://ngreifer.github.io/WeightIt/reference/method_optweight.md) | Stable balancing weights |
| [`"super"`](https://ngreifer.github.io/WeightIt/reference/method_super.md) | Propensity score weighting using SuperLearner |
| [`"bart"`](https://ngreifer.github.io/WeightIt/reference/method_bart.md) | Propensity score weighting using Bayesian additive regression trees (BART) |
| [`"energy"`](https://ngreifer.github.io/WeightIt/reference/method_energy.md) | Energy balancing |
| [`"cfd"`](https://ngreifer.github.io/WeightIt/reference/method_cfd.md) | Characteristic function distance balancing |

`method` can also be supplied as a user-defined function; see
[`method_user`](https://ngreifer.github.io/WeightIt/reference/method_user.md)
for instructions and examples. Setting `method = NULL` computes unit
weights.

### Empty model formulas

The right hand side of `formula` may be empty, as in `A ~ 1`, requesting
a marginal model in which the treatment (or censoring) is taken to be
independent of the covariates. With no covariates there is nothing for
any method to model or balance, and every method's target is met by the
same weights: the inverse of the marginal treatment probability, or
\\1/P(C = 0)\\ for a censoring model. Those weights are computed by
fitting an intercept-only generalized linear model whatever `method` is
supplied, which is simply the easiest way to get them; any
method-specific arguments are ignored, since none of them can apply to
covariates that do not exist. This is invisible: `method` is reported as
supplied, the package for the requested method need not be installed,
and the weights are what that method would have produced. Every method
therefore accepts an empty formula.

For a continuous treatment the conditional density of the treatment is
its marginal density, so all the weights are exactly 1 and nothing is
estimated; no M-estimation components are produced. For the other
treatment types the marginal probability is estimated, so M-estimation
is available as usual.

`estimand`, `focal`, `by`, `s.weights`, and `stabilize` are unaffected,
and the arguments are still checked against the requested `method`, so a
method that cannot handle the treatment type at all (e.g., `"npcbps"`
with a censoring model) still produces an error.

A formula whose only terms are *lme4*-style random effects, such as
`A ~ (1 | school)`, is *not* empty in this sense and is fit as the
multilevel model it describes.

### Censoring weights (IPCW)

Wrapping the left side of `formula` in
[`.cens()`](https://ngreifer.github.io/WeightIt/reference/dot-cens.md)
requests inverse probability of censoring weights instead of treatment
weights, as in `weightit(.cens(C) ~ x1 + x2, data = d, method = "glm")`.
Censoring is treated as its own treatment type, distinct from binary,
multi-category, and continuous treatments; the indicator must be 0 or
`FALSE` for units still under observation and 1 or `TRUE` for units that
are censored.

Weights are estimated only for the units still under observation, and
are those that make their covariate distribution resemble that of the
full at-risk sample. Writing \\e(X) = P(C = 1 \| X)\\, the weights are
\\1 / (1 - e(X))\\ for units with \\C = 0\\ and exactly 0 for units with
\\C = 1\\. Because only one group is weighted, the estimation problem is
smaller and better conditioned than the corresponding binary-treatment
problem, which would additionally solve for weights among the censored
units; this matters most when few units are censored.

`estimand`, `focal`, and `subclass` do not apply and are rejected or
ignored. `by` and `stabilize` can be used; stabilization multiplies the
weights by \\P(C = 0 \| V)\\ from a second censoring model (marginal
when `stabilize` is `TRUE`, otherwise fit with the predictors in the
supplied formula), giving \\P(C = 0 \| V) / P(C = 0 \| X)\\ for the
units still under observation and leaving the censored units at exactly
0. `ps` is the predicted probability of *being censored*. Not all
methods support estimating censoring weights; see the `treat_type`
component of
[`.weightit_methods`](https://ngreifer.github.io/WeightIt/reference/dot-weightit_methods.md).

As with a treatment model, the right side of the formula may be empty,
as in `.cens(C) ~ 1`, which requests a marginal censoring model that
assumes censoring is independent of the covariates. The resulting
weights are \\1/P(C = 0)\\ for the units still under observation and 0
for the censored units. The rest of the censoring machinery is
unaffected, so such a model can be combined with `by`, `s.weights`,
`stabilize` (which then makes all nonzero weights exactly 1), and
M-estimation, and it can be interleaved with covariate-dependent
censoring models in
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
See *Empty model formulas* above for how `method` is handled.

Because censored units receive a weight of exactly 0, they contribute
nothing to a weighted outcome model, and
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
and friends tolerate missing values in the model variables for those
units, including a missing event time in the `Surv()` response of a
[`coxph_weightit()`](https://ngreifer.github.io/WeightIt/reference/coxph_weightit.md)
model. Missing values in units with a nonzero weight still produce an
error. See
[`.cens()`](https://ngreifer.github.io/WeightIt/reference/dot-cens.md)
for how to assess balance, which requires a little care.

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
or with both treatment and censoring indicators for use in estimating
marginal structural models (MSMs).

[`weightit.fit()`](https://ngreifer.github.io/WeightIt/reference/weightit.fit.md),
which is a lower-level dispatcher function that accepts a matrix of
covariates and a vector of treatment statuses rather than a formula and
data frame and performs minimal argument checking and processing. It may
be useful for speeding up simulation studies for which the correct
arguments are known. In general, `weightit()` should be used.

[`summary.weightit()`](https://ngreifer.github.io/WeightIt/reference/summary.weightit.md)
for summarizing the distribution of the weights.

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
#>             5     4   3     2     1
#>  treated    1     1   1     1     1
#>           411   595 269   409   296
#>  control 1.33 1.437 1.5 1.637 2.044
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.    0.       0.         0
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
#>          Min                                  Max
#> black  1.397 |-----------|                 13.517
#> hispan 1.201 |---------------------------| 28.417
#> white  0.817 |-|                            3.949
#> 
#> - Units with the 5 most extreme weights by group:
#>                                           
#>            226    244    485    181    182
#>   black  6.371  6.441   7.09  8.983 13.517
#>            392    564    269    345    371
#>  hispan 17.451 21.573 22.447 23.067 28.417
#>             68    457    599    589    531
#>   white  3.513  3.537   3.58  3.643  3.949
#> 
#> - Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.59  0.413   0.131       0
#> hispan       0.609 0.44    0.163       0
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
#>         485    481    482    484    483
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
#> married   Binary        0
#> nodegree  Binary       -0
#> re74     Contin.        0
#> 
#> Effective sample sizes
#>             Total
#> Unadjusted 614.  
#> Adjusted   197.36
```
