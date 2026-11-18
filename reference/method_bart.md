# Propensity Score Weighting Using BART

This page explains the details of estimating weights from Bayesian
additive regression trees (BART)-based propensity scores by setting
`method = "bart"` in the call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
This method can be used with binary, multi-category, and continuous
treatments, as well as for estimating censoring weights.

In general, this method relies on estimating propensity scores using
BART and then converting those propensity scores into weights using a
formula that depends on the desired estimand. This method relies on
[`bartisan::bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.html)
from the [bartisan](https://CRAN.R-project.org/package=bartisan)
package, which fits BART models with soft decision rules and with
likelihoods beyond the normal using the sampler of Linero (2025).

### Binary Treatments

For binary treatments, this method estimates the propensity scores using
[`bartisan::bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.html)
with a binomial likelihood, by default with a logit link (see `link` in
*Additional Arguments* below). The following estimands are allowed: ATE,
ATT, ATC, ATO, ATM, and ATOS. Weights can also be computed using
marginal mean weighting through stratification for the ATE, ATT, and
ATC. See
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
for details.

### Multi-Category Treatments

For multi-category treatments, the propensity scores are estimated using
a single call to
[`bartisan::bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.html)
with a multinomial likelihood
([`bartisan::multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.html)
), by default with a logit link (see `link` in *Additional Arguments*
below). One forest is fit for each treatment level, and the predicted
probabilities for each unit sum to 1. The following estimands are
allowed: ATE, ATT, ATC, ATO, and ATM. The weights for each estimand are
computed using the standard formulas. Weights can also be computed using
marginal mean weighting through stratification for the ATE, ATT, and
ATC. See
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
for details.

### Continuous Treatments

For continuous treatments, weights are estimated as \\w_i = f_A(a_i) /
f\_{A\|X}(a_i)\\, where \\f_A(a_i)\\ (known as the stabilization factor)
is the unconditional density of treatment evaluated at the observed
treatment value and \\f\_{A\|X}(a_i)\\ (known as the generalized
propensity score) is the conditional density of treatment given the
covariates evaluated at the observed treatment value. The conditional
density is estimated directly by a BART model for the treatment whose
form is set by the `family` argument (see *Additional Arguments* below).
By default, this is
[`bartisan::dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.html)
, a model whose mean is a sum of trees and whose error distribution is a
Dirichlet process mixture of normals estimated from the data (George et
al., 2019), so no parametric form for the density needs to be assumed.
The error distribution is the same for all units, though, so when the
variance of the treatment may depend on the covariates, a different
`family` may be preferable, such as `family = "gaussian_ls"`, which fits
a normal model in which both the mean and the logarithm of the standard
deviation are sums of trees. \\f_A(.)\\ is estimated by marginalizing
over \\f\_{A\|X}(.)\\.

### Multilevel Treatment Models

When the model `formula` contains
[lme4](https://CRAN.R-project.org/package=lme4)-style random intercept
terms (e.g., `treat ~ x1 + x2 + (1 | school)`), the random intercepts
are included in the BART model directly; see *Random Intercepts* at
[`bartisan::bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.html)
. This can improve balance when units are clustered (e.g., patients
within hospitals or students within schools). Multiple grouping factors
(e.g., `(1 | school) + (1 | district)`) and nesting
(`(1 | district/school)`) are allowed, but random slopes are not; a
variable whose effect varies by group belongs in the fixed part of the
formula, where a tree can split on it and on the group together. The
grouping variables are taken from `data`. The estimated propensity
scores (or conditional densities for continuous treatments) include the
estimated random effects; i.e., they are cluster-specific.

### Censoring Weights

For censoring weights, requested by wrapping the censoring indicator in
[`.cens()`](https://ngreifer.github.io/WeightIt/reference/dot-cens.md),
a single BART model of the probability of being censored is fit, and the
weights are `1/P(C = 0 | X)` for the units still under observation and 0
for the censored units. The returned propensity score is the probability
of *being censored*.

### Longitudinal Treatments

For longitudinal treatments, the weights are the product of the weights
estimated at each time point.

### Sampling Weights

Sampling weights are supported through `s.weights`, except with the
default `family` for continuous treatments. Before being supplied to
[`bartisan::bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.html)
as prior weights, they are scaled to have a mean of 1 among the units in
the model, so that the fit does not depend on their scale. For
continuous treatments, they are also used in computing the marginal
density of the treatment in the numerator of the weights. The default
`family` for continuous treatments,
[`bartisan::dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.html)
, does not accept prior weights, so supplying sampling weights that vary
across units with it is an error. To use sampling weights with a
continuous treatment, supply a family that accepts them, such as
`family = gaussian()` or `family = "gaussian_ls"`.

### Missing Data

In the presence of missing data, the following value(s) for `missing`
are allowed:

- `"mia"` (default):

  The missing values are left in the covariates and handled by the trees
  themselves using missingness incorporated in attributes; see *Missing
  Predictor Values* at
  [`bartisan::bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.html)
  . Each splitting rule on a covariate with missing values also
  determines which side of the split the units missing that covariate go
  to, or splits on whether the covariate is missing at all, and these
  choices are made by the sampler along with the rest of the tree. The
  estimated propensity scores are then conditional on the observed
  covariate values and on which covariates are missing. With
  `use.offset = TRUE`, the generalized linear model that supplies the
  offset uses missingness indicators as described for `"ind"` below. The
  covariates output in the resulting `weightit` object will be the
  original covariates with the `NA`s.

- `"ind"`:

  First, for each variable with missingness, a new missingness indicator
  variable is created which takes the value 1 if the original covariate
  is `NA` and 0 otherwise. The missingness indicators are added to the
  model formula as main effects. The missing values in the covariates
  are then replaced with the covariate medians. For a categorical
  covariate, the missing values are instead given a category of their
  own. The weight estimation then proceeds with this new formula and set
  of covariates. The covariates output in the resulting `weightit`
  object will be the original covariates with the `NA`s.

### M-estimation

M-estimation is not supported.

## Details

BART works by fitting a sum-of-trees model for the treatment or
probability of treatment. The number of trees is determined by the
`num_trees` argument. Bayesian priors are used for the hyperparameters,
so the result is a posterior distribution of predicted values for each
unit. The mean of these for each unit is taken for use in computing the
(generalized) propensity score. Although the hyperparameters governing
the priors can be modified by supplying arguments to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
that are passed to the BART fitting function, the default values tend to
work well and require little modification. Unlike many other machine
learning methods, no loss function is optimized and the hyperparameters
do not need to be tuned (e.g., using cross-validation), though
performance can benefit from tuning. BART tends to balance sparseness
with flexibility by using very weak learners as the trees, which makes
it suitable for capturing complex functions without specifying a
particular functional form and without overfitting.

The covariates are supplied to
[`bartisan::bartisan()`](https://ngreifer.github.io/bartisan/reference/bartisan.html)
as they appear in the model formula rather than as a model matrix, so a
categorical covariate remains a single variable whose levels a tree can
split into any two groups (see `categorical` at
[`bartisan::bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.html)
). Transformed terms (e.g., `I(x^2)`) are included as given, and an
interaction term is represented by its component variables, since the
trees can find interactions among them without one being specified.

By default, the decision rules in each tree are soft, as in the SoftBart
model of Linero and Yang (2018): an observation reaches every leaf with
some weight, so the fitted propensity score is a smooth function of the
covariates rather than a step function. Smoother propensity scores tend
to be better calibrated and to yield less extreme weights than those
from standard BART. Soft rules are slower than hard rules; the standard
BART sampler can be requested with `gate = "hard"`. See
[`bartisan::bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.html)
for details.

### Reproducibility

BART has a random component. bartisan uses R's random number generator,
so calling [`set.seed()`](https://rdrr.io/r/base/Random.html) before
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
makes the results reproducible, including when the chains run in
parallel, since each chain is given its own random number stream. The
`seed` argument accepted by earlier versions of this method, which used
dbarts, is ignored with a warning.

### Parallel Processing

When `chains` is greater than 1 and the
[future.apply](https://CRAN.R-project.org/package=future.apply) package
is installed, the chains run in parallel under whatever
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
has been set. See *Parallelization with future* at
[`bartisan::bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.html)
.

## Additional Arguments

All arguments to
[`bartisan::bartisan_control()`](https://ngreifer.github.io/bartisan/reference/bartisan_control.html)
can be passed through
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
to control the sampler and the priors, including `num_trees`,
`num_draws`, `num_burn`, `chains`, `gate`, and `sparsity`. Their
defaults are those of `bartisan_control()`. For compatibility with
earlier versions of this method, which used dbarts, the following
arguments are translated to their `bartisan_control()` equivalents when
supplied: `n.trees` to `num_trees`, `n.samples` to `num_draws`, `n.burn`
to `num_burn`, and `n.chains` to `chains`. A natively named argument
takes precedence over its translation. `n.threads` is ignored (see
*Parallel Processing* below), and `seed` is ignored with a warning (see
*Reproducibility* below).

For binary and multi-category treatments and for censoring weights, only
the link of the model described above can be changed, using `link`;
`family` is ignored. For continuous treatments, the model is set by
`family`, described below.

The following additional arguments may be supplied:

- `use.offset`:

  `logical`; whether to use the linear predictor resulting from a
  generalized linear model as an offset to the BART model. If `TRUE`,
  this fits a binomial regression model with the link specified by
  `link` (for binary treatments) or a linear regression model (for
  continuous treatments) and supplies the linear predictor to the
  `offset` argument of `bartisan()`, overriding any supplied argument to
  `offset`. With `missing = "mia"`, this model uses missingness
  indicators for the covariates with missing values. Default is `FALSE`
  to omit this offset. Only allowed for binary and continuous
  treatments.

For binary and multi-category treatments, the following arguments may be
supplied:

- `link`:

  for binary treatments and censoring weights, the link function of the
  binomial model. Allowable options include `"logit"` (the default),
  `"probit"`, and `"cloglog"`. For multi-category treatments, the link
  function of the multinomial model, which can be `"logit"` (the
  default) or `"probit"`; with `"probit"`, the probabilities have no
  closed form and are computed by simulation (see
  [`bartisan::multinomial()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.html)
  ).

- `subclass`:

  `integer`; the number of subclasses to use for computing weights using
  marginal mean weighting through stratification (MMWS). If `NULL`,
  standard inverse probability weights (and their extensions) will be
  computed; if a number greater than 1, subclasses will be formed and
  weights will be computed based on subclass membership. See
  [`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
  for details and references.

For continuous treatments, the following arguments may be supplied:

- `family`:

  the model for the conditional density of the treatment, given as a
  family object or the name of one; see
  [`bartisan-families`](https://ngreifer.github.io/bartisan/reference/bartisan-families.html).
  Default is
  [`bartisan::dpm()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.html)
  , which estimates the shape of the error distribution from the data
  but does not accept sampling weights. Other options include
  [`bartisan::gaussian_ls()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.html)
  for a normal density whose mean and standard deviation both depend on
  the covariates and [`gaussian()`](https://rdrr.io/r/stats/family.html)
  for a normal density with a single standard deviation; both accept
  sampling weights. `gaussian_ls()` can capture a variance that differs
  across units but takes several times as long to fit as `dpm()`, and,
  like [`gaussian()`](https://rdrr.io/r/stats/family.html), it assumes
  the conditional density is normal, which can yield extreme weights
  when the treatment is skewed. With `family = "ordinal"`
  ([`bartisan::ordinal()`](https://ngreifer.github.io/bartisan/reference/bartisan-families.html)
  ), each observed value of the treatment is a category of an ordinal
  model, and the weights are the ratio of the proportion of units with
  each unit's treatment value (weighted by the sampling weights) to the
  probability of that value given the covariates, as for a discrete
  treatment.

The `density`, `bw`, `adjust`, `kernel`, and `n` arguments that other
methods use to specify the conditional density are not used by this
method and are ignored with a warning, since the density is estimated by
the model specified by `family`.

## Additional Outputs

- `obj`:

  When `include.obj = TRUE`, the `bartisan_fit` object used to generate
  the predicted values. For multi-category treatments, this is a single
  fit, since one multinomial model is fit for all treatment levels. The
  predicted probabilities used to compute the propensity scores can be
  extracted using
  [`bartisan::predict.bartisan_fit()`](https://ngreifer.github.io/bartisan/reference/predict.bartisan_fit.html)
  with `type = "prob"`, and the conditional density of a continuous
  treatment with `type = "density"`.

## References

Chipman, H. A., George, E. I., & McCulloch, R. E. (2010). BART: Bayesian
additive regression trees. *The Annals of Applied Statistics*, 4(1),
266–298. [doi:10.1214/09-AOAS285](https://doi.org/10.1214/09-AOAS285)

Hill, J., Weiss, C., & Zhai, F. (2011). Challenges With Propensity Score
Strategies in a High-Dimensional Setting and a Potential Alternative.
*Multivariate Behavioral Research*, 46(3), 477–513.
[doi:10.1080/00273171.2011.570161](https://doi.org/10.1080/00273171.2011.570161)

Linero, A. R. (2025). Generalized Bayesian additive regression trees
models: beyond conditional conjugacy. *Journal of the American
Statistical Association*, 120(549), 356–369.
[doi:10.1080/01621459.2024.2337156](https://doi.org/10.1080/01621459.2024.2337156)

Linero, A. R., & Yang, Y. (2018). Bayesian regression tree ensembles
that adapt to smoothness and sparsity. *Journal of the Royal Statistical
Society Series B*, 80(5), 1087–1110.
[doi:10.1111/rssb.12293](https://doi.org/10.1111/rssb.12293)

- Continuous treatments

George, E., Laud, P., Logan, B., McCulloch, R., & Sparapani, R. (2019).
Fully Nonparametric Bayesian Additive Regression Trees. In *Topics in
Identification, Limited Dependent Variables, Partial Observability,
Experimentation, and Flexible Modeling: Part B* (Advances in
Econometrics, Vol. 40B, pp. 89–110). Emerald Publishing Limited.
[doi:10.1108/S0731-90532019000040B006](https://doi.org/10.1108/S0731-90532019000040B006)

Note that many references that deal with BART for causal inference focus
on estimating potential outcomes with BART, not the propensity scores,
and so are not directly relevant when using BART to estimate propensity
scores for weights.

See
[`method_glm`](https://ngreifer.github.io/WeightIt/reference/method_glm.md)
for additional references on propensity score weighting more generally.

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md),
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)

[`method_super`](https://ngreifer.github.io/WeightIt/reference/method_super.md)
for stacking predictions from several machine learning methods,
including BART.

[`vignette("weighting-methods")`](https://ngreifer.github.io/WeightIt/articles/weighting-methods.md)
for descriptions of all the weighting methods and guidance on choosing
among them.

## Examples

``` r
data("lalonde", package = "cobalt")

set.seed(123)

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "bart", estimand = "ATT"))
#> A weightit object
#>  - method: "bart" (propensity score weighting with BART)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74

summary(W1)
#>                   Summary of weights
#> 
#> ─ Weight ranges:
#> 
#>           Min                                  Max
#> Treated 1.        │                          1.   
#> Control 0.006 ╞═══════════════════════════╡ 11.793
#> 
#> ─ Units with the 5 most extreme weights by group:
#>                                      
#>              5    4    3     2      1
#>  Treated     1    1    1     1      1
#>            592  409  585   374    608
#>  Control 3.239 3.39 3.56 4.387 11.793
#> 
#> ─ Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> Treated       0.    0.      0.          0
#> Control       2.152 0.972   0.855       0
#> 
#> ─ Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted     76.35     185

cobalt::bal.tab(W1)
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.3887
#> age         Contin.  -0.0103
#> educ        Contin.  -0.0315
#> married      Binary  -0.0215
#> nodegree     Binary   0.0374
#> re74        Contin.  -0.0518
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted     76.35     185

#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "bart", estimand = "ATE"))
#> A weightit object
#>  - method: "bart" (propensity score weighting with BART)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 3-category (black, hispan, white)
#>  - estimand: ATE
#>  - covariates: age, educ, married, nodegree, re74

summary(W2)
#>                   Summary of weights
#> 
#> ─ Weight ranges:
#> 
#>          Min                                  Max
#> black  1.235 ╞═══════════════╡              8.539
#> hispan 4.045        ╞════════════════════╡ 13.25 
#> white  1.118 ╞═══════════════╡              8.232
#> 
#> ─ Units with the 5 most extreme weights by group:
#>                                          
#>            181    209    231    244   423
#>   black  7.373  8.243  8.381  8.388 8.539
#>            426    564     28    180   570
#>  hispan 12.561 12.645 13.109 13.145 13.25
#>             68     60     23     76   140
#>   white  5.598  5.665  5.894   7.93 8.232
#> 
#> ─ Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.619 0.402   0.145       0
#> hispan       0.339 0.292   0.055       0
#> white        0.476 0.324   0.088       0
#> 
#> ─ Effective Sample Sizes:
#> 
#>             black hispan  white
#> Unweighted 243.    72.   299.  
#> Weighted   175.81  64.66 243.95

cobalt::bal.tab(W2)
#> 
#> Balance summary across all treatment pairs
#>             Type Max.Diff.Adj
#> age      Contin.       0.1917
#> educ     Contin.       0.2670
#> married   Binary       0.0391
#> nodegree  Binary       0.0428
#> re74     Contin.       0.0733
#> Effective sample sizes
#>             black hispan  white
#> Unadjusted 243.    72.   299.  
#> Adjusted   175.81  64.66 243.95

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "bart"))
#> A weightit object
#>  - method: "bart" (propensity score weighting with BART)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: continuous
#>  - covariates: age, educ, married, nodegree, re74

summary(W3)
#>                   Summary of weights
#> 
#> ─ Weight ranges:
#> 
#>       Min                                 Max
#> All 0.508 ╞═══════════════════════════╡ 1.928
#> 
#> ─ Units with the 5 most extreme weights:
#>                                  
#>       485   396   126   409   128
#>  All 1.66 1.681 1.851 1.915 1.928
#> 
#> ─ Weight statistics:
#> 
#>     Coef of Var   MAD Entropy # Zeros
#> All       0.199 0.137   0.018       0
#> 
#> ─ Effective Sample Sizes:
#> 
#>             Total
#> Unweighted 614.  
#> Weighted   590.69

cobalt::bal.tab(W3)
#> Balance Measures
#>             Type Corr.Adj
#> age      Contin.   0.1362
#> educ     Contin.   0.0122
#> married   Binary   0.3186
#> nodegree  Binary  -0.0677
#> re74     Contin.   0.4874
#> 
#> Effective sample sizes
#>             Total
#> Unadjusted 614.  
#> Adjusted   590.69

#Using a normal model whose variance depends on the
#covariates for the conditional density
(W4 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "bart",
                family = "gaussian_ls"))
#> Warning: The leaf scale settled 9 times above its prior median, which usually means the
#> response is close to separable by the predictors.
#> ℹ The additive predictor is then only weakly identified; fix the scale with
#>   `update_sigma_mu = FALSE` in `bartisan_control()` if the draws look unstable.
#> A weightit object
#>  - method: "bart" (propensity score weighting with BART)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: continuous
#>  - covariates: age, educ, married, nodegree, re74

summary(W4)
#>                   Summary of weights
#> 
#> ─ Weight ranges:
#> 
#>       Min                                   Max
#> All 0.007 ╞═══════════════════════════╡ 979.914
#> 
#> ─ Units with the 5 most extreme weights:
#>                                             
#>          490     487       8     484     563
#>  All 634.937 639.779 735.771 757.083 979.914
#> 
#> ─ Weight statistics:
#> 
#>     Coef of Var   MAD Entropy # Zeros
#> All       2.606 1.479   1.691       0
#> 
#> ─ Effective Sample Sizes:
#> 
#>             Total
#> Unweighted 614.  
#> Weighted    78.93

cobalt::bal.tab(W4)
#> Balance Measures
#>             Type Corr.Adj
#> age      Contin.  -0.0005
#> educ     Contin.  -0.0014
#> married   Binary  -0.0009
#> nodegree  Binary  -0.0005
#> re74     Contin.   0.0000
#> 
#> Effective sample sizes
#>             Total
#> Unadjusted 614.  
#> Adjusted    78.93
```
