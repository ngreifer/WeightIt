# Propensity Score Weighting Using Generalized Linear Models

This page explains the details of estimating weights from generalized
linear model-based propensity scores by setting `method = "glm"` in the
call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
This method can be used with binary, multi-category, and continuous
treatments, as well as for estimating censoring weights.

In general, this method relies on estimating propensity scores with a
parametric generalized linear model and then converting those propensity
scores into weights using a formula that depends on the desired
estimand. For binary and multi-category treatments, a binomial or
multinomial regression model (logistic, by default) is used to estimate
the propensity scores as the predicted probability of being in each
treatment given the covariates. For ordinal treatments, an ordinal
regression model is used to estimate generalized propensity scores. For
continuous treatments, a generalized linear model is used to estimate
generalized propensity scores as the conditional density of treatment
given the covariates.

### Binary Treatments

For binary treatments, this method estimates the propensity scores using
[`glm()`](https://rdrr.io/r/stats/glm.html). The following estimands are
allowed: ATE, ATT, ATC, ATO, ATM, and ATOS. Weights can also be computed
using marginal mean weighting through stratification for the ATE, ATT,
and ATC. See
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
for details.

### Multi-Category Treatments

For multi-category treatments, the propensity scores are estimated using
multinomial regression from one of a few functions depending on the
argument supplied to `multi.method` (see Additional Arguments below).
The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The
weights for each estimand are computed using the standard formulas.
Weights can also be computed using marginal mean weighting through
stratification for the ATE, ATT, and ATC. See
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
for details. Ordinal treatments are treated exactly the same as
non-ordinal multi-category treatments except that additional models are
available to estimate the generalized propensity score (e.g., ordinal
logistic regression).

### Continuous Treatments

For continuous treatments, weights are estimated as \\w_i = f_A(a_i) /
f\_{A\|X}(a_i)\\, where \\f_A(a_i)\\ (known as the stabilization factor)
is the unconditional density of treatment evaluated the observed
treatment value and \\f\_{A\|X}(a_i)\\ (known as the generalized
propensity score) is the conditional density of treatment given the
covariates evaluated at the observed treatment value. The shape of
\\f\_{A\|X}(.)\\ is controlled by the `density` argument described below
(normal distribution by default), and the predicted values used for the
mean of the conditional density are estimated using linear regression.
\\f_A(.)\\ is estimated by marginalizing over \\f\_{A\|X}(.)\\. Kernel
density estimation can be used instead of assuming a specific density
for the denominator by setting `density = "kernel"`. Other arguments to
[`density()`](https://rdrr.io/r/stats/density.html) can be specified to
refine the density estimation parameters.

### Multilevel Treatment Models

When the model `formula` contains
[lme4](https://CRAN.R-project.org/package=lme4)-style random effects
terms (e.g., `treat ~ x1 + x2 + (1 | school)`), a multilevel
(mixed-effects) model is used to estimate the propensity scores. This
can improve balance and overlap when units are clustered (e.g., patients
within hospitals or students within schools). The grouping (and any
random-slope) variables are taken from `data`. For binary treatments,
the model is fit using
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) with the
requested `link` (`"logit"`, `"probit"`, or `"cloglog"`); for continuous
treatments, using
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) or
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) ; and for
multi-category treatments, using
[`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html)
with its `random` argument (i.e., `multi.method` is set to `"mclogit"`).
The propensity scores are the predicted probabilities (or conditional
densities) that include the estimated random effects, i.e.,
cluster-specific predictions. M-estimation is not available for these
models; robust (`HC0`) or bootstrap standard errors should be used
instead when estimating treatment effects (see the M-estimation section
below and
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)).

### Censoring Weights

For censoring weights, requested by wrapping the censoring indicator in
[`.cens()`](https://ngreifer.github.io/WeightIt/reference/dot-cens.md),
a single binomial regression model of the probability of being censored
is fit, and the weights are `1/P(C = 0 | X)` for the units still under
observation and 0 for the censored units. All of the `link` options,
`missing` options, and multilevel models described above apply
unchanged, and the returned propensity score is the probability of
*being censored*.

### Longitudinal Treatments

For longitudinal treatments, the weights are the product of the weights
estimated at each time point.

### Sampling Weights

Sampling weights are supported through `s.weights` in all scenarios
except for multi-category treatments with `multi.method = "mnp"` and for
binary and continuous treatments with `missing = "saem"` (see below).
Warning messages may appear otherwise about non-integer successes, and
these can be ignored.

### Missing Data

In the presence of missing data, the following value(s) for `missing`
are allowed:

- `"ind"` (default):

  First, for each variable with missingness, a new missingness indicator
  variable is created which takes the value 1 if the original covariate
  is `NA` and 0 otherwise. The missingness indicators are added to the
  model formula as main effects. The missing values in the covariates
  are then replaced with the covariate medians (this value is arbitrary
  and does not affect estimation). The weight estimation then proceeds
  with this new formula and set of covariates. The covariates output in
  the resulting `weightit` object will be the original covariates with
  the `NA`s.

- `"saem"`:

  For binary treatments with `link = "logit"` or continuous treatments,
  a stochastic approximation version of the EM algorithm (SAEM) is used
  via the [misaem](https://CRAN.R-project.org/package=misaem) package.
  No additional covariates are created. See Jiang et al. (2019) for
  information on this method. In some cases, this is a suitable
  alternative to multiple imputation.

### M-estimation

For binary treatments, M-estimation is supported when `link` is neither
`"flic"` nor `"flac"` (see below). For multi-category treatments,
M-estimation is supported when `multi.method` is `"weightit"` (the
default) or `"glm"`. M-estimation is not supported when `subclass` is
specified. For continuous treatments, M-estimation is supported when
`density` is not `"kernel"`. The conditional treatment variance and
unconditional treatment mean and variance are included as parameters to
estimate, as these all go into calculation of the weights. For all
treatment types, M-estimation is not supported when `missing = "saem"`
or when the model `formula` includes random effects terms (i.e., a
multilevel model is fit; see *Multilevel Treatment Models* above). See
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
and
[`vignette("estimating-effects")`](https://ngreifer.github.io/WeightIt/articles/estimating-effects.md)
for details. For longitudinal treatments, M-estimation is supported
whenever the underlying methods are.

## Additional Arguments

For binary treatments, the following additional argument can be
specified:

- `link`:

  the link used in the generalized linear model for the propensity
  scores. `link` can be any of those allowed by
  [`binomial()`](https://rdrr.io/r/stats/family.html) as well as
  `"loglog"` and `"clog"`. A `br.` prefix can be added (e.g.,
  `"br.logit"`); this changes the fitting method to the bias-corrected
  generalized linear models implemented in the
  [brglm2](https://CRAN.R-project.org/package=brglm2) package. `link`
  can also be either `"flic"` or `"flac"` to fit the corresponding Firth
  corrected logistic regression models implemented in the
  [logistf](https://CRAN.R-project.org/package=logistf) package.

- `subclass`:

  `integer`; the number of subclasses to use for computing weights using
  marginal mean weighting through stratification (MMWS). If `NULL`,
  standard inverse probability weights (and their extensions) will be
  computed; if a number greater than 1, subclasses will be formed and
  weights will be computed based on subclass membership. See
  [`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
  for details and references.

For multi-category treatments, the following additional arguments can be
specified:

- `multi.method`:

  the method used to estimate the generalized propensity scores.
  Allowable options include `"weightit"` (the default) to use
  multinomial logistic regression implemented in WeightIt, `"glm"` to
  use a series of binomial models using
  [`glm()`](https://rdrr.io/r/stats/glm.html), `"mclogit"` to use
  multinomial logistic regression as implemented in
  [`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html)
  , and `"mnp"` to use Bayesian multinomial probit regression as
  implemented in [`MNP::MNP()`](https://rdrr.io/pkg/MNP/man/mnp.html) .
  `"weightit"` and `"mclogit"` should give near-identical results, the
  main difference being increased robustness and customizability when
  using `"mclogit"` at the expense of not being able to use M-estimation
  to compute standard errors after weighting. For ordered treatments,
  allowable options include `"weightit"` (the default) to use ordinal
  regression implemented in WeightIt or `"polr"` to use ordinal
  regression implemented in
  [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) . Ignored
  when `missing = "saem"`. Using the defaults allows for the use of
  M-estimation and requires no additional dependencies, but other
  packages may provide benefits such as speed and flexibility.

  Bias-reduced models are requested by adding a `br.` prefix to `link`
  rather than by supplying a `multi.method` (see below). Because
  [`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html) cannot fit
  bias-reduced models, supplying a `br.` link with
  `multi.method = "polr"` uses the WeightIt implementation instead,
  which also makes M-estimation available.

- `link`:

  The link used in the multinomial, binomial, or ordered regression
  model for the generalized propensity scores depending on the argument
  supplied to `multi.method`. When `multi.method = "glm"`, `link` can be
  any of those allowed by
  [`binomial()`](https://rdrr.io/r/stats/family.html) as well as
  `"loglog"` and `"clog"`. When treatment is ordered and `multi.method`
  is `"weightit"` or `"polr"`, `link` can be any of `"logit"`,
  `"probit"`, `"loglog"`, `"cloglog"`, or `"cauchit"`. Otherwise, `link`
  should be `"logit"` or not specified.

  For unordered treatments with `multi.method = "weightit"` and for
  ordered treatments with `multi.method` of `"weightit"` or `"polr"`, a
  `br.` prefix can be added (e.g., `"br.logit"`) to request mean bias
  reduction, i.e., to solve the bias-reducing adjusted score equations
  of Firth (1993) rather than the score equations. See
  [`multinom_weightit()`](https://ngreifer.github.io/WeightIt/reference/multinom_weightit.md)
  and
  [`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/ordinal_weightit.md)
  for details. The estimates are always finite, even when the maximum
  likelihood estimates are not, which can help when the treatment groups
  are nearly separated. M-estimation remains available.

- `subclass`:

  `integer`; the number of subclasses to use for computing weights using
  marginal mean weighting through stratification (MMWS). If `NULL`,
  standard inverse probability weights (and their extensions) will be
  computed; if a number greater than 1, subclasses will be formed and
  weights will be computed based on subclass membership. See
  [`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)
  for details and references.

For continuous treatments, the following additional arguments may be
supplied:

- `density`:

  A function corresponding the conditional density of the treatment. The
  standardized residuals of the treatment model will be fed through this
  function to produce the denominator of the generalized propensity
  score weights. If blank,
  [`dnorm()`](https://rdrr.io/r/stats/Normal.html) is used as
  recommended by Robins et al. (2000). This can also be supplied as a
  string containing the name of the function to be called. If the string
  contains underscores, the call will be split by the underscores and
  the latter splits will be supplied as arguments to the second argument
  and beyond. For example, if `density = "dt_2"` is specified, the
  density used will be that of a t-distribution with 2 degrees of
  freedom. Using a t-distribution can be useful when extreme outcome
  values are observed (Naimi et al., 2014).

  Can also be `"kernel"` to use kernel density estimation, which calls
  [`density()`](https://rdrr.io/r/stats/density.html) to estimate the
  denominator density for the weights. (This used to be requested by
  setting `use.kernel = TRUE`, which is now deprecated.)

- `bw`, `adjust`, `kernel`, `n`:

  If `density = "kernel"`, the arguments to
  [`density()`](https://rdrr.io/r/stats/density.html). The defaults are
  the same as those in
  [`density()`](https://rdrr.io/r/stats/density.html).

- `link`:

  The link used to fit the linear model for the generalized propensity
  score. Can be any allowed by
  [`gaussian()`](https://rdrr.io/r/stats/family.html).

Additional arguments to [`glm()`](https://rdrr.io/r/stats/glm.html) can
be specified as well when it is used for fitting. The `method` argument
in [`glm()`](https://rdrr.io/r/stats/glm.html) is renamed to
`glm.method`. This can be used to supply alternative fitting functions,
such as those implemented in the
[glm2](https://CRAN.R-project.org/package=glm2) package. Other arguments
to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
are passed to `...` in [`glm()`](https://rdrr.io/r/stats/glm.html). In
the presence of missing data with `link = "logit"` and
`missing = "saem"`, additional arguments are passed to
[`misaem::miss.glm()`](https://rdrr.io/pkg/misaem/man/miss.glm.html) and
[`misaem::predict.miss.glm()`](https://rdrr.io/pkg/misaem/man/predict.miss.glm.html)
, except the `method` argument in
[`misaem::predict.miss.glm()`](https://rdrr.io/pkg/misaem/man/predict.miss.glm.html)
is replaced with `saem.method`.

For continuous treatments in the presence of missing data with
`missing = "saem"`, additional arguments are passed to
[`misaem::miss.lm()`](https://rdrr.io/pkg/misaem/man/miss.lm.html) and
[`misaem::predict.miss.lm()`](https://rdrr.io/pkg/misaem/man/predict.miss.lm.html)
.

When the model `formula` includes random effects terms (see *Multilevel
Treatment Models* above), additional arguments are passed to the
corresponding fitting function:
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) for binary
treatments,
[`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html)
for multi-category treatments, and
[`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) or
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) for
continuous treatments.

## Additional Outputs

- `obj`:

  When `include.obj = TRUE`, the (generalized) propensity score model
  fit. For binary treatments, the output of the call to
  [`glm()`](https://rdrr.io/r/stats/glm.html) or the requested fitting
  function. For multi-category treatments, the output of the call to the
  fitting function (or a list thereof if `multi.method = "glm"`). For
  continuous treatments, the output of the call to
  [`glm()`](https://rdrr.io/r/stats/glm.html) for the predicted values
  in the denominator density. When the model `formula` includes random
  effects terms, the output of the call to
  [`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) ,
  [`lme4::lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) , or
  [`mclogit::mblogit()`](https://melff.github.io/mclogit/reference/mblogit.html)
  .

## References

### Binary treatments

- `estimand = "ATO"`

Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates
via propensity score weighting. *Journal of the American Statistical
Association*, 113(521), 390–400.
[doi:10.1080/01621459.2016.1260466](https://doi.org/10.1080/01621459.2016.1260466)

- `estimand = "ATM"`

Li, L., & Greene, T. (2013). A Weighting Analogue to Pair Matching in
Propensity Score Analysis. *The International Journal of Biostatistics*,
9(2). [doi:10.1515/ijb-2012-0030](https://doi.org/10.1515/ijb-2012-0030)

- `estimand = "ATOS"`

Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009).
Dealing with limited overlap in estimation of average treatment effects.
*Biometrika*, 96(1), 187–199.
[doi:10.1093/biomet/asn055](https://doi.org/10.1093/biomet/asn055)

- Other estimands

Austin, P. C. (2011). An Introduction to Propensity Score Methods for
Reducing the Effects of Confounding in Observational Studies.
*Multivariate Behavioral Research*, 46(3), 399–424.
[doi:10.1080/00273171.2011.568786](https://doi.org/10.1080/00273171.2011.568786)

- Marginal mean weighting through stratification

Hong, G. (2010). Marginal mean weighting through stratification:
Adjustment for selection bias in multilevel data. *Journal of
Educational and Behavioral Statistics*, 35(5), 499–531.
[doi:10.3102/1076998609359785](https://doi.org/10.3102/1076998609359785)

- Bias-reduced regression

Firth, D. (1993). Bias reduction of maximum likelihood estimates.
*Biometrika*, 80(1), 27–38.
[doi:10.1093/biomet/80.1.27](https://doi.org/10.1093/biomet/80.1.27)

For binary treatments, see also the references for the
[brglm2](https://CRAN.R-project.org/package=brglm2) package, which does
the fitting. For multi-category treatments, the fitting is done by
WeightIt; see
[`multinom_weightit()`](https://ngreifer.github.io/WeightIt/reference/multinom_weightit.md)
and
[`ordinal_weightit()`](https://ngreifer.github.io/WeightIt/reference/ordinal_weightit.md),
which document the adjustments used.

- Firth corrected logistic regression

Puhr, R., Heinze, G., Nold, M., Lusa, L., & Geroldinger, A. (2017).
Firth’s logistic regression with rare events: Accurate effect estimates
and predictions? *Statistics in Medicine*, 36(14), 2302–2317.
[doi:10.1002/sim.7273](https://doi.org/10.1002/sim.7273)

- SAEM logistic regression for missing data

Jiang, W., Josse, J., & Lavielle, M. (2019). Logistic regression with
missing covariates — Parameter estimation, model selection and
prediction within a joint-modeling framework. *Computational Statistics
& Data Analysis*, 106907.
[doi:10.1016/j.csda.2019.106907](https://doi.org/10.1016/j.csda.2019.106907)

### Multi-Category Treatments

- `estimand = "ATO"`

Li, F., & Li, F. (2019). Propensity score weighting for causal inference
with multiple treatments. *The Annals of Applied Statistics*, 13(4),
2389–2415.
[doi:10.1214/19-AOAS1282](https://doi.org/10.1214/19-AOAS1282)

- `estimand = "ATM"`

Yoshida, K., Hernández-Díaz, S., Solomon, D. H., Jackson, J. W., Gagne,
J. J., Glynn, R. J., & Franklin, J. M. (2017). Matching weights to
simultaneously compare three treatment groups: Comparison to three-way
matching. *Epidemiology* (Cambridge, Mass.), 28(3), 387–395.
[doi:10.1097/EDE.0000000000000627](https://doi.org/10.1097/EDE.0000000000000627)

- Other estimands

McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E.,
Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score
Estimation for Multiple Treatments Using Generalized Boosted Models.
*Statistics in Medicine*, 32(19), 3388–3414.
[doi:10.1002/sim.5753](https://doi.org/10.1002/sim.5753)

- Marginal mean weighting through stratification

Hong, G. (2012). Marginal mean weighting through stratification: A
generalized method for evaluating multivalued and multiple treatments
with nonexperimental data. *Psychological Methods*, 17(1), 44–60.
[doi:10.1037/a0024918](https://doi.org/10.1037/a0024918)

### Continuous treatments

Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural
Models and Causal Inference in Epidemiology. *Epidemiology*, 11(5),
550–560.

- Using non-normal conditional densities

Naimi, A. I., Moodie, E. E. M., Auger, N., & Kaufman, J. S. (2014).
Constructing Inverse Probability Weights for Continuous Exposures: A
Comparison of Methods. *Epidemiology*, 25(2), 292–299.
[doi:10.1097/EDE.0000000000000053](https://doi.org/10.1097/EDE.0000000000000053)

- SAEM linear regression for missing data

Jiang, W., Josse, J., & Lavielle, M. (2019). Logistic regression with
missing covariates — Parameter estimation, model selection and
prediction within a joint-modeling framework. *Computational Statistics
& Data Analysis*, 106907.
[doi:10.1016/j.csda.2019.106907](https://doi.org/10.1016/j.csda.2019.106907)

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md),
[`get_w_from_ps()`](https://ngreifer.github.io/WeightIt/reference/get_w_from_ps.md)

## Examples

``` r
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "glm", estimand = "ATT",
                link = "probit"))
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
#> treated 1.                    ||            1.   
#> control 0.018 |---------------------------| 1.834
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>              5     4     3     2     1
#>  treated     1     1     1     1     1
#>            612   595   269   409   296
#>  control 1.278 1.351 1.412 1.518 1.834
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.    0.      0.          0
#> control       0.804 0.691   0.322       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    260.83     185

bal.tab(W1)
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.0252
#> age         Contin.   0.0716
#> educ        Contin.  -0.0565
#> married      Binary   0.0058
#> nodegree     Binary   0.0128
#> re74        Contin.  -0.0507
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    260.83     185

#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "glm", estimand = "ATE"))
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
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
#> black  1.453 |---------------------------| 30.807
#> hispan 1.799  |---------------------|      25.953
#> white  1.112 |-|                            3.977
#> 
#> - Units with the 5 most extreme weights by group:
#>                                           
#>            226    231    485    181    182
#>   black  7.532  7.624  8.525 12.351 30.807
#>            392    564    345    269    371
#>  hispan 15.992 16.752 18.784 23.411 25.953
#>            432    589    437    404    599
#>   white  3.686  3.686  3.712  3.784  3.977
#> 
#> - Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.89  0.426   0.192       0
#> hispan       0.541 0.404   0.132       0
#> white        0.382 0.317   0.068       0
#> 
#> - Effective Sample Sizes:
#> 
#>            black hispan  white
#> Unweighted 243.   72.   299.  
#> Weighted   135.8  55.86 261.02

bal.tab(W2)
#> Balance summary across all treatment pairs
#>             Type Max.Diff.Adj
#> age      Contin.       0.0419
#> educ     Contin.       0.1276
#> married   Binary       0.0500
#> nodegree  Binary       0.0605
#> re74     Contin.       0.2023
#> 
#> Effective sample sizes
#>            black hispan  white
#> Unadjusted 243.   72.   299.  
#> Adjusted   135.8  55.86 261.02

#Balancing covariates with respect to re75 (continuous)
#with kernel density estimate
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "glm", density = "kernel"))
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: continuous
#>  - covariates: age, educ, married, nodegree, re74

summary(W3)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>       Min                                  Max
#> all 0.043 |---------------------------| 44.225
#> 
#> - Units with the 5 most extreme weights:
#>                                        
#>         482    481    484    483    485
#>  all 30.586 30.632 37.608 41.195 44.225
#> 
#> - Weight statistics:
#> 
#>     Coef of Var   MAD Entropy # Zeros
#> all       2.405 0.933    0.92       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Total
#> Unweighted 614. 
#> Weighted    90.6

bal.tab(W3)
#> Balance Measures
#>             Type Corr.Adj
#> age      Contin.  -0.0712
#> educ     Contin.   0.0100
#> married   Binary  -0.0834
#> nodegree  Binary  -0.0186
#> re74     Contin.  -0.1429
#> 
#> Effective sample sizes
#>            Total
#> Unadjusted 614. 
#> Adjusted    90.6
```
