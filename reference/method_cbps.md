# Covariate Balancing Propensity Score Weighting

This page explains the details of estimating weights from covariate
balancing propensity scores by setting `method = "cbps"` in the call to
[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md)
or
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md).
This method can be used with binary, multi-category, and continuous
treatments.

In general, this method relies on estimating propensity scores using
generalized method of moments and then converting those propensity
scores into weights using a formula that depends on the desired
estimand. This method relies on code written for WeightIt using
[`optim()`](https://rdrr.io/r/stats/optim.html).

### Binary Treatments

For binary treatments, this method estimates the propensity scores and
weights using [`optim()`](https://rdrr.io/r/stats/optim.html) using
formulas described by Imai and Ratkovic (2014). The following estimands
are allowed: ATE, ATT, ATC, and ATO.

### Multi-Category Treatments

For multi-category treatments, this method estimates the generalized
propensity scores and weights using
[`optim()`](https://rdrr.io/r/stats/optim.html) using formulas described
by Imai and Ratkovic (2014). The following estimands are allowed: ATE
and ATT.

### Continuous Treatments

For continuous treatments, this method estimates the generalized
propensity scores and weights using
[`optim()`](https://rdrr.io/r/stats/optim.html) using a modification of
the formulas described by Fong, Hazlett, and Imai (2018). See Details.

### Longitudinal Treatments

For longitudinal treatments, the weights are computed using methods
similar to those described by Huffman and van Gameren (2018). This
involves specifying moment conditions for the models at each time point
as with single-time point treatments but using the product of the
time-specific weights as the weights that are used in the balance moment
conditions. This yields weights that balance the covariate at each time
point. This is not the same implementation as is implemented in
[`CBPS::CBMSM()`](https://rdrr.io/pkg/CBPS/man/CBMSM.html), and results
should not be expected to align between the two methods. Any combination
of treatment types is supported.

For the over-identified version (i.e., setting `over = TRUE`), the
empirical variance is used in the objective function, whereas the
expected variance averaging over the treatment is used with binary and
multi-category point treatments.

### Sampling Weights

Sampling weights are supported through `s.weights` in all scenarios.

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

### M-estimation

M-estimation is supported for the just-identified CBPS (the default,
setting `over = FALSE`) for binary and multi-category treatments.
Otherwise (i.e., for continuous or longitudinal treatments or when
`over = TRUE`), M-estimation is not supported. See
[`glm_weightit()`](https://ngreifer.github.io/WeightIt/reference/glm_weightit.md)
and
[`vignette("estimating-effects")`](https://ngreifer.github.io/WeightIt/articles/estimating-effects.md)
for details.

## Details

CBPS estimates the coefficients of a generalized linear model (for
binary treatments), multinomial logistic regression model (for
multi-category treatments), or linear regression model (for continuous
treatments) that is used to compute (generalized) propensity scores,
from which the weights are computed. It involves replacing (or
augmenting, in the case of the over-identified version) the standard
maximum likelihood score equations with the balance constraints in a
generalized method of moments estimation. The idea is to nudge the
estimation of the coefficients toward those that produce balance in the
weighted sample. The just-identified version (with `over = FALSE`) does
away with the maximum likelihood score equations for the coefficients so
that only the balance constraints are used, which will therefore produce
superior balance on the means (i.e., corresponding to the balance
constraints) for binary and multi-category treatments and linear terms
for continuous treatments than will the over-identified version.

Just-identified CBPS is very similar to entropy balancing and inverse
probability tilting. For the ATT, all three methods will yield identical
estimates. For other estimands, the results will differ.

Note that WeightIt provides different functionality from the CBPS
package in terms of the versions of CBPS available; for extensions to
CBPS (e.g., optimal CBPS and CBPS for instrumental variables), the CBPS
package may be preferred. Note that for longitudinal treatments,
[`CBPS::CBMSM()`](https://rdrr.io/pkg/CBPS/man/CBMSM.html) uses
different methods and produces different results from
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)
called with `method = "cbps"`.

## Note

This method used to rely on functionality in the CBPS package, but no
longer does. Slight differences may be found between the two packages in
some cases due to numerical imprecision (or, for continuous and
longitudinal treatments, due to a difference in the estimator). WeightIt
supports arbitrary numbers of groups for the multi-category CBPS and any
estimand, whereas CBPS only supports up to four groups and only the ATE.
The implementation of the just-identified CBPS for continuous treatments
also differs from that of CBPS, and departs slightly from that described
by Fong et al. (2018). The treatment mean and treatment variance are
treated as random parameters to be estimated and are included in the
balance moment conditions. In Fong et al. (2018), the treatment mean and
variance are fixed to their empirical counterparts. For continuous
treatments with the over-identified CBPS, WeightIt and CBPS use
different methods of specifying the GMM variance matrix, which may lead
to differing results.

Note that the default method differs between the two implementations; by
default WeightIt uses the just-identified CBPS, which is faster to fit,
yields better balance, and is compatible with M-estimation for
estimating the standard error of the treatment effect, whereas CBPS uses
the over-identified CBPS by default. However, both the just-identified
and over-identified versions are available in both packages.

When the rootSolve package is installed, the optimization process will
be slightly faster and more accurate because starting values are
provided by an initial call to
[`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html)
. However, the package is not required.

## Additional Arguments

The following additional arguments can be specified:

- `over`:

  `logical`; whether to request the over-identified CBPS, which combines
  the generalized linear model regression score equations (for binary
  treatments), multinomial logistic regression score equations (for
  multi-category treatments), or linear regression score equations (for
  continuous treatments) to the balance moment conditions. Default is
  `FALSE` to use the just-identified CBPS.

- `twostep`:

  `logical`; when `over = TRUE`, whether to use the two-step
  approximation to the generalized method of moments variance. Default
  is `TRUE`. Setting to `FALSE` increases computation time but may
  improve estimation. Ignored with a warning when `over = FALSE`.

- `link`:

  the link used in the generalized linear model for the propensity
  scores when treatment is binary. Default is `"logit"` for logistic
  regression, which is used in the original description of the method by
  Imai and Ratkovic (2014), but others are allowed, including
  `"probit"`, `"cauchit"`, `"cloglog"`, `"loglog"`, `"log"`, `"clog"`,
  and `"identity"`. Note that negative weights are possible with these
  last three and they should be used with caution. An object of class
  `"link-glm"` can also be supplied. The argument is passed to
  [`quasibinomial()`](https://rdrr.io/r/stats/family.html). Ignored for
  multi-category, continuous, and longitudinal treatments.

- `reltol`:

  the relative tolerance for convergence of the optimization. Passed to
  the `control` argument of
  [`optim()`](https://rdrr.io/r/stats/optim.html). Default is `1e-10`.

- `maxit`:

  the maximum number of iterations for convergence of the optimization.
  Passed to the `control` argument of
  [`optim()`](https://rdrr.io/r/stats/optim.html). Default is 1000 for
  binary and multi-category treatments and 10000 for continuous and
  longitudinal treatments.

- `solver`:

  the solver to use to estimate the parameters of the just-identified
  CBPS. Allowable options include `"multiroot"` to use
  [`rootSolve::multiroot()`](https://rdrr.io/pkg/rootSolve/man/multiroot.html)
  and `"optim"` to use
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html). `"multiroot"`
  is the default when rootSolve is installed, as it tends to be much
  faster and more accurate; otherwise, `"optim"` is the default and
  requires no dependencies. Regardless of `solver`, the output of
  [`optim()`](https://rdrr.io/r/stats/optim.html) is returned when
  `include.obj = TRUE` (see below). When `over = TRUE`, the parameter
  estimates of the just-identified CBPS are used as starting values for
  the over-identified CBPS.

- `moments`:

  `integer`; the highest power of each covariate to be balanced. For
  example, if `moments = 3`, each covariate, its square, and its cube
  will be balanced. Can also be a named vector with a value for each
  covariate (e.g., `moments = c(x1 = 2, x2 = 4)`). Values greater than 1
  for categorical covariates are ignored. Default is 1 to balance
  covariate means.

- `int`:

  `logical`; whether first-order interactions of the covariates are to
  be balanced. Default is `FALSE`.

- `quantile`:

  a named list of quantiles (values between 0 and 1) for each continuous
  covariate, which are used to create additional variables that when
  balanced ensure balance on the corresponding quantile of the variable.
  For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures
  the 25th, 50th, and 75th percentiles of `x1` in each treatment group
  will be balanced in the weighted sample. Can also be a single number
  (e.g., `.5`) or a vector (e.g., `c(.25, .5, .75)`) to request the same
  quantile(s) for all continuous covariates. Only allowed with binary
  and multi-category treatments.

## Additional Outputs

- `obj`:

  When `include.obj = TRUE`, the output of the final call to
  [`optim()`](https://rdrr.io/r/stats/optim.html) used to produce the
  model parameters. Note that because of variable transformations, the
  resulting parameter estimates may not be interpretable.

## References

### Binary treatments

Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score.
*Journal of the Royal Statistical Society: Series B (Statistical
Methodology)*, 76(1), 243–263.

### Multi-Category treatments

Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score.
*Journal of the Royal Statistical Society: Series B (Statistical
Methodology)*, 76(1), 243–263.

### Continuous treatments

Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity
score for a continuous treatment: Application to the efficacy of
political advertisements. *The Annals of Applied Statistics*, 12(1),
156–177. [doi:10.1214/17-AOAS1101](https://doi.org/10.1214/17-AOAS1101)

### Longitudinal treatments

Huffman, C., & van Gameren, E. (2018). Covariate Balancing Inverse
Probability Weights for Time-Varying Continuous Interventions. *Journal
of Causal Inference*, 6(2).
[doi:10.1515/jci-2017-0002](https://doi.org/10.1515/jci-2017-0002)

Note: one should not cite Imai & Ratkovic (2015) when using CBPS for
longitudinal treatments.

Some of the code was inspired by the source code of the CBPS package.

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`weightitMSM()`](https://ngreifer.github.io/WeightIt/reference/weightitMSM.md)

[method_ebal](https://ngreifer.github.io/WeightIt/reference/method_ebal.md)
and
[method_ipt](https://ngreifer.github.io/WeightIt/reference/method_ipt.md)
for entropy balancing and inverse probability tilting, which work
similarly.

## Examples

``` r
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1a <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cbps", estimand = "ATT"))
#> A weightit object
#>  - method: "cbps" (covariate balancing propensity score weighting)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74

summary(W1a)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                 Max
#> treated 1.                 ||               1.   
#> control 0.017 |---------------------------| 2.263
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>              1     2     3     4     5
#>  treated     1     1     1     1     1
#>            410   404   224   111    84
#>  control 1.464 1.485 1.576 1.743 2.263
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.000 0.000   0.000       0
#> control       0.839 0.707   0.341       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    252.12     185

cobalt::bal.tab(W1a)
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.0164
#> age         Contin.  -0.0000
#> educ        Contin.   0.0000
#> married      Binary  -0.0000
#> nodegree     Binary   0.0000
#> re74        Contin.  -0.0000
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    252.12     185

#Balancing covariates between treatment groups (binary)
#using over-identified CBPS with probit link
(W1b <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cbps", estimand = "ATT",
                over = TRUE, link = "probit"))
#> A weightit object
#>  - method: "cbps" (covariate balancing propensity score weighting)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, re74

summary(W1b)
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                 Max
#> treated 1.                  ||              1.   
#> control 0.012 |---------------------------| 2.053
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>              1     2     3     4     5
#>  treated     1     1     1     1     1
#>            410   404   224   111    84
#>  control 1.368 1.378 1.472 1.607 2.053
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated        0.00 0.000   0.000       0
#> control        0.81 0.693   0.326       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  429.       185
#> Weighted    259.24     185

cobalt::bal.tab(W1b)
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.0334
#> age         Contin.  -0.0021
#> educ        Contin.   0.0028
#> married      Binary   0.0011
#> nodegree     Binary   0.0041
#> re74        Contin.  -0.0284
#> 
#> Effective sample sizes
#>            Control Treated
#> Unadjusted  429.       185
#> Adjusted    259.24     185

#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cbps", estimand = "ATE"))
#> A weightit object
#>  - method: "cbps" (covariate balancing propensity score weighting)
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
#> black  1.501  |------------------|         17.966
#> hispan 1.631  |--------------------------| 24.561
#> white  1.131 |--|                           4.134
#> 
#> - Units with the 5 most extreme weights by group:
#>                                           
#>            203    164    163    153    152
#>   black  6.799  6.838  7.267  9.897 17.966
#>             67     43     39     36     28
#>  hispan 16.762 19.853 22.019 23.878 24.561
#>            270    195    190    172    169
#>   white  3.688  3.781  3.848  3.895  4.134
#> 
#> - Weight statistics:
#> 
#>        Coef of Var   MAD Entropy # Zeros
#> black        0.635 0.387   0.133       0
#> hispan       0.582 0.447   0.155       0
#> white        0.389 0.327   0.071       0
#> 
#> - Effective Sample Sizes:
#> 
#>             black hispan  white
#> Unweighted 243.    72.   299.  
#> Weighted   173.37  53.95 259.76

cobalt::bal.tab(W2)
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
#> Adjusted   173.37  53.95 259.76

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

cobalt::bal.tab(W3)
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
# \donttest{
#Longitudinal treatments
data("msmdata")
(W4 <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                        A_2 ~ X1_1 + X2_1 +
                          A_1 + X1_0 + X2_0),
                   data = msmdata,
                   method = "cbps"))
#> A weightitMSM object
#>  - method: "cbps" (covariate balancing propensity score weighting)
#>  - number of obs.: 7500
#>  - sampling weights: none
#>  - number of time points: 2 (A_1, A_2)
#>  - treatment:
#>     + time 1: 2-category
#>     + time 2: 2-category
#>  - covariates:
#>     + baseline: X1_0, X2_0
#>     + after time 1: X1_1, X2_1, A_1, X1_0, X2_0

summary(W4)
#>                         Time 1                        
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                   Max
#> treated 1.05  |---------------------------| 110.88 
#> control 1.239 |---------------------|        86.951
#> 
#> - Units with the 5 most extreme weights by group:
#>                                            
#>            2115   2177   2058   1611     90
#>  treated 50.114 50.114 56.079  96.45 110.88
#>            2778   1532   2516    832    586
#>  control 52.431 55.212 55.212 58.593 86.951
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       1.107 0.552   0.291       0
#> control       1.033 0.600   0.315       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted 3306.   4194.  
#> Weighted   1598.88 1884.92
#> 
#>                         Time 2                        
#>                   Summary of weights
#> 
#> - Weight ranges:
#> 
#>           Min                                  Max
#> treated 1.05  |-----------------------|      96.45
#> control 1.239 |---------------------------| 110.88
#> 
#> - Units with the 5 most extreme weights by group:
#>                                            
#>            2004   1954    460   1933     82
#>  treated 42.387 42.822 50.114 50.114  96.45
#>            1753   1689   1386    911    664
#>  control 55.212 56.079 58.593 86.951 110.88
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.986 0.588   0.292       0
#> control       1.165 0.583   0.329       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted 3701.      3799
#> Weighted   1570.46    1926
#> 

cobalt::bal.tab(W4)
#> Balance summary across all time points
#>      Times    Type Max.Diff.Adj
#> X1_0  1, 2 Contin.            0
#> X2_0  1, 2  Binary            0
#> X1_1     2 Contin.            0
#> X2_1     2  Binary            0
#> A_1      2  Binary            0
#> 
#> Effective sample sizes
#>  - Time 1
#>            Control Treated
#> Unadjusted 3306.   4194.  
#> Adjusted   1598.88 1884.92
#>  - Time 2
#>            Control Treated
#> Unadjusted 3701.      3799
#> Adjusted   1570.46    1926
# }
```
