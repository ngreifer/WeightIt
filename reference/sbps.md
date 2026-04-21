# Subgroup Balancing Propensity Score

Implements the subgroup balancing propensity score (SBPS), which is an
algorithm that attempts to achieve balance in subgroups by sharing
information from the overall sample and subgroups (Dong, Zhang, Zeng, &
Li, 2020; DZZL). Each subgroup can use either weights estimated using
the whole sample, weights estimated using just that subgroup, or a
combination of the two. The optimal combination is chosen as that which
minimizes an imbalance criterion that includes subgroup as well as
overall balance.

## Usage

``` r
sbps(
  obj,
  obj2 = NULL,
  moderator = NULL,
  formula = NULL,
  data = NULL,
  smooth = FALSE,
  full.search
)
```

## Arguments

- obj:

  a `weightit` object containing weights estimated in the overall
  sample.

- obj2:

  a `weightit` object containing weights estimated in the subgroups.
  Typically this has been estimated by including `by` in the call to
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md).
  Either `obj2` or `moderator` must be specified.

- moderator:

  optional; a string containing the name of the variable in `data` for
  which weighting is to be done within subgroups or a one-sided formula
  with the subgrouping variable on the right-hand side. This argument is
  analogous to the `by` argument in
  [`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
  and in fact it is passed on to `by`. Either `obj2` or `moderator` must
  be specified.

- formula:

  an optional formula with the covariates for which balance is to be
  optimized. If not specified, the formula in `obj$call` will be used.

- data:

  an optional data set in the form of a data frame that contains the
  variables in `formula` or `moderator`.

- smooth:

  `logical`; whether the smooth version of the SBPS should be used. This
  is only compatible with `weightit` methods that return a propensity
  score.

- full.search:

  `logical`; when `smooth = FALSE`, whether every combination of
  subgroup and overall weights should be evaluated. If `FALSE`, a
  stochastic search as described in DZZL will be used instead. If
  `TRUE`, all \\2^R\\ combinations will be checked, where \\R\\ is the
  number of subgroups, which can take a long time with many subgroups.
  If unspecified, will default to `TRUE` if \\R \le 8\\ and `FALSE`
  otherwise.

## Value

A `weightit.sbps` object, which inherits from `weightit`. This contains
all the information in `obj` with the weights, propensity scores, call,
and possibly covariates updated from `sbps()`. In addition, the
`prop.subgroup` component contains the values of the coefficients \\C\\
for the subgroups (which are either 0 or 1 for the standard SBPS), and
the `moderator` component contains a data.frame with the moderator.

This object has its own summary method and is compatible with cobalt
functions. The `cluster` argument should be used with cobalt functions
to accurately reflect the performance of the weights in balancing the
subgroups.

## Details

The SBPS relies on two sets of weights: one estimated in the overall
sample and one estimated within each subgroup. The algorithm decides
whether each subgroup should use the weights estimated in the overall
sample or those estimated in the subgroup. There are \\2^R\\
permutations of overall and subgroup weights, where \\R\\ is the number
of subgroups. The optimal permutation is chosen as that which minimizes
a balance criterion as described in DZZL. The balance criterion used
here is, for binary and multi-category treatments, the sum of the
squared standardized mean differences within subgroups and overall,
which are computed using
[`cobalt::col_w_smd()`](https://ngreifer.github.io/cobalt/reference/balance-summary.html),
and for continuous treatments, the sum of the squared correlations
between each covariate and treatment within subgroups and overall, which
are computed using
[`cobalt::col_w_corr()`](https://ngreifer.github.io/cobalt/reference/balance-summary.html).

The smooth version estimates weights that determine the relative
contribution of the overall and subgroup propensity scores to a weighted
average propensity score for each subgroup. If \\P_O\\ are the
propensity scores estimated in the overall sample and \\P_S\\ are the
propensity scores estimated in each subgroup, the smooth SBPS finds
\\R\\ coefficients \\C\\ so that for each subgroup, the ultimate
propensity score is \\C\*P_S + (1-C)\*P_O\\, and weights are computed
from this propensity score. The coefficients are estimated using
[`optim()`](https://rdrr.io/r/stats/optim.html) with
`method = "L-BFGS-B"`. When \\C\\ is estimated to be 1 or 0 for each
subgroup, the smooth SBPS coincides with the standard SBPS.

If `obj2` is not specified and `moderator` is, `sbps()` will attempt to
refit the model specified in `obj` with the `moderator` in the `by`
argument. This relies on the environment in which `obj` was created to
be intact and can take some time if `obj` was hard to fit. It's safer to
estimate `obj` and `obj2` (the latter simply by including the moderator
in the `by` argument) and supply these to `sbps()`.

## References

Dong, J., Zhang, J. L., Zeng, S., & Li, F. (2020). Subgroup balancing
propensity score. *Statistical Methods in Medical Research*, 29(3),
659–676.
[doi:10.1177/0962280219870836](https://doi.org/10.1177/0962280219870836)

## See also

[`weightit()`](https://ngreifer.github.io/WeightIt/reference/weightit.md),
[`summary.weightit()`](https://ngreifer.github.io/WeightIt/reference/summary.weightit.md)

## Examples

``` r
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups within races
(W1 <- weightit(treat ~ age + educ + married +
                nodegree + race + re74, data = lalonde,
                method = "glm", estimand = "ATT"))
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, race, re74

(W2 <- weightit(treat ~ age + educ + married +
                nodegree + race + re74, data = lalonde,
                method = "glm", estimand = "ATT",
                by = "race"))
#> A weightit object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, race, re74
#>  - by: race

S <- sbps(W1, W2)

print(S)
#> A weightit.sbps object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, race, re74
#>  - moderator: race (3 subgroups)

summary(S)
#>                   Summary of weights
#> 
#> - Overall vs. subgroup proportion contribution:
#> 
#>          race = "black" race = "hispan" race = "white"
#> Overall               0               0              0
#> Subgroup              1               1              1
#> 
#>                     Subgroup: race = "black"                    
#> 
#> - Weight ranges:
#> 
#>           Min                                Max
#> treated 1.         ||                       1.  
#> control 0.466 |---------------------------| 3.59
#> 
#> - Units with the 5 most extreme weights by group:
#>                                      
#>              6     5     4     3    1
#>  treated     1     1     1     1    1
#>            559   381   573   303  411
#>  control 2.949 2.949 3.006 3.064 3.59
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.    0.      0.          0
#> control       0.425 0.366   0.093       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted   87.       156
#> Weighted     73.82     156
#> 
#>                     Subgroup: race = "hispan"                    
#> 
#> - Weight ranges:
#> 
#>           Min                                   Max
#> treated 1.                                 || 1.   
#> control 0.021   |------------|                0.505
#> 
#> - Units with the 5 most extreme weights by group:
#>                                       
#>            100    87    44    28     2
#>  treated     1     1     1     1     1
#>            425   456   480   424   412
#>  control 0.412 0.477 0.484 0.497 0.505
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.    0.      0.          0
#> control       0.714 0.579   0.246       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted   61.        11
#> Weighted     40.62      11
#> 
#>                     Subgroup: race = "white"                    
#> 
#> - Weight ranges:
#> 
#>         Min                                   Max
#> treated   1                              || 1.   
#> control   0   |---------|                   0.385
#> 
#> - Units with the 5 most extreme weights by group:
#>                                      
#>             60   42    23    22    10
#>  treated     1    1     1     1     1
#>            592  589   595   546   580
#>  control 0.239 0.27 0.294 0.296 0.385
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.    0.      0.          0
#> control       1.154 0.952   0.619       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  281.        18
#> Weighted    120.78      18

bal.tab(S, cluster = "race")
#> Balance by cluster
#> 
#>  - - - Cluster: black - - - 
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.0016
#> age         Contin.   0.0126
#> educ        Contin.  -0.0332
#> married      Binary   0.0030
#> nodegree     Binary   0.0062
#> re74        Contin.  -0.0826
#> 
#> Effective sample sizes
#>                0   1
#> Unadjusted 87.   156
#> Adjusted   73.82 156
#> 
#>  - - - Cluster: hispan - - - 
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance  -0.2678
#> age         Contin.   0.1196
#> educ        Contin.  -0.0756
#> married      Binary   0.0217
#> nodegree     Binary   0.0018
#> re74        Contin.   0.0114
#> 
#> Effective sample sizes
#>                0  1
#> Unadjusted 61.   11
#> Adjusted   40.62 11
#> 
#>  - - - Cluster: white - - - 
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.0652
#> age         Contin.   0.0191
#> educ        Contin.  -0.0185
#> married      Binary  -0.0015
#> nodegree     Binary   0.0039
#> re74        Contin.  -0.0117
#> 
#> Effective sample sizes
#>                 0  1
#> Unadjusted 281.   18
#> Adjusted   120.78 18
#>  - - - - - - - - - - - - - - 
#> 

#Could also have run
#S <- sbps(W1, moderator = "race")

S_ <- sbps(W1, W2, smooth = TRUE)

print(S_)
#> A weightit.sbps object
#>  - method: "glm" (propensity score weighting with GLM)
#>  - number of obs.: 614
#>  - sampling weights: none
#>  - treatment: 2-category
#>  - estimand: ATT (focal: 1)
#>  - covariates: age, educ, married, nodegree, race, re74
#>  - moderator: race (3 subgroups)

summary(S_)
#>                   Summary of weights
#> 
#> - Overall vs. subgroup proportion contribution:
#> 
#>          race = "black" race = "hispan" race = "white"
#> Overall            0.17            0.25              0
#> Subgroup           0.83            0.75              1
#> 
#>                     Subgroup: race = "black"                    
#> 
#> - Weight ranges:
#> 
#>           Min                                Max
#> treated 1.         ||                       1.  
#> control 0.465 |---------------------------| 3.57
#> 
#> - Units with the 5 most extreme weights by group:
#>                                     
#>              6     5     4    3    1
#>  treated     1     1     1    1    1
#>            559   381   573  303  411
#>  control 2.979 2.979 3.034 3.09 3.57
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.    0.      0.          0
#> control       0.426 0.366   0.094       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted   87.       156
#> Weighted     73.74     156
#> 
#>                     Subgroup: race = "hispan"                    
#> 
#> - Weight ranges:
#> 
#>           Min                                   Max
#> treated 1.                                 || 1.   
#> control 0.025   |-----------|                 0.474
#> 
#> - Units with the 5 most extreme weights by group:
#>                                     
#>            100   87    44   28     2
#>  treated     1    1     1    1     1
#>            269  456   480  424   412
#>  control 0.391 0.45 0.456 0.47 0.474
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.    0.      0.          0
#> control       0.679 0.553   0.225       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted   61.        11
#> Weighted     41.95      11
#> 
#>                     Subgroup: race = "white"                    
#> 
#> - Weight ranges:
#> 
#>         Min                                   Max
#> treated   1                              || 1.   
#> control   0   |---------|                   0.385
#> 
#> - Units with the 5 most extreme weights by group:
#>                                      
#>             60   42    23    22    10
#>  treated     1    1     1     1     1
#>            592  589   595   546   580
#>  control 0.239 0.27 0.294 0.296 0.385
#> 
#> - Weight statistics:
#> 
#>         Coef of Var   MAD Entropy # Zeros
#> treated       0.    0.      0.          0
#> control       1.154 0.952   0.619       0
#> 
#> - Effective Sample Sizes:
#> 
#>            Control Treated
#> Unweighted  281.        18
#> Weighted    120.78      18

bal.tab(S_, cluster = "race")
#> Balance by cluster
#> 
#>  - - - Cluster: black - - - 
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.0019
#> age         Contin.   0.0388
#> educ        Contin.  -0.0305
#> married      Binary   0.0096
#> nodegree     Binary   0.0086
#> re74        Contin.  -0.0561
#> 
#> Effective sample sizes
#>                0   1
#> Unadjusted 87.   156
#> Adjusted   73.74 156
#> 
#>  - - - Cluster: hispan - - - 
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance  -0.1909
#> age         Contin.   0.0167
#> educ        Contin.  -0.0654
#> married      Binary   0.0314
#> nodegree     Binary   0.0084
#> re74        Contin.  -0.0175
#> 
#> Effective sample sizes
#>                0  1
#> Unadjusted 61.   11
#> Adjusted   41.95 11
#> 
#>  - - - Cluster: white - - - 
#> Balance Measures
#>                Type Diff.Adj
#> prop.score Distance   0.0652
#> age         Contin.   0.0191
#> educ        Contin.  -0.0185
#> married      Binary  -0.0015
#> nodegree     Binary   0.0039
#> re74        Contin.  -0.0117
#> 
#> Effective sample sizes
#>                 0  1
#> Unadjusted 281.   18
#> Adjusted   120.78 18
#>  - - - - - - - - - - - - - - 
#> 
```
