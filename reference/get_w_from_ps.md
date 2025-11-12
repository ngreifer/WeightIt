# Compute weights from propensity scores

Given a vector or matrix of propensity scores, outputs a vector of
weights that target the provided estimand.

## Usage

``` r
get_w_from_ps(
  ps,
  treat,
  estimand = "ATE",
  focal = NULL,
  treated = NULL,
  subclass = NULL,
  stabilize = FALSE
)
```

## Arguments

- ps:

  a vector, matrix, or data frame of propensity scores. See Details.

- treat:

  a vector of treatment status for each individual. See Details.

- estimand:

  the desired estimand that the weights should target. Current options
  include `"ATE"` (average treatment effect), `"ATT"` (average treatment
  effect on the treated), `"ATC"` (average treatment effect on the
  control), `"ATO"` (average treatment effect in the overlap), `"ATM"`
  (average treatment effect in the matched sample), and `"ATOS"`
  (average treatment effect in the optimal subset). See Details.

- focal:

  when `estimand` is `"ATT"` or `"ATC"`, which group should be consider
  the (focal) "treated" or "control" group, respectively. If not `NULL`
  and `estimand` is not `"ATT"` or `"ATC"`, `estimand` will
  automatically be set to `"ATT"`.

- treated:

  when treatment is binary, the value of `treat` that is considered the
  "treated" group (i.e., the group for which the propensity scores are
  the probability of being in). If `NULL`, `get_w_from_ps()` will
  attempt to figure it out on its own using some heuristics. This really
  only matters when `treat` has values other than 0 and 1 and when `ps`
  is given as a vector or an unnamed single-column matrix or data frame.

- subclass:

  `numeric`; the number of subclasses to use when computing weights
  using marginal mean weighting through stratification (MMWS; also known
  as fine stratification). If `NULL`, standard inverse probability
  weights (and their extensions) will be computed; if a number greater
  than 1, subclasses will be formed and weights will be computed based
  on subclass membership. `estimand` must be `"ATE"`, `"ATT"`, or
  `"ATC"` if `subclass` is non-`NULL`. See Details.

- stabilize:

  `logical`; whether to compute stabilized weights or not. This simply
  involves multiplying each unit's weight by the proportion of units in
  their treatment group. For saturated outcome models and in balance
  checking, this won't make a difference; otherwise, this can improve
  performance.

## Value

A vector of weights. When `subclass` is not `NULL`, the subclasses are
returned as the `"subclass"` attribute. When `estimand = "ATOS"`, the
chosen value of `alpha` (the smallest propensity score allowed to remain
in the sample) is returned in the `"alpha"` attribute.

## Details

`get_w_from_ps()` applies the formula for computing weights from
propensity scores for the desired estimand. The formula for each
estimand is below, with \\A_i\\ the treatment value for unit \\i\\
taking on values \\\mathcal{A} = (1, \ldots, g)\\, \\p\_{a, i}\\ the
probability of receiving treatment level \\a\\ for unit \\i\\, and \\f\\
is the focal group (the treated group for the ATT and the control group
for the ATC):

\$\$ \begin{aligned} w^{ATE}\_i &= 1 / p\_{A_i, i} \\ w^{ATT}\_i &=
w^{ATE}\_i \times p\_{f, i} \\ w^{ATO}\_i &= w^{ATE}\_i / \sum\_{a \in
\mathcal{A}}{1/p\_{a, i}} \\ w^{ATM}\_i &= w^{ATE}\_i \times \min(p\_{1,
i}, \ldots, p\_{g, i}) \\ w^{ATOS}\_i &= w^{ATE}\_i \times
\mathbb{1}\left(\alpha \< p\_{2, i} \< 1 - \alpha\right) \end{aligned}
\$\$

`get_w_from_ps()` can only be used with binary and multi-category
treatments.

### Supplying the `ps` argument

The `ps` argument can be entered in two ways:

- A numeric matrix with a row for each unit and a (named) column for
  each treatment level, with each cell corresponding to the probability
  of receiving the corresponding treatment level

- A numeric vector with a value for each unit corresponding to the
  probability of being "treated" (only allowed for binary treatments)

When supplied as a vector, `get_w_from_ps()` has to know which value of
`treat` corresponds to the "treated" group. For 0/1 variables, 1 will be
considered treated. For other types of variables, `get_w_from_ps()` will
try to figure it out using heuristics, but it's safer to supply an
argument to `treated`. When `estimand` is `"ATT"` or `"ATC"`, supplying
a value to `focal` is sufficient (for ATT, `focal` is the treated group,
and for ATC, `focal` is the control group).

When supplied as a matrix, the columns must be named with the levels of
the treatment, and it is assumed that each column corresponds to the
probability of being in that treatment group. This is the safest way to
supply `ps` unless `treat` is a 0/1 variable. When `estimand` is `"ATT"`
or `"ATC"`, a value for `focal` must be specified.

### Marginal mean weighting through stratification (MMWS)

When `subclass` is not `NULL`, MMWS weights are computed. The
implementation differs slightly from that described in Hong (2010,
2012). First, subclasses are formed by finding the quantiles of the
propensity scores in the target group (for the ATE, all units; for the
ATT or ATC, just the units in the focal group). Any subclasses lacking
members of a treatment group will be filled in with them from
neighboring subclasses so each subclass will always have at least one
member of each treatment group. A new subclass-propensity score matrix
is formed, where each unit's subclass-propensity score for each
treatment value is computed as the proportion of units with that
treatment value in the unit's subclass. For example, if a subclass had
10 treated units and 90 control units in it, the subclass-propensity
score for being treated would be .1 and the subclass-propensity score
for being control would be .9 for all units in the subclass.

For multi-category treatments, the propensity scores for each treatment
are stratified separately as described in Hong (2012); for binary
treatments, only one set of propensity scores are stratified and the
subclass-propensity scores for the other treatment are computed as the
complement of the propensity scores for the stratified treatment.

After the subclass-propensity scores have been computed, the standard
propensity score weighting formulas are used to compute the unstabilized
MMWS weights. To estimate MMWS weights equivalent to those described in
Hong (2010, 2012), `stabilize` must be set to `TRUE`, but, as with
standard propensity score weights, this is optional. Note that MMWS
weights are also known as fine stratification weights and described by
Desai et al. (2017).

## References

### Binary treatments

- `estimand = "ATO"`

Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates
via propensity score weighting. Journal of the American Statistical
Association, 113(521), 390–400.
[doi:10.1080/01621459.2016.1260466](https://doi.org/10.1080/01621459.2016.1260466)

- `estimand = "ATM"`

Li, L., & Greene, T. (2013). A Weighting Analogue to Pair Matching in
Propensity Score Analysis. The International Journal of Biostatistics,
9(2). [doi:10.1515/ijb-2012-0030](https://doi.org/10.1515/ijb-2012-0030)

- `estimand = "ATOS"`

Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009).
Dealing with limited overlap in estimation of average treatment effects.
Biometrika, 96(1), 187–199.
[doi:10.1093/biomet/asn055](https://doi.org/10.1093/biomet/asn055)

- Other estimands

Austin, P. C. (2011). An Introduction to Propensity Score Methods for
Reducing the Effects of Confounding in Observational Studies.
Multivariate Behavioral Research, 46(3), 399–424.
[doi:10.1080/00273171.2011.568786](https://doi.org/10.1080/00273171.2011.568786)

- Marginal mean weighting through stratification (MMWS)

Hong, G. (2010). Marginal mean weighting through stratification:
Adjustment for selection bias in multilevel data. Journal of Educational
and Behavioral Statistics, 35(5), 499–531.
[doi:10.3102/1076998609359785](https://doi.org/10.3102/1076998609359785)

Desai, R. J., Rothman, K. J., Bateman, B. . T., Hernandez-Diaz, S., &
Huybrechts, K. F. (2017). A Propensity-score-based Fine Stratification
Approach for Confounding Adjustment When Exposure Is Infrequent:
Epidemiology, 28(2), 249–257.
[doi:10.1097/EDE.0000000000000595](https://doi.org/10.1097/EDE.0000000000000595)

### Multi-Category Treatments

- `estimand = "ATO"`

Li, F., & Li, F. (2019). Propensity score weighting for causal inference
with multiple treatments. The Annals of Applied Statistics, 13(4),
2389–2415.
[doi:10.1214/19-AOAS1282](https://doi.org/10.1214/19-AOAS1282)

- `estimand = "ATM"`

Yoshida, K., Hernández-Díaz, S., Solomon, D. H., Jackson, J. W., Gagne,
J. J., Glynn, R. J., & Franklin, J. M. (2017). Matching weights to
simultaneously compare three treatment groups: Comparison to three-way
matching. Epidemiology (Cambridge, Mass.), 28(3), 387–395.
[doi:10.1097/EDE.0000000000000627](https://doi.org/10.1097/EDE.0000000000000627)

- Other estimands

McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E.,
Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score
Estimation for Multiple Treatments Using Generalized Boosted Models.
Statistics in Medicine, 32(19), 3388–3414.
[doi:10.1002/sim.5753](https://doi.org/10.1002/sim.5753)

- Marginal mean weighting through stratification (MMWS)

Hong, G. (2012). Marginal mean weighting through stratification: A
generalized method for evaluating multivalued and multiple treatments
with nonexperimental data. *Psychological Methods*, 17(1), 44–60.
[doi:10.1037/a0024918](https://doi.org/10.1037/a0024918)

## See also

[`method_glm`](https://ngreifer.github.io/WeightIt/reference/method_glm.md)

## Examples

``` r
library("cobalt")
data("lalonde", package = "cobalt")

ps.fit <- glm(treat ~ age + educ + race + married +
                nodegree + re74 + re75, data = lalonde,
              family = binomial)
ps <- ps.fit$fitted.values

w1 <- get_w_from_ps(ps, treat = lalonde$treat,
                    estimand = "ATT")

treatAB <- factor(ifelse(lalonde$treat == 1, "A", "B"))
w2 <- get_w_from_ps(ps, treat = treatAB,
                    estimand = "ATT", focal = "A")
all.equal(w1, w2)
#> [1] TRUE
w3 <- get_w_from_ps(ps, treat = treatAB,
                    estimand = "ATT", treated = "A")
all.equal(w1, w3)
#> [1] TRUE

# Using MMWS
w4 <- get_w_from_ps(ps, treat = lalonde$treat,
                    estimand = "ATE", subclass = 20,
                    stabilize = TRUE)

# A multi-category example using predicted probabilities
# from multinomial logistic regression
T3 <- factor(sample(c("A", "B", "C"), nrow(lalonde),
                    replace = TRUE))

multi.fit <- multinom_weightit(
  T3 ~ age + educ + race + married +
    nodegree + re74 + re75, data = lalonde,
  vcov = "none"
)

ps.multi <- fitted(multi.fit)
head(ps.multi)
#>           A         B         C
#> 1 0.3976058 0.2609169 0.3414773
#> 2 0.2793092 0.2324854 0.4882054
#> 3 0.3728859 0.2908927 0.3362214
#> 4 0.3711736 0.2665370 0.3622894
#> 5 0.3310437 0.2942148 0.3747415
#> 6 0.3485804 0.3109687 0.3404509

w5 <- get_w_from_ps(ps.multi, treat = T3,
                    estimand = "ATE")
```
