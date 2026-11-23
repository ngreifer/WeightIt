# Families for continuous treatments with `method = "glm"` (2026-11-20)

Requested: other conditional distributions for continuous treatments with
`method = "glm"`, like `family` for `method = "bart"`; Poisson, negative
binomial, Gamma, and a Gaussian location-scale model, with bias-reduced versions
where they make sense. Beta omitted at the user's request. CBPS left for later
(no published non-normal continuous CBPS that I know of; see the conversation of
2026-11-20). SuperLearner left alone (its learners cannot take other families).

## Interface

- `family`: `"gaussian"` (default, unchanged), `"gaussian_ls"`, `"poisson"`,
  `"negbin"`, `"Gamma"` (case-insensitive), or a `gaussian()`, `poisson()`, or
  `Gamma()` family object, whose link is used when `link` is not supplied.
  Ignored for binary and multi-category treatments, as with BART.
- `link`: the link for the mean. Gaussian and location-scale: `"identity"`
  (default), `"log"`, `"inverse"`. Poisson and negative binomial: `"log"`
  (default), `"identity"`, `"sqrt"`. Gamma: `"log"` (default), `"inverse"`,
  `"identity"`. No bias reduction (see below).
- `density` and the kernel arguments apply to the Gaussian families only (for
  `"gaussian_ls"`, the residuals standardized by each unit's own SD go through
  `density`); they are ignored with a warning for the others.
- Support is checked: nonnegative integers for Poisson and negative binomial,
  positive values for Gamma.

## Models and estimating equations

All with sampling weights normalized to mean 1.

- Poisson: score for the coefficients.
- Negative binomial: score for the coefficients given theta and the ML score for
  log(theta). Fit by `MASS::glm.nb()` (ML).
- Gamma: score for the coefficients (does not involve the shape) and the ML score
  for log(shape) given the means: `log k + 1 - digamma(k) + log(A/mu) - A/mu`.
  Fit by `glm()`, then the shape by solving its score (not `glm()`'s Pearson
  dispersion, so that the density used for the weights is the ML one).
- Gaussian location-scale: mean linear on the `link` scale, log variance linear
  in the same covariates (a double GLM). Scores: mean
  `X * w * mu.eta * (A - mu) / s2_i`; log variance
  `X * w * ((A - mu)^2 / s2_i - 1) / 2`. Fit by alternating a weighted `glm()`
  for the mean (weights `w / s2_i`) with a Gamma log-link `glm()` on the squared
  residuals (whose score is exactly the normal score for the variance
  coefficients) until the coefficients converge.
- Bias reduction was implemented (a `br.` prefix, with the adjusted scores of
  Firth 1993 / Kosmidis and Firth 2009 solved by WeightIt's own quasi-Fisher
  scoring, the dispersion by ML) and removed on 2026-11-20 at the user's
  request: it was WeightIt's own implementation, which the user could not vet
  and did not want to maintain, and `br.` was not allowed with continuous
  treatments before. It is still not.

## Weights

`w = f(A) / f(A | X)`, with `f(A)` the sampling-weighted average over units of
the conditional density (as for the Gaussian family now and for BART).

- Discrete families: exact, evaluated at each distinct treatment value.
- Gamma: on a grid of 200 points evenly spaced in `log(A)`, with the log density
  interpolated against `log(A)`. A uniform grid on the treatment's own scale was
  poor (worst relative error in the weights 25% at shape 3 with 200 points, and
  far worse for shapes below 1, on 2000 units); the log grid keeps it below 5e-4
  for shapes from 0.3 to 200 (checked against the exact O(n^2) average).
- Gaussian location-scale: exact, at each distinct treatment value
  (`.get_w_from_gps_internal_cont_ls()`, called from
  `.get_w_from_gps_internal_cont()` when the SD varies). The marginal density is
  a sum of one narrow peak per small-SD unit, which no grid captures: with log
  SDs with an SD of 1, a uniform or quantile grid of 200 points, and even a
  1000-point grid at 10,000 units, was off by 13% to over 100% for some units.
  Exact costs O(n^2): 0.4 s at 5000 units, 1.5 s at 10,000. Above 20,000 distinct
  treatment values, a 1000-point grid with log interpolation is used instead (a
  known limitation). The constant-SD Gaussian path is untouched (bit-identical).

## M-estimation

`Mparts` for every family without random effects, `missing = "saem"`, or
`density = "kernel"`: `psi_treat` as above, `wfun` recomputing the weights. For
the location-scale model with the normal density, `dw_dBtreat` is given in closed
form (`.dw_dB_gaussian_ls()`), since a numeric derivative recomputes the O(n^2)
numerator for every parameter and step (`glm_weightit()` at 4000 units: 23 s
numeric, 2.3 s closed form). The others use the numeric derivative.

## Other restrictions

- Random effects: Gaussian (as now), Poisson and Gamma (`lme4::glmer()`),
  negative binomial (`lme4::glmer.nb()`, log link only); not location-scale,
  since lme4 cannot model the variance. For Gamma, the shape is the ML estimate
  given the fitted conditional means: `1 / sigma(fit)^2` from `glmer()` is a
  moment-type estimate that ran low for small shapes (.64 for a true .8 on 2000
  units; ML given the means gave .78, and 3.03 for a true 3).
- `missing = "saem"`: Gaussian only.
- Stabilization fits the numerator with `method = "glm"` and the user's other
  arguments, so with another method (e.g., BART, whose `family` is a *bartisan*
  family) `family` is dropped from the numerator's arguments.

## Verification

- `tests/testthat/test-glm_cont_families.R`: every family, with and without
  sampling weights, passes `expect_M_parts_okay()`
  (the estimates solve `psi_treat`, `wfun` reproduces the weights, and the
  analytic `dw_dBtreat` matches the numeric one); exact marginal densities; ML
  dispersion parameters; errors and warnings; stabilization, `weightitMSM()`,
  and random effects.
- Calibration: `_dev/glm-cont-families-sim.R` (results below).

## Calibration

`_dev/glm-cont-families-sim.R`, results in `_dev/glm-cont-families-sim.rds`
(2026-11-20; 200 replications per cell, 2400 cells, no failures other than the
negative binomial errors noted below). Slope of the MSM, true value .5.

- The estimated models track the oracle in every setting: with strong
  confounding (k = 1), the oracle itself is biased (Gamma +.090, location-scale
  +.13, Poisson +.044, negative binomial +.015) and the correct family matches
  it (+.091, +.13, +.047, +.017). That bias is a property of inverse probability
  weighting with heavy-tailed weights in that design, not of the models.
- With k = .3 and n = 4000, the correct family is unbiased (Poisson +.001,
  negative binomial -.0002, Gamma -.00001, location-scale +.013, the oracle's
  value).
- Against the default Gaussian model: heteroscedastic treatment, bias -.38 (k = 1,
  ESS 14 vs. 408) and -.04 (k = .3) vs. the location-scale model's +.13 (=
  oracle) and +.005/+.013; Gamma treatment, -.019/.013 vs. .008/.000 at k = .3.
  For counts with a mean near 2.7, the Gaussian model is nearly as good.
- Misspecified families: Poisson on negative binomial data is biased (-.08) with
  poor coverage; the location-scale model on Gamma data gives extreme weights
  (ESS 6 to 73) and large bias, hence the note on `?method_glm` that it does not
  accommodate skewness.
- Coverage of the M-estimation intervals at k = .3: Poisson .93/.94, negative
  binomial .94/.93, Gamma .93/.92, location-scale .91/.89 (n = 1000/4000); HC0
  intervals, which ignore estimation of the weights, cover .95 to .99. The
  oracle, whose weights are known and whose HC0 interval is the right one, also
  covers only .91 for Gamma and location-scale at n = 4000, so the shortfall is
  the slow convergence of IPW with heavy-tailed weights rather than an error in
  the estimating functions (which `expect_M_parts_okay()` checks, including the
  analytic weight derivative). With k = 1, every interval undercovers.
- The negative binomial model on Poisson data fails with the intended error in
  33 of 600 fits (no overdispersion to estimate).

## Name of the location-scale family

Checked 2026-11-20 against the installed packages: *bartisan* `gaussian_ls()`
(and `Gamma_ls()`); *mgcv* `gaulss()`, documented as "Gaussian location-scale
model family" (it models the log SD); *VGAM* `uninormal()` with `zero = NULL`;
*gamlss* `NO()` with `sigma.formula`; *brms* and *glmmTMB* `gaussian()` with a
separate formula for the dispersion. `gaussian_ls` matches *bartisan* and is the
same phrase *mgcv* abbreviates; the others have no family name for it.

## Numerator

`_dev/glm-cont-families-numerator.R` (results `.rds`, 2026-11-20; 200
replications, n = 1000) keeps each family's estimated denominator and compares
numerators: the average of the conditional densities (`avg_cond`, used now), a
KDE of the treatment (`kde`), the empirical frequencies for counts (`freq`), and
the same family fit with no covariates (`marginal`, the classic stabilized
numerator). With k = .3:

| | avg_cond | marginal | kde | freq |
|---|---|---|---|---|
| Gamma: SD of estimate / ESS | .045 / 914 | .043 / 917 | .055 / 833 | |
| Location-scale | .076 / 827 | .074 / 829 | .085 / 776 | |
| Poisson | .032 / 900 | .030 / 912 | .038 / 856 | .038 / 856 |
| Negative binomial | .016 / 959 | .016 / 960 | .023 / 864 | .024 / 846 |

The model-free numerators (KDE, frequencies) are the worst on every measure:
SDs 15% to 50% larger, lower ESS, roughly twice the largest weighted
treatment-covariate correlation, and more bias for Gamma (.023 vs. .008). The
intercept-only fit matched or slightly beat the averaged conditional density
with k = .3, but only because those designs barely made it wrong: with weak
confounding, the marginal of the treatment is nearly a member of the same
family, and the two numerators gave almost the same estimate in every data set
(correlations .975 to .997). With k = 1 the paired differences reached .23 in
single data sets (three times the estimator's SD), and the intercept-only ESS
fell to .14 of the averaged numerator's.

Stress test: `_dev/glm-cont-families-numerator-stress.R` (results `.rds`,
2026-11-20; 300 replications, n = 1000). A binary covariate shifts the log mean
(or, for the location-scale design, the log SD) by 1 or 1.8, making the
marginal bimodal or a heavy-tailed scale mixture; the treatment model is
correct. Two outcomes per data set: one linear in A (the linear MSM is correct,
slope .5) and one nonlinear (g(a) = 2 sqrt(a), or a + .2 a^2 for the
location-scale design), fit with a linear MSM, whose target is the projection
under the true marginal of A (by Monte Carlo on 2 million draws). Oracle
weights (true conditional density, sample average of the true conditionals as
the numerator) are the reference.

- The averaged conditional density is within about 2% of the oracle numerator
  (median absolute log ratio .014 to .023) in every design; the intercept-only
  fit is off by .08 to .51, and by 1.49 (a factor of 4.5) for the strongly
  bimodal Poisson design.
- Correct MSM: every numerator is consistent, but the intercept-only fit loses
  precision: RMSE ratio to the averaged numerator 1.0 to 1.2, and 1.94 for the
  strongly bimodal Poisson design (RMSE .083 vs. .043, mean ESS 47 vs. 103); in
  single data sets, its ESS fell to .03 (Gamma) and .11 (Poisson) of the
  averaged numerator's.
- Misspecified MSM: the averaged numerator's bias equals the oracle's in every
  design (e.g., Gamma .052 vs. .051, location-scale -.015 vs. -.020; the oracle
  itself is biased in finite samples with these weights); the intercept-only
  fit shifts it, targeting the projection under the wrong marginal: Gamma
  moderate .081, negative binomial .069/.040 (vs. .050/.028), location-scale
  strong -.107 (vs. -.015), with RMSE ratios up to 1.42.
- It is not uniformly worse: for the moderately bimodal Poisson design, it gave
  a higher ESS (239 vs. 196) and better balance (.12 vs. .19).

Decision (the user's, 2026-11-20): keep the averaged conditional density. The intercept-only fit is
simpler, but it can lose half the precision when the marginal departs from the
family, and with a misspecified MSM it estimates a different estimand.

## Not done

- CBPS for continuous treatments with other families (see the conversation of
  2026-11-20: no published version known; the just-identified balance
  conditions generalize, but standardizing the treatment does not, the
  numerator must be the averaged conditional density, and solver behavior is
  unknown).
- `method = "super"` and `"gbm"` still use the residual-density approach; a
  count-family density could be used with their predicted means later.
- *bartisan*'s grid for the marginal density (`weightit2bart.cont()`) is uniform
  with 50 points, which the Gamma and location-scale checks above suggest can be
  inaccurate for right-skewed or heteroscedastic treatments. Worth checking.
