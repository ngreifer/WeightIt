# Plan: Weights for Modified Treatment Policies (2026-11-18)

This plan covers how *WeightIt* could estimate weights for shift interventions and other modified treatment policies (MTPs), how balance would be assessed, and how effects would be estimated, including with *adrftools*' shift curves. No package code has been changed. The numerical claims marked "checked" come from `_dev/mtp-checks.R`.

**Recommendation.** Start with additive shifts, $q(a, x) = a + \delta(x)$, where $\delta$ is a constant or a per-unit vector, requested with a marker on the left side of the formula, `weightit(.shift(A, 0.5) ~ X1 + X2, ...)`, the way `.cens()` requests censoring weights. Every continuous method can support this with a small, local change, because an MTP changes only *what* the weights target, not how they are estimated:

- Density-ratio methods (`glm`, `gbm`, `super`, `bart`) change the numerator of the weight from the marginal density at $A$ to the conditional density at $A - \delta$.
- Balancing methods (`ebal`, `optweight`, `npcbps`, `energy`, `cbps`) change their targets from the moments of an independent $(A, X)$ to the moments of the shifted sample $(A + \delta, X)$.

Balance should be assessed against the shifted sample, which needs a new *cobalt* class modeled on `bal.tab.cens`; stacking the two samples works with *cobalt* as it is in the meantime (checked). Effects come from weighted g-computation with `glm_weightit()`, which is doubly robust for this estimand.

MTP weights are fit for one shift value at a time. Correctly estimated ADRF weights do not substitute for them, because the shift estimand keeps each unit's covariates with its own treatment (checked). So *WeightIt*'s workflow is a single shift and its contrast. Shift curves and the ADE stay in *adrftools* as outcome-model g-computation, which needs nothing new from *WeightIt*. Weights for the average derivative effect (ADE) are possible, including with energy balancing, but they are signed, and with linear balance constraints they reproduce regression exactly. If they are offered, it should be outside `weightit()`.

## Sources

- Haneuse and Rotnitzky (2013), `haneuseEstimationEffectInterventions2013`. Equations were read from the rendered pages, not the text layer.
- Jiang and Huling (2023), `jiangEnhancingModifiedTreatment2023`, arXiv 2310.11620. Both versions were read from the LaTeX source. Version 2 is retitled "Modified treatment policy effect estimation with weighted energy distance" and replaces the bootstrap threshold for choosing the shift size with a split-sample permutation test.
- *adrftools* at commit 53bc809: `R/shift_curve.R`, `R/ade.R`, `R/adrf.R`, `R/aux_functions.R`, `_dev/plan-shift-effects.md`, `_dev/plan-doubly-robust.md`, and `_dev/variance-theory.md`.

Items in the Zotero library that look relevant but were **not read** for this plan: Hines, Diaz-Ordaz, and Vansteelandt (2026), `hinesParameterizingEffectContinuous2026`, on ADEs; Williams, Hines, and Rudolph (2025), `williamsRieszRepresentersRest2025`; Pal and Huling (2026), `palHeterogeneousEffectsContinuous2026`, on conditional MTPs; and Sengupta and Imai (2026), `senguptaBalancingCovariatesContinuous2026`. Two papers that the design leans on are **not in the library** and should be added before they are cited: Díaz et al.'s *lmtp* paper (*JASA*, 2023), the source of the classification density-ratio weights that Jiang and Huling use as a comparator, and Díaz Muñoz and van der Laan (2012) on stochastic interventions.

## Notation and Estimand

A policy $q(a, x)$ maps each unit's observed treatment and covariates to a new treatment value, $Q = q(A, X)$. The estimand is $\mu^q = E[Y(Q)]$, usually contrasted with the observed mean $E[Y]$. Under consistency, positivity in the form "$(q(a, x), x)$ is in the support whenever $(a, x)$ is", and exchangeability between the units that received $a$ and those that received $q(a, x)$ (both papers' assumptions; Haneuse and Rotnitzky call these C0 to C2), the estimand is identified by the g-formula,

$$
\mu^q = E\big[m(q(A, X), X)\big], \qquad m(a, x) = E[Y \mid A = a, X = x],
$$

and by weighting, $\mu^q = E[\lambda(A, X) Y] / E[\lambda(A, X)]$. For a policy that is piecewise smooth and invertible in $a$, with pieces $j$ and inverses $h_j$,

$$
\lambda(a, x) = \sum_j I\{h_j(a, x) \in \mathcal{I}_j(x)\}\, \frac{f(h_j(a, x) \mid x)}{f(a \mid x)}\, \big|h_j'(a, x)\big|.
$$

Both papers print $h_j'$ without the absolute value. Every example in either paper is increasing in $a$, so it never matters there, but an implementation should take the absolute value. For an additive shift the weight is $f(a - \delta \mid x) / f(a \mid x)$.

Haneuse and Rotnitzky normalize the weights (Hájek). They also warn (p. 5271) that, as with the ATT, standard errors that treat the weights as known are not guaranteed to be conservative for $\mu^q$, so the variance has to account for estimating the weights.

A policy that sends positive mass to a single point, such as Taubman's "raise everyone below 30 to 30", has no density ratio. The g-formula still covers it, but weighting does not, so *WeightIt* should not accept it.

## The Distribution the Weights Target

Jiang and Huling's main point is that MTP weights should make the weighted sample $\{(A_i, X_i), w_i\}$ look like the *unweighted shifted sample* $\{(Q_i, X_i)\}$. Their supplement shows that the target distribution in the identification formula is the joint distribution of $(Q, X)$. The covariate distribution is unchanged; only the treatment moves.

This gives one frame for every continuous-treatment weight in *WeightIt*. Weights for the ADRF target the product of the marginals, $f(a) f(x)$. Weights for an MTP target the law of $(q(A, X), X)$. In moment form, the condition is

$$
\sum_i w_i\, g(A_i, X_i) = \sum_i g(q(A_i, X_i), X_i)
$$

for basis functions $g$. If $m(a, x)$ is linear in those functions, exact balance on them removes all confounding bias. The ADRF condition is the same equation with $(A, X)$ on the right side drawn independently.

For an additive shift with constant $\delta$, the consequences are specific:

- The covariate means and distributions should match those of the unweighted sample.
- The weighted mean of $A$ should be $\bar{A} + \delta$, and the weighted distribution of $A$ should be the observed one translated by $\delta$.
- The weighted covariance of $A$ with each covariate should equal the *unweighted* covariance, since $\mathrm{Cov}(A + \delta, X) = \mathrm{Cov}(A, X)$.

The last point is the one most likely to mislead. The current continuous-treatment balance statistics in *cobalt*, the treatment–covariate correlations, would report a correctly weighted MTP sample as badly imbalanced, because the target correlations are the observed ones, not zero.

Under a normal linear model for $A \mid X$ with residual SD $\sigma$, the weight is an exponential tilt of the residual $r = A - \mu(X)$:

$$
\lambda = \exp\!\left(\frac{\delta r}{\sigma^2} - \frac{\delta^2}{2\sigma^2}\right).
$$

So the log weight is linear in $A$ and $X$ (checked: regressing $\log \lambda$ on $A$, $X_1$, $X_2$ returns a coefficient on $A$ equal to $\delta/\sigma^2$). Two consequences:

- These weights lie in the same exponential family as entropy balancing weights on the means of $A$ and $X$. The relationship between the normal-GPS weights and `ebal` is therefore the same as that between logistic regression and `ebal` for the ATT.
- $E[\lambda^2] = \exp(\delta^2/\sigma^2)$, so the effective sample size is about $n \exp(-\delta^2/\sigma^2)$ (checked: 0.790 vs. 0.790 at $\delta/\sigma = 0.49$). Shifts of 0.5, 1, and 2 conditional SDs keep about 78%, 37%, and 2% of the sample.

What matters is the shift in units of the **conditional** SD of the treatment, not its marginal SD. *adrftools* currently expresses shifts in marginal SD units, so this is worth surfacing there too.

## ADRF Weights Used for a Shift

Correctly estimated ADRF weights do **not** in general give valid shift effects through an unadjusted weighted outcome model. The weighted sample resembles a trial in which $A$ is assigned independently of $X$, so the unadjusted weighted regression estimates the ADRF, $\theta(a) = E[m(a, X)]$, which averages over the *marginal* distribution of $X$. The shift estimand averages $m(q(a, x), x)$ over the *joint* distribution of $(A, X)$, keeping each unit's covariates with its own treatment. Haneuse and Rotnitzky call their estimand "ATT-like" for this reason. The relationship is the one between the ATE and the ATT: ATE weights followed by an unadjusted comparison estimate the ATE, which differs from the ATT when the effect varies with covariates that predict treatment.

Checked with $n = 200{,}000$, $A \mid X \sim N(0.8X, 1)$, $\delta = 0.5$, true ADRF weights, and a quartic weighted regression of $Y$ on $A$ alone:

| Outcome model $m(a, x)$ | True contrast | ADRF weights | MTP weights |
|---|---|---|---|
| $a + 0.5a^2 + x$ (no effect modification) | 0.625 | 0.637 | 0.632 |
| $a + ax + x$ (modification linear in $a$) | 0.500 | 0.524 | 0.506 |
| $a + 0.5a^2 x + x$ (modification nonlinear in $a$) | 0.903 | 0.562 | 0.912 |

In the third case the gap is $\delta\,\mathrm{Cov}(A, X) = 0.4$, as the algebra predicts. In the second case the ADRF-weighted contrast is right only if computed as $\hat\theta(A + \delta) - \hat\theta(A)$. Computed against the observed mean of $Y$ it gives $-0.25$, because the *level* $E[\theta(A + \delta)]$ is off by $\mathrm{Cov}(A, X)$ and the observed mean is not.

The ADRF-weighted approach therefore works when no covariate associated with the treatment modifies the dose–response, or when that modification is linear in $a$ and the contrast is computed against $\hat\theta(A)$. A covariate-adjusted outcome model fit with ADRF weights is a different matter: predicting each unit at $q(A_i, X_i)$ with its own $X_i$ is consistent if the outcome model is right. It is not doubly robust, though, since the weights do not target the shifted distribution. ADRF weights also need positivity over the whole range of $A$ at every $x$, which is the assumption MTPs are designed to avoid, so they are usually much more variable than the MTP weights for a modest shift.

## Weight Estimation

### Density-Ratio Methods

All of `glm`, `gbm`, and `super` compute the weight in `.get_w_from_gps_internal_cont()` (`R/functions_for_processing.R:1665`). It evaluates the density of the standardized residual at the observed treatment for the denominator, and marginalizes it over a grid for the numerator. For a shift, the numerator becomes the same density evaluated at the residual of the inverse-shifted treatment, $(h(A_i, X_i) - \mu_i)/\sigma$, times the Jacobian. Nothing is marginalized, so the grid and spline steps are skipped. This is a new argument to that function, not a new function.

M-estimation comes almost free for `glm`. The `psi_treat` in `weightit2glm.cont()` (`R/weightit2glm.R:1336`) is unchanged, since the treatment model is the same, and `wfun` (`:1346`) only needs to close over the shift. The same holds wherever a continuous method already returns `Mparts`.

Points to check:

- **Kernel densities.** With `density = "kernel"`, `.get_make_dens_fun()` builds the density on `[min(grid), max(grid)]`, and `approxfun()` returns `NA` outside it, which `.check_density()` turns into an error. The grid must cover the shifted residuals.
- **`bart`.** The conditional density comes from the *bartisan* family (`dpm()` by default), not from `.get_w_from_gps_internal_cont()`. Whether *bartisan* can evaluate the fitted conditional density at an arbitrary treatment value for each unit is an open question. Without that, `bart` is out of the first phase.
- **`gbm` tuning.** `gbm` chooses the number of trees by a balance criterion from `cobalt::available.stats("continuous")` (`p.mean` by default; `R/weightit2gbm.R`, `criterion`). Those statistics measure independence, which is the wrong target here. Censoring already solved the same problem by switching to *cobalt*'s "target" statistics, so the MTP version needs the *cobalt* support described under *Balance Assessment* before `gbm` can be tuned correctly. `super` with a balance-based `SL.method` has the same dependency.

### Balancing Methods

Balancing methods need only the forward map $q$, never its inverse or Jacobian, because the target is computed from the vector $Q_i$. That makes them the natural route for complicated policies.

- **`ebal`** (`weightit2ebal.cont()`, `R/weightit2ebal.R:651`) already builds the constraint matrix from columns $A$, $X$, and $A \cdot X$, standardized so that the targets are the means and zero cross-moments. For an MTP, the targets become the sample means of $Q$, $X$, and $Q \cdot X$, standardized with $A$'s centering and scale. `d.moments` carries over (higher moments of $Q$). Because `ebal` returns `Mparts` when `tols` is 0, the M-estimation part needs the constant targets moved into its estimating equations.
- **`optweight`**: this is a one-sample calibration problem, so it goes through *optweight*'s target interface, `optweight.svy()`, not its continuous-treatment interface. The "covariates" are the basis columns $A$, $X$, $A \cdot X$ (and $A^2$ with `moments = 2`), and `targets` are the sample means of the same columns computed with $Q$ in place of $A$. Nothing in *optweight* has to change. A `shift` argument there would only be a convenience that builds those columns and targets. That is reasonable if *optweight* is meant to be used for MTPs on its own, but *WeightIt* does not need it. The same view applies to `ebal`, which is solving the same target problem with a different objective.
- **`npcbps`**: empirical likelihood with the same moment conditions. This is lower priority.
- **`energy`**: Jiang and Huling's weighted energy distance is the two-sample energy distance between the weighted observed points $(X_i, A_i)$ and the unit-weighted shifted points $(X_i, Q_i)$, on standardized $(X, A)$ jointly. $Q$ must be scaled with $A$'s mean and SD, not standardized on its own. This is structurally the binary ATT problem in `weightit2energy()` (weighted controls to a fixed treated sample), so that code can be reused rather than the distance-covariance objective of `weightit2energy.cont()`. Two notes from the paper:
  - The printed quadratic form has the wrong sign: it should be $-w^\top D w + 2 w^\top b$, not $+w^\top D w$.
  - The paper recommends adding the ridge penalty $\lambda n^{-2} \sum_i w_i^2$; *WeightIt*'s existing `lambda` argument plays this role.
- **`cbps`**: covariate balancing GPS with MTP balance conditions in place of the independence ones. This needs its own derivation, so it belongs in a later phase.

### Stacked Classification

There is a third route that needs no new estimation code. Stack the observed sample (label 0) with the shifted sample (label 1), fit any *binary* method with `estimand = "ATT"`, and keep the weights of the label-0 rows. The ATT weights $p/(1 - p)$ for the label-0 rows estimate the density ratio of $(Q, X)$ to $(A, X)$, which is $\lambda$. This is the classification trick Jiang and Huling compare against (from Díaz et al.).

Checked with $n = 4000$ and a normal linear treatment model:

- Logistic regression on $A + X_1 + X_2$ in the stacked data gives weights perfectly correlated with the normal-GPS density ratio (the log density ratio is linear in $A$ and $X$, so logistic regression is correctly specified).
- The Hájek estimates are 5.434 for both, against a truth of 5.418, and 4.385 unweighted.
- Entropy balancing in the stacked data on $A$, $X$, $A^2$, and $A \cdot X$ gives 5.412.
- Energy balancing in the stacked data ($n = 1000$) gives 5.852 against a truth of 5.784; that is one draw, not a comparison.

Stacking works for any $q$ with only its forward map, and with every binary method, including `gbm`, `super`, and `bart`. Its drawbacks are the doubled sample size (an energy distance matrix of $(2n)^2$) and that the classifier must learn that $X$ alone carries no information. Its M-estimation would need the two rows of each unit summed. I would use it as the reference implementation in tests, and possibly as the internal fallback for density-ratio methods when $q$ has no supplied inverse, rather than as the main route.

### Method Support

| Method | Route | Needs inverse of $q$ | M-estimation |
|---|---|---|---|
| `glm` | density ratio | yes | yes (same `psi_treat`, new `wfun`) |
| `gbm`, `super` | density ratio | yes | no (as now) |
| `bart` | density ratio | yes | no; depends on *bartisan* |
| `ebal` | moment targets | no | yes, with targets in the equations |
| `optweight` | moment targets via `optweight.svy()` | no | no (as now) |
| `npcbps` | moment targets | no | no |
| `energy` | two-sample energy distance | no | no |
| `cbps` | GPS with MTP conditions | yes | later |

### Policies Beyond Additive Shifts

| Policy | Inverse | Notes |
|---|---|---|
| Additive, $a + \delta$ | $a - \delta$, Jacobian 1 | Phase 1. |
| Covariate-dependent additive, $a + \delta(x)$ | $a - \delta(x)$, Jacobian 1 | Phase 1, as a per-unit vector of $\delta$. |
| Multiplicative, $a(1 + \delta)$ | $a/(1 + \delta)$, Jacobian $1/(1 + \delta)$ | Matches *adrftools*' `"multiplicative"`; needs $A \ge 0$. |
| Thresholded shift (Haneuse and Rotnitzky's eq. 3) | two pieces | A unit can stand in for two pre-images, giving weights like $1 + f(a + \delta \mid x)/f(a \mid x)$. Needs a lower bound $b(x)$ from the user. |
| User function | user-supplied inverse, or none | Balancing methods need only the forward map. Density methods need the inverse and its derivative, or fall back to stacking. |

## Interface

```r
W <- weightit(.shift(A, 0.5) ~ X1 + X2 + X3, data = d,
              method = "glm")

W2 <- weightit(.shift(A, d$delta_i) ~ X1 + X2 + X3, data = d,
               method = "ebal", moments = 2)
```

`.shift(x, by, type = "additive")` would be defined in *cobalt* and re-exported by *WeightIt*, exactly as `.cens()` is, so that `bal.tab()` recognizes the marker too. The sign convention should match *adrftools*' `shift_fun`: additive is $a + \delta$, multiplicative is $a(1 + \delta)$. Haneuse and Rotnitzky write $a - \delta$, so the documentation should be explicit about the direction. Later arguments would be a `fun` (and `inverse`) for user policies and `lower`/`upper` for thresholded shifts.

A marker is preferable to a new `shift` argument to `weightit()` for three reasons:

- *cobalt* needs the information to assess balance, and the marker carries it into `bal.tab(formula)` with no extra arguments.
- In `weightitMSM()`, each entry of `formula.list` can carry its own shift. That is the longitudinal MTP setting of *lmtp*, and a later extension.
- It follows the precedent already set by `.cens()`.

The `weightit` object would add a `shift` component holding the policy specification and the vector $Q$. *cobalt* needs $Q$ to assess balance. `print()` names the policy; `summary()` adds the share of $Q$ outside the observed range of $A$, alongside the usual ESS and weight range. `stabilize` does not apply to MTPs and should be rejected: the weights already have a numerator. `estimand` and `focal` do not apply. `by` works unchanged. Sampling weights enter as in Jiang and Huling's setup: they weight the shifted target sample and multiply the estimated weights.

## Balance Assessment

What to compare, between the weighted observed sample and the shifted sample:

- each covariate's mean and distribution (SMD and KS against the unweighted sample, which *cobalt* already computes as its "target" statistics);
- the treatment's weighted mean, SD, and KS statistic against those of $Q$;
- the cross-moments $A \cdot X_j$ against $Q \cdot X_j$, standardized by the target's SD; the equivalent in correlation form is the weighted $\mathrm{corr}(A, X_j)$ against $\mathrm{corr}(Q, X_j)$, but the cross-moments are what the balancing methods target;
- a joint summary, Jiang and Huling's energy distance on $(X, A)$.

**Today, with no changes.** Stacking reproduces all of this in *cobalt* as it is (checked):

```r
Q  <- d$A + delta
st <- rbind(transform(d, S = 0), transform(d, A = Q, S = 1))

bal.tab(S ~ A + X1 + X2 + A:X1 + A:X2, data = st,
        weights = c(W$weights, rep(1, nrow(d))),
        estimand = "ATT", s.d.denom = "treated",
        stats = c("m", "ks"), un = TRUE)

bal.compute(st[c("A", "X1", "X2")], treat = st$S,
            weights = c(W$weights, rep(1, nrow(d))),
            stat = "energy.dist")
```

The "Control" group is the weighted observed sample and the "Treated" group the shifted target, so `Diff.Adj` is what the weights should drive to zero. In the check, the unweighted SMD for $A$ was 0.376 and both weightings brought it under 0.01. The covariate rows were exactly 0 before weighting, as they should be, and the $A \cdot X_1$ row was small before weighting because the covariance does not change under the shift. The energy distance fell from 0.056 unweighted to 0.0007 (density ratio) and 0.0003 (energy).

**In *cobalt*.** A `bal.tab.shift` class, modeled on `bal.tab.cens` (`cobalt/R/class-bal.tab.cens.R`), would do the stacking internally when the treatment carries a `.shift()` marker or the input is a `weightit` object with a `shift` component. Rows for the covariates, the treatment, and the treatment–covariate products, plus the energy distance in `bal.compute()`, and `available.stats("shift")` for `gbm` and `super` to tune on. `bal.plot()` would overlay the weighted density of $A$ on the density of $Q$, and the weighted and target scatterplots of $A$ against a covariate. `love.plot()` should work as is on the new class.

**Choosing the shift size.** Jiang and Huling's diagnostic plots the minimized energy distance against the shift size, with a threshold from a split-sample permutation test (version 2; version 1 used a bootstrap). The feasible shifts are those below the crossing point. A plain sweep over shift values, reporting ESS and the energy distance, is enough for a first version and is a few lines of user code. The normal-model formula $\mathrm{ESS} \approx n \exp(-\delta^2/\sigma^2)$ gives a quick prior expectation. The permutation test could live in *cobalt* later. The ESS curve over $\delta$ is something neither paper has, and it is the natural companion to *adrftools*' shift curves.

## Effect Estimation

**Hájek IPW.** `lm_weightit(Y ~ 1, data = d, weightit = W)` estimates $\mu^q$, with an M-estimation standard error that accounts for estimating the weights whenever the method returns `Mparts`. The contrast with $E[Y]$ is not available this way, because the observed mean is unweighted and `lm_weightit()` fits one weighted model. The influence function of the contrast is simple (Haneuse and Rotnitzky, p. 5269), so a contrast belongs either in a small *WeightIt* helper or in *adrftools*. I lean toward *adrftools*.

**Weighted g-computation.** Fit the outcome model with the MTP weights, then average its predictions at the shifted treatment:

```r
fit <- glm_weightit(Y ~ splines::ns(A, 3) * (X1 + X2), data = d,
                    weightit = W)

marginaleffects::avg_predictions(fit, newdata = transform(d, A = A + delta))
```

This is doubly robust for $\mu^q$. A canonical-link model with an intercept sets $\sum_i w_i (Y_i - \hat m_i) = 0$, so the augmented estimator, $n^{-1} \sum_i \hat m(Q_i, X_i) + \sum_i w_i (Y_i - \hat m_i) / \sum_i w_i$ (Jiang and Huling's eq. 15, Hájek form), collapses to the plain g-computation average. Checked: with a misspecified linear outcome model, unweighted g-computation gave 5.291 and MTP-weighted g-computation 5.426, against a truth of 5.418.

One subtlety for the contrast: `avg_comparisons()` at $A$ and $A + \delta$ from this fit is doubly robust only for its first term. The average prediction at the observed $A$ estimates $E[Y]$ only if the outcome model is right, because the weighted fit does not force the unweighted residuals to sum to zero. The doubly robust contrast uses the sample mean of $Y$ for $E[Y]$.

**Augmented estimator with one outcome model.** For a grid of shifts, refitting the outcome model with each shift's weights is wasteful. The augmented form uses one outcome model, fit however the user likes, and one weight vector per shift value:

$$
\hat\mu^{q_\delta} = \frac{1}{n} \sum_i \hat m(q_\delta(A_i, X_i), X_i) + \frac{\sum_i \hat\lambda_{\delta, i} (Y_i - \hat m(A_i, X_i))}{\sum_i \hat\lambda_{\delta, i}}.
$$

This is the doubly robust shift curve in *adrftools*' `_dev/plan-doubly-robust.md`. With a density model, $\hat\lambda_\delta$ for every $\delta$ comes from one fit. This plan leaves it to *adrftools* to take that treatment model directly, rather than having *WeightIt* supply a weight per shift value (see *Integration with adrftools*).

## Integration with *adrftools*

**Current state.** `shift_curve()` is g-computation from one outcome model, averaged with one vector of target weights `wts` for every shift value. It never uses the observed $Y$. *WeightIt* enters only through `estfun()` and `bread()` of `glm_weightit` fits and through re-calling `weightit()` in the bootstrap. There is no slot for weights that vary with the shift.

**What the weights add to a curve.** A shift curve from `shift_curve()` is valid whenever the outcome model is right, whatever weights it was fit with. ADRF weights do not make it valid without that (see *ADRF Weights Used for a Shift*), and MTP weights apply to one shift value each. So the weights add robustness to one point of the curve at a time, never to the curve as a whole.

That makes the division of labor simple:

- **A single shift** is a *WeightIt* workflow, with no *adrftools* needed: `weightit(.shift(A, delta) ~ X)`, then balance against the shifted sample, then `glm_weightit()` and the contrast at $A + \delta$. The estimate is doubly robust and the standard error accounts for the weights.
- **A curve over shift values, or the ADE** is an *adrftools* workflow built on the outcome model, as now.
- If doubly robust curves are wanted later, the cheapest route is for *adrftools* to take a *treatment model* directly (e.g., `lm(A ~ X)`), from which $\lambda_\delta$ for every $\delta$ is a closed-form density ratio. That needs nothing from *WeightIt*. A weight matrix with one column per shift value is dropped from this plan. Haneuse and Rotnitzky's marginal structural model over $\delta$ (their eqs. 17 to 21) does fit weights at each of ten or more shift values, so there is a precedent if the question comes back.

*WeightIt* still owes *adrftools* two small things. Once a call carries `.shift()`, *adrftools*' bootstrap path, which re-calls `weightit()`, works unchanged. Its default grid in marginal SD units should also be reported in conditional SD units, given the ESS formula above.

## The Average Derivative Effect

The ADE, $\tau = E[\partial_a m(A, X)]$, is the derivative of $\mu^{q_\delta}$ at $\delta = 0$. A weighting estimator is $\hat\tau = n^{-1} \sum_i \omega_i Y_i$, where the population weight is $\omega = -\partial_a \log f(A \mid X)$ (integrate by parts, assuming the density vanishes at the ends of the support). Under a normal linear model this is $(A - \mu(X))/\sigma^2$. Checked: $\partial \lambda_\delta / \partial \delta$ at 0 matches it.

### Balancing Constraints

For an outcome model in the span of basis functions $g_k(a, x)$, exact balance means

$$
\frac{1}{n} \sum_i \omega_i\, g_k(A_i, X_i) = \frac{1}{n} \sum_i \partial_a g_k(A_i, X_i) \quad \text{for each } k,
$$

so that $n^{-1} \sum_i \omega_i \hat m(A_i, X_i)$ equals the average derivative of $\hat m$ for every $m$ in the span. For the basis $(1, X, A, A X, A^2)$ the constraints are:

| Basis function | Constraint on $n^{-1} \sum_i \omega_i g(A_i, X_i)$ | Meaning |
|---|---|---|
| $1$ | $= 0$ | positive and negative weights carry equal mass |
| $X_j$ | $= 0$ | the positive and negative groups have the same covariate means |
| $A$ | $= 1$ | the comparison is per unit of treatment |
| $A X_j$ | $= \bar X_j$ | effect modification by $X_j$ is averaged over the sample's $X_j$ |
| $A^2$ | $= 2 \bar A$ | curvature in $A$ is averaged over the sample's $A$ |

The weights are necessarily signed: $\sum_i \omega_i = 0$ together with $\sum_i \omega_i A_i = n$ cannot hold with nonnegative weights. Write $M$ for the total positive mass and $p_i = |\omega_i| / M$ for the normalized weights within each sign group. The estimator is then a Wald-type ratio:

$$
\hat\tau = \frac{E_{p^+}[Y] - E_{p^-}[Y]}{E_{p^+}[A] - E_{p^-}[A]},
$$

the difference in mean outcome between two covariate-balanced groups, divided by their difference in mean treatment. That reading gives a natural balance check: compare the positive and negative groups on the covariates with *cobalt*, as for a binary treatment, and report the difference in mean treatment.

### Example

*adrftools*' `nhanes3lead`, with $A$ = `logBLL`, $Y$ = `Math`, and $X$ = age, sex, race, PIR, food security, smoking in the home, smoking in pregnancy, and NICU stay ($n = 2521$). The basis is $(1, X, A, A X, A^2)$. Three weight vectors:

- `gps`: $(A - \hat\mu(X))/\hat\sigma^2$ from `lm(logBLL ~ X)`;
- `balancing`: minimum $\sum_i \omega_i^2$ subject to the constraints, which has a closed form, $\omega = G (G^\top G)^{-1} \sum_i \partial_a g(A_i, X_i)$;
- `calibrated`: the `gps` weights moved as little as possible to satisfy the constraints.

Selected constraints:

| Basis function | Target | `gps` | `balancing` | `calibrated` |
|---|---|---|---|---|
| $1$ | 0 | 0 | 0 | 0 |
| Age | 0 | 0 | 0 | 0 |
| PIR | 0 | 0 | 0 | 0 |
| $A$ | 1 | 1 | 1 | 1 |
| $A \times$ Age | 9.020 | 8.969 | 9.020 | 9.020 |
| $A \times$ PIR | 1.628 | 1.442 | 1.628 | 1.628 |
| $A^2$ | 1.941 | 2.030 | 1.941 | 1.941 |

The `gps` weights balance the means of $X$ and fix the scale through $A$ exactly, because they are OLS residuals, but they miss the interaction and curvature constraints. Their estimate, $-0.636$, is the OLS coefficient on $A$ in `Math ~ logBLL + X`, by the Frisch–Waugh theorem. The `balancing` weights give $-0.618$, which equals the regression ADE from `lm(Math ~ G)`.

As a two-group comparison:

- The positive group (52% of units) has a weighted mean `logBLL` of 1.746 and the negative group 0.135, a difference of 1.611.
- The difference in mean `Math` is $-0.995$, so the ratio is $-0.618$.
- Every covariate mean differs by exactly 0 between the groups.
- The effective sample sizes are 815 and 775.

Two equalities in the example hold in general and are worth knowing before building anything:

- **Minimum-norm balancing weights reproduce regression.** With $\omega = G (G^\top G)^{-1} t$, $n^{-1} \omega^\top Y = n^{-1} t^\top \hat\beta_{\text{OLS}}$. This is the ADE analog of the known result that unconstrained-sign balancing weights reproduce regression imputation.
- **Calibrating linear-GPS weights changes nothing beyond that.** The `calibrated` weights and the augmented estimator, $\overline{\partial_a \hat m} + n^{-1} \sum_i \omega^{\text{gps}}_i (Y_i - \hat m_i)$, both returned $-0.618$. The normal-linear GPS score is linear in $A$ and $X$, so it already lies in the span of the basis. Double robustness adds something only when the treatment model has terms the outcome basis lacks, such as a richer mean model or a covariate-dependent variance.

So ADE weights are only worth offering in a form that regression does not already give: a kernel-balancing version, or GPS-score weights from a flexible treatment model.

### Energy and Kernel Balancing

Energy balancing works for the ADE, not only for substantial shifts. Kernel balancing of the derivative functional is a quadratic problem with the same structure as Jiang and Huling's:

- For a kernel $k$ on standardized $(A, X)$, the worst-case imbalance over the unit ball of $k$'s native space is $n^{-2}(\omega^\top K \omega - 2\omega^\top b)$, up to a constant.
- Here $b_i = \sum_j \partial_{a_j} k(Z_i, Z_j)$. For the energy kernel $k(z, z') = -\lVert z - z' \rVert$ this is $b_i = \sum_{j \ne i} (A_i - A_j) / \lVert Z_i - Z_j \rVert$, rescaled to the units of $A$.
- With a ridge term and $\sum_i \omega_i = 0$, the solution is one linear solve. There is no sign constraint, so no QP is needed.

Checked, in a simulation with three covariates, a nonlinear treatment mean, and an outcome with $a^2 x_1$, $a x_3$, and $\sin(2x_2)$ terms ($n = 500$, 100 replications, $\lambda = 1$ untuned):

| ADE weights | Bias | SD | RMSE |
|---|---|---|---|
| Linear-GPS score (misspecified mean) | 0.014 | 0.097 | 0.097 |
| Linear-GPS score (correct mean) | 0.015 | 0.075 | 0.076 |
| Moment balancing on $(1, X, A, AX, A^2)$ | -0.005 | 0.071 | 0.071 |
| Gaussian kernel, median bandwidth | -0.120 | 0.051 | 0.131 |
| Energy kernel | -0.037 | 0.053 | 0.065 |

The energy-kernel weights are also the limit of Jiang and Huling's energy *shift* weights as the shift shrinks. Checked: fitting their QP at $\delta$ = 1, 0.5, 0.25, 0.1, and 0.05 conditional SDs and taking $(\bar Y_w - \bar Y)/\delta$ gave mean errors of -0.029, -0.035, -0.038, -0.040, and -0.040 over 12 replications. The energy ADE weights gave -0.040, with SD 0.082 for both. The ESS rose from 228 to 499 of 500 as $\delta$ shrank. At small shifts, $w_\delta \approx 1 + \delta \omega$, so the nonnegativity constraint stops binding and the problem becomes the linear one above. The constant contribution of each unit's distance to its own shifted copy cancels under $\sum_i w_i = n$.

One theoretical caveat. The derivative functional is not bounded on the energy kernel's native space, which contains functions with kinks, so the energy distance gives no finite bias bound for the ADE the way Jiang and Huling's eq. 11 does for a shift. The weights are well defined and behaved well here with a smooth outcome. A Gaussian or Matérn kernel has a bounded derivative functional, but it needs a bandwidth, and the median heuristic was badly biased in this simulation.

### Where ADE Weights Would Live

Signed weights break ESS, `trim()`, `summary()`, *cobalt*'s weighted statistics, and `glm_weightit()` (a weighted likelihood with negative weights is meaningless), so the ADE should not be an `estimand` of `weightit()`. Two options:

- A separate function returning an object of its own class, with a print method showing the constraint table, the two-group balance, and the ratio form above, and an estimator for $n^{-1} \sum_i \omega_i Y_i$ with an M-estimation standard error.
- Leaving ADE weights to *adrftools*' `ade()`, where they would augment the g-computation estimate.

Williams, Hines, and Rudolph (2025) is in the library and probably relevant, but I have not read it. The variance-weighted ADE of Hines et al. needs only $E[A \mid X]$ and $E[Y \mid X]$, so it needs nothing from *WeightIt* at all.

## Phases

1. **Additive shifts, constant or per unit.**
   - `.shift()` in *cobalt*, re-exported by *WeightIt*.
   - The numerator change in `.get_w_from_gps_internal_cont()`, with the kernel-grid fix, for `glm`, `gbm`, and `super`, and M-estimation for `glm`.
   - Target changes for `ebal` (with M-estimation) and `energy` (reusing the binary ATT code).
   - `summary()` and `print()` changes, and an error for unsupported methods.
   - Tests against the stacked-classification reference and the closed-form normal case.
   - A section in `vignette("estimating-effects")` showing weighted g-computation, and the stacking code for balance until *cobalt* catches up.
2. **Balance in *cobalt*.**
   - `bal.tab.shift`, `bal.compute()` statistics for shifts, and `bal.plot()`.
   - Then `gbm` and `super` tuning on those statistics.
3. **ADE weights (optional).**
   - Decide between a separate *WeightIt* function with its own class and leaving ADE weights to *adrftools*' `ade()`.
   - If they go in *WeightIt*: moment-balancing weights (closed form), energy-kernel weights (one linear solve), and GPS-score weights from a flexible treatment model, with the two-group balance display.
   - Doubly robust shift curves in *adrftools* from a treatment model, if wanted, are *adrftools*' own item.
4. **Other policies and methods.**
   - Multiplicative and thresholded shifts, then user policies (with inverse for density methods, forward map only for balancing methods), with stacking as the possible fallback.
   - `optweight`, `npcbps`, `cbps`, and `bart`.
   - Shifts per time point in `weightitMSM()`.

## Open Questions

- Whether `.shift()` should live in *cobalt*, as `.cens()` does. I assume yes.
- Whether ADE weights belong in *WeightIt* at all, given that with linear constraints they reproduce regression, and if so, how to choose the kernel and ridge penalty for the kernel versions.
- Whether *bartisan* can evaluate its fitted conditional density at arbitrary treatment values.
- Whether *optweight* should get a `shift` convenience on top of `optweight.svy()`.
- Whether the stacked-classification route should be exposed to users (e.g., for arbitrary policies with `gbm`) or kept internal.
- How to handle positivity: report the share of shifted values outside the range of $A$ only, or also something based on the conditional density, as *adrftools*' `_dev/feature-ideas.md` suggests.
- Whether to implement Jiang and Huling's permutation test for the shift size, and where.
