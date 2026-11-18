# `hatvalues()` methods for `ordinal_weightit` and `multinom_weightit` (2026-11-16)

*adrftools* (development version) now computes its default unconditional variance with the HC2 adjustment, dividing each unit's influence on the model coefficients by $\sqrt{1 - h_i}$, where $h_i$ is the unit's leverage. It reads $h_i$ from the QR decomposition for `lm()` and `glm()` fits (which covers `glm_weightit()` and `lm_weightit()`) and calls `hatvalues()` for anything else. `ordinal_weightit()` and `multinom_weightit()` fits have no `hatvalues()` method, so *adrftools* falls back to HC0 for them with a message. Both models can have one, and adding it would let *adrftools* (and `sandwich::vcovHC(type = "HC2")`-style code generally) use the adjustment with no change on its side.

## The generalized leverage

For a model fit by maximizing $\sum_i w_i \ell_i(\theta)$, the natural generalization of the hat value is each unit's share of the information:

$$h_i = \mathrm{tr}\left(\mathcal{I}^{-1} \mathcal{I}_i\right), \qquad \mathcal{I} = \sum_i \mathcal{I}_i,$$

where $\mathcal{I}_i$ is unit $i$'s contribution to the expected (Fisher) information. These sum to the number of coefficients, $\sum_i h_i = p$, and for a GLM fit by iteratively reweighted least squares, $\mathcal{I}_i = w_i x_i x_i'$ with $w_i$ the working weight, so $h_i = w_i x_i' (X'WX)^{-1} x_i$, the usual hat value, for any link. Use the expected information rather than the observed: the observed contributions can be negative, and the expected ones reproduce `hatvalues()` for GLMs.

Units with zero weight get $h_i = 0$; aliased coefficients are dropped from $\theta$ before computing $\mathcal{I}$, as in `.get_hess_multinom()` and `.get_hess_ordinal()`. The returned vector should have one value per row of `estfun()`, in the same order, since *adrftools* checks that the lengths match.

## `multinom_weightit`

With $K$ non-reference categories, fitted probabilities $p_{ik}$, and coefficients $\beta_1, \dots, \beta_K$ (one block of length $q$ per category), unit $i$'s information block for categories $k$ and $l$ is

$$\mathcal{I}_{i,kl} = w_i \left(\delta_{kl}\, p_{ik} - p_{ik} p_{il}\right) x_i x_i'.$$

This is the per-unit term inside the `crossprod()` calls in `.get_hess_multinom()` (the multinomial logit uses the canonical link, so observed and expected information agree). With $\mathcal{I}^{-1}$ partitioned into $q \times q$ blocks $V_{lk}$,

$$h_i = w_i \sum_{k, l} \left(\delta_{kl}\, p_{ik} - p_{ik} p_{il}\right) x_i' V_{lk}\, x_i,$$

which vectorizes as `rowSums((x_ %*% V[l_ind, k_ind]) * x_)` for each pair $(k, l)$, multiplied by the probability term and summed: $O(n q^2 K^2)$, with $\mathcal{I}$ itself being the negative of the Hessian that `.get_hess_multinom()` already returns.

## `ordinal_weightit`

With linear predictor $\eta_i = x_i'\beta + \mathrm{offset}_i$, thresholds $\alpha_1 < \dots < \alpha_{m-1}$ (with $\alpha_0 = -\infty$ and $\alpha_m = \infty$), inverse link $F$ with density $f$, and category probabilities $\pi_{ic} = F(\alpha_c - \eta_i) - F(\alpha_{c-1} - \eta_i)$, the derivatives of $\pi_{ic}$ with respect to $\theta = (\beta, \alpha)$ are

$$\frac{\partial \pi_{ic}}{\partial \beta} = -x_i \left[f(\alpha_c - \eta_i) - f(\alpha_{c-1} - \eta_i)\right], \qquad \frac{\partial \pi_{ic}}{\partial \alpha_c} = f(\alpha_c - \eta_i), \qquad \frac{\partial \pi_{ic}}{\partial \alpha_{c-1}} = -f(\alpha_{c-1} - \eta_i),$$

and the unit's expected information is $\mathcal{I}_i = w_i \sum_{c=1}^{m} d_{ic} d_{ic}' / \pi_{ic}$ with $d_{ic} = \partial \pi_{ic} / \partial \theta$. Then

$$h_i = w_i \sum_{c=1}^{m} \frac{d_{ic}' \mathcal{I}^{-1} d_{ic}}{\pi_{ic}},$$

computed for each category as `rowSums((D_c %*% Iinv) * D_c) / pi_c` with `D_c` the $n \times p$ matrix of the $d_{ic}$, at $O(n p^2 m)$. This should use the natural (not the cumulative-sum) parameterization of the thresholds, and $\mathcal{I}$ should be the expected information computed the same way, not the `optimHess()` Hessian of `.get_hess_ordinal()`, which is observed and on the internal scale. `fam$linkinv` and `fam$mu.eta` already give $F$ and $f$.

## Notes

- **HC2 with a scalar leverage is an approximation for these models.** Each unit has a vector of residuals (one per category), and the exact analogue of HC2 adjusts that vector by $(I - H_{ii})^{-1/2}$, with $H_{ii}$ the unit's diagonal block of the generalized hat matrix. A scalar $h_i = \mathrm{tr}(H_{ii})$ scales the whole unit by $1/\sqrt{1 - h_i}$, which removes the average own-observation bias but not its distribution across categories. That is still a clear improvement over HC0, and it is what *adrftools* can use without knowing the model's structure.
- **Register them** as `S3method(stats::hatvalues, multinom_weightit)` and `S3method(stats::hatvalues, ordinal_weightit)`.
- **The multinomial formula was checked** (2026-11-16) on `nhanes3lead` from *adrftools*: with a two-category outcome, it matched `hatvalues()` of the equivalent `glm()` to $8 \times 10^{-9}$ (the fitting tolerance), and with two and three categories the leverages summed to the number of coefficients (4 and 8). The ordinal formula was checked on 2026-11-18 (see below).
- **Implemented** (2026-11-18) in `R/glm_weightit-methods.R`, reusing `.multinom_info()` and `.ordinal_info()` for $\mathcal{I}$. Both were checked against `VGAM::hatvalues()` on `vglm()` fits with `multinomial(refLevel = 1)` and `cumulative(parallel = TRUE)` (logit, probit, cloglog), summing its $n \times M$ matrix across columns: agreement to $5 \times 10^{-9}$ (multinomial) and $10^{-7}$ (ordinal, the level at which the two fits' coefficients differ). Both also match a brute-force $\mathrm{tr}(\mathcal{I}^{-1} \mathcal{I}_i)$ built from finite differences of the fitted probabilities to $10^{-11}$, and `glm()` for two categories to $10^{-9}$. These checks are in `tests/testthat/test-hatvalues.R` (without *VGAM*, which is not in Suggests).
- **The parameterization does not matter.** $\mathrm{tr}(\mathcal{I}^{-1} \mathcal{I}_i)$ is invariant to any smooth reparameterization ($\mathcal{I} \to J' \mathcal{I} J$), so the cumulative-sum thresholds or the standardized covariates would give the same values. Only the use of expected rather than observed information matters.
- **$h_i$ can exceed 1.** The unit's block $H_{ii}$ has eigenvalues in $[0, 1]$ and dimension $m - 1$, so $0 \le h_i \le m - 1$. *adrftools*' `.hat_values()` rejects any vector with an element above 1 and falls back to HC0, which is right for the scalar adjustment ($1 / \sqrt{1 - h_i}$ is undefined there) but means a single extreme unit turns HC2 off for the whole fit. With moderate $n$ every $h_i$ is far below 1 (max .08 on `lalonde` with 8 coefficients).
- **Check them** against the identity $\sum_i h_i = p$, against `hatvalues()` of an equivalent `glm()` when the outcome has two categories (both models then reduce to logistic regression), and against a brute-force $\mathrm{tr}(\mathcal{I}^{-1} \mathcal{I}_i)$ with $\mathcal{I}_i$ built unit by unit on a small dataset.
