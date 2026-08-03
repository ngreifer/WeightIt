#' Mark a censoring indicator in a model formula
#'
#' @description
#' `.cens()` marks a variable as a censoring indicator rather than a treatment, requesting inverse probability of censoring weights (IPCW). It is most often used on the left side of a formula supplied to [weightit()] or in the `formula.list` of [weightitMSM()], as in `.cens(C) ~ x1 + x2`. It can also be called directly to tag a variable for use with the lower-level interfaces [weightit.fit()] and [get_w_from_ps()].
#'
#' @param x a censoring indicator: a numeric variable taking only the values 0 (still under observation) and 1 (censored), a logical variable, or a factor with levels `0`/`1` or `FALSE`/`TRUE`.
#'
#' @details
#' Censoring is treated as its own treatment type, distinct from binary, multi-category, and continuous treatments. Weights are estimated only for the units still under observation, and are those that make the covariate distribution of the uncensored units match that of the full at-risk sample. Formally, for \eqn{e(X) = P(C = 1 | X)}, the weights are \eqn{(1 - e(X))^{-1}} for units with \eqn{C = 0} and exactly 0 for units with \eqn{C = 1}.
#'
#' Because only one group is weighted, the estimation problem is smaller and better conditioned than the corresponding binary-treatment problem, which would additionally solve for weights among the censored units. This matters most when few units are censored.
#'
#' Because censored units receive a weight of exactly 0, they contribute nothing to a weighted outcome model, and [glm_weightit()], [lm_weightit()], [multinom_weightit()], [ordinal_weightit()], and [coxph_weightit()] all tolerate missing values in the model variables for those units. This is what makes it possible to fit an outcome model whose outcome is unobserved after censoring. See the *Censoring weights (IPCW)* section of [weightit()] for details.
#'
#' No `estimand` or `focal` argument applies to censoring models; supplying one produces a warning and it is ignored. `subclass` cannot be used. Not all methods support censoring weights; see the `treat_type` component of [`.weightit_methods`] to check.
#'
#' The right side of the formula may be empty, as in `.cens(C) ~ 1`, requesting a marginal censoring model in which censoring is assumed independent of the covariates. The weights are then \eqn{1/P(C = 0)} for the units still under observation and 0 for the censored units, whatever `method` is supplied; see *Empty model formulas* in Details at [weightit()]. Nothing else changes, so a marginal censoring model still composes with `by`, `stabilize`, M-estimation, and, in [weightitMSM()], the risk sets and the missing values permitted after censoring.
#'
#' ## Use outside a formula
#'
#' Calling `.cens()` directly returns a tagged indicator that routes the lower-level interfaces to the censoring path:
#'
#' ```
#' # Equivalent to weightit(.cens(C) ~ x1 + x2, ...)
#' weightit.fit(covs, treat = .cens(C), method = "glm")
#'
#' # Censoring weights from a set of censoring probabilities
#' get_w_from_ps(ps, treat = .cens(C))
#' ```
#'
#' Tagging is the only way to request censoring weights from these interfaces; an untagged 0/1 vector is treated as an ordinary binary treatment.
#'
#' @section Note on the coding convention:
#' The survival convention is used: `1` means the unit is *censored* (drops out of observation) and `0` means it remains under observation. This is the opposite of an "observed" or "event" indicator. When a censoring model is supplied to [weightit()], the returned propensity score is \eqn{P(C = 1 | X)}, the probability of *being censored*.
#'
#' @section Assessing balance:
#' [cobalt::bal.tab()] cannot be used directly on the output of a censoring model. For a point-treatment `weightit` object it errors because every censored unit has a weight of 0, leaving one "treatment group" with no weight; for a `weightitMSM` object it errors because censoring leaves missing values in the later treatments, which `bal.tab()` does not allow.
#'
#' Because the target of a censoring model is the full at-risk sample rather than another treatment group, balance is assessed by comparing the weighted uncensored units against that sample. This can be done by stacking the two:
#'
#' ```
#' u <- which(W$treat == 0)
#' bal.tab(rbind(covs[u, ], covs),
#'         treat = c(rep(0L, length(u)), rep(1L, nrow(covs))),
#'         weights = c(W$weights[u], rep(1, nrow(covs))),
#'         estimand = "ATT", s.d.denom = "treated")
#' ```
#'
#' The "Control" group is then the weighted uncensored sample and the "Treated" group the full at-risk sample, so `Diff.Adj` is the quantity the weights are designed to zero out. For a `weightitMSM` object, use the `at.risk` component, which has one column per time point, to restrict each model to the units that were under observation when it was fit:
#'
#' ```
#' ar <- W$at.risk[, "A_3"]
#' bal.tab(A_3 ~ X1_2 + X2_2 + A_2, data = data[ar, ],
#'         weights = W$weights[ar], s.weights = W$s.weights[ar])
#' ```
#'
#' @returns
#' `x` coerced to a 0/1 numeric vector of class `treat` with a `"treat.type"` attribute of `"censoring"`. Missing values are preserved; any value other than 0 or 1 throws an error.
#'
#' Inside a formula the marker is stripped before the formula is processed, so `.cens()` is not actually evaluated there and the treatment name remains that of the indicator itself (e.g., `C` rather than `.cens(C)`).
#'
#' @seealso
#' [weightit()] and [weightitMSM()] for estimating censoring weights; [`.weightit_methods`] for which methods support them.
#'
#' @examples
#' data("msmdata")
#'
#' # A censoring indicator: 1 = lost to follow-up after time 2
#' set.seed(1234)
#' msmdata$C_2 <- rbinom(nrow(msmdata), 1,
#'                       prob = plogis(-2 + 0.2 * msmdata$X1_1))
#'
#' # Inverse probability of censoring weights
#' (Wc <- weightit(.cens(C_2) ~ X1_1 + X2_1 + A_1,
#'                 data = msmdata, method = "glm"))
#'
#' # Censored units have a weight of 0
#' all(Wc$weights[msmdata$C_2 == 1] == 0)
#'
#' # The weights are the inverse probability of remaining
#' # under observation
#' all.equal(Wc$weights[msmdata$C_2 == 0],
#'           1 / (1 - Wc$ps[msmdata$C_2 == 0]))
#'
#' # `.cens()` can also tag a variable directly for use with
#' # the lower-level interfaces
#' C <- .cens(msmdata$C_2)
#'
#' all.equal(get_w_from_ps(Wc$ps, treat = C),
#'           Wc$weights, check.attributes = FALSE)

#' @export
.cens <- function(x) {
  arg::arg_supplied(x)

  treat.name <- deparse1(substitute(x))

  #`.make_cens_treat()` validates the 0/1 domain here, at the call site, rather
  #than deep inside a fitting function. The `"treat"` class is what makes the
  #tag survive subsetting (see `[.treat`), which the `by` and longitudinal
  #risk-set paths rely on.
  out <- as.treat(.make_cens_treat(x), censoring = TRUE)

  attr(out, "treat.name") <- treat.name

  out
}
