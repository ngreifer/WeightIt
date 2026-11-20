#' Weighting from Supplied Propensity Scores
#' @name method_ps
#' @usage NULL
#'
#' @description
#' This page explains the details of computing weights from propensity scores
#' that have already been estimated, requested by supplying them to the `ps`
#' argument of [weightit()] or [weightitMSM()] instead of setting `method`. This
#' can be used with binary, multi-category, and continuous treatments, as well as
#' with censoring.
#'
#' No model is fit; the supplied scores are converted to weights exactly as
#' [get_w_from_ps()] would, and the result is a full `weightit` object rather than
#' a bare vector of weights. `formula` must still contain the treatment variable,
#' but the covariates on the right hand side play no role in the computation. Any
#' `method` supplied alongside `ps` is ignored unless it is a user-supplied
#' function.
#'
#' Because nothing is estimated, no M-estimation components are produced, and
#' [glm_weightit()] will treat the weights as fixed.
#'
#' ## Binary Treatments
#'
#' `ps` may be a numeric vector of the probability of being in the treated group,
#' a one-column matrix or data frame of the same, or a two-column matrix or data
#' frame with one column per treatment level. When two columns are supplied and
#' their names match the treatment levels, those names determine which column is
#' which; otherwise the columns are assumed to be in the order of the treatment
#' levels. All estimands allowed by [get_w_from_ps()] are available: `"ATE"`,
#' `"ATT"`, `"ATC"`, `"ATO"`, `"ATM"`, and `"ATOS"`.
#'
#' ## Multi-Category Treatments
#'
#' `ps` may be a matrix or data frame with one column per treatment level (when
#' the column names match the treatment levels, those names determine which
#' column is which; otherwise the columns are assumed to be in the order of the
#' treatment levels), or a numeric vector (or one-column matrix) giving each unit's probability of being
#' in the treatment group it is actually in. The estimands `"ATE"`, `"ATT"`,
#' `"ATC"`, `"ATO"`, and `"ATM"` are available.
#'
#' ## Continuous Treatments
#'
#' `ps` is interpreted as the conditional mean of the treatment (i.e., the fitted
#' values of a treatment model) rather than as a probability. The generalized
#' propensity score is formed by feeding the standardized residuals through the
#' requested `density`, exactly as in [`method_glm`].
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, `ps` is supplied as a list with one entry per
#' time point, and the weights are the product of the time-specific weights.
#'
#' ## Censoring Weights
#'
#' Wrapping the censoring indicator in [.cens()] requests inverse probability of
#' censoring weights. `ps` is then the probability of *being censored*; the
#' weights are `1/P(C = 0 | X)` for the units still under observation and exactly
#' 0 for the censored units.
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported and are applied to the resulting weights, but
#' they play no part in computing them, since no model is fit, except for
#' continuous treatments, where they are used to estimate the residual standard
#' deviation and the marginal density of the treatment.
#'
#' @section Additional Arguments:
#' For binary and multi-category treatments, the following additional argument
#' can be specified:
#' \describe{
#'   \item{`subclass`}{`integer`; the number of subclasses to use for computing weights using marginal mean weighting through stratification (MMWS). If `NULL`, standard inverse probability weights (and their extensions) will be computed; if a number greater than 1, subclasses will be formed and weights will be computed based on subclass membership. See [get_w_from_ps()] for details and references.}
#' }
#'
#' For continuous treatments, the following additional arguments may be supplied:
#' \describe{
#'   \item{`density`}{A function corresponding to the conditional density of the treatment. The standardized residuals of the treatment will be fed through this function to produce the denominator of the generalized propensity score weights. If blank, [dnorm()] is used as recommended by Robins et al. (2000). This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if `density = "dt_2"` is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme treatment values are observed (Naimi et al., 2014).
#'
#' Can also be `"kernel"` to use kernel density estimation, which calls [density()] to estimate the denominator density for the weights. (This used to be requested by setting `use.kernel = TRUE`, which is now deprecated.)}
#'   \item{`bw`, `adjust`, `kernel`, `n`}{If `density = "kernel"`, the arguments to [density()]. The defaults are the same as those in `density()`.}
#' }
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{No fit object is produced, since no model is fit; `include.obj` has no effect.}
#' }
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' [get_w_from_ps()], which performs the same computation but returns only the
#' weights
#'
#' [`method_glm`], for estimating the propensity scores within `weightit()`
#' instead of supplying them
#'
#' @references
#' See [get_w_from_ps()] for references on each estimand and on marginal mean
#' weighting through stratification.
#'
#' Naimi, A. I., Moodie, E. E. M., Auger, N., & Kaufman, J. S. (2014).
#' Constructing Inverse Probability Weights for Continuous Exposures: A
#' Comparison of Methods. *Epidemiology*, 25(2), 292–299.
#' \doi{10.1097/EDE.0000000000000053}
#'
#' Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural
#' Models and Causal Inference in Epidemiology. *Epidemiology*, 11(5), 550–560.
#'
#' @examples
#' data("lalonde", package = "cobalt")
#'
#' # Estimate the propensity score separately, then supply it
#' fit <- glm(treat ~ age + educ + married + nodegree + re74,
#'            data = lalonde, family = binomial)
#'
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 ps = fitted(fit), estimand = "ATT"))
#'
#' summary(W1)
#'
#' # Marginal mean weighting through stratification
#' (W2 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 ps = fitted(fit), estimand = "ATT",
#'                 subclass = 20))
#'
#' summary(W2)
NULL

weightit2ps <- function(covs, treat, s.weights, subset, estimand, focal,
                        missing, ps, .data, verbose, ...) {

  fit.obj <- NULL

  n <- length(treat)
  p.score <- NULL
  treat_sub <- factor(treat[subset])

  t.lev <- .get_treated_level(treat, estimand, focal)
  c.lev <- setdiff(levels(treat_sub), t.lev)

  if (is.matrix(ps) || is.data.frame(ps)) {
    if (nrow(ps) == n) {
      if (ncol(ps) == 1L) {

        ps <- data.frame(ps[subset, 1L], 1 - ps[subset, 1L])

        names(ps) <- c(t.lev, c.lev)

        p.score <- ps[[t.lev]]
      }
      else if (ncol(ps) == 2L) {

        if (all(colnames(ps) %in% levels(treat_sub))) {
          ps <- as.data.frame(ps[subset, , drop = FALSE])
        }
        else {
          ps <- as.data.frame(ps[subset, , drop = FALSE])
          names(ps) <- levels(treat_sub)
        }

        p.score <- ps[[t.lev]]
      }
    }
  }
  else if (is.numeric(ps) && length(ps) == n) {
    ps <- data.frame(ps[subset], 1 - ps[subset])

    names(ps) <- c(t.lev, c.lev)

    p.score <- ps[[t.lev]]
  }

  if (is_null(p.score)) {
    arg::err("{.arg ps} must be a numeric vector with a propensity score for each unit")
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_bin(ps = p.score,
                                   treat = as.numeric(treat_sub == t.lev), estimand,
                                   subclass = ...get("subclass"))

  list(w = w, ps = p.score, fit.obj = fit.obj)
}

weightit2ps.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                              missing, ps, .data, verbose, ...) {

  n <- length(treat)
  treat <- factor(treat)
  treat_sub <- factor(treat[subset])

  bad.ps <- FALSE
  if (is.matrix(ps) || is.data.frame(ps)) {
    if (all(dim(ps) == c(n, nunique(treat)))) {
      #Columns named for the treatment levels are matched by name, as for binary
      #treatments; otherwise they are assumed to be in the order of the levels
      cn <- colnames(ps)
      ps <- as.data.frame(ps)

      ps <- {
        if (is_not_null(cn) && setequal(cn, levels(treat))) ps[levels(treat)]
        else setNames(ps, levels(treat))
      }

      ps <- ps[subset, , drop = FALSE]
    }
    else if (nrow(ps) == n && ncol(ps) == 1L) {
      ps <- setNames(list2DF(lapply(levels(treat), function(x) {
        p_ <- rep_with(1, treat)
        p_[treat == x] <- ps[treat == x, 1L]
        p_
      })), levels(treat))[subset, , drop = FALSE]
    }
    else {
      bad.ps <- TRUE
    }
  }
  else if (is.numeric(ps)) {
    if (length(ps) == n) {
      ps <- setNames(list2DF(lapply(levels(treat), function(x) {
        p_ <- rep_with(1, treat)
        p_[treat == x] <- ps[treat == x]
        p_
      })), levels(treat))[subset, , drop = FALSE]
    }
    else {
      bad.ps <- TRUE
    }
  }
  else {
    bad.ps <- TRUE
  }

  if (bad.ps) {
    arg::err("{.arg ps} must be a numeric vector with a propensity score for each unit or a matrix with the probability of being in each treatment for each unit")
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_multi(ps = ps, treat = treat_sub, estimand, focal = focal,
                                     subclass = ...get("subclass"))

  list(w = w)
}

weightit2ps.cens <- function(covs, treat, s.weights, subset, missing, ps, verbose,
                             estimand = NULL, focal = NULL, ...) {

  C <- .make_cens_treat(treat)

  out <- .cens_degenerate_out(C[subset])

  if (is_not_null(out)) {
    return(out)
  }

  #`ps` is the probability of being censored, so the censored units are the focal
  #("treated") group of the equivalent ATT problem. Delegating this way inherits
  #all of `weightit2ps()`'s input parsing.
  out <- weightit2ps(covs = covs, treat = C, s.weights = s.weights,
                     subset = subset, estimand = "ATT", focal = 1,
                     missing = missing, ps = ps,
                     verbose = verbose, ...)

  .att_out_to_cens(out, C[subset])
}

weightit2ps.cont <- function(covs, treat, s.weights, subset, missing, ps, verbose, ...) {

  treat <- treat[subset]
  s.weights <- s.weights[subset]

  # Process density params
  make_dens_fun <- .get_make_dens_fun(density = ...get("density"),
                                      bw = ...get("bw"),
                                      adjust = ...get("adjust"),
                                      kernel = ...get("kernel"),
                                      n = ...get("n"),
                                      use.kernel = ...get("use.kernel"))

  #Get weights
  w <- .get_w_from_gps_internal_cont(mu = ps[subset], treat = treat,
                                     s.weights = s.weights,
                                     make_dens_fun = make_dens_fun)

  list(w = w)
}
