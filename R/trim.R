#' Trim (Winsorize) Large Weights
#' @name trim
#'
#' @description Trims (i.e., winsorizes) large weights by setting all weights
#' higher than that at a given quantile to the weight at the quantile or to 0.
#' This can be useful in controlling extreme weights, which can reduce effective
#' sample size by enlarging the variability of the weights. Note that by
#' default, no observations are fully discarded when using `trim()`, which may
#' differ from the some uses of the word "trim" (see the `drop` argument below).
#'
#' @param x A `weightit` object or a vector of weights.
#' @param at `numeric`; either the quantile of the weights above which weights
#'   are to be trimmed. A single number between .5 and 1, or the number of
#'   weights to be trimmed (e.g., `at = 3` for the top 3 weights to be set to
#'   the 4th largest weight).
#' @param lower `logical`; whether also to trim at the lower quantile (e.g., for
#'   `at = .9`, trimming at both .1 and .9, or for `at = 3`, trimming the top
#'   and bottom 3 weights). Default is `FALSE` to only trim the higher weights.
#' @param treat A vector of treatment status for each unit. This should always
#'   be included when `x` is numeric, but you can get away with leaving it out
#'   if the treatment is continuous or the estimand is the ATE for binary or
#'   multi-category treatments.
#' @param drop `logical`; whether to set the weights of the trimmed units to 0
#'   or not. Default is `FALSE` to retain all trimmed units. Setting to `TRUE`
#'   may change the original targeted estimand when not the ATT or ATC.
#' @param \dots Not used.
#'
#' @returns If the input is a `weightit` object, the output will be a `weightit`
#' object with the weights replaced by the trimmed weights (or 0) and will have
#' an additional attribute, `"trim"`, equal to the quantile of trimming.
#'
#' If the input is a numeric vector of weights, the output will be a numeric
#' vector of the trimmed weights, again with the aforementioned attribute.
#'
#' @details `trim()` takes in a `weightit` object (the output of a call to
#' [weightit()] or [weightitMSM()]) or a numeric vector of weights and trims
#' (winsorizes) them to the specified quantile. All weights above that quantile
#' are set to the weight at that quantile unless `drop = TRUE`, in which case
#' they are set to 0. If `lower = TRUE`, all weights below 1 minus the quantile
#' are trimmed. In general, trimming weights can decrease balance but also
#' decreases the variability of the weights, improving precision at the
#' potential expense of unbiasedness (Cole & Hernán, 2008). See Lee, Lessler,
#' and Stuart (2011) and Thoemmes and Ong (2015) for discussions and simulation
#' results of trimming weights at various quantiles. Note that trimming weights
#' can also change the target population and therefore the estimand.
#'
#' When using `trim()` on a numeric vector of weights, it is helpful to include
#' the treatment vector as well. The helps determine the type of treatment and
#' estimand, which are used to specify how trimming is performed. In particular,
#' if the estimand is determined to be the ATT or ATC, the weights of the target
#' (i.e., focal) group are ignored, since they should all be equal to 1.
#' Otherwise, if the estimand is the ATE or the treatment is continuous, all
#' weights are considered for trimming. In general, weights for any group for
#' which all the weights are the same will not be considered in the trimming.
#'
#' @seealso [weightit()], [weightitMSM()]
#'
#' @references Cole, S. R., & Hernán, M. Á. (2008). Constructing Inverse
#' Probability Weights for Marginal Structural Models. American Journal of
#' Epidemiology, 168(6), 656–664.
#'
#' Lee, B. K., Lessler, J., & Stuart, E. A. (2011). Weight Trimming and
#' Propensity Score Weighting. PLoS ONE, 6(3), e18174.
#'
#' Thoemmes, F., & Ong, A. D. (2016). A Primer on Inverse Probability of
#' Treatment Weighting and Marginal Structural Models. Emerging Adulthood, 4(1),
#' 40–59.
#'
#' @examples
#'
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' (W <- weightit(treat ~ age + educ + married +
#'                  nodegree + re74, data = lalonde,
#'                method = "glm", estimand = "ATT"))
#' summary(W)
#'
#' #Trimming the top and bottom 5 weights
#' trim(W, at = 5, lower = TRUE)
#'
#' #Trimming at 90th percentile
#' (W.trim <- trim(W, at = .9))
#'
#' summary(W.trim)
#' #Note that only the control weights were trimmed
#'
#' #Trimming a numeric vector of weights
#' all.equal(trim(W$weights, at = .9, treat = lalonde$treat),
#'           W.trim$weights)
#'
#' #Dropping trimmed units
#' (W.trim <- trim(W, at = .9, drop = TRUE))
#'
#' summary(W.trim)
#' #Note that we now have zeros in the control group
#'
#' #Using made up data and as.weightit()
#' treat <- rbinom(500, 1, .3)
#' weights <- rchisq(500, df = 2)
#' W <- as.weightit(weights, treat = treat,
#'                  estimand = "ATE")
#' summary(W)
#' summary(trim(W, at = .95))

#' @export
trim <- function(x, ...) {
  UseMethod("trim")
}

#' @exportS3Method trim weightit
#' @rdname trim
trim.weightit <- function(x, at = 0, lower = FALSE, drop = FALSE, ...) {

  if (is_null(at)) {
    return(x)
  }

  chk::chk_number(at)
  chk::chk_gte(at, 0)

  if (check_if_zero(at)) {
    return(x)
  }

  chk::chk_flag(lower)
  chk::chk_flag(drop)

  x[["weights"]] <- .trim_weights(x[["weights"]],
                                 at = at,
                                 treat = x[["treat"]],
                                 groups.not.to.trim = x[["focal"]],
                                 lower = lower,
                                 drop = drop)
  x
}

#' @exportS3Method trim default
#' @rdname trim
trim.default <- function(x, at = 0, lower = FALSE, treat = NULL, drop = FALSE, ...) {

  if (!is.numeric(x) && !inherits(x, "weightit")) {
    .err("`x` must be a numeric vector or `weightit` object")
  }

  if (is_null(at)) {
    return(x)
  }

  chk::chk_number(at)
  chk::chk_gte(at, 0)

  if (check_if_zero(at)) {
    return(x)
  }

  if (is_not_null(treat)) {
    if (!has_treat_type(treat)) {
      treat <- assign_treat_type(treat)
    }
    treat.type <- get_treat_type(treat)
  }

  groups.not.to.trim <- NULL

  if (is_not_null(treat) && treat.type != "continuous") {
    w_all_same <- tapply(x, treat, all_the_same)
    if (all(w_all_same)) {
      .wrn("weights are all the same in each treatment group and will not be trimmed")
      return(x)
    }

    if (any(w_all_same)) {
      groups.not.to.trim <- names(w_all_same)[w_all_same]
    }
  }
  else if (all_the_same(x)) {
    .wrn("weights are all the same and will not be trimmed")
    return(x)
  }

  chk::chk_flag(lower)
  chk::chk_flag(drop)

  .trim_weights(x, at = at,
               treat = treat,
               groups.not.to.trim = groups.not.to.trim,
               lower = lower,
               drop = drop)
}

#Internal function that does the work
.trim_weights <- function(weights, at = 0, treat = NULL, groups.not.to.trim = NULL,
                          lower = FALSE, drop = FALSE) {

  if (at < 1) {
    at <- max(at, 1 - at)

    trim.q <- {
      if (lower) c(1 - at, at)
      else c(0, at)
    }

    if (is_not_null(groups.not.to.trim)) {
      to.be.trimmed <- treat %nin% groups.not.to.trim
      trim.w <- quantile(weights[to.be.trimmed], probs = trim.q, type = 3)

      if (drop) {
        weights[to.be.trimmed & (weights < trim.w[1L] | weights > trim.w[2L])] <- 0
        .msg(sprintf("setting weights beyond %s where treat is not %s to 0",
                     word_list(paste0(round(100 * trim.q[c(lower, TRUE)], 2L), "%")),
                     word_list(groups.not.to.trim, and.or = "or")))
      }
      else {
        weights[to.be.trimmed & weights < trim.w[1L]] <- trim.w[1L]
        weights[to.be.trimmed & weights > trim.w[2L]] <- trim.w[2L]
        .msg(sprintf("trimming weights where treat is not %s to %s",
                     word_list(groups.not.to.trim, and.or = "or"),
                     word_list(paste0(round(100 * trim.q[c(lower, TRUE)], 2L), "%"))))
      }

    }
    else {
      trim.w <- quantile(weights, probs = trim.q, type = 3)

      if (drop) {
        weights[weights < trim.w[1L] | weights > trim.w[2L]] <- 0
        .msg(sprintf("Setting weights beyond %s to 0",
                     word_list(paste0(round(100 * trim.q[c(lower, TRUE)], 2L), "%"))))
      }
      else {
        weights[weights < trim.w[1L]] <- trim.w[1L]
        weights[weights > trim.w[2L]] <- trim.w[2L]

        .msg(sprintf("Trimming weights to %s",
                     word_list(paste0(round(100 * trim.q[c(lower, TRUE)], 2L), "%"))))
      }

    }
  }
  else {
    if (is_not_null(groups.not.to.trim)) {
      to.be.trimmed <- treat %nin% groups.not.to.trim
      if (at >= sum(to.be.trimmed)) {
        .wrn(sprintf("`at` must be less than %s, the number of units for which treat is not %s. Weights will not be trimmed",
                     sum(to.be.trimmed),
                     word_list(groups.not.to.trim, and.or = "or")))
        return(weights)
      }

      at <- as.integer(min(at, sum(to.be.trimmed) - at))

      trim.top <- {
        if (lower) c(at + 1, sum(to.be.trimmed) - at)
        else c(1, sum(to.be.trimmed) - at)
      }

      weights.text <- {
        if (at > 1) sprintf("%s weights", at)
        else if (lower) "weights"
        else "weight"
      }

      trim.w <- sort(weights[to.be.trimmed], partial = trim.top)[trim.top]

      if (drop) {
        weights[to.be.trimmed & (weights < trim.w[1L] | weights > trim.w[2L])] <- 0

        .msg(sprintf("setting the %s %s where treat is not %s to 0",
                     word_list(c("top", "bottom")[c(TRUE, lower)]),
                     weights.text,
                     word_list(groups.not.to.trim, and.or = "or")))
      }
      else {
        weights[to.be.trimmed & weights < trim.w[1L]] <- trim.w[1L]
        weights[to.be.trimmed & weights > trim.w[2L]] <- trim.w[2L]

        .msg(sprintf("trimming the %s %s where treat is not %s",
                     word_list(c("top", "bottom")[c(TRUE, lower)]),
                     weights.text,
                     word_list(groups.not.to.trim, and.or = "or")))
      }
    }
    else {
      if (at >= length(weights)) {
        .wrn(sprintf("`at` must be less than %s, the number of units. Weights will not be trimmed",
                     length(weights)))
        return(weights)
      }

      at <- as.integer(min(at, length(weights) - at))

      trim.top <- {
        if (lower) c(at + 1, length(weights) - at)
        else c(1, length(weights) - at)
      }

      weights.text <- {
        if (at > 1) sprintf("%s weights", at)
        else if (lower) "weights"
        else "weight"
      }

      trim.w <- sort(weights)[trim.top]

      if (drop) {
        weights[weights < trim.w[1L] | weights > trim.w[2L]] <- 0

        .msg(sprintf("setting the %s %s to 0",
                     word_list(c("top", "bottom")[c(TRUE, lower)]),
                     weights.text))
      }
      else {
        weights[weights < trim.w[1L]] <- trim.w[1L]
        weights[weights > trim.w[2L]] <- trim.w[2L]

        .msg(sprintf("trimming the %s %s",
                     word_list(c("top", "bottom")[c(TRUE, lower)]),
                     weights.text))
      }
    }
  }

  attr(weights, "trim") <- at
  attr(weights, "trim.lower") <- lower

  weights
}
