#' Subgroup Balancing Propensity Score
#'
#' @description
#' Implements the subgroup balancing propensity score (SBPS), which is an
#' algorithm that attempts to achieve balance in subgroups by sharing
#' information from the overall sample and subgroups (Dong, Zhang, Zeng, & Li,
#' 2020; DZZL). Each subgroup can use either weights estimated using the whole
#' sample, weights estimated using just that subgroup, or a combination of the
#' two. The optimal combination is chosen as that which minimizes an imbalance
#' criterion that includes subgroup as well as overall balance.
#'
#' @param obj a `weightit` object containing weights estimated in the
#' overall sample.
#' @param obj2 a `weightit` object containing weights estimated in the
#' subgroups. Typically this has been estimated by including `by` in the
#' call to [weightit()]. Either `obj2` or `moderator` must be
#' specified.
#' @param moderator optional; a string containing the name of the variable in
#' `data` for which weighting is to be done within subgroups or a
#' one-sided formula with the subgrouping variable on the right-hand side. This
#' argument is analogous to the `by` argument in `weightit()`, and in
#' fact it is passed on to `by`. Either `obj2` or `moderator`
#' must be specified.
#' @param formula an optional formula with the covariates for which balance is
#' to be optimized. If not specified, the formula in `obj$call` will be
#' used.
#' @param data an optional data set in the form of a data frame that contains
#' the variables in `formula` or `moderator`.
#' @param smooth `logical`; whether the smooth version of the SBPS should
#' be used. This is only compatible with `weightit` methods that return a
#' propensity score.
#' @param full.search `logical`; when `smooth = FALSE`, whether every
#' combination of subgroup and overall weights should be evaluated. If
#' `FALSE`, a stochastic search as described in DZZL will be used instead.
#' If `TRUE`, all \eqn{2^R} combinations will be checked, where \eqn{R} is
#' the number of subgroups, which can take a long time with many subgroups. If
#' unspecified, will default to `TRUE` if \eqn{R <= 8} and `FALSE`
#' otherwise.
#'
#' @returns
#' A `weightit.sbps` object, which inherits from `weightit`.
#' This contains all the information in `obj` with the weights, propensity
#' scores, call, and possibly covariates updated from `sbps()`. In
#' addition, the `prop.subgroup` component contains the values of the
#' coefficients C for the subgroups (which are either 0 or 1 for the standard
#' SBPS), and the `moderator` component contains a data.frame with the
#' moderator.
#'
#' This object has its own summary method and is compatible with \pkg{cobalt}
#' functions. The `cluster` argument should be used with \pkg{cobalt}
#' functions to accurately reflect the performance of the weights in balancing
#' the subgroups.
#'
#' @details
#' The SBPS relies on two sets of weights: one estimated in the overall sample
#' and one estimated within each subgroup. The algorithm decides whether each
#' subgroup should use the weights estimated in the overall sample or those
#' estimated in the subgroup. There are 2^R permutations of overall and
#' subgroup weights, where R is the number of subgroups. The optimal
#' permutation is chosen as that which minimizes a balance criterion as
#' described in DZZL. The balance criterion used here is, for binary and
#' multi-category treatments, the sum of the squared standardized mean differences
#' within subgroups and overall, which are computed using
#' [cobalt::col_w_smd()], and for continuous treatments, the
#' sum of the squared correlations between each covariate and treatment within
#' subgroups and overall, which are computed using [cobalt::col_w_corr()].
#'
#' The smooth version estimates weights that determine the relative
#' contribution of the overall and subgroup propensity scores to a weighted
#' average propensity score for each subgroup. If P_O are the propensity scores
#' estimated in the overall sample and P_S are the propensity scores estimated
#' in each subgroup, the smooth SBPS finds R coefficients C so that for each
#' subgroup, the ultimate propensity score is \eqn{C*P_S + (1-C)*P_O}, and
#' weights are computed from this propensity score. The coefficients are
#' estimated using [optim()] with `method = "L-BFGS-B"`. When C is
#' estimated to be 1 or 0 for each subgroup, the smooth SBPS coincides with the
#' standard SBPS.
#'
#' If `obj2` is not specified and `moderator` is, `sbps()` will
#' attempt to refit the model specified in `obj` with the `moderator`
#' in the `by` argument. This relies on the environment in which
#' `obj` was created to be intact and can take some time if `obj` was
#' hard to fit. It's safer to estimate `obj` and `obj2` (the latter
#' simply by including the moderator in the `by` argument) and supply
#' these to `sbps()`.
#'
#' @seealso
#' [weightit()], [summary.weightit()]
#'
#' @references
#' Dong, J., Zhang, J. L., Zeng, S., & Li, F. (2020). Subgroup
#' balancing propensity score. Statistical Methods in Medical Research, 29(3),
#' 659–676. \doi{10.1177/0962280219870836}
#'
#' @examples
#'
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups within races
#' (W1 <- weightit(treat ~ age + educ + married +
#'                 nodegree + race + re74, data = lalonde,
#'                 method = "glm", estimand = "ATT"))
#'
#' (W2 <- weightit(treat ~ age + educ + married +
#'                 nodegree + race + re74, data = lalonde,
#'                 method = "glm", estimand = "ATT",
#'                 by = "race"))
#' S <- sbps(W1, W2)
#' print(S)
#' summary(S)
#' bal.tab(S, cluster = "race")
#'
#' #Could also have run
#' #  sbps(W1, moderator = "race")
#'
#' S_ <- sbps(W1, W2, smooth = TRUE)
#' print(S_)
#' summary(S_)
#' bal.tab(S_, cluster = "race")

#' @export
sbps <- function(obj, obj2 = NULL, moderator = NULL, formula = NULL, data = NULL, smooth = FALSE, full.search) {

  if (is_null(obj2) && is_null(moderator)) {
    .err("either `obj2` or `moderator` must be specified")
  }

  treat <- obj[["treat"]]
  treat.type <- get_treat_type(treat)

  focal <- obj[["focal"]]
  estimand <- obj[["estimand"]]

  data.list <- list(data, obj2[["covs"]], obj[["covs"]])
  combined.data <- do.call("data.frame", clear_null(data.list))
  processed.moderator <- process.by(moderator, data = clear_null(combined.data),
                                    treat = obj[["treat"]], treat.name = NULL,
                                    by.arg = "moderator")
  moderator.factor <- attr(processed.moderator, "by.factor")

  if (is_not_null(obj2)) {
    if (!inherits(obj2, "weightit")) {
      .err("`obj2` must be a `weightit` object, ideally with a 'by' component")
    }
    else if (is_not_null(obj2[["by"]])) {
      if (is_not_null(obj[["by"]])) {
        if (is_null(processed.moderator)) .err("cannot figure out moderator. Please supply a value to `moderator`")
      }
      else {
        processed.moderator <- obj2[["by"]]

        moderator.factor <- attr(processed.moderator, "by.factor")
      }
    }
    else if (is_null(processed.moderator)) {
      .err("no moderator was specified")
    }
  }
  else {
    call <- obj[["call"]]

    if (is_not_null(obj[["by"]])) {
      call[["by"]] <- setNames(data.frame(factor(paste(processed.moderator[[1]],
                                                       obj[["by"]][[1]], sep = " | "))),
                               paste(names(processed.moderator), names(obj[["by"]]), sep = " | "))

    }
    else {
      call[["by"]] <- processed.moderator
    }

    obj2 <- eval(call)
  }

  if ((is_null(obj[["ps"]]) || is_null(obj2[["ps"]])) && smooth) {
    .err("smooth SBPS can only be used with methods that produce a propensity score")
  }

  call <- obj[["call"]]
  if (is_null(formula) && "formula" %in% names(call)) {
    formula <- eval(call[["formula"]])
  }
  formula <- delete.response(terms(formula))

  t.c <- get_covs_and_treat_from_formula(formula, combined.data)
  if (is_null(t.c[["reported.covs"]])) {
    .err("No covariates were found")
  }
  covs <- t.c[["model.covs"]]
  s.weights <- obj[["s.weights"]]

  mod.split <- cobalt::splitfactor(moderator.factor, drop.first = "if2")
  same.as.moderator <- apply(covs, 2, function(c) {
    any(vapply(mod.split, function(x) equivalent.factors(x, c), logical(1L)))
  })
  covs <- covs[, !same.as.moderator, drop = FALSE]

  bin.vars <- is_binary_col(covs)
  s.d.denom <- get.s.d.denom.weightit(estimand = obj[["estimand"]], weights = obj[["weights"]],
                                      treat = treat)

  R <- levels(moderator.factor)

  if (smooth) {
    if (!missing(full.search)) {
      .wrn("`full.search` is ignored when `smooth = TRUE`")
    }

    ps_o <- obj[["ps"]]
    ps_s <- obj2[["ps"]]

    get_w_smooth <- function(coefs, moderator.factor, treat, ps_o, ps_s, estimand) {
      ind.coefs <- coefs[moderator.factor] #Gives each unit the coef for their subgroup
      ps_ <- (1 - ind.coefs) * ps_o + ind.coefs * ps_s

      get_w_from_ps(ps_, treat, estimand)
    }

    get_F_smooth <- function(ps_o, ps_s, treat.type, ...) {
      coefs <- unlist(list(...))
      w_ <- get_w_smooth(coefs, moderator.factor, treat, ps_o, ps_s, estimand = obj[["estimand"]])

      if (treat.type == "binary") {
        F0_o <- cobalt::col_w_smd(covs, treat, w_, std = TRUE, s.d.denom = s.d.denom,
                                  abs = TRUE, s.weights = s.weights, bin.vars = bin.vars)
        F0_s <- unlist(lapply(R, function(g) cobalt::col_w_smd(covs[moderator.factor == g, , drop = FALSE],
                                                               treat[moderator.factor == g], w_[moderator.factor == g],
                                                               std = TRUE, s.d.denom = s.d.denom,
                                                               abs = TRUE, s.weights = s.weights[moderator.factor == g],
                                                               bin.vars = bin.vars)))
      }
      else if (treat.type == "multinomial") {
        if (is_not_null(focal)) {
          bin.treat <- as.numeric(treat == focal)
          s.d.denom <- switch(estimand, ATT = "treated", ATC = "control", "all")
          F0_o <- unlist(lapply(levels(treat)[levels(treat) != focal], function(t) {
            cobalt::col_w_smd(covs, bin.treat, w_, std = TRUE, s.d.denom = s.d.denom,
                              abs = TRUE, s.weights = s.weights, bin.vars = bin.vars,
                              subset = treat %in% c(t, focal))
          }))
          F0_s <- unlist(lapply(levels(treat)[levels(treat) != focal], function(t) {
            unlist(lapply(R, function(g) cobalt::col_w_smd(covs[moderator.factor == g, , drop = FALSE],
                                                           bin.treat[moderator.factor == g], w_[moderator.factor == g],
                                                           std = TRUE, s.d.denom = s.d.denom,
                                                           abs = TRUE, s.weights = s.weights[moderator.factor == g],
                                                           bin.vars = bin.vars,
                                                           subset = treat[moderator.factor == g] %in% c(t, focal))))
          }))

        }
        else {
          F0_o <- unlist(lapply(levels(treat), function(t) {
            covs_i <- rbind(covs, covs[treat == t, , drop = FALSE])
            treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat == t)))
            w_i <- c(rep(1, nrow(covs)), w_[treat == t])
            if (is_not_null(s.weights)) s.weights_i <- c(s.weights, s.weights[treat == t])
            else s.weights_i <- NULL
            cobalt::col_w_smd(covs_i, treat_i, w_i, std = TRUE, s.d.denom = "treated",
                              abs = TRUE, s.weights = s.weights_i, bin.vars = bin.vars)
          }))
          F0_s <- unlist(lapply(levels(treat), function(t) {
            covs_i <- rbind(covs, covs[treat == t, , drop = FALSE])
            treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat == t)))
            w_i <- c(rep(1, nrow(covs)), w_[treat == t])
            moderator.factor_i <- c(moderator.factor, moderator.factor[treat == t])
            if (is_not_null(s.weights)) s.weights_i <- c(s.weights, s.weights[treat == t])
            else s.weights_i <- NULL
            unlist(lapply(R, function(g) cobalt::col_w_smd(covs_i[moderator.factor_i == g, , drop = FALSE],
                                                           treat_i[moderator.factor_i == g], w_i[moderator.factor_i == g],
                                                           std = TRUE, s.d.denom = "treated",
                                                           abs = TRUE, s.weights = s.weights_i[moderator.factor_i == g],
                                                           bin.vars = bin.vars)))
          }))
        }
      }
      else if (treat.type == "continuous") {
        F0_o <- cobalt::col_w_corr(covs, treat, w_, abs = TRUE, s.weights = s.weights, bin.vars = bin.vars)
        F0_s <- unlist(lapply(R, function(g) cobalt::col_w_corr(covs[moderator.factor == g, , drop = FALSE],
                                                                treat[moderator.factor == g], w_[moderator.factor == g],
                                                                abs = TRUE, s.weights = s.weights[moderator.factor == g],
                                                                bin.vars = bin.vars)))
      }

      # F0_g <- cobalt::col_w_smd(cobalt::splitfactor(moderator.factor, drop.first = FALSE),
      #                           treat, w_, std = FALSE,
      #                           abs = TRUE, s.weights = s.weights,
      #                           bin.vars = rep(TRUE, length(R)))

      sum(F0_o^2) + sum(F0_s^2) #+ sum(F0_g^2)
    }

    opt.out <- optim(rep(.5, length(R)), fn = get_F_smooth,
                     ps_o = ps_o, ps_s = ps_s, treat.type = treat.type,
                     lower = 0, upper = 1,
                     method = "L-BFGS-B")

    s_min <- setNames(opt.out$par, R) #coef is proportion subgroup vs. overall
    weights <- get_w_smooth(s_min, moderator.factor, treat, ps_o, ps_s, estimand = obj[["estimand"]])
    ps <- (1 - s_min[moderator.factor]) * ps_o + s_min[moderator.factor] * ps_s
  }
  else {
    w_o <- obj[["weights"]]
    w_s <- obj2[["weights"]]

    if (missing(full.search)) {
      full.search <- (length(R) <= 8)
    }
    else {
      chk::chk_flag(full.search)
    }

    get_w <- function(s, moderator.factor, w_o, w_s) {
      #Get weights for given permutation of "O" and "S"
      w_ <- numeric(length(moderator.factor))
      for (g in levels(moderator.factor)) {
        if (s[g] == 0) w_[moderator.factor == g] <- w_o[moderator.factor == g]
        else if (s[g] == 1) w_[moderator.factor == g] <- w_s[moderator.factor == g]
      }

      w_
    }

    get_F <- function(s, moderator.factor, w_o, w_s, treat.type) {
      #Get value of loss function for given permutation of "O" and "S"
      w_ <- get_w(s, moderator.factor, w_o, w_s)

      if (treat.type == "binary") {
        F0_o <- cobalt::col_w_smd(covs, treat, w_, std = TRUE, s.d.denom = s.d.denom,
                                  abs = TRUE, s.weights = s.weights, bin.vars = bin.vars)
        F0_s <- unlist(lapply(R, function(g) cobalt::col_w_smd(covs[moderator.factor == g, , drop = FALSE],
                                                               treat[moderator.factor == g], w_[moderator.factor == g],
                                                               std = TRUE, s.d.denom = s.d.denom,
                                                               abs = TRUE, s.weights = s.weights[moderator.factor == g],
                                                               bin.vars = bin.vars)))
      }
      else if (treat.type == "multinomial") {
        if (is_not_null(focal)) {
          bin.treat <- as.numeric(treat == focal)
          s.d.denom <- switch(estimand, ATT = "treated", ATC = "control", "all")
          F0_o <- unlist(lapply(levels(treat)[levels(treat) != focal], function(t) {
            cobalt::col_w_smd(covs, bin.treat, w_, std = TRUE, s.d.denom = s.d.denom,
                              abs = TRUE, s.weights = s.weights, bin.vars = bin.vars,
                              subset = treat %in% c(t, focal))
          }))
          F0_s <- unlist(lapply(levels(treat)[levels(treat) != focal], function(t) {
            unlist(lapply(R, function(g) cobalt::col_w_smd(covs[moderator.factor == g, , drop = FALSE],
                                                           bin.treat[moderator.factor == g], w_[moderator.factor == g],
                                                           std = TRUE, s.d.denom = s.d.denom,
                                                           abs = TRUE, s.weights = s.weights[moderator.factor == g],
                                                           bin.vars = bin.vars,
                                                           subset = treat[moderator.factor == g] %in% c(t, focal))))
          }))

        }
        else {
          F0_o <- unlist(lapply(levels(treat), function(t) {
            covs_i <- rbind(covs, covs[treat == t, , drop = FALSE])
            treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat == t)))
            w_i <- c(rep(1, nrow(covs)), w_[treat == t])
            if (is_not_null(s.weights)) s.weights_i <- c(s.weights, s.weights[treat == t])
            else s.weights_i <- NULL
            cobalt::col_w_smd(covs_i, treat_i, w_i, std = TRUE, s.d.denom = "treated",
                              abs = TRUE, s.weights = s.weights_i, bin.vars = bin.vars)
          }))
          F0_s <- unlist(lapply(levels(treat), function(t) {
            covs_i <- rbind(covs, covs[treat == t, , drop = FALSE])
            treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat == t)))
            w_i <- c(rep(1, nrow(covs)), w_[treat == t])
            moderator.factor_i <- c(moderator.factor, moderator.factor[treat == t])
            if (is_not_null(s.weights)) s.weights_i <- c(s.weights, s.weights[treat == t])
            else s.weights_i <- NULL
            unlist(lapply(R, function(g) cobalt::col_w_smd(covs_i[moderator.factor_i == g, , drop = FALSE],
                                                           treat_i[moderator.factor_i == g], w_i[moderator.factor_i == g],
                                                           std = TRUE, s.d.denom = "treated",
                                                           abs = TRUE, s.weights = s.weights_i[moderator.factor_i == g],
                                                           bin.vars = bin.vars)))
          }))
        }
      }
      else if (treat.type == "continuous") {
        F0_o <- cobalt::col_w_corr(covs, treat, w_, abs = TRUE, s.weights = s.weights, bin.vars = bin.vars)
        F0_s <- unlist(lapply(R, function(g) cobalt::col_w_corr(covs[moderator.factor == g, , drop = FALSE],
                                                                treat[moderator.factor == g], w_[moderator.factor == g],
                                                                abs = TRUE, s.weights = s.weights[moderator.factor == g],
                                                                bin.vars = bin.vars)))
      }

      # F0_g <- cobalt::col_w_smd(cobalt::splitfactor(moderator.factor, drop.first = FALSE),
      #                           treat, w_, std = FALSE,
      #                           abs = TRUE, s.weights = s.weights,
      #                           bin.vars = rep(TRUE, length(R)))

      sum(F0_o^2) + sum(F0_s^2) #+ sum(F0_g^2)
    }

    if (full.search) {
      S <- as.matrix(setNames(do.call("expand.grid", replicate(length(R), c(0, 1), simplify = FALSE)),
                              R))

      F_min <- Inf

      for (i in seq_row(S)) {
        s_try <- S[i,]
        F_try <- get_F(s_try, moderator.factor, w_o, w_s, treat.type)

        if (F_try < F_min) {
          F_min <- F_try
          s_min <- s_try
        }
      }
    }
    else {
      #Stochastic search described by Dong et al (2019)

      s_min <- setNames(rep(0, length(R), replace = TRUE), R)
      F_min <- get_F(s_min, moderator.factor, w_o, w_s, treat.type)

      L1 <- 25
      L2 <- 10

      k <- 0
      iters_since_change <- 0

      while (k < L1 || iters_since_change < L2) {
        s_try <- setNames(sample(c(0, 1), length(R), replace = TRUE), R)
        F_try <- get_F(s_try, moderator.factor, w_o, w_s, treat.type)


        Ar <- sample(R)
        #Optimize s_try for given Ar
        repeat {
          s_try_prev <- s_try
          for (i in Ar) {
            s_alt <- s_try; s_alt[i] <- ifelse(s_try[i] == 0, 1, 0)
            F_alt <- get_F(s_alt, moderator.factor, w_o, w_s, treat.type)
            if (F_alt < F_try) {
              s_try <- s_alt
              F_try <- F_alt
            }
          }
          if (identical(s_try_prev, s_try)) break
        }

        if (F_try < F_min) {
          F_min <- F_try
          s_min <- s_try
          iters_since_change <- 0
        }
        else {
          iters_since_change <- iters_since_change + 1
        }

        k <- k + 1
      }
    }

    weights <- get_w(s_min, moderator.factor, w_o, w_s)
    if (is_not_null(obj[["ps"]]) && is_not_null(obj2[["ps"]])) {
      ps <- get_w(s_min, moderator.factor, obj[["ps"]], obj2[["ps"]])
    }
    else ps <- NULL
  }

  out <- obj
  out[["covs"]] <- t.c[["reported.covs"]]
  out[["weights"]] <- weights
  out[["ps"]] <- ps
  out[["moderator"]] <- processed.moderator
  out[["prop.subgroup"]] <- s_min
  out[["call"]] <- match.call()

  out <- clear_null(out)

  class(out) <- c("weightit.sbps", "weightit")

  out
}

#' @exportS3Method print weightit.sbps
print.weightit.sbps <- function(x, ...) {
  treat.type <- get_treat_type(x[["treat"]])
  trim <- attr(x[["weights"]], "trim")

  cat("A weightit.sbps object\n")
  if (is_not_null(x[["method"]])) cat(paste0(" - method: \"", attr(x[["method"]], "name"), "\" (", method.to.phrase(x[["method"]]), ")\n"))
  cat(paste0(" - number of obs.: ", length(x[["weights"]]), "\n"))
  cat(paste0(" - sampling weights: ", if (is_null(x[["s.weights"]]) || all_the_same(x[["s.weights"]])) "none" else "present", "\n"))
  cat(paste0(" - treatment: ", ifelse(treat.type == "continuous", "continuous", paste0(nunique(x[["treat"]]), "-category", ifelse(treat.type == "multinomial", paste0(" (", paste(levels(x[["treat"]]), collapse = ", "), ")"), ""))), "\n"))
  if (is_not_null(x[["estimand"]])) cat(paste0(" - estimand: ", x[["estimand"]], ifelse(is_not_null(x[["focal"]]), paste0(" (focal: ", x[["focal"]], ")"), ""), "\n"))
  if (is_not_null(x[["covs"]])) cat(paste0(" - covariates: ", ifelse(length(names(x[["covs"]])) > 60, "too many to name", paste(names(x[["covs"]]), collapse = ", ")), "\n"))
  if (is_not_null(x[["by"]])) {
    cat(paste0(" - by: ", paste(names(x[["by"]]), collapse = ", "), "\n"))
  }
  if (is_not_null(x[["moderator"]])) {
    nsubgroups <- nlevels(attr(x[["moderator"]], "by.factor"))
    cat(paste0(" - moderator: ", paste(names(x[["moderator"]]), collapse = ", ")," (", nsubgroups, " subgroups)\n"))
  }
  if (is_not_null(trim)) {
    if (trim < 1) {
      if (attr(x[["weights"]], "trim.lower")) trim <- c(1 - trim, trim)
      cat(paste(" - weights trimmed at", word_list(paste0(round(100*trim, 2), "%")), "\n"))
    }
    else {
      if (attr(x[["weights"]], "trim.lower")) t.b <- "top and bottom" else t.b <- "top"
      cat(paste(" - weights trimmed at the", t.b, trim, "\n"))
    }
  }
  invisible(x)
}

#' @exportS3Method summary weightit.sbps
summary.weightit.sbps <- function(object, top = 5, ignore.s.weights = FALSE, ...) {

  if (ignore.s.weights || is_null(object$s.weights)) sw_ <- rep(1, length(object$weights))
  else sw_ <- object$s.weights
  w_ <- object$weights*sw_
  t_ <- object$treat
  mod <- object$moderator
  mod_factor <- attr(mod, "by.factor")
  mod_levels <- levels(mod_factor)
  treat.type <- get_treat_type(object[["treat"]])

  out.list <- lapply(mod_levels, function(i) {
    outnames <- c("weight.range", "weight.top","weight.ratio",
                  "coef.of.var",
                  "effective.sample.size")
    out <- make_list(outnames)

    in.subgroup <- mod_factor == i
    w <- w_[in.subgroup]
    sw <- sw_[in.subgroup]
    t <- t_[in.subgroup]
    if (treat.type == "continuous") {
      out$weight.range <- list(all = c(min(w[w > 0]),
                                       max(w[w > 0])))
      out$weight.ratio <- c(all = out$weight.range[["all"]][2]/out$weight.range[["all"]][1])
      top.weights <- sort(w, decreasing = TRUE)[seq_len(top)]
      out$weight.top <- list(all = sort(setNames(top.weights, which(w %in% top.weights)[seq_len(top)])))
      out$coef.of.var <- c(all = sd(w)/mean_fast(w))

      nn <- make_df("Total", c("Unweighted", "Weighted"))
      nn["Unweighted", ] <- ESS(sw)
      nn["Weighted", ] <- ESS(w)

    }
    else if (treat.type == "binary") {
      top0 <- c(treated = min(top, sum(t == 1)),
                control = min(top, sum(t == 0)))
      out$weight.range <- list(treated = c(min(w[w > 0 & t == 1]),
                                           max(w[w > 0 & t == 1])),
                               control = c(min(w[w > 0 & t == 0]),
                                           max(w[w > 0 & t == 0])))
      out$weight.ratio <- c(treated = out$weight.range$treated[2]/out$weight.range$treated[1],
                            control = out$weight.range$control[2]/out$weight.range$control[1],
                            overall = max(unlist(out$weight.range)/min(unlist(out$weight.range))))
      top.weights <- list(treated = sort(w[t == 1], decreasing = TRUE)[seq_len(top0["treated"])],
                          control = sort(w[t == 0], decreasing = TRUE)[seq_len(top0["control"])])
      out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w %in% top.weights[[x]] & t == {if (x == "control") 0 else 1})[seq_len(top0[x])]))),
                                 names(top.weights))

      out$coef.of.var <- c(treated = sd(w[t==1])/mean_fast(w[t==1]),
                           control = sd(w[t==0])/mean_fast(w[t==0]),
                           overall = sd(w)/mean_fast(w))

      #dc <- weightit$discarded

      nn <- make_df(c("Control", "Treated"), c("Unweighted", "Weighted"))
      nn["Unweighted", ] <- c(ESS(sw[t==0]),
                              ESS(sw[t==1]))
      nn["Weighted", ] <- c(ESS(w[t==0]),
                            ESS(w[t==1]))
    }
    else if (treat.type == "multinomial") {
      out$weight.range <- setNames(lapply(levels(t), function(x) c(min(w[w > 0 & t == x]),
                                                                   max(w[w > 0 & t == x]))),
                                   levels(t))
      out$weight.ratio <- setNames(c(vapply(out$weight.range, function(x) x[2]/x[1], numeric(1L)),
                                     max(unlist(out$weight.range)/min(unlist(out$weight.range)))),
                                   c(levels(t), "overall"))
      top.weights <- setNames(lapply(levels(t), function(x) sort(w[t == x], decreasing = TRUE)[seq_len(top)]),
                              levels(t))
      out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w %in% top.weights[[x]] & t == x)[seq_len(top)]))),
                                 names(top.weights))
      out$coef.of.var <- c(vapply(levels(t), function(x) sd(w[t==x])/mean(w[t==x]), numeric(1L)),
                           overall = sd(w)/mean(w))

      nn <- make_df(levels(t), c("Unweighted", "Weighted"))
      for (i in levels(t)) {
        nn["Unweighted", i] <- ESS(sw[t==i])
        nn["Weighted", i] <- ESS(w[t==i])
      }
    }

    out$effective.sample.size <- nn

    if (is_not_null(object$focal)) {
      w <- w[t != object$focal]
      attr(w, "focal") <- object$focal
    }
    attr(out, "weights") <- w

    out
  })

  attr(out.list, "prop.subgroup") <- matrix(c(1 - object$prop.subgroup,
                                              object$prop.subgroup),
                                            nrow = 2, byrow = TRUE,
                                            dimnames = list(c("Overall", "Subgroup"),
                                                            names(object$prop.subgroup)))
  names(out.list) <- mod_levels
  class(out.list) <- "summary.weightit.sbps"
  out.list
}

#' @exportS3Method print summary.weightit.sbps
print.summary.weightit.sbps <- function(x, ...) {
  cat("Summary of weights:\n")
  cat("\n - Overall vs. subgroup proportion contribution:\n")
  print.data.frame(round_df_char(attr(x, "prop.subgroup"), 2))

  for (g in seq_along(x)) {
    cat(paste("\n - - - - - - - Subgroup", names(x)[g], "- - - - - - -\n"))
    top <- max(lengths(x[[g]]$weight.top))
    cat("- Weight ranges:\n")
    print.data.frame(round_df_char(text_box_plot(x[[g]]$weight.range, 28), 4), ...)
    df <- setNames(data.frame(do.call("c", lapply(names(x[[g]]$weight.top), function(j) c(" ", j))),
                              matrix(do.call("c", lapply(x[[g]]$weight.top, function(j) c(names(j), rep("", top - length(j)), round(j, 4), rep("", top - length(j))))),
                                     byrow = TRUE, nrow = 2*length(x[[g]]$weight.top))),
                   rep("", 1 + top))
    cat(paste("\n- Units with", top, "greatest weights by group:\n"))
    print.data.frame(df, row.names = FALSE)
    cat("\n")
    print.data.frame(round_df_char(as.data.frame(matrix(c(x[[g]]$weight.ratio, x[[g]]$coef.of.var), ncol = 2,
                                                        dimnames = list(names(x[[g]]$weight.ratio),
                                                                        c("Ratio", "Coef of Var")))), 4))
    cat("\n- Effective Sample Sizes:\n")
    print.data.frame(round_df_char(x[[g]]$effective.sample.size, 3))
  }
  invisible(x)
}
