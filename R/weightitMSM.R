#' Generate Balancing Weights for Longitudinal Treatments
#'
#' @description `weightitMSM()` allows for the easy generation of balancing
#' weights for marginal structural models for time-varying treatments using a
#' variety of available methods for binary, continuous, and multi-category
#' treatments. Many of these methods exist in other packages, which [weightit()]
#' calls; these packages must be installed to use the desired method.
#'
#' @inheritParams weightit
#' @param formula.list a list of formulas corresponding to each time point with
#'   the time-specific treatment variable on the left hand side and
#'   pre-treatment covariates to be balanced on the right hand side. The
#'   formulas must be in temporal order, and must contain all covariates to be
#'   balanced at that time point (i.e., treatments and covariates featured in
#'   early formulas should appear in later ones). Interactions and functions of
#'   covariates are allowed.
#' @param data an optional data set in the form of a data frame that contains
#'   the variables in the formulas in `formula.list`. This must be a wide data
#'   set with exactly one row per unit.
#' @param method a string of length 1 containing the name of the method that
#'   will be used to estimate weights. See [weightit()] for allowable options.
#'   The default is `"glm"`, which estimates the weights using generalized
#'   linear models.
#' @param stabilize `logical`; whether or not to stabilize the weights.
#'   Stabilizing the weights involves fitting a model predicting treatment at
#'   each time point from treatment status at prior time points. If `TRUE`, a
#'   fully saturated model will be fit (i.e., all interactions between all
#'   treatments up to each time point), essentially using the observed treatment
#'   probabilities in the numerator (for binary and multi-category treatments).
#'   This may yield an error if some combinations are not observed. Default is
#'   `FALSE`. To manually specify stabilization model formulas, e.g., to specify
#'   non-saturated models, use `num.formula`. With many time points, saturated
#'   models may be time-consuming or impossible to fit.
#' @param num.formula optional; a one-sided formula with the stabilization
#'   factors (other than the previous treatments) on the right hand side, which
#'   adds, for each time point, the stabilization factors to a model saturated
#'   with previous treatments. See Cole & Hernán (2008) for a discussion of how
#'   to specify this model; including stabilization factors can change the
#'   estimand without proper adjustment, and should be done with caution. Can
#'   also be a list of one-sided formulas, one for each time point. Unless you
#'   know what you are doing, we recommend setting `stabilize = TRUE` and
#'   ignoring `num.formula`.
#' @param include.obj whether to include in the output a list of the fit objects
#'   created in the process of estimating the weights at each time point. For
#'   example, with `method = "glm"`, a list of the `glm` objects containing the
#'   propensity score models at each time point will be included. See the help
#'   pages for each method for information on what object will be included if
#'   `TRUE`.
#' @param is.MSM.method whether the method estimates weights for multiple time
#'   points all at once (`TRUE`) or by estimating weights at each time point and
#'   then multiplying them together (`FALSE`). This is only relevant for
#'   user-specified functions.
#' @param weightit.force several methods are not valid for estimating weights
#'   with longitudinal treatments, and will produce an error message if
#'   attempted. Set to `TRUE` to bypass this error message.
#' @param ...  other arguments for functions called by `weightit()` that control
#'   aspects of fitting that are not covered by the above arguments. See Details
#'   at [weightit()].
#'
#' @returns A `weightitMSM` object with the following elements:
#' \item{weights}{The estimated weights, one for each unit.} \item{treat.list}{A
#' list of the values of the time-varying treatment variables.}
#' \item{covs.list}{A list of the covariates used in the fitting at each time
#' point. Only includes the raw covariates, which may have been altered in the
#' fitting process.} \item{data}{The data.frame originally entered to
#' `weightitMSM()`.} \item{estimand}{"ATE", currently the only estimand for MSMs
#' with binary or multi-category treatments.} \item{method}{The weight
#' estimation method specified.} \item{ps.list}{A list of the estimated
#' propensity scores (if any) at each time point.} \item{s.weights}{The provided
#' sampling weights.} \item{by}{A data.frame containing the `by` variable when
#' specified.} \item{stabilization}{The stabilization factors, if any.}
#'
#' When `keep.mparts` is `TRUE` (the default) and the chosen method is
#' compatible with M-estimation, the components related to M-estimation for use
#' in [glm_weightit()] are stored in the `"Mparts.list"` attribute. When `by` is
#' specified, `keep.mparts` is set to `FALSE`.
#'
#' @details Currently only "wide" data sets, where each row corresponds to a
#' unit's entire variable history, are supported. You can use [reshape()] or
#' other functions to transform your data into this format; see example below.
#'
#' In general, `weightitMSM()` works by separating the estimation of weights
#' into separate procedures for each time period based on the formulas provided.
#' For each formula, `weightitMSM()` simply calls `weightit()` to that formula,
#' collects the weights for each time period, and multiplies them together to
#' arrive at longitudinal balancing weights.
#'
#' Each formula should contain all the covariates to be balanced on. For
#' example, the formula corresponding to the second time period should contain
#' all the baseline covariates, the treatment variable at the first time period,
#' and the time-varying covariates that took on values after the first treatment
#' and before the second. Currently, only wide data sets are supported, where
#' each unit is represented by exactly one row that contains the covariate and
#' treatment history encoded in separate variables.
#'
#' The `"cbps"` method, which calls `CBPS()` in \pkg{CBPS}, will yield different
#' results from `CBMSM()` in \pkg{CBPS} because `CBMSM()` takes a different
#' approach to generating weights than simply estimating several time-specific
#' models.
#'
#' @seealso [weightit()] for information on the allowable methods
#'
#' [summary.weightitMSM()] for summarizing the weights
#'
#' @references Cole, S. R., & Hernán, M. A. (2008). Constructing Inverse
#' Probability Weights for Marginal Structural Models. American Journal of
#' Epidemiology, 168(6), 656–664. \doi{10.1093/aje/kwn164}
#'
#' @examples
#'
#' library("cobalt")
#'
#' data("msmdata")
#' (W1 <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
#'                         A_2 ~ X1_1 + X2_1 +
#'                           A_1 + X1_0 + X2_0,
#'                         A_3 ~ X1_2 + X2_2 +
#'                           A_2 + X1_1 + X2_1 +
#'                           A_1 + X1_0 + X2_0),
#'                    data = msmdata,
#'                    method = "glm"))
#' summary(W1)
#' bal.tab(W1)
#'
#' #Using stabilization factors
#' W2 <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
#'                         A_2 ~ X1_1 + X2_1 +
#'                           A_1 + X1_0 + X2_0,
#'                         A_3 ~ X1_2 + X2_2 +
#'                           A_2 + X1_1 + X2_1 +
#'                           A_1 + X1_0 + X2_0),
#'                    data = msmdata,
#'                    method = "glm",
#'                    stabilize = TRUE,
#'                    num.formula = list(~ 1,
#'                                       ~ A_1,
#'                                       ~ A_1 + A_2))
#'
#' #Same as above but with fully saturated stabilization factors
#' #(i.e., making the last entry in 'num.formula' A_1*A_2)
#' W3 <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
#'                         A_2 ~ X1_1 + X2_1 +
#'                           A_1 + X1_0 + X2_0,
#'                         A_3 ~ X1_2 + X2_2 +
#'                           A_2 + X1_1 + X2_1 +
#'                           A_1 + X1_0 + X2_0),
#'                    data = msmdata,
#'                    method = "glm",
#'                    stabilize = TRUE)

#' @export
weightitMSM <- function(formula.list, data = NULL, method = "glm",
                        stabilize = FALSE, by = NULL,
                        s.weights = NULL, num.formula = NULL, moments = NULL,
                        int = FALSE, missing = NULL, verbose = FALSE,
                        include.obj = FALSE, keep.mparts = TRUE,
                        is.MSM.method, weightit.force = FALSE, ...) {

  A <- list(...)

  call <- match.call()

  ## Checks and processing ----

  #Checks

  ##Process method
  .check_acceptable_method(method, msm = TRUE, force = weightit.force)

  if (is_null(method)) {
    method <- NULL
    is.MSM.method <- TRUE
  }
  else if (is.character(method)) {
    method <- .method_to_proper_method(method)
    attr(method, "name") <- method
    if (missing(is.MSM.method)) is.MSM.method <- NULL
    is.MSM.method <- .process_MSM_method(is.MSM.method, method)
  }
  else { #function
    method.name <- paste(deparse(substitute(method)))
    .check_user_method(method)
    if (missing(is.MSM.method)) is.MSM.method <- NULL
    is.MSM.method <- .process_MSM_method(is.MSM.method, method)
    attr(method, "name") <- method.name
  }

  ##Process by
  if (is_not_null(A[["exact"]])) {
    .wrn("`by` has replaced `exact` in the `weightit()` syntax, but `exact` will always work")
    by <- A[["exact"]]
    by.arg <- "exact"
  }
  else {
    by.arg <- "by"
  }

  reported.covs.list <- covs.list <- treat.list <- w.list <- ps.list <-
    stabout <- sw.list <- Mparts.list <- stab.Mparts.list <- make_list(length(formula.list))

  if (is_null(formula.list) || !is.list(formula.list) ||
      !all_apply(formula.list, rlang::is_formula, lhs = TRUE)) {
    .err("`formula.list` must be a list of formulas")
  }

  for (i in seq_along(formula.list)) {

    #Process treat and covs from formula and data
    t.c <- get_covs_and_treat_from_formula(formula.list[[i]], data)
    reported.covs.list[[i]] <- t.c[["reported.covs"]]
    covs.list[[i]] <- t.c[["model.covs"]]
    treat.list[[i]] <- t.c[["treat"]]
    treat.name <- t.c[["treat.name"]]
    names(treat.list)[i] <- treat.name
    names(reported.covs.list)[i] <- treat.name

    if (is_null(covs.list[[i]])) {
      .err(sprintf("no covariates were specified in the %s formula", ordinal(i)))
    }

    if (is_null(treat.list[[i]])) {
      .err(sprintf("no treatment variable was specified in the %s formula", ordinal(i)))
    }

    n <- length(treat.list[[i]])

    if (nrow(covs.list[[i]]) != n) {
      .err("treatment and covariates must have the same number of units")
    }

    if (anyNA(treat.list[[i]])) {
      .err(sprintf("no missing values are allowed in the treatment variable. Missing values found in %s", treat.name))
    }

    treat.list[[i]] <- assign_treat_type(treat.list[[i]])

    if (!is.MSM.method) {
      .check_method_treat.type(method, get_treat_type(treat.list[[i]]))
    }

    #By is processed each for each time to check, but only last time is used for by.factor.
    processed.by <- .process_by(by, data = data,
                                treat = treat.list[[i]],
                                treat.name = treat.name,
                                by.arg = by.arg)
  }

  #Process missing
  missing <- {
    if (is_null(method) || !any_apply(reported.covs.list, anyNA)) ""
    else .process_missing(missing, method)
  }

  #Process s.weights
  s.weights <- .process.s.weights(s.weights, data)

  if (is_null(s.weights)) s.weights <- rep.int(1, n)
  else .check_method_s.weights(method, s.weights)

  if (is_null(method)) {
    num.formula <- NULL
    stabilize <- FALSE
  }
  else if (is_not_null(num.formula)) {
    if (!isTRUE(stabilize)) {
      .msg("setting `stabilize` to `TRUE` based on `num.formula` input")
    }
    stabilize <- TRUE
  }

  if (stabilize) {
    if (!is.function(method) && !.weightit_methods[[method]]$stabilize_ok) {
      .wrn(sprintf("`stabilize` cannot be used with %s and will be ignored",
                   .method_to_phrase(method)))
      stabilize <- FALSE
      num.formula <- NULL
    }

    if (is_not_null(num.formula)) {
      if (rlang::is_formula(num.formula)) {
        if (!rlang::is_formula(num.formula, lhs = FALSE)) {
          .err("the argument to `num.formula` must have right hand side variables but not a response variable (e.g., ~ V1 + V2)")
        }

        rhs.vars.mentioned.lang <- attr(terms(num.formula), "variables")[-1L]
        rhs.vars.mentioned <- vapply(rhs.vars.mentioned.lang, deparse1, character(1L))
        rhs.vars.failed <- vapply(rhs.vars.mentioned.lang, function(v) {
          null_or_error(try(eval(v, c(data, .GlobalEnv)), silent = TRUE))
        }, logical(1L))

        if (any(rhs.vars.failed)) {
          .err(paste0(c("All variables in `num.formula` must be variables in `data` or objects in the global environment.\nMissing variables: ",
                        word_list(rhs.vars.mentioned[rhs.vars.failed], and.or = FALSE))), tidy = FALSE)
        }
      }
      else if (is.list(num.formula)) {
        if (length(num.formula) != length(formula.list)) {
          .err("when supplied as a list, `num.formula` must have as many entries as `formula.list`")
        }

        if (!all_apply(num.formula, rlang::is_formula, lhs = FALSE)) {
          .err("`num.formula` must be a single formula with no response variable and with the stabilization factors on the right hand side or a list thereof")
        }

        rhs.vars.mentioned.lang.list <- lapply(num.formula, function(nf) attr(terms(nf), "variables")[-1L])
        rhs.vars.mentioned <- unique(unlist(lapply(rhs.vars.mentioned.lang.list,
                                                   function(r) vapply(r, deparse1, character(1L)))))
        rhs.vars.failed <- vapply(rhs.vars.mentioned, function(v) {
          null_or_error(try(eval(parse(text = v), c(data, .GlobalEnv)), silent = TRUE))
        }, logical(1L))

        if (any(rhs.vars.failed)) {
          .err(paste0(c("All variables in `num.formula` must be variables in `data` or objects in the global environment.\nMissing variables: ",
                        word_list(rhs.vars.mentioned[rhs.vars.failed], and.or = FALSE))), tidy = FALSE)
        }
      }
      else {
        .err("`num.formula` must be a single formula with no response variable and with the stabilization factors on the right hand side or a list thereof")
      }
    }
  }

  #Process moments and int
  m.i.q <- .process_moments_int_quantile(moments = moments,
                                         int = int,
                                         quantile = A[["quantile"]],
                                         method = method)

  A["s.weights"] <- list(s.weights)
  A["by.factor"] <- list(attr(processed.by, "by.factor"))
  A["method"] <- list(method)
  A[c("moments", "int", "quantile")] <- m.i.q[c("moments", "int", "quantile")]
  A["subclass"] <- list(numeric())
  A["missing"] <- list(missing)
  A["verbose"] <- list(verbose)
  A["include.obj"] <- list(include.obj)

  if (is.MSM.method) {
    #Returns weights (w)

    A["covs.list"] <- list(covs.list)
    A["treat.list"] <- list(treat.list)
    A["stabilize"] <- list(stabilize)

    obj <- do.call("weightitMSM.fit", A)

    w <- obj[["weights"]]
    stabout <- NULL
    obj.list <- obj[["fit.obj"]]
    Mparts.list <- attr(obj, "Mparts")
  }
  else {

    if (is_not_null(A[["link"]])) {
      if (length(A[["link"]]) == 1L) {
        A[["link"]] <- rep.int(A[["link"]], length(formula.list))
      }
      else if (length(A[["link"]]) != length(formula.list)) {
        .err(sprintf("the argument to `link` must have length 1 or %s", length(formula.list)))
      }
    }

    obj.list <- make_list(length(formula.list))

    A["estimand"] <- list("ATE")
    A["focal"] <- list(character())
    A["stabilize"] <- list(FALSE)
    A["ps"] <- list(numeric())
    A["is.MSM.method"] <- list(FALSE)

    for (i in seq_along(formula.list)) {
      A_i <- A
      if (length(A[["link"]]) == length(formula.list)) {
        A_i["link"] <- list(A[["link"]][[i]])
      }

      A_i["covs"] <- list(covs.list[[i]])
      A_i["treat"] <- list(treat.list[[i]])
      A_i["treat.type"] <- list(get_treat_type(treat.list[[i]]))
      A_i[".data"] <- list(data)
      A_i[".covs"] <- list(reported.covs.list[[i]])

      ## Running models ----

      #Returns weights (w) and propensty score (ps)
      obj <- do.call("weightit.fit", A_i)

      w.list[i] <- list(obj[["weights"]])
      ps.list[i] <- list(obj[["ps"]])
      obj.list[i] <- list(obj[["fit.obj"]])
      Mparts.list[i] <- list(attr(obj, "Mparts"))

      if (stabilize) {
        #Process stabilization formulas and get stab weights
        if (rlang::is_formula(num.formula)) {
          if (i == 1L) {
            stab.f <- update.formula(as.formula(paste(names(treat.list)[i], "~ 1")),
                                     as.formula(paste(paste(num.formula, collapse = ""), "+ .")))
          }
          else {
            stab.f <- update.formula(as.formula(paste(names(treat.list)[i], "~",
                                                      paste(names(treat.list)[seq_along(names(treat.list)) < i],
                                                            collapse = " * "))),
                                     as.formula(paste(num.formula, "+ .")))
          }
        }
        else if (is.list(num.formula)) {
          stab.f <- update.formula(as.formula(paste(names(treat.list)[i], "~ 1")),
                                   as.formula(paste(paste(num.formula[[i]], collapse = ""), "+ .")))
        }
        else {
          if (i == 1L) {
            stab.f <- as.formula(paste(names(treat.list)[i], "~ 1"))
          }
          else {
            stab.f <- as.formula(paste(names(treat.list)[i], "~",
                                       paste(names(treat.list)[seq_along(names(treat.list)) < i],
                                             collapse = " * ")))
          }
        }

        stab.t.c_i <- get_covs_and_treat_from_formula(stab.f, data)

        A_i["covs"] <- list(stab.t.c_i[["model.covs"]])
        A_i["method"] <- list("glm")
        A_i["moments"] <- list(integer())
        A_i["int"] <- list(FALSE)
        A_i["quantile"] <- list(list())

        sw_obj <- do.call("weightit.fit", A_i)

        sw.list[[i]] <- 1 / sw_obj[["weights"]]
        stabout[[i]] <- stab.f[-2L]

        stab.Mparts.list[i] <- list(attr(sw_obj, "Mparts"))

        if (is_not_null(stab.Mparts.list[[i]])) {
          #Invert wfun and compute derivative of inverted wfun
          .wfun <- stab.Mparts.list[[i]]$wfun
          stab.Mparts.list[[i]]$wfun <- Invert(.wfun)

          if (is_not_null(stab.Mparts.list[[i]]$dw_dBtreat)) {
            .dw_dBtreat <- stab.Mparts.list[[i]]$dw_dBtreat
            stab.Mparts.list[[i]]$dw_dBtreat <- function(Btreat, Xtreat, A, SW) {
              -.dw_dBtreat(Btreat, Xtreat, A, SW) / .wfun(Btreat, Xtreat, A)^2
            }
          }
        }
      }
    }

    w <- Reduce("*", w.list, init = 1)

    if (stabilize) {
      w <-  Reduce("*", sw.list, init = w)

      unique.stabout <- unique(stabout)

      if (length(unique.stabout) <= 1L) stabout <- unique.stabout
    }
    else {
      stabout <- NULL
    }

    if (include.obj) names(obj.list) <- names(treat.list)
  }

  if (is_not_null(method) && all_the_same(w)) {
    .wrn(sprintf("all weights are %s, possibly indicating an estimation failure", w[1L]))
  }

  ## Assemble output object----
  out <- list(weights = w,
              treat.list = treat.list,
              covs.list = reported.covs.list,
              estimand = "ATE",
              method = method,
              s.weights = s.weights,
              by = processed.by,
              call = call,
              formula.list = formula.list,
              stabilization = stabout,
              missing = if (nzchar(missing)) missing else NULL,
              env = parent.frame(),
              obj = obj.list
  )

  out <- clear_null(out)

  if (keep.mparts && all(lengths(Mparts.list) > 0L)) {
    attr(out, "Mparts.list") <- clear_null(c(Mparts.list, stab.Mparts.list))
  }

  class(out) <- c("weightitMSM", "weightit")

  out
}

#' @exportS3Method print weightitMSM
print.weightitMSM <- function(x, ...) {
  treat.types <- vapply(x[["treat.list"]], get_treat_type, character(1L))

  cat(sprintf("A %s object\n", italic(class(x)[1L])))

  if (is_not_null(x[["method"]])) {
    cat(sprintf(" - method: %s (%s)\n",
                add_quotes(attr(x[["method"]], "name")),
                .method_to_phrase(x[["method"]])))
  }
  else if (all_the_same(x[["weights"]])) {
    cat(" - method: no weighting\n")
  }

  cat(sprintf(" - number of obs.: %s\n",
              nobs(x)))

  cat(sprintf(" - sampling weights: %s\n",
              if (is_null(x[["s.weights"]]) || all_the_same(x[["s.weights"]])) "none" else "present"))

  cat(sprintf(" - number of time points: %s (%s)\n",
              length(x[["treat.list"]]),
              word_list(names(x[["treat.list"]]), and.or = FALSE)))

  cat(" - treatment:\n")
  for (i in seq_along(x[["treat.list"]])) {
    cat(sprintf("    + time %s: %s\n",
                i,
                switch(treat.types[i],
                       "continuous" = "continuous",
                       "multi-category" = sprintf("%s-category (%s)",
                                                  nunique(x[["treat.list"]][[i]]),
                                                  word_list(levels(x[["treat.list"]][[i]]), and.or = FALSE)),
                       "binary" = "2-category")))
  }

  if (is_not_null(x[["covs.list"]])) {
    cat(" - covariates:\n")
    for (i in seq_along(x[["covs.list"]])) {
      if (i == 1L) {
        cat(sprintf("    + baseline: %s\n",
                    if (is_null(x$covs.list[[i]])) "(none)"
                    else word_list(names(x$covs.list[[i]]), and.or = FALSE)))
      }
      else {
        cat(sprintf("    + after time %s: %s\n",
                    i - 1L,
                    word_list(names(x$covs.list[[i]]), and.or = FALSE)))
      }
    }
  }

  if (is_not_null(x[["missing"]]) && !identical(x[["missing"]], "")) {
    cat(sprintf(" - missingness method: %s\n",
                .missing_to_phrase(x[["missing"]])))
  }

  if (is_not_null(x[["by"]])) {
    cat(sprintf(" - by: %s\n", word_list(names(x[["by"]]), and.or = FALSE)))
  }

  if (is_not_null(x$stabilization)) {
    cat(" - stabilized")
    if (any_apply(x$stabilization, function(s) is_not_null(all.vars(s)))) {
      cat(paste0("; stabilization factors:\n",
                 if (length(x$stabilization) == 1L) {
                   sprintf("      %s", word_list(attr(terms(x[["stabilization"]][[1L]]), "term.labels"),
                                                 and.or = FALSE))
                 }
                 else {
                   paste(vapply(seq_along(x$stabilization), function(i) {
                     if (i == 1L) {
                       sprintf("    + baseline: %s",
                               if (is_null(attr(terms(x[["stabilization"]][[i]]), "term.labels"))) "(none)"
                               else word_list(attr(terms(x[["stabilization"]][[i]]), "term.labels"), and.or = FALSE))
                     }
                     else {
                       sprintf("    + after time %s: %s",
                               i - 1L,
                               word_list(attr(terms(x[["stabilization"]][[i]]), "term.labels"), and.or = FALSE))
                     }
                   }, character(1L)), collapse = "\n")
                 }))
    }
  }

  trim <- attr(x[["weights"]], "trim")
  if (is_not_null(trim)) {
    if (trim < 1) {
      if (attr(x[["weights"]], "trim.lower")) trim <- c(1 - trim, trim)
      cat(sprintf(" - weights trimmed at %s\n", word_list(paste0(round(100 * trim, 2L), "%"))))
    }
    else {
      cat(sprintf(" - weights trimmed at the %s %s\n",
                  if (attr(x[["weights"]], "trim.lower")) "top and bottom"
                  else "top",
                  trim))
    }
  }

  invisible(x)
}
