#' Create a `weightit` object manually
#'
#' @description
#' This function allows users to get the benefits of a `weightit`
#' object when using weights not estimated with [weightit()] or [weightitMSM()].
#' These benefits include diagnostics, plots, and direct compatibility with
#' \pkg{cobalt} for assessing balance.
#'
#' @param x required; a `numeric` vector of weights, one for each unit, or a
#'   `weightit.fit` object from [weightit.fit()].
#' @param treat a vector of treatment statuses, one for each unit. Required when
#'   `x` is a vector of weights.
#' @param covs an optional `data.frame` of covariates. For using \pkg{WeightIt}
#'   functions, this is not necessary, but for use with \pkg{cobalt} it is. Note
#'   that when using with a `weightit.fit` object, this should not be the matrix
#'   supplied to the `covs` argument of `weightit.fit()` unless there are no
#'   factor/character variables in it. Ideally this is the original, unprocessed
#'   covariate data frame with factor variables included.
#' @param estimand an optional `character` of length 1 giving the estimand. The
#'   text is not checked.
#' @param s.weights an optional `numeric` vector of sampling weights, one for
#'   each unit.
#' @param ps an optional `numeric` vector of propensity scores, one for each
#'   unit.
#' @param treat.list a list of treatment statuses at each time point.
#' @param covs.list an optional list of `data.frame`s of covariates of
#'   covariates at each time point. For using \pkg{WeightIt} functions, this is
#'   not necessary, but for use with \pkg{cobalt} it is.
#' @param ps.list an optional list of `numeric` vectors of propensity scores at
#'   each time point.
#' @param \dots additional arguments. These must be named. They will be included
#'   in the output object.
#'
#' @returns
#' An object of class `weightit` (for `as.weightit()`) or `weightitMSM` (for `as.weightitMSM()`).
#'
#' @examples
#' treat <- rbinom(500, 1, .3)
#' weights <- rchisq(500, df = 2)
#'
#' W <- as.weightit(weights, treat = treat, estimand = "ATE")
#' summary(W)
#'
#' # See ?weightit.fit for using as.weightit() with a
#' # weightit.fit object.
#'

#' @export
as.weightit <- function(x, ...) {
  UseMethod("as.weightit")
}

#' @exportS3Method as.weightit weightit.fit
#' @rdname as.weightit
as.weightit.weightit.fit <- function(x, covs = NULL, ...) {
  if (is_not_null(covs)) {
    if (is.matrix(covs)) {
      covs <- as.data.frame.matrix(covs)
    }
    else if (!is.data.frame(covs)) {
      .err("`covs` must be a data.frame of covariates")
    }

    if (nrow(covs) != length(x$weights)) {
      .err("`covs` must have as many rows as there are units in the original call to `weightit.fit()`")
    }

    x$covs <- covs
  }

  x <- clear_null(x)

  names(x)[names(x) == "weights"] <- "x"
  names(x)[names(x) == "fit.obj"] <- "obj"

  do.call("as.weightit.default", c(x, list(...)))
}

#' @exportS3Method as.weightit default
#' @rdname as.weightit
as.weightit.default <- function(x, treat, covs = NULL, estimand = NULL,
                                s.weights = NULL, ps = NULL, ...) {

  if (!is.numeric(x) || length(dim(x)) > 1L) {
    .err("`x` must be a numeric vector of weights")
  }

  chk::chk_not_missing(treat, "`treat`")
  .chk_basic_vector(treat)

  if (length(x) != length(treat)) {
    .err("`treat` and `x` must be the same length")
  }

  treat <- as.treat(treat, process = TRUE)

  if (is_not_null(covs)) {
    if (is.matrix(covs)) {
      covs <- as.data.frame.matrix(covs)
    }
    else if (!is.data.frame(covs)) {
      .err("`covs` must be a data.frame of covariates")
    }

    if (nrow(covs) != length(treat)) {
      .err("`covs` and `treat` must be the same length")
    }
  }

  chk::chk_null_or(estimand, vld = chk::vld_string)
  chk::chk_null_or(s.weights, vld = chk::vld_numeric)

  if (is_not_null(s.weights) && length(s.weights) != length(treat)) {
    .err("`s.weights` and `treat` must be the same length")
  }

  chk::chk_null_or(ps, vld = chk::vld_numeric)

  if (is_not_null(ps) && length(ps) != length(treat)) {
    .err("`ps` and `treat` must be the same length")
  }

  w.list <- list(weights = x,
                 treat = treat,
                 covs = covs,
                 estimand = estimand,
                 s.weights = s.weights,
                 ps = ps)

  if (...length() > 0L) {
    nm <- ...names()
    if (is_null(nm) || !all(nzchar(nm))) {
      .err("all arguments in `...` must be named")
    }

    for (i in seq_along(nm)) {
      w.list[[nm[i]]] <- ...elt(i)
    }
  }

  class(w.list) <- "weightit"

  w.list
}

# @exportS3Method as.weightit optweight
# @rdname as.weightit
.as.weightit.optweight <- function(x, ...) {
  names(x)[names(x) == "weights"] <- "x"

  x$method <- "optweight"
  attr(x$method, "name") <- x$method
  attr(x$method, "package") <- "optweight"

  w <- do.call("as.weightit.default", x, quote = TRUE)

  if (is_null(w$info)) {
    w$info <- list(duals = w$duals)
  }
  else {
    w$info$duals <- w$duals
  }

  w$duals <- NULL

  w
}

#' @export
#' @rdname as.weightit
as.weightitMSM <- function(x, ...) {
  UseMethod("as.weightitMSM")
}

#' @exportS3Method as.weightitMSM default
#' @rdname as.weightit
as.weightitMSM.default <- function(x, treat.list, covs.list = NULL, estimand = NULL,
                                   s.weights = NULL, ps.list = NULL, ...) {

  if (!is.numeric(x) || length(dim(x)) > 1L) {
    .err("`x` must be a numeric vector of weights")
  }

  chk::chk_not_missing(treat.list, "`treat.list`")
  chk::chk_list(treat.list)

  if (!all_apply(treat.list, is.atomic) ||
      any_apply(treat.list, function(z) length(dim(z)) > 1L)) {
    .err("`treat.list` must be a list of atomic vectors (i.e., numeric, logical, or character) or factors")
  }

  if (!all_the_same(lengths(treat.list))) {
    .err("each component of `treat.list` must have the same length")
  }

  if (length(x) != length(treat.list[[1L]])) {
    .err("`x` and each component of `treat.list` must be the same length")
  }

  for (i in seq_along(treat.list)) {
    treat.list[[i]] <- as.treat(treat.list[[i]], process = TRUE)
  }

  chk::chk_null_or(covs.list, vld = chk::vld_list)
  if (is_not_null(covs.list)) {
    if (!all_apply(covs.list, is.data.frame)) {
      .err("`covs.list` must be a list of data.frames for each time point")
    }
    if (length(covs.list) != length(treat.list)) {
      .err("`covs.list` must have the same number of time points as `treat.list`")
    }
    if (!all_the_same(vapply(covs.list, nrow, numeric(1L)))) {
      .err("each component of `covs.list` must have the same number of rows")
    }
    if (length(weights) != nrow(covs.list[[1L]])) {
      .err("`x` and each component of `covs.list` must be the same length")
    }
  }

  chk::chk_null_or(estimand, vld = chk::vld_string)
  chk::chk_null_or(s.weights, vld = chk::vld_numeric)

  if (is_not_null(s.weights) && length(s.weights) != length(weights)) {
    .err("`s.weights` and `x` must be the same length")
  }

  chk::chk_null_or(ps.list, vld = chk::vld_list)
  if (is_not_null(ps.list)) {
    if (length(ps.list) != length(treat.list)) {
      .err("`ps.list` must have the same number of time points as `treat.list`")
    }

    if (!all_apply(ps.list, is.numeric) ||
        any_apply(ps.list, function(z) length(dim(z)) > 1L)) {
      .err("`ps.list` must be a list of numeric vectors")
    }

    if (!all_the_same(lengths(ps.list))) {
      .err("each component of `ps.list` must have the same length")
    }

    if (length(weights) != length(ps.list[[1L]])) {
      .err("`x` and each component of `ps.list` must be the same length")
    }
  }

  w.list <- list(weights = weights,
                 treat.list = treat.list,
                 covs.list = covs.list,
                 estimand = estimand,
                 s.weights = s.weights,
                 ps.list = ps.list)


  if (...length() > 0L) {
    nm <- ...names()
    if (is_null(nm) || !all(nzchar(nm))) {
      .err("all arguments in `...` must be named")
    }

    for (i in seq_along(nm)) {
      w.list[[nm[i]]] <- ...elt(i)
    }
  }

  class(w.list) <- c("weightitMSM", "weightit")

  w.list
}
