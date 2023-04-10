as.weightit <- function(...) {
  UseMethod("as.weightit")
}
as.weightit.default <- function(weights, treat, covs = NULL, estimand = NULL, s.weights = NULL, ps = NULL, ...) {

  chk::chk_not_missing(weights, "`weights`")
  chk::chk_numeric(weights)

  chk::chk_not_missing(treat, "`treat`")
  chk::chk_atomic(treat)

  if (length(weights) != length(treat)) {
    .err("`weights` and `treat` must be the same length")
  }

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)

  if (is_not_null(covs)) {
    if (!is.data.frame(covs) && !is.matrix(covs)) {
      .err("`covs` must be a data.frame of covariates")
    }
    else if (is.matrix(covs)) {
      covs <- as.data.frame.matrix(covs)
    }

    if (nrow(covs) != length(weights)) {
      .err("`covs` and `weights` must be the same length")
    }
  }

  .chk_null_or(estimand, chk::chk_string)
  .chk_null_or(s.weights, chk::chk_numeric)

  if (is_not_null(s.weights) && length(s.weights) != length(weights)) {
    .err("`s.weights` and `weights` must be the same length")
  }

  .chk_null_or(ps, chk::chk_numeric)

  if (is_not_null(ps) && length(ps) != length(weights)) {
    .err("`ps` and `weights` must be the same length")
  }

  w.list <- list(weights = weights,
                 treat = assign.treat.type(treat),
                 covs = covs,
                 estimand = estimand,
                 s.weights = s.weights,
                 ps = ps)

  if (...length() > 0L) {
    A <- list(...)
    if (is_null(names(A)) || any(names(A) == "")) {
      .err("all arguments in `...` must be named")
    }
    w.list <- c(w.list, A)
  }

  class(w.list) <- "weightit"

  w.list
}

as.weightit.CBPS <- function(cbps, s.weights = NULL, ...) {

  treat <- model.response(model.frame(cbps$terms, cbps$data))
  treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)
  covs <- cbps$data[names(cbps$data) %in% attributes(terms(cbps))$term.labels]
  weights <- cbps$weights
  c.data <- cbps$data

  if (is_not_null(s.weights)) {
    if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
      .err("the argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights")
    }
    if (is.character(s.weights) && length(s.weights)==1) {
      if (any(names(data) == s.weights)) {
        s.weights <- data[[s.weights]]
      }
      else if (any(names(c.data) == s.weights)) {
        s.weights <- data[[s.weights]]
      }
      else {
        .err("the name supplied to `s.weights` is not the name of a variable in `data`")
      }
    }
  }
  else s.weights <- rep(1, length(treat))

  weights <- weights/s.weights #Because CBPS weights contain s.weights in them

  if (is_binary(treat)) {
    if (all_the_same(weights[treat==1 & !check_if_zero(weights)]) &&
        !all_the_same(weights[treat==0 & !check_if_zero(weights)])
    ) { #if treated weights are the same and control weights differ; ATT
      estimand <- "att"
    }
    else if (all_the_same(weights[treat==0 & !check_if_zero(weights)]) &&
             !all_the_same(weights[treat==1 & !check_if_zero(weights)])
    ) { #if control weights are the same and treated weights differ; ATC
      estimand <- "ATC"
    }
    else {
      estimand <- "ATE"
    }
  }
  else {
    estimand <- "ATE"
  }

  method <- "cbps"
  attr(method, "name") <- "cbps"

  if (treat.type != "continuous") {
    ps <- cbps$fitted.values
  }
  else ps <- NULL

  call <- cbps$call

  w.list <- list(weights = weights,
                 treat = treat,
                 covs = covs,
                 estimand = estimand,
                 method = method,
                 ps = ps,
                 s.weights = s.weights,
                 call = call)
  class(w.list) <- "weightit"

  w.list
}

as.weightitMSM <- function(...) {
  UseMethod("as.weightitMSM")
}
as.weightitMSM.default <- function(weights, treat.list, covs.list = NULL, estimand = NULL, s.weights = NULL, ps.list = NULL, ...) {

  chk::chk_not_missing(weights, "`weights`")
  chk::chk_numeric(weights)

  chk::chk_not_missing(treat.list, "`treat.list`")
  chk::chk_list(treat.list)

  if (!all(vapply(treat.list, is.atomic, logical(1L)))) {
    .err("`treat.list` must be a list of atomic vectors (i.e., numeric, logical, or character) or factors")
  }

  if (!all_the_same(lengths(treat.list))) {
    .err("each component of `treat.list` must have the same length")
  }

  if (length(weights) != length(treat.list[[1]])) {
    .err("`weights` and each component of `treat.list` must be the same length")
  }

  for (i in seq_along(treat.list)) {
    if (!has.treat.type(treat.list[[i]])) {
      treat.list[[i]] <- assign.treat.type(treat.list[[i]])
    }
  }

  .chk_null_or(covs.list, chk::chk_list)
  if (is_not_null(covs.list)) {
    if (!all(vapply(covs.list, is.data.frame, logical(1L)))) {
      .err("`covs.list` must be a list of data.frames for each time point")
    }
    if (length(covs.list) != length(treat.list)) {
      .err("`covs.list` must have the same number of time points as `treat.list`")
    }
    if (!all_the_same(vapply(covs.list, nrow, numeric(1L)))) {
      .err("each component of `covs.list` must have the same number of rows")
    }
    if (length(weights) != nrow(covs.list[[1]])) {
      .err("`weights` and each component of `covs.list` must be the same length")
    }
  }

  .chk_null_or(estimand, chk::chk_string)
  .chk_null_or(s.weights, chk::chk_numeric)

  if (is_not_null(s.weights) && length(s.weights) != length(weights)) {
    .err("`s.weights` and `weights` must be the same length")
  }

  .chk_null_or(ps.list, chk::chk_list)
  if (is_not_null(ps.list)) {
    if (length(ps.list) != length(treat.list)) {
      .err("`ps.list` must have the same number of time points as `treat.list`")
    }
    if (!all(vapply(ps.list, is.vector, logical(1L), "numeric"))) {
      .err("`ps.list` must be a list of numeric vectors")
    }
    if (!all_the_same(lengths(ps.list))) {
      .err("each component of `ps.list` must have the same length")
    }
    if (length(weights) != length(ps.list[[1]])) {
      .err("`weights` and each component of `ps.list` must be the same length")
    }
  }

  w.list <- list(weights = weights,
                 treat.list = treat.list,
                 covs.list = covs.list,
                 estimand = estimand,
                 s.weights = s.weights,
                 ps.list = ps.list)


  if (...length() > 0L) {
    A <- list(...)
    if (is_null(names(A)) || any(names(A) == "")) {
      .err("all arguments in `...` must be named")
    }
    w.list <- c(w.list, A)
  }

  class(w.list) <- c("weightitMSM", "weightit")

  w.list
}