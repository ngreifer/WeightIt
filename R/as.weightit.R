as.weightit <- function(...) {
  UseMethod("as.weightit")
}
as.weightit.default <- function(weights, treat, covs = NULL, estimand = NULL, s.weights = NULL, ps = NULL, ...) {

  if (missing(weights)) stop("'weights' must be supplied to as.weightit().", call. = FALSE)
  if (missing(treat)) stop("'treat' must be supplied to as.weightit().", call. = FALSE)

  if (!is.vector(weights, "numeric")) stop("'weights' must be a numeric vector.", call. = FALSE)

  if (!is.atomic(treat) && !is.factor(treat)) stop("'treat' must be an atomic vector (i.e., numeric, logical, or character) or a factor.", call. = FALSE)
  if (length(weights) != length(treat)) stop("'weights' and 'treat' must be the same length.", call. = FALSE)
  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)

  if (is_not_null(covs)) {
    if (!is.data.frame(covs) && !is.matrix(covs)) stop("'covs' must be a data.frame of covariates.", call. = FALSE)
    else if (is.matrix(covs)) covs <- as.data.frame.matrix(covs)

    if (nrow(covs) != length(weights)) stop("'covs' and 'weights' must be the same length.", call. = FALSE)
  }
  if (is_not_null(estimand)) {
    if (!is.vector(estimand, "character") || length(estimand) != 1L) stop("'estimand' must be a character vector of length 1.", call. = FALSE)
  }

  if (is_not_null(s.weights)) {
    if (!is.vector(s.weights, "numeric")) stop("'s.weights' must be a numeric vector.", call. = FALSE)
    if (length(s.weights) != length(weights)) stop("'s.weights' and 'weights' must be the same length.", call. = FALSE)
  }
  if (is_not_null(ps)) {
    if (!is.vector(ps, "numeric")) stop("'ps' must be a numeric vector.", call. = FALSE)
    if (length(ps) != length(ps)) stop("'ps' and 'weights' must be the same length.", call. = FALSE)
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
      stop("All arguments in ... must be named.")
    }
    w.list <- c(w.list, A)
  }

  class(w.list) <- "weightit"

  return(w.list)
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
      stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
    }
    if (is.character(s.weights) && length(s.weights)==1) {
      if (any(names(data) == s.weights)) {
        s.weights <- data[[s.weights]]
      }
      else if (any(names(c.data) == s.weights)) {
        s.weights <- data[[s.weights]]
      }
      else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
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
  return(w.list)
}

as.weightitMSM <- function(...) {
  UseMethod("as.weightitMSM")
}
as.weightitMSM.default <- function(weights, treat.list, covs.list = NULL, estimand = NULL, s.weights = NULL, ps.list = NULL, ...) {

  if (missing(weights)) stop("weights must be supplied to as.weightitMSM.", call. = FALSE)
  if (missing(treat.list)) stop("treat.list must be supplied to as.weightitMSM.", call. = FALSE)

  if (!is.vector(weights, "numeric")) stop("weights must be a numeric vector.", call. = FALSE)

  if (!is.vector(treat.list, "list")) stop("treat.list must be a list of treatment statuses for each individual at each time point.", call. = FALSE)
  if (any(vapply(treat.list, function(x) !is.atomic(x) && !is.factor(x), logical(1L)))) stop("treat.list must be a list of atomic vectors (i.e., numeric, logical, or character) or factors.", call. = FALSE)
  if (!all_the_same(lengths(treat.list))) stop("Each component of treat.list must have the same length.", call. = FALSE)
  if (length(weights) != length(treat.list[[1]])) stop("weights and each component of treat.list must be the same length.", call. = FALSE)
  treat.list <- sapply(treat.list, function(t) if (!has.treat.type(t)) assign.treat.type(t) else t, simplify = FALSE)

  if (is_not_null(covs.list)) {
    if (!is.vector(covs.list, "list") || any(vapply(covs.list, function(x) !is.data.frame(x), logical(1L)))) stop("covs.list must be a list of data.frames for each time point.", call. = FALSE)
    if (length(covs.list) != length(treat.list)) stop("covs.list must have the same number of time points as treat.list.", call. = FALSE)
    if (!all_the_same(vapply(covs.list, nrow, numeric(1L)))) stop("Each component of covs.list must have the same number of rows.", call. = FALSE)
    if (length(weights) != nrow(covs.list[[1]])) stop("weights and each component of covs.list must be the same length.", call. = FALSE)
  }
  if (is_not_null(estimand)) {
    if (!is.vector(estimand, "character") || length(estimand) != 1L) stop("estimand must be a character vector of length 1.", call. = FALSE)
  }
  if (is_not_null(s.weights)) {
    if (!is.vector(s.weights, "numeric")) stop("s.weights must be a numeric vector.", call. = FALSE)
    if (length(s.weights) != length(weights)) stop("s.weights and weights must be the same length.", call. = FALSE)
  }
  if (is_not_null(ps.list)) {
    if (!is.vector(ps.list, "list")) stop("ps.list must be a list of propensity scores for each individual at each time point.", call. = FALSE)
    if (length(ps.list) != length(treat.list)) stop("ps.list must have the same number of time points as treat.list.", call. = FALSE)
    if (any(vapply(ps.list, function(x) !is.vector(x, "numeric"), logical(1L)))) stop("ps.list must be a list of numeric vectors.", call. = FALSE)
    if (!all_the_same(lengths(ps.list))) stop("Each component of ps.list must have the same length.", call. = FALSE)
    if (length(weights) != length(ps.list[[1]])) stop("weights and each component of ps.list must be the same length.", call. = FALSE)
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
      stop("All arguments in ... must be named.")
    }
    w.list <- c(w.list, A)
  }

  class(w.list) <- c("weightitMSM", "weightit")

  return(w.list)
}