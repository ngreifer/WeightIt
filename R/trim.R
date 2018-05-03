trim <- function(x, ...) {
  UseMethod("trim")
}

trim.weightit <- function(x, at = .99, lower = FALSE, ...) {
  x[["weights"]] <- trim.weights(x[["weights"]],
                                 at = at,
                                 treat = x[["treat"]],
                                 estimand = x[["estimand"]],
                                 focal = x[["focal"]],
                                 treat.type = attr(x[["treat"]], "treat.type"),
                                 lower = lower)
  return(x)
}
trim.numeric <- function(x, at = .99, lower = FALSE, treat = NULL, ...) {
  if (length(treat) > 0 && length(attr(treat, "treat.type")) == 0) {
    nunique.treat <- nunique(treat)
    if (nunique.treat == 2) {
      treat.type <- "binary"
    }
    else if (nunique.treat < 2) {
      stop("The treatment must have at least two unique values.", call. = FALSE)
    }
    else if (is.factor(treat) || is.character(treat)) {
      treat.type <- "multinomial"
      treat <- factor(treat)
    }
    else {
      treat.type <- "continuous"
    }
  }
  else treat.type <- "continuous"

  if (treat.type != "continuous" && length(treat) > 0) {
    nunique.w.gt.1 <- tapply(x, treat, function(y) (max(y) - min(y))^2 > .Machine$double.eps)
    if (all(nunique.w.gt.1)) {
      estimand <- "ATE"
      focal <- NULL
    }
    else if (sum(!nunique.w.gt.1) == 1) {
      estimand <- "ATT"
      focal <- names(nunique.w.gt.1)[!nunique.w.gt.1]
    }
    else {
      stop("It appears there is more than one focal group, which is a no-no.", call. = FALSE)
    }
  }
  else {
    estimand <- "ATE"
    focal <- NULL
  }

  w <- trim.weights(x, trim = at,
                    treat = treat,
                    estimand = estimand,
                    focal = focal,
                    treat.type = treat.type,
                    lower = lower)
  return(w)
}

trim.weights <- function(weights, at, treat, estimand, focal, treat.type = NULL, lower) {
  estimand <- toupper(estimand)
  if (length(estimand) != 1 ||
      !is.character(estimand) ||
      !estimand %in% c("ATT", "ATC", "ATE", "ATO", "ATM")) {
    stop("estimand must be a character vector of length 1 with an acceptable estimand value (e.g., ATT, ATC, ATE).", call. = FALSE)
  }

  f.e.r <- process.focal.and.estimand(focal, estimand, treat, treat.type)
  focal <- f.e.r[["focal"]]
  estimand <- f.e.r[["estimand"]]
  reported.estimand <- f.e.r[["reported.estimand"]]

  if (length(at) == 0) {
    at <- 1
  }
  else if (length(at) > 1 || !is.numeric(at)) {
    warning("\"at\" must be a single number between 0 and 1. Weights will not be trimmed.", call. = FALSE)
    at <- 1
  }
  else  if (at <= 0 || at >= 1) {
    warning("\"at\" must be between 0 and 1. Weights will not be trimmed.", call. = FALSE)
    at <- 1
  }
  else {
    if (lower) at <- sort(c(at, 1 - at))
    else at <- c(0, max(at, 1 - at))
  }

  if (all(at != 1)) {
    if (toupper(estimand) == "ATT") {
      trim.w <- quantile(weights[treat != focal], probs = at)
      weights[treat != focal & weights < trim.w] <- trim.w[1]
      weights[treat != focal & weights > trim.w] <- trim.w[2]
      message(paste0("Trimming weights where treat ≠ ", focal, " to ", round(100*at, 2), "%."))
    }
    else {
      trim.w <- quantile(weights, probs = at)
      weights[weights < trim.w] <- trim.w[1]
      weights[weights > trim.w] <- trim.w[2]
      if (sum(abs(weights - 1) < sqrt(.Machine$double.eps)) > 10) {
        warning("Several weights are equal to 1. You should enter the treatment variable as an argument to treat in trim().", call. = FALSE)
      }
      message(paste0("Trimming weights to ", word.list(paste0(round(100*at[c(lower, TRUE)], 2), "%")), "."))
    }
  }

  attr(weights, "trim") <- at
  return(weights)
}
