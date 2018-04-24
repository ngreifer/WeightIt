trim <- function(x, at = .99, ...) {
  UseMethod("trim")
}

trim.weightit <- function(x, at = .99, ...) {
  x[["weights"]] <- trim.weights(x[["weights"]],
                                     trim = at,
                                     treat = x[["treat"]],
                                     estimand = x[["estimand"]],
                                     focal = x[["focal"]],
                                     treat.type = attr(x[["treat"]], "treat.type"))
  return(x)
}
trim.numeric <- function(x, at = .99, treat = NULL, ...) {
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
                        treat.type = treat.type)
  return(w)
}

trim.weights <- function(weights, trim, treat, estimand, focal, treat.type = NULL) {
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

  if (length(trim) == 0) {
    trim <- 1
  }
  else if (length(trim) > 1 || !is.numeric(trim)) {
    warning("trim must be a single number between 0 and 1. Weights will not be trimmed.", call. = FALSE)
    trim <- 1
  }
  else  {
    if (trim <= 0 || trim >= 1) {
      warning("trim must be between 0 and 1. Weights will not be trimmed.", call. = FALSE)
      trim <- 1
    }
    else {
      trim <- max(trim, 1 - trim)
    }
  }

  if (trim != 1) {
    if (toupper(estimand) == "ATT") {
      max.trim.w <- quantile(weights[treat != focal], probs = trim)
      weights[treat != focal & weights > max.trim.w] <- max.trim.w
      message(paste0("Trimming weights where treat ≠ ", focal, " to ", round(100*trim, 2), "%."))
    }
    else {
      max.trim.w <- quantile(weights, probs = trim)
      weights[weights > max.trim.w] <- max.trim.w
      if (sum(abs(weights - 1) < sqrt(.Machine$double.eps)) > 10) {
        warning("Several weights are equal to 1. You should enter the treatment variable as an argument to treat in trim().", call. = FALSE)
      }
      message(paste0("Trimming weights to ", round(100*trim, 2), "%."))
    }
  }

  attr(weights, "trim") <- trim
  return(weights)
}
