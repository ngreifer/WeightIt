trim <- function(x, ...) {
  UseMethod("trim")
}

trim.weightit <- function(x, at = .99, lower = FALSE, ...) {
  x[["weights"]] <- trim_weights(x[["weights"]],
                                 at = at,
                                 treat = x[["treat"]],
                                 estimand = x[["estimand"]],
                                 focal = x[["focal"]],
                                 treat.type = get.treat.type(x[["treat"]]),
                                 lower = lower)
  return(x)
}
trim.numeric <- function(x, at = .99, lower = FALSE, treat = NULL, ...) {
  if (is_not_null(treat)) {
    if (!has.treat.type(treat)) {
      treat <- assign.treat.type(treat)
    }
    treat.type <- get.treat.type(treat)
  }
  else treat.type <- "continuous"

  if (treat.type != "continuous" && is_not_null(treat)) {
    w_all_same <- tapply(x, treat, all_the_same)
    if (!any(w_all_same)) {
      estimand <- "ATE"
      focal <- NULL
    }
    else if (sum(w_all_same) == 1) {
      estimand <- "ATT"
      focal <- names(w_all_same)[w_all_same]
    }
    else {
      stop("It appears there is more than one focal group, which is a no-no.", call. = FALSE)
    }
  }
  else {
    estimand <- "ATE"
    focal <- NULL
  }

  w <- trim_weights(x, at = at,
                    treat = treat,
                    estimand = estimand,
                    focal = focal,
                    treat.type = treat.type,
                    lower = lower)
  return(w)
}

trim_weights <- function(weights, at, treat, estimand, focal, treat.type = NULL, lower) {
  estimand <- toupper(estimand)
  if (treat.type != "continuous" && (length(estimand) != 1 ||
      !is.character(estimand) ||
      !estimand %in% c("ATT", "ATC", "ATE", "ATO", "ATM"))) {
    stop("estimand must be a character vector of length 1 with an acceptable estimand value (e.g., ATT, ATC, ATE).", call. = FALSE)
  }

  f.e.r <- process.focal.and.estimand(focal, estimand, treat, treat.type)
  focal <- f.e.r[["focal"]]
  estimand <- f.e.r[["estimand"]]
  #reported.estimand <- f.e.r[["reported.estimand"]]

  if (is_null(at) || isTRUE(at == 0)) {
    at <- NULL
  }
  else if (length(at) > 1 || !is.numeric(at) || at < 0) {
    warning("\"at\" must be a single positive number. Weights will not be trimmed.", call. = FALSE)
    at <- NULL
  }
  else {
    if (at < 1) {
      at <- max(at, 1 - at)

      if (lower) trim.q <- c(1 - at, at)
      else trim.q <- c(0, at)

      if (treat.type != "continuous" && toupper(estimand) == "ATT") {
        trim.w <- quantile(weights[treat != focal], probs = trim.q, type = 3)
        weights[treat != focal & weights < trim.w[1]] <- trim.w[1]
        weights[treat != focal & weights > trim.w[2]] <- trim.w[2]
        message(paste0("Trimming weights where treat is not ", focal, " to ", word_list(paste0(round(100*trim.q[c(lower, TRUE)], 2), "%")), "."))
      }
      else {
        trim.w <- quantile(weights, probs = trim.q, type = 3)
        weights[weights < trim.w[1]] <- trim.w[1]
        weights[weights > trim.w[2]] <- trim.w[2]
        if (sum(check_if_zero(weights - 1)) > 10) {
          warning("Several weights are equal to 1. You should enter the treatment variable as an argument to treat in trim().", call. = FALSE)
        }
        message(paste0("Trimming weights to ", word_list(paste0(round(100*trim.q[c(lower, TRUE)], 2), "%")), "."))
      }
    }
    else {
      if (treat.type != "continuous" && toupper(estimand) == "ATT") {
        if (at >= sum(treat != focal)) {
          warning(paste0("'at' must be less than ", sum(treat != focal), ", the number of units not in the focal treatment. Weights will not be trimmed."), call. = FALSE)
          at <- NULL
        }
        else {
          at <- as.integer(min(at, sum(treat != focal) - at))

          if (lower) trim.top <- c(at + 1, sum(treat != focal) - at)
          else trim.top <- c(1, sum(treat != focal) - at)

          trim.w <- sort(weights[treat != focal])[trim.top]
          weights[treat != focal & weights < trim.w[1]] <- trim.w[1]
          weights[treat != focal & weights > trim.w[2]] <- trim.w[2]
          if (at == 1) {
            if (lower) weights.text <- "weights"
            else weights.text <- "weight"
          }
          else weights.text <- paste(at, "weights")
          message(paste0("Trimming the ", word_list(c("top", "bottom")[c(TRUE, lower)]), " ", weights.text, " where treat \u2260 ", focal, "."))
        }
      }
      else {
        if (at >= length(treat)) {
          warning(paste0("'at' must be less than ", length(treat), ", the number of units. Weights will not be trimmed."), call. = FALSE)
          at <- NULL
        }
        else {
          at <- as.integer(min(at, length(treat) - at))

          if (lower) trim.top <- c(at + 1, length(treat) - at)
          else trim.top <- c(1, length(treat) - at)

          trim.w <- sort(weights)[trim.top]
          weights[weights < trim.w[1]] <- trim.w[1]
          weights[weights > trim.w[2]] <- trim.w[2]
          if (at == 1) {
            if (lower) weights.text <- "weights"
            else weights.text <- "weight"
          }
          else weights.text <- paste(at, "weights")
          message(paste0("Trimming the ", word_list(c("top", "bottom")[c(TRUE, lower)]), " ", weights.text, "."))
        }
      }
    }
  }

  attr(weights, "trim") <- at
  if (is_not_null(at)) attr(weights, "trim.lower") <- lower
  return(weights)
}
