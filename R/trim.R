trim <- function(w, ...) {
  UseMethod("trim")
}

trim.weightit <- function(w, at = 0, lower = FALSE, ...) {
  w[["weights"]] <- trim_weights(w[["weights"]],
                                 at = at,
                                 treat = w[["treat"]],
                                 groups.not.to.trim = w[["focal"]],
                                 lower = lower)
  return(w)
}
trim.numeric <- function(w, at = 0, lower = FALSE, treat = NULL, ...) {
  if (is_not_null(treat)) {
    if (!has.treat.type(treat)) {
      treat <- assign.treat.type(treat)
    }
    treat.type <- get.treat.type(treat)
  }

  groups.not.to.trim <- NULL
  if (is_not_null(treat) && treat.type != "continuous") {
    w_all_same <- tapply(w, treat, all_the_same)
    if (all(w_all_same)) {
      warning("Weights are all the same in each treatment group and will not be trimmed.", call. = FALSE)
      at <- NULL
    }
    else if (any(w_all_same)) {
      groups.not.to.trim <- names(w_all_same)[w_all_same]
    }
  }
  else if (all_the_same(w)) {
    warning("Weights are all the same and will not be trimmed.", call. = FALSE)
    at <- NULL
  }

  w <- trim_weights(w, at = at,
                    treat = treat,
                    groups.not.to.trim = groups.not.to.trim,
                    lower = lower)
  return(w)
}

trim_weights <- function(weights, at = 0, treat = NULL, groups.not.to.trim = NULL, lower = FALSE) {

  if (is_null(at) || isTRUE(at == 0)) {
    at <- NULL
  }
  else if (length(at) > 1 || !is.numeric(at) || is.na(at) || at < 0) {
    warning("'at' must be a single positive number. Weights will not be trimmed.", call. = FALSE)
    at <- NULL
  }
  else {
    if (at < 1) {
      at <- max(at, 1 - at)

      if (lower) trim.q <- c(1 - at, at)
      else trim.q <- c(0, at)

      if (is_not_null(groups.not.to.trim)) {
        to.be.trimmed <- treat %nin% groups.not.to.trim
        trim.w <- quantile(weights[to.be.trimmed], probs = trim.q, type = 3)
        weights[to.be.trimmed & weights < trim.w[1]] <- trim.w[1]
        weights[to.be.trimmed & weights > trim.w[2]] <- trim.w[2]
        message(paste0("Trimming weights where treat is not ", word_list(groups.not.to.trim, and.or = "or"), " to ",
                       word_list(paste0(round(100*trim.q[c(lower, TRUE)], 2), "%")), "."))
      }
      else {
        trim.w <- quantile(weights, probs = trim.q, type = 3)
        weights[weights < trim.w[1]] <- trim.w[1]
        weights[weights > trim.w[2]] <- trim.w[2]
        # if (sum(check_if_zero(weights - 1)) > 10) {
        #   warning("Several weights are equal to 1. You should enter the treatment variable as an argument to 'treat' in trim().", call. = FALSE)
        # }
        message(paste0("Trimming weights to ", word_list(paste0(round(100*trim.q[c(lower, TRUE)], 2), "%")), "."))
      }
    }
    else {
      if (is_not_null(groups.not.to.trim)) {
        to.be.trimmed <- treat %nin% groups.not.to.trim
        if (at >= sum(to.be.trimmed)) {
          warning(paste0("'at' must be less than ", sum(to.be.trimmed), ", the number of units for which treat is not ",
                         word_list(groups.not.to.trim, and.or = "or"), ". Weights will not be trimmed."), call. = FALSE)
          at <- NULL
        }
        else {
          at <- as.integer(min(at, sum(to.be.trimmed) - at))

          if (lower) trim.top <- c(at + 1, sum(to.be.trimmed) - at)
          else trim.top <- c(1, sum(to.be.trimmed) - at)

          trim.w <- sort(weights[to.be.trimmed], partial = trim.top)[trim.top]
          weights[to.be.trimmed & weights < trim.w[1]] <- trim.w[1]
          weights[to.be.trimmed & weights > trim.w[2]] <- trim.w[2]
          if (at == 1) {
            if (lower) weights.text <- "weights"
            else weights.text <- "weight"
          }
          else weights.text <- paste(at, "weights")
          message(paste0("Trimming the ", word_list(c("top", "bottom")[c(TRUE, lower)]), " ", weights.text,
                         " where treat is not ", word_list(groups.not.to.trim, and.or = "or"), "."))
        }
      }
      else {
        if (at >= length(weights)) {
          warning(paste0("'at' must be less than ", length(weights), ", the number of units. Weights will not be trimmed."), call. = FALSE)
          at <- NULL
        }
        else {
          at <- as.integer(min(at, length(weights) - at))

          if (lower) trim.top <- c(at + 1, length(weights) - at)
          else trim.top <- c(1, length(weights) - at)

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
