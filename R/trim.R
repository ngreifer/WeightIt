trim <- function(w, ...) {
  UseMethod("trim")
}

trim.weightit <- function(w, at = 0, lower = FALSE, ...) {

  if (is_null(at)) return(w)

  chk::chk_number(at)
  chk::chk_gte(at, 0)

  if (check_if_zero(at)) return(w)

  chk::chk_flag(lower)

  w[["weights"]] <- trim_weights(w[["weights"]],
                                 at = at,
                                 treat = w[["treat"]],
                                 groups.not.to.trim = w[["focal"]],
                                 lower = lower)
  w
}

trim.numeric <- function(w, at = 0, lower = FALSE, treat = NULL, ...) {
  if (is_not_null(treat)) {
    if (!has_treat_type(treat)) {
      treat <- assign_treat_type(treat)
    }
    treat.type <- get_treat_type(treat)
  }

  groups.not.to.trim <- NULL
  if (is_not_null(treat) && treat.type != "continuous") {
    w_all_same <- tapply(w, treat, all_the_same)
    if (all(w_all_same)) {
      .wrn("weights are all the same in each treatment group and will not be trimmed")
      return(w)
    }

    if (any(w_all_same)) {
      groups.not.to.trim <- names(w_all_same)[w_all_same]
    }
  }
  else if (all_the_same(w)) {
    .wrn("weights are all the same and will not be trimmed")
    return(w)
  }


  if (is_null(at)) return(w)

  chk::chk_number(at)
  chk::chk_gte(at, 0)

  if (check_if_zero(at)) return(w)

  chk::chk_flag(lower)

  trim_weights(w, at = at,
               treat = treat,
               groups.not.to.trim = groups.not.to.trim,
               lower = lower)
}

trim_weights <- function(weights, at = 0, treat = NULL, groups.not.to.trim = NULL, lower = FALSE) {

  if (at < 1) {
    at <- max(at, 1 - at)

    trim.q <- {
      if (lower) c(1 - at, at)
      else c(0, at)
    }

    if (is_not_null(groups.not.to.trim)) {
      to.be.trimmed <- treat %nin% groups.not.to.trim
      trim.w <- quantile(weights[to.be.trimmed], probs = trim.q, type = 3)
      weights[to.be.trimmed & weights < trim.w[1]] <- trim.w[1]
      weights[to.be.trimmed & weights > trim.w[2]] <- trim.w[2]
      .msg(sprintf("Trimming weights where treat is not %s to %s",
                   word_list(groups.not.to.trim, and.or = "or"),
                   word_list(paste0(round(100*trim.q[c(lower, TRUE)], 2), "%"))))
    }
    else {
      trim.w <- quantile(weights, probs = trim.q, type = 3)
      weights[weights < trim.w[1]] <- trim.w[1]
      weights[weights > trim.w[2]] <- trim.w[2]
      # if (sum(check_if_zero(weights - 1)) > 10) {
      #   warning("Several weights are equal to 1. You should enter the treatment variable as an argument to 'treat' in trim().", call. = FALSE)
      # }
      .msg(sprintf("Trimming weights to %s",
                   word_list(paste0(round(100*trim.q[c(lower, TRUE)], 2), "%"))))
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

      trim.w <- sort(weights[to.be.trimmed], partial = trim.top)[trim.top]
      weights[to.be.trimmed & weights < trim.w[1]] <- trim.w[1]
      weights[to.be.trimmed & weights > trim.w[2]] <- trim.w[2]
      if (at == 1) {
        if (lower) weights.text <- "weights"
        else weights.text <- "weight"
      }
      else weights.text <- paste(at, "weights")

      .msg(sprintf("trimming the %s %s where treat is not %s",
                   word_list(c("top", "bottom")[c(TRUE, lower)]),
                   weights.text,
                   word_list(groups.not.to.trim, and.or = "or")))

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

      trim.w <- sort(weights)[trim.top]
      weights[weights < trim.w[1]] <- trim.w[1]
      weights[weights > trim.w[2]] <- trim.w[2]

      if (at == 1) {
        if (lower) weights.text <- "weights"
        else weights.text <- "weight"
      }
      else weights.text <- paste(at, "weights")
      .msg(sprintf("trimming the %s %s",
                   word_list(c("top", "bottom")[c(TRUE, lower)]),
                   weights.text))

    }
  }

  attr(weights, "trim") <- at
  attr(weights, "trim.lower") <- lower

  weights
}
