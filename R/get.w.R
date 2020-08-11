#Extract weights from output objects; duplicated from cobalt

get.w <- function(x, ...) UseMethod("get.w")
get.w.ps <- function(x, stop.method = NULL, estimand = NULL, s.weights = FALSE, ...) {
  ps <- x
  estimand <- tolower(estimand)
  if (is_not_null(stop.method)) {
    if (any(is.character(stop.method))) {
      rule1 <- names(ps$w)[pmatch(tolower(names(ps$w)), tolower(stop.method), 0L)]
      if (is_null(rule1)) {
        message(paste0("Warning: stop.method should be ", word_list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
        rule1 <- names(ps$w)
      }
    }
    else if (is.numeric(stop.method) && any(stop.method %in% seq_along(names(ps$w)))) {
      if (any(!stop.method %in% seq_along(names(ps$w)))) {
        message(paste0("Warning: There are ", length(names(ps$w)), " stop methods available, but you requested ",
                       word_list(stop.method[!stop.method %in% seq_along(names(ps$w))], and.or = "and"),"."))
      }
      rule1 <- names(ps$w)[stop.method %in% seq_along(names(ps$w))]
    }
    else {
      warning("stop.method should be ", word_list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
      rule1 <- names(ps$w)
    }
  }
  else {
    rule1 <- names(ps$w)
  }

  s <- names(ps$w)[match(tolower(rule1), tolower(names(ps$w)))]
  criterion <- substr(tolower(s), 1, nchar(s)-4)

  if (is_null(estimand)) estimand <- setNames(substr(toupper(s), nchar(s)-2, nchar(s)), s)
  else if (!all(toupper(estimand) %in% c("ATT", "ATE", "ATC"))) {
    stop('estimand must be "ATT", "ATE", or "ATC".', call. = FALSE)
  }
  else {
    if (length(estimand) == 1) estimand <- setNames(toupper(rep(estimand, length(s))), s)
    else if (length(estimand) >= length(s)) estimand <- setNames(toupper(estimand[seq_along(s)]), s)
    else stop("estimand must be the same length as the number of sets of weights requested.", call. = FALSE)
  }

  w <- make_df(s, nrow(ps$ps))
  for (p in s) {
    if (estimand[p] == "ATT") w[[p]] <- ps$treat + (1-ps$treat)*ps$ps[,p]/(1-ps$ps[,p])
    else if (estimand[p] == "ATE") w[[p]] <- ps$treat/ps$ps[,p] + (1-ps$treat)/(1-ps$ps[,p])
    else if (estimand[p] == "ATC") w[[p]] <- (1-ps$treat) + ps$treat*ps$ps[,p]/(1-ps$ps[,p])
    else w[[p]] <- ps$w[,p]
    if (s.weights) w[[p]] <- w[[p]] * ps$sampw
  }

  names(w) <- ifelse(toupper(substr(s, nchar(s)-2, nchar(s))) == estimand, criterion, paste0(criterion, " (", estimand, ")"))
  if (ncol(w) == 1) w <- w[[1]]

  return(w)
}
get.w.mnps <- function(x, stop.method = NULL, s.weights = FALSE, ...) {
  mnps <- x
  if (is_not_null(stop.method)) {
    if (is.character(stop.method)) {
      rule1 <- mnps$stopMethods[pmatch(tolower(stop.method), tolower(mnps$stopMethods), nomatch = 0L)]
      if (is_null(rule1)) {
        message(paste0("Warning: stop.method should be ", word_list(mnps$stopMethods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
        rule1 <- mnps$stopMethods
      }
    }
    else if (is.numeric(stop.method) && any(stop.method %in% seq_along(mnps$stopMethods))) {
      if (any(!stop.method %in% seq_along(mnps$stopMethods))) {
        message(paste0("Warning: There are ", length(mnps$stopMethods), " stop methods available, but you requested ",
                       word_list(stop.method[!stop.method %in% seq_along(mnps$stopMethods)], and.or = "and"),"."))
      }
      rule1 <- mnps$stopMethods[stop.method %in% seq_along(mnps$stopMethods)]
    }
    else {
      warning("stop.method should be ", word_list(mnps$stopMethods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
      rule1 <- mnps$stopMethods
    }
  }
  else {
    rule1 <- mnps$stopMethods
  }

  s <- paste.(mnps$stopMethods[match(tolower(rule1), tolower(mnps$stopMethods))],
              mnps$estimand)

  estimand <- mnps$estimand
  criterion <- mnps$stopMethods[match(tolower(rule1), tolower(mnps$stopMethods))]

  w <- make_df(criterion, length(mnps$treatVar))

  if (estimand == "ATT") {
    for (i in mnps$levExceptTreatATT) {
      if (length(s) > 1) {
        w[mnps$treatVar == i, criterion] <- get.w.ps(mnps$psList[[i]])[mnps$psList[[i]]$treat == FALSE, criterion]
      }
      else {
        w[mnps$treatVar == i, criterion] <- get.w.ps(mnps$psList[[i]])[mnps$psList[[i]]$treat == FALSE]
      }
    }
  }
  else if (estimand == "ATE") {
    for (i in mnps$treatLev) {
      if (length(s) > 1) {
        w[mnps$treatVar == i, criterion] <- get.w.ps(mnps$psList[[i]])[mnps$psList[[i]]$treat == TRUE, criterion]
      }
      else {
        w[mnps$treatVar == i, criterion] <- get.w.ps(mnps$psList[[i]])[mnps$psList[[i]]$treat == TRUE]
      }
    }
  }

  if (s.weights) {
    w <- w * mnps$sampw
  }

  names(w) <- ifelse(toupper(substr(s, nchar(s)-2, nchar(s))) == estimand, criterion, paste0(criterion, " (", estimand, ")"))

  if (ncol(w) == 1) w <- w[[1]]

  return(w)
}
get.w.ps.cont <- function(x, stop.method = NULL, s.weights = FALSE, ...) {
  ps.cont <- x
  if (is_not_null(stop.method)) {
    if (any(is.character(stop.method))) {
      rule1 <- names(ps.cont$w)[pmatch(tolower(names(ps.cont$w)), tolower(stop.method), 0L)]
      if (is_null(rule1)) {
        message(paste0("Warning: stop.method should be ", word_list(names(ps.cont$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
        rule1 <- names(ps.cont$w)
      }

    }
    else if (is.numeric(stop.method) && any(stop.method %in% seq_along(names(ps.cont$w)))) {
      if (any(!stop.method %in% seq_along(names(ps.cont$w)))) {
        message(paste0("Warning: There are ", length(names(ps.cont$w)), " stop methods available, but you requested ",
                       word_list(stop.method[!stop.method %in% seq_along(names(ps.cont$w))], and.or = "and"),"."))
      }
      rule1 <- names(ps.cont$w)[stop.method %in% seq_along(names(ps.cont$w))]
    }
    else {
      warning("stop.method should be ", word_list(names(ps.cont$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
      rule1 <- names(ps.cont$w)
    }
  }
  else {
    rule1 <- names(ps.cont$w)
  }

  s <- names(ps.cont$w)[match(tolower(rule1), tolower(names(ps.cont$w)))]

  w <- make_df(s, nrow(ps.cont$w))

  for (p in s) {
    w[[p]] <- ps.cont$w[[p]]
    if (!s.weights && is_not_null(ps.cont$sampw)) w[[p]] <- w[[p]] / ps.cont$sampw
  }

  if (ncol(w) == 1) w <- w[[1]]

  return(w)
}
get.w.iptw <- function(x, stop.method = NULL, s.weights = FALSE, ...) {
  iptw <- x
  if (is_not_null(stop.method)) {
    if (any(is.character(stop.method))) {
      rule1 <- names(iptw$psList[[1]]$ps)[pmatch(tolower(names(iptw$psList[[1]]$ps)), tolower(stop.method), 0L)]
      if (is_null(rule1)) {
        message(paste0("Warning: stop.method should be ", word_list(names(iptw$psList[[1]]$ps), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
        rule1 <- names(iptw$psList[[1]]$ps)
      }
    }
    else if (is.numeric(stop.method) && any(stop.method %in% seq_along(names(iptw$psList[[1]]$ps)))) {
      if (any(!stop.method %in% seq_along(names(iptw$psList[[1]]$ps)))) {
        message(paste0("Warning: There are ", length(names(iptw$psList[[1]]$ps)), " stop methods available, but you requested ",
                       word_list(stop.method[!stop.method %in% seq_along(names(iptw$psList[[1]]$ps))], and.or = "and"),"."))
      }
      rule1 <- names(iptw$psList[[1]]$ps)[stop.method %in% seq_along(names(iptw$psList[[1]]$ps))]
    }
    else {
      warning("stop.method should be ", word_list(names(iptw$psList[[1]]$ps), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
      rule1 <- names(iptw$psList[[1]]$ps)
    }
  }
  else {
    rule1 <- names(iptw$psList[[1]]$ps)
  }

  w <- make_df(rule1, nrow(iptw$psList[[1]]$ps))

  for (i in rule1) {
    w[i] <- Reduce("*", lapply(iptw$psList, function(x) get.w.ps(x, stop.method = i)))
  }

  if (s.weights) {
    w <- w * iptw$psList[[1]]$sampw
  }

  return(w)
}
get.w.CBPS <- function(x, estimand = NULL, ...) {
  c <- x
  A <- list(...)
  if (is_null(A$use.weights)) use.weights <- TRUE
  else use.weights <- A$use.weights

  estimand <- tolower(estimand)

  if ("CBPSContinuous" %in% class(c) || is.factor(c$y)) { #continuous
    return(c$weights)
  }
  else {
    if (!use.weights) {
      ps <- c$fitted.values
      t <- c$y
      if (is_null(estimand)) {
        if (nunique.gt(c$weights[t == 1], 1)) {
          estimand <- "ate"
        }
        else estimand <- "att"
      }

      estimand <- match_arg(tolower(estimand), c("att", "atc", "ate"))
      if (estimand == "att") {
        return(ifelse(t == 1, 1, ps/(1-ps)))
      }
      if (estimand == "atc") {
        return(ifelse(t == 1, (1-ps)/ps, 1))
      }
      else if (estimand == "ate") {
        return(ifelse(t == 1, 1/ps, 1/(1-ps)))
      }
    }
    else {
      return(c$weights)
    }

  }
}
get.w.npCBPS <- function(x, ...) {
  return(x$weights)
}
get.w.CBMSM <- function(x, ...) {
  return(x$weights)
}
get.w.ebalance <- function(x, treat, ...) {
  if (missing(treat)) stop("treat must be specified.", call. = FALSE)

  weights <- rep(1, length(treat))

  if (length(x$w) != sum(treat == 0)) {
    stop("There are more control units in treat than weights in the ebalance object.", call. = FALSE)
  }
  weights[treat == 0] <- x$w
  return(weights)
}
get.w.ebalance.trim <- get.w.ebalance
get.w.weightit <- function(x, s.weights = FALSE, ...) {
  if (s.weights) return(x$weights * x$s.weights)
  else return(x$weights)
}
