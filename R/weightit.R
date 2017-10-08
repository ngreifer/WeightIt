weightit <- function(formula, data, method, estimand = "ATE", stabilize = FALSE, focal = NULL,
                     exact = NULL, s.weights = NULL, ps = NULL, trim = c(0, Inf), trim.q = c(0, 1),
                     truncate = c(0, Inf), truncate.q = c(0, 1), verbose = FALSE, ...) {

  #Checks
  if (length(ps) == 0 && (length(formula) == 0 || length(class(formula)) == 0)) {
    stop("formula must be a formula relating treatment to covariates.", call. = FALSE)
  }
  if (missing(data)) {
    stop("Data must be specified.", call. = FALSE)}
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)}

  if (length(ps) > 0) {
    if (!(is.character(ps) && length(ps) == 1) && !is.numeric(ps)) {
      stop("The argument to ps must be a vector or data frame of propensity scores or the (quoted) names of variables in data that contain propensity scores.", call. = FALSE)
    }
    if (is.character(ps) && length(ps)==1) {
      if (ps %in% names(data)) {
        ps <- data[, ps]
      }
      else stop("The name supplied to ps is not the name of a variable in data.", call. = FALSE)
    }
  }



  ##Process method
  bad.method <- FALSE
  acceptable.methods <- c("ps",
                          "gbm", "twang", "gbr",
                          "cbps",
                          "ebal", "entropy", "ebalance",
                          "sbw",
                          "ebcw", "ate")
  if (missing(method) || length(ps) > 0) method <- "ps"
  else if (!is.character(method)) bad.method <- TRUE
  else if (length(method) != 1) bad.method <- TRUE
  else if (!tolower(method) %in% acceptable.methods) bad.method <- TRUE

  if (bad.method) stop("method must be a string of length 1 containing the name of an acceptable weighting method.", call. = FALSE)
  else method <- method.to.proper.method(tolower(method))

  #Process treat and covs from formula and data
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.na(match(rownames(attr(tt, "factors"))[1], names(data)))) {
    stop(paste0("The given response variable, \"", rownames(attr(tt, "factors"))[1], "\", is not a variable in data."))
  }
  m.try <- try({mf <- model.frame(tt, data)}, TRUE)
  if (class(m.try) == "try-error") {
    stop(paste0(c("All variables of formula must be variables in data.\nVariables not in data: ",
                  paste(attr(tt, "term.labels")[is.na(match(attr(tt, "term.labels"), names(data)))], collapse=", "))), call. = FALSE)}
  treat <- model.response(mf)
  covs <- data[, !is.na(match(names(data), attr(tt, "term.labels"))), drop = FALSE]

  nunique.treat <- nunique(treat)
  if (nunique.treat == 2) {
    treat.type = "binary"
  }
  else if (nunique.treat < 2) {
    stop("treatment must have at least two unique values.", call. = FALSE)
  }
  else if (is.factor(treat) || is.character(treat)) {
    treat.type = "multi"
  }
  else {
    treat.type = "continuous"
  }

#kgkg
  if (length(s.weights) > 0) {
    if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
      stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
    }
    if (is.character(s.weights) && length(s.weights)==1) {
      if (s.weights %in% names(data)) {
        s.weights <- data[, s.weights]
      }
      else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
    }
  }
  else s.weights <- rep(1, length(treat))

  n <- nrow(data)

  ##Process exact
  bad.exact <- FALSE
  acceptable.exacts <- names(data)
  if (missing(exact)) exact.factor <- factor(rep(1, n))
  else if (!is.atomic(exact)) bad.exact <- TRUE
  else if (is.character(exact) && all(exact %in% acceptable.exacts)) {
    exact.factor <- factor(apply(data[, exact, drop = FALSE], 1, paste, collapse = "|"))
  }
  else if (length(exact) == n) exact.factor <- factor(exact)
  else bad.exact <- TRUE

  if (bad.exact) stop("exact must be the quoted names of variables in data for which weighting is to occur within strata or the variable itself.", call. = FALSE)

  if (any(sapply(levels(exact), function(x) nunique(treat) != nunique(treat[exact == x])))) {
    stop("Not all the groups formed by exact contain all treatment levels. Consider coarsening exact.", call. = FALSE)
  }

  ##Process s.weights
  if (length(s.weights) == 0) s.weights <- rep(1, nrow(data))

  w <- p.score <- rep(NA_real_, n)

  for (i in levels(exact.factor)) {
    #Run method
    if (method == "ps") {
      if (treat.type == "binary") {
        estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE", "ATO"), method)
        obj <- weightit2ps(formula = formula,
                           data = data,
                           s.weights = s.weights,
                           subset = exact.factor == i,
                           estimand = estimand,
                           stabilize = stabilize,
                           ps = ps,
                           ...)
      }
      else if (treat.type == "multi") {
        estimand <- process.estimand(estimand, c("ATT", "ATE"), method)
        obj <- weightit2ps.multi(formula = formula,
                           data = data,
                           s.weights = s.weights,
                           subset = exact.factor == i,
                           estimand = estimand,
                           focal = focal,
                           stabilize = stabilize,
                           ps = ps,
                           ...)
      }
      else if (treat.type == "continuous") {
        obj <- weightit2ps.cont(formula = formula,
                                 data = data,
                                 s.weights = s.weights,
                                 subset = exact.factor == i,
                                 stabilize = stabilize,
                                 ps = ps,
                                 ...)
      }

      p.score[exact.factor == i] <- obj$ps
    }
    else if (method == "gbm") {
      if (treat.type == "binary") {
        estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method)
        obj <- weightit2gbm(formula = formula,
                            data = data,
                            s.weights = s.weights,
                            estimand = estimand,
                            subset = exact.factor == i,
                            verbose = verbose, ...)
      }
      else if (treat.type == "multi") {
        estimand <- process.estimand(estimand, c("ATT", "ATE"), method)
        obj <- weightit2gbm.multi(formula = formula,
                            data = data,
                            s.weights = s.weights,
                            estimand = estimand,
                            focal = focal,
                            subset = exact.factor == i,
                            verbose = verbose, ...)
      }
      else stop("Generalized boosted modeling is not compatible with continuous treatments.", call. = FALSE)

      p.score[exact.factor == i] <- obj$ps
    }
    else if (method == "cbps") {
      if (treat.type == "binary") {
        estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method)
        obj <- weightit2cbps(formula = formula,
                             data = data,
                             subset = exact.factor == i,
                             estimand = estimand,
                             verbose = verbose,
                             ...)
      }
      else if (treat.type == "multi") {
        estimand <- process.estimand(estimand, c("ATE"), method)
        obj <- weightit2cbps.multi(formula = formula,
                             data = data,
                             subset = exact.factor == i,
                             verbose = verbose,
                             ...)
      }
      else if (treat.type == "multi") {
        obj <- weightit2cbps.cont(formula = formula,
                                   data = data,
                                   subset = exact.factor == i,
                                   verbose = verbose,
                                   ...)
      }

      p.score[exact.factor == i] <- obj$ps
    }
    else if (method == "nbcbps") {
      if (treat.type == "binary") {
        estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method)
        obj <- weightit2nbcbps(formula = formula,
                             data = data,
                             subset = exact.factor == i,
                             estimand = estimand,
                             verbose = verbose,
                             ...)
      }
      else if (treat.type == "multi") {
        estimand <- process.estimand(estimand, c("ATE"), method)
        obj <- weightit2nbcbps.multi(formula = formula,
                                   data = data,
                                   subset = exact.factor == i,
                                   verbose = verbose,
                                   ...)
      }
      else if (treat.type == "multi") {
        obj <- weightit2nbcbps.cont(formula = formula,
                                  data = data,
                                  subset = exact.factor == i,
                                  verbose = verbose,
                                  ...)
      }

      p.score[exact.factor == i] <- obj$ps
    }
    else if (method == "ebal") {
      if (treat.type == "binary") {
        estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method)
        obj <- weightit2ebal(formula = formula,
                             data = data,
                             s.weights = s.weights,
                             subset = exact.factor == i,
                             estimand = estimand,
                             stabilize = stabilize,
                             verbose = verbose,
                             ...)
      }
      else if (treat.type == "multi") {
        estimand <- process.estimand(estimand, c("ATT", "ATE"), method)
        obj <- weightit2ebal.multi(formula = formula,
                             data = data,
                             s.weights = s.weights,
                             subset = exact.factor == i,
                             estimand = estimand,
                             focal = focal,
                             stabilize = stabilize,
                             verbose = verbose,
                             ...)
      }
      else stop("Entropy balancing is not compatible with continuous treatments.", call. = FALSE)

    }
    else if (method == "sbw") {
      stop("Stable Balancing Weights are not yet available.", call. = FALSE)
    }
    else if (method == "ebcw") {
      if (treat.type == "binary") {
        estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method)
        obj <- weightit2ebcw(formula = formula,
                             data = data,
                             #s.weights = s.weights,
                             subset = exact.factor == i,
                             estimand = estimand,
                             #stabilize = stabilize,
                             ...)
      }
      else if (treat.type == "multi") {
        estimand <- process.estimand(estimand, c("ATE"), method)
        obj <- weightit2ebcw(formula = formula,
                             data = data,
                             #s.weights = s.weights,
                             subset = exact.factor == i,
                             estimand = estimand,
                             #stabilize = stabilize,
                             ...)
      }
      else if (treat.type == "continuous") {
        stop("Empirical balancing calibration weights are not compatible with continuous treatments.", call. = FALSE)
      }

    }

    #Extract weights
    w[exact.factor == i] <- obj$w
  }

  #Trim/Truncate

  if (length(trim.q) == 2) trim <- quantile(w, trim.q)

  if (length(trim) == 2) {
    if (length(ps) == 0 || all(is.na(ps))) {
      to.trim <- !between(w, trim, inclusive = TRUE)
    }
    else {
      to.trim <- !between(ps, trim, inclusive = TRUE)
    }

    w[to.trim] <- 0

    if (nunique(treat[!to.trim]) != nunique(treat)) {
      warning("Trimming will remove an entire treatment group.", call. = FALSE)
    }
  }

  if (length(truncate.q) == 2) truncate <- quantile(w, truncate.q)
  if (length(truncate) == 2) {
    if (length(ps) == 0 || all(is.na(ps))) {
      w[w < min(truncate)] <- min(truncate)
      w[w > max(truncate)] <- max(truncate)
    }
    else {
      w[ps < min(truncate)] <- w[ps == min(truncate)][1]
      w[ps > max(truncate)] <- w[ps == max(truncate)][1]
    }
  }

  #Assemble output object
  out <- list(weights = w,
             treat = treat,
             covs = covs,
             data = data,
             estimand = estimand,
             method = method,
             ps = p.score,
             s.weights = s.weights,
             discarded = NULL,
             treat.type = treat.type)
  class(out) <- "weightit"

  return(out)
}

summary.weightit <- function(weightit, ...) {

}
