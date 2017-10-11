weightit <- function(formula, data, method, estimand = "ATE", stabilize = FALSE, focal = NULL,
                     exact = NULL, s.weights = NULL, ps = NULL, trim = c(0, Inf), trim.q = c(0, 1),
                     truncate = c(0, Inf), truncate.q = c(0, 1), verbose = FALSE, ...) {

  ## Checks and processing ----

  #Checks
  if (length(ps) == 0) {
    if (length(formula) == 0 || length(class(formula)) == 0) {
      stop("formula must be a formula relating treatment to covariates.", call. = FALSE)
    }
    if (missing(data)) {
      stop("Data must be specified.", call. = FALSE)}
    if (!is.data.frame(data)) {
      stop("Data must be a data.frame.", call. = FALSE)}
  }
  else {
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

  n <- nrow(data)

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

  #Process s.weights
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
  else s.weights <- rep(1, n)

  ##Process exact
  bad.exact <- FALSE
  acceptable.exacts <- names(data)
  exact.vars <- character(0)

  if (missing(exact)) exact.factor <- factor(rep(1, n))
  else if (!is.atomic(exact)) bad.exact <- TRUE
  else if (is.character(exact) && all(exact %in% acceptable.exacts)) {
    exact.factor <- factor(apply(data[, exact, drop = FALSE], 1, paste, collapse = "|"))
    exact.vars <- exact
  }
  else if (length(exact) == n) {
    exact.factor <- factor(exact)
    exact.vars <- acceptable.exacts[sapply(acceptable.exacts, function(x) equivalent.factors(exact, data[,x]))]
  }
  else bad.exact <- TRUE

  if (bad.exact) stop("exact must be the quoted names of variables in data for which weighting is to occur within strata or the variable itself.", call. = FALSE)

  if (any(sapply(levels(exact), function(x) nunique(treat) != nunique(treat[exact == x])))) {
    stop("Not all the groups formed by exact contain all treatment levels. Consider coarsening exact.", call. = FALSE)
  }

  #Check to ensure formula makes sense with levels
  if (length(exact.vars) > 0) {
    formula <- update.formula(formula, as.formula(paste("~ . -", paste(exact.vars, collapse = " - "))))
  }

  ## Running models ----
  w <- p.score <- rep(NA_real_, n)

  if (verbose) {
    eval.verbose <- base::eval
  }
  else eval.verbose <- utils::capture.output

  eval.verbose({
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

  })

  ## Trim/Truncate ----

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
  ####----
}

print.weightit <- function(weightit, ...) {
  cat("A weightit object\n")
  cat(paste0(" - method: \"cbps\" (", method.to.phrase(weightit$method), ")\n"))
  cat(paste0(" - number of obs.: ", length(weightit$weights), "\n"))
  cat(paste0(" - estimand: ", weightit$estimand, "\n"))
  cat(paste0(" - treatment: ", ifelse(weightit$treat.type == "continuous", "continuous", paste0(nunique(weightit$treat), "-category")), "\n"))
  cat(paste0(" - covariates: ", ifelse(length(names(weightit$covs)) > 60, "too many to name", paste(names(weightit$covs), collapse = ", "))))
}
summary.weightit <- function(weightit, top = 5, ...) {
  outnames <- c("weight.range", "weight.top","weight.ratio",
                "effective.sample.size", "coef.of.var")
  out <- setNames(vector("list", length(outnames)), outnames)

  sw <- weightit$s.weights
  w <- weightit$weights*sw
  t <- weightit$treat

  if (weightit$treat.type == "continuous") {
    out$weight.range <- list(all = c(min(w[w > 0]),
                          max(w[w > 0])))
    out$weight.ratio <- list(all = out$weight.range[2]/out$weight.range[1])
    top.weights <- sort(w, decreasing = TRUE)[seq_len(top)]
    out$weight.top <- list(all = sort(setNames(top.weights, which(w %in% top.weights)[seq_len(top)])))
    out$coef.of.var <- c(all = sd(w)/mean(w))

    nn <- as.data.frame(matrix(0, ncol = 1, nrow = 2))
    nn[1, ] <- (sum(sw)^2)/sum(sw^2)
    nn[2, ] <- (sum(w)^2)/sum((w)^2)
    dimnames(nn) <- list(c("Unweighted", "Weighted"),
                         c("Total"))

  }
  else if (weightit$treat.type == "binary") {
    out$weight.range <- list(treated = c(min(w[w > 0 & t == 1]),
                                         max(w[w > 0 & t == 1])),
                             control = c(min(w[w > 0 & t == 0]),
                                         max(w[w > 0 & t == 0])))
    out$weight.ratio <- c(treated = out$weight.range$treated[2]/out$weight.range$treated[1],
                          control = out$weight.range$control[2]/out$weight.range$control[1],
                          overall = max(unlist(out$weight.range)/min(unlist(out$weight.range))))
    top.weights <- list(treated = sort(w[t == 1], decreasing = TRUE)[seq_len(top)],
                        control = sort(w[t == 0], decreasing = TRUE)[seq_len(top)])
    out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w[t == ifelse(x == "control", 0, 1)] %in% top.weights[[x]])[seq_len(top)]))),
                               names(top.weights))
    out$coef.of.var <- c(treated = sd(w[t==1])/mean(w[t==1]),
                         control = sd(w[t==0])/mean(w[t==0]),
                         overall = sd(w)/mean(w))

    #dc <- weightit$discarded

    nn <- as.data.frame(matrix(0, ncol = 2, nrow = 2))
    nn[1, ] <- c((sum(sw[t==0])^2)/sum(sw[t==0]^2),
                 (sum(sw[t==1])^2)/sum(sw[t==1]^2))
    nn[2, ] <- c((sum(w[t==0])^2)/sum((w[t==0])^2),
                 (sum(w[t==1])^2)/sum((w[t==1])^2))
    # nn[3, ] <- c(sum(t==0 & dc==1), #Discarded
    #              sum(t==1 & dc==1))
    dimnames(nn) <- list(c("Unweighted", "Weighted"),
                         c("Control", "Treated"))
  }
  else if (weightit$treat.type == "multi") {
    out$weight.range <- lapply(levels(t), function(x) c(min(w[w > 0 & t == x]),
                                                        max(w[w > 0 & t == x])))
    out$weight.ratio <- sapply(out$weight.range, function(x) x[2]/x[1])
    top.weights <- lapply(levels(t), function(x) sort(w[t == x], decreasing = TRUE)[seq_len(top)])
    out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w[t == x] %in% top.weights[[x]])[seq_len(top)]))),
                               names(top.weights))
    out$coef.of.var <- c(sapply(levels(t), function(x) sd(w[t==x])/mean(w[t==x])),
                         overall = sd(w)/mean(w))

    nn <- as.data.frame(matrix(0, ncol = 2, nrow = nunique(t)))
    for (i in seq_len(nunique(t))) {
      nn[1, i] <- (sum(sw[t==levels(t)[i]])^2)/sum(sw[t==levels(t)[i]]^2)
      nn[2, i] <- (sum(w[t==levels(t)[i]])^2)/sum((w[t==levels(t)[i]])^2)
      # nn[3, i] <- sum(t==levels(t)[i] & dc==1) #Discarded
    }

    dimnames(nn) <- list(c("Unweighted", "Weighted"),
                         levels(t))

  }

  out$effective.sample.size <- nn


  class(out) <- "summary.weightit"
  return(out)
}
print.summary.weightit <- function(summary.weightit, ...) {
  cat("Summary of weights:\n\n")
  cat("- Weight ranges:\n")
  print.data.frame(round_df(text.box.plot(summary.weightit$weight.range, 28), 4))
  cat(paste("\n- Units with", length(summary.weightit$weight.top[[1]]), "greatest weights by group:\n"))
  for (i in seq_along(summary.weightit$weight.top)) {
    print(c(names(summary.weightit$weight.top)[i], setNames(as.character(round(summary.weightit$weight.top[[i]], 4)), names(summary.weightit$weight.top[[i]]))), quote = FALSE)
  }
  cat("\n")
  print.data.frame(as.data.frame(matrix(round(c(summary.weightit$weight.ratio, summary.weightit$coef.of.var), 4), ncol = 2,
                                        dimnames = list(names(summary.weightit$weight.ratio),
                                                        c("Ratio", "Coef of Var")))))
  cat("\n- Effective Sample Sizes:\n")
  print.data.frame(round(summary.weightit$effective.sample.size, 3))
}
