weightit <- function(formula, data, method = "ps", estimand = "ATE", stabilize = FALSE, focal = NULL,
                     exact = NULL, s.weights = NULL, ps = NULL, moments = 1L, int = FALSE,
                     #trim = c(0, Inf), trim.q = c(0, 1), truncate = c(0, Inf), truncate.q = c(0, 1),
                     verbose = FALSE, ...) {

  ## Checks and processing ----

  #Checks
  if (length(ps) == 0) {
    if (length(formula) == 0 || length(class(formula)) == 0) {
      stop("formula must be a formula relating treatment to covariates.", call. = FALSE)
    }
    if (missing(data)) {
      data <- environment(formula)
      #stop("Data must be specified.", call. = FALSE)
    }
    if (!is.data.frame(data)) {
      stop("Data must be a data.frame.", call. = FALSE)}
  }
  else {
    if (!(is.character(ps) && length(ps) == 1) && !is.numeric(ps)) {
      stop("The argument to ps must be a vector or data frame of propensity scores or the (quoted) names of variables in data that contain propensity scores.", call. = FALSE)
    }
    if (is.character(ps) && length(ps)==1) {
      if (ps %in% names(data)) {
        ps <- data[[ps]]
      }
      else stop("The name supplied to ps is not the name of a variable in data.", call. = FALSE)
    }
  }

  ##Process method
  bad.method <- FALSE
  acceptable.methods <- c("ps",
                          "gbm", "twang", "gbr",
                          "cbps",
                          "npcbps",
                          "ebal", "entropy", "ebalance",
                          "sbw",
                          "ebcw", "ate")
  if (missing(method) || length(ps) > 0) method <- "ps"
  else if (!is.character(method)) bad.method <- TRUE
  else if (length(method) != 1) bad.method <- TRUE
  else if (!tolower(method) %in% acceptable.methods) bad.method <- TRUE

  if (bad.method) stop("method must be a string of length 1 containing the name of an acceptable weighting method.", call. = FALSE)
  method <- method.to.proper.method(tolower(method))

  #Process treat and covs from formula and data
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  tt <- terms(formula)

#mf0 <<- mf

#covs <- mf[,-1, drop = FALSE]


  #tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.na(match(all.vars(tt[[2]]), names(data)))) {
    stop(paste0("The given response variable, \"", all.vars(tt[[2]]), "\", is not a variable in data."), call. = FALSE)
  }
  vars.mentioned <- all.vars(tt[[3]])

  tryCatch({mf <- eval(mf, parent.frame())}, error = function(e) {
    stop(paste0(c("All variables in formula must be variables in data or objects in the global environment.\nMissing variables: ",
                  paste(vars.mentioned[is.na(match(vars.mentioned, c(names(data), names(.GlobalEnv))))], collapse=", "))), call. = FALSE)})

  covs <- mf[,-1, drop = FALSE]
  #colnames(covs) <- paste0("`", colnames(covs), "`")
  treat <- model.response(mf)
  #covs <- data[!is.na(match(names(data), vars.mentioned[vars.mentioned != all.vars(tt[[2]])]))]

  n <- nrow(data)

  if (any(is.na(covs)) || nrow(covs) != n) {
    stop("No missing values are allowed in the covariates.", call. = FALSE)
  }
  if (any(is.na(treat)) || length(treat) != n) {
    stop("No missing values are allowed in the treatment variable.", call. = FALSE)
  }
  nunique.treat <- nunique(treat)
  if (nunique.treat == 2) {
    treat.type = "binary"
  }
  else if (nunique.treat < 2) {
    stop("The treatment must have at least two unique values.", call. = FALSE)
  }
  else if (is.factor(treat) || is.character(treat)) {
    treat.type = "multinomial"
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
        s.weights <- data[[s.weights]]
      }
      else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
    }
    s.weights.specified <- TRUE
  }
  else {
    s.weights.specified <- FALSE
    s.weights <- rep(1, n)
  }

  ##Process exact
  bad.exact <- FALSE
  acceptable.exacts <- names(data)
  exact.vars <- character(0)

  if (missing(exact) || length(exact) == 0) exact.factor <- factor(rep(1, n))
  else if (!is.atomic(exact)) bad.exact <- TRUE
  else if (is.character(exact) && all(exact %in% acceptable.exacts)) {
    exact.factor <- factor(apply(data[exact], 1, paste, collapse = "|"))
    exact.vars <- exact
  }
  else if (length(exact) == n) {
    exact.factor <- factor(exact)
    exact.vars <- acceptable.exacts[sapply(acceptable.exacts, function(x) equivalent.factors(exact, data[,x]))]
  }
  else bad.exact <- TRUE

  if (bad.exact) stop("exact must be the quoted names of variables in data for which weighting is to occur within strata or the variable itself.", call. = FALSE)

  if (any(sapply(levels(exact.factor), function(x) nunique(treat) != nunique(treat[exact.factor == x])))) {
    stop("Not all the groups formed by exact contain all treatment levels. Consider coarsening exact.", call. = FALSE)
  }

  #Check to ensure formula makes sense with levels
  if (length(exact.vars) > 0) {
    formula <- update.formula(formula, as.formula(paste("~ . -", paste(exact.vars, collapse = " - "))))
  }

  #Process moments and int
  moments.int <- check.moments.int(method, moments, int)
  moments <- moments.int["moments"]; int <- moments.int["int"]

  call <- match.call()
  args <- list(...)

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
          estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE", "ATO"), method, treat.type)
          obj <- weightit2ps(formula = formula,
                             data = data,
                             s.weights = s.weights,
                             subset = exact.factor == i,
                             estimand = estimand,
                             stabilize = stabilize,
                             ps = ps,
                             ...)

        }
        else if (treat.type == "multinomial") {
          estimand <- process.estimand(estimand, c("ATT", "ATE", "ATO"), method, treat.type)
          process.focal(focal, estimand, treat)
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

      }
      else if (method == "gbm") {
        if (treat.type == "binary") {
          estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method, treat.type)
          obj <- weightit2gbm(formula = formula,
                              data = data,
                              s.weights = s.weights,
                              estimand = estimand,
                              subset = exact.factor == i,
                              stabilize = stabilize,
                              ...)
        }
        else if (treat.type == "multinomial") {
          estimand <- process.estimand(estimand, c("ATT", "ATE"), method, treat.type)
          process.focal(focal, estimand, treat)
          obj <- weightit2gbm.multi(formula = formula,
                                    data = data,
                                    s.weights = s.weights,
                                    estimand = estimand,
                                    focal = focal,
                                    subset = exact.factor == i,
                                    stabilize = stabilize,
                                    ...)
        }
        else stop("Generalized boosted modeling is not compatible with continuous treatments.", call. = FALSE)

      }
      else if (method == "cbps") {
        if (treat.type == "binary") {
          estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method, treat.type)
          obj <- weightit2cbps(formula = formula,
                               data = data,
                               subset = exact.factor == i,
                               estimand = estimand,
                               s.weights = s.weights,
                               stabilize = stabilize,
                               ...)

        }
        else if (treat.type == "multinomial") {
          estimand <- process.estimand(estimand, c("ATE"), method, treat.type)
          process.focal(focal, estimand, treat)
          obj <- weightit2cbps.multi(formula = formula,
                                     data = data,
                                     subset = exact.factor == i,
                                     s.weights = s.weights,
                                     stabilize = stabilize,
                                     ...)
        }
        else if (treat.type == "continuous") {
          obj <- weightit2cbps.cont(formula = formula,
                                    data = data,
                                    subset = exact.factor == i,
                                    s.weights = s.weights,
                                    #stabilize = stabilize,
                                    ...)

        }

      }
      else if (method == "npcbps") {
        if (s.weights.specified) stop(paste0("Sampling weights cannot be used with ", method.to.phrase(method), "."),
                                      call. = FALSE)
        if (treat.type == "binary") {
          estimand <- process.estimand(estimand, c("ATE"), method, treat.type)
          obj <- weightit2npcbps(formula = formula,
                                 data = data,
                                 subset = exact.factor == i,
                                 ...)
        }
        else if (treat.type == "multinomial") {
          estimand <- process.estimand(estimand, c("ATE"), method, treat.type)
          process.focal(focal, estimand, treat)
          obj <- weightit2npcbps.multi(formula = formula,
                                       data = data,
                                       subset = exact.factor == i,
                                       ...)
        }
        else if (treat.type == "continuous") {
          obj <- weightit2npcbps.cont(formula = formula,
                                      data = data,
                                      subset = exact.factor == i,
                                      ...)
        }

      }
      else if (method == "ebal") {
        if (treat.type == "binary") {
          estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method, treat.type)
          obj <- weightit2ebal(formula = formula,
                               data = data,
                               s.weights = s.weights,
                               subset = exact.factor == i,
                               estimand = estimand,
                               stabilize = stabilize,
                               moments = moments,
                               int = int,
                               ...)
        }
        else if (treat.type == "multinomial") {
          estimand <- process.estimand(estimand, c("ATT", "ATE"), method, treat.type)
          process.focal(focal, estimand, treat)
          obj <- weightit2ebal.multi(formula = formula,
                                     data = data,
                                     s.weights = s.weights,
                                     subset = exact.factor == i,
                                     estimand = estimand,
                                     focal = focal,
                                     stabilize = stabilize,
                                     moments = moments,
                                     int = int,
                                     ...)
        }
        else stop("Entropy balancing is not compatible with continuous treatments.", call. = FALSE)

      }
      else if (method == "sbw") {
        if (treat.type == "binary") {
          estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method, treat.type)
          obj <- weightit2sbw(formula = formula,
                               data = data,
                               s.weights = s.weights,
                               subset = exact.factor == i,
                               estimand = estimand,
                              moments = moments,
                              int = int,
                              ...)
        }
        else if (treat.type == "multinomial") {
          estimand <- process.estimand(estimand, c("ATT", "ATE"), method, treat.type)
          process.focal(focal, estimand, treat)
          obj <- weightit2sbw.multi(formula = formula,
                              data = data,
                              s.weights = s.weights,
                              subset = exact.factor == i,
                              estimand = estimand,
                              focal = focal,
                              moments = moments,
                              int = int,
                              ...)
        }
        else {
          stop("Stable balancing weights are not compatible with continuous treatments.", call. = FALSE)
        }
      }
      else if (method == "ebcw") {
        if (treat.type == "binary") {
          estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method, treat.type)
          obj <- weightit2ebcw(formula = formula,
                               data = data,
                               s.weights = s.weights,
                               subset = exact.factor == i,
                               estimand = estimand,
                               #stabilize = stabilize,
                               moments = moments,
                               int = int,
                               ...)
        }
        else if (treat.type == "multinomial") {
          estimand <- process.estimand(estimand, c("ATE", "ATT"), method, treat.type)
          process.focal(focal, estimand, treat)
          obj <- weightit2ebcw.multi(formula = formula,
                               data = data,
                               s.weights = s.weights,
                               subset = exact.factor == i,
                               estimand = estimand,
                               focal = focal,
                               #stabilize = stabilize,
                               moments = moments,
                               int = int,
                               ...)
        }
        else {
          stop("Empirical balancing calibration weights are not compatible with continuous treatments.", call. = FALSE)
        }

      }

      #Extract weights
      w[exact.factor == i] <- obj$w
      if (length(obj$ps) > 0) p.score[exact.factor == i] <- obj$ps
    }

  })

  ## Trim/Truncate ----

  # if (length(trim.q) == 2 && all(trim == c(0, Inf))) trim <- quantile(w, trim.q)
  #
  # if (length(trim) == 2) {
  #   # if (length(ps) == 0 || all(is.na(ps))) {
  #     to.trim <- !between(w, trim, inclusive = TRUE)
  #   # }
  #   # else {
  #   #   to.trim <- !between(ps, trim, inclusive = TRUE)
  #   # }
  #
  #   w[to.trim] <- 0
  #
  #   if (treat.type %in% c("binary", "multinomial") && nunique(treat[!to.trim]) != nunique(treat)) {
  #     warning("Trimming will remove an entire treatment group.", call. = FALSE)
  #   }
  # }
  #
  # if (length(truncate.q) == 2 && all(truncate == c(0, Inf))) truncate <- quantile(w, truncate.q)
  # if (length(truncate) == 2) {
  #   # if (length(ps) == 0 || all(is.na(ps))) {
  #     w[w < min(truncate)] <- min(truncate)
  #     w[w > max(truncate)] <- max(truncate)
  #   # }
  #   # else {
  #   #   w[ps < min(truncate)] <- w[ps == min(truncate)][1]
  #   #   w[ps > max(truncate)] <- w[ps == max(truncate)][1]
  #   # }
  # }

  if (sum(w) == 0) stop("All weights are 0.", call. = FALSE)

  ## Assemble output object----
  out <- list(weights = w,
              treat = treat,
              covs = covs,
              data = data,
              estimand = if (treat.type == "continuous") NULL else estimand,
              method = method,
              ps = if (all(is.na(p.score))) NULL else p.score,
              s.weights = s.weights,
              #discarded = NULL,
              treat.type = treat.type,
              focal = focal,
              call = call)
  class(out) <- "weightit"

  return(out)
  ####----
}

print.weightit <- function(x, ...) {
  cat("A weightit object\n")
  cat(paste0(" - method: \"", x$method, "\" (", method.to.phrase(x$method), ")\n"))
  cat(paste0(" - number of obs.: ", length(x$weights), "\n"))
  cat(paste0(" - sampling weights: ", ifelse(max(x$s.weights) - min(x$s.weights) < sqrt(.Machine$double.eps),
                                             "none", "present"), "\n"))
  cat(paste0(" - treatment: ", ifelse(x$treat.type == "continuous", "continuous", paste0(nunique(x$treat), "-category", ifelse(x$treat.type == "multinomial", paste0(" (", paste(levels(x$treat), collapse = ", "), ")"), ""))), "\n"))
  if (length(x$estimand) > 0) cat(paste0(" - estimand: ", x$estimand, ifelse(length(x$focal)>0, paste0(" (focal: ", x$focal, ")"), ""), "\n"))
  cat(paste0(" - covariates: ", ifelse(length(names(x$covs)) > 60, "too many to name", paste(names(x$covs), collapse = ", ")), "\n"))
  invisible(x)
}
summary.weightit <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  outnames <- c("weight.range", "weight.top","weight.ratio",
                "coef.of.var",
                "effective.sample.size")
  out <- setNames(vector("list", length(outnames)), outnames)

  if (ignore.s.weights) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights
  w <- object$weights*sw
  t <- object$treat

  if (object$treat.type == "continuous") {
    out$weight.range <- list(all = c(min(w[w > 0]),
                          max(w[w > 0])))
    out$weight.ratio <- c(all = out$weight.range[["all"]][2]/out$weight.range[["all"]][1])
    top.weights <- sort(w, decreasing = TRUE)[seq_len(top)]
    out$weight.top <- list(all = sort(setNames(top.weights, which(w %in% top.weights)[seq_len(top)])))
    out$coef.of.var <- c(all = sd(w)/mean(w))

    nn <- as.data.frame(matrix(0, ncol = 1, nrow = 2))
    nn[1, ] <- (sum(sw)^2)/sum(sw^2)
    nn[2, ] <- (sum(w)^2)/sum((w)^2)
    dimnames(nn) <- list(c("Unweighted", "Weighted"),
                         c("Total"))

  }
  else if (object$treat.type == "binary") {
    top0 <- c(treated = min(top, sum(t == 1)),
              control = min(top, sum(t == 0)))
    out$weight.range <- list(treated = c(min(w[w > 0 & t == 1]),
                                         max(w[w > 0 & t == 1])),
                             control = c(min(w[w > 0 & t == 0]),
                                         max(w[w > 0 & t == 0])))
    out$weight.ratio <- c(treated = out$weight.range$treated[2]/out$weight.range$treated[1],
                          control = out$weight.range$control[2]/out$weight.range$control[1],
                          overall = max(unlist(out$weight.range)/min(unlist(out$weight.range))))
    top.weights <- list(treated = sort(w[t == 1], decreasing = TRUE)[seq_len(top0["treated"])],
                        control = sort(w[t == 0], decreasing = TRUE)[seq_len(top0["control"])])
    out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w %in% top.weights[[x]] & t == {if (x == "control") 0 else 1})[seq_len(top0[x])]))),
                               names(top.weights))
    top.weights <- list(treated = sort(w[t == 1], decreasing = TRUE)[seq_len(top0["treated"])],
                        control = sort(w[t == 0], decreasing = TRUE)[seq_len(top0["control"])])

    out$coef.of.var <- c(treated = sd(w[t==1])/mean(w[t==1]),
                         control = sd(w[t==0])/mean(w[t==0]),
                         overall = sd(w)/mean(w))

    #dc <- weightit$discarded

    nn <- as.data.frame(matrix(0, nrow = 2, ncol = 2))
    nn[1, ] <- c((sum(sw[t==0])^2)/sum(sw[t==0]^2),
                 (sum(sw[t==1])^2)/sum(sw[t==1]^2))
    nn[2, ] <- c((sum(w[t==0])^2)/sum((w[t==0])^2),
                 (sum(w[t==1])^2)/sum((w[t==1])^2))
    # nn[3, ] <- c(sum(t==0 & dc==1), #Discarded
    #              sum(t==1 & dc==1))
    dimnames(nn) <- list(c("Unweighted", "Weighted"),
                         c("Control", "Treated"))
  }
  else if (object$treat.type == "multinomial") {
    out$weight.range <- setNames(lapply(levels(t), function(x) c(min(w[w > 0 & t == x]),
                                                        max(w[w > 0 & t == x]))),
                                 levels(t))
    out$weight.ratio <- setNames(c(sapply(out$weight.range, function(x) x[2]/x[1]),
                          max(unlist(out$weight.range)/min(unlist(out$weight.range)))),
                          c(levels(t), "overall"))
    top.weights <- setNames(lapply(levels(t), function(x) sort(w[t == x], decreasing = TRUE)[seq_len(top)]),
                            levels(t))
    out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w %in% top.weights[[x]] & t == x)[seq_len(top)]))),
                               names(top.weights))
    out$coef.of.var <- c(sapply(levels(t), function(x) sd(w[t==x])/mean(w[t==x])),
                         overall = sd(w)/mean(w))

    nn <- as.data.frame(matrix(0, nrow = 2, ncol = nunique(t)))
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
print.summary.weightit <- function(x, ...) {
  top <- max(lengths(x$weight.top))
  cat("Summary of weights:\n\n")
  cat("- Weight ranges:\n")
  print.data.frame(round_df_char(text.box.plot(x$weight.range, 28), 4), ...)
  df <- setNames(data.frame(do.call("c", lapply(names(x$weight.top), function(x) c(" ", x))),
                   matrix(do.call("c", lapply(x$weight.top, function(x) c(names(x), rep("", top - length(x)), round(x, 4), rep("", top - length(x))))),
               byrow = TRUE, nrow = 2*length(x$weight.top))),
               rep("", 1 + top))
  cat(paste("\n- Units with", top, "greatest weights by group:\n"))
  print.data.frame(df, row.names = FALSE)
  cat("\n")
  print.data.frame(round_df_char(as.data.frame(matrix(c(x$weight.ratio, x$coef.of.var), ncol = 2,
                                        dimnames = list(names(x$weight.ratio),
                                                        c("Ratio", "Coef of Var")))), 4))
  cat("\n- Effective Sample Sizes:\n")
  print.data.frame(round_df_char(x$effective.sample.size, 3))
  invisible(x)
}
