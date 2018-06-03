weightit <- function(formula, data = NULL, method = "ps", estimand = "ATE", stabilize = FALSE, focal = NULL,
                     exact = NULL, s.weights = NULL, ps = NULL, moments = 1L, int = FALSE,
                     verbose = FALSE, ...) {

  ## Checks and processing ----

  #Checks
  if (is_null(ps)) {
    if (is_null(formula) || is_null(class(formula))) {
      stop("formula must be a formula relating treatment to covariates.", call. = FALSE)
    }
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
  if (missing(method) || is_not_null(ps)) method <- "ps"
  else if (!is.character(method)) bad.method <- TRUE
  else if (is_null(method) || length(method) > 1) bad.method <- TRUE
  else if (tolower(method) %nin% acceptable.methods) bad.method <- TRUE

  if (bad.method) stop("method must be a string of length 1 containing the name of an acceptable weighting method.", call. = FALSE)
  method <- method.to.proper.method(tolower(method))

  #Process treat and covs from formula and data
  t.c <- get.covs.and.treat.from.formula(formula, data)
  reported.covs <- t.c[["reported.covs"]]
  covs <- t.c[["model.covs"]]
  treat <- t.c[["treat"]]
  treat.name <- t.c[["treat.name"]]

  if (is_null(covs)) stop("No covariates were specified.", call. = FALSE)
  if (is_null(treat)) stop("No treatment variable was specified.", call. = FALSE)

  n <- length(treat)

  if (any(is.na(reported.covs)) || nrow(reported.covs) != n) {
    warning("Missing values are present in the covariates. See ?weightit for information on how these are handled.", call. = FALSE)
    #stop("No missing values are allowed in the covariates.", call. = FALSE)
  }
  if (any(is.na(treat))) {
    stop("No missing values are allowed in the treatment variable.", call. = FALSE)
  }

  #Get treat type
  treat <- get.treat.type(treat)
  treat.type <- attr(treat, "treat.type")


  #Process estimand and focal
  estimand <- process.estimand(estimand, method, treat.type)
  f.e.r <- process.focal.and.estimand(focal, estimand, treat, treat.type)
  focal <- f.e.r[["focal"]]
  estimand <- f.e.r[["estimand"]]
  reported.estimand <- f.e.r[["reported.estimand"]]

  #Process s.weights
  s.weights <- process.s.weights(s.weights, data)

  if (is_null(s.weights)) s.weights <- rep(1, n)

  ##Process exact
  processed.exact <- process.exact(exact = exact, data= data,
                                   treat = treat)

  # #Recreate data and formula
  # w.data <- data.frame(treat, covs)
  # w.formula <- formula(w.data)

  #Process moments and int
  moments.int <- check.moments.int(method, moments, int)
  moments <- moments.int["moments"]; int <- moments.int["int"]

  call <- match.call()
  args <- list(...)

  ## Running models ----

  if (verbose) eval.verbose <- base::eval
  else eval.verbose <- utils::capture.output

  eval.verbose({
    #Returns weights (w) and propensty score (ps)
      obj <- weightit.fit(treat = treat,
                          covs = covs,
                          treat.type = treat.type,
                          s.weights = s.weights,
                          exact.factor = processed.exact[["exact.factor"]],
                          estimand = estimand,
                          focal = focal,
                          stabilize = stabilize,
                          method = method,
                          moments = moments,
                          int = int,
                          ps = ps,
                          ...)
  })

  if (all_the_same(obj$w)) stop(paste0("All weights are ", obj$w[1], "."), call. = FALSE)

  warn <- FALSE
  test.w <- obj$w*s.weights
  if (treat.type == "continuous") {if (sd(test.w)/mean(test.w) > 4) warn <- TRUE}
  else {if (any(sapply(unique(treat), function(x) sd(test.w[treat == x])/mean(test.w[treat == x]) > 4))) warn <- TRUE}
  if (warn) warning("Some extreme weights were generated. Examine them with summary() and maybe trim them with trim().", call. = FALSE)
  # #Create new data set
  # #treat, covs, data (not in treat or covs), exact
  # treat.in.data <- treat; attr(treat.in.data, "treat.type") <- NULL
  # data.list <- list(treat.in.data, reported.covs)
  # o.data <- setNames(do.call("data.frame", data.list[sapply(data.list, is_not_null)]),
  #                    c(treat.name, names(reported.covs)))
  # o.data <- data.frame(o.data, data[names(data) %nin% names(o.data)])

  ## Assemble output object----
  out <- list(weights = obj$w,
              treat = treat,
              covs = reported.covs,
              #data = o.data,
              estimand = if (treat.type == "continuous") NULL else reported.estimand,
              method = method,
              ps = if (is_null(obj$ps) || all(is.na(obj$ps))) NULL else obj$ps,
              s.weights = s.weights,
              #discarded = NULL,
              focal = if (reported.estimand == "ATT") focal else NULL,
              exact = processed.exact[["exact.components"]],
              call = call)
  class(out) <- "weightit"

  return(out)
  ####----
}

print.weightit <- function(x, ...) {
  treat.type <- attr(x[["treat"]], "treat.type")
  trim <- attr(x[["weights"]], "trim")

  cat("A weightit object\n")
  cat(paste0(" - method: \"", x[["method"]], "\" (", method.to.phrase(x[["method"]]), ")\n"))
  cat(paste0(" - number of obs.: ", length(x[["weights"]]), "\n"))
  cat(paste0(" - sampling weights: ", ifelse(all_the_same(x[["s.weights"]]),"none", "present"), "\n"))
  cat(paste0(" - treatment: ", ifelse(treat.type == "continuous", "continuous", paste0(nunique(x[["treat"]]), "-category", ifelse(treat.type == "multinomial", paste0(" (", paste(levels(x[["treat"]]), collapse = ", "), ")"), ""))), "\n"))
  if (is_not_null(x[["estimand"]])) cat(paste0(" - estimand: ", x[["estimand"]], ifelse(is_not_null(x[["focal"]]), paste0(" (focal: ", x[["focal"]], ")"), ""), "\n"))
  cat(paste0(" - covariates: ", ifelse(length(names(x[["covs"]])) > 60, "too many to name", paste(names(x[["covs"]]), collapse = ", ")), "\n"))
  if (is_not_null(x[["exact"]])) {
    cat(paste0(" - exact: ", paste(names(x[["exact"]]), collapse = ", "), "\n"))
  }
  if (is_not_null(trim)) {
    if (trim < 1) {
      if (attr(x[["weights"]], "trim.lower")) trim <- c(1 - trim, trim)
      cat(paste(" - weights trimmed at", word.list(paste0(round(100*trim, 2), "%")), "\n"))
    }
    else {
      if (attr(x[["weights"]], "trim.lower")) t.b <- "top and bottom" else t.b <- "top"
      cat(paste(" - weights trimmed at the", t.b, trim, "\n"))
    }
  }
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
  treat.type <- attr(object[["treat"]], "treat.type")

  if (treat.type == "continuous") {
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
  else if (treat.type == "binary") {
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
  else if (treat.type == "multinomial") {
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
  else if (treat.type == "ordinal") {
    stop("Sneaky, sneaky! Ordinal coming soon :)", call. = FALSE)
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
