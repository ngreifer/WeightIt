weightit <- function(formula, data = NULL, method = "ps", estimand = "ATE", stabilize = FALSE, focal = NULL,
                     by = NULL, s.weights = NULL, ps = NULL, moments = 1L, int = FALSE,
                     verbose = FALSE, include.obj = FALSE, ...) {

  ## Checks and processing ----

  A <- list(...)

  #Checks
  if (is_null(ps)) {
    if (is_null(formula) || is_null(class(formula)) || !is.formula(formula, 2)) {
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
    method <- "ps"
  }

  ##Process method
  check.acceptable.method(method, msm = FALSE, force = FALSE)

  if (is.character(method)) {
    method <- method.to.proper.method(method)
    attr(method, "name") <- method
  }
  else if (is.function(method)) {
    method.name <- paste(deparse(substitute(method)))
    check.user.method(method)
    attr(method, "name") <- method.name
  }

  #Process treat and covs from formula and data
  t.c <- get.covs.and.treat.from.formula(formula, data)
  reported.covs <- t.c[["reported.covs"]]
  covs <- t.c[["model.covs"]]
  treat <- t.c[["treat"]]
  # treat.name <- t.c[["treat.name"]]

  if (is_null(covs)) stop("No covariates were specified.", call. = FALSE)
  if (is_null(treat)) stop("No treatment variable was specified.", call. = FALSE)

  n <- length(treat)

  if (anyNA(reported.covs) || nrow(reported.covs) != n) {
    warning("Missing values are present in the covariates. See ?weightit for information on how these are handled.", call. = FALSE)
    #stop("No missing values are allowed in the covariates.", call. = FALSE)
  }
  if (anyNA(treat)) {
    stop("No missing values are allowed in the treatment variable.", call. = FALSE)
  }

  #Get treat type
  treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  #Process estimand and focal
  estimand <- process.estimand(estimand, method, treat.type)
  f.e.r <- process.focal.and.estimand(focal, estimand, treat, treat.type)
  focal <- f.e.r[["focal"]]
  estimand <- f.e.r[["estimand"]]
  reported.estimand <- f.e.r[["reported.estimand"]]

  #Process s.weights
  s.weights <- process.s.weights(s.weights, data)

  if (is_null(s.weights)) s.weights <- rep(1, n)

  ##Process by
  if (is_not_null(A[["exact"]])) {
    message("'by' has replaced 'exact' in the weightit() syntax, but 'exact' will always work.")
    # by.name <- deparse(A[["exact"]])
    by <- A[["exact"]]
    by.arg <- "exact"
  }
  else by.arg <- "by"

  # processed.by <- process.by(by.name, data = data, treat = treat)
  processed.by <- process.by(by, data = data, treat = treat, by.arg = by.arg)

  #Process moments and int
  moments.int <- process.moments.int(moments, int, method)
  moments <- moments.int["moments"]; int <- moments.int["int"]

  call <- match.call()
  # args <- list(...)

  ## Running models ----

  if (verbose) eval.verbose <- base::eval
  else eval.verbose <- utils::capture.output

  eval.verbose({
    #Returns weights (w) and propensty score (ps)
      obj <- weightit.fit(treat = treat,
                          covs = covs,
                          treat.type = treat.type,
                          s.weights = s.weights,
                          by.factor = attr(processed.by, "by.factor"),
                          estimand = estimand,
                          focal = focal,
                          stabilize = stabilize,
                          method = method,
                          moments = moments,
                          int = int,
                          ps = ps,
                          include.obj = include.obj,
                          ...)
  })

  if (all_the_same(obj$w)) stop(paste0("All weights are ", obj$w[1], "."), call. = FALSE)

  warn <- FALSE
  test.w <- obj$w*s.weights
  if (treat.type == "continuous") {if (sd(test.w, na.rm = TRUE)/mean(test.w, na.rm = TRUE) > 4) warn <- TRUE}
  else {if (any(sapply(unique(treat), function(x) sd(test.w[treat == x], na.rm = TRUE)/mean(test.w[treat == x], na.rm = TRUE) > 4))) warn <- TRUE}
  if (warn) warning("Some extreme weights were generated. Examine them with summary() and maybe trim them with trim().", call. = FALSE)
  # #Create new data set
  # #treat, covs, data (not in treat or covs), by
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
              by = processed.by,
              call = call,
              obj = obj$fit.obj)

  out <- clear_null(out)

  class(out) <- "weightit"

  return(out)
  ####----
}

print.weightit <- function(x, ...) {
  treat.type <- get.treat.type(x[["treat"]])
  trim <- attr(x[["weights"]], "trim")

  cat("A weightit object\n")
  if (is_not_null(x[["method"]])) cat(paste0(" - method: \"", attr(x[["method"]], "name"), "\" (", method.to.phrase(x[["method"]]), ")\n"))
  cat(paste0(" - number of obs.: ", length(x[["weights"]]), "\n"))
  cat(paste0(" - sampling weights: ", if (is_null(x[["s.weights"]]) || all_the_same(x[["s.weights"]])) "none" else "present", "\n"))
  cat(paste0(" - treatment: ", ifelse(treat.type == "continuous", "continuous", paste0(nunique(x[["treat"]]), "-category", ifelse(treat.type == "multinomial", paste0(" (", paste(levels(x[["treat"]]), collapse = ", "), ")"), ""))), "\n"))
  if (is_not_null(x[["estimand"]])) cat(paste0(" - estimand: ", x[["estimand"]], ifelse(is_not_null(x[["focal"]]), paste0(" (focal: ", x[["focal"]], ")"), ""), "\n"))
  if (is_not_null(x[["covs"]])) cat(paste0(" - covariates: ", ifelse(length(names(x[["covs"]])) > 60, "too many to name", paste(names(x[["covs"]]), collapse = ", ")), "\n"))
  if (is_not_null(x[["by"]])) {
    cat(paste0(" - by: ", paste(names(x[["by"]]), collapse = ", "), "\n"))
  }
  if (is_not_null(trim)) {
    if (trim < 1) {
      if (attr(x[["weights"]], "trim.lower")) trim <- c(1 - trim, trim)
      cat(paste(" - weights trimmed at", word_list(paste0(round(100*trim, 2), "%")), "\n"))
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

  if (ignore.s.weights  || is_null(object$s.weights)) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights
  w <- object$weights*sw
  t <- object$treat
  treat.type <- get.treat.type(object[["treat"]])

  if (treat.type == "continuous") {
    out$weight.range <- list(all = c(min(w[w > 0]),
                          max(w[w > 0])))
    out$weight.ratio <- c(all = out$weight.range[["all"]][2]/out$weight.range[["all"]][1])
    top.weights <- sort(w, decreasing = TRUE)[seq_len(top)]
    out$weight.top <- list(all = sort(setNames(top.weights, which(w %in% top.weights)[seq_len(top)])))
    out$coef.of.var <- c(all = sd(w)/mean(w))

    nn <- as.data.frame(matrix(0, ncol = 1, nrow = 2))
    nn[1, ] <- ESS(sw)
    nn[2, ] <- ESS(w)
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
    nn[1, ] <- c(ESS(sw[t==0]),
                 ESS(sw[t==1]))
    nn[2, ] <- c(ESS(w[t==0]),
                 ESS(w[t==1]))
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
      nn[1, i] <- ESS(sw[t==levels(t)[i]])
      nn[2, i] <- ESS(w[t==levels(t)[i]])
      # nn[3, i] <- sum(t==levels(t)[i] & dc==1) #Discarded
    }
    dimnames(nn) <- list(c("Unweighted", "Weighted"),
                         levels(t))
  }
  else if (treat.type == "ordinal") {
    stop("Sneaky, sneaky! Ordinal coming soon :)", call. = FALSE)
  }

  out$effective.sample.size <- nn

  if (is_not_null(object$focal)) {
    w <- w[t != object$focal]
    attr(w, "focal") <- object$focal
  }
  attr(out, "weights") <- w
  class(out) <- "summary.weightit"
  return(out)
}
print.summary.weightit <- function(x, ...) {
  top <- max(lengths(x$weight.top))
  cat("Summary of weights:\n\n")
  cat("- Weight ranges:\n")
  print.data.frame(round_df_char(text_box_plot(x$weight.range, 28), 4), ...)
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
plot.summary.weightit <- function(x, ...) {
  w <- attr(x, "weights")
  focal <- attr(w, "focal")

  if (is_not_null(focal)) subtitle <- paste0("For Units Not in Treatment Group \"", focal, "\"")
  else subtitle <- NULL

  p <- ggplot(data = data.frame(w), mapping = aes(x = w)) +
    geom_histogram(breaks = hist(w, plot = FALSE,
                                 ...)$breaks,
                   color = "black",
                   fill = "gray", alpha = .8) +
    scale_y_continuous(expand = expand_scale(c(0, .05))) +
    geom_vline(xintercept = mean(w), linetype = "12") +
    labs(x = "Weight", y = "Count", title = "Distribution of Weights",
         subtitle = subtitle) +
    theme_bw()
  p
}