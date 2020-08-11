weightit <- function(formula, data = NULL, method = "ps", estimand = "ATE", stabilize = FALSE, focal = NULL,
                     by = NULL, s.weights = NULL, ps = NULL, moments = NULL, int = FALSE, subclass = NULL,
                     missing = NULL, verbose = FALSE, include.obj = FALSE, ...) {

  ## Checks and processing ----

  A <- list(...)

  #Checks
  if (is_null(ps)) {
    if (is_null(formula) || is_null(class(formula)) || !is.formula(formula, 2)) {
      stop("'formula' must be a formula relating treatment to covariates.", call. = FALSE)
    }
  }
  else {
    if (!(is.character(ps) && length(ps) == 1) && !is.numeric(ps)) {
      stop("The argument to 'ps' must be a vector or data frame of propensity scores or the (quoted) names of variables in 'data' that contain propensity scores.", call. = FALSE)
    }
    if (is.character(ps) && length(ps)==1) {
      if (ps %in% names(data)) {
        ps <- data[[ps]]
      }
      else stop("The name supplied to 'ps' is not the name of a variable in 'data'.", call. = FALSE)
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
  if (length(treat) != nrow(covs)) {
    stop("Treatment and covariates must have the same number of units.", call. = FALSE)
  }

  n <- length(treat)


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

  #Process missing
  if (anyNA(reported.covs)) {
    missing <- process.missing(missing, method, treat.type)
  }
  else missing <- ""

  #Check subclass
  if(is_not_null(subclass)) check.subclass(method, treat.type)

  #Process s.weights
  s.weights <- process.s.weights(s.weights, data)

  if (is_null(s.weights)) s.weights <- rep(1, n)

  ##Process by
  if (is_not_null(A[["exact"]])) {
    message("'by' has replaced 'exact' in the weightit() syntax, but 'exact' will always work.")
    by <- A[["exact"]]
    by.arg <- "exact"
  }
  else by.arg <- "by"

  # processed.by <- process.by(by.name, data = data, treat = treat)
  processed.by <- process.by(by, data = data, treat = treat, by.arg = by.arg)

  #Process moments and int
  moments.int <- process.moments.int(moments, int, method)
  moments <- moments.int[["moments"]]; int <- moments.int[["int"]]

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
                          subclass = subclass,
                          ps = ps,
                          missing = missing,
                          include.obj = include.obj,
                          ...)
  })

  if (all_the_same(obj$w)) stop(paste0("All weights are ", obj$w[1], "."), call. = FALSE)

  warn <- FALSE
  test.w <- obj$w*s.weights
  if (treat.type == "continuous") {if (sd(test.w, na.rm = TRUE)/mean(test.w, na.rm = TRUE) > 4) warn <- TRUE}
  else {if (any(sapply(unique(treat), function(x) sd(test.w[treat == x], na.rm = TRUE)/mean(test.w[treat == x], na.rm = TRUE) > 4))) warn <- TRUE}
  if (warn) warning("Some extreme weights were generated. Examine them with summary() and maybe trim them with trim().", call. = FALSE)

  ## Assemble output object----
  out <- list(weights = obj$weights,
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
              info = obj$info,
              obj = obj$fit.obj)

  out <- clear_null(out)

  class(out) <- "weightit"

  return(out)
  ####----
}

print.weightit <- function(x, ...) {
  treat.type <- get.treat.type(x[["treat"]])
  trim <- attr(x[["weights"]], "trim")

  cat("A " %+% italic("weightit") %+% " object\n")
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
  outnames <- c("weight.range", "weight.top",
                "coef.of.var", "scaled.mad", "negative.entropy",
                "effective.sample.size")
  out <- make_list(outnames)

  if (ignore.s.weights  || is_null(object$s.weights)) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights
  w <- object$weights*sw
  t <- object$treat
  treat.type <- get.treat.type(object[["treat"]])

  if (treat.type == "continuous") {
    out$weight.range <- list(all = c(min(w[w > 0]),
                          max(w[w > 0])))
    top.weights <- sort(w, decreasing = TRUE)[seq_len(top)]
    out$weight.top <- list(all = sort(setNames(top.weights, which(w %in% top.weights)[seq_len(top)])))
    out$coef.of.var <- c(all = sd(w)/mean_fast(w))
    out$scaled.mad <- c(all = mean.abs.dev(w)/mean_fast(w))
    out$negative.entropy <- c(all = sum(w[w>0]*log(w[w>0]))/sum(w[w>0]))
    out$num.zeros <- c(overall = sum(check_if_zero(w)))

    nn <- make_df("Total", c("Unweighted", "Weighted"))
    nn["Unweighted", ] <- ESS(sw)
    nn["Weighted", ] <- ESS(w)

  }
  else if (treat.type == "binary") {
    top0 <- c(treated = min(top, sum(t == 1)),
              control = min(top, sum(t == 0)))
    out$weight.range <- list(treated = c(min(w[w > 0 & t == 1]),
                                         max(w[w > 0 & t == 1])),
                             control = c(min(w[w > 0 & t == 0]),
                                         max(w[w > 0 & t == 0])))
    top.weights <- list(treated = sort(w[t == 1], decreasing = TRUE)[seq_len(top0["treated"])],
                        control = sort(w[t == 0], decreasing = TRUE)[seq_len(top0["control"])])
    out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w %in% top.weights[[x]] & t == {if (x == "control") 0 else 1})[seq_len(top0[x])]))),
                               names(top.weights))
    top.weights <- list(treated = sort(w[t == 1], decreasing = TRUE)[seq_len(top0["treated"])],
                        control = sort(w[t == 0], decreasing = TRUE)[seq_len(top0["control"])])

    out$coef.of.var <- c(treated = sd(w[t==1])/mean_fast(w[t==1]),
                         control = sd(w[t==0])/mean_fast(w[t==0]),
                         overall = sd(w)/mean_fast(w))
    out$scaled.mad <- c(treated = mean.abs.dev(w[t==1])/mean_fast(w[t==1]),
                        control = mean.abs.dev(w[t==0])/mean_fast(w[t==0]),
                        overall = mean.abs.dev(w)/mean_fast(w))
    out$negative.entropy <- c(treated = sum(w[t==1 & w>0]*log(w[t==1 & w>0]))/sum(w[t==1 & w>0]),
                              control = sum(w[t==0 & w>0]*log(w[t==0 & w>0]))/sum(w[t==0 & w>0]),
                              overall = sum(w[w>0]*log(w[w>0]))/sum(w[w>0]))
    out$num.zeros <- c(treated = sum(check_if_zero(w[t==1])),
                       control = sum(check_if_zero(w[t==0])),
                       overall = sum(check_if_zero(w)))

    #dc <- weightit$discarded

    nn <- make_df(c("Control", "Treated"), c("Unweighted", "Weighted"))
    nn["Unweighted", ] <- c(ESS(sw[t==0]),
                            ESS(sw[t==1]))
    nn["Weighted", ] <- c(ESS(w[t==0]),
                          ESS(w[t==1]))
  }
  else if (treat.type == "multinomial") {
    out$weight.range <- setNames(lapply(levels(t), function(x) c(min(w[w > 0 & t == x]),
                                                        max(w[w > 0 & t == x]))),
                                 levels(t))
    top.weights <- setNames(lapply(levels(t), function(x) sort(w[t == x], decreasing = TRUE)[seq_len(top)]),
                            levels(t))
    out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w %in% top.weights[[x]] & t == x)[seq_len(top)]))),
                               names(top.weights))
    out$coef.of.var <- c(vapply(levels(t), function(x) sd(w[t==x])/mean_fast(w[t==x]), numeric(1L)),
                         overall = sd(w)/mean_fast(w))
    out$scaled.mad <- c(vapply(levels(t), function(x) mean.abs.dev(w[t==x])/mean_fast(w[t==x]), numeric(1L)),
                        overall = mean.abs.dev(w)/mean_fast(w))
    out$negative.entropy <- c(vapply(levels(t), function(x) sum(w[t==x & w>0]*log(w[t==x & w>0]))/sum(w[t==x & w>0]), numeric(1L)),
                         overall = sum(w[w>0]*log(w[w>0]))/sum(w[w>0]))
    out$num.zeros <- c(vapply(levels(t), function(x) sum(check_if_zero(w[t==x])), numeric(1L)),
                       overall = sum(check_if_zero(w)))

    nn <- make_df(levels(t), c("Unweighted", "Weighted"))
    for (i in levels(t)) {
      nn["Unweighted", i] <- ESS(sw[t==i])
      nn["Weighted", i] <- ESS(w[t==i])
    }
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
  cat(paste(rep(" ", 17), collapse = "") %+% underline("Summary of weights") %+% "\n\n")

  tryCatch({
    cat("- " %+% italic("Weight ranges") %+% ":\n\n")
    print.data.frame(round_df_char(text_box_plot(x$weight.range, 28), 4), ...)
  })
  df <- setNames(data.frame(do.call("c", lapply(names(x$weight.top), function(x) c(" ", x))),
                   matrix(do.call("c", lapply(x$weight.top, function(x) c(names(x), rep("", top - length(x)), round(x, 4), rep("", top - length(x))))),
               byrow = TRUE, nrow = 2*length(x$weight.top))),
               rep("", 1 + top))
  cat("\n- " %+% italic("Units with", top, "greatest weights by group") %+% ":\n")
  print.data.frame(df, row.names = FALSE)
  cat("\n- " %+% italic("Weight statistics") %+% ":\n\n")
  print.data.frame(round_df_char(setNames(as.data.frame(cbind(x$coef.of.var,
                                                              x$scaled.mad,
                                                              x$negative.entropy,
                                                              x$num.zeros)),
                                        c("Coef of Var", "MAD", "Entropy", "# Zeros")), 3))
  cat("\n- " %+% italic("Effective Sample Sizes") %+% ":\n\n")
  print.data.frame(round_df_char(x$effective.sample.size, 2, pad = " "))
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