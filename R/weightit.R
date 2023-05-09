weightit <- function(formula, data = NULL, method = "glm", estimand = "ATE", stabilize = FALSE, focal = NULL,
                     by = NULL, s.weights = NULL, ps = NULL, moments = NULL, int = FALSE, subclass = NULL,
                     missing = NULL, verbose = FALSE, include.obj = FALSE, ...) {

  ## Checks and processing ----

  A <- list(...)

  #Checks
  if (is_null(formula) || !rlang::is_formula(formula, lhs = TRUE)) {
    .err("`formula` must be a formula relating treatment to covariates")
  }

  #Process treat and covs from formula and data
  t.c <- get_covs_and_treat_from_formula(formula, data)
  reported.covs <- t.c[["reported.covs"]]
  covs <- t.c[["model.covs"]]
  treat <- t.c[["treat"]]
  # treat.name <- t.c[["treat.name"]]

  if (is_null(covs)) {
    .err("no covariates were specified")
  }
  if (is_null(treat)) {
    .err("no treatment variable was specified")
  }
  if (length(treat) != nrow(covs)) {
    .err("the treatment and covariates must have the same number of units")
  }

  n <- length(treat)

  if (anyNA(treat)) {
    .err("no missing values are allowed in the treatment variable")
  }

  #Get treat type
  treat <- assign_treat_type(treat)
  treat.type <- get_treat_type(treat)

  #Process ps
  ps <- process.ps(ps, data, treat)
  if (is_not_null(ps)) {
    method <- "glm"
  }

  ##Process method
  check.acceptable.method(method, msm = FALSE, force = FALSE)

  if (is.character(method)) {
    method <- method.to.proper.method(method)
    attr(method, "name") <- method
  }
  else if (is.function(method)) {
    method.name <- deparse1(substitute(method))
    check.user.method(method)
    attr(method, "name") <- method.name
  }

  #Process estimand and focal
  estimand <- process.estimand(estimand, method, treat.type)
  f.e.r <- process.focal.and.estimand(focal, estimand, treat)
  focal <- f.e.r[["focal"]]
  estimand <- f.e.r[["estimand"]]
  reported.estimand <- f.e.r[["reported.estimand"]]

  #Process missing
  if (anyNA(reported.covs)) {
    missing <- process.missing(missing, method, treat.type)
  }
  else missing <- ""

  #Check subclass
  if (is_not_null(subclass)) check.subclass(method, treat.type)

  #Process s.weights
  s.weights <- process.s.weights(s.weights, data)

  if (is_null(s.weights)) s.weights <- rep(1, n)

  ##Process by
  if (is_not_null(A[["exact"]])) {
    .msg("`by` has replaced `exact` in the `weightit()` syntax, but `exact` will always work")
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

  #Returns weights (weights) and propensity score (ps)
  A[["treat"]] <- treat
  A[["covs"]] <- covs
  A[["s.weights"]] <- s.weights
  A[["by.factor"]] <- attr(processed.by, "by.factor")
  A[["estimand"]] <- estimand
  A[["focal"]] <- focal
  A[["stabilize"]] <- stabilize
  A[["method"]] <- method
  A[["moments"]] <- moments
  A[["int"]] <- int
  A[["subclass"]] <- subclass
  A[["ps"]] <- ps
  A[["missing"]] <- missing
  A[["verbose"]] <- verbose
  A[["include.obj"]] <- include.obj
  A[[".data"]] <- data
  A[[".covs"]] <- reported.covs

  obj <- do.call("weightit.fit", A)

  check_estimated_weights(obj$weights, treat, treat.type, s.weights)

  ## Assemble output object----
  out <- list(weights = obj$weights,
              treat = treat,
              covs = reported.covs,
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

  out
  ####----
}

print.weightit <- function(x, ...) {
  treat.type <- get_treat_type(x[["treat"]])
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
  outnames <- c("weight.range", "weight.top", "weight.mean",
                "coef.of.var", "scaled.mad", "negative.entropy",
                "effective.sample.size")
  out <- make_list(outnames)

  if (ignore.s.weights  || is_null(object$s.weights)) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights
  w <- setNames(object$weights*sw, seq_along(sw))
  t <- object$treat
  treat.type <- get_treat_type(object[["treat"]])
  stabilized <- is_not_null(object[["stabilization"]])

  attr(out, "weights") <- w
  attr(out, "treat") <- t

  if (treat.type == "continuous") {
    out$weight.range <- list(all = c(min(w[w != 0]),
                                     max(w[w != 0])))
    out$weight.top <- list(all = rev(w[order(abs(w), decreasing = TRUE)][seq_len(top)]))
    out$coef.of.var <- c(all = sd(w)/mean_fast(w))
    out$scaled.mad <- c(all = mean_abs_dev(w/mean_fast(w)))
    out$negative.entropy <- c(all = neg_ent(w))
    out$num.zeros <- c(overall = sum(check_if_zero(w)))
    out$weight.mean <- if (stabilized) mean_fast(w) else NULL

    nn <- make_df("Total", c("Unweighted", "Weighted"))
    nn["Unweighted", ] <- ESS(sw)
    nn["Weighted", ] <- ESS(w)

  }
  else if (treat.type == "binary" && !is_(t, c("factor", "character"))) {
    treated <- get_treated_level(t)
    t <- as.integer(t == treated)

    top0 <- c(treated = min(top, sum(t == 1)),
              control = min(top, sum(t == 0)))
    out$weight.range <- list(treated = c(min(w[w != 0 & t == 1]),
                                         max(w[w != 0 & t == 1])),
                             control = c(min(w[w != 0 & t == 0]),
                                         max(w[w != 0 & t == 0])))
    out$weight.top <- list(treated = rev(w[t == 1][order(abs(w[t == 1]), decreasing = TRUE)][seq_len(top0["treated"])]),
                           control = rev(w[t == 0][order(abs(w[t == 0]), decreasing = TRUE)][seq_len(top0["control"])]))
    out$coef.of.var <- c(treated = sd(w[t==1])/mean_fast(w[t==1]),
                         control = sd(w[t==0])/mean_fast(w[t==0]))
    out$scaled.mad <- c(treated = mean_abs_dev(w[t==1]/mean_fast(w[t==1])),
                        control = mean_abs_dev(w[t==0]/mean_fast(w[t==0])))
    out$negative.entropy <- c(treated = neg_ent(w[t==1]),
                              control = neg_ent(w[t==0]))
    out$num.zeros <- c(treated = sum(check_if_zero(w[t==1])),
                       control = sum(check_if_zero(w[t==0])))
    out$weight.mean <- if (stabilized) mean_fast(w) else NULL

    #dc <- weightit$discarded

    nn <- make_df(c("Control", "Treated"), c("Unweighted", "Weighted"))
    nn["Unweighted", ] <- c(ESS(sw[t==0]),
                            ESS(sw[t==1]))
    nn["Weighted", ] <- c(ESS(w[t==0]),
                          ESS(w[t==1]))
  }
  else if (treat.type == "multinomial" || is_(t, c("factor", "character"))) {
    t <- as.factor(t)

    top0 <- setNames(lapply(levels(t), function(x) min(top, sum(t == x))), levels(t))
    out$weight.range <- setNames(lapply(levels(t), function(x) c(min(w[w != 0 & t == x]),
                                                                 max(w[w != 0 & t == x]))),
                                 levels(t))
    out$weight.top <- setNames(lapply(levels(t), function(x) rev(w[t == x][order(abs(w[t == x]), decreasing = TRUE)][seq_len(top0[[x]])])),
                            levels(t))
    out$coef.of.var <- c(vapply(levels(t), function(x) sd(w[t==x])/mean_fast(w[t==x]), numeric(1L)))
    out$scaled.mad <- c(vapply(levels(t), function(x) mean_abs_dev(w[t==x])/mean_fast(w[t==x]), numeric(1L)))
    out$negative.entropy <- c(vapply(levels(t), function(x) neg_ent(w[t==x]), numeric(1L)))
    out$num.zeros <- c(vapply(levels(t), function(x) sum(check_if_zero(w[t==x])), numeric(1L)))
    out$weight.mean <- if (stabilized) mean_fast(w) else NULL

    nn <- make_df(levels(t), c("Unweighted", "Weighted"))
    for (i in levels(t)) {
      nn["Unweighted", i] <- ESS(sw[t==i])
      nn["Weighted", i] <- ESS(w[t==i])
    }
  }
  else if (treat.type == "ordinal") {
    .err("Sneaky, sneaky! Ordinal coming one day :)", tidy = FALSE)
  }

  out$effective.sample.size <- nn

  if (is_not_null(object$focal)) {
    attr(w, "focal") <- object$focal
  }

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
  cat("\n- " %+% italic(sprintf("Units with the %s most extreme weights%s",
                        top, ngettext(length(x$weight.top), "",
                                      " by group"))) %+% ":\n")
  print.data.frame(df, row.names = FALSE)
  cat("\n- " %+% italic("Weight statistics") %+% ":\n\n")
  print.data.frame(round_df_char(setNames(as.data.frame(cbind(x$coef.of.var,
                                                              x$scaled.mad,
                                                              x$negative.entropy,
                                                              x$num.zeros)),
                                          c("Coef of Var", "MAD", "Entropy", "# Zeros")), 3))
  if (is_not_null(x$weight.mean)) cat("\n- " %+% italic("Mean of Weights") %+% " = " %+% round(x$weight.mean, 2) %+% "\n")

  cat("\n- " %+% italic("Effective Sample Sizes") %+% ":\n\n")
  print.data.frame(round_df_char(x$effective.sample.size, 2, pad = " "))
  invisible(x)
}
plot.summary.weightit <- function(x, binwidth = NULL, bins = NULL, ...) {
  w <- attr(x, "weights")
  t <- attr(x, "treat")
  focal <- attr(w, "focal")
  treat.type <- get_treat_type(t)

  A <- list(...)
  if (is_not_null(A[["breaks"]])) {
    breaks <- hist(w, breaks = A[["breaks"]], plot = FALSE)$breaks
    bins <- binwidth <- NULL
  }
  else {
    breaks <- NULL
    if (is_null(bins)) bins <- 20
  }

  if (is_not_null(focal)) subtitle <- paste0("For Units Not in Treatment Group \"", focal, "\"")
  else subtitle <- NULL

  if (treat.type == "continuous") {
  p <- ggplot(data = data.frame(w), mapping = aes(x = w)) +
    geom_histogram(binwidth = binwidth,
                   bins = bins,
                   breaks = breaks,
                   center = mean(w),
                   color = "gray70",
                   fill = "gray70", alpha = 1) +
    scale_y_continuous(expand = expansion(c(0, .05))) +
    geom_vline(xintercept = mean(w), linetype = "12", color = "blue", size = .75) +
    labs(x = "Weight", y = "Count", title = "Distribution of Weights") +
    theme_bw()
  }
  else {
    d <- data.frame(w, t = factor(t))
    if (is_not_null(focal)) d <- d[t != focal,]

    levels(d$t) <- paste("Treat =", levels(d$t))
    w_means <- aggregate(w ~ t, data = d, FUN = mean)

    p <- ggplot(data = d, mapping = aes(x = w)) +
      geom_histogram(binwidth = binwidth,
                     bins = bins,
                     breaks = breaks,
                     # center = mean(w),
                     color = "gray70",
                     fill = "gray70", alpha = 1) +
      scale_y_continuous(expand = expansion(c(0, .05))) +
      geom_vline(data = w_means, aes(xintercept = w), linetype = "12", color = "red") +
      labs(x = "Weight", y = "Count", title = "Distribution of Weights") +
      theme_bw() + facet_wrap(vars(t), ncol = 1, scales = "free") +
      theme(panel.background = element_blank(), panel.border = element_rect(fill = NA, color = "black", size = .25))
  }
  p
}
