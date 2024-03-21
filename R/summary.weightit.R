#' Print and Summarize Output
#' @name summary.weightit
#'
#' @description
#' `summary()` generates a summary of the `weightit` or
#' `weightitMSM` object to evaluate the properties of the estimated
#' weights. `plot()` plots the distribution of the weights.
#'
#' @param object a `weightit` or `weightitMSM` object; the output of
#' a call to [weightit()] or [weightitMSM()].
#' @param top how many of the largest and smallest weights to display. Default
#' is 5.
#' @param ignore.s.weights whether or not to ignore sampling weights when
#' computing the weight summary. If `FALSE`, the default, the estimated
#' weights will be multiplied by the sampling weights (if any) before values
#' are computed.
#' @param binwidth,bins arguments passed to [ggplot2::geom_histogram()] to
#' control the size and/or number of bins.
#' @param x a `summary.weightit` or `summary.weightitMSM` object; the
#' output of a call to `summary.weightit()` or
#' `summary.weightitMSM()`.
#' @param time `numeric`; the time point for which to display the distribution of weights. Default is to plot the distribution for the first time points.
#' @param ... For `plot()`, additional arguments passed to [graphics::hist()] to
#' determine the number of bins, though [ggplot2::geom_histogram()] is actually used to create the plot.
#'
#' @returns
#' For point treatments (i.e., `weightit` objects), a `summary.weightit` object with the following elements:
#' \item{weight.range}{The range (minimum and maximum) weight for each treatment group.}
#' \item{weight.top}{The units with the greatest weights in each treatment group; how many are included is determined by `top`.}
#' \item{coef.of.var (Coef of Var)}{The coefficient of variation (standard deviation divided by mean) of the weights in each treatment group and overall.}
#' \item{scaled.mad (MAD)}{The mean absolute deviation of the weights in each treatment group and overall divided by the mean of the weights in the corresponding group.}
#' \item{negative entropy (Entropy)}{The negative entropy (\eqn{\sum w log(w)}) of the weights in each treatment group and overall divided by the mean of the weights in the corresponding group.}
#' \item{num.zeros}{The number of weights equal to zero.}
#' \item{effective.sample.size}{The effective sample size for each treatment group before and after weighting. See [ESS()].}
#'
#' For longitudinal treatments (i.e., `weightitMSM` objects), a list of
#' the above elements for each treatment period.
#'
#' `plot()` returns a `ggplot` object with a histogram displaying the
#' distribution of the estimated weights. If the estimand is the ATT or ATC,
#' only the weights for the non-focal group(s) will be displayed (since the
#' weights for the focal group are all 1). A dotted line is displayed at the
#' mean of the weights.
#'
#' @seealso
#' [weightit()], [weightitMSM()], [summary()]
#'
#' @examples
#'
#' # See example at ?weightit or ?weightitMSM
#'

#' @exportS3Method summary weightit
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

  out$effective.sample.size <- nn

  if (is_not_null(object$focal)) {
    attr(w, "focal") <- object$focal
  }

  class(out) <- "summary.weightit"

  out
}

#' @exportS3Method print summary.weightit
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

#' @exportS3Method plot summary.weightit
#' @rdname summary.weightit
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

  subtitle <- {
    if (is_not_null(focal))  paste0("For Units Not in Treatment Group \"", focal, "\"")
    else NULL
  }

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

.summary.weightitMSM <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  outnames <- c("weight.range", "weight.top",
                "coef.of.var", "scaled.mad", "negative.entropy",
                "weight.mean",
                "effective.sample.size")

  out.list <- make_list(names(object$treat.list))

  if (ignore.s.weights || is_null(object$s.weights)) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights
  w <- setNames(object$weights*sw, seq_along(sw))
  treat.types <- vapply(object[["treat.list"]], get_treat_type, character(1L))
  stabilized <- is_not_null(object[["stabilization"]])

  for (ti in seq_along(object$treat.list)) {
    out <- make_list(outnames)
    if (treat.types[ti] == "continuous") {
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

      out$effective.sample.size <- nn

      out.list[[ti]] <- out

    }
    else if (treat.types[ti] == "binary") {
      t <- object$treat.list[[ti]]
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

      nn <- make_df(c("Control", "Treated"), c("Unweighted", "Weighted"))
      nn["Unweighted", ] <- c(ESS(sw[t==0]),
                              ESS(sw[t==1]))
      nn["Weighted", ] <- c(ESS(w[t==0]),
                            ESS(w[t==1]))

      out$effective.sample.size <- nn
      out.list[[ti]] <- out

    }
    else if (treat.types[ti] == "multinomial") {
      t <- object$treat.list[[ti]]
      top0 <- setNames(lapply(levels(t), function(x) min(top, sum(t == x))), levels(t))
      out$weight.range <- setNames(lapply(levels(t), function(x) c(min(w[w != 0 & t == x]),
                                                                   max(w[w != 0 & t == x]))),
                                   levels(t))
      out$weight.top <- setNames(lapply(levels(t), function(x) rev(w[t == x][order(abs(w[t == x]), decreasing = TRUE)][seq_len(top0[[x]])])),
                                 levels(t))
      out$coef.of.var <- c(vapply(levels(t), function(x) sd(w[t==x])/mean_fast(w[t==x]), numeric(1L)),
                           overall = sd(w)/mean_fast(w))
      out$scaled.mad <- c(vapply(levels(t), function(x) mean_abs_dev(w[t==x]/mean_fast(w[t==x])), numeric(1L)),
                          overall = mean_abs_dev(w)/mean_fast(w))
      out$negative.entropy <- c(vapply(levels(t), function(x) neg_ent(w[t==x]), numeric(1L)),
                                overall = sum(w[w>0]*log(w[w>0]))/sum(w[w>0]))
      out$num.zeros <- c(vapply(levels(t), function(x) sum(check_if_zero(w[t==x])), numeric(1L)),
                         overall = sum(check_if_zero(w)))
      out$weight.mean <- if (stabilized) mean_fast(w) else NULL

      nn <- make_df(levels(t), c("Unweighted", "Weighted"))
      for (i in levels(t)) {
        nn["Unweighted", i] <- ESS(sw[t==i])
        nn["Weighted", i] <- ESS(w[t==i])
      }

      out$effective.sample.size <- nn
      out.list[[ti]] <- out
    }
    else if (treat.types[ti] == "ordinal") {
      .wrn("Sneaky, sneaky! Ordinal coming soon :)", tidy = FALSE)
    }
  }

  class(out.list) <- "summary.weightitMSM"
  attr(out.list, "weights") <- w

  out.list
}

#' @exportS3Method summary weightitMSM
#' @rdname summary.weightit
summary.weightitMSM <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  out.list <- make_list(names(object$treat.list))

  if (ignore.s.weights || is_null(object$s.weights)) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights

  for (ti in seq_along(object$treat.list)) {
    obj <- as.weightit(object$weights, treat = object$treat.list[[ti]],
                       s.weights = sw, stabilization = object$stabilization)
    out.list[[ti]] <- summary.weightit(obj, top = top, ignore.s.weights = ignore.s.weights, ...)
  }

  class(out.list) <- "summary.weightitMSM"

  out.list
}

#' @exportS3Method plot summary.weightitMSM
print.summary.weightitMSM <- function(x, ...) {
  only.one <- all(vapply(x, function(y) isTRUE(all.equal(x[[1]], y)), logical(1L)))

  cat(paste(rep(" ", 17), collapse = "") %+% underline("Summary of weights") %+% "\n\n")
  for (ti in seq_along(x)) {
    if (!only.one) cat(strikethrough(paste(rep(" ", 22), collapse = "")) %+% italic(" Time " %+% ti %+% " ") %+% strikethrough(paste(rep(" ", 22), collapse = "")) %+% "\n")
    print(x[[ti]])
    cat("\n")
    if (only.one) break
  }

  invisible(x)
}

#' @exportS3Method plot summary.weightitMSM
#' @rdname summary.weightit
plot.summary.weightitMSM <- function(x, binwidth = NULL, bins = NULL, time = 1, ...) {
  if (!is.numeric(time) || length(time) != 1 || time %nin% seq_along(x)) {
    .err("`time` must be a number corresponding to the time point for which to display the distribution of weights")
  }

  plot.summary.weightit(x[[time]], binwidth = binwidth, bins = bins, ...) +
    labs(subtitle = paste0("For Time ", time))
}
