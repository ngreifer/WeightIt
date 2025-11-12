#' Print and Summarize Output
#'
#' @description
#' `summary()` generates a summary of the `weightit` or
#' `weightitMSM` object to evaluate the properties of the estimated weights.
#' `plot()` plots the distribution of the weights. `nobs()` extracts the number
#' of observations.
#'
#' @param object a `weightit` or `weightitMSM` object; the output of a call to
#'   [weightit()] or [weightitMSM()].
#' @param top how many of the largest and smallest weights to display. Default
#'   is 5. Ignored when `weight.range = FALSE`.
#' @param ignore.s.weights `logical`; whether or not to ignore sampling weights when
#'   computing the weight summary. If `FALSE`, the default, the estimated
#'   weights will be multiplied by the sampling weights (if any) before values
#'   are computed.
#' @param weight.range `logical`; whether to display statistics about the range of weights and the highest and lowest weights for each group. Default is `TRUE`.
#' @param binwidth,bins arguments passed to [ggplot2::geom_histogram()] to
#'   control the size and/or number of bins.
#' @param x a `summary.weightit` or `summary.weightitMSM` object; the output of
#'   a call to `summary.weightit()` or `summary.weightitMSM()`.
#' @param time `numeric`; the time point for which to display the distribution
#'   of weights. Default is to plot the distribution for the first time points.
#' @param ... For `plot()`, additional arguments passed to [graphics::hist()] to
#'   determine the number of bins, though [ggplot2::geom_histogram()] is
#'   actually used to create the plot.
#'
#' @returns
#' For point treatments (i.e., `weightit` objects), `summary()` returns
#' a `summary.weightit` object with the following elements:
#'
#' \item{weight.range}{The range (minimum and maximum) weight for each treatment group.}
#' \item{weight.top}{The units with the greatest weights in each treatment group; how many are included is determined by `top`.}
#' \item{coef.of.var (Coef of Var)}{The coefficient of variation (standard deviation divided by mean) of the weights in each treatment group and overall.}
#' \item{scaled.mad (MAD)}{The mean absolute deviation of the weights in each treatment group and overall divided by the mean of the weights in the corresponding group.}
#' \item{negative entropy (Entropy)}{The negative entropy (\eqn{\sum w log(w)}) of the weights in each treatment group and overall divided by the mean of the weights in the corresponding group.}
#' \item{num.zeros}{The number of weights equal to zero.}
#' \item{effective.sample.size}{The effective sample size for each treatment group before and after weighting. See [ESS()].}
#'
#' For longitudinal treatments (i.e., `weightitMSM` objects), `summary()`
#' returns a list of the above elements for each treatment period.
#'
#' `plot()` returns a `ggplot` object with a histogram displaying the
#' distribution of the estimated weights. If the estimand is the ATT or ATC,
#' only the weights for the non-focal group(s) will be displayed (since the
#' weights for the focal group are all 1). A dotted line is displayed at the
#' mean of the weights.
#'
#' `nobs()` returns a single number. Note that even units with `weights` or
#' `s.weights` of 0 are included.
#'
#' @seealso [weightit()], [weightitMSM()], [summary()]
#'
#' @examples
#' # See example at ?weightit or ?weightitMSM

#' @exportS3Method summary weightit
summary.weightit <- function(object, top = 5L, ignore.s.weights = FALSE, weight.range = TRUE, ...) {

  chk::chk_count(top)
  chk::chk_flag(ignore.s.weights)
  chk::chk_flag(weight.range)

  outnames <- c("weight.range", "weight.top", "weight.mean",
                "coef.of.var", "scaled.mad", "negative.entropy",
                "effective.sample.size")
  out <- make_list(outnames)

  sw <- {
    if (ignore.s.weights || is_null(object$s.weights)) rep.int(1.0, nobs(object))
    else object$s.weights
  }

  w <- object$weights
  t <- object$treat

  treat.type <- get_treat_type(object[["treat"]])
  stabilized <- is_not_null(object[["stabilization"]])

  treat.type[treat.type == "multinomial"] <- "multi-category"

  ww <- setNames(w * sw, seq_along(sw))

  attr(out, "weights") <- ww
  attr(out, "treat") <- t

  if (treat.type == "binary") {
    treated <- get_treated_level(t, object$estimand, object$focal)

    tx <- list(treated = which(t == treated),
               control = which(t != treated))
  }
  else if (treat.type == "multi-category") {
    tx <- lapply(levels(t), function(i) which(t == i)) |>
      setNames(levels(t))
  }
  else {
    tx <- list(all = seq_along(w))
  }

  if (weight.range) {
    out$weight.range <- lapply(tx, function(ti) c(min(ww[ti]), max(w[ti]))) |>
      setNames(names(tx))

    top.weights <- lapply(tx, function(ti) sort(ww[ti], decreasing = TRUE)[seq_len(top)]) |>
      setNames(names(tx))

    out$weight.top <- lapply(names(tx), function(i) {
      sort(setNames(top.weights[[i]], which(ww[tx[[i]]] %in% top.weights[[i]])[seq_len(top)]))
    }) |>
      setNames(names(tx))
  }

  out$coef.of.var <- vapply(tx, function(ti) sd(ww[ti]) / mean_fast(ww[ti]), numeric(1L))
  out$scaled.mad <- vapply(tx, function(ti) mean_abs_dev(ww[ti] / mean_fast(ww[ti])), numeric(1L))
  out$negative.entropy <- vapply(tx, function(ti) neg_ent(ww[ti]), numeric(1L))
  out$num.zeros <- vapply(tx, function(ti) sum(check_if_zero(ww[ti], tol = 1e-10)), numeric(1L))

  if (stabilized) {
    out$weight.mean <- vapply(tx, function(ti) mean_fast(ww[ti]), numeric(1L))
  }

  if (treat.type == "binary") {
    nn <- make_df(c("Control", "Treated"),
                  c("Unweighted", "Weighted"))

    nn[["Control"]] <- c(ESS(sw[tx$control]), ESS(ww[tx$control]))
    nn[["Treated"]] <- c(ESS(sw[tx$treated]), ESS(ww[tx$treated]))
  }
  else if (treat.type == "multi-category") {
    nn <- make_df(levels(t),
                  c("Unweighted", "Weighted"))

    for (i in levels(t)) {
      nn[[i]] <- c(ESS(sw[tx[[i]]]), ESS(ww[tx[[i]]]))
    }
  }
  else {
    nn <- make_df("Total",
                  c("Unweighted", "Weighted"))
    nn[["Total"]] <- c(ESS(sw), ESS(ww))
  }

  out$effective.sample.size <- nn

  attr(w, "focal") <- object$focal %or% NULL

  class(out) <- "summary.weightit"

  out
}

#' @exportS3Method print summary.weightit
print.summary.weightit <- function(x, digits = 3L, ...) {
  cat0(space(18L), .ul("Summary of weights"), "\n\n")

  if (is_not_null(x$weight.range)) {
    cat0("- ", .it("Weight ranges"), ":\n\n")

    x$weight.range |>
      text_box_plot(width = 28L) |>
      round_df_char(digits = digits, pad = " ") |>
      print.data.frame()

    top <- max(lengths(x$weight.top))

    cat0("\n- ", .it(sprintf("Units with the %s most extreme weights%s",
                             top,
                             ngettext(length(x$weight.top), "", " by group"))),
         ":\n")

    data.frame(unlist(lapply(names(x$weight.top), function(y) c(" ", y))),
               matrix(unlist(lapply(x$weight.top, function(y) c(names(y), character(top - length(y)),
                                                                round(y, digits), character(top - length(y))))),
                      byrow = TRUE, nrow = 2 * length(x$weight.top))) |>
      setNames(character(1L + top)) |>
      print.data.frame(row.names = FALSE, digits = digits)

    cat("\n")
  }

  cat0("- ", .it("Weight statistics"), ":\n\n")
  cbind(x$coef.of.var,
        x$scaled.mad,
        x$negative.entropy,
        x$num.zeros) |>
    as.data.frame() |>
    setNames(c("Coef of Var", "MAD", "Entropy", "# Zeros")) |>
    round_df_char(digits = 3L) |>
    print.data.frame()

  if (is_not_null(x$weight.mean)) {
    cat0("\n- ", .it("Mean of Weights"), " = ", round(x$weight.mean, 2L), "\n")
  }

  cat0("\n- ", .it("Effective Sample Sizes"), ":\n\n")
  x$effective.sample.size |>
    round_df_char(digits = 2L, pad = " ") |>
    print.data.frame()

  invisible(x)
}

#' @exportS3Method plot summary.weightit
#' @rdname summary.weightit
plot.summary.weightit <- function(x, binwidth = NULL, bins = NULL, ...) {
  w <- .attr(x, "weights")
  t <- .attr(x, "treat")
  focal <- .attr(w, "focal")

  treat.type <- get_treat_type(t)

  breaks <- ...get("breaks")
  if (is_null(breaks)) {
    bins <- 20
  }
  else {
    breaks <- hist(w, breaks = breaks, plot = FALSE)[["breaks"]]
    bins <- binwidth <- NULL
  }

  subtitle <- {
    if (is_not_null(focal)) sprintf("For Units Not in Treatment Group %s", add_quotes(focal))
    else NULL
  }

  if (treat.type == "continuous") {
    p <- ggplot(data = data.frame(w), mapping = aes(x = .data$w)) +
      geom_histogram(binwidth = binwidth,
                     bins = bins,
                     breaks = breaks,
                     center = mean(w),
                     color = "gray70",
                     fill = "gray70", alpha = 1) +
      scale_y_continuous(expand = expansion(c(0, .05))) +
      geom_vline(xintercept = mean(w), linetype = "12", color = "blue", size = .75) +
      labs(x = "Weight", y = "Count", title = "Distribution of Weights",
           subtitle = subtitle) +
      theme_bw()
  }
  else {
    d <- data.frame(w, t)

    if (is_not_null(focal)) {
      d <- d[t != focal, , drop = FALSE]
    }

    d$t <- factor(d$t)

    levels(d$t) <- sprintf("Treat = %s", levels(d$t))
    w_means <- aggregate(w ~ t, data = d, FUN = mean)

    p <- ggplot(data = d, mapping = aes(x = .data$w)) +
      geom_histogram(binwidth = binwidth,
                     bins = bins,
                     breaks = breaks,
                     color = "gray70",
                     fill = "gray70", alpha = 1) +
      scale_y_continuous(expand = expansion(c(0, .05))) +
      geom_vline(data = w_means, aes(xintercept = .data$w), linetype = "12", color = "red") +
      labs(x = "Weight", y = "Count", title = "Distribution of Weights",
           subtitle = subtitle) +
      theme_bw() +
      facet_wrap(vars(.data$t), ncol = 1L, scales = "free") +
      theme(panel.background = element_blank(),
            panel.border = element_rect(fill = NA, color = "black", size = .25))
  }

  p
}

#' @exportS3Method summary weightitMSM
#' @rdname summary.weightit
summary.weightitMSM <- function(object, top = 5L, ignore.s.weights = FALSE, weight.range = TRUE, ...) {

  chk::chk_count(top)
  chk::chk_flag(ignore.s.weights)
  chk::chk_flag(weight.range)

  out.list <- make_list(names(object$treat.list))

  sw <- {
    if (ignore.s.weights || is_null(object$s.weights)) rep.int(1, nobs(object))
    else object$s.weights
  }

  for (ti in seq_along(object$treat.list)) {
    obj <- as.weightit(object$weights, treat = object$treat.list[[ti]],
                       s.weights = sw, stabilization = object$stabilization)
    out.list[[ti]] <- summary.weightit(obj, top = top, ignore.s.weights = ignore.s.weights, ...)
  }

  class(out.list) <- "summary.weightitMSM"

  out.list
}

#' @exportS3Method print summary.weightitMSM
print.summary.weightitMSM <- function(x, digits = 3L, ...) {
  chk::chk_whole_number(digits)

  only.one <- length(x) == 1L || all_apply(x, function(y) isTRUE(all.equal(x[[1L]], y)))

  for (ti in seq_along(x)) {
    if (!only.one) {
      cat0(.st(space(23L)),
           .it(sprintf(" Time %s ", ti)),
           .st(space(23L)), "\n")
    }

    print(x[[ti]])
    cat("\n")

    if (only.one) {
      break
    }
  }

  invisible(x)
}

#' @exportS3Method plot summary.weightitMSM
#' @rdname summary.weightit
plot.summary.weightitMSM <- function(x, binwidth = NULL, bins = NULL, time = 1L, ...) {
  if (!chk::vld_count(time) || time %nin% seq_along(x)) {
    .err("`time` must be a number corresponding to the time point for which to display the distribution of weights")
  }

  plot.summary.weightit(x[[time]], binwidth = binwidth, bins = bins, ...) +
    labs(subtitle = sprintf("For Time %s", time))
}

#' @exportS3Method nobs weightit
nobs.weightit <- function(object, ...) {
  length(object[["weights"]])
}
