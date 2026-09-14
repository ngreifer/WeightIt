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
#' @param which.time for `summary()`, which models to summarize, given as a vector
#'   of positions in `formula.list` or of treatment or censoring variable names;
#'   omit it to summarize all of them. For `plot()`, which single model to display
#'   the distribution of weights for, given as one position or one name; default
#'   is the first. This mirrors the argument of the same name in
#'   \pkgfun{cobalt}{bal.tab}, except that a value matching no model is an error
#'   rather than a warning, since here it decides what is computed. Note that when
#'   censoring is modeled, the censoring models occupy positions of their own, so
#'   positions do not count treatment time points alone; naming the variable
#'   avoids having to count.
#' @param time for `plot()`, the former name of `which.time`, which still works.
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
#' \item{negative entropy (Entropy)}{The negative entropy (\eqn{\frac{1}{n}\sum w \log(w)}) of the weights in each treatment group and overall, after dividing the weights by their mean in the corresponding group.}
#' \item{weight.mean (Mean of Weights)}{The mean of the weights in each treatment group and overall. Only included when the weights are stabilized.}
#' \item{num.zeros}{The number of weights equal to zero.}
#' \item{effective.sample.size}{The effective sample size for each treatment group before and after weighting. See [ESS()].}
#'
#' For longitudinal treatments (i.e., `weightitMSM` objects), `summary()`
#' returns a list of the above elements for each treatment period. When censoring
#' is modeled (see [.cens()]), each censoring model gets an entry of its own,
#' placed among the treatment entries in the order the models were fit and named
#' for its censoring indicator. Every entry summarizes the same weights -- the
#' product across all the models -- and differs only in the sample it summarizes
#' them over: a treatment entry splits by treatment group, and a censoring entry
#' covers the units still under observation when that model was fit, since the
#' censored units have a weight of exactly 0 and are not part of the weighted
#' sample.
#'
#' `plot()` returns a `ggplot` object with a histogram displaying the
#' distribution of the estimated weights. For a censoring model, only the weights
#' of the units still under observation are displayed, matching what `summary()`
#' reports. If the estimand is the ATT or ATC,
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

  arg::arg_count(top)
  arg::arg_flag(ignore.s.weights)
  arg::arg_flag(weight.range)

  outnames <- c("weight.range", "weight.top", "weight.mean",
                "coef.of.var", "scaled.mad", "negative.entropy",
                "num.zeros", "effective.sample.size")
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

  ww <- setNames(w * sw, names(t) %or% seq_along(t))

  if (treat.type == "binary") {
    treated <- .get_treated_level(t, object$estimand, object$focal)

    tx <- list("Treated" = which(t == treated),
               "Control" = which(t != treated))
  }
  else if (treat.type == "multi-category") {
    tx <- lapply(levels(t), function(i) which(t == i)) |>
      setNames(levels(t))
  }
  else if (treat.type == "censoring") {
    #Censored units have a weight of exactly 0 and are not part of the weighted
    #sample, so the weight distribution is summarized over the units still under
    #observation.
    C <- .make_cens_treat(t)

    tx <- list(All = which(C == 0))

    #The effective sample size, though, is measured against the units this model was
    #fit on -- those under observation entering it, censored here or not -- so that it
    #shows what this censoring cost. A unit censored at an earlier time point was never
    #eligible and has no indicator here; counting it would understate the weighted
    #sample's efficiency relative to a sample it was never drawn from.
    at.risk <- which(!is.na(C))
  }
  else {
    tx <- list(All = seq_along(w))
  }

  if (weight.range) {
    out$weight.range <- lapply(tx, function(ti) c(min(ww[ti]), max(ww[ti]))) |>
      setNames(names(tx))

    out$weight.top <- lapply(tx, function(ti) rev(sort(ww[ti], decreasing = TRUE)[seq_len(top)])) |>
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

    nn[["Control"]] <- c(ESS(sw[tx$Control]), ESS(ww[tx$Control]))
    nn[["Treated"]] <- c(ESS(sw[tx$Treated]), ESS(ww[tx$Treated]))
  }
  else if (treat.type == "multi-category") {
    nn <- make_df(levels(t),
                  c("Unweighted", "Weighted"))

    for (i in levels(t)) {
      nn[[i]] <- c(ESS(sw[tx[[i]]]), ESS(ww[tx[[i]]]))
    }
  }
  else if (treat.type == "censoring") {
    nn <- make_df("Total",
                  c("Unweighted", "Weighted"))
    nn[["Total"]] <- c(ESS(sw[at.risk]), ESS(ww[at.risk]))
  }
  else {
    nn <- make_df("Total",
                  c("Unweighted", "Weighted"))
    nn[["Total"]] <- c(ESS(sw), ESS(ww))
  }

  out$effective.sample.size <- nn

  attr(ww, "focal") <- object$focal %or% NULL
  attr(out, "weights") <- ww
  attr(out, "treat") <- t

  class(out) <- "summary.weightit"

  out
}

#' @exportS3Method print summary.weightit
print.summary.weightit <- function(x, digits = 3L, ...) {
  arg::arg_whole_number(digits)

  cli::cat_line(space(18L), .ul("Summary of weights"), "\n")

  .print_summary_weightit_internal(x, digits = digits, ...)
}

.print_summary_weightit_internal <- function(x, digits = 3L, ...) {
  bullet <- "line"

  if (is_not_null(x$weight.range)) {
    cli::cat_bullet(.it("Weight ranges"), ":\n", bullet = bullet)
    x$weight.range |>
      .text_box_plot(width = 28L) |>
      round_df_char(digits = digits, pad = " ") |>
      print.data.frame()
    cat("\n")

    top <- max(lengths(x$weight.top))
    label <- sprintf("Units with the %s most extreme weights%s",
                      top,
                      ngettext(length(x$weight.top), "", " by group"))

    cli::cat_bullet(.it(label), ":", bullet = bullet)

    data.frame(unlist(lapply(names(x$weight.top), function(y) c(" ", y))),
               matrix(unlist(lapply(x$weight.top, function(y) c(names(y), character(top - length(y)),
                                                                round(y, digits), character(top - length(y))))),
                      byrow = TRUE, nrow = 2 * length(x$weight.top))) |>
      setNames(character(1L + top)) |>
      print.data.frame(row.names = FALSE, digits = digits)

    cat("\n")
  }

  cli::cat_bullet(.it("Weight statistics"), ":\n", bullet = bullet)
  cbind(x$coef.of.var,
        x$scaled.mad,
        x$negative.entropy,
        x$num.zeros) |>
    as.data.frame() |>
    setNames(c("Coef of Var", "MAD", "Entropy", "# Zeros")) |>
    round_df_char(digits = 3L, pad = " ") |>
    print.data.frame()
  cat("\n")

  if (is_not_null(x$weight.mean)) {
    cli::cat_bullet(.it("Mean of Weights"), ":", bullet = bullet)
    x$weight.mean |>
      as.data.frame() |>
      setNames("") |>
      round_df_char(digits = digits, pad = " ") |>
      print.data.frame()
    cat("\n")
  }

  cli::cat_bullet(.it("Effective Sample Sizes"), ":\n", bullet = bullet)
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
    #Only fall back to a bin count when the user has specified neither `bins` nor
    #`binwidth`; supplying both to `geom_histogram()` makes it ignore one of them.
    if (is_null(binwidth)) {
      bins <- bins %or% 20
    }
  }
  else {
    breaks <- hist(w, breaks = breaks, plot = FALSE)[["breaks"]]
    bins <- binwidth <- NULL
  }

  subtitle <- if (is_not_null(focal)) sprintf("For Units Not in Treatment Group %s", add_quotes(focal))

  if (treat.type == "censoring") {
    #The censored units have a weight of exactly 0 and are not part of the weighted
    #sample, so the distribution shown is that of the units still under observation,
    #matching what `summary()` reports. Faceting by the indicator would instead show a
    #spike at 0 next to the weights of interest, labeled as treatment groups.
    w <- w[which(.make_cens_treat(t) == 0)]

    subtitle <- "For Units Still Under Observation"

    treat.type <- "continuous"
  }

  if (treat.type == "continuous") {
    p <- ggplot(data = data.frame(w), mapping = aes(x = .data$w)) +
      geom_histogram(binwidth = binwidth,
                     bins = bins,
                     breaks = breaks,
                     boundary = 0,
                     color = "gray70",
                     fill = "gray70", alpha = 1) +
      scale_y_continuous(expand = expansion(c(0, .05))) +
      geom_vline(xintercept = mean(w), linetype = "12", color = "blue",
                 linewidth = .75) +
      labs(x = "Weight", y = "Count", title = "Distribution of Weights",
           subtitle = subtitle) +
      theme_bw()
  }
  else {
    d <- data.frame(w, t)

    if (is_not_null(focal)) {
      d <- d[which(t != focal), , drop = FALSE]
    }

    #A unit censored at an earlier time point has no treatment here and is not part of
    #the weighted sample; `summary()` leaves it out, and without this it would get a
    #facet of its own holding nothing but weights of 0.
    d <- d[!is.na(d$t), , drop = FALSE]

    d$t <- factor(d$t)

    levels(d$t) <- sprintf("Treat = %s", levels(d$t))
    w_means <- aggregate(w ~ t, data = d, FUN = mean)

    p <- ggplot(data = d, mapping = aes(x = .data$w)) +
      geom_histogram(binwidth = binwidth,
                     bins = bins,
                     breaks = breaks,
                     boundary = 0,
                     color = "gray70",
                     fill = "gray70", alpha = 1) +
      scale_y_continuous(expand = expansion(c(0, .05))) +
      geom_vline(data = w_means, aes(xintercept = .data$w), linetype = "12", color = "red") +
      labs(x = "Weight", y = "Count", title = "Distribution of Weights",
           subtitle = subtitle) +
      theme_bw() +
      facet_wrap(vars(.data$t), ncol = 1L, scales = "free") +
      theme(panel.background = element_blank(),
            panel.border = element_rect(fill = NA, color = "black",
                                        linewidth = .25))
  }

  p
}

#' @exportS3Method summary weightitMSM
#' @rdname summary.weightit
summary.weightitMSM <- function(object, top = 5L, ignore.s.weights = FALSE, weight.range = TRUE,
                                which.time, ...) {

  arg::arg_count(top)
  arg::arg_flag(ignore.s.weights)
  arg::arg_flag(weight.range)

  #Treatment and censoring models in the order they were fit. The weights are the
  #product across all of them, so every entry summarizes the same weights and differs
  #only in the sample it summarizes them over: a treatment entry splits by treatment
  #group, a censoring entry covers the units still under observation when that model
  #was fit.
  treat.list <- object[["treat.list"]]

  n.models <- length(treat.list)

  #Computed over every model and subset afterwards, so that a treatment keeps the time
  #number it has in the whole sequence however few of them are asked for.
  labels <- .msm_labels(names(treat.list), is_cens_treat(treat.list))

  keep <- {
    if (missing(which.time)) seq_along(treat.list)
    else .process_which.time(which.time, names(treat.list))
  }

  treat.list <- treat.list[keep]
  labels <- labels[keep]

  out.list <- make_list(names(treat.list))

  sw <- {
    if (ignore.s.weights || is_null(object$s.weights)) rep.int(1, nobs(object))
    else object$s.weights
  }

  for (ti in seq_along(treat.list)) {
    obj <- as.weightit(object$weights, treat = treat.list[[ti]],
                       s.weights = sw, stabilization = object$stabilization)
    out.list[[ti]] <- summary.weightit(obj, top = top, ignore.s.weights = ignore.s.weights,
                                       weight.range = weight.range, ...)
  }

  attr(out.list, "labels") <- labels

  #A subset is always labeled when printed: which time point is being shown cannot be
  #inferred from a summary standing on its own.
  attr(out.list, "subset") <- length(keep) < n.models

  class(out.list) <- "summary.weightitMSM"

  out.list
}

#`which.time` as `cobalt::bal.tab()` understands it, less its `.all` and `.none`:
#requesting no models has no use here, and requesting all of them is what leaving the
#argument out already does. Returns the positions to keep. Unlike `bal.tab()`, where
#`which.time` only chooses what to display and an unusable value costs nothing, here it
#decides what is computed, so a value matching no model is an error rather than a
#warning.
.process_which.time <- function(which.time, nm) {
  if (is_null(which.time)) {
    return(seq_along(nm))
  }

  if (is.numeric(which.time)) {
    out <- seq_along(nm)[seq_along(nm) %in% which.time]

    if (is_null(out)) {
      arg::err("no numbers in {.arg which.time} correspond to a model; there {?is/are} only {length(nm)}")
    }

    return(out)
  }

  if (is.character(which.time)) {
    out <- seq_along(nm)[nm %in% which.time]

    if (is_null(out)) {
      arg::err(c("no names in {.arg which.time} correspond to a model.",
                 "i" = "Available: {.or {.val {nm}}}"))
    }

    return(out)
  }

  arg::err("{.arg which.time} must be a vector of model numbers or of treatment or censoring variable names")
}

#What each entry of a `summary.weightitMSM` is called: its position in `formula.list`,
#whether it is a treatment or a censoring model, and the variable it is about. Numbered
#by position rather than by treatment time point so that the number is the one
#`which.time` takes, and so that a censoring model has a number too.
.msm_labels <- function(nm, is.cens) {
  sprintf("%s. %s: %s",
          seq_along(nm),
          ifelse(is.cens, "Censoring", "Treatment"),
          nm)
}

#The labels a `summary.weightitMSM` was built with, which record each entry's position
#in the whole sequence rather than in whatever subset is in hand. Derived from the
#summaries themselves for an object made before they were stored, where nothing was
#subset and so the positions are the same either way.
.summary_msm_labels <- function(x) {
  .attr(x, "labels") %or% {
    .msm_labels(names(x), is_cens_treat(lapply(x, .attr, "treat")))
  }
}

#' @exportS3Method print summary.weightitMSM
print.summary.weightitMSM <- function(x, digits = 3L, ...) {
  arg::arg_whole_number(digits)

  cli::cat_line(space(18L), .ul("Summary of weights"), "\n")

  .cat_segments(.summary_msm_labels(x),
                function(i) .print_summary_weightit_internal(x[[i]], digits = digits, ...))

  invisible(x)
}

#' @exportS3Method plot summary.weightitMSM
#' @rdname summary.weightit
plot.summary.weightitMSM <- function(x, binwidth = NULL, bins = NULL, which.time = 1L, ...,
                                     time) {
  #`time` was this argument's name before it was aligned with `which.time` in
  #`summary()` and in `cobalt::bal.tab()`. It still works, silently.
  if (!missing(time)) {
    which.time <- time
  }

  labels <- .summary_msm_labels(x)

  #A censoring model gets an entry of its own, so the entries are no longer numbered by
  #treatment time point alone; naming the treatment or censoring variable picks one out
  #without having to count.
  if (rlang::is_string(which.time)) {
    i <- match(which.time, names(x))

    if (is.na(i)) {
      arg::err("{.arg which.time} must be {.or {.val {names(x)}}} or the position of one of them")
    }

    which.time <- i
  }
  else if (!is_number(which.time) || which.time %nin% seq_along(x)) {
    arg::err("{.arg which.time} must be a single number corresponding to the time point for which to display the distribution of weights, or the name of a treatment or censoring variable")
  }

  plot.summary.weightit(x[[which.time]], binwidth = binwidth, bins = bins, ...) +
    labs(subtitle = labels[which.time])
}

#' @exportS3Method nobs weightit
nobs.weightit <- function(object, ...) {
  length(object[["weights"]])
}

#The narrowest a segment divider is ever drawn, whatever sits under it.
.RULE_WIDTH <- 15L

.cat_rule <- function(label, width = 0L) {
  width <- max(.RULE_WIDTH, width, nchar(label) + 5L)
  pad <- width - nchar(label) - 2L

  cli::cat_line(strrep(cli::symbol$line, 3L),
                " ",
                .it(.bd(label)),
                " ",
                strrep(cli::symbol$line, pad - 3L), "\n")
}
.capture_printed <- function(expr) {
  rlang::with_options(utils::capture.output(invisible(expr)),
                      cli.num_colors = cli::num_ansi_colors())
}
.printed_width <- function(out) {
  max(0L, cli::ansi_nchar(out))
}
.cat_segments <- function(labels, print.one, .close = FALSE) {
  out <- lapply(seq_along(labels), function(i) .capture_printed(print.one(i)))

  #Settled here rather than left to `.cat_rule()` so that the closing rule below is drawn
  #to the same width as the dividers above it.
  width <- max(.RULE_WIDTH, vapply(out, .printed_width, numeric(1L)), nchar(labels) + 5L)

  for (i in seq_along(labels)) {
    cli::cat_line()
    .cat_rule(labels[i], width)
    writeLines(out[[i]])
  }

  if (.close) {
    cli::cat_line(strrep(cli::symbol$line, width))
  }

  invisible(width)
}
.text_box_plot <- function(range.list, width = 12L) {
  #The span has a cli symbol and so degrades on its own wherever Unicode is unusable
  #(a non-UTF-8 console, `cli.unicode = FALSE`, which is what `testthat` sets). The
  #brackets and the point marker have none, so they follow the same test rather than
  #staying Unicode in output that has gone ASCII everywhere else.
  dashchar <- cli::symbol$double_line

  utf8 <- cli::is_utf8_output()

  lborder <- if (utf8) "\u255E" else "|"
  rborder <- if (utf8) "\u2561" else "|"
  middle <- if (utf8) "\u2502" else "|"

  full.range <- range(unlist(range.list))
  if (all_the_same(full.range)) {
    for (i in seq_along(range.list)) {
      range.list[[i]][1L] <- range.list[[i]][1L] - 1e-6
      range.list[[i]][2L] <- range.list[[i]][2L] + 1e-6
    }
    full.range <- range(unlist(range.list))
  }
  ratio <- diff(full.range) / (width + 1)
  rescaled.range.list <- lapply(range.list, function(x) round(x / ratio))
  rescaled.full.range <- round(full.range / ratio)
  d <- make_df(c("Min", space(width + 1L), "Max"),
               names(range.list),
               "character")
  d[["Min"]] <- vapply(range.list, function(x) x[1L], numeric(1L))
  d[["Max"]] <- vapply(range.list, function(x) x[2L], numeric(1L))
  for (i in seq_row(d)) {
    spaces1 <- rescaled.range.list[[i]][1L] - rescaled.full.range[1L]
    if (diff(rescaled.range.list[[i]]) == 0) {
      dashes <- max(c(0L, diff(rescaled.range.list[[i]]) - 2L))
      spaces2 <- max(c(0L, diff(rescaled.full.range) - (spaces1 + 3L)))

      d[i, 2L] <- sprintf("%s%s%s",
                          space(spaces1),
                          middle,
                          space(spaces2))
    }
    else {
      dashes <- max(c(0L, diff(rescaled.range.list[[i]]) - 2L))
      spaces2 <- max(c(0L, diff(rescaled.full.range) - (spaces1 + 1L + dashes + 1L)))

      d[i, 2L] <- sprintf("%s%s%s%s%s",
                          space(spaces1),
                          lborder,
                          strrep(dashchar, dashes),
                          rborder,
                          space(spaces2))
    }
  }

  d
}
