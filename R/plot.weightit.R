#' Plot information about the weight estimation process
#' @name plot.weightit
#'
#' @description
#' `plot.weightit()` plots information about the weights depending on how they were estimated. Currently, only weighting using `method = "gbm"` is supported. To plot the distribution of weights, see [plot.summary.weightit()].
#'
#' @param x a `weightit` object; the output of
#' a call to [weightit()].
#' @param ... Unused.
#'
#' @returns
#' A `ggplot` object.
#'
#' @details
#'
#' ## `method = "gbm"`
#'
#' After weighting with generalized boosted modeling, `plot()` displays the results of the tuning process used to find the optimal number of trees (and tuning parameter values, if modified) that are used in the final weights. The plot produced has the number of trees on the x-axis and the value of the criterion on the y axis with a diamond at the optimal point. When multiple parameters are selected by tuning, a separate line is displayed on the plot for each combination of tuning parameters. When `by` is used in the call to `weightit()`, the plot is faceted by the `by` variable. See [`method_gbm`] for more information on selecting tuning parameters.
#'
#'
#' @seealso
#' [weightit()], [plot.summary.weightit()]
#'
#' @examples
#'
#' # See example at the corresponding methods page
#'

#' @exportS3Method plot weightit
plot.weightit <- function(x, ...) {
  if (inherits(x, "weightitMSM")) {
    .err("`plot(.)` can currently only be used with `weightit()` output objects. To view the distribution of weights, use `plot(summary(.))`")
  }

  if (x$method != "gbm") {
    .err("`plot(.)` can currently only be used with `method = \"gbm\"`. To view the distribution of weights, use `plot(summary(.))`")
  }

  .plot_tune_gbm(x$info, x$by)
}