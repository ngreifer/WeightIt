#' Plot information about the weight estimation process
#' @name plot.weightit
#'
#' @description
#' `plot.weightit()` plots information about the weights depending
#' on how they were estimated. Currently, only weighting using `method = "gbm"`
#' or `"optweight"` are supported. To plot the distribution of weights, see
#' [plot.summary.weightit()].
#'
#' @param x a `weightit` object; the output of a call to [weightit()].
#' @param ... unused.
#'
#' @returns
#' A `ggplot` object.
#'
#' @details
#'
#' ## `method = "gbm"`
#'
#' After weighting with generalized boosted modeling, `plot()` displays the
#' results of the tuning process used to find the optimal number of trees (and
#' tuning parameter values, if modified) that are used in the final weights. The
#' plot produced has the number of trees on the x-axis and the value of the
#' criterion on the y-axis with a diamond at the optimal point. When multiple
#' parameters are selected by tuning, a separate line is displayed on the plot
#' for each combination of tuning parameters. When `by` is used in the call to
#' `weightit()`, the plot is faceted by the `by` variable. See [`method_gbm`]
#' for more information on selecting tuning parameters.
#'
#' ## `method = "optweight"`
#'
#' After estimating stable balancing weights, `plot()` displays the values of
#' the dual variables for each balance constraint in a bar graph. Large values
#' of the dual variables indicate the covariates for which the balance
#' constraint is causing increases in the variability of the weights, i.e., the
#' covariates for which relaxing the imbalance tolerance would yield the
#' greatest gains in effective sample size. For continuous treatments, the dual
#' variables are split into those for the target (i.e., ensuring the mean of
#' each covariate after weighting is equal to its unweighted mean) and those for
#' balance (i.e., ensuring the treatment-covariate correlations are no larger
#' than the imbalance tolerance). This is essentially a wrapper for
#' \pkgfun{optweight}{plot.optweight}. See [`method_optweight`] for details.
#'
#' @seealso [weightit()], [plot.summary.weightit()]
#'
#' @examples
#' # See example at the corresponding methods page

#' @exportS3Method plot weightit
plot.weightit <- function(x, ...) {
  if (inherits(x, "weightitMSM")) {
    .err("`plot(.)` can currently only be used with `weightit()` output objects. To view the distribution of weights, use `plot(summary(.))`")
  }

  if (!chk::vld_string(x$method) || !.weightit_methods[[x$method]]$plot.weightit_ok) {
    .err(sprintf("`plot(.)` cannot be used with %s. To view the distribution of weights, use `plot(summary(.))`",
                 .method_to_phrase(x$method)))
  }

  switch(x$method,
         gbm = .plot_tune_gbm(x$info, x$by),
         optweight = .plot_duals_optweight(x$info, x$by))
}
