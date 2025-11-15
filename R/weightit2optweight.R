#' Stable Balancing Weights
#' @name method_optweight
#' @aliases method_sbw
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating stable balancing
#' weights (also known as optimization-based weights) by setting `method = "optweight"` in the call to [weightit()]. This method can be used with binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating weights by solving a quadratic
#' programming problem subject to approximate or exact balance constraints. This
#' method relies on \pkgfun{optweight}{optweight.fit} from the \CRANpkg{optweight}
#' package.
#'
#' Because \pkgfun{optweight}{optweight} offers finer control and uses the same syntax as
#' `weightit()`, it is recommended that `optweight()` be used
#' instead of `weightit()` with `method = "optweight"`.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the weights using
#' \pkgfun{optweight}{optweight.fit}. The following estimands are allowed: ATE, ATT,
#' and ATC. The weights are taken from the output of the `optweight.fit` fit object.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the weights using
#' \pkgfun{optweight}{optweight.fit}. The following estimands are allowed: ATE and
#' ATT. The weights are taken from the output of the `optweight.fit` fit object.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, this method estimates the weights using
#' \pkgfun{optweight}{optweight.fit}. The weights are taken from the output of the
#' `optweight.fit` fit object.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights
#' estimated at each time point. This method is not guaranteed to yield exact
#' balance at each time point. **NOTE: the use of stable balancing weights with
#' longitudinal treatments has not been validated and should not be done!**
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios, but only
#' for versions of \pkg{optweight} greater than or equal to 1.0.0.
#'
#' ## Missing Data
#'
#' In the presence of missing data, the following value(s) for `missing` are
#' allowed:
#'     \describe{
#'       \item{`"ind"` (default)}{First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is `NA` and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with the covariate medians (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting `weightit` object will be the original covariates with the `NA`s.
#'       }
#'     }
#'
#' ## M-estimation
#'
#' M-estimation is not supported.
#'
#' @section Additional Arguments:
#'
#' \describe{
#'   \item{`moments`}{`integer`; the highest power of each covariate to be balanced. For example, if `moments = 3`, each covariate, its square, and its cube will be balanced. Can also be a named vector with a value for each covariate (e.g., `moments = c(x1 = 2, x2 = 4)`). Values greater than 1 for categorical covariates are ignored. Default is 1 to balance covariate means.
#'     }
#'     \item{`int`}{`logical`; whether first-order interactions of the covariates are to be balanced. Default is `FALSE`.
#'     }
#'     \item{`quantile`}{a named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures the 25th, 50th, and 75th percentiles of `x1` in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., `.5`) or a vector (e.g., `c(.25, .5, .75)`) to request the same quantile(s) for all continuous covariates. Only allowed with binary and multi-category treatments.
#'     }
#' }
#'
#' All arguments to `optweight.fit()` can be passed through `weightit()` or `weightitMSM()`, with the following exception:
#'
#' * `targets` cannot be used and is ignored.
#'
#' All arguments take on the defaults of those in `optweight.fit()`.
#'
#' @section Additional Outputs:
#'  \describe{
#'  \item{`info`}{
#'  A list with one entry:
#'      \describe{
#'        \item{`duals`}{A data frame of dual variables for each balance constraint.}
#'      }
#'  }
#'  \item{`obj`}{When `include.obj = TRUE`, the output of the call to \pkgfun{optweight}{optweight.fit}.}
#'  }
#'
#' @details
#' Stable balancing weights are weights that solve a constrained
#' optimization problem, where the constraints correspond to covariate balance
#' and the loss function is the variance (or other norm) of the weights. These
#' weights maximize the effective sample size of the weighted sample subject to
#' user-supplied balance constraints. An advantage of this method over entropy
#' balancing is the ability to allow approximate, rather than exact, balance
#' through the `tols` argument, which can increase precision even for slight
#' relaxations of the constraints.
#'
#' The function of the weights that is optimized can be changed using the `norm` argument. The default `norm = "l2"`, minimizes the variance of the weights (i.e., maximizes the ESS). `norm = "entropy"` minimizes the negative entropy of the weights and is equivalent to entropy balancing, though in this implementation, inexact balance is allowed. `norm = "log"` minimizes the sum of the negative logs of the weights and is equivalent to nonparametric covariate balancing propensity score weighting (npCBPS). See \pkgfun{optweight}{optweight.fit} for the other allowed options to `norm` and other arguments.
#'
#' `plot()` can be used on the output of `weightit()` with `method = "optweight"`
#' to display the dual variables; see Examples and [plot.weightit()] for more
#' details.
#'
#' @note
#' The specification of `tols` differs between `weightit()` and
#' `optweight()`. In `weightit()`, one tolerance value should be included per
#' level of each factor variable, whereas in `optweight()`, all levels of a
#' factor are given the same tolerance, and only one value needs to be supplied
#' for a factor variable. Because of the potential for confusion and ambiguity,
#' it is recommended to only supply one value for `tols` in `weightit()` that
#' applies to all variables. For finer control, use `optweight()` directly.
#'
#' Seriously, just use \pkgfun{optweight}{optweight}. The syntax is almost
#' identical and it's compatible with \pkg{cobalt}, too.
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' \pkgfun{optweight}{optweight.fit} for the fitting function.
#'
#' [`method_entropy`] for entropy balancing, which is a special case of stable balancing weights.
#'
#' [`method_npcbps`] for npCBPS weighting, which is also a special case of stable balancing weights.
#'
#' @references
#' ## Binary treatments
#'
#' Wang, Y., & Zubizarreta, J. R. (2020). Minimal dispersion approximately
#' balancing weights: Asymptotic properties and practical considerations.
#' *Biometrika*, 107(1), 93–105. \doi{10.1093/biomet/asz050}
#'
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for
#' Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' ## Multi-Category Treatments
#'
#' de los Angeles Resa, M., & Zubizarreta, J. R. (2020). Direct and Stable Weight Adjustment in Non-Experimental Studies With Multivalued Treatments: Analysis of the Effect of an Earthquake on Post-Traumatic Stress. *Journal of the Royal Statistical Society Series A: Statistics in Society*, 183(4), 1387–1410. \doi{10.1111/rssa.12561}
#'
#' ## Continuous treatments
#'
#' Greifer, N. (2020). *Estimating Balancing Weights for Continuous Treatments Using Constrained Optimization*. \doi{10.17615/DYSS-B342}
#'
#' @examplesIf rlang::is_installed("optweight (>= 1.0.0)")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + race +
#'                   nodegree + re74, data = lalonde,
#'                 method = "optweight", estimand = "ATT",
#'                 tols = 0))
#'
#' summary(W1)
#'
#' cobalt::bal.tab(W1)
#'
#' plot(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "optweight", estimand = "ATE",
#'                 tols = .01))
#'
#' summary(W2)
#'
#' cobalt::bal.tab(W2)
#'
#' plot(W2)
#'
#' #Balancing covariates with respect to re75 (continuous)
#' (W3 <- weightit(re75 ~ age + educ + race +
#'                   nodegree + re74, data = lalonde,
#'                 method = "optweight", tols = .02))
#'
#' summary(W3)
#'
#' cobalt::bal.tab(W3)
#'
#' plot(W3)
NULL

weightit2optweight <- function(covs, treat, s.weights, subset, estimand, focal, missing,
                               verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in c("b.weights", "base.weights", "base.weight")) {
    bw <- ...get(i)

    if (is_not_null(bw)) {
      if (!is.numeric(bw) || length(bw) != length(treat)) {
        .err(sprintf("the argument to `%s` must be a numeric vector with length equal to the number of units",
                     i))
      }

      break
    }
  }

  if (is_null(bw)) {
    bw <- rep_with(1, treat)
  }

  treat <- factor(treat[subset])

  bw <- bw[subset]

  covs <- .apply_moments_int_quantile(covs,
                                      moments = ...get("moments"),
                                      int = ...get("int"),
                                      quantile = ...get("quantile"),
                                      s.weights = s.weights, focal = focal,
                                      treat = treat)

  A <- list(...)

  A[["covs"]] <- covs
  A[["treat"]] <- treat
  A[["estimand"]] <- estimand
  A[["s.weights"]] <- s.weights
  A[["b.weights"]] <- bw
  A[["focal"]] <- focal
  A[["verbose"]] <- TRUE

  if (is_not_null(A[["targets"]])) {
    .wrn("`targets` cannot be used through `WeightIt` and will be ignored")
    A[["targets"]] <- NULL
  }

  verbosely({
    out <- do.call(optweight::optweight.fit, A, quote = TRUE)
  }, verbose = verbose)

  list(w = out[["w"]],
       info = list(duals = .process_duals(out[["duals"]], covs)),
       fit.obj = out)
}

weightit2optweight.multi <- weightit2optweight

weightit2optweight.cont <- function(covs, treat, s.weights, subset, missing, verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in c("b.weights", "base.weights", "base.weight")) {
    bw <- ...get(i)

    if (is_not_null(bw)) {
      if (!is.numeric(bw) || length(bw) != length(treat)) {
        .err(sprintf("the argument to `%s` must be a numeric vector with length equal to the number of units",
                     i))
      }

      break
    }
  }

  if (is_null(bw)) {
    bw <- rep_with(1, treat)
  }

  treat <- treat[subset]

  bw <- bw[subset]

  covs <- covs |>
    .apply_moments_int_quantile(moments = ...get("moments"),
                                int = ...get("int")) |>
    .make_covs_closer_to_1()


  A <- list(...)

  A[["covs"]] <- covs
  A[["treat"]] <- treat
  A[["estimand"]] <- "ATE"
  A[["s.weights"]] <- s.weights
  A[["b.weights"]] <- bw
  A[["verbose"]] <- TRUE

  if (is_not_null(A[["targets"]])) {
    .wrn("`targets` cannot be used through `WeightIt` and will be ignored")
    A[["targets"]] <- NULL
  }

  verbosely({
    out <- do.call(optweight::optweight.fit, A, quote = TRUE)
  }, verbose = verbose)

  list(w = out[["w"]],
       info = list(duals = .process_duals(out[["duals"]], covs)),
       fit.obj = out)
}

.plot_duals_optweight <- function(info, by = NULL) {

  use.by <- is_not_null(by)

  if (!use.by) {
    info <- list("z" = info)
  }

  d <- do.call("rbind", lapply(seq_along(info), function(i) {
    cbind(info[[i]]$duals, by = names(info)[i])
  }))

  d$by <- factor(d$by, levels = names(info))

  d$cov <- factor(d$cov, levels = rev(unique(d$cov)))

  constraint_types <- unique(d$constraint, nmax = 2L)

  d$constraint <- factor(d$constraint,
                         levels = constraint_types,
                         labels = paste("Constraint:", constraint_types))

  title <- "Dual Variables for Constraints"

  ggplot(d) +
    geom_col(aes(y = .data$cov, x = .data$dual)) +
    scale_x_continuous(expand = expansion(c(0, .05))) +
    facet_grid(rows = vars(.data$constraint),
               cols = if (use.by) vars(.data$by),
               scales = "free_y", space = "free") +
    labs(x = "Absolute Dual Variable",
         y = "Variable",
         title = title) +
    theme_bw()
}

.process_duals <- function(d, covs) {

  if (is_null(d$dual)) {
    return(NULL)
  }

  na.cov <- is.na(d$cov)

  if (is_not_null(d$treat) && any(na.cov)) {
    d$cov[na.cov] <- d$treat[na.cov]
  }

  d$dual <- ave(d$dual, d$constraint, d$cov, FUN = sum) #Total effect of constraint on obj. fun. is sum of abs(duals)

  d$treat <- NULL

  d <- unique(d)

  rownames(d) <- NULL

  d
}
