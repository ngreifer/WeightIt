#' Nonparametric Covariate Balancing Propensity Score Weighting
#' @name method_npcbps
#' @aliases method_npcbps
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights from nonparametric covariate balancing propensity scores by setting `method = "npcbps"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating weights by maximizing the empirical likelihood of the data subject to balance constraints. This method relies on \pkgfun{CBPS}{npCBPS} from the \CRANpkg{CBPS} package.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the weights using \pkgfun{CBPS}{npCBPS}. The ATE is the only estimand allowed. The weights are taken from the output of the `npCBPS` fit object.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the weights using \pkgfun{CBPS}{npCBPS}. The ATE is the only estimand allowed. The weights are taken from the output of the `npCBPS` fit object.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, this method estimates the weights using \pkgfun{CBPS}{npCBPS}. The weights are taken from the output of the `npCBPS` fit object.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights estimated at each time point. This is not how \pkgfun{CBPS}{CBMSM} estimates weights for longitudinal treatments.
#'
#' ## Sampling Weights
#'
#' Sampling weights are \bold{not} supported with `method = "npcbps"`.
#'
#' ## Missing Data
#'
#' In the presence of missing data, the following value(s) for `missing` are allowed:
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
#' `moments` and `int` are accepted. See [weightit()] for details.
#'
#' \describe{
#'   \item{`quantile`}{
#'     A named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures the 25th, 50th, and 75th percentiles of `x1` in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., `.5`) or an unnamed list of length 1 (e.g., `list(c(.25, .5, .75))`) to request the same quantile(s) for all continuous covariates, or a named vector (e.g., `c(x1 = .5, x2 = .75`) to request one quantile for each covariate. Only allowed with binary and multi-category treatments.
#'   }
#' }
#'
#' All arguments to `npCBPS()` can be passed through `weightit()` or `weightitMSM()`.
#'
#' All arguments take on the defaults of those in `npCBPS()`.
#'
#' @section Additional Outputs:
#' \describe{
#'   \item{`obj`}{When `include.obj = TRUE`, the nonparametric CB(G)PS model fit. The output of the call to \pkgfun{CBPS}{npCBPS}.
#'   }
#' }
#'
#' @details
#' Nonparametric CBPS involves the specification of a constrained optimization problem over the weights. The constraints correspond to covariate balance, and the loss function is the empirical likelihood of the data given the weights. npCBPS is similar to \link[=method_ebal]{entropy balancing} and will generally produce similar results. Because the optimization problem of npCBPS is not convex it can be slow to converge or not converge at all, so approximate balance is allowed instead using the `cor.prior` argument, which controls the average deviation from zero correlation between the treatment and covariates allowed.
#'
#' @note
#' When sampling weights are used with `CBPS::CBPS()`, the estimated weights already incorporate the sampling weights. When `weightit()` is used with `method = "cbps"`, the estimated weights are separated from the sampling weights, as they are with all other methods.
#'
#' @seealso
#' [weightit()], [weightitMSM()], [`method_cbps`]
#'
#' \pkgfun{CBPS}{npCBPS} for the fitting function
#'
#' @references
#' Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity score for a continuous treatment: Application to the efficacy of political advertisements. *The Annals of Applied Statistics*, 12(1), 156–177. \doi{10.1214/17-AOAS1101}
#'
#' @examplesIf requireNamespace("CBPS", quietly = TRUE)
#' # Examples take a long time to run
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#' \donttest{
#'   #Balancing covariates between treatment groups (binary)
#'   (W1 <- weightit(treat ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "npcbps", estimand = "ATE"))
#'   summary(W1)
#'   bal.tab(W1)
#'
#'   #Balancing covariates with respect to race (multi-category)
#'   (W2 <- weightit(race ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "npcbps", estimand = "ATE"))
#'   summary(W2)
#'   bal.tab(W2)
#' }
NULL

weightit2npcbps <- function(covs, treat, s.weights, subset, missing, moments, int, verbose, ...) {
  rlang::check_installed("CBPS")

  A <- list(...)

  if (!all_the_same(s.weights)) {
    .err("sampling weights cannot be used with `method = \"npcbps\"`")
  }

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- cbind(covs, .int_poly_f(covs, poly = moments, int = int))

  covs <- cbind(covs, .quantile_f(covs, qu = A[["quantile"]], s.weights = s.weights))

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  tryCatch({verbosely({
    fit <- do.call(CBPS::npCBPS, c(list(formula(new.data), data = new.data, print.level = 1), A),
                   quote = TRUE)
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    .err("(from `CBPS::npCBPS()`) ", e., tidy = FALSE)
  })

  w <- fit$weights

  for (i in levels(treat)) w[treat == i] <- w[treat == i]/mean(w[treat == i])

  list(w = w, fit.obj = fit)
}

weightit2npcbps.multi <- weightit2npcbps

weightit2npcbps.cont <- function(covs, treat, s.weights, subset, missing, moments, int, verbose, ...) {
  rlang::check_installed("CBPS")

  A <- list(...)

  if (!all_the_same(s.weights)) {
    .err("sampling weights cannot be used with `method = \"npcbps\"`")
  }

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  covs <- cbind(covs, .int_poly_f(covs, poly = moments, int = int))

  colinear.covs.to.remove <- colnames(covs)[colnames(covs) %nin% colnames(make_full_rank(covs))]
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  tryCatch({verbosely({
    fit <- do.call(CBPS::npCBPS, c(list(formula(new.data), data = new.data, print.level = 1), A),
                   quote = TRUE)
  }, verbose = verbose)},
  error = function(e) {
    e. <- conditionMessage(e)
    .err("(from `CBPS::npCBPS()`) ", e., tidy = FALSE)
  })

  w <- fit$weights

  w <- w/mean(w)

  list(w = w, fit.obj = fit)
}
