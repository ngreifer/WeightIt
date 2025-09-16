#' Nonparametric Covariate Balancing Propensity Score Weighting
#' @name method_npcbps
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights from
#' nonparametric covariate balancing propensity scores by setting `method = "npcbps"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating weights by maximizing the
#' empirical likelihood of the data subject to balance constraints. This method
#' relies on \pkgfun{CBPS}{npCBPS} from the \CRANpkg{CBPS} package.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the weights using
#' \pkgfun{CBPS}{npCBPS}. The ATE is the only estimand allowed. The weights are
#' taken from the output of the `npCBPS` fit object.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the weights using
#' \pkgfun{CBPS}{npCBPS}. The ATE is the only estimand allowed. The weights are
#' taken from the output of the `npCBPS` fit object.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, this method estimates the weights using
#' \pkgfun{CBPS}{npCBPS}. The weights are taken from the output of the `npCBPS`
#' fit object.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights
#' estimated at each time point. **NOTE: the use of npCBPS with longitudinal treatments has not been validated!**
#'
#' ## Sampling Weights
#'
#' Sampling weights are \bold{not} supported with `method = "npcbps"`.
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
#' Nonparametric CBPS involves the specification of a constrained
#' optimization problem over the weights. The constraints correspond to
#' covariate balance, and the loss function is the empirical likelihood of the
#' data given the weights. npCBPS is similar to \link[=method_ebal]{entropy balancing} and will generally produce similar results. Because the optimization problem of npCBPS is not convex it can be slow to converge or
#' not converge at all, so approximate balance is allowed instead using the
#' `cor.prior` argument, which controls the average deviation from zero
#' correlation between the treatment and covariates allowed.
#'
#' @seealso
#' [weightit()], [weightitMSM()], [`method_cbps`]
#'
#' [`method_optweight`], which can also be used to perform npCBPS by setting `norm = "log"`. In generally, this `"optweight"` implementation is more stable and flexible.
#'
#' \pkgfun{CBPS}{npCBPS} for the fitting function
#'
#' @references
#' Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing
#' propensity score for a continuous treatment: Application to the efficacy of
#' political advertisements. *The Annals of Applied Statistics*, 12(1), 156–177.
#' \doi{10.1214/17-AOAS1101}
#'
#' @examplesIf rlang::is_installed("CBPS")
#' # Examples take a long time to run
#' data("lalonde", package = "cobalt")
#' \donttest{
#'   #Balancing covariates between treatment groups (binary)
#'   (W1 <- weightit(treat ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "npcbps", estimand = "ATE"))
#'
#'   summary(W1)
#'
#'   cobalt::bal.tab(W1)
#'
#'   #Balancing covariates with respect to race (multi-category)
#'   (W2 <- weightit(race ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "npcbps", estimand = "ATE"))
#'
#'   summary(W2)
#'
#'   cobalt::bal.tab(W2)
#' }
NULL

weightit2npcbps <- function(covs, treat, s.weights, subset, missing, verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  covs <- .apply_moments_int_quantile(covs,
                                      moments = ...get("moments"),
                                      int = ...get("int"),
                                      quantile = ...get("quantile"),
                                      s.weights = s.weights,
                                      treat = treat)

  for (i in seq_col(covs)) {
    covs[, i] <- .make_closer_to_1(covs[, i])
  }

  colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  corprior <- ...get("corprior", .01)

  tryCatch({verbosely({
    fit <- CBPS::npCBPS(formula(new.data),
                        data = new.data,
                        corprior = corprior,
                        print.level = 1)
  }, verbose = verbose)},
  error = function(e) {
    .err(sprintf("(from `CBPS::npCBPS()`): %s",
                 conditionMessage(e)),
         tidy = FALSE)
  })

  w <- fit$weights

  for (i in levels(treat)) {
    w[treat == i] <- w[treat == i] / mean(w[treat == i])
  }

  list(w = w, fit.obj = fit)
}

weightit2npcbps.multi <- weightit2npcbps

weightit2npcbps.cont <- function(covs, treat, s.weights, subset, missing, verbose, ...) {

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]

  missing <- .process_missing2(missing, covs)

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  for (i in seq_col(covs)) {
    covs[, i] <- .make_closer_to_1(covs[, i])
  }

  covs <- .apply_moments_int_quantile(covs,
                                      moments = ...get("moments"),
                                      int = ...get("int"))

  colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))
  covs <- covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]

  new.data <- data.frame(treat = treat, covs)

  corprior <- ...get("corprior", .01)

  tryCatch({verbosely({
    fit <- CBPS::npCBPS(formula(new.data),
                        data = new.data,
                        corprior = corprior,
                        print.level = 1)
  }, verbose = verbose)},
  error = function(e) {
    .err(sprintf("(from `CBPS::npCBPS()`): %s",
                 conditionMessage(e)),
         tidy = FALSE)
  })

  w <- fit$weights

  w <- w / mean(w)

  list(w = w, fit.obj = fit)
}
