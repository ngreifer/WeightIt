#' Optimization-Based Weighting
#' @name method_optweight
#' @aliases method_optweight
#' @usage NULL
#'
#' @description
#'This page explains the details of estimating optimization-based weights (also known as stable balancing weights) by setting `method = "optweight"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating weights by solving a quadratic programming problem subject to approximate or exact balance constraints. This method relies on \pkgfun{optweight}{optweight} from the \CRANpkg{optweight} package.
#'
#' Because `optweight()` offers finer control and uses the same syntax as `weightit()`, it is recommended that \pkgfun{optweight}{optweight} be used instead of `weightit()` with `method = "optweight"`.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the weights using \pkgfun{optweight}{optweight}. The following estimands are allowed: ATE, ATT, and ATC. The weights are taken from the output of the `optweight` fit object.
#'
#' ## Multi-Category Treatments
#'
#' For multi-category treatments, this method estimates the weights using \pkgfun{optweight}{optweight}. The following estimands are allowed: ATE and ATT. The weights are taken from the output of the `optweight` fit object.
#'
#' ## Continuous Treatments
#'
#' For binary treatments, this method estimates the weights using \pkgfun{optweight}{optweight}. The weights are taken from the output of the `optweight` fit object.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, `optweight()` estimates weights that simultaneously satisfy balance constraints at all time points, so only one model is fit to obtain the weights. Using `method = "optweight"` in `weightitMSM()` causes `is.MSM.method` to be set to `TRUE` by default. Setting it to `FALSE` will run one model for each time point and multiply the weights together, a method that is not recommended. NOTE: neither use of optimization-based weights with longitudinal treatments has been validated!
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios.
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
#'     A named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting `quantile = list(x1 = c(.25, .5. , .75))` ensures the 25th, 50th, and 75th percentiles of `x1` in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., `.5`) or an unnamed list of length 1 (e.g., `list(c(.25, .5, .75))`) to request the same quantile(s) for all continuous covariates, or a named vector (e.g., `c(x1 = .5, x2 = .75)` to request one quantile for each covariate. Only allowed with binary and multi-category treatments.
#'   }
#' }
#'
#' All arguments to `optweight()` can be passed through `weightit()` or `weightitMSM()`, with the following exception:
#'
#' * `targets` cannot be used and is ignored.
#'
#' All arguments take on the defaults of those in `optweight()`.
#'
#' @section Additional Outputs:
#' \describe{
#' \item{`info`}{
#' A list with one entry:
#'     \describe{
#'       \item{`duals`}{A data frame of dual variables for each balance constraint.}
#'     }
#' }
#' \item{`obj`}{When `include.obj = TRUE`, the output of the call to \pkgfun{optweight}{optweight}.}
#' }
#'
#' @details
#' Stable balancing weights are weights that solve a constrained optimization problem, where the constraints correspond to covariate balance and the loss function is the variance (or other norm) of the weights. These weights maximize the effective sample size of the weighted sample subject to user-supplied balance constraints. An advantage of this method over entropy balancing is the ability to allow approximate, rather than exact, balance through the `tols` argument, which can increase precision even for slight relaxations of the constraints.
#'
#' `plot()` can be used on the output of `weightit()` with `method = "optweight"` to display the dual variables; see Examples and [plot.weightit()] for more details.
#'
#' @note
#' The specification of `tols` differs between `weightit()` and `optweight()`. In `weightit()`, one tolerance value should be included per level of each factor variable, whereas in `optweight()`, all levels of a factor are given the same tolerance, and only one value needs to be supplied for a factor variable. Because of the potential for confusion and ambiguity, it is recommended to only supply one value for `tols` in `weightit()` that applies to all variables. For finer control, use `optweight()` directly.
#'
#' Seriously, just use \pkgfun{optweight}{optweight}. The syntax is almost identical and it's compatible with \pkg{cobalt}, too.
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' \pkgfun{optweight}{optweight} for the fitting function
#'
#' @references
#' ## Binary treatments
#'
#' Wang, Y., & Zubizarreta, J. R. (2020). Minimal dispersion approximately balancing weights: Asymptotic properties and practical considerations. *Biometrika*, 107(1), 93–105. \doi{10.1093/biomet/asz050}
#'
#' Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. *Journal of the American Statistical Association*, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}
#'
#' ## Multi-Category Treatments
#'
#' de los Angeles Resa, M., & Zubizarreta, J. R. (2020). Direct and stable weight adjustment in non-experimental studies with multivalued treatments: Analysis of the effect of an earthquake on post-traumatic stress. *Journal of the Royal Statistical Society: Series A (Statistics in Society)*, n/a(n/a). \doi{10.1111/rssa.12561}
#'
#' ## Continuous treatments
#'
#' Greifer, N. (2020). *Estimating Balancing Weights for Continuous Treatments Using Constrained Optimization*. \doi{10.17615/DYSS-B342}
#'
#' @examplesIf all(sapply(c("optweight", "osqp"), requireNamespace, quietly = TRUE))
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "optweight", estimand = "ATT",
#'                 tols = 0))
#' summary(W1)
#' cobalt::bal.tab(W1)
#' plot(W1)
#'
#' #Balancing covariates with respect to race (multi-category)
#' (W2 <- weightit(race ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "optweight", estimand = "ATE",
#'                 tols = .01))
#' summary(W2)
#' cobalt::bal.tab(W2)
#' plot(W2)
#'
#' #Balancing covariates with respect to re75 (continuous)
# (W3 <- weightit(re75 ~ age + educ + married +
#                   nodegree + re74, data = lalonde,
#                 method = "optweight", tols = .05))
# summary(W3)
# cobalt::bal.tab(W3)
# plot(W3)
NULL

weightit2optweight <- function(covs, treat, s.weights, subset, estimand, focal, missing,
                               moments, int, verbose, ...) {
  A <- list(...)

  rlang::check_installed("optweight")

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  covs <- cbind(.int_poly_f(covs, poly = moments, int = int, center = TRUE),
                .quantile_f(covs, qu = A[["quantile"]], s.weights = s.weights,
                            focal = focal, treat = treat))

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  new.data <- data.frame(treat, covs)
  new.formula <- formula(new.data)

  for (f in names(formals(optweight::optweight))) {
    if (is_null(A[[f]])) A[[f]] <- formals(optweight::optweight)[[f]]
  }
  A[names(A) %in% names(formals(weightit2optweight))] <- NULL

  if ("tols" %in% names(A)) A[["tols"]] <- optweight::check.tols(new.formula, new.data, A[["tols"]], stop = TRUE)
  if ("targets" %in% names(A)) {
    .wrn("`targets` cannot be used through WeightIt and will be ignored")
    A[["targets"]] <- NULL
  }

  A[["formula"]] <- new.formula
  A[["data"]] <- new.data
  A[["estimand"]] <- estimand
  A[["s.weights"]] <- s.weights
  A[["focal"]] <- focal
  A[["verbose"]] <- TRUE

  verbosely({
    out <- do.call(optweight::optweight, A, quote = TRUE)
  }, verbose = verbose)

  list(w = out[["weights"]], info = list(duals = out$duals), fit.obj = out)
}

weightit2optweight.multi <- weightit2optweight

weightit2optweight.cont <- function(covs, treat, s.weights, subset, missing, moments, int, verbose, ...) {
  A <- list(...)
  rlang::check_installed("optweight")

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  missing <- .process_missing2(missing, covs)

  covs <- .int_poly_f(covs, poly = moments, int = int)

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs)
  }

  new.data <- data.frame(treat, covs)
  new.formula <- formula(new.data)

  for (f in names(formals(optweight::optweight))) {
    if (is_null(A[[f]])) A[[f]] <- formals(optweight::optweight)[[f]]
  }
  A[names(A) %in% names(formals(weightit2optweight.cont))] <- NULL

  if ("tols" %in% names(A)) {
    A[["tols"]] <- optweight::check.tols(new.formula, new.data, A[["tols"]], stop = TRUE)
  }

  if ("targets" %in% names(A)) {
    .wrn("`targets` cannot be used through WeightIt and will be ignored")
    A[["targets"]] <- NULL
  }

  A[["formula"]] <- new.formula
  A[["data"]] <- new.data
  A[["s.weights"]] <- s.weights
  A[["verbose"]] <- TRUE

  verbosely({
    out <- do.call(optweight::optweight, A, quote = TRUE)
  }, verbose = verbose)

  list(w = out[["weights"]], info = list(duals = out$duals), fit.obj = out)
}

weightitMSM2optweight <- function(covs.list, treat.list, s.weights, subset, missing, moments, int, verbose, ...) {
  A <- list(...)
  rlang::check_installed("optweight")

  s.weights <- s.weights[subset]
  treat.types <- character(length(treat.list))

  for (i in seq_along(treat.list)) {
    treat.list[[i]] <- treat.list[[i]][subset]

    if (!has_treat_type(treat.list[[i]])) {
      treat.list[[i]] <- assign_treat_type(treat.list[[i]])
    }
    treat.types[i] <- get_treat_type(treat.list[[i]])

    if (get_treat_type(treat.list[[i]]) != "continuous") {
      treat.list[[i]] <- factor(treat.list[[i]])
    }

    covs.list[[i]] <- covs.list[[i]][subset, , drop = FALSE]

    if (.process_missing2(missing, covs.list[[i]]) == "ind") {
      covs.list[[i]] <- add_missing_indicators(covs.list[[i]])
    }

    covs.list[[i]] <- cbind(covs.list[[i]],
                            .int_poly_f(covs.list[[i]], poly = moments, int = int))

    if (treat.types[i] %in% c("binary", "multi-category")) {
      covs.list[[i]] <- cbind(.int_poly_f(covs.list[[i]], poly = moments, int = int, center = TRUE),
                              .quantile_f(covs.list[[i]], qu = A[["quantile"]], s.weights = s.weights,
                                          treat = treat.list[[i]]))
    }
    else {
      covs.list[[i]] <- cbind(covs.list[[i]], .int_poly_f(covs.list[[i]], poly = moments,
                                                          int = int, center = TRUE))
    }

    for (j in seq_col(covs.list[[i]])) {
      covs.list[[i]][,j] <- .make_closer_to_1(covs.list[[i]][,j])
    }
  }

  baseline.data <- data.frame(treat.list[[1]], covs.list[[1]])
  baseline.formula <- formula(baseline.data)
  if ("tols" %in% names(A)) {
    A[["tols"]] <- optweight::check.tols(baseline.formula, baseline.data, A[["tols"]], stop = TRUE)
  }

  if ("targets" %in% names(A)) {
    .wrn("`targets` cannot be used through WeightIt and will be ignored")
    A[["targets"]] <- NULL
  }

  verbosely({
    out <- do.call(optweight::optweight.fit,
                   c(list(treat = treat.list,
                          covs = covs.list,
                          s.weights = s.weights,
                          verbose = TRUE),
                     A),
                   quote = TRUE)
  }, verbose = verbose)

  list(w = out$w, fit.obj = out)
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

  title <- "Dual Variables for Constraints"
  # }
  d$cov <- factor(d$cov, levels = rev(unique(d$cov)))
  d$constraint <- factor(d$constraint, levels = unique(d$constraint, nmax = 2),
                         labels = paste("Constraint:", unique(d$constraint, nmax = 2)))

  p <- ggplot(d, aes(y = .data$cov, x = .data$dual)) +
    geom_col() +
    geom_vline(xintercept = 0) +
    labs(x = "Absolute Dual Variable", y = "Covariate", title = title) +
    scale_x_continuous(expand = expansion(c(0, 0.05))) +
    facet_grid(rows = vars(.data$constraint),
               if (use.by) vars(.data$by) else NULL) +
    theme_bw()
  p
}
