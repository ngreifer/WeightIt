#' User-Defined Functions for Estimating Weights
#' @name method_user
#' @aliases method_user
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights using a user-defined function. The function must take in arguments that are passed to it by [weightit()] or [weightitMSM()] and return a vector of weights or a list containing the weights.
#'
#' To supply a user-defined function, the function object should be entered directly to `method`; for example, for a function `fun`, `method = fun`.
#'
#' ## Point Treatments
#'
#' The following arguments are automatically passed to the user-defined function, which should have named parameters corresponding to them:
#'   \itemize{
#'     \item{`treat`: a vector of treatment status for each unit. This comes directly from the left hand side of the formula passed to `weightit()` and so will have it's type (e.g., numeric, factor, etc.), which may need to be converted.}
#' \item{`covs`: a data frame of covariate values for each unit. This comes directly from the right hand side of the formula passed to `weightit()`. The covariates are processed so that all columns are numeric; all factor variables are split into dummies and all interactions are evaluated. All levels of factor variables are given dummies, so the matrix of the covariates is not full rank. Users can use [make_full_rank()], which accepts a numeric matrix or data frame and removes columns to make it full rank, if a full rank covariate matrix is desired.}
#' \item{`s.weights`: a numeric vector of sampling weights, one for each unit.}
#' \item{`ps`: a numeric vector of propensity scores.}
#' \item{`subset`: a logical vector the same length as `treat` that is `TRUE` for units to be included in the estimation and `FALSE` otherwise. This is used to subset the input objects when `exact` is used. `treat`, `covs`, `s.weights`, and `ps`, if supplied, will already have been subsetted by `subset`.}
#' \item{`estimand`: a character vector of length 1 containing the desired estimand. The characters will have been converted to uppercase. If "ATC" was supplied to estimand, `weightit()` sets `focal` to the control level (usually 0 or the lowest level of `treat`) and sets `estimand` to "ATT".}
#' \item{`focal`: a character vector of length 1 containing the focal level of the treatment when the estimand is the ATT (or the ATC as detailed above). `weightit()` ensures the value of focal is a level of `treat`.}
#' \item{`stabilize`: a logical vector of length 1. It is not processed by `weightit()` before it reaches the fitting function.}
#' \item{`moments`: a numeric vector of length 1. It is not processed by `weightit()` before it reaches the fitting function except that `as.integer()` is applied to it. This is used in other methods to determine whether polynomials of the entered covariates are to be used in the weight estimation.}
#' \item{`int`: a logical vector of length 1. It is not processed by `weightit()` before it reaches the fitting function. This is used in other methods to determine whether interactions of the entered covariates are to be used in the weight estimation.}
#' }
#' None of these parameters are required to be in the fitting function. These are simply those that are automatically available.
#'
#' In addition, any additional arguments supplied to `weightit()` will be passed on to the fitting function. `weightit()` ensures the arguments correspond to the parameters of the fitting function and throws an error if an incorrectly named argument is supplied and the fitting function doesn't include `\dots` as a parameter.
#'
#' The fitting function must output either a numeric vector of weights or a list (or list-like object) with an entry named wither "w" or "weights". If a list, the list can contain other named entries, but only entries named "w", "weights", "ps", and "fit.obj" will be processed. "ps" is a vector of propensity scores and "fit.obj" should be an object used in the fitting process that a user may want to examine and that is included in the `weightit` output object as "obj" when `include.obj = TRUE`. The "ps" and "fit.obj" components are optional, but "weights" or "w" is required.
#'
#' ## Longitudinal Treatments
#'
#' Longitudinal treatments can be handled either by running the fitting function for point treatments for each time point and multiplying the resulting weights together or by running a method that accommodates multiple time points and outputs a single set of weights. For the former, `weightitMSM()` can be used with the user-defined function just as it is with `weightit()`. The latter method is not yet accommodated by `weightitMSM()`, but will be someday, maybe.
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' @examples
#'
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #A user-defined version of method = "ps"
#' my.ps <- function(treat, covs, estimand, focal = NULL) {
#'   covs <- make_full_rank(covs)
#'   d <- data.frame(treat, covs)
#'   f <- formula(d)
#'   ps <- glm(f, data = d, family = "binomial")$fitted
#'   w <- get_w_from_ps(ps, treat = treat, estimand = estimand,
#'                      focal = focal)
#'
#'   list(w = w, ps = ps)
#' }
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = my.ps, estimand = "ATT"))
#' summary(W1)
#' bal.tab(W1)
#'
#' data("msmdata")
#' (W2 <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
#'                         A_2 ~ X1_1 + X2_1 +
#'                           A_1 + X1_0 + X2_0,
#'                         A_3 ~ X1_2 + X2_2 +
#'                           A_2 + X1_1 + X2_1 +
#'                           A_1 + X1_0 + X2_0),
#'                    data = msmdata,
#'                    method = my.ps))
#'
#' summary(W2)
#' bal.tab(W2)
#'
#' # Kernel balancing using the `kbal` package, available
#' # using `pak::pak_install("chadhazlett/KBAL")`.
#' # Only the ATT and ATC are available.
#'
#' \dontrun{
#'   kbal.fun <- function(treat, covs, estimand, focal, verbose, ...) {
#'     args <- list(...)
#'
#'     if (!estimand %in% c("ATT", "ATC"))
#'       stop("`estimand` must be \"ATT\" or \"ATC\".", call. = FALSE)
#'
#'     treat <- as.numeric(treat == focal)
#'
#'     args <- args[names(args) %in% names(formals(kbal::kbal))]
#'     args$allx <- covs
#'     args$treatment <- treat
#'     args$printprogress <- verbose
#'
#'     cat_cols <- apply(covs, 2, function(x) length(unique(x)) <= 2)
#'
#'     if (all(cat_cols)) {
#'       args$cat_data <- TRUE
#'       args$mixed_data <- FALSE
#'       args$scale_data <- FALSE
#'       args$linkernel <- FALSE
#'       args$drop_MC <- FALSE
#'     }
#'     else if (any(cat_cols)) {
#'       args$cat_data <- FALSE
#'       args$mixed_data <- TRUE
#'       args$cat_columns <- colnames(covs)[cat_cols]
#'       args$allx[,!cat_cols] <- scale(args$allx[,!cat_cols])
#'       args$cont_scale <- 1
#'     }
#'     else {
#'       args$cat_data <- FALSE
#'       args$mixed_data <- FALSE
#'     }
#'
#'     k.out <- do.call(kbal::kbal, args)
#'     w <- k.out$w
#'
#'     list(w = w, fit.obj = k.out)
#'   }
#'
#'   (Wk <- weightit(treat ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = kbal.fun, estimand = "ATT",
#'                   include.obj = TRUE))
#'   summary(Wk)
#'   bal.tab(Wk, stats = c("m", "ks"))
#' }
#'
NULL

weightit2user <- function(Fun, covs, treat, s.weights, subset, estimand, focal,
                          stabilize, subclass, missing, ps, moments, int, verbose, ...) {
  A <- list(...)

  if (is_not_null(covs)) {
    covs <- covs[subset, , drop = FALSE]
  }

  if (is_not_null(treat)) {
    treat <- treat[subset]
  }

  if (is_not_null(s.weights)) {
    s.weights <- s.weights[subset]
  }

  if (is_not_null(ps)) {
    ps <- ps[subset]
  }

  #Get a list of function args for the user-defined function Fun
  Fun_formal <- as.list(formals(Fun))
  if (has_dots <- ("..." %in% names(Fun_formal))) {
    Fun_formal[["..."]] <- NULL
  }

  fun_args <- Fun_formal
  for (i in names(fun_args)) {
    if (exists(i, inherits = FALSE)) {
      fun_args[i] <- list(get0(i, inherits = FALSE))
    }
    else if (i %in% names(A)) {
      fun_args[i] <- A[i]
      A[[i]] <- NULL
    }
    #else just use Fun default
  }

  if (has_dots) fun_args <- c(fun_args, A)

  obj <- do.call(Fun, fun_args)

  if (is.numeric(obj)) {
    obj <- list(w = obj)
  }
  else if (is.list(obj) && any(c("w", "weights") %in% names(obj))) {
    names(obj)[names(obj) == "weights"] <- "w"
  }
  else {
    .err('the output of the user-provided function must be a list with an entry named "w" or "weights" containing the estimated weights')
  }

  if (is_null(obj[["w"]]))
    .err("no weights were estimated")

  if (!is.numeric(obj[["w"]]) || is_not_null(dim(obj[["w"]])))
    .err("the \"w\" or \"weights\" entry of the output of the user-provided function must be a numeric vector of weights")

  if (all(is.na(obj[["w"]])))
    .err("all weights were generated as `NA`")

  if (length(obj[["w"]]) != length(treat))
    .err(sprintf("%s weights were estimated, but there are %s units",
                 length(obj[["w"]]), length(treat)))

  obj
}

weightitMSM2user <- function(Fun, covs.list, treat.list, s.weights, subset, stabilize,
                             missing, moments, int, verbose, ...) {
  A <- list(...)

  if (is_not_null(covs.list)) {
    for (i in seq_along(covs.list)) {
      covs.list[[i]] <- covs.list[[i]][subset, , drop = FALSE]
    }

    covs <- covs.list
  }

  if (is_not_null(treat.list)) {
    for (i in seq_along(treat.list)) {
      treat.list[[i]] <- treat.list[[i]][subset]
    }

    treat <- treat.list
  }

  if (is_not_null(s.weights)) {
    s.weights <- s.weights[subset]
  }

  #Get a list of function args for the user-defined function Fun
  Fun_formal <- as.list(formals(Fun))
  if (has_dots <- ("..." %in% names(Fun_formal))) {
    Fun_formal[["..."]] <- NULL
  }

  fun_args <- Fun_formal
  for (i in names(fun_args)) {
    if (exists(i, inherits = FALSE)) {
      fun_args[i] <- list(get0(i, inherits = FALSE))
    }
    else if (is_not_null(A[[i]])) {
      fun_args[[i]] <- A[[i]]
      A[[i]] <- NULL
    }
    #else just use Fun default
  }

  if (has_dots) fun_args <- c(fun_args, A)

  obj <- do.call(Fun, fun_args)

  if (is.numeric(obj)) {
    obj <- list(w = obj)
  }
  else if (is.list(obj) && any(c("w", "weights") %in% names(obj))) {
    names(obj)[names(obj) == "weights"] <- "w"
  }
  else {
    .err('the output of the user-provided function must be a list with an entry named "w" or "weights" containing the estimated weights')
  }

  if (is_null(obj[["w"]]))
    .err("no weights were estimated")

  if (!is.numeric(obj[["w"]]) || is_not_null(dim(obj[["w"]])))
    .err("the \"w\" or \"weights\" entry of the output of the user-provided function must be a numeric vector of weights")

  if (all(is.na(obj[["w"]])))
    .err("all weights were generated as `NA`")

  if (length(obj[["w"]]) != length(treat.list[[1]]))
    .err(sprintf("%s weights were estimated, but there are %s units",
                 length(obj[["w"]]), length(treat.list[[1]])))

  obj
}