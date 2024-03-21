#' Propensity Score Weighting Using Generalized Boosted Models
#' @name method_gbm
#' @aliases method_gbm
#' @usage NULL
#'
#' @description
#' This page explains the details of estimating weights from generalized boosted model-based propensity scores by setting `method = "gbm"` in the call to [weightit()] or [weightitMSM()]. This method can be used with binary, multi-category, and continuous treatments.
#'
#' In general, this method relies on estimating propensity scores using generalized boosted modeling and then converting those propensity scores into weights using a formula that depends on the desired estimand. The algorithm involves using a balance-based or prediction-based criterion to optimize in choosing the value of tuning parameters (the number of trees and possibly others). The method relies on the \CRANpkg{gbm} package.
#'
#' This method mimics the functionality of functions in the \pkg{twang} package, but has improved performance and more flexible options. See Details section for more details.
#'
#' ## Binary Treatments
#'
#' For binary treatments, this method estimates the propensity scores using \pkgfun{gbm}{gbm.fit} and then selects the optimal tuning parameter values using the method specified in the `criterion` argument. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights are computed from the estimated propensity scores using [get_w_from_ps()], which implements the standard formulas. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See [get_w_from_ps()] for details.
#'
#' ## Multi-Category Treatments
#'
#' For binary treatments, this method estimates the propensity scores using \pkgfun{gbm}{gbm.fit} and then selects the optimal tuning parameter values using the method specified in the `criterion` argument. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights are computed from the estimated propensity scores using [get_w_from_ps()], which implements the standard formulas. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See [get_w_from_ps()] for details.
#'
#' ## Continuous Treatments
#'
#' For continuous treatments, this method estimates the generalized propensity score using \pkgfun{gbm}{gbm.fit} and then selects the optimal tuning parameter values using the method specified in the `criterion` argument.
#'
#' ## Longitudinal Treatments
#'
#' For longitudinal treatments, the weights are the product of the weights estimated at each time point.
#'
#' ## Sampling Weights
#'
#' Sampling weights are supported through `s.weights` in all scenarios.
#'
#' ## Missing Data
#'
#' In the presence of missing data, the following value(s) for `missing` are allowed:
#'     \describe{
#'       \item{`"ind"` (default)}{First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is `NA` and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with the covariate medians (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting `weightit` object will be the original covariates with the `NA`s.}
#'       \item{`"surr"`}{Surrogate splitting is used to process `NA`s. No missingness indicators are created. Nodes are split using only the non-missing values of each variable. To generate predicted values for each unit, a non-missing variable that operates similarly to the variable with missingness is used as a surrogate. Missing values are ignored when calculating balance statistics to choose the optimal tree.}
#'     }
#'
#' ## M-estimation
#'
#' M-estimation is not supported.
#'
#' @section Additional Arguments:
#' The following additional arguments can be specified:
#'   \describe{
#'     \item{`criterion`}{A string describing the balance criterion used to select the best weights. See \pkgfun{cobalt}{bal.compute} for allowable options for each treatment type. In addition, to optimize the cross-validation error instead of balance, `criterion` can be set as `"cv{#}`", where `{#}` is replaced by a number representing the number of cross-validation folds used (e.g., `"cv5"` for 5-fold cross-validation). For binary and multi-category treatments, the default is `"smd.mean"`, which minimizes the average absolute standard mean difference among the covariates between treatment groups. For continuous treatments, the default is `"p.mean"`, which minimizes the average absolute Pearson correlation between the treatment and covariates.
#'     }
#'       \item{`trim.at`}{A number supplied to `at` in [trim()] which trims the weights from all the trees before choosing the best tree. This can be valuable when some weights are extreme, which occurs especially with continuous treatments. The default is 0 (i.e., no trimming).
#'       }
#'       \item{`distribution`}{A string with the distribution used in the loss function of the boosted model. This is supplied to the `distribution` argument in \pkgfun{gbm}{gbm.fit}. For binary treatments, `"bernoulli"` and `"adaboost"` are available, with `"bernoulli"` the default. For multi-category treatments, only `"multinomial"` is allowed. For continuous treatments `"gaussian"`, `"laplace"`, and `"tdist"` are available, with `"gaussian"` the default. This argument is tunable.
#'       }
#'       \item{`n.trees`}{The maximum number of trees used. This is passed onto the `n.trees` argument in `gbm.fit()`. The default is 10000 for binary and multi-category treatments and 20000 for continuous treatments.
#'       }
#'       \item{`start.tree`}{The tree at which to start balance checking. If you know the best balance isn't in the first 100 trees, for example, you can set `start.tree = 101` so that balance statistics are not computed on the first 100 trees. This can save some time since balance checking takes up the bulk of the run time for some balance-based stopping methods, and is especially useful when running the same model adding more and more trees. The default is 1, i.e., to start from the very first tree in assessing balance.
#'       }
#'       \item{`interaction.depth`}{The depth of the trees. This is passed onto the `interaction.depth` argument in `gbm.fit()`. Higher values indicate better ability to capture nonlinear and nonadditive relationships. The default is 3 for binary and multi-category treatments and 4 for continuous treatments. This argument is tunable.
#'       }
#'       \item{`shrinkage`}{The shrinkage parameter applied to the trees. This is passed onto the `shrinkage` argument in `gbm.fit()`. The default is .01 for binary and multi-category treatments and .0005 for continuous treatments. The lower this value is, the more trees one may have to include to reach the optimum. This argument is tunable.
#'       }
#'       \item{`bag.fraction`}{The fraction of the units randomly selected to propose the next tree in the expansion. This is passed onto the `bag.fraction` argument in `gbm.fit()`. The default is 1, but smaller values should be tried. For values less then 1, subsequent runs with the same parameters will yield different results due to random sampling; be sure to seed the seed using [set.seed()] to ensure replicability of results.
#'        }
#' }
#'
#' All other arguments take on the defaults of those in \pkgfun{gbm}{gbm.fit}, and some are not used at all.
#'
#' The `w` argument in `gbm.fit()` is ignored because sampling weights are passed using `s.weights`.
#'
#' For continuous treatments only, the following arguments may be supplied:
#' \describe{
#'       \item{`density`}{A function corresponding to the conditional density of the treatment. The standardized residuals of the treatment model will be fed through this function to produce the numerator and denominator of the generalized propensity score weights. If blank, [dnorm()] is used as recommended by Robins et al. (2000). This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if `density = "dt_2"` is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme outcome values are observed (Naimi et al., 2014). Ignored if `use.kernel = TRUE` (described below).
#'       }
#'       \item{`use.kernel`}{If `TRUE`, uses kernel density estimation through the [density()] function to estimate the numerator and denominator densities for the weights. If `FALSE` (the default), the argument to the `density` parameter is used instead.
#'       }
#'       \item{`bw`, `adjust`, `kernel`, `n`}{If `use.kernel = TRUE`, the arguments to [density()]. The defaults are the same as those in `density` except that `n` is 10 times the number of units in the sample.
#'       }
#'       \item{`plot`}{If `use.kernel = TRUE` with continuous treatments, whether to plot the estimated density.
#'       }
#' }
#'
#' For tunable arguments, multiple entries may be supplied, and `weightit()` will choose the best value by optimizing the criterion specified in `criterion`. See below for additional outputs that are included when arguments are supplied to be tuned. See Examples for an example of tuning.
#'
#' @section Additional Outputs:
#' \describe{
#' \item{`info`}{
#'   A list with the following entries:
#'     \describe{
#'       \item{`best.tree`}{
#'         The number of trees at the optimum. If this is close to `n.trees`, `weightit()` should be rerun with a larger value for `n.trees`, and `start.tree` can be set to just below `best.tree`. When other parameters are tuned, this is the best tree value in the best combination of tuned parameters. See example.
#'       }
#'       \item{`tree.val`}{
#'         A data frame with two columns: the first is the number of trees and the second is the value of the criterion corresponding to that tree. Running [plot()] on this object will plot the criterion by the number of trees and is a good way to see patterns in the relationship between them and to determine if more trees are needed. When other parameters are tuned, these are the number of trees and the criterion values in the best combination of tuned parameters. See example.
#'       }
#'     }
#'   If any arguments are to be tuned (i.e., they have been supplied more than one value), the following two additional components are included in `info`:
#'     \describe{
#'       \item{`tune`}{
#'         A data frame with a column for each argument being tuned, the best value of the balance criterion for the given combination of parameters, and the number of trees at which the best value was reached.
#'       }
#'       \item{`best.tune`}{
#'         A one-row data frame containing the values of the arguments being tuned that were ultimately selected to estimate the returned weights.
#'       }
#'     }
#' }
#' \item{`obj`}{
#'   When `include.obj = TRUE`, the `gbm` fit used to generate the predicted values.
#' }
#' }
#'
#' @details
#' Generalized boosted modeling (GBM, also known as gradient boosting machines) is a machine learning method that generates predicted values from a flexible regression of the treatment on the covariates, which are treated as propensity scores and used to compute weights. It does this by building a series of regression trees, each fit to the residuals of the last, minimizing a loss function that depends on the distribution chosen. The optimal number of trees is a tuning parameter that must be chosen; McCaffrey et al. (2004) were innovative in using covariate balance to select this value rather than traditional machine learning performance metrics such as cross-validation accuracy. GBM is particularly effective for fitting nonlinear treatment models characterized by curves and interactions, but performs worse for simpler treatment models. It is unclear which balance measure should be used to select the number of trees, though research has indicated that balance measures tend to perform better than cross-validation accuracy for estimating effective propensity score weights.
#'
#' \pkg{WeightIt} offers almost identical functionality to \pkg{twang}, the first package to implement this method. Compared to the current version of \pkg{twang}, \pkg{WeightIt} offers more options for the measure of balance used to select the number of trees, improved performance, tuning of hyperparameters, more estimands, and support for continuous treatments. \pkg{WeightIt} computes weights for multi-category treatments differently from how \pkg{twang} does; rather than fitting a separate binary GBM for each pair of treatments, \pkg{WeightIt} fits a single multi-class GBM model and uses balance measures appropriate for multi-category treatments.
#'
#' @note
#' The `criterion` argument used to be called `stop.method`, which is its name in \pkg{twang}. `stop.method` still works for backward compatibility. Additionally, the criteria formerly named as `es.mean`, `es.max`, and `es.rms` have been renamed to `smd.mean`, `smd.max`, and `smd.rms`. The former are used in \pkg{twang} and will still work with `weightit()` for backward compatibility.
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' \pkgfun{gbm}{gbm.fit} for the fitting function.
#'
#' @references
#' ## Binary treatments
#'
#' McCaffrey, D. F., Ridgeway, G., & Morral, A. R. (2004). Propensity Score Estimation With Boosted Regression for Evaluating Causal Effects in Observational Studies. Psychological Methods, 9(4), 403–425. \doi{10.1037/1082-989X.9.4.403}
#'
#' ## Multi-Category Treatments
#'
#' McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. Statistics in Medicine, 32(19), 3388–3414. \doi{10.1002/sim.5753}
#'
#'
#' ## Continuous treatments
#'
#' Zhu, Y., Coffman, D. L., & Ghosh, D. (2015). A Boosting Algorithm for Estimating Generalized Propensity Scores with Continuous Treatments. Journal of Causal Inference, 3(1). \doi{10.1515/jci-2014-0022}
#'
#' @examplesIf requireNamespace("gbm", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "gbm", estimand = "ATE",
#'                 criterion = "smd.max"))
#' summary(W1)
#' bal.tab(W1)
#'
#' \dontrun{
#'   #Balancing covariates with respect to race (multi-category)
#'   (W2 <- weightit(race ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "gbm", estimand = "ATT",
#'                   focal = "hispan", criterion = "ks.mean"))
#'   summary(W2)
#'   bal.tab(W2, stats = c("m", "ks"))
#'
#'   #Balancing covariates with respect to re75 (continuous)
#'   (W3 <- weightit(re75 ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "gbm", use.kernel = TRUE,
#'                   criterion = "p.rms", trim.at = .97))
#'   summary(W3)
#'   bal.tab(W3)
#'
#'   #Using a t(3) density and illustrating the search for
#'   #more trees.
#'   W4a <- weightit(re75 ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "gbm", density = "dt_3",
#'                   criterion = "p.max",
#'                   n.trees = 10000)
#'
#'   W4a$info$best.tree #10000; optimum hasn't been found
#'   plot(W4a$info$tree.val, type = "l") #decreasing at right edge
#'
#'   W4b <- weightit(re75 ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "gbm", density = "dt_3",
#'                   criterion = "p.max",
#'                   start.tree = 10000,
#'                   n.trees = 20000)
#'
#'   W4b$info$best.tree #13417; optimum has been found
#'   plot(W4b$info$tree.val, type = "l") #increasing at right edge
#'
#'   bal.tab(W4b)
#'
#'   #Tuning hyperparameters
#'   (W5 <- weightit(treat ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "gbm", estimand = "ATT",
#'                   criterion = "ks.max",
#'                   interaction.depth = 2:4,
#'                   distribution = c("bernoulli", "adaboost")))
#'
#'   W5$info$tune
#'
#'   W5$info$best.tune #Best values of tuned parameters
#'
#'   bal.tab(W5, stats = c("m", "ks"))
#' }
NULL

weightit2gbm <- function(covs, treat, s.weights, estimand, focal, subset,
                         stabilize, subclass, missing, verbose, ...) {

  rlang::check_installed("gbm")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (!has_treat_type(treat)) treat <- assign_treat_type(treat)
  treat.type <- get_treat_type(treat)

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs, replace_with = NA)
  }

  if (is_null(A[["criterion"]])) {
    A[["criterion"]] <- A[["stop.method"]]
  }

  if (is_null(A[["criterion"]])) {
    .wrn("no `criterion` was provided. Using \"smd.mean\"")
    A[["criterion"]] <- "smd.mean"
  }
  else if (length(A[["criterion"]]) > 1) {
    .wrn("only one `criterion` is allowed at a time. Using just the first `criterion`")
    A[["criterion"]] <- A[["criterion"]][1]
  }

  available.criteria <- cobalt::available.stats(treat.type)

  if (is.character(A[["criterion"]]) &&
      startsWith(A[["criterion"]], "es.")) {
    subbed.crit <- sub("es.", "smd.", A[["criterion"]], fixed = TRUE)
    subbed.match <- charmatch(subbed.crit, available.criteria)
    if (!anyNA(subbed.match) && subbed.match != 0L) {
      A[["criterion"]] <- subbed.crit
    }
  }

  cv <- 0

  s.m.matches <- charmatch(A[["criterion"]], available.criteria)
  if (anyNA(s.m.matches) || s.m.matches == 0L) {
    if (startsWith(A[["criterion"]], "cv") &&
        can_str2num(numcv <- substr(A[["criterion"]], 3, nchar(A[["criterion"]])))) {
      cv <- round(str2num(numcv))
      if (cv < 2) .err("at least 2 CV-folds must be specified in `criterion`")
    }
    else {
      .err(sprintf("`criterion` must be one of %s.",
                   word_list(c(available.criteria, "cv{#}"), "or", quotes = TRUE)))
    }
  }
  else criterion <- available.criteria[s.m.matches]

  tunable <- c("interaction.depth", "shrinkage", "distribution")

  trim.at <- if_null_then(A[["trim.at"]], 0)

  for (f in names(formals(gbm::gbm.fit))) {
    if (is_null(A[[f]])) {
      if (f %in% c("x", "y", "misc", "w", "verbose", "var.names",
                   "response.name", "group", "distribution")) A[f] <- list(NULL)
      else A[f] <- list(switch(f, n.trees = 1e4,
                               interaction.depth = 3,
                               shrinkage = .01,
                               bag.fraction = 1,
                               keep.data = FALSE,
                               formals(gbm::gbm.fit)[[f]]))
    }
  }

  chk::chk_count(A[["n.trees"]], "`n.trees`")
  chk::chk_gt(A[["n.trees"]], 1, "`n.trees`")

  if (treat.type == "binary")  {
    available.distributions <- c("bernoulli", "adaboost")
    t.lev <- get_treated_level(treat)
    treat <- binarize(treat, one = focal)
  }
  else {
    available.distributions <- "multinomial"
    treat <- factor(treat)
  }

  if (cv == 0) {
    start.tree <- if_null_then(A[["start.tree"]], 1)
    if (is_null(A[["n.grid"]])) n.grid <- round(1 + sqrt(2 * (A[["n.trees"]] - start.tree + 1)))
    else if (!is_(A[["n.grid"]], "numeric") || length(A[["n.grid"]]) > 1 ||
             !between(A[["n.grid"]], c(2, A[["n.trees"]]))) {
      .err("`n.grid` must be a numeric value between 2 and `n.trees`")
    }
    else n.grid <- round(A[["n.grid"]])

    init <- cobalt::bal.init(
      if (!anyNA(covs)) covs
      else if (missing == "surr") add_missing_indicators(covs)
      else replace_na_with(covs),
      treat = treat, stat = criterion,
      estimand = estimand, s.weights = s.weights,
      focal = focal, ...)
  }

  A[["x"]] <- covs
  A[["y"]] <- treat
  A[["distribution"]] <- if (is_null(distribution <- A[["distribution"]])) {
    available.distributions[1]} else {
      match_arg(distribution, available.distributions, several.ok = TRUE)}
  A[["w"]] <- s.weights
  A[["verbose"]] <- FALSE

  tune <- do.call("expand.grid", c(A[names(A) %in% tunable],
                                   list(stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)))

  current.best.loss <- Inf

  for (i in seq_row(tune)) {

    A[["distribution"]] <- list(name = tune[["distribution"]][i])
    tune_args <- as.list(tune[i, setdiff(tunable, "distribution")])

    gbm.call <- as.call(c(list(quote(gbm::gbm.fit)),
                          A[names(A) %in% setdiff(names(formals(gbm::gbm.fit)), names(tune_args))],
                          tune_args))
    verbosely({
      fit <- eval(gbm.call)
    }, verbose = verbose)

    if (cv == 0) {

      n.trees <- fit[["n.trees"]]
      iters <- 1:n.trees
      iters.grid <- round(seq(start.tree, n.trees, length.out = n.grid))

      if (is_null(iters.grid) || anyNA(iters.grid) || any(iters.grid > n.trees)) {
        .err("a problem has occurred")
      }

      ps <- gbm::predict.gbm(fit, n.trees = iters.grid, type = "response", newdata = covs)
      w <- .get_w_from_ps_internal_array(ps, treat = treat, estimand = estimand,
                                         focal = focal, stabilize = stabilize, subclass = subclass)
      if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

      iter.grid.balance <- apply(w, 2, cobalt::bal.compute, x = init)

      if (n.grid == n.trees) {
        best.tree.index <- which.min(iter.grid.balance)
        best.loss <- iter.grid.balance[best.tree.index]
        best.tree <- iters.grid[best.tree.index]
        tree.val <- setNames(data.frame(iters.grid,
                                        iter.grid.balance),
                             c("tree", criterion))
      }
      else {
        it <- which.min(iter.grid.balance) + c(-1, 1)
        it[1] <- iters.grid[max(1, it[1])]
        it[2] <- iters.grid[min(length(iters.grid), it[2])]
        iters.to.check <- iters[between(iters, iters[it])]

        if (is_null(iters.to.check) || anyNA(iters.to.check) || any(iters.to.check > n.trees)) {
          .err("a problem has occurred")
        }

        ps <- gbm::predict.gbm(fit, n.trees = iters.to.check, type = "response", newdata = covs)
        w <- .get_w_from_ps_internal_array(ps, treat = treat, estimand = estimand,
                                           focal = focal, stabilize = stabilize, subclass = subclass)
        if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

        iter.grid.balance.fine <- apply(w, 2, cobalt::bal.compute, x = init)

        best.tree.index <- which.min(iter.grid.balance.fine)
        best.loss <- iter.grid.balance.fine[best.tree.index]
        best.tree <- iters.to.check[best.tree.index]
        tree.val <- setNames(data.frame(c(iters.grid, iters.to.check),
                                        c(iter.grid.balance, iter.grid.balance.fine)),
                             c("tree", criterion))
      }

      tree.val <- unique(tree.val[order(tree.val$tree),])
      w <- w[,best.tree.index]
      ps <- if (treat.type == "binary") ps[,best.tree.index] else NULL

      tune[[paste.("best", criterion)]][i] <- best.loss
      tune[["best.tree"]][i] <- best.tree

      if (best.loss < current.best.loss) {
        best.fit <- fit
        best.w <- w
        best.ps <- ps
        current.best.loss <- best.loss
        best.tune.index <- i

        info <- list(best.tree = best.tree,
                     tree.val = tree.val)
      }
    }
    else {
      A["data"] <- list(data.frame(treat, covs))
      A[["cv.folds"]] <- cv
      A["n.cores"] <- list(A[["n.cores"]])
      A["var.names"] <- list(A[["var.names"]])
      A["offset"] <- list(NULL)
      A[["nTrain"]] <- floor(nrow(covs))
      A[["class.stratify.cv"]] <- FALSE
      A[["y"]] <- treat
      A[["x"]] <- covs
      A[["distribution"]] <- list(name = tune[["distribution"]][i])
      A[["w"]] <- s.weights

      tune_args <- as.list(tune[i, setdiff(tunable, "distribution")])

      gbmCrossVal.call <- as.call(c(list(quote(gbm::gbmCrossVal)),
                                    A[names(A) %in% setdiff(names(formals(gbm::gbmCrossVal)), names(tune_args))],
                                    tune_args))

      verbosely({
        cv.results <- eval(gbmCrossVal.call)
      }, verbose = verbose)

      best.tree.index <- which.min(cv.results$error)
      best.loss <- cv.results$error[best.tree.index]
      best.tree <- best.tree.index

      tune[[paste.("best", names(fit$name))]][i] <- best.loss
      tune[["best.tree"]][i] <- best.tree

      if (best.loss < current.best.loss) {
        best.fit <- fit
        best.ps <- gbm::predict.gbm(fit, n.trees = best.tree, type = "response", newdata = covs)
        best.w <- drop(.get_w_from_ps_internal_array(best.ps, treat = treat, estimand = estimand,
                                                     focal = focal, stabilize = stabilize, subclass = subclass))
        # if (trim.at != 0) best.w <- suppressMessages(trim(best.w, at = trim.at, treat = treat))
        current.best.loss <- best.loss
        best.tune.index <- i

        tree.val <- data.frame(tree = seq_along(cv.results$error),
                               error = cv.results$error)

        info <- list(best.tree = best.tree,
                     tree.val = tree.val)

        if (treat.type == "multinomial") best.ps <- NULL
      }
    }

    if (treat.type == "multinomial") ps <- NULL
  }

  tune[tunable[lengths(A[tunable]) == 1]] <- NULL

  if (ncol(tune) > 2) {
    info[["tune"]] <- tune
    info[["best.tune"]] <- tune[best.tune.index,]
  }

  if (is_not_null(best.ps) && is_not_null(focal) && focal != t.lev) {
    best.ps <- 1 - best.ps
  }

  list(w = best.w, ps = best.ps, info = info, fit.obj = best.fit)
}

weightit2gbm.multi <- weightit2gbm

weightit2gbm.cont <- function(covs, treat, s.weights, estimand, focal, subset,
                              stabilize, subclass, missing, verbose, ...) {

  rlang::check_installed("gbm")

  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs, replace_with = NA)
  }

  if (is_null(A[["criterion"]])) {
    A[["criterion"]] <- A[["stop.method"]]
  }
  if (is_null(A[["criterion"]])) {
    .wrn("no `criterion` was provided. Using \"p.mean\"")
    A[["criterion"]] <- "p.mean"
  }
  else if (length(A[["criterion"]]) > 1) {
    .wrn("only one `criterion` is allowed at a time. Using just the first `criterion`")
    A[["criterion"]] <- A[["criterion"]][1]
  }

  available.criteria <- cobalt::available.stats("continuous")

  cv <- 0

  s.m.matches <- charmatch(A[["criterion"]], available.criteria)
  if (anyNA(s.m.matches) || s.m.matches == 0L) {
    if (startsWith(A[["criterion"]], "cv") &&
        can_str2num(numcv <- substr(A[["criterion"]], 3, nchar(A[["criterion"]])))) {
      cv <- round(str2num(numcv))
      if (cv < 2) .err("at least 2 CV-folds must be specified in `criterion`")
    }
    else .err(sprintf("`criterion` must be one of %s",
                      word_list(c(available.criteria, "cv{#}"), "or", quotes = TRUE)))
  }
  else criterion <- available.criteria[s.m.matches]

  tunable <- c("interaction.depth", "shrinkage", "distribution")

  trim.at <- if_null_then(A[["trim.at"]], 0)

  for (f in names(formals(gbm::gbm.fit))) {
    if (is_null(A[[f]])) {
      if (f %in% c("x", "y", "misc", "w", "verbose", "var.names",
                   "response.name", "group", "distribution")) A[f] <- list(NULL)
      else A[f] <- list(switch(f, n.trees = 2e4,
                               interaction.depth = 4,
                               shrinkage = 0.0005,
                               bag.fraction = 1,
                               formals(gbm::gbm.fit)[[f]]))
    }
  }

  chk::chk_count(A[["n.trees"]], "`n.trees`")
  chk::chk_gt(A[["n.trees"]], 1, "`n.trees`")

  available.distributions <- c("gaussian", "laplace", "tdist", "poisson")

  if (cv == 0) {
    start.tree <- if_null_then(A[["start.tree"]], 1)
    if (is_null(A[["n.grid"]])) n.grid <- round(1 + sqrt(2 * (A[["n.trees"]] - start.tree + 1)))
    else if (!is_(A[["n.grid"]], "numeric") || length(A[["n.grid"]]) > 1 ||
             !between(A[["n.grid"]], c(2, A[["n.trees"]]))) {
      .err("`n.grid` must be a numeric value between 2 and `n.trees`")
    }
    else n.grid <- round(A[["n.grid"]])

    init <- cobalt::bal.init(
      if (!anyNA(covs)) covs
      else if (missing == "surr") add_missing_indicators(covs)
      else replace_na_with(covs),
      treat = treat, stat = criterion,
      s.weights = s.weights, ...)
  }

  A[["x"]] <- covs
  A[["y"]] <- treat
  A[["distribution"]] <- {
    if (is_null(distribution <- A[["distribution"]])) {
      available.distributions[1]
    }
    else {
      match_arg(distribution, available.distributions, several.ok = TRUE)
    }
  }
  A[["w"]] <- s.weights
  A[["verbose"]] <- FALSE

  tune <- do.call("expand.grid", c(A[names(A) %in% tunable],
                                   list(stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)))

  #Process density params
  densfun <- get_dens_fun(use.kernel = isTRUE(A[["use.kernel"]]), bw = A[["bw"]],
                          adjust = A[["adjust"]], kernel = A[["kernel"]],
                          n = A[["n"]], treat = treat, density = A[["density"]])

  #Stabilization - get dens.num
  dens.num <- densfun(treat - mean(treat), s.weights)

  current.best.loss <- Inf

  for (i in seq_row(tune)) {
    gbm.call <- as.call(c(list(quote(gbm::gbm.fit)),
                          A[names(A) %in% setdiff(names(formals(gbm::gbm.fit)), tunable)],
                          tune[i, tunable[tunable %in% names(formals(gbm::gbm.fit))]]))
    verbosely({
      fit <- eval(gbm.call)
    }, verbose = verbose)

    if (cv == 0) {

      n.trees <- fit[["n.trees"]]
      iters <- 1:n.trees
      iters.grid <- round(seq(start.tree, n.trees, length.out = n.grid))

      if (is_null(iters.grid) || anyNA(iters.grid) || any(iters.grid > n.trees)) {
        .err("a problem has occurred")
      }

      gps <- gbm::predict.gbm(fit, n.trees = iters.grid, newdata = covs)

      w <- apply(gps, 2, function(gp.score) {
        dens.num/densfun(treat - gp.score, s.weights)
      })

      if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

      iter.grid.balance <- apply(w, 2, cobalt::bal.compute, x = init)

      if (n.grid == n.trees) {
        best.tree.index <- which.min(iter.grid.balance)
        best.loss <- iter.grid.balance[best.tree.index]
        best.tree <- iters.grid[best.tree.index]
        tree.val <- setNames(data.frame(iters.grid,
                                        iter.grid.balance),
                             c("tree", criterion))
      }
      else {
        it <- which.min(iter.grid.balance) + c(-1, 1)
        it[1] <- iters.grid[max(1, it[1])]
        it[2] <- iters.grid[min(length(iters.grid), it[2])]
        iters.to.check <- iters[between(iters, iters[it])]

        if (is_null(iters.to.check) || anyNA(iters.to.check) || any(iters.to.check > n.trees)) {
          .err("a problem has occurred")
        }

        gps <- gbm::predict.gbm(fit, n.trees = iters.to.check, newdata = covs)
        w <- apply(gps, 2, function(gp.score) {
          dens.num/densfun(treat - gp.score, s.weights)
        })

        if (trim.at != 0) w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))

        iter.grid.balance.fine <- apply(w, 2, cobalt::bal.compute, x = init)

        best.tree.index <- which.min(iter.grid.balance.fine)
        best.loss <- iter.grid.balance.fine[best.tree.index]
        best.tree <- iters.to.check[best.tree.index]
        tree.val <- setNames(data.frame(c(iters.grid, iters.to.check),
                                        c(iter.grid.balance, iter.grid.balance.fine)),
                             c("tree", criterion))
      }

      tree.val <- unique(tree.val[order(tree.val$tree),])
      w <- w[,best.tree.index]
      gps <- gps[,as.character(best.tree)]

      tune[[paste.("best", criterion)]][i] <- best.loss
      tune[["best.tree"]][i] <- best.tree

      if (best.loss < current.best.loss) {
        best.fit <- fit
        best.w <- w
        best.gps <- gps
        current.best.loss <- best.loss
        best.tune.index <- i

        info <- list(best.tree = best.tree,
                     tree.val = tree.val)
      }
    }
    else {
      A["data"] <- list(data.frame(treat, covs))
      A[["cv.folds"]] <- cv
      A["n.cores"] <- list(A[["n.cores"]])
      A["var.names"] <- list(A[["var.names"]])
      A["offset"] <- list(NULL)
      A[["nTrain"]] <- floor(nrow(covs))
      A[["class.stratify.cv"]] <- FALSE
      A[["y"]] <- treat
      A[["x"]] <- covs
      A[["distribution"]] <- list(name = A[["distribution"]])
      A[["w"]] <- s.weights

      gbmCrossVal.call <- as.call(c(list(quote(gbm::gbmCrossVal)),
                                    A[names(A) %in% setdiff(names(formals(gbm::gbmCrossVal)), tunable)],
                                    tune[i, tunable[tunable %in% names(formals(gbm::gbmCrossVal))]]))

      verbosely({
        cv.results <- eval(gbmCrossVal.call)
      }, verbose = verbose)

      best.tree.index <- which.min(cv.results$error)
      best.loss <- cv.results$error[best.tree.index]
      best.tree <- best.tree.index

      tune[[paste.("best", "error")]][i] <- best.loss
      tune[["best.tree"]][i] <- best.tree

      if (best.loss < current.best.loss) {
        best.fit <- fit
        best.gps <- gbm::predict.gbm(fit, n.trees = best.tree, newdata = covs)

        dens.denom <- densfun(treat - best.gps, s.weights)
        best.w <- dens.num/dens.denom

        # if (trim.at != 0) best.w <- suppressMessages(trim(best.w, at = trim.at, treat = treat))
        current.best.loss <- best.loss
        best.tune.index <- i

        tree.val <- data.frame(tree = seq_along(cv.results$error),
                               error = cv.results$error)

        info <- list(best.tree = best.tree,
                     tree.val = tree.val)
      }
    }

  }

  if (isTRUE(A[["use.kernel"]]) && isTRUE(A[["plot"]])) {
    d.n <- attr(dens.num, "density")
    dens.denom <- densfun(treat - best.gps, s.weights)
    d.d <- attr(dens.denom, "density")
    plot_density(d.n, d.d)
  }

  tune[tunable[lengths(A[tunable]) == 1]] <- NULL

  if (ncol(tune) > 2) {
    info[["tune"]] <- tune
    info[["best.tune"]] <- tune[best.tune.index,]
  }

  list(w = best.w, info = info, fit.obj = best.fit)
}
