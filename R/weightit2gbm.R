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
#'        \item{`use.offset`}{`logical`; whether to use the linear predictor resulting from a generalized linear model as an offset to the GBM model. If `TRUE`, this fits a logistic regression model (for binary treatments) or a linear regression model (for continuous treatments) and supplies the linear predict to the `offset` argument of `gbm.fit()`. This often improves performance generally but especially when the true propensity score model is well approximated by a GLM, and this yields uniformly superior performance over `method = "glm"` with respect to `criterion`. Default is `FALSE` to omit the offset. Only allowed for binary and continuous treatments. This argument is tunable.
#'        }
#' }
#'
#' All other arguments take on the defaults of those in \pkgfun{gbm}{gbm.fit}, and some are not used at all. For binary and multi-category treatments with a with cross-validation used as the criterion, `class.stratify.cv` is set to `TRUE` by default.
#'
#' The `w` argument in `gbm.fit()` is ignored because sampling weights are passed using `s.weights`.
#'
#' For continuous treatments only, the following arguments may be supplied:
#' \describe{
#'       \item{`density`}{A function corresponding to the conditional density of the treatment. The standardized residuals of the treatment model will be fed through this function to produce the numerator and denominator of the generalized propensity score weights. This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if `density = "dt_2"` is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme outcome values are observed (Naimi et al., 2014).
#'
#' Can also be `"kernel"` to use kernel density estimation, which calls [density()] to estimate the numerator and denominator densities for the weights. (This used to be requested by setting `use.kernel = TRUE`, which is now deprecated.)
#'
#' If unspecified, a density corresponding to the argument passed to `distribution`. If `"gaussian"` (the default), [dnorm()] is used. If `"tdist"`, a t-distribution with 4 degrees of freedom is used. If `"laplace"`, a laplace distribution is used.}
#'       \item{`bw`, `adjust`, `kernel`, `n`}{If `density = "kernel"`, the arguments to [density()]. The defaults are the same as those in `density()` except that `n` is 10 times the number of units in the sample.}
#'       \item{`plot`}{If `density = "kernel"`, whether to plot the estimated densities.}
#' }
#'
#' For tunable arguments, multiple entries may be supplied, and `weightit()` will choose the best value by optimizing the criterion specified in `criterion`. See below for additional outputs that are included when arguments are supplied to be tuned. See Examples for an example of tuning. The same seed is used for every run to ensure any variation in performance across tuning parameters is due to the specification and not to using a random seed. This only matters when `bag.fraction` differs from 1 (its default) or cross-validation is used as the criterion; otherwise, there are no random components in the model.
#'
#' @section Additional Outputs:
#' \describe{
#' \item{`info`}{
#'   A list with the following entries:
#'     \describe{
#'       \item{`best.tree`}{
#'         The number of trees at the optimum. If this is close to `n.trees`, `weightit()` should be rerun with a larger value for `n.trees`, and `start.tree` can be set to just below `best.tree`. When other parameters are tuned, this is the best tree value in the best combination of tuned parameters. See example.}
#'       \item{`tree.val`}{
#'         A data frame with two columns: the first is the number of trees and the second is the value of the criterion corresponding to that tree. Running [plot()] on this object will plot the criterion by the number of trees and is a good way to see patterns in the relationship between them and to determine if more trees are needed. When other parameters are tuned, these are the number of trees and the criterion values in the best combination of tuned parameters. See example.}
#'     }
#'   If any arguments are to be tuned (i.e., they have been supplied more than one value), the following two additional components are included in `info`:
#'     \describe{
#'       \item{`tune`}{
#'         A data frame with a column for each argument being tuned, the best value of the balance criterion for the given combination of parameters, and the number of trees at which the best value was reached.}
#'       \item{`best.tune`}{
#'         A one-row data frame containing the values of the arguments being tuned that were ultimately selected to estimate the returned weights.}
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
#' `plot()` can be used on the output of `weightit()` with `method = "gbm"` to display the results of the tuning process; see Examples and [plot.weightit()] for more details.
#'
#' @note
#' The `criterion` argument used to be called `stop.method`, which is its name in \pkg{twang}. `stop.method` still works for backward compatibility. Additionally, the criteria formerly named as `"es.mean"`, `"es.max"`, and `"es.rms"` have been renamed to `"smd.mean"`, `"smd.max"`, and `"smd.rms"`. The former are used in \pkg{twang} and will still work with `weightit()` for backward compatibility.
#'
#' Estimated propensity scores are trimmed to \eqn{10^{-8}} and \eqn{1 - 10^{-8}} to ensure balance statistics can be computed.
#'
#' @seealso
#' [weightit()], [weightitMSM()]
#'
#' \pkgfun{gbm}{gbm.fit} for the fitting function.
#'
#' @references
#' ## Binary treatments
#'
#' McCaffrey, D. F., Ridgeway, G., & Morral, A. R. (2004). Propensity Score Estimation With Boosted Regression for Evaluating Causal Effects in Observational Studies. *Psychological Methods*, 9(4), 403–425. \doi{10.1037/1082-989X.9.4.403}
#'
#' ## Multi-Category Treatments
#'
#' McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. *Statistics in Medicine*, 32(19), 3388–3414. \doi{10.1002/sim.5753}
#'
#' ## Continuous treatments
#'
#' Zhu, Y., Coffman, D. L., & Ghosh, D. (2015). A Boosting Algorithm for Estimating Generalized Propensity Scores with Continuous Treatments. *Journal of Causal Inference*, 3(1). \doi{10.1515/jci-2014-0022}
#'
#' @examplesIf requireNamespace("gbm", quietly = TRUE)
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "gbm", estimand = "ATE",
#'                 criterion = "smd.max",
#'                 use.offset = TRUE))
#' summary(W1)
#' bal.tab(W1)
#'
#' # View information about the fitting process
#' W1$info$best.tree #best tree
#' plot(W1) #plot of criterion value against number of trees
#'
#' \donttest{
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
#'                   method = "gbm", density = "kernel",
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
#'   plot(W4a) #decreasing at right edge
#'
#'   W4b <- weightit(re75 ~ age + educ + married +
#'                     nodegree + re74, data = lalonde,
#'                   method = "gbm", density = "dt_3",
#'                   criterion = "p.max",
#'                   start.tree = 10000,
#'                   n.trees = 20000)
#'
#'   W4b$info$best.tree #13417; optimum has been found
#'   plot(W4b) #increasing at right edge
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
#'   plot(W5) #plot criterion values against number of trees
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

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs, replace_with = NA)
  }

  criterion <- A[["criterion"]]
  if (is_null(criterion)) {
    criterion <- A[["stop.method"]]
  }

  if (is_null(criterion)) {
    .wrn("no `criterion` was provided. Using \"smd.mean\"")
    criterion <- "smd.mean"
  }
  else {
    chk::chk_string(criterion)
  }

  available.criteria <- cobalt::available.stats(switch(treat.type, "multi-category" = "multi", treat.type))

  if (startsWith(criterion, "es.")) {
    subbed.crit <- sub("es.", "smd.", criterion, fixed = TRUE)
    subbed.match <- charmatch(subbed.crit, available.criteria)
    if (!anyNA(subbed.match) && subbed.match != 0L) {
      criterion <- subbed.crit
    }
  }

  cv <- 0

  s.m.matches <- charmatch(criterion, available.criteria)
  if (anyNA(s.m.matches) || s.m.matches == 0L) {
    if (!startsWith(criterion, "cv") ||
        !can_str2num(numcv <- substr(criterion, 3, nchar(criterion)))) {
      .err(sprintf("`criterion` must be one of %s",
                   word_list(c(available.criteria, "cv{#}"), "or", quotes = TRUE)))
    }

    cv <- round(str2num(numcv))

    if (cv < 2)
      .err("at least 2 CV-folds must be specified in `criterion`")
  }
  else {
    criterion <- available.criteria[s.m.matches]
  }

  tunable <- c("interaction.depth", "shrinkage", "distribution", "use.offset")

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
                               class.stratify.cv = if (cv > 0) TRUE else NULL,
                               formals(gbm::gbm.fit)[[f]]))
    }
  }

  n.trees <- A[["n.trees"]]
  chk::chk_count(n.trees)
  chk::chk_gt(n.trees, 1)

  if (treat.type == "binary")  {
    available.distributions <- c("bernoulli", "adaboost")
    t.lev <- get_treated_level(treat)
    treat <- binarize(treat, one = t.lev)
  }
  else {
    available.distributions <- "multinomial"
    treat <- factor(treat)
  }

  if (cv == 0) {
    start.tree <- if_null_then(A[["start.tree"]], 1)
    chk::chk_count(start.tree)
    chk::chk_range(start.tree, c(1, n.trees))

    n.grid <- if_null_then(A[["n.grid"]],
                           round(1 + sqrt(2 * (n.trees - start.tree + 1))))
    chk::chk_count(n.grid)
    chk::chk_range(n.grid, c(2, n.trees))

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
  A[["distribution"]] <- {
    if (is_null(distribution <- A[["distribution"]])) available.distributions[1]
    else match_arg(distribution, available.distributions, several.ok = TRUE)
  }
  A[["w"]] <- s.weights
  A[["verbose"]] <- FALSE
  A[["n.trees"]] <- n.trees

  # Offset
  if (is_not_null(A[["use.offset"]])) {
    chk::chk_logical(A[["use.offset"]], "`use.offset`")
    chk::chk_not_any_na(A[["use.offset"]], "`use.offset`")
  }

  if (any(A[["use.offset"]])) {
    if (treat.type == "multi-category") {
      .err("`use.offset` cannot be used with multi-category treatments")
    }

    if (!identical(A[["distribution"]], "bernoulli")) {
      .err("`use.offset` can only be used with `distribution = \"bernoulli\"`")
    }

    fit <- glm.fit(x = as.matrix(cbind(1, covs)), y = treat,
                   weights = s.weights, family = quasibinomial())
    offset <- fit$linear.predictors
  }
  else {
    A[["use.offset"]] <- FALSE
    offset <- NULL
  }

  tune <- do.call("expand.grid", c(A[names(A) %in% tunable],
                                   list(stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)))

  tune <- unique(tune)

  info <- list()

  current.best.loss <- Inf

  # Maintain seed across tunable params
  genv <- globalenv()
  if (is_null(genv$.Random.seed)) runif(1)
  curr_seed <- genv$.Random.seed

  for (i in seq_row(tune)) {
    assign(".Random.seed", value = curr_seed, envir = genv)

    use.offset <- tune[["use.offset"]][i]
    A[["offset"]] <- if (use.offset) offset else NULL
    A[["distribution"]] <- list(name = tune[["distribution"]][i])
    tune_args <- as.list(tune[i, setdiff(tunable, c("distribution", "use.offset"))])

    gbm.call <- as.call(c(list(quote(gbm::gbm.fit)),
                          A[names(A) %in% setdiff(names(formals(gbm::gbm.fit)), names(tune_args))],
                          tune_args))
    verbosely({
      fit <- eval(gbm.call)
    }, verbose = verbose)

    if (cv == 0) {

      n.trees <- fit[["n.trees"]]
      iters <- seq_len(n.trees)
      iters.grid <- unique(round(seq(start.tree, n.trees, length.out = n.grid)))

      if (is_null(iters.grid) || anyNA(iters.grid) || any(iters.grid > n.trees)) {
        .err("a problem has occurred")
      }

      ps <- {
        if (use.offset) plogis(A[["offset"]] + gbm::predict.gbm(fit, n.trees = iters.grid,
                                                                type = "link", newdata = covs))
        else gbm::predict.gbm(fit, n.trees = iters.grid,
                              type = "response", newdata = covs)
      }

      w <- .get_w_from_ps_internal_array(ps, treat = treat, estimand = estimand,
                                         focal = focal, stabilize = stabilize, subclass = subclass)
      if (trim.at != 0) {
        w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))
      }

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

        ps <- {
          if (use.offset) plogis(A[["offset"]] + gbm::predict.gbm(fit, n.trees = iters.to.check,
                                                                  type = "link", newdata = covs))
          else gbm::predict.gbm(fit, n.trees = iters.to.check,
                                type = "response", newdata = covs)
        }

        w <- .get_w_from_ps_internal_array(ps, treat = treat, estimand = estimand,
                                           focal = focal, stabilize = stabilize, subclass = subclass)
        if (trim.at != 0) {
          w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))
        }

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
      }

      info <- list(best.tree = c(info$best.tree, setNames(best.tree, i)),
                   tree.val = rbind(info$tree.val, cbind(tune = i, tree.val)))
    }
    else {
      if (i == 1) {
        A["data"] <- list(data.frame(treat, covs))
        A[["cv.folds"]] <- cv
        A["n.cores"] <- list(A[["n.cores"]])
        A["var.names"] <- list(A[["var.names"]])
        A[["nTrain"]] <- nrow(covs)
      }

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

        best.ps <- {
          if (use.offset) plogis(A[["offset"]] + gbm::predict.gbm(best.fit, n.trees = best.tree,
                                                                  type = "link", newdata = covs))
          else gbm::predict.gbm(best.fit, n.trees = best.tree,
                                type = "response", newdata = covs)
        }

        best.w <- drop(.get_w_from_ps_internal_array(best.ps, treat = treat, estimand = estimand,
                                                     focal = focal, stabilize = stabilize, subclass = subclass))
        # if (trim.at != 0) best.w <- suppressMessages(trim(best.w, at = trim.at, treat = treat))
        current.best.loss <- best.loss
        best.tune.index <- i

        tree.val <- data.frame(tree = seq_along(cv.results$error),
                               error = cv.results$error)

        info <- list(best.tree = best.tree,
                     tree.val = tree.val)

        if (treat.type == "multi-category") best.ps <- NULL
      }
    }

    if (treat.type == "multi-category") ps <- NULL
  }

  if (nrow(tune) > 1) {
    info[["tune"]] <- tune
    info[["best.tune"]] <- tune[best.tune.index,]
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

  for (i in seq_col(covs)) covs[,i] <- .make_closer_to_1(covs[,i])

  if (missing == "ind") {
    covs <- add_missing_indicators(covs, replace_with = NA)
  }

  criterion <- A[["criterion"]]
  if (is_null(criterion)) {
    criterion <- A[["stop.method"]]
  }

  if (is_null(criterion)) {
    .wrn("no `criterion` was provided. Using \"p.mean\"")
    criterion <- "p.mean"
  }
  else {
    chk::chk_string(criterion)
  }

  available.criteria <- cobalt::available.stats("continuous")

  cv <- 0

  s.m.matches <- charmatch(criterion, available.criteria)
  if (anyNA(s.m.matches) || s.m.matches == 0L) {
    if (!startsWith(criterion, "cv") ||
        !can_str2num(numcv <- substr(criterion, 3, nchar(criterion)))) {
      .err(sprintf("`criterion` must be one of %s",
                   word_list(c(available.criteria, "cv{#}"), "or", quotes = TRUE)))
    }

    cv <- round(str2num(numcv))

    if (cv < 2)
      .err("at least 2 CV-folds must be specified in `criterion`")
  }
  else {
    criterion <- available.criteria[s.m.matches]
  }

  tunable <- c("interaction.depth", "shrinkage", "use.offset", "distribution")

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

  n.trees <- A[["n.trees"]]
  chk::chk_count(n.trees)
  chk::chk_gt(n.trees, 1)

  available.distributions <- c("gaussian", "laplace", "tdist")

  if (cv == 0) {
    start.tree <- if_null_then(A[["start.tree"]], 1)
    chk::chk_count(start.tree)
    chk::chk_range(start.tree, c(1, n.trees))

    n.grid <- if_null_then(A[["n.grid"]],
                           round(1 + sqrt(2 * (n.trees - start.tree + 1))))
    chk::chk_count(n.grid)
    chk::chk_range(n.grid, c(2, n.trees))

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
  A[["n.trees"]] <- n.trees

  # Offset
  if (is_not_null(A[["use.offset"]])) {
    chk::chk_logical(A[["use.offset"]], "`use.offset`")
    chk::chk_not_any_na(A[["use.offset"]], "`use.offset`")
  }

  if (any(A[["use.offset"]])) {
    fit <- lm.wfit(x = as.matrix(cbind(1, covs)), y = treat,
                   w = s.weights)
    offset <- fit$fitted.values
  }
  else {
    A[["use.offset"]] <- FALSE
    offset <- NULL
  }

  tune <- do.call("expand.grid", c(A[names(A) %in% tunable],
                                   list(stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)))

  tune <- unique(tune)

  null_density <- !isTRUE(A[["use.kernel"]]) && is_null(A[["density"]])

  info <- list()
  current.best.loss <- Inf

  # Maintain seed across tunable params
  genv <- globalenv()
  if (is_null(genv$.Random.seed)) runif(1)
  curr_seed <- genv$.Random.seed

  for (i in seq_row(tune)) {
    assign(".Random.seed", value = curr_seed, envir = genv)

    use.offset <- tune[["use.offset"]][i]
    A["offset"] <- list(if (use.offset) offset else NULL)
    A[["distribution"]] <- list(name = tune[["distribution"]][i])
    tune_args <- as.list(tune[i, setdiff(tunable, c("distribution", "use.offset"))])

    if (null_density) {
      A[["use.kernel"]] <- FALSE
      A[["density"]] <- switch(A[["distribution"]]$name,
                               "gaussian" = "dnorm",
                               "tdist" = "dt_4",
                               "laplace" = function(x) exp(-abs(x))/2)
    }

    if (i == 1 || (null_density && !identical(tune[["distribution"]][i], tune[["distribution"]][i - 1]))) {
      #Process density params
      densfun <- .get_dens_fun(use.kernel = isTRUE(A[["use.kernel"]]), bw = A[["bw"]],
                               adjust = A[["adjust"]], kernel = A[["kernel"]],
                               n = A[["n"]], treat = treat, density = A[["density"]],
                               weights = s.weights)

      #Stabilization - get dens.num
      dens.num <- densfun(scale_w(treat, s.weights))
    }

    gbm.call <- as.call(c(list(quote(gbm::gbm.fit)),
                          A[names(A) %in% setdiff(names(formals(gbm::gbm.fit)), names(tune_args))],
                          tune_args))
    verbosely({
      fit <- eval(gbm.call)
    }, verbose = verbose)

    if (cv == 0) {

      n.trees <- fit[["n.trees"]]
      iters <- seq_len(n.trees)
      iters.grid <- unique(round(seq(start.tree, n.trees, length.out = n.grid)))

      if (is_null(iters.grid) || anyNA(iters.grid) || any(iters.grid > n.trees)) {
        .err("a problem has occurred")
      }

      gps <- gbm::predict.gbm(fit, n.trees = iters.grid, newdata = covs)
      if (use.offset) gps <- gps + A[["offset"]]

      w <- apply(gps, 2, function(p) {
        r <- treat - p
        dens.num / densfun(r / sqrt(col.w.v(r, s.weights)))
      })

      if (trim.at != 0) {
        w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))
      }

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
        if (use.offset) gps <- gps + A[["offset"]]

        w <- apply(gps, 2, function(p) {
          r <- treat - p
          dens.num / densfun(r / sqrt(col.w.v(r, s.weights)))
        })

        if (trim.at != 0) {
          w <- suppressMessages(apply(w, 2, trim, at = trim.at, treat = treat))
        }

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
      }

    }
    else {
      A["data"] <- list(data.frame(treat, covs))
      A[["cv.folds"]] <- cv
      A["n.cores"] <- list(A[["n.cores"]])
      A["var.names"] <- list(A[["var.names"]])
      A[["nTrain"]] <- floor(nrow(covs))
      A[["class.stratify.cv"]] <- FALSE

      gbmCrossVal.call <- as.call(c(list(quote(gbm::gbmCrossVal)),
                                    A[names(A) %in% setdiff(names(formals(gbm::gbmCrossVal)), names(tune_args))],
                                    tune_args))

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
        if (use.offset) best.gps <- best.gps + A[["offset"]]

        r <- treat - best.gps
        dens.denom <- densfun(r / sqrt(col.w.v(r, s.weights)))
        best.w <- dens.num / dens.denom

        current.best.loss <- best.loss
        best.tune.index <- i
      }

      tree.val <- data.frame(tree = seq_along(cv.results$error),
                             error = cv.results$error)
    }

    info <- list(best.tree = c(info$best.tree, setNames(best.tree, i)),
                 tree.val = rbind(info$tree.val, cbind(tune = i, tree.val)))
  }

  if (isTRUE(A[["plot"]])) {
    d.n <- attr(dens.num, "density")
    r <- treat - best.gps
    dens.denom <- densfun(r / sqrt(col.w.v(r, s.weights)))
    d.d <- attr(dens.denom, "density")
    plot_density(d.n, d.d)
  }

  if (nrow(tune) > 1) {
    info[["tune"]] <- tune
    info[["best.tune"]] <- tune[best.tune.index,]
  }

  list(w = best.w, info = info, fit.obj = best.fit)
}

.plot_tune_gbm <- function(info, by = NULL) {

  use.by <- is_not_null(by)
  use.tune <- {
    if (use.by) is_not_null(info[[1]]$tune)
    else is_not_null(info$tune)
  }

  shape <- "diamond"
  size <- 3
  subsample <- 5000

  if (use.by) {
    d <- do.call("rbind", lapply(names(info), function(i) {
      cbind(info[[i]]$tree.val, by = i)
    }))

    d$by <- factor(d$by, levels = names(info))

    criterion <- names(d)[3]

    if (use.tune) {
      best <- do.call("rbind", lapply(names(info), function(i) {
        data.frame(tune = factor(seq_along(info[[i]]$best.tree)),
                   best.tree = info[[i]]$best.tree,
                   y = info[[i]]$tune[[paste.("best", criterion)]],
                   by = i)
      }))

      d$tune <- factor(d$tune)

      #Subsample if too big
      ind <- unlist(lapply(split(seq_row(d), d[c("by", "tune")]), function(i) {
        if (length(i) <= subsample) return(i)
        b <- d$by[i][1]
        t <- d$tune[i][1]

        trees <- round(seq(min(d$tree[i]), max(d$tree[i]), length.out = round(subsample * .8)))
        trees <- c(trees, d$tree[i][d$tree[i] >= best$best.tree[best$by == b & best$tune == t] - subsample * .1 &
                                      d$tree[i] <= best$best.tree[best$by == b & best$tune == t] + subsample * .1])

        i[d$tree[i] %in% trees]
      }))

      d <- d[sort(ind),]
    }
    else {
      best <- do.call("rbind", lapply(names(info), function(i) {
        data.frame(best.tree = info[[i]]$best.tree,
                   y = info[[i]]$tree.val[[criterion]][info[[i]]$tree.val$tree == info[[i]]$best.tree],
                   by = i)
      }))

      #Subsample if too big
      ind <- unlist(lapply(split(seq_row(d), d["by"]), function(i) {
        if (length(i) <= subsample) return(i)
        b <- d$by[i][1]

        trees <- round(seq(min(d$tree[i]), max(d$tree[i]), length.out = round(subsample * .8)))
        trees <- c(trees, d$tree[i][d$tree[i] >= best$best.tree[best$by == b] - subsample * .1 &
                                      d$tree[i] <= best$best.tree[best$by == b] + subsample * .1])

        i[d$tree[i] %in% trees]
      }))

      d <- d[sort(ind),]
    }

    best$by <- factor(best$by, levels = names(info))
  }
  else {
    d <- info$tree.val
    criterion <- names(d)[3]

    if (use.tune) {
      best <- data.frame(tune = factor(seq_along(info$best.tree)),
                         best.tree = info$best.tree,
                         y = info$tune[[paste.("best", criterion)]])

      d$tune <- factor(d$tune)

      #Subsample if too big
      ind <- unlist(lapply(split(seq_row(d), d["tune"]), function(i) {
        if (length(i) <= subsample) return(i)
        t <- d$tune[i][1]

        trees <- round(seq(min(d$tree[i]), max(d$tree[i]), length.out = round(subsample * .8)))
        trees <- c(trees, d$tree[i][d$tree[i] >= best$best.tree[best$tune == t] - subsample * .1 &
                                      d$tree[i] <= best$best.tree[best$tune == t] + subsample * .1])

        i[d$tree[i] %in% trees]
      }))

      d <- d[sort(ind),]
    }
    else {
      best <- data.frame(best.tree = info$best.tree,
                         y = d[[criterion]][d$tree == info$best.tree])

      #Subsample if too big
      if (nrow(d) > subsample) {
        trees <- round(seq(min(d$tree), max(d$tree), length.out = round(subsample * .8)))
        trees <- c(trees, d$tree[d$tree >= best$best.tree - subsample * .1 &
                                   d$tree <= best$best.tree + subsample * .1])

        d <- d[d$tree %in% trees,]
      }
    }
  }

  p <- ggplot() +
    labs(x = "Tree", y = criterion)

  if (use.tune) {
    tune <- {
      if (use.by) info[[1]]$tune
      else info$tune
    }

    tune_args <- setdiff(names(tune), c(paste.("best", criterion), "best.tree"))
    tune_args <- tune_args[!vapply(tune[tune_args], all_the_same, logical(1L))]

    levels(d$tune) <- levels(best$tune) <- vapply(seq_row(tune), function(i) {
      do.call("paste", c(lapply(tune_args, function(a) {
        sprintf("%s = %s", a, add_quotes(tune[[a]][i], is.character(tune[[a]][i])))
      }), list(sep = ", ")))
    }, character(1L))

    p <- p +
      geom_line(data = d, aes(x = .data$tree, y = .data[[criterion]],
                              color = .data$tune)) +
      geom_point(data = best, aes(y = .data$y, x = .data$best.tree,
                                  color = .data$tune),
                 shape = shape, size = size) +
      labs(color = "Parameters") +
      guides(color = guide_legend(position = "bottom", ncol = 1, direction = "vertical"))
  }
  else {
    p <- p +
      geom_line(data = d, aes(x = .data$tree, y = .data[[criterion]])) +
      geom_point(data = best, aes(y = .data$y, x = .data$best.tree),
                 shape = shape, size = size)
  }

  if (use.by) {
    p <- p + facet_wrap(vars(.data$by))
  }

  p +
    theme_bw()
}

