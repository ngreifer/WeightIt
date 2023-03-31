ps.multi <- function(formula, data, n.trees = 20000, interaction.depth = 4, shrinkage = 0.001, bag.fraction = 1,
                     verbose = FALSE, estimand = "ATE", stop.method = "es.mean", sampw = NULL, treatATT = NULL) {

  terms <- match.call()

  #Process stop.method
  allowable.stop.methods <- expand.grid_string(c("es", "ks"), ".", c("mean", "max", "l2"))
  stop.method <- tolower(stop.method)
  if (stop.method %nin% allowable.stop.methods) stop(paste0("Stop method must be one of ", word.list(allowable.stop.methods, "or"), "."), call. = FALSE)
  else {
    split.stop.method <- strsplit(stop.method, ".", fixed = TRUE)[[1]]
    stop.stat <- split.stop.method[1]
    stop.fun <- switch(split.stop.method[2], "l2" = function(x) sqrt(mean(x^2)),
                       get(split.stop.method[2]))
  }

  if (stop.stat == "es") {
    bal <- function(w, stop.fun) {
      all.bal <- do.call(cobalt::bal.tab, list(formula, data = data,
                                 weights = as.data.frame(w), method = "w",
                                 estimand = estimand, binary = "std", focal = treatATT,
                                 quick = TRUE))[["Balance.Across.Pairs"]]
      all.es <- all.bal[,startsWith(names(all.bal), "Max.Diff.")][,-1]
      return(apply(all.es, 2, stop.fun))
    }
  }
  else if (stop.stat == "ks") {
    bal <- function(w, stop.fun) {
      all.bal <- cobalt::bal.tab(formula, data = data,
                                 weights = as.data.frame(w), method = "w",
                                 estimand = estimand, continuous = "raw", focal = treatATT,
                                 quick = TRUE, disp.ks = TRUE)[["Balance.Across.Pairs"]]
      all.ks <- all.bal[,startsWith(names(all.bal), "Max.KS.")][,-1]
      all.es <- all.bal[,startsWith(names(all.bal), "Max.Diff.")][,-1]
      all.ks[all.bal[["Type"]] == "Binary",] <- all.es[all.bal[["Type"]] == "Binary",]

      return(apply(all.ks, 2, stop.fun))
    }
  }

  t.c <- get.covs.and.treat.from.formula(formula, data)

  covs <- t.c[["model.covs"]]
  treat <- t.c[["treat"]]

  if (!nunique.gt(treat, 2)) stop("treat must have more than 2 unique values.", call. = FALSE)

  fit <- do.call(gbm::gbm.fit, list(x = covs,
                                    y = setNames(treat, 1:length(treat)),
                                    distribution = "multinomial",
                                    n.trees = n.trees,
                                    interaction.depth = interaction.depth,
                                    n.minobsinnode = 10,
                                    shrinkage = shrinkage,
                                    bag.fraction = bag.fraction,
                                    verbose = verbose,
                                    response.name = "treat",
                                    keep.data = FALSE))

  iters <- 1:n.trees
  iters.grid <- round(seq(1, n.trees, length.out = round(1+sqrt(2*n.trees))))

  if (any(is.na(iters.grid)) || length(iters.grid) == 0 || any(iters.grid > n.trees)) stop("A problem has occurred")

  ps <- gbm::predict.gbm(fit, newdata = covs, n.trees = iters.grid, type = "response")

  w.mat <- apply(ps, 3, get_w_from_ps, treat = treat, estimand = "ATE", focal = NULL)

  iter.grid.balance <- bal(w.mat, stop.fun)

  it <- which.min(iter.grid.balance) + c(-1, 1)
  it[1] <- iters.grid[max(1, it[1])]
  it[2] <- iters.grid[min(length(iters.grid), it[2])]
  iters.to.check <- iters[iters <= iters[it[2]] & iters >= iters[it[1]]]

  if (any(is.na(iters.to.check)) || length(iters.to.check) == 0 || any(iters.to.check > n.trees)) stop("A problem has occurred")

  ps <- gbm::predict.gbm(fit, newdata = covs, n.trees = iters.to.check, type = "response")

  w.mat <- apply(ps, 3, get_w_from_ps, treat = treat, estimand = "ATE", focal = NULL)

  iter.grid.balance.fine <- bal(w.mat, stop.fun)

  best.tree <- iters.to.check[which.min(iter.grid.balance.fine)]

  w <- w.mat[,as.character(best.tree)]

  out <- list(weights = w,
              formula = formula,
              data = data,
              estimand = estimand,
              sampw = sampw,
              best.tree = best.tree,
              treatATT = treatATT,
              call = terms)
  class(out) <- "ps.multi"

  return(out)

}