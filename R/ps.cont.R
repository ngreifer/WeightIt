#GBM for continuous treatments adapted from Zhu, Coffman, & Ghosh (2015) code
ps.cont <- function(formula, data, n.trees = 20000, interaction.depth = 4, shrinkage = 0.0005, bag.fraction = 1,
                    print.level = 0, verbose = FALSE, stop.method, sampw = NULL, optimize = 1, use.kernel = FALSE, ...) {
  #Program checks

  A <- list(...)
  terms <- match.call()

  #Set up subfunctions
  # cor2z <- function(x) {return(.5 * log((1+x)/(1-x)))}
  F.aac.w <- function(i, data, t, covs, ps.model, ps.num, corr.type, mean.max,
                      # z.trans,
                      s.weights) {
    GBM.fitted <- predict(ps.model, newdata = data, n.trees = floor(i),
                          type = "response")
    ps.den <- dnorm((t - GBM.fitted)/sd(t - GBM.fitted), 0, 1)
    wt <- ps.num/ps.den

    if (corr.type == "spearman") {
      ranked.t <- rank(t)
      corr_ <- col.w.r(apply(covs, 2, rank), ranked.t, w = wt, s.weights = s.weights)
      # corr_ <- apply(covs, 2, function(c) w.r(rank(c), y = ranked.t, w = wt, s.weights = s.weights))
    }
    else if (corr.type == "pearson") {
      corr_ <- col.w.r(covs, t, w = wt, s.weights = s.weights)

      # corr_ <- apply(covs, 2, function(c) {
      # w.r(c, y = t, w = wt, s.weights = s.weights)})
    }
    else stop("stop.method is not correctly specified.", call. = FALSE)

    # if (z.trans) corr_ <- cor2z(corr_)

    return(mean.max(abs(corr_)))
  }

  desc.wts.cont <- function(t, covs, weights, which.tree) {
    desc <- setNames(vector("list", 10),
                     c("ess", "n", "max.p.cor", "mean.p.cor", "rms.p.cor", "max.s.cor", "mean.s.cor", "rms.s.cor", "bal.tab", "n.trees"))
    desc[["bal.tab"]][["results"]] <- data.frame(p.cor = col.w.r(covs, t, w = weights, s.weights = s.weights),
                                                 # p.cor.z = NA,
                                                 s.cor = col.w.r(apply(covs, 2, rank), rank(t), w = weights, s.weights = s.weights),
                                                 # s.cor.z = NA,
                                                 row.names = colnames(covs))
    # desc[["bal.tab"]][["results"]][["p.cor.z"]] <- cor2z(desc[["bal.tab"]][["results"]][["p.cor"]])
    # desc[["bal.tab"]][["results"]][["s.cor.z"]] <- cor2z(desc[["bal.tab"]][["results"]][["s.cor"]])
    desc[["ess"]] <- ESS(weights)
    desc[["n"]] <- length(t)
    desc[["max.p.cor"]] <- max(abs(desc[["bal.tab"]][["results"]][["p.cor"]]))
    desc[["mean.p.cor"]] <- mean(abs(desc[["bal.tab"]][["results"]][["p.cor"]]))
    desc[["rms.p.cor"]] <- sqrt(mean(desc[["bal.tab"]][["results"]][["p.cor"]]^2))
    desc[["max.s.cor"]] <- max(abs(desc[["bal.tab"]][["results"]][["s.cor"]]))
    desc[["mean.s.cor"]] <- mean(abs(desc[["bal.tab"]][["results"]][["s.cor"]]))
    desc[["rms.s.cor"]] <- sqrt(mean(desc[["bal.tab"]][["results"]][["s.cor"]]^2))
    desc[["n.trees"]] <- which.tree

    return(desc)
  }

  check.package("gbm")

  if (missing(stop.method)) {
    warning("No stop.method was entered. Using \"p.mean\", the mean of the absolute Pearson correlations.", call. = FALSE, immediate. = TRUE)
    stop.method <- "p.mean"
  }
  else {
    stop.method <- tryCatch(match.arg(tolower(stop.method), apply(expand.grid(c("s", "p"),
                                                                              c(".mean", ".max", ".rms"),
                                                                              c(""
                                                                                # , ".z"
                                                                                )), 1, paste, collapse = ""),
                                      several.ok = TRUE),
                            error = function(e) {
                              stop("The entered stop.method is not one of the accepted values.\n\tSee ?ps.cont for the accepted values of stop.method for continuous treatments.",
                                      call. = FALSE)
                            })
  }

  stop.method.split <- strsplit(stop.method, ".", fixed = TRUE)
  corr.type <- sapply(stop.method.split, function(x) switch(x[1], s = "spearman", p = "pearson", k = "kendall"))
  mean.max <- lapply(stop.method.split, function(x) switch(x[2],
                                                           mean = function(y) mean(abs(y), na.rm = TRUE),
                                                           max = function(y) max(abs(y), na.rm = TRUE),
                                                           rms = function(y) sqrt(mean(y^2, na.rm = TRUE))))
  # z.trans <- sapply(stop.method.split, function(x) length(x) == 3 && x[3] == "z")

  if (missing(formula)) stop("A formula must be supplied.", call. = FALSE)
  if (missing(data)) data <- NULL

  t.c <- get.covs.and.treat.from.formula(formula, data)
  t <- t.c[["treat"]]
  covs <- t.c[["model.covs"]]
  #covs <- apply(covs, 2, function(x) if (is.factor(x) || is.character(x) || !nunique.gt(x, 2)) factor(x))
  new.data <- data.frame(t = t, t.c[["reported.covs"]])

  if (any(is.na(t))) stop("Missingness is not allowed in the treatment variable.", call. = FALSE)

  if (is_null(sampw)) {
    if (is_not_null(A[["s.weights"]])) s.weights <- A[["s.weights"]]
    else s.weights <- rep(1, length(t))
  }
  else s.weights <- sampw

  if (verbose) cat("Estimating marginal density of the treatment ")
  num.fit <- glm(t ~ 1,
                 data = new.data,
                 weights = s.weights)
  p.num <- t - num.fit$fitted.values

  if (isTRUE(use.kernel)) {
    if (verbose) cat("using kernel density estimation\n")
    if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
    if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
    if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
    if (is_null(A[["n"]])) A[["n"]] <- 10*length(t)

    d.n <- density(p.num, n = A[["n"]],
                   weights = s.weights/sum(s.weights), give.Rkern = FALSE,
                   bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
    ps.num <- with(d.n, approxfun(x = x, y = y))(p.num)
    #if anyNA(ps.num) stop("NAs were estimated in the numerator densty")
    if (isTRUE(A[["plot"]])) {
      d.n_ <- cbind(as.data.frame(d.n[c("x", "y")]), dens = "Numerator Density", stringsAsfactors = FALSE)
      d.all <- rbind(d.n_)
      d.all$dens <- factor(d.all$dens, levels = c("Numerator Density"))
      pl <- ggplot(d.all, aes(x=x,y=y)) + geom_line() + labs(title = "Weight Component Densities", x = "E[Treat|X]", y = "Density") +
        facet_grid(rows = vars(dens)) + theme(panel.background = element_rect(fill = "white"),
                                              panel.border = element_rect(fill = NA, color = "black"),
                                              axis.text.x = element_text(color = "black"),
                                              axis.text.y = element_text(color = "black"),
                                              panel.grid.major = element_blank(),
                                              panel.grid.minor = element_blank()
        )
      print(pl)
    }
  }
  else {
    if (verbose) cat("using a normal approximation\n")
    ps.num <- dnorm(p.num, 0, sqrt(summary(num.fit)$dispersion))
  }

  if (verbose) cat("Fitting gbm model\n")
  model.den <- gbm::gbm(formula(new.data),
                        data = new.data,
                        shrinkage = shrinkage,
                        interaction.depth = interaction.depth,
                        distribution = "gaussian",
                        n.trees = n.trees,
                        bag.fraction = bag.fraction,
                        n.minobsinnode = 10,
                        train.fraction = 1,
                        keep.data = FALSE,
                        verbose = FALSE,
                        weights = s.weights)

  desc <- setNames(vector("list", 1 + length(stop.method)),
                   c("unw", stop.method))

  if (verbose) cat("Diagnosis of unweighted analysis\n")
  desc[["unw"]] <- desc.wts.cont(t, covs, s.weights, NA)

  if (optimize == 1) {
    w <- ps <- setNames(as.data.frame(matrix(0, nrow = nrow(covs), ncol = length(stop.method))), stop.method)
    best.tree <- setNames(numeric(length(stop.method)), stop.method)
    nintervals <- round(1+sqrt(2*n.trees))
    iters <- round(seq(1, n.trees, length = nintervals))

    balance <- matrix(NA, ncol = length(stop.method), nrow = nintervals,
                      dimnames = list(iters, stop.method))

    for (s in seq_along(stop.method)) {

      sm <- stop.method[s]

      if (verbose) cat("Optimizing with", sm,"stopping rule\n")

      # get optimal number of iterations
      # Step #1: evaluate at (round(1+sqrt(2*n.trees))) equally spaced points
      for (j in 1:length(iters)) {

        balance[j, sm] <- F.aac.w(iters[j], data = new.data, t = t, covs = covs,
                                  ps.model = model.den,
                                  ps.num = ps.num, corr.type = corr.type[s], mean.max = mean.max[[s]],
                                  # z.trans = z.trans[s],
                                  s.weights = s.weights)
      }

      # Step #2: find the interval containing the approximate minimum
      interval <- which.min(balance[,sm]) + c(-1,1)
      interval[1] <- max(1, interval[1])
      interval[2] <- min(length(iters), interval[2])

      # Step #3: refine the minimum by searching with the identified interval

      opt <- optimize(F.aac.w, interval = iters[interval], data = new.data, t = t, covs = covs,
                      ps.model = model.den,
                      ps.num = ps.num, corr.type = corr.type[s], mean.max = mean.max[[s]],
                      # z.trans = z.trans[s],
                      s.weights = s.weights, tol = .Machine$double.eps)

      best.tree[sm] <- floor(opt$minimum)

      # compute propensity score weights
      GBM.fitted <- predict(model.den, newdata = new.data, n.trees = floor(best.tree[sm]),
                            type = "response")
      ps[[sm]] <- dnorm((t - GBM.fitted)/sd(t - GBM.fitted), 0, 1)
      w[[sm]] <- ps.num/ps[[s]]
    }
  }
  else if (optimize == 2) {
    w <- ps <- setNames(as.data.frame(matrix(0, nrow = nrow(covs), ncol = length(stop.method))), stop.method)
    best.tree <- setNames(numeric(length(stop.method)), stop.method)
    for (s in seq_along(stop.method)) {
      if (verbose) cat("Optimizing with", stop.method[s],"stopping rule\n")

      opt <- optimize(F.aac.w, interval = c(1, n.trees), data = new.data, t = t, covs = covs,
                      ps.model = model.den,
                      ps.num = ps.num, corr.type = corr.type[s], mean.max = mean.max[[s]],
                      # z.trans = z.trans[s],
                      s.weights = s.weights, tol = .Machine$double.eps)
      best.tree[s] <- floor(opt$minimum)

      GBM.fitted <- predict(model.den, newdata = new.data, n.trees = floor(best.tree[s]),
                            type = "response")
      ps[[s]] <- dnorm((t - GBM.fitted)/sd(t - GBM.fitted), 0, 1)
      w[[s]] <- ps.num/ps[[s]]
    }
    balance <- NULL
    iters <- NULL

  }
  else {
    wt <- ps.den <- vector("list", n.trees)
    iters <- 1:n.trees
    balance <- matrix(NA, ncol = length(stop.method), nrow = n.trees,
                      dimnames = list(iters, stop.method))
    corr_ <- array(NA, dim = c(n.trees, ncol(covs), length(stop.method)),
                    dimnames = list(iters, colnames(covs), stop.method))
    ps <- predict(model.den, newdata = data,
                               n.trees = iters, type = "response")
    for (i in iters) {
      # Calculate the inverse probability weights
      model.den$fitted = ps[,i]
      ps.den[[i]] = dnorm((t - model.den$fitted)/sd(t - model.den$fitted), 0, 1)
      wt[[i]] <- ps.num/ps.den[[i]]

      #bal[[i]] <- setNames(vector('list', length(stop.method)), stop.method)
      for (s in seq_along(stop.method)) {
        if (s > 1 && corr.type[s] %in% corr.type[1:(s-1)]) {
          corr_[i, , s] <- corr_[i, , which(corr.type == corr.type[s])[1]]
        }
        else {
          if (corr.type[s] == "spearman") {
            ranked.t <- rank(t)
            corr_[i, , s] <- col.w.r(apply(covs, 2, rank), ranked.t, w = wt[[i]], s.weights = s.weights)
            # corr_[i, , s] <- apply(covs, 2, function(c) {
            #   w.r(rank(c), y = ranked.t, w = wt[[i]], s.weights = s.weights)
            #   })
          }
          else if (corr.type[s] == "pearson") {
            corr_[i, , s] <- col.w.r(covs, t, w = wt[[i]], s.weights = s.weights)
            # corr_[i, , s] <- apply(covs, 2, function(c) {
            #   w.r(c, y = t, w = wt[[i]], s.weights = s.weights)})
          }
          else stop("stop.method is not correctly specified.", call. = FALSE)

          #bal[[i]][[s]] <- setNames(corr_, colnames(covs))

          # if (z.trans[s]) {
          #   balance[i, s] <- mean.max[[s]](cor2z(corr_[i, , s]))
          # }
          # else {
            balance[i, s] <- mean.max[[s]](corr_[i, , s])
          # }
        }
      }

    }

    w <- ps <- setNames(as.data.frame(matrix(0, nrow = nrow(covs), ncol = length(stop.method))), stop.method)
    best.tree <- setNames(numeric(length(stop.method)), stop.method)
    for (s in seq_along(stop.method)) {
      if (verbose) cat("Optimizing with", stop.method[s],"stopping rule\n")
      best.tree[s] <- floor(which.min(balance[,s]))

      ps[[s]] <- ps.den[[best.tree[s]]]
      w[[s]] <- wt[[best.tree[s]]]
    }

  }

  if (any(n.trees - best.tree < 100)) warning("Optimal number of iterations is close to the specified n.trees. n.trees is likely set too small and better balance might be obtainable by setting n.trees to be larger.", call. = FALSE)

  for (sm in stop.method) {
    if (verbose) cat("Diagnosis of",sm,"weights\n")

    w[[sm]] <- w[[sm]]*s.weights
    desc[[sm]] <- desc.wts.cont(t, covs, w[[sm]], best.tree[sm])
  }

  out <- list(treat = t, desc = desc, ps = ps, w = w,
              sampw = sampw, estimand = NULL, datestamp = date(),
              parameters = terms, alerts = NULL, iters = iters, balance = balance,
              n.trees = n.trees, data = data, gbm.obj = model.den)

  class(out) <- c("ps.cont")

  return(out)

}

summary.ps.cont <- function (object, ...) {
  summary.tab <- NULL
  typ <- NULL
  n.tp <- length(object$desc)
  for (i.tp in 1:n.tp) {
    desc.temp <- object$desc[[i.tp]]
    iter <- desc.temp$n.trees
    tp <- names(object$desc)[i.tp]
    summary.tab <- rbind(summary.tab, with(desc.temp, c(n, ess, max.p.cor, mean.p.cor, rms.p.cor,
                                                        max.s.cor, mean.s.cor, rms.s.cor, iter)))
    typ <- c(typ, tp)
  }
  summary.tab <- matrix(summary.tab, nrow = n.tp)
  rownames(summary.tab) <- typ
  colnames(summary.tab) <- c("n", "ess", "max.p.cor", "mean.p.cor", "rms.p.cor",
                             "max.s.cor", "mean.s.cor", "rms.s.cor", "iter")
  class(summary.tab) <- "summary.ps"
  return(summary.tab)
}

plot.ps.cont <- function(x, ...) {
  class(x) <- "ps"
  plot(x, ...)
}

boxplot.ps.cont <- function(x, ...) {
  class(x) <- "ps"
  boxplot(x, ...)
}