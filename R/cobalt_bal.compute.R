# A version of cobalt::bal.init() and cobalt::bal.compute() so that WeightIt
# can be a standalone package before cobalt is updated

bal.compute <- function(x,
                        weights = NULL,
                        ...) {
  fun <- attr(x, "fun")

  fun(init = x, weights = weights)
}

bal.init <- function(x,
                     treat,
                     stat,
                     s.weights = NULL,
                     ...) {
  chk::chk_not_missing(x, "`x`")
  chk::chk_not_missing(treat, "`treat`")
  chk::chk_not_missing(stat, "`stat`")

  chk::chk_string(stat)
  chk::chk_vector(treat)

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  stat <- match_arg(stat, available.stats(treat.type))

  .chk_null_or(s.weights, chk = chk::chk_numeric)

  init <- bal_criterion(treat.type, stat)

  out <- init$init(x = x, treat = treat, s.weights = s.weights, ...)

  attr(out, "fun") <- init$fun
  attr(out, "treat.type") <- attr(init, "treat.type")
  attr(out, "stat") <- attr(init, "stat")

  class(out) <- c("bal.init", class(out))

  out
}

available.stats <- function(treat.type = "binary") {
  chk::chk_string(treat.type)
  treat.type <- match_arg(treat.type, c("binary", "multinomial", "continuous"))

  criteria <- switch(
    treat.type,
    binary = c(
      "smd.mean",
      "smd.max",
      "smd.rms",
      "ks.mean",
      "ks.max",
      "ks.rms",
      "ovl.mean",
      "ovl.max",
      "ovl.rms",
      "mahalanobis",
      "energy.dist",
      "kernel.dist",
      "l1.med",
      "r2",
      "r2.2",
      "r2.3"
    ),
    multinomial = c(
      "smd.mean",
      "smd.max",
      "smd.rms",
      "ks.mean",
      "ks.max",
      "ks.rms",
      "ovl.mean",
      "ovl.max",
      "ovl.rms",
      "energy.dist",
      "l1.med"
    ),
    continuous = c(
      "p.mean",
      "p.max",
      "p.rms",
      "s.mean",
      "s.max",
      "s.rms",
      "r2",
      "r2.2",
      "r2.3",
      "distance.cov"
    )
  )

  criteria
}

process_init_covs <- function(covs) {
  nm <- deparse1(substitute(covs))
  needs.splitting <- FALSE
  if (!is.matrix(covs)) {
    if (is.data.frame(covs)) {
      if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
        needs.splitting <- TRUE
      }
      else covs <- as.matrix(covs)
    }
    else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
    else .err(sprintf("`%s` must be a data.frame or numeric matrix", nm))
  }
  else if (!is.numeric(covs)) .err(sprintf("`%s` must be a data.frame or numeric matrix", nm))

  bin.vars <- {
    if (is_null(attr(covs, "bin"))) process.bin.vars(mat = covs)
    else attr(covs, "bin")
  }

  if (needs.splitting) {
    bin.vars[to.split] <- TRUE
    covs <- do.call(cobalt::splitfactor, list(covs, drop.first ="if2",
                                              split.with = bin.vars))
    bin.vars <- attr(covs, "split.with")[[1]]
  }

  attr(covs, "bin") <- bin.vars
  covs
}

#init functions
init_smd <- function(x, treat, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE, ...) {
  chk::chk_flag(pairwise)

  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("binary", "multinomial")) {
    .err("`treat` must be a binary or multi-category variable")
  }

  .chk_null_or(estimand, chk::chk_string)
  if (is_null(estimand)) estimand <- "ATE"

  f.e <- process.focal.and.estimand(focal, estimand, treat)
  focal <- f.e[["focal"]]
  estimand <- f.e[["estimand"]]

  unique.treats <- unique(treat)

  if (treat.type == "multinomial") {
    if (is_null(focal) && !pairwise) {
      treat.all <- last(make.unique(unique.treats, "All"))
      treat <- factor(c(as.character(treat), rep(treat.all, length(treat))),
                      levels = c(unique.treats, treat.all))
      x <- rbind(x, x)
      s.weights <- rep(s.weights, 2)
      focal <- treat.all
    }

    if (is_null(focal) || pairwise) {
      treatment.pairs <- combn(unique.treats, 2, simplify = FALSE)
    }
    else {
      treatment.pairs <- lapply(setdiff(unique.treats, focal), c, focal)
    }
  }
  else {
    treatment.pairs <- list(unique.treats)
    pairwise <- TRUE
  }

  s.d.denom <- get.s.d.denom.weightit(estimand = estimand, treat = treat, focal = focal)

  denoms <- compute_s.d.denom(x, treat = treat,
                              s.d.denom = s.d.denom, s.weights = s.weights,
                              bin.vars = bin.vars)

  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              s.d.denom = denoms,
              focal = focal,
              pairwise = pairwise,
              treatment.pairs = treatment.pairs)
  class(out) <- "init_smd"
  out
}
init_ks <- function(x, treat, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE, ...) {
  chk::chk_not_missing(treat)
  chk::chk_atomic(treat)

  chk::chk_flag(pairwise)

  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("binary", "multinomial")) {
    .err("`treat` must be a binary or multi-category variable")
  }

  .chk_null_or(estimand, chk::chk_string)
  if (is_null(estimand)) estimand <- "ATE"

  f.e <- process.focal.and.estimand(focal, estimand, treat)
  focal <- f.e[["focal"]]
  estimand <- f.e[["estimand"]]

  unique.treats <- unique(treat)

  if (treat.type == "multinomial") {
    if (is_null(focal) && !pairwise) {
      treat.all <- last(make.unique(unique.treats, "All"))
      treat <- factor(c(as.character(treat), rep(treat.all, length(treat))),
                      levels = c(unique.treats, treat.all))
      x <- rbind(x, x)
      s.weights <- rep(s.weights, 2)
      focal <- treat.all
    }

    if (is_null(focal) || pairwise) {
      treatment.pairs <- combn(unique.treats, 2, simplify = FALSE)
    }
    else {
      treatment.pairs <- lapply(setdiff(unique.treats, focal), c, focal)
    }
  }
  else {
    treatment.pairs <- list(unique.treats)
    pairwise <- TRUE
  }

  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              focal = focal,
              pairwise = pairwise,
              treatment.pairs = treatment.pairs)
  class(out) <- "init_ks"
  out
}
init_ovl <- function(x, treat, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE,
                     integrate = FALSE, ...) {
  chk::chk_not_missing(treat)
  chk::chk_atomic(treat)

  chk::chk_flag(pairwise)

  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")

  chk::chk_flag(integrate)

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("binary", "multinomial")) {
    .err("`treat` must be a binary or multi-category variable")
  }

  .chk_null_or(estimand, chk::chk_string)
  if (is_null(estimand)) estimand <- "ATE"

  f.e <- process.focal.and.estimand(focal, estimand, treat)
  focal <- f.e[["focal"]]
  estimand <- f.e[["estimand"]]

  unique.treats <- unique(treat)

  if (treat.type == "multinomial") {
    if (is_null(focal) && !pairwise) {
      treat.all <- last(make.unique(unique.treats, "All"))
      treat <- factor(c(as.character(treat), rep(treat.all, length(treat))),
                      levels = c(unique.treats, treat.all))
      x <- rbind(x, x)
      s.weights <- rep(s.weights, 2)
      focal <- treat.all
    }

    if (is_null(focal) || pairwise) {
      treatment.pairs <- combn(unique.treats, 2, simplify = FALSE)
    }
    else {
      treatment.pairs <- lapply(setdiff(unique.treats, focal), c, focal)
    }
  }
  else {
    treatment.pairs <- list(unique.treats)
    pairwise <- TRUE
  }

  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              focal = focal,
              pairwise = pairwise,
              treatment.pairs = treatment.pairs,
              integrate = integrate)
  class(out) <- "init_ovl"
  out
}
init_ent <- function(x, treat, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE,
                     integrate = FALSE, ...) {
  chk::chk_not_missing(treat)
  chk::chk_atomic(treat)

  chk::chk_flag(pairwise)

  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")

  chk::chk_flag(integrate)

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("binary", "multinomial")) {
    .err("`treat` must be a binary or multi-category variable")
  }

  .chk_null_or(estimand, chk::chk_string)
  if (is_null(estimand)) estimand <- "ATE"

  f.e <- process.focal.and.estimand(focal, estimand, treat)
  focal <- f.e[["focal"]]
  estimand <- f.e[["estimand"]]

  unique.treats <- unique(treat)

  if (treat.type == "multinomial") {
    if (is_null(focal) && !pairwise) {
      treat.all <- last(make.unique(unique.treats, "All"))
      treat <- factor(c(as.character(treat), rep(treat.all, length(treat))),
                      levels = c(unique.treats, treat.all))
      x <- rbind(x, x)
      s.weights <- rep(s.weights, 2)
      focal <- treat.all
    }

    if (is_null(focal) || pairwise) {
      treatment.pairs <- combn(unique.treats, 2, simplify = FALSE)
    }
    else {
      treatment.pairs <- lapply(setdiff(unique.treats, focal), c, focal)
    }
  }
  else {
    treatment.pairs <- list(unique.treats)
    pairwise <- TRUE
  }

  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              focal = focal,
              pairwise = pairwise,
              treatment.pairs = treatment.pairs,
              integrate = integrate)
  class(out) <- "init_ent"
  out
}
init_mahalanobis <- function(x, treat, s.weights = NULL, estimand = NULL, focal = NULL, ...) {
  chk::chk_not_missing(treat)
  chk::chk_atomic(treat)

  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("binary", "multinomial")) {
    .err("`treat` must be a binary or multi-category variable")
  }

  .chk_null_or(estimand, chk::chk_string)
  if (is_null(estimand)) estimand <- "ATE"

  f.e <- process.focal.and.estimand(focal, estimand, treat)
  focal <- f.e[["focal"]]
  estimand <- f.e[["estimand"]]

  s.d.denom <- get.s.d.denom.weightit(estimand = estimand, treat = treat, focal = focal)

  if (any(!bin.vars)) x[,!bin.vars] <- scale(x[,!bin.vars])

  if (s.d.denom %in% as.character(treat)) {
    sigma <- cov.wt(x[treat == s.d.denom,,drop = FALSE], s.weights[treat == s.d.denom])$cov
    if (any(zeros <- diag(sigma) == 0)) {
      sigma_all <- cov.wt(x, s.weights)$cov
      sigma[zeros,] <- sigma_all[zeros,]
      sigma[,zeros] <- sigma_all[,zeros]
      if (any(zeros <- diag(sigma) == 0)) diag(sigma)[zeros]  <- 1
    }
  }
  else if (s.d.denom == "pooled") {
    sigma <- .5 * (cov.wt(x[treat == treat[1],,drop = FALSE], s.weights[treat == treat[1]])$cov +
                     cov.wt(x[treat != treat[1],,drop = FALSE], s.weights[treat != treat[1]])$cov)
    if (any(zeros <- diag(sigma) == 0)) diag(sigma)[zeros] <- 1
  }
  else {
    sigma <- cov.wt(x, s.weights)$cov
    if (any(zeros <- diag(sigma) == 0)) diag(sigma)[zeros] <- 1
  }

  #MASS::ginv
  sigmasvd <- svd(sigma)
  pos <- sigmasvd$d > max(1e-8 * sigmasvd$d[1L], 0)
  sigma_inv <- sigmasvd$v[, pos, drop = FALSE] %*% ((1/sigmasvd$d[pos]) *
                                                      t(sigmasvd$u[, pos, drop = FALSE]))

  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              s.d.denom = s.d.denom,
              sigma_inv = sigma_inv)
  class(out) <- "init_mahalanobis"
  out
}
init_energy.dist <- function(x, treat, s.weights = NULL, estimand = NULL, focal = NULL, improved = TRUE, ...) {
  chk::chk_not_missing(treat)
  chk::chk_atomic(treat)

  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("binary", "multinomial")) {
    .err("`treat` must be a binary or multi-category variable")
  }

  treat <- factor(treat)

  if (is_not_null(estimand)) {
    chk::chk_string(estimand)
    f.e <- process.focal.and.estimand(focal, estimand, treat)
    focal <- f.e[["focal"]]
    estimand <- f.e[["estimand"]]
  }

  dist.covs <- scale(x, scale = sqrt(col.w.v(x, s.weights, bin.vars)))

  d <- unname(as.matrix(dist(dist.covs)))

  n <- length(treat)
  unique.treats <- levels(treat)

  for (t in unique.treats) {
    s.weights[treat == t] <- s.weights[treat == t]/mean_fast(s.weights[treat == t])
  }

  treat_t <- vapply(unique.treats, function(t) treat == t, logical(n))
  n_t <- colSums(treat_t)

  s.weights_n_t <- setNames(lapply(unique.treats, function(t) treat_t[,t] * s.weights / n_t[t]),
                            unique.treats)

  if (is_null(estimand)) {
    all_pairs <- combn(unique.treats, 2, simplify = FALSE)
    P <- - d * Reduce("+", lapply(all_pairs, function(p) {
      tcrossprod(s.weights_n_t[[p[1]]] - s.weights_n_t[[p[2]]])
    }))
    q <- rep(0, n)
  }
  else if (is_null(focal)) {

    P <- -d * Reduce("+", lapply(s.weights_n_t, tcrossprod))

    q <- ((s.weights * 2 / n) %*% d) * Reduce("+", s.weights_n_t)

    if (improved) {
      all_pairs <- combn(unique.treats, 2, simplify = FALSE)
      P <- P - d * Reduce("+", lapply(all_pairs, function(p) {
        tcrossprod(s.weights_n_t[[p[1]]] - s.weights_n_t[[p[2]]])
      }))
    }
  }
  else {
    non_focal <- setdiff(unique.treats, focal)
    in_focal <- treat == focal

    P <- -d[!in_focal, !in_focal] *
      Reduce("+", lapply(s.weights_n_t[non_focal], function(s) tcrossprod(s[!in_focal])))

    q <- 2 * (s.weights_n_t[[focal]][in_focal] %*% d[in_focal, !in_focal]) *
      Reduce("+", lapply(s.weights_n_t[non_focal], function(s) s[!in_focal]))
  }

  out <- list(q = q,
              P = P,
              s.weights = s.weights,
              treat = treat,
              unique.treats = unique.treats,
              focal = focal)
  class(out) <- "init_energy.dist"
  out
}
init_kernel.dist <- function(x, treat, s.weights = NULL, estimand = NULL, focal = NULL, ...) {
  chk::chk_not_missing(treat)
  chk::chk_atomic(treat)

  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("binary")) {
    .err("`treat` must be a binary variable")
  }

  treat <- as.numeric(treat == treat[1])

  for (t in 0:1) {
    s.weights[treat == t] <- s.weights[treat == t]/mean_fast(s.weights[treat == t])
  }

  dist.covs <- scale(x, scale = sqrt(col.w.v(x, s.weights, bin.vars)))

  d <- unname(as.matrix(dist(dist.covs)))

  K <- exp(-(d^2)/median(d))

  T_star <- numeric(length(treat))
  T_star[treat == 1] <- 1/sum(treat == 1)
  T_star[treat == 0] <- -1/sum(treat == 0)

  out <- list(K = K,
              T_star = T_star,
              s.weights = s.weights,
              treat = treat)
  class(out) <- "init_kernel.dist"
  out
}
init_p <- function(x, treat, s.weights = NULL, ...) {
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")

  chk::chk_not_missing(treat)
  chk::chk_atomic(treat)

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("continuous")) {
    .err("`treat` must be a continuous (numeric) variable")
  }

  s.d.denom <- get.s.d.denom.cont.weightit()

  denoms <- compute_s.d.denom(x, treat = treat,
                              s.d.denom = s.d.denom, s.weights = s.weights,
                              bin.vars = bin.vars)

  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              s.d.denom = denoms)
  class(out) <- "init_p"
  out
}
init_s <- function(x, treat, s.weights = NULL, ...) {
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")

  chk::chk_not_missing(treat)
  chk::chk_atomic(treat)

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("continuous")) {
    .err("`treat` must be a continuous (numeric) variable")
  }

  for (i in seq_len(ncol(x))[!bin.vars[i]]) {
    x[,i] <- rank(x[,i], na.last = "keep")
  }
  treat <- rank(treat, na.last = "keep")

  s.d.denom <- get.s.d.denom.cont.weightit()

  denoms <- compute_s.d.denom(x, treat = treat,
                              s.d.denom = s.d.denom, s.weights = s.weights,
                              bin.vars = bin.vars)

  out <- list(treat = treat,
              covs = x,
              bin.vars = bin.vars,
              s.weights = s.weights,
              s.d.denom = denoms)
  class(out) <- "init_s"
  out
}
init_r2 <- function(x, treat, s.weights = NULL, poly = 1, int = FALSE, ...) {
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")

  chk::chk_not_missing(treat)
  chk::chk_atomic(treat)

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("binary", "continuous")) {
    .err("`treat` must be a binary or continuous (numeric) variable")
  }

  if (treat.type == "binary") treat <- as.numeric(treat == treat[1])

  x <- cbind(`(Intercept)` = 1, x, int.poly.f(x, poly = poly, int = int))

  out <- list(treat = treat,
              x = x,
              s.weights = s.weights)
  class(out) <- "init_r2"
  out
}
init_distance.cov <- function(x, treat, s.weights = NULL, ...) {
  x <- process_init_covs(x)
  bin.vars <- attr(x, "bin")

  chk::chk_not_missing(treat)
  chk::chk_atomic(treat)

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("continuous")) {
    .err("`treat` must be a continuous (numeric) variable")
  }

  dist.covs <- scale(x, scale = sqrt(col.w.v(x, s.weights, bin.vars)))

  Xdist <- unname(as.matrix(dist(dist.covs)))

  n <- length(treat)

  Adist <- unname(as.matrix(dist(treat/sqrt(col.w.v(treat, s.weights)))))

  Xmeans <- colMeans(Xdist)
  Xgrand_mean <- mean(Xmeans)
  XA <- Xdist + Xgrand_mean - outer(Xmeans, Xmeans, "+")

  Ameans <- colMeans(Adist)
  Agrand_mean <- mean(Ameans)
  AA <- Adist + Agrand_mean - outer(Ameans, Ameans, "+")

  P <- XA * AA/n^2

  out <- list(P = P,
              s.weights = s.weights,
              treat = treat)
  class(out) <- "init_distance.cov"
  out
}
init_l1.med <- function(x, treat, s.weights = NULL, estimand = NULL, focal = NULL,
                        .covs = NULL, l1.min.bin = 2, l1.max.bin = 12, l1.n = 101, ...) {

  if (is_not_null(.covs)) x <- .covs
  if (!is.data.frame(x)) {
    if (is.atomic(x) && is_null(dim(x))) x <- data.frame(x)
    else if (!is.matrix(x)) .err("`x` must be a data.frame or matrix.")
  }
  x <- as.data.frame(x)

  chk::chk_not_missing(treat)
  chk::chk_atomic(treat)

  check_arg_lengths(x, treat, s.weights)

  if (is_null(s.weights)) s.weights <- rep(1, NROW(x))

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (treat.type %nin% c("binary", "multinomial")) {
    .err("`treat` must be a binary or multi-category variable")
  }

  coarsen <- function(covs, cutpoints = NULL, grouping = NULL) {
    is.numeric.cov <- setNames(vapply(covs, is.numeric, logical(1L)), names(covs))
    for (i in names(cutpoints)) {
      if (cutpoints[[i]] == 0) is.numeric.cov[i] <- FALSE #Will not be binned
    }

    #Process grouping
    if (!is.null(grouping) && !is.null(names(grouping))) {
      covs[names(grouping)] <- lapply(names(grouping), function(g) {
        x <- covs[[g]]
        groups <- grouping[[g]]

        for (i in seq_along(groups)) {
          x[x %in% groups[[i]]] <- groups[[i]][1]
        }
        x
      })
      cutpoints[names(cutpoints) %in% names(grouping)] <- NULL
    }

    #Create bins for numeric variables
    for (i in names(covs)[is.numeric.cov]) {
      bins <- cutpoints[[i]]

      #cutpoints is number of bins, unlike in cem
      breaks <- seq(min(covs[[i]]), max(covs[[i]]), length = bins + 1)
      breaks[c(1, bins + 1)] <- c(-Inf, Inf)

      covs[[i]] <- findInterval(covs[[i]], breaks)
    }

    #Reduce to strata
    factor(do.call("paste", c(covs, sep = " | ")))
  }

  .chk_null_or(estimand, chk::chk_string)
  if (is_null(estimand)) estimand <- "ATE"

  f.e <- process.focal.and.estimand(focal, estimand, treat)
  focal <- f.e[["focal"]]
  estimand <- f.e[["estimand"]]

  unique.treats <- unique(treat)
  for (t in unique.treats) s.weights[treat == t] <- s.weights[treat == t]/sum(s.weights[treat == t])

  is.numeric.cov <- setNames(vapply(x, is.numeric, logical(1L)), names(x))
  nunique.covs <- vapply(x, nunique, integer(1L))

  coarsenings <- lapply(1:l1.n, function(i) {
    cutpoints <- setNames(lapply(nunique.covs[is.numeric.cov], function(nu) {
      sample(seq(min(l1.min.bin, nu), min(l1.max.bin, nu)), 1)
    }), names(x)[is.numeric.cov])
    grouping <- setNames(lapply(seq_along(x)[!is.numeric.cov], function(i) {
      nu <- nunique.covs[i]
      u <- unique(x[[i]], nmax = nu)

      #Randomly select number of bins
      nbins <- sample(seq(min(l1.min.bin, nu), min(l1.max.bin, nu)), 1)

      #Randomly assign bin numbers to levels of covariate
      bin.assignments <- sample(seq_len(nbins), nu, replace = TRUE)

      #Group levels with same bin number
      lapply(unique(bin.assignments, nmax = nbins),
             function(b) u[bin.assignments == b])

    }), names(x)[!is.numeric.cov])
    list(cutpoints = cutpoints, grouping = grouping, treat_cutpoints = NULL)
  })

  l1s <- unlist(lapply(coarsenings, function(co) {
    x <- coarsen(x, cutpoints = co[["cutpoints"]], grouping = co[["grouping"]])

    if (treat.type == "binary" || is_null(focal)) {
      sum(vapply(levels(x), function(l) {
        in_l <- which(x == l)
        abs(diff(range(vapply(unique.treats, function(t) sum(s.weights[in_l][treat[in_l] == t]), numeric(1L)))))
      }, numeric(1L))) / length(unique.treats)
    }
    else {
      sum(vapply(levels(x), function(l) {
        in_l <- which(x == l)
        sum.s.weights.focal <- sum(s.weights[in_l][treat[in_l] == focal])
        max(abs(vapply(unique.treats[unique.treats != focal],
                       function(t) sum(s.weights[in_l][treat[in_l] == t]) - sum.s.weights.focal, numeric(1L))))
      }, numeric(1L))) / length(unique.treats)
    }
  }))

  l1.med <- sort(l1s, partial = ceiling(l1.n/2))[ceiling(l1.n/2)]

  l1.med.coarsening <- coarsenings[[which(l1s == l1.med)[1]]]

  out <- list(coarsened.covs = coarsen(x,
                                       cutpoints = l1.med.coarsening[["cutpoints"]],
                                       grouping = l1.med.coarsening[["grouping"]]),
              s.weights = s.weights,
              treat = treat,
              unique.treats = unique.treats,
              focal = focal)
  class(out) <- "init_l1.med"
  out

}

#Statistics
smd.binary <- function(init, weights = NULL) {
  check_init(init, "init_smd")
  cobalt::col_w_smd(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
                    bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE)
}
ks.binary <- function(init, weights = NULL) {
  check_init(init, "init_ks")
  cobalt::col_w_ks(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
                   bin.vars = init$bin.vars)
}
ovl.binary <- function(init, weights = NULL) {
  check_init(init, "init_ovl")
  cobalt::col_w_ovl(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
                    bin.vars = init$bin.vars, integrate = init$integrate)
}
mahalanobis.binary <- function(init, weights = NULL) {
  check_init(init, "init_mahalanobis")
  mean.diffs <- cobalt::col_w_smd(init$covs, init$treat, weights, s.weights = init$s.weights,
                                  bin.vars = init$bin.vars, std = FALSE)
  drop(sqrt(t(mean.diffs) %*% init$sigma_inv %*% mean.diffs))
}
energy.dist.binary <- function(init, weights = NULL) {
  check_init(init, "init_energy.dist")

  if (is_null(weights)) weights <- init[["s.weights"]]
  else weights <- weights * init[["s.weights"]]

  for (t in init[["unique.treats"]]) {
    weights[init[["treat"]] == t] <- weights[init[["treat"]] == t]/mean_fast(weights[init[["treat"]] == t])
  }

  if (is_not_null(init[["focal"]])) {
    weights <- weights[init[["treat"]] != init[["focal"]]]
  }

  return(drop(weights %*% init[["P"]] %*% weights + init[["q"]] %*% weights))
}
kernel.dist.binary <- function(init, weights = NULL) {
  check_init(init, "init_kernel.dist")

  if (is_null(weights)) weights <- init[["s.weights"]]
  else weights <- weights * init[["s.weights"]]

  for (t in 0:1) {
    weights[init[["treat"]] == t] <- weights[init[["treat"]] == t]/mean_fast(weights[init[["treat"]] == t])
  }

  T_weights <- init[["T_star"]] * weights

  drop(sqrt(T_weights %*% init[["K"]] %*% T_weights))
}
r2.binary <- function(init, weights = NULL) {
  check_init(init, "init_r2")

  if (is_null(weights)) weights <- init[["s.weights"]]
  else weights <- weights * init[["s.weights"]]

  fit <- glm.fit(init$x, init$treat, weights, family = quasibinomial())

  wmtreat <- sum(weights*fit$linear.predictors)/sum(weights)

  SSmodel <- sum(weights * (fit$linear.predictors - wmtreat)^2)

  SSmodel / (sum(weights) * pi^2/3 + SSmodel)
}
l1.med.binary <- function(init, weights = NULL) {
  check_init(init, "init_l1.med")

  if (is_null(weights)) weights <- init[["s.weights"]]
  else weights <- weights * init[["s.weights"]]

  for (t in init[["unique.treats"]]) {
    weights[init[["treat"]] == t] <- weights[init[["treat"]] == t]/sum(weights[init[["treat"]] == t])
  }

  x <- init[["coarsened.covs"]]

  sum(vapply(levels(x), function(l) {
    in_l <- which(x == l)
    abs(diff(vapply(init[["unique.treats"]], function(t) sum(weights[in_l][init[["treat"]][in_l] == t]), numeric(1L))))
  }, numeric(1L))) / 2

}
smd.multinomial <- function(init, weights = NULL) {
  check_init(init, "init_smd")

  if (!init$pairwise) {
    weights <- c(weights, rep(1, length(weights)))
  }

  unlist(lapply(init$treatment.pairs, function(x) {
    cobalt::col_w_smd(init$covs[init$treat %in% x,,drop = FALSE],
                      treat = init$treat[init$treat %in% x],
                      weights = weights[init$treat %in% x],
                      s.weights = init$s.weights[init$treat %in% x],
                      bin.vars = init$bin.vars,
                      s.d.denom = init$s.d.denom, abs = TRUE)
  }))
}
ks.multinomial <- function(init, weights = NULL) {
  check_init(init, "init_ks")

  if (!init$pairwise) {
    weights <- c(weights, rep(1, length(weights)))
  }

  unlist(lapply(init$treatment.pairs, function(x) {
    cobalt::col_w_ks(init$covs[init$treat %in% x,,drop = FALSE],
                     treat = init$treat[init$treat %in% x],
                     weights = weights[init$treat %in% x],
                     s.weights = init$s.weights[init$treat %in% x],
                     bin.vars = init$bin.vars)
  }))
}
ovl.multinomial <- function(init, weights = NULL) {
  check_init(init, "init_ovl")

  if (!init$pairwise) {
    weights <- c(weights, rep(1, length(weights)))
  }

  unlist(lapply(init$treatment.pairs, function(x) {
    cobalt::col_w_ovl(init$covs[init$treat %in% x,,drop = FALSE],
                      treat = init$treat[init$treat %in% x],
                      weights = weights[init$treat %in% x],
                      s.weights = init$s.weights[init$treat %in% x],
                      bin.vars = init$bin.vars,
                      integrate = init$integrate)
  }))
}
energy.dist.multinomial <- function(init, weights = NULL) {
  energy.dist.binary(init, weights)
}
l1.med.multinomial <- function(init, weights = NULL) {
  check_init(init, "init_l1.med")

  if (is_null(weights)) weights <- init[["s.weights"]]
  else weights <- weights * init[["s.weights"]]

  for (t in init[["unique.treats"]]) {
    weights[init[["treat"]] == t] <- weights[init[["treat"]] == t] / sum(weights[init[["treat"]] == t])
  }

  x <- init[["coarsened.covs"]]

  if (is_null(init[["focal"]])) {
    sum(vapply(levels(x), function(l) {
      in_l <- which(x == l)
      abs(diff(range(vapply(init[["unique.treats"]], function(t) {
        sum(weights[in_l][init[["treat"]][in_l] == t])
      }, numeric(1L)))))
    }, numeric(1L))) / length(init[["unique.treats"]])
  }
  else {
    sum(vapply(levels(x), function(l) {
      in_l <- which(x == l)
      sum.weights.focal <- sum(weights[in_l][init[["treat"]][in_l] == init[["focal"]]])
      max(abs(vapply(init[["unique.treats"]][init[["unique.treats"]] != init[["focal"]]],
                     function(t) {
                       sum(weights[in_l][init[["treat"]][in_l] == t]) - sum.weights.focal
                     }, numeric(1L))))
    }, numeric(1L))) / length(init[["unique.treats"]])
  }

}
pearson.corr.continuous <- function(init, weights = NULL) {
  check_init(init, "init_p")
  cobalt::col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
                    bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
                    std = TRUE)
}
spearman.corr.continuous <- function(init, weights = NULL) {
  check_init(init, "init_s")
  cobalt::col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
                    bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
                    std = TRUE)
}
r2.continuous <- function(init, weights = NULL) {
  check_init(init, "init_r2")

  if (is_null(weights)) weights <- init[["s.weights"]]
  else weights <- weights * init[["s.weights"]]

  fit <- lm.wfit(init$x, init$treat, weights)

  SSresid <- sum(weights * (fit$residuals - w.m(fit$residuals, weights))^2)
  SStreat <- sum(weights * (init$treat - w.m(init$treat, weights))^2)

  1 - SSresid/SStreat
}
distance.cov.continuous <- function(init, weights = NULL) {
  check_init(init, "init_distance.cov")

  if (is_null(weights)) weights <- init[["s.weights"]]
  else weights <- weights * init[["s.weights"]]

  weights <- weights/mean_fast(weights)

  return(drop(t(weights) %*% init[["P"]] %*% weights))
}

bal_criterion <- function(treat.type, criterion) {
  chk::chk_not_missing(criterion)
  chk::chk_not_missing(treat.type)

  bal.obj <- switch(
    treat.type,
    binary = switch(criterion,
                    smd.mean = list(
                      fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_smd(covs, treat, s.weights, estimand)
                        }
                        mean_fast(smd.binary(init, weights))
                      },
                      init = init_smd
                    ),
                    smd.max = list(
                      fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_smd(covs, treat, s.weights, estimand)
                        }
                        max(smd.binary(init, weights))
                      },
                      init = init_smd
                    ),
                    smd.rms = list(
                      fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_smd(covs, treat, s.weights, estimand)
                        }
                        rms(smd.binary(init, weights))
                      },
                      init = init_smd
                    ),
                    ks.mean = list(
                      fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_ks(covs, treat, s.weights)
                        }
                        mean_fast(ks.binary(init, weights))
                      },
                      init = init_ks
                    ),
                    ks.max = list(
                      fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_ks(covs, treat, s.weights)
                        }
                        max(ks.binary(init, weights))
                      },
                      init = init_ks
                    ),
                    ks.rms = list(
                      fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_ks(covs, treat, s.weights)
                        }
                        rms(ks.binary(init, weights))
                      },
                      init = init_ks
                    ),
                    ovl.mean = list(
                      fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, integrate = FALSE,
                                     init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_ovl(covs, treat, s.weights, integrate = integrate)
                        }
                        mean_fast(ovl.binary(init, weights))
                      },
                      init = init_ovl
                    ),
                    ovl.max = list(
                      fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, integrate = FALSE,
                                     init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_ovl(covs, treat, s.weights, integrate = integrate)
                        }
                        max(ovl.binary(init, weights))
                      },
                      init = init_ovl
                    ),
                    ovl.rms = list(
                      fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, integrate = FALSE,
                                     init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_ovl(covs, treat, s.weights, integrate = integrate)
                        }
                        rms(ovl.binary(init, weights))
                      },
                      init = init_ovl
                    ),
                    mahalanobis = list(
                      fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_mahalanobis(covs, treat, s.weights, estimand)
                        }
                        mahalanobis.binary(init, weights)
                      },
                      init = init_mahalanobis
                    ),
                    energy.dist = list(
                      fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = NULL, focal = NULL, improved = TRUE, init = NULL, ...) {

                        if (is_null(init)) {
                          init <- init_energy.dist(covs, treat, s.weights, estimand, focal, improved)
                        }
                        energy.dist.binary(init, weights)
                      },
                      init = init_energy.dist
                    ),
                    kernel.dist = list(
                      fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = NULL, focal = NULL, init = NULL, ...) {

                        if (is_null(init)) {
                          init <- init_kernel.dist(covs, treat, s.weights, estimand, focal)
                        }
                        kernel.dist.binary(init, weights)
                      },
                      init = init_kernel.dist
                    ),
                    l1.med = list(
                      fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = NULL, focal = NULL, init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_l1.med(covs, treat, s.weights, estimand, focal, ...)
                        }
                        l1.med.binary(init, weights)
                      },
                      init = init_l1.med
                    ),
                    r2 = list(
                      fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_r2(covs, treat, s.weights)
                        }
                        r2.binary(init, weights)
                      },
                      init = init_r2
                    ),
                    r2.2 = list(
                      fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_r2(covs, treat, s.weights, poly = 2)
                        }
                        r2.binary(init, weights)
                      },
                      init = function(...) init_r2(..., poly = 2)
                    ),
                    r2.3 = list(
                      fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                        if (is_null(init)) {
                          init <- init_r2(covs, treat, s.weights, poly = 3)
                        }
                        r2.binary(init, weights)
                      },
                      init = function(...) init_r2(..., poly = 3)
                    )
    ),
    multinomial = switch(criterion,
                         smd.mean = list(
                           fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE, init = NULL, ...) {
                             if (is_null(init)) {
                               init <- init_smd(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                             }
                             mean_fast(smd.multinomial(init, weights))
                           },
                           init = init_smd
                         ),
                         smd.max = list(
                           fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE, init = NULL, ...) {
                             if (is_null(init)) {
                               init <- init_smd(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                             }
                             max(smd.multinomial(init, weights))
                           },
                           init = init_smd
                         ),
                         smd.rms = list(
                           fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = NULL, focal = NULL, pairwise = TRUE, init = NULL, ...) {
                             if (is_null(init)) {
                               init <- init_smd(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                             }
                             rms(smd.multinomial(init, weights))
                           },
                           init = init_smd
                         ),
                         ks.mean = list(
                           fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                                          pairwise = TRUE, init = NULL, ...) {
                             if (is_null(init)) {
                               init <- init_ks(covs, treat, s.weights, focal = focal, pairwise = pairwise)
                             }
                             mean_fast(ks.multinomial(init, weights))
                           },
                           init = init_ks
                         ),
                         ks.max = list(
                           fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                                          pairwise = TRUE, init = NULL, ...) {
                             if (is_null(init)) {
                               init <- init_ks(covs, treat, s.weights, focal = focal, pairwise = pairwise)
                             }
                             max(ks.multinomial(init, weights))
                           },
                           init = init_ks
                         ),
                         ks.rms = list(
                           fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                                          pairwise = TRUE, init = NULL, ...) {
                             if (is_null(init)) {
                               init <- init_ks(covs, treat, s.weights, focal = focal, pairwise = pairwise)
                             }
                             rms(ks.multinomial(init, weights))
                           },
                           init = init_ks
                         ),
                         ovl.mean = list(
                           fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                                          pairwise = TRUE, integrate = FALSE, init = NULL, ...) {
                             if (is_null(init)) {
                               init <- init_ovl(covs, treat, s.weights, focal = focal, pairwise = pairwise,
                                                integrate = integrate)
                             }
                             mean_fast(ovl.multinomial(init, weights))
                           },
                           init = init_ovl
                         ),
                         ovl.max = list(
                           fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                                          pairwise = TRUE, integrate = FALSE, init = NULL, ...) {
                             if (is_null(init)) {
                               init <- init_ovl(covs, treat, s.weights, focal = focal, pairwise = pairwise,
                                                integrate = integrate)
                             }
                             max(ovl.multinomial(init, weights))
                           },
                           init = init_ovl
                         ),
                         ovl.rms = list(
                           fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL,
                                          pairwise = TRUE, integrate = FALSE, init = NULL, ...) {
                             if (is_null(init)) {
                               init <- init_ovl(covs, treat, s.weights, focal = focal, pairwise = pairwise,
                                                integrate = integrate)
                             }
                             rms(ovl.multinomial(init, weights))
                           },
                           init = init_ovl
                         ),
                         energy.dist = list(
                           fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = NULL, focal = NULL, improved = TRUE, init = NULL, ...) {

                             if (is_null(init)) {
                               init <- init_energy.dist(covs, treat, s.weights, estimand, focal, improved)
                             }
                             energy.dist.multinomial(init, weights)
                           },
                           init = init_energy.dist
                         ),
                         l1.med = list(
                           fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = NULL, focal = NULL, init = NULL, ...) {
                             if (is_null(init)) {
                               init <- init_l1.med(covs, treat, s.weights, estimand, focal, ...)
                             }
                             l1.med.multinomial(init, weights)
                           },
                           init = init_l1.med
                         )
    ),

    continuous = switch(criterion,
                        p.mean = list(
                          fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                            if (is_null(init)) {
                              init <- init_p(covs, treat, s.weights)
                            }
                            mean_fast(pearson.corr.continuous(init, weights))
                          },
                          init = init_p
                        ),
                        p.max = list(
                          fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                            if (is_null(init)) {
                              init <- init_p(covs, treat, s.weights)
                            }
                            max(pearson.corr.continuous(init, weights))
                          },
                          init = init_p
                        ),
                        p.rms = list(
                          fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                            if (is_null(init)) {
                              init <- init_p(covs, treat, s.weights)
                            }
                            rms(pearson.corr.continuous(init, weights))
                          },
                          init = init_p
                        ),
                        s.mean = list(
                          fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                            if (is_null(init)) {
                              init <- init_s(covs, treat, s.weights)
                            }
                            mean_fast(spearman.corr.continuous(init, weights))
                          },
                          init = init_s
                        ),
                        s.max = list(
                          fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                            if (is_null(init)) {
                              init <- init_s(covs, treat, s.weights)
                            }
                            max(spearman.corr.continuous(init, weights))
                          },
                          init = init_s
                        ),
                        s.rms = list(
                          fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                            if (is_null(init)) {
                              init <- init_s(covs, treat, s.weights)
                            }
                            rms(spearman.corr.continuous(init, weights))
                          },
                          init = init_s
                        ),
                        r2 = list(
                          fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                            if (is_null(init)) {
                              init <- init_r2(covs, treat, s.weights)
                            }
                            r2.continuous(init, weights)
                          },
                          init = init_r2
                        ),
                        r2.2 = list(
                          fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                            if (is_null(init)) {
                              init <- init_r2(covs, treat, s.weights, poly = 2)
                            }
                            r2.continuous(init, weights)
                          },
                          init = function(...) init_r2(..., poly = 2)
                        ),
                        r2.3 = list(
                          fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                            if (is_null(init)) {
                              init <- init_r2(covs, treat, s.weights, poly = 3)
                            }
                            r2.continuous(init, weights)
                          },
                          init = function(...) init_r2(..., poly = 3)
                        ),
                        distance.cov = list(
                          fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                            if (is_null(init)) {
                              init <- init_distance.cov(covs, treat, s.weights)
                            }
                            distance.cov.continuous(init, weights)
                          },
                          init = init_distance.cov
                        )
    )
  )

  bal.obj
}

check_init <- function(init, init_class) {
  chk::chk_not_missing(init)
  chk::chk_not_missing(init_class)
  chk::chk_is(init, init_class)
}

check_arg_lengths <- function(...) {
  dots_names <- vapply(match.call(expand.dots = FALSE)$..., deparse1,
                       character(1L))
  lengths <- setNames(vapply(list(...), len, integer(1L)),
                      dots_names)
  supplied <- lengths > 0
  if (!all_the_same(lengths[supplied])) {
    .err(sprintf("%s must have the same number of units",
                 word_list(dots_names[supplied], quotes = "`")))
  }
}

compute_s.d.denom <- function(mat, treat, s.d.denom = "pooled", s.weights = NULL,
                              bin.vars = NULL, subset = NULL, weighted.weights = NULL,
                              to.sd = rep(TRUE, ncol(mat)), na.rm = TRUE) {
  denoms <- setNames(rep(1, ncol(mat)), colnames(mat))
  if (is.character(s.d.denom) && length(s.d.denom) == 1L) {
    if (is_null(bin.vars)) {
      bin.vars <- rep(FALSE, ncol(mat))
      bin.vars[to.sd] <- is_binary_col(mat[subset, to.sd,drop = FALSE])
    }
    else if (!is.atomic(bin.vars) || length(bin.vars) != ncol(mat) ||
             anyNA(as.logical(bin.vars))) {
      stop("'bin.vars' must be a logical vector with length equal to the number of columns of 'mat'.")
    }

    possibly.supplied <- c("mat", "treat", "weighted.weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
      stop(paste(word_list(possibly.supplied[supplied], quotes = 1), "must have the same number of units."))
    }

    if (lengths["weighted.weights"] == 0) weighted.weights <- rep(1, NROW(mat))
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("'subset' must be a logical vector.")

    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    cont.treat <- get.treat.type(treat) == "continuous"

    if (!cont.treat) {
      treat <- as.character(treat)
      unique.treats <- unique(treat)
    }
    else unique.treats <- NULL

    if (s.d.denom %in% unique.treats)
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(col.w.v(mat[treat == s.d.denom, , drop = FALSE],
                     w = s.weights[treat == s.d.denom],
                     bin.vars = bin.vars, na.rm = na.rm))
      }

    else if (s.d.denom == "pooled")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(Reduce("+", lapply(unique.treats,
                                function(t) col.w.v(mat[treat == t, , drop = FALSE],
                                                    w = s.weights[treat == t],
                                                    bin.vars = bin.vars, na.rm = na.rm))) / length(unique.treats))
      }
    else if (s.d.denom == "all")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(col.w.v(mat, w = s.weights, bin.vars = bin.vars, na.rm = na.rm))
      }
    else if (s.d.denom == "weighted")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(col.w.v(mat, w = weighted.weights * s.weights, bin.vars = bin.vars, na.rm = na.rm))
      }
    else if (s.d.denom == "hedges")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        (1 - 3/(4*length(treat) - 9))^-1 * sqrt(Reduce("+", lapply(unique.treats,
                                                                   function(t) (sum(treat == t) - 1) * col.w.v(mat[treat == t, , drop = FALSE],
                                                                                                               w = s.weights[treat == t],
                                                                                                               bin.vars = bin.vars, na.rm = na.rm))) / (length(treat) - 2))
      }
    else stop("s.d.denom is not an allowed value.")

    denoms[to.sd] <- denom.fun(mat = mat[, to.sd, drop = FALSE], treat = treat, s.weights = s.weights,
                               weighted.weights = weighted.weights, bin.vars = bin.vars[to.sd],
                               unique.treats = unique.treats, na.rm = na.rm)

    if (any(zero_sds <- check_if_zero(denoms[to.sd]))) {
      denoms[to.sd][zero_sds] <- sqrt(col.w.v(mat[, to.sd, drop = FALSE][, zero_sds, drop = FALSE],
                                              w = s.weights,
                                              bin.vars = bin.vars[to.sd][zero_sds], na.rm = na.rm))
    }

    if (cont.treat) {
      treat.sd <- denom.fun(mat = treat, s.weights = s.weights,
                            weighted.weights = weighted.weights, bin.vars = FALSE,
                            na.rm = na.rm)
      denoms[to.sd] <- denoms[to.sd]*treat.sd
    }
  }
  else {
    if (is.numeric(s.d.denom)) {
      if (is_not_null(names(s.d.denom)) && any(colnames(mat) %in% names(s.d.denom))) {
        denoms[colnames(mat)[colnames(mat) %in% names(s.d.denom)]] <- s.d.denom[names(s.d.denom)[names(s.d.denom) %in% colnames(mat)]]
      }
      else if (length(s.d.denom) == sum(to.sd)) {
        denoms[to.sd] <- s.d.denom
      }
      else if (length(s.d.denom) == ncol(mat)) {
        denoms[] <- s.d.denom
      }
      else {
        stop("'s.d.denom' must be an allowable value or a numeric vector of with length equal to the number of columns of 'mat'. See ?cobalt::col_w_smd for allowable values.")
      }
    }
    else {
      stop("'s.d.denom' must be an allowable value or a numeric vector of with length equal to the number of columns of 'mat'. See ?cobalt::col_w_smd for allowable values.")
    }
  }
  return(denoms)
}