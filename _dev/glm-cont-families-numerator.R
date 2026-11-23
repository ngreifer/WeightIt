# Comparison of numerators f(A) for the weights f(A) / f(A | X) with the families
# for continuous treatments with method = "glm" (see
# _dev/glm-continuous-families-notes.md). The denominator is always the one the
# correct family estimates; only the numerator changes:
#
#   avg_cond:   the sampling-weighted average over units of the estimated
#               conditional densities (what the families use)
#   kde:        a kernel density estimate of the treatment, `density()` with its
#               default bandwidth (what `density = "kernel"` would suggest);
#               for count treatments, it is still computed, as a continuous KDE
#               evaluated at the counts
#   freq:       for count treatments, the empirical proportion of units with
#               each treatment value
#   marginal:   the same family fit with no covariates (the classic stabilized
#               numerator of Robins et al., 2000), by maximum likelihood
#
# Any positive function of A gives weights under which A is independent of X in
# the weighted sample, so the numerator does not affect consistency when the
# marginal structural model is correctly specified; it affects the variability
# of the weights and so the precision of the estimate. Data-generating processes
# as in _dev/glm-cont-families-sim.R, with n = 1000 and k = .3 and 1.
#
# Run from the package root with the development version installed in the
# library given by WEIGHTIT_LIB:
#   WEIGHTIT_LIB=<lib> Rscript _dev/glm-cont-families-numerator.R
# Results: _dev/glm-cont-families-numerator.rds.

source("~/.claude/skills/live-progress/assets/progress.R")

lib <- Sys.getenv("WEIGHTIT_LIB")
reps <- as.integer(Sys.getenv("SIM_REPS", "200"))
workers <- as.integer(Sys.getenv("SIM_WORKERS", "2"))
out_file <- "_dev/glm-cont-families-numerator.rds"

set.seed(20261121)

tau <- 0.5

dgps <- c(normal_ls = "gaussian_ls", poisson = "poisson", negbin = "negbin",
          gamma = "Gamma")

sim_data <- function(dgp, n, k) {
  X1 <- rnorm(n)
  X2 <- rbinom(n, 1, .5)
  X3 <- runif(n, -1, 1)
  lp <- k * (.5 * X1 - .4 * X2 + .4 * X3)

  A <- switch(dgp,
              normal_ls = lp + exp(-.2 + k * .6 * X1) * rnorm(n),
              poisson = rpois(n, exp(1 + lp)),
              negbin = rnbinom(n, size = 1.5, mu = exp(1 + lp)),
              gamma = rgamma(n, shape = 2, rate = 2 / exp(.5 + lp)))

  Y <- 1 + tau * A + X1 - X2 + .5 * X3 + rnorm(n)

  data.frame(A, Y, X1, X2, X3)
}

w_cor <- function(x, y, w) {
  w <- w / sum(w)
  mx <- sum(w * x)
  my <- sum(w * y)
  sum(w * (x - mx) * (y - my)) / sqrt(sum(w * (x - mx)^2) * sum(w * (y - my)^2))
}

#Marginal (intercept-only) density of the same family, by maximum likelihood
marginal_density <- function(family, A) {
  switch(family,
         gaussian_ls = dnorm(A, mean(A), sqrt(mean((A - mean(A))^2))),
         poisson = dpois(A, mean(A)),
         negbin = {
           theta <- MASS::theta.ml(A, rep(mean(A), length(A)), limit = 100)
           dnbinom(A, size = theta, mu = mean(A))
         },
         Gamma = {
           k <- WeightIt:::.gamma_shape_ml(A, rep(mean(A), length(A)), rep(1, length(A)))
           dgamma(A, shape = k, rate = k / mean(A))
         })
}

cells <- expand.grid(rep = seq_len(reps), dgp = names(dgps), k = c(1, .3),
                     stringsAsFactors = FALSE)
cells$n <- 1000L

one_cell <- function(i) {
  suppressPackageStartupMessages(library(WeightIt, lib.loc = lib))

  cell <- cells[i, ]
  family <- dgps[[cell$dgp]]
  d <- sim_data(cell$dgp, cell$n, cell$k)
  A <- d$A

  W <- suppressWarnings(WeightIt::weightit(A ~ X1 + X2 + X3, data = d,
                                           method = "glm", family = family,
                                           include.obj = TRUE))

  #Conditional density at each unit's treatment (the denominator) and the
  #current numerator, recovered from the weights
  M <- attr(W, "Mparts")

  denom <- {
    if (family == "gaussian_ls") dnorm(A, W$obj$fitted.values, W$obj$sd)
    else {
      disp <- if (family %in% c("negbin", "Gamma")) exp(M$btreat[[1L]]) else NULL
      B <- if (is.null(disp)) M$btreat else M$btreat[-1L]
      mu <- exp(drop(M$Xtreat %*% B))
      switch(family,
             poisson = dpois(A, mu),
             negbin = dnbinom(A, size = disp, mu = mu),
             Gamma = dgamma(A, shape = disp, rate = disp / mu))
    }
  }

  nums <- list(avg_cond = W$weights * denom)

  kde <- density(A)
  nums$kde <- approx(kde$x, kde$y, xout = A)$y

  if (family %in% c("poisson", "negbin")) {
    nums$freq <- as.numeric(table(A)[as.character(A)]) / length(A)
  }

  nums$marginal <- marginal_density(family, A)

  rows <- lapply(names(nums), function(nm) {
    w <- nums[[nm]] / denom

    fit <- lm(Y ~ A, data = d, weights = w)

    data.frame(dgp = cell$dgp, k = cell$k, n = cell$n, rep = cell$rep,
               numerator = nm,
               est = unname(coef(fit)["A"]),
               ess = WeightIt::ESS(w),
               max_cor = max(abs(vapply(d[c("X1", "X2", "X3")], w_cor,
                                        numeric(1), y = A, w = w))),
               max_w = max(w / mean(w)))
  })

  do.call(rbind, rows)
}

pr <- prog_init(total = nrow(cells), title = "glm families: numerators",
                unit = "replicate", kind = "simulation", workers = workers)
on.exit(prog_end(pr, "failed"), add = TRUE)

res <- prog_parallel(pr, seq_len(nrow(cells)), one_cell, workers = workers,
                     label = function(i) sprintf("%s k=%.1f #%d", cells$dgp[i], cells$k[i], cells$rep[i]),
                     checkpoint = "_dev/glm-cont-families-numerator-sim-parts/sim")

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, res)

summ <- do.call(rbind, lapply(split(res, list(res$dgp, res$k, res$numerator), drop = TRUE), function(r) {
  data.frame(dgp = r$dgp[1], k = r$k[1], numerator = r$numerator[1],
             bias = mean(r$est) - tau,
             mcse_bias = sd(r$est) / sqrt(nrow(r)),
             sd_est = sd(r$est),
             rmse = sqrt(mean((r$est - tau)^2)),
             ess = mean(r$ess),
             max_cor = mean(r$max_cor),
             max_w = median(r$max_w))
}))

summ <- summ[order(summ$dgp, -summ$k, summ$numerator), ]
rownames(summ) <- NULL

saveRDS(list(results = res, summary = summ, tau = tau, reps = reps,
             complete = TRUE), out_file)

print(summ, digits = 3)
