# Stress test of the numerator f(A) for the families for continuous treatments
# with method = "glm" (see _dev/glm-continuous-families-notes.md and
# _dev/glm-cont-families-numerator.R, whose designs turned out too mild: with
# weak confounding, the marginal distribution of the treatment is nearly a
# member of the same family, so an intercept-only fit of that family is nearly
# right).
#
# Here a binary covariate X2 shifts the mean of the treatment (or, for the
# location-scale design, its SD) by `sep` on the log scale, so the marginal
# distribution of the treatment is a two-component mixture: bimodal for the
# Poisson, negative binomial, and Gamma designs, and a heavy-tailed scale mixture
# for the location-scale design. The treatment model given X1, X2, and X3 is
# correctly specified in every design.
#
# The denominator is always the conditional density estimated by the correct
# family. The numerators compared:
#
#   avg_cond:   the sampling-weighted average of the estimated conditional
#               densities (what the families use)
#   marginal:   the same family fit with no covariates, by maximum likelihood
#   kde:        a kernel density estimate of the treatment
#   oracle:     not a numerator alone; the oracle weights, from the true
#               conditional density with the sample average of the true
#               conditional densities as the numerator
#
# Two outcomes are computed on each data set:
#
#   linear:     Y = 1 + .5 A + X1 - X2 + .5 X3 + e. The linear marginal
#               structural model is correct, so every numerator gives a
#               consistent slope (.5); the numerator affects precision only.
#   nonlinear:  Y = 1 + g(A) + X1 - X2 + .5 X3 + e, with g(a) = 2 sqrt(a) for
#               the count and Gamma designs and g(a) = a + .2 a^2 for the
#               location-scale design, still fit with a linear model. The
#               target is the slope of the linear projection of E[Y(a)] onto a
#               under the true marginal distribution of the treatment, computed
#               by Monte Carlo on 2 million draws. A numerator that does not
#               estimate the true marginal targets a different projection.
#
# Reported: bias and RMSE of the slope for each outcome, ESS, the largest
# weighted treatment-covariate correlation, and the median absolute log ratio
# of each numerator to the oracle numerator (how far it is from the true
# marginal density at the observed treatment values).
#
# Run from the package root with the development version installed in the
# library given by WEIGHTIT_LIB:
#   WEIGHTIT_LIB=<lib> Rscript _dev/glm-cont-families-numerator-stress.R
# Results: _dev/glm-cont-families-numerator-stress.rds.

source("~/.claude/skills/live-progress/assets/progress.R")

lib <- Sys.getenv("WEIGHTIT_LIB")
reps <- as.integer(Sys.getenv("SIM_REPS", "300"))
workers <- as.integer(Sys.getenv("SIM_WORKERS", "2"))
n <- as.integer(Sys.getenv("SIM_N", "1000"))
out_file <- "_dev/glm-cont-families-numerator-stress.rds"

set.seed(20261122)

dgps <- c(normal_ls = "gaussian_ls", poisson = "poisson", negbin = "negbin",
          gamma = "Gamma")

seps <- c(moderate = 1, strong = 1.8)

#Covariates, treatment, and the true conditional density (a function of a
#treatment value returning the density for every unit)
sim_treat <- function(dgp, n, sep) {
  X1 <- rnorm(n)
  X2 <- rbinom(n, 1, .5)
  X3 <- runif(n, -1, 1)
  lp <- .3 * X1 + .3 * X3

  switch(dgp,
         normal_ls = {
           mu <- lp + .5 * X2
           sd <- exp(-.3 + sep * (X2 - .5))
           A <- mu + sd * rnorm(n)
           f <- function(a) dnorm(a, mu, sd)
         },
         poisson = {
           mu <- exp(.7 + sep * X2 + lp)
           A <- rpois(n, mu)
           f <- function(a) dpois(a, mu)
         },
         negbin = {
           mu <- exp(.7 + sep * X2 + lp)
           A <- rnbinom(n, size = 3, mu = mu)
           f <- function(a) dnbinom(a, size = 3, mu = mu)
         },
         gamma = {
           mu <- exp(.3 + sep * X2 + lp)
           A <- rgamma(n, shape = 4, rate = 4 / mu)
           f <- function(a) dgamma(a, shape = 4, rate = 4 / mu)
         })

  list(X1 = X1, X2 = X2, X3 = X3, A = A, f = f)
}

g_nonlinear <- function(dgp, a) {
  if (dgp == "normal_ls") a + .2 * a^2 else 2 * sqrt(a)
}

#Slope of the linear projection of g(A) on A under the true marginal
#distribution of the treatment, by Monte Carlo
truth <- expand.grid(dgp = names(dgps), sep = names(seps), stringsAsFactors = FALSE)
truth$nonlinear <- mapply(function(dgp, sep) {
  s <- sim_treat(dgp, 2e6, seps[[sep]])
  unname(lm.fit(cbind(1, s$A), g_nonlinear(dgp, s$A))$coefficients[2L])
}, truth$dgp, truth$sep)
truth$linear <- .5

w_cor <- function(x, y, w) {
  w <- w / sum(w)
  mx <- sum(w * x)
  my <- sum(w * y)
  sum(w * (x - mx) * (y - my)) / sqrt(sum(w * (x - mx)^2) * sum(w * (y - my)^2))
}

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

#Sample average of the conditional densities at each treatment value
avg_density <- function(A, f) {
  ua <- unique(A)
  vapply(ua, function(a) mean(f(a)), numeric(1L))[match(A, ua)]
}

cells <- expand.grid(rep = seq_len(reps), dgp = names(dgps), sep = names(seps),
                     stringsAsFactors = FALSE)

one_cell <- function(i) {
  suppressPackageStartupMessages(library(WeightIt, lib.loc = lib))

  cell <- cells[i, ]
  family <- dgps[[cell$dgp]]
  s <- sim_treat(cell$dgp, n, seps[[cell$sep]])
  A <- s$A

  base <- 1 + s$X1 - s$X2 + .5 * s$X3 + rnorm(n)
  d <- data.frame(A, X1 = s$X1, X2 = s$X2, X3 = s$X3,
                  Y_lin = base + .5 * A,
                  Y_nl = base + g_nonlinear(cell$dgp, A))

  W <- suppressWarnings(WeightIt::weightit(A ~ X1 + X2 + X3, data = d,
                                           method = "glm", family = family,
                                           include.obj = TRUE))

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

  oracle_num <- avg_density(A, s$f)
  oracle_denom <- vapply(seq_along(A), function(i) s$f(A[i])[i], numeric(1L))

  kde <- density(A)

  weights <- list(avg_cond = W$weights,
                  marginal = marginal_density(family, A) / denom,
                  kde = approx(kde$x, kde$y, xout = A)$y / denom,
                  oracle = oracle_num / oracle_denom)

  nums <- list(avg_cond = W$weights * denom,
               marginal = marginal_density(family, A),
               kde = approx(kde$x, kde$y, xout = A)$y,
               oracle = oracle_num)

  tr <- truth[truth$dgp == cell$dgp & truth$sep == cell$sep, ]

  rows <- lapply(names(weights), function(nm) {
    w <- weights[[nm]]

    est_lin <- unname(coef(lm(Y_lin ~ A, data = d, weights = w))["A"])
    est_nl <- unname(coef(lm(Y_nl ~ A, data = d, weights = w))["A"])

    data.frame(dgp = cell$dgp, sep = cell$sep, rep = cell$rep, numerator = nm,
               err_lin = est_lin - tr$linear,
               err_nl = est_nl - tr$nonlinear,
               ess = WeightIt::ESS(w),
               max_cor = max(abs(vapply(d[c("X1", "X2", "X3")], w_cor,
                                        numeric(1), y = A, w = w))),
               max_w = max(w / mean(w)),
               num_misfit = median(abs(log(nums[[nm]] / oracle_num))))
  })

  do.call(rbind, rows)
}

pr <- prog_init(total = nrow(cells), title = "glm families: numerator stress test",
                unit = "replicate", kind = "simulation", workers = workers)
on.exit(prog_end(pr, "failed"), add = TRUE)

res <- prog_parallel(pr, seq_len(nrow(cells)), one_cell, workers = workers,
                     label = function(i) sprintf("%s %s #%d", cells$dgp[i], cells$sep[i], cells$rep[i]),
                     checkpoint = "_dev/glm-cont-families-numerator-stress-sim-parts/sim")

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, res)

summ <- do.call(rbind, lapply(split(res, list(res$dgp, res$sep, res$numerator), drop = TRUE), function(r) {
  data.frame(dgp = r$dgp[1], sep = r$sep[1], numerator = r$numerator[1],
             bias_lin = mean(r$err_lin),
             rmse_lin = sqrt(mean(r$err_lin^2)),
             bias_nl = mean(r$err_nl),
             mcse_nl = sd(r$err_nl) / sqrt(nrow(r)),
             rmse_nl = sqrt(mean(r$err_nl^2)),
             ess = mean(r$ess),
             max_cor = mean(r$max_cor),
             max_w = median(r$max_w),
             num_misfit = median(r$num_misfit))
}))

summ <- summ[order(summ$dgp, summ$sep, summ$numerator), ]
rownames(summ) <- NULL

saveRDS(list(results = res, summary = summ, truth = truth, seps = seps, n = n,
             reps = reps, complete = TRUE), out_file)

print(truth)
print(summ, digits = 3)
