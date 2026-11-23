# Calibration of the families for continuous treatments with method = "glm"
# (see _dev/glm-continuous-families-notes.md).
#
# For four treatment distributions, the treatment depends on three covariates and
# the outcome is linear in the treatment and the covariates, so the slope of the
# marginal structural model E[Y(a)] = b0 + b1 a is known (b1 = 0.5). Each
# distribution is fit with the default Gaussian model and with the families meant
# for it, and the weighted regression of Y on A is fit with glm_weightit(), whose
# standard errors use M-estimation (HC0 standard errors, which treat the weights as
# fixed, are recorded too). The same regression with oracle weights, computed from
# the true conditional density and its average over the sample, separates the
# behavior of inverse probability weighting in each design from that of the
# estimated models.
#
# Settings: confounding strength k (the coefficients of the treatment model, and
# the slope of the log SD for the location-scale design, are multiplied by k) and
# sample size n: (k = 1, n = 1000), (k = .3, n = 1000), and (k = .3, n = 4000).
# With k = 1 the weights are heavy-tailed, and even the oracle weights give a
# biased slope with poor coverage, so the implementation is judged against the
# oracle there and against the truth at k = .3.
#
# Reported: bias and RMSE of the slope, coverage of the 95% Wald interval, ESS, and
# the largest weighted treatment-covariate correlation.
#
# Run from the package root with the development version installed in the
# library given by WEIGHTIT_LIB:
#   WEIGHTIT_LIB=<lib> Rscript _dev/glm-cont-families-sim.R
# Results: _dev/glm-cont-families-sim.rds (checkpointed while running). The stored
# results come from a version that also had bias-reduced fits (`br.` links), since
# removed; those rows can be ignored.

source("~/.claude/skills/live-progress/assets/progress.R")

lib <- Sys.getenv("WEIGHTIT_LIB")
reps <- as.integer(Sys.getenv("SIM_REPS", "200"))
workers <- as.integer(Sys.getenv("SIM_WORKERS", "3"))
out_file <- "_dev/glm-cont-families-sim.rds"

set.seed(20261120)

tau <- 0.5

dgps <- c("normal_ls", "poisson", "negbin", "gamma")

methods <- list(
  normal_ls = list(gaussian = list(family = "gaussian"),
                   gaussian_ls = list(family = "gaussian_ls")),
  poisson = list(gaussian = list(family = "gaussian"),
                 poisson = list(family = "poisson"),
                 negbin = list(family = "negbin")),
  negbin = list(gaussian = list(family = "gaussian"),
                poisson = list(family = "poisson"),
                negbin = list(family = "negbin")),
  gamma = list(gaussian = list(family = "gaussian"),
               gaussian_ls = list(family = "gaussian_ls"),
               Gamma = list(family = "Gamma"))
)

#Returns the data and the true conditional density of the treatment, as a
#function of a treatment value giving the density for every unit
sim_data <- function(dgp, n, k) {
  X1 <- rnorm(n)
  X2 <- rbinom(n, 1, .5)
  X3 <- runif(n, -1, 1)
  lp <- k * (.5 * X1 - .4 * X2 + .4 * X3)

  switch(dgp,
         normal_ls = {
           sd <- exp(-.2 + k * .6 * X1)
           A <- lp + sd * rnorm(n)
           f <- function(a) dnorm(a, lp, sd)
         },
         poisson = {
           mu <- exp(1 + lp)
           A <- rpois(n, mu)
           f <- function(a) dpois(a, mu)
         },
         negbin = {
           mu <- exp(1 + lp)
           A <- rnbinom(n, size = 1.5, mu = mu)
           f <- function(a) dnbinom(a, size = 1.5, mu = mu)
         },
         gamma = {
           mu <- exp(.5 + lp)
           A <- rgamma(n, shape = 2, rate = 2 / mu)
           f <- function(a) dgamma(a, shape = 2, rate = 2 / mu)
         })

  Y <- 1 + tau * A + X1 - X2 + .5 * X3 + rnorm(n)

  list(data = data.frame(A, Y, X1, X2, X3), f = f)
}

#Oracle weights: the true conditional density, with its average over the sample
#as the numerator
oracle_weights <- function(A, f) {
  ua <- unique(A)
  num <- vapply(ua, function(a) mean(f(a)), numeric(1L))[match(A, ua)]

  num / vapply(seq_along(A), function(i) f(A[i])[i], numeric(1L))
}

w_cor <- function(x, y, w) {
  w <- w / sum(w)
  mx <- sum(w * x)
  my <- sum(w * y)
  sum(w * (x - mx) * (y - my)) / sqrt(sum(w * (x - mx)^2) * sum(w * (y - my)^2))
}

settings <- data.frame(k = c(1, .3, .3), n = c(1000L, 1000L, 4000L))

cells <- merge(expand.grid(rep = seq_len(reps), dgp = dgps, setting = seq_len(nrow(settings)),
                           stringsAsFactors = FALSE),
               cbind(setting = seq_len(nrow(settings)), settings))
cells <- cells[order(cells$setting, cells$dgp, cells$rep), ]
rownames(cells) <- NULL

one_cell <- function(i) {
  suppressPackageStartupMessages(library(WeightIt, lib.loc = lib))

  cell <- cells[i, ]
  sim <- sim_data(cell$dgp, cell$n, cell$k)
  d <- sim$data

  summarize_w <- function(w, fit) {
    data.frame(est = unname(coef(fit)["A"]),
               se = unname(sqrt(diag(vcov(fit)))["A"]),
               se_hc0 = unname(sqrt(diag(vcov(fit, vcov = "HC0")))["A"]),
               ess = WeightIt::ESS(w),
               max_cor = max(abs(vapply(d[c("X1", "X2", "X3")], w_cor,
                                        numeric(1), y = d$A, w = w))),
               max_w = max(w / mean(w)),
               error = NA_character_)
  }

  failed <- function(e) {
    data.frame(est = NA_real_, se = NA_real_, se_hc0 = NA_real_, ess = NA_real_,
               max_cor = NA_real_, max_w = NA_real_, error = conditionMessage(e))
  }

  rows <- lapply(names(methods[[cell$dgp]]), function(m) {
    args <- methods[[cell$dgp]][[m]]

    res <- tryCatch({
      W <- suppressWarnings(do.call(WeightIt::weightit,
                                    c(list(A ~ X1 + X2 + X3, data = d, method = "glm"),
                                      args)))
      fit <- WeightIt::glm_weightit(Y ~ A, data = d, weightit = W)

      summarize_w(W$weights, fit)
    }, error = failed)

    cbind(dgp = cell$dgp, k = cell$k, n = cell$n, rep = cell$rep, method = m, res)
  })

  #Oracle weights; no weights are estimated, so the HC0 standard error is the
  #appropriate one and is also recorded as `se`
  w_or <- oracle_weights(d$A, sim$f)
  res_or <- tryCatch({
    fit <- WeightIt::glm_weightit(Y ~ A, data = d, weights = w_or, vcov = "HC0")
    summarize_w(w_or, fit)
  }, error = failed)

  rows <- c(rows, list(cbind(dgp = cell$dgp, k = cell$k, n = cell$n, rep = cell$rep,
                             method = "oracle", res_or)))

  do.call(rbind, rows)
}

pr <- prog_init(total = nrow(cells), title = "glm families: calibration",
                unit = "replicate", kind = "simulation", workers = workers)
on.exit(prog_end(pr, "failed"), add = TRUE)

res <- prog_parallel(pr, seq_len(nrow(cells)), one_cell, workers = workers,
                     label = function(i) sprintf("%s k=%.1f n=%d #%d", cells$dgp[i], cells$k[i], cells$n[i], cells$rep[i]),
                     checkpoint = "_dev/glm-cont-families-sim-parts/sim")

on.exit()
prog_end(pr, "done")

res <- do.call(rbind, res)

summ <- do.call(rbind, lapply(split(res, list(res$dgp, res$k, res$n, res$method), drop = TRUE), function(r) {
  ok <- !is.na(r$est)
  data.frame(dgp = r$dgp[1], k = r$k[1], n = r$n[1], method = r$method[1],
             fits = sum(ok), errors = sum(!ok),
             bias = mean(r$est[ok]) - tau,
             mcse_bias = sd(r$est[ok]) / sqrt(sum(ok)),
             rmse = sqrt(mean((r$est[ok] - tau)^2)),
             sd_est = sd(r$est[ok]),
             mean_se = mean(r$se[ok]),
             mean_se_hc0 = mean(r$se_hc0[ok]),
             coverage = mean(abs(r$est[ok] - tau) <= qnorm(.975) * r$se[ok]),
             coverage_hc0 = mean(abs(r$est[ok] - tau) <= qnorm(.975) * r$se_hc0[ok]),
             ess = mean(r$ess[ok]),
             max_cor = mean(r$max_cor[ok]),
             max_w = median(r$max_w[ok]))
}))

summ <- summ[order(summ$dgp, -summ$k, summ$n, summ$method), ]
rownames(summ) <- NULL

saveRDS(list(results = res, summary = summ, tau = tau, reps = reps, settings = settings,
             complete = TRUE), out_file)

print(summ, digits = 3)
