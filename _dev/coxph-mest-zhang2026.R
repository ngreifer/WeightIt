# How much did the Cox M-estimation bug fixed in 2.2.0 affect the standard errors
# reported by Zhang, Chen, Du, and Wang (2026, arXiv:2608.06774), Table 2 (Type I
# error and variance estimation for the marginal HR, n = 400)?
#
# The bug: `.compute_vcov()` formed the cross-derivative block of the M-estimation
# Jacobian from Cox score residuals evaluated with every estimated weight set to 1,
# i.e., with unweighted risk sets, rather than at the weights used in the fit.
#
# The same data sets are analyzed by CRAN WeightIt 2.1.0 (installed library) and
# by the development version (2.2.0, via `pkgload::load_all()`), so the two sets
# of standard errors are paired. Run once per version:
#
#   script -q /dev/null Rscript _dev/coxph-mest-zhang2026.R cran 2000 < /dev/null
#   script -q /dev/null Rscript _dev/coxph-mest-zhang2026.R dev 2000 < /dev/null
#
# Data-generating process from the paper's Section 5.1 and Supplement (Data
# Generation Process Details), with beta_A = 0 and beta_Z = 9 as in the
# inference study. EB is `method = "ebal"`, IPW is `method = "glm"` (logistic
# regression), both with `estimand = "ATE"`, which balance each arm to the pooled
# sample as the paper's calibration constraints do.

args <- commandArgs(trailingOnly = TRUE)
version <- args[1L]
n_sim <- as.integer(args[2L])
stopifnot(version %in% c("cran", "dev", "v1709"), !is.na(n_sim))

root <- "/Users/NoahGreifer/Dropbox/Research/R/WeightIt"
source("/Users/NoahGreifer/.claude/skills/live-progress/assets/progress.R")

gen_data <- function(n = 400, beta_A = 0, beta_Z = 9) {
  D <- diag(c(.3, .25, .1))
  R <- matrix(c(1, .245, .239,
                .245, 1, -.036,
                .239, -.036, 1), 3L)
  L <- MASS::mvrnorm(n, mu = c(4, 3.5, 4.8), Sigma = D %*% R %*% D)

  d <- data.frame(LAge = L[, 1L], LHDL = L[, 2L], LSBP = L[, 3L],
                  DM = rbinom(n, 1, .2), Smoke = rbinom(n, 1, .15),
                  N1 = rnorm(n), N2 = rnorm(n), N3 = rnorm(n))

  g <- with(d, 2.5 * LAge + .7 * DM - .8 * LHDL + 2.5 * LSBP + .5 * Smoke - 19.415)
  d$Z <- 1 - .9^exp(g)

  d$A <- rbinom(n, 1, .5)

  T_ <- rexp(n, rate = .05 * exp(beta_A * d$A + beta_Z * d$Z))
  C_ <- runif(n, 1, 4)

  d$time <- pmin(T_, C_)
  d$status <- as.integer(T_ <= C_)

  d
}

specs <- list(
  oracle = A ~ Z,
  over = A ~ LAge + DM + LHDL + LSBP + Smoke + N1 + N2 + N3,
  omit = A ~ LAge + DM + LHDL + Smoke + N1,
  wrong = A ~ I(LAge^2) + DM + LHDL + LSBP + Smoke
)

methods <- c(EB = "ebal", IPW = "glm")

one_rep <- function(i) {
  if (!"WeightIt" %in% loadedNamespaces()) {
    if (version == "dev") {
      suppressMessages(pkgload::load_all(root, quiet = TRUE, export_all = FALSE))
    }
    else if (version == "v1709") {
      # 1.7.0.9000, the version the authors cite (commit 1a8324c), installed with
      #   git archive 1a8324c | tar -x -C <src>; R CMD INSTALL -l <lib1709> <src>
      library(WeightIt, lib.loc = Sys.getenv("WEIGHTIT_LIB1709"))
    }
    else {
      library(WeightIt)
    }
  }

  d <- gen_data()

  out <- vector("list", length(specs) * length(methods))
  k <- 0L

  for (s in names(specs)) {
    for (m in names(methods)) {
      k <- k + 1L

      W <- WeightIt::weightit(specs[[s]], data = d, method = methods[[m]],
                              estimand = "ATE")

      fit <- WeightIt::coxph_weightit(survival::Surv(time, status) ~ A,
                                      data = d, weightit = W)

      out[[k]] <- data.frame(rep = i, spec = s, method = m,
                             est = unname(coef(fit)["A"]),
                             se = sqrt(vcov(fit)["A", "A"]),
                             vcov_type = fit$vcov_type %||% NA_character_,
                             n_events = sum(d$status))
    }
  }

  do.call("rbind", out)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

set.seed(20261123)

pr <- prog_init(total = n_sim,
                title = sprintf("Zhang et al. Cox M-estimation SEs (%s)", version),
                unit = "replicate", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

res <- prog_parallel(pr, seq_len(n_sim), one_rep,
                     label = function(i) paste0("rep ", i),
                     workers = 7L, schedule = "dynamic", chunk_size = 10L,
                     checkpoint = file.path(root, "_dev", "coxph-mest-zhang2026-parts",
                                            paste0(version, "-", n_sim)),
                     packages = c("MASS", "survival"))

on.exit()
prog_end(pr, "done")

res <- do.call("rbind", res)
res$version <- version

saveRDS(res, file.path(root, "_dev", sprintf("coxph-mest-zhang2026-%s.rds", version)))
