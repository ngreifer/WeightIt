# Compare three ways of handling sampling weights with method = "bart" and a
# continuous treatment: balance, effective sample size, and time.
#
#   - "dpm": the default, dpm(), which refuses prior weights, so the sampling
#     weights enter only the marginal density in the numerator of the weights
#   - "dpm + SW covariate": dpm() with the sampling weights also added to the
#     model as a covariate
#   - "gaussian": gaussian(), which takes the sampling weights as prior weights in
#     the fit, as well as in the marginal density
#
# Balance is assessed in the sample weighted by the sampling weights, which is the
# target. Run from the package root:
#
#   Rscript _dev/bart-sampling-weights.R
#
# Set BARTISAN_LIB to a library holding the bartisan version to test (e.g., the
# CRAN release); otherwise the installed one is used. Set QUICK=1 for a smoke test
# with tiny settings. Results are checkpointed after every fit to
# `_dev/bart-sampling-weights.rds` (or `..._quick.rds`).
#
# The results saved on 2026-11-18 came from the version of the method in which
# dpm() used varying sampling weights only in the marginal density. Since then,
# that is an error, so the two dpm() variants no longer run as written.

if (nzchar(Sys.getenv("BARTISAN_LIB"))) {
  .libPaths(c(Sys.getenv("BARTISAN_LIB"), .libPaths()))
}

quick <- nzchar(Sys.getenv("QUICK"))

source(file.path(Sys.getenv("HOME"), ".claude/skills/live-progress/assets/progress.R"))

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("cobalt")
})

cat("bartisan", as.character(packageVersion("bartisan")), "from", find.package("bartisan"), "\n")

out_file <- if (quick) "_dev/bart-sampling-weights_quick.rds" else "_dev/bart-sampling-weights.rds"

# --- Data-generating processes --------------------------------------------------
# Six covariates; the sampling weights depend on X2 and X4 in the first three, and
# also on the treatment itself in the fourth, as when selection into the sample
# depends on the treatment.

gen_data <- function(dgp, n) {
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rnorm(n)
  X4 <- rbinom(n, 1, .5)
  X5 <- runif(n, -1, 1)
  X6 <- rbinom(n, 1, plogis(X1))

  m_lin <- .5 * X1 + .5 * X2 - .4 * X3 + .6 * X4 + .3 * X5 + .2 * X6
  m_nl <- .5 * X1 + .6 * (X2^2 - 1) + sin(2 * X3) + .6 * X4 * X1 + .3 * X5 + .4 * X6 * X2
  s_het <- exp(.4 * X1 - .3 * X4 + .4 * X5)

  A <- switch(dgp,
              "linear, normal; SW on X" = m_lin + rnorm(n),
              "nonlinear, heteroskedastic; SW on X" = m_nl + s_het * rnorm(n),
              "nonlinear, skewed; SW on X" = m_nl + (rchisq(n, 3) - 3) / sqrt(6),
              "nonlinear, normal; SW on X and A" = m_nl + rnorm(n))

  SW <- {
    if (dgp == "nonlinear, normal; SW on X and A") exp(.4 * X2 + .5 * (A - mean(A)) / sd(A))
    else exp(.6 * X2 - .5 * X4)
  }

  data.frame(A, X1, X2, X3, X4, X5, X6, SW = SW / mean(SW))
}

# Covariates whose correlation with the treatment is assessed: the main terms and
# the nonlinear terms the DGPs use
sim_covs <- function(d) {
  with(d, cbind(X1, X2, X3, X4, X5, X6,
                `X1^2` = X1^2, `X2^2` = X2^2, `X3^2` = X3^2, `sin(2*X3)` = sin(2 * X3),
                `X1*X4` = X1 * X4, `X2*X6` = X2 * X6))
}

dgps <- c("linear, normal; SW on X", "nonlinear, heteroskedastic; SW on X",
          "nonlinear, skewed; SW on X", "nonlinear, normal; SW on X and A")

# --- Lalonde, with synthetic sampling weights -----------------------------------

data("lalonde", package = "cobalt")
lalonde$black <- as.numeric(lalonde$race == "black")
lalonde$hispan <- as.numeric(lalonde$race == "hispan")
lalonde$SW <- with(lalonde, exp(.5 * (educ - mean(educ)) / sd(educ) - .6 * black + .4 * married))
lalonde$SW <- lalonde$SW / mean(lalonde$SW)

lalonde_cases <- list(
  "lalonde: age" = c("educ", "black", "hispan", "married", "nodegree", "re74", "re75"),
  "lalonde: re75" = c("age", "educ", "black", "hispan", "married", "nodegree", "re74")
)

lalonde_covs <- function(vars) {
  X <- as.matrix(lalonde[vars])
  cont <- vars[vapply(lalonde[vars], function(x) length(unique(x)) > 2L, logical(1L))]
  sq <- X[, cont, drop = FALSE]^2
  colnames(sq) <- paste0(cont, "^2")
  cbind(X, sq)
}

# --- Assessment -----------------------------------------------------------------

assess <- function(W, covs, n_main, s.weights) {
  r_un <- abs(col_w_corr(covs, W$treat, s.weights = s.weights))
  r_adj <- abs(col_w_corr(covs, W$treat, weights = W$weights, s.weights = s.weights))

  w <- W$weights * s.weights

  c(max_r_un = max(r_un[seq_len(n_main)]),
    max_r = max(r_adj[seq_len(n_main)]),
    mean_r = mean(r_adj[seq_len(n_main)]),
    max_r_nl = max(r_adj),
    ess = sum(w)^2 / sum(w^2),
    max_w = max(W$weights) / mean(W$weights))
}

# --- The design -----------------------------------------------------------------

sampler <- {
  if (quick) list(num_trees = 10L, num_draws = 30L, num_burn = 20L)
  else list()
}

variants <- c("dpm", "dpm + SW covariate", "gaussian")

reps <- if (quick) 1L else 3L
n_sim <- if (quick) 300L else 1000L

design <- rbind(
  expand.grid(study = "dgp", case = dgps, n = n_sim, rep = seq_len(reps),
              variant = variants, stringsAsFactors = FALSE),
  expand.grid(study = "lalonde", case = names(lalonde_cases), n = nrow(lalonde),
              rep = seq_len(reps), variant = variants, stringsAsFactors = FALSE)
)

# --- Run ------------------------------------------------------------------------

pr <- prog_init(total = nrow(design),
                title = if (quick) "BART sampling weights comparison (quick)"
                        else "BART sampling weights comparison",
                unit = "fit", kind = "benchmark")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- vector("list", nrow(design))

for (k in seq_len(nrow(design))) {
  cfg <- design[k, ]
  label <- sprintf("%s | n = %d | rep %d | %s", cfg$case, cfg$n, cfg$rep, cfg$variant)

  res <- prog_do(pr, i = k, label = label, f = function(i) {
    # The same data for every variant within a case and replication
    set.seed(1000L * cfg$rep + match(cfg$case, c(dgps, names(lalonde_cases))))

    if (cfg$study == "lalonde") {
      vars <- lalonde_cases[[cfg$case]]
      d <- lalonde
      treat <- sub("lalonde: ", "", cfg$case, fixed = TRUE)
      covs <- lalonde_covs(vars)
      nm <- length(vars)
    }
    else {
      d <- gen_data(cfg$case, cfg$n)
      treat <- "A"
      vars <- paste0("X", 1:6)
      covs <- sim_covs(d)
      nm <- length(vars)
    }

    if (cfg$variant == "dpm + SW covariate") {
      vars <- c(vars, "SW")
    }

    f <- reformulate(vars, response = treat)
    family <- if (cfg$variant == "gaussian") "gaussian" else "dpm"

    # The same sampler seed for every variant
    set.seed(cfg$rep)
    t <- system.time({
      W <- do.call("weightit", c(list(f, data = quote(d), method = "bart",
                                      s.weights = d$SW, family = family),
                                 sampler))
    })[["elapsed"]]

    c(assess(W, covs, nm, d$SW), seconds = t)
  })

  if (!is.null(res)) {
    rows[[k]] <- cbind(cfg, as.data.frame(as.list(res)))
  }

  # Checkpoint after every fit
  saveRDS(list(results = do.call("rbind", rows), complete = FALSE,
               done = k, total = nrow(design)), out_file)
}

results <- do.call("rbind", rows)
saveRDS(list(results = results, complete = TRUE, done = nrow(design),
             total = nrow(design)), out_file)

on.exit()
prog_end(pr, "done")

# --- Summaries ------------------------------------------------------------------

summ <- aggregate(cbind(max_r_un, max_r, mean_r, max_r_nl, ess, max_w, seconds) ~
                    study + case + variant,
                  data = results, FUN = mean)
summ <- summ[order(summ$study, summ$case, summ$variant), ]
num <- vapply(summ, is.numeric, logical(1L))
summ[num] <- lapply(summ[num], function(x) round(x, 3))
options(width = 200)
print(summ, row.names = FALSE)
