# Compare gaussian_ls() and dpm() as the model for continuous treatments with
# method = "bart": balance, effective sample size, and time, on simulated data from
# several data-generating processes and on lalonde. Run from the package root:
#
#   Rscript _dev/bart-continuous-family.R
#
# Set BARTISAN_LIB to a library holding the bartisan version to test (e.g., the
# CRAN release); otherwise the installed one is used. Set QUICK=1 for a smoke test
# with tiny settings. Results are checkpointed after every fit to
# `_dev/bart-continuous-family.rds` (or `..._quick.rds`); set RESUME=1 to run only
# the configurations missing from that file (e.g., ones that failed) and add them.
#
# The results saved on 2026-11-18 came from an earlier version of this script,
# written before `family` existed, which swapped the continuous method in the
# namespace. Its `gaussian_ls()` variant computed both densities in closed form from
# `predict(type = "link", draws = TRUE)` rather than with `predict(type =
# "density")` on a grid, so its times there are a few seconds shorter than this
# script gives, and its `dpm()` variant interpolated the marginal density itself
# rather than its log.

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

out_file <- if (quick) "_dev/bart-continuous-family_quick.rds" else "_dev/bart-continuous-family.rds"

# --- Data-generating processes --------------------------------------------------
# Six covariates; the treatment's mean is linear or nonlinear in them, and its
# error is normal, heteroskedastic normal, skewed, or bimodal. The last DGP adds
# sampling weights that depend on the covariates.

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

  e <- switch(dgp,
              "linear, normal" = ,
              "nonlinear, normal" = rnorm(n),
              "nonlinear, heteroskedastic" = ,
              "heteroskedastic, s.weights" = s_het * rnorm(n),
              "nonlinear, skewed" = (rchisq(n, 3) - 3) / sqrt(6),
              "nonlinear, bimodal" = ifelse(rbinom(n, 1, .5) == 1, -1, 1) * 1.2 + rnorm(n, sd = .5))

  m <- if (dgp == "linear, normal") m_lin else m_nl

  d <- data.frame(A = m + e, X1, X2, X3, X4, X5, X6)
  d$SW <- {
    if (dgp == "heteroskedastic, s.weights") exp(.6 * X2 - .5 * X4)
    else 1
  }

  d
}

# Covariates whose correlation with the treatment is assessed: the main terms and
# the nonlinear terms the DGPs use
sim_covs <- function(d) {
  with(d, cbind(X1, X2, X3, X4, X5, X6,
                `X1^2` = X1^2, `X2^2` = X2^2, `X3^2` = X3^2, `sin(2*X3)` = sin(2 * X3),
                `X1*X4` = X1 * X4, `X2*X6` = X2 * X6))
}

n_main <- 6L

# --- Lalonde --------------------------------------------------------------------

data("lalonde", package = "cobalt")
lalonde$black <- as.numeric(lalonde$race == "black")
lalonde$hispan <- as.numeric(lalonde$race == "hispan")

lalonde_cases <- list(
  "lalonde: re75" = list(f = re75 ~ age + educ + black + hispan + married + nodegree + re74,
                         vars = c("age", "educ", "black", "hispan", "married", "nodegree", "re74")),
  "lalonde: age" = list(f = age ~ educ + black + hispan + married + nodegree + re74 + re75,
                        vars = c("educ", "black", "hispan", "married", "nodegree", "re74", "re75"))
)

lalonde_covs <- function(vars) {
  X <- as.matrix(lalonde[vars])
  cont <- vars[vapply(lalonde[vars], function(x) length(unique(x)) > 2L, logical(1L))]
  cbind(X, setNames(as.data.frame(X[, cont, drop = FALSE]^2), paste0(cont, "^2")) |> as.matrix())
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

families <- c("gaussian_ls", "dpm")

dgps <- c("linear, normal", "nonlinear, normal", "nonlinear, heteroskedastic",
          "nonlinear, skewed", "nonlinear, bimodal", "heteroskedastic, s.weights")

reps <- if (quick) 1L else 3L
n_sim <- if (quick) 300L else 1000L
n_scale <- if (quick) 600L else c(500L, 2000L, 4000L)

design <- rbind(
  expand.grid(study = "dgp", case = dgps, n = n_sim, rep = seq_len(reps),
              family = families, stringsAsFactors = FALSE),
  expand.grid(study = "lalonde", case = names(lalonde_cases), n = nrow(lalonde),
              rep = seq_len(reps), family = families, stringsAsFactors = FALSE),
  expand.grid(study = "scaling", case = "nonlinear, heteroskedastic", n = n_scale,
              rep = 1L, family = families, stringsAsFactors = FALSE)
)

# --- Run ------------------------------------------------------------------------

key <- function(x) paste(x$study, x$case, x$n, x$rep, x$family, sep = "|")

previous <- NULL
if (nzchar(Sys.getenv("RESUME")) && file.exists(out_file)) {
  previous <- readRDS(out_file)[["results"]]
  design <- design[key(design) %nin% key(previous), , drop = FALSE]
  cat("Resuming:", nrow(design), "configurations to run\n")
}

pr <- prog_init(total = nrow(design),
                title = if (quick) "BART continuous family comparison (quick)"
                        else "BART continuous family comparison",
                unit = "fit", kind = "benchmark")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- vector("list", nrow(design))

for (k in seq_len(nrow(design))) {
  cfg <- design[k, ]
  label <- sprintf("%s | n = %d | rep %d | %s", cfg$case, cfg$n, cfg$rep, cfg$family)

  res <- prog_do(pr, i = k, label = label, f = function(i) {
    # The same data for both families within a case and replication
    set.seed(1000L * cfg$rep + cfg$n + match(cfg$case, c(dgps, names(lalonde_cases))))

    if (cfg$study == "lalonde") {
      lc <- lalonde_cases[[cfg$case]]
      d <- lalonde
      f <- lc$f
      covs <- lalonde_covs(lc$vars)
      nm <- length(lc$vars)
      sw <- rep(1, nrow(d))
    }
    else {
      d <- gen_data(cfg$case, cfg$n)
      f <- A ~ X1 + X2 + X3 + X4 + X5 + X6
      covs <- sim_covs(d)
      nm <- n_main
      sw <- d$SW
    }

    # The same sampler seed for both families
    set.seed(cfg$rep)
    t <- system.time({
      W <- do.call("weightit", c(list(f, data = quote(d), method = "bart",
                                      s.weights = sw, family = cfg$family),
                                 sampler))
    })[["elapsed"]]

    c(assess(W, covs, nm, sw), seconds = t)
  })

  if (!is.null(res)) {
    rows[[k]] <- cbind(cfg, as.data.frame(as.list(res)))
  }

  # Checkpoint after every fit
  saveRDS(list(results = do.call("rbind", c(list(previous), rows)), complete = FALSE,
               done = k, total = nrow(design)), out_file)
}

results <- do.call("rbind", c(list(previous), rows))
saveRDS(list(results = results, complete = TRUE, done = nrow(design),
             total = nrow(design)), out_file)

on.exit()
prog_end(pr, "done")

# --- Summaries ------------------------------------------------------------------

summ <- aggregate(cbind(max_r_un, max_r, mean_r, max_r_nl, ess, max_w, seconds) ~
                    study + case + n + family,
                  data = results, FUN = mean)
summ <- summ[order(summ$study, summ$case, summ$n, summ$family), ]
num <- vapply(summ, is.numeric, logical(1L))
summ[num] <- lapply(summ[num], function(x) round(x, 3))
options(width = 200)
print(summ, row.names = FALSE)
