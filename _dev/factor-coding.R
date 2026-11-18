# How the tree-based methods should receive a categorical covariate: as a factor
# (split on subsets of its levels), as a dummy variable for each level, or as
# integer codes (split as if ordered). Each coding is fit with the factor's levels
# in their original order and in a permuted order, to check whether the weights
# depend on how the factor is coded, and is assessed on balance and ESS. Run from
# the package root:
#
#   Rscript _dev/factor-coding.R
#
# Set BARTISAN_LIB to a library holding the bartisan version to test; QUICK=1 for
# a smoke test. Results are checkpointed after every fit to
# `_dev/factor-coding.rds` (or `..._quick.rds`).

if (nzchar(Sys.getenv("BARTISAN_LIB"))) {
  .libPaths(c(Sys.getenv("BARTISAN_LIB"), .libPaths()))
}

quick <- nzchar(Sys.getenv("QUICK"))

source(file.path(Sys.getenv("HOME"), ".claude/skills/live-progress/assets/progress.R"))

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
  library("cobalt")
})

out_file <- if (quick) "_dev/factor-coding_quick.rds" else "_dev/factor-coding.rds"

# --- Data -------------------------------------------------------------------------
# A six-level factor whose effect on treatment either alternates across the levels
# in their original order (so no ordering of the codes captures it) or rises with
# that order

effects <- list(
  "non-monotone" = c(a = 1, b = -1, c = 1, d = -1, e = .5, f = -.5),
  "monotone" = c(a = -1, b = -.6, c = -.2, d = .2, e = .6, f = 1)
)

gen_data <- function(dgp, n) {
  F <- factor(sample(letters[1:6], n, replace = TRUE), levels = letters[1:6])
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  A <- rbinom(n, 1, plogis(-.3 + .6 * X1 - .4 * X2 + effects[[dgp]][as.character(F)]))
  data.frame(A, X1, X2, F)
}

# Balance on X1, X2, and an indicator for every level of the factor, which is the
# same whatever the coding
assess <- function(w, d) {
  covs <- cbind(X1 = d$X1, X2 = d$X2, model.matrix(~ F - 1, data = d))
  smd <- abs(col_w_smd(covs, d$A, weights = w, s.d.denom = "pooled"))
  c(max_smd = max(smd), mean_smd = mean(smd),
    max_smd_F = max(smd[startsWith(colnames(covs), "F")]),
    ess0 = ESS(w[d$A == 0]), ess1 = ESS(w[d$A == 1]))
}

# --- Fits -------------------------------------------------------------------------

gbm_args <- {
  if (quick) list(criterion = "smd.mean", n.trees = 200)
  else list(criterion = "smd.mean", n.trees = 3000)
}

bart_args <- {
  if (quick) list(num_trees = 10L, num_draws = 30L, num_burn = 20L)
  else list()
}

fit_weights <- function(method, coding, d) {
  args <- if (method == "gbm") gbm_args else bart_args

  switch(coding,
         factor = do.call("weightit", c(list(A ~ X1 + X2 + F, data = d, method = method),
                                        args))$weights,
         onehot = do.call("weightit", c(list(A ~ X1 + X2 + F, data = d, method = method,
                                             categorical = "onehot"), args))$weights,
         dummies = do.call("weightit.fit",
                           c(list(cbind(X1 = d$X1, X2 = d$X2, model.matrix(~ F - 1, data = d)),
                                  treat = d$A, method = method), args))$weights,
         codes = do.call("weightit", c(list(A ~ X1 + X2 + Fc, data = transform(d, Fc = as.integer(F)),
                                            method = method), args))$weights)
}

codings <- list(gbm = c("factor", "dummies", "codes"),
                bart = c("factor", "onehot", "dummies"))

reps <- if (quick) 1L else 10L
n <- if (quick) 300L else 1000L

design <- do.call("rbind", lapply(names(codings), function(m) {
  expand.grid(method = m, dgp = names(effects), rep = seq_len(reps),
              coding = codings[[m]], order = c("original", "permuted"),
              stringsAsFactors = FALSE)
}))

pr <- prog_init(total = nrow(design),
                title = if (quick) "Factor coding (quick)" else "Factor coding for gbm and bart",
                unit = "fit", kind = "simulation")
on.exit(prog_end(pr, "failed"), add = TRUE)

rows <- vector("list", nrow(design))
weights <- vector("list", nrow(design))

for (k in seq_len(nrow(design))) {
  cfg <- design[k, ]
  label <- sprintf("%s | %s | rep %d | %s | %s", cfg$method, cfg$dgp, cfg$rep, cfg$coding, cfg$order)

  res <- prog_do(pr, i = k, label = label, f = function(i) {
    set.seed(100L * cfg$rep + match(cfg$dgp, names(effects)))
    d <- gen_data(cfg$dgp, n)

    if (cfg$order == "permuted") {
      set.seed(cfg$rep)
      d$F <- factor(d$F, levels = sample(levels(d$F)))
    }

    set.seed(cfg$rep)
    t <- system.time(w <- suppressWarnings(fit_weights(cfg$method, cfg$coding, d)))[["elapsed"]]

    list(stats = c(assess(w, d), seconds = t), w = unname(w))
  })

  if (!is.null(res)) {
    rows[[k]] <- cbind(cfg, as.data.frame(as.list(res$stats)))
    weights[[k]] <- res$w
  }

  saveRDS(list(results = do.call("rbind", rows), weights = weights, design = design,
               complete = FALSE, done = k, total = nrow(design)), out_file)
}

results <- do.call("rbind", rows)
saveRDS(list(results = results, weights = weights, design = design, complete = TRUE,
             done = nrow(design), total = nrow(design)), out_file)

on.exit()
prog_end(pr, "done")

# --- Summaries --------------------------------------------------------------------

options(width = 200)

summ <- aggregate(cbind(max_smd, mean_smd, max_smd_F, ess0, ess1, seconds) ~
                    method + dgp + coding + order, data = results, FUN = mean)
num <- vapply(summ, is.numeric, logical(1L))
summ[num] <- lapply(summ[num], function(x) round(x, 3))
print(summ[order(summ$method, summ$dgp, summ$coding, summ$order), ], row.names = FALSE)

# Whether the weights change when the levels are reordered: the largest absolute
# difference between the original-order and permuted-order weights, by rep
key <- with(design, paste(method, dgp, rep, coding))
diffs <- tapply(seq_len(nrow(design)), key, function(ii) {
  w <- weights[ii]
  if (length(w) == 2L && !is.null(w[[1L]]) && !is.null(w[[2L]])) max(abs(w[[1L]] - w[[2L]]))
  else NA_real_
})
inv <- data.frame(key = names(diffs), max_abs_diff = as.numeric(diffs))
inv$method <- sub(" .*", "", inv$key)
inv$coding <- sub(".* ", "", inv$key)
print(aggregate(max_abs_diff ~ method + coding, data = inv,
                FUN = function(x) signif(max(x), 3)), row.names = FALSE)
