# Tests for the bartisan-based `method = "bart"` in `_dev/weightit2bart_bartisan.R`.
#
# Run from the package root with
#
#   Rscript _dev/test-weightit2bart_bartisan.R
#
# The script loads the development version of WeightIt, replaces the dbarts-based
# `weightit2bart*()` functions in the namespace with the bartisan-based ones, patches
# the `.weightit_methods` entry for `"bart"` to match, and then runs testthat
# expectations against `weightit()` and `weightitMSM()`. The final section compares
# the two implementations on the lalonde data. Once the new file replaces
# `R/weightit2bart.R`, these expectations belong in `tests/testthat/test-method_bart.R`;
# see `_dev/bartisan-migration.md`.

suppressPackageStartupMessages({
  library("testthat")
  library("cobalt")
})

source("_dev/bartisan-inject.R")

# --- Helpers ---------------------------------------------------------------------

test_data <- readRDS("tests/testthat/fixtures/test_data.rds")

max_smd <- function(W) {
  b <- bal.tab(W, un = TRUE, s.d.denom = "pooled")

  # Binary treatment: one balance table
  if (!is.null(b$Balance) && !is.null(b$Balance$Diff.Un)) {
    return(c(un = max(abs(b$Balance$Diff.Un), na.rm = TRUE),
             adj = max(abs(b$Balance$Diff.Adj), na.rm = TRUE)))
  }

  # Multi-category treatment: the summary across pairs
  if (!is.null(b$Balance.Across.Pairs)) {
    return(c(un = max(abs(b$Balance.Across.Pairs$Max.Diff.Un), na.rm = TRUE),
             adj = max(abs(b$Balance.Across.Pairs$Max.Diff.Adj), na.rm = TRUE)))
  }

  # Otherwise, aggregate over the pairwise tables
  un <- adj <- numeric(0)
  for (p in b$Pair.Balance) {
    un <- c(un, max(abs(p$Balance$Diff.Un), na.rm = TRUE))
    adj <- c(adj, max(abs(p$Balance$Diff.Adj), na.rm = TRUE))
  }
  c(un = max(un), adj = max(adj))
}

max_cor <- function(W) {
  c(un = max(abs(col_w_cov(W$covs, W$treat, std = TRUE, s.weights = W$s.weights))),
    adj = max(abs(col_w_cov(W$covs, W$treat, weights = W$weights, std = TRUE,
                            s.weights = W$s.weights))))
}

fast <- list(num_draws = 150L, num_burn = 100L)

# --- Binary treatments -------------------------------------------------------------

test_that("Binary treatment: basic fit, outputs, reproducibility", {
  set.seed(123)
  expect_no_condition({
    W0 <- do.call(weightit, c(list(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                                   data = test_data, method = "bart",
                                   estimand = "ATE", include.obj = TRUE), fast))
  })

  expect_s3_class(W0$obj, "bartisan_fit")
  expect_true(is.numeric(W0$ps))
  expect_true(all(W0$ps > 0 & W0$ps < 1))
  expect_true(all(is.finite(W0$weights) & W0$weights >= 0))
  expect_null(attr(W0, "Mparts", exact = TRUE))

  b <- max_smd(W0)
  expect_lt(b[["adj"]], b[["un"]])

  # set.seed() reproduces exactly
  set.seed(123)
  W0b <- do.call(weightit, c(list(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                                  data = test_data, method = "bart",
                                  estimand = "ATE"), fast))
  expect_equal(W0$weights, W0b$weights)

  # Different seed, different weights
  set.seed(124)
  W0c <- do.call(weightit, c(list(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                                  data = test_data, method = "bart",
                                  estimand = "ATE"), fast))
  expect_false(isTRUE(all.equal(W0$weights, W0c$weights)))
})

test_that("Binary treatment: dbarts-style arguments are translated and `seed` warns", {
  set.seed(1)
  W.native <- weightit(A ~ X1 + X2 + X3, data = test_data, method = "bart",
                       num_trees = 10L, num_draws = 60L, num_burn = 40L,
                       include.obj = TRUE)
  set.seed(1)
  W.old <- weightit(A ~ X1 + X2 + X3, data = test_data, method = "bart",
                    n.trees = 10L, n.samples = 60L, n.burn = 40L, n.chains = 1L,
                    n.threads = 2L, include.obj = TRUE)

  expect_equal(W.native$weights, W.old$weights)
  expect_equal(W.old$obj$control$num_trees, 10L)
  expect_equal(nrow(W.old$obj$eta[[1L]]), 60L)

  # A native name wins over its translation
  set.seed(1)
  W.both <- weightit(A ~ X1 + X2 + X3, data = test_data, method = "bart",
                     n.trees = 30L, num_trees = 10L, num_draws = 60L, num_burn = 40L,
                     include.obj = TRUE)
  expect_equal(W.both$obj$control$num_trees, 10L)

  expect_warning({
    weightit(A ~ X1 + X2 + X3, data = test_data, method = "bart",
             seed = 5L, num_draws = 30L, num_burn = 20L)
  }, "seed.*ignored")
})

test_that("Binary treatment: sampling weights are accepted and scaled", {
  set.seed(123)
  W.sw <- do.call(weightit, c(list(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                                   data = test_data, method = "bart",
                                   estimand = "ATE", s.weights = "SW",
                                   include.obj = TRUE), fast))
  expect_true(all(is.finite(W.sw$weights) & W.sw$weights >= 0))

  # The weights passed to bartisan have a mean of 1
  expect_equal(mean(W.sw$obj$prior_weights), 1)
  expect_equal(unname(W.sw$obj$prior_weights), test_data$SW / mean(test_data$SW))

  # Scaling the sampling weights by a constant does not change the fit
  set.seed(123)
  W.sw2 <- do.call(weightit, c(list(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                                    data = test_data, method = "bart",
                                    estimand = "ATE", s.weights = 10 * test_data$SW),
                               fast))
  expect_equal(W.sw$weights, W.sw2$weights)

  # Varying sampling weights change the weights relative to none
  set.seed(123)
  W.none <- do.call(weightit, c(list(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                                     data = test_data, method = "bart",
                                     estimand = "ATE"), fast))
  expect_false(isTRUE(all.equal(W.sw$weights, W.none$weights)))

  # Balance assessed with the sampling weights improves
  b <- max_smd(W.sw)
  expect_lt(b[["adj"]], b[["un"]])
})

test_that("Binary treatment: estimands, link, use.offset, subclass, missing data", {
  for (estimand in c("ATT", "ATC", "ATO", "ATM", "ATOS")) {
    set.seed(123)
    expect_no_condition({
      W <- do.call(weightit, c(list(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                                    method = "bart", estimand = estimand), fast))
    }, message = estimand)
    expect_true(all(is.finite(W$weights) & W$weights >= 0), label = estimand)
    if (estimand == "ATT") expect_equal(unname(W$weights[W$treat == 1]), rep(1, sum(W$treat == 1)))
    if (estimand == "ATC") expect_equal(unname(W$weights[W$treat == 0]), rep(1, sum(W$treat == 0)))
  }

  set.seed(123)
  expect_no_condition({
    W.probit <- do.call(weightit, c(list(A ~ X1 + X2 + X3, data = test_data,
                                         method = "bart", link = "probit",
                                         include.obj = TRUE), fast))
  })
  expect_equal(W.probit$obj$family$link, "probit")

  expect_error({
    weightit(A ~ X1 + X2 + X3, data = test_data, method = "bart", link = "cauchit")
  })

  set.seed(123)
  expect_no_condition({
    W.off <- do.call(weightit, c(list(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                                      method = "bart", use.offset = TRUE,
                                      include.obj = TRUE), fast))
  })
  expect_true(isTRUE(W.off$obj$has_offset))
  set.seed(123)
  W.nooff <- do.call(weightit, c(list(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                                      method = "bart"), fast))
  expect_false(isTRUE(all.equal(W.off$weights, W.nooff$weights)))

  set.seed(123)
  expect_no_condition({
    W.sub <- do.call(weightit, c(list(A ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                                      method = "bart", subclass = 10L), fast))
  })
  expect_lt(length(unique(W.sub$weights)), nrow(test_data) / 10)

  data_na <- test_data
  set.seed(4321)
  is.na(data_na$X1) <- sample(nrow(data_na), round(.05 * nrow(data_na)))
  set.seed(123)
  expect_no_condition({
    W.na <- do.call(weightit, c(list(A ~ X1 + X2 + X3, data = data_na,
                                     method = "bart", missing = "ind"), fast))
  })
  expect_true(anyNA(W.na$covs$X1))
  expect_true(all(W.na$ps > 0 & W.na$ps < 1))
})

# --- Multi-category treatments -----------------------------------------------------

test_that("Multi-category treatment: single multinomial fit, estimands, weights", {
  set.seed(123)
  expect_no_condition({
    W0 <- do.call(weightit, c(list(Am ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                                   method = "bart", estimand = "ATE",
                                   include.obj = TRUE), fast))
  })

  expect_null(W0$ps)
  expect_s3_class(W0$obj, "bartisan_fit")
  expect_equal(W0$obj$family$family, "multinomial")

  p <- predict(W0$obj, type = "prob")
  expect_equal(dim(p), c(nrow(test_data), 3L))
  expect_equal(unname(rowSums(p)), rep(1, nrow(test_data)))
  expect_true(all(is.finite(W0$weights) & W0$weights > 0))

  b <- max_smd(W0)
  expect_lt(b[["adj"]], b[["un"]])

  set.seed(123)
  W.att <- do.call(weightit, c(list(Am ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                                    method = "bart", estimand = "ATT", focal = "T"),
                               fast))
  expect_equal(unname(W.att$weights[W.att$treat == "T"]), rep(1, sum(W.att$treat == "T")))

  for (estimand in c("ATO", "ATM")) {
    set.seed(123)
    expect_no_condition({
      W <- do.call(weightit, c(list(Am ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                                    method = "bart", estimand = estimand), fast))
    }, message = estimand)
    expect_true(all(is.finite(W$weights) & W$weights >= 0), label = estimand)
  }

  expect_error({
    weightit(Am ~ X1 + X2 + X3, data = test_data, method = "bart", estimand = "ATOS")
  }, "not an allowable estimand", ignore.case = TRUE)

  expect_error({
    weightit(Am ~ X1 + X2 + X3, data = test_data, method = "bart", use.offset = TRUE)
  }, "use.offset")

  # Sampling weights, reference level, probit link
  set.seed(123)
  expect_no_condition({
    W.sw <- do.call(weightit, c(list(Am ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                                     method = "bart", s.weights = "SW",
                                     include.obj = TRUE), fast))
  })
  expect_equal(mean(W.sw$obj$prior_weights), 1)

  set.seed(123)
  expect_no_condition({
    W.ref <- do.call(weightit, c(list(Am ~ X1 + X2 + X3, data = test_data,
                                      method = "bart", reference = "T",
                                      include.obj = TRUE), fast))
  })
  expect_equal(length(W.ref$obj$eta), 2L)

  set.seed(123)
  expect_no_condition({
    W.probit <- weightit(Am ~ X1 + X2 + X3, data = test_data, method = "bart",
                         link = "probit", replicates = 50L,
                         num_draws = 40L, num_burn = 40L, include.obj = TRUE)
  })
  expect_equal(W.probit$obj$family$link, "probit")
  expect_true(all(is.finite(W.probit$weights) & W.probit$weights > 0))
})

# --- Continuous treatments ---------------------------------------------------------

test_that("Continuous treatment: DPM density, families, sampling weights", {
  set.seed(123)
  expect_no_condition({
    W0 <- do.call(weightit, c(list(Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                                   data = test_data, method = "bart",
                                   include.obj = TRUE), fast))
  })

  expect_s3_class(W0$obj, "bartisan_fit")
  expect_equal(W0$obj$family$family, "dpm")
  expect_true(all(is.finite(W0$weights) & W0$weights > 0))

  r <- max_cor(W0)
  expect_lt(r[["adj"]], r[["un"]])

  # Weights have a mean near 1, as stabilized GPS weights should
  expect_lt(abs(mean(W0$weights) - 1), .25)

  # The density arguments of the other methods are ignored with a warning
  set.seed(123)
  expect_warning({
    weightit(Ac ~ X1 + X2 + X3, data = test_data, method = "bart",
             density = "kernel", num_draws = 40L, num_burn = 40L)
  }, "ignored")

  # gaussian family uses the sampling weights in the fit
  set.seed(123)
  expect_no_condition({
    W.g <- do.call(weightit, c(list(Ac ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                                    method = "bart", family = gaussian(),
                                    s.weights = "SW", include.obj = TRUE), fast))
  })
  expect_equal(W.g$obj$family$family, "gaussian")
  expect_equal(mean(W.g$obj$prior_weights), 1)
  expect_false(all(W.g$obj$prior_weights == 1))
  expect_true(all(is.finite(W.g$weights) & W.g$weights > 0))
  expect_false(is.null(W.g$obj$call$weights))

  # dpm with sampling weights: fit without them, but they enter the numerator
  set.seed(123)
  expect_no_condition({
    W.dsw <- do.call(weightit, c(list(Ac ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                                      method = "bart", s.weights = "SW",
                                      include.obj = TRUE), fast))
  })
  expect_true(is.null(W.dsw$obj$call$weights))
  set.seed(123)
  W.d <- do.call(weightit, c(list(Ac ~ X1 + X2 + X3 + X4 + X5, data = test_data,
                                  method = "bart"), fast))
  expect_false(isTRUE(all.equal(W.dsw$weights, W.d$weights)))

  # gaussian_ls and a family given by name
  set.seed(123)
  expect_no_condition({
    W.ls <- do.call(weightit, c(list(Ac ~ X1 + X2 + X3, data = test_data,
                                     method = "bart", family = bartisan::gaussian_ls(),
                                     include.obj = TRUE), fast))
  })
  expect_true(all(is.finite(W.ls$weights) & W.ls$weights > 0))

  set.seed(123)
  expect_no_condition({
    W.nm <- do.call(weightit, c(list(Ac ~ X1 + X2 + X3, data = test_data,
                                     method = "bart", family = "gaussian"), fast))
  })
  expect_true(all(is.finite(W.nm$weights) & W.nm$weights > 0))

  expect_error({
    weightit(Ac ~ X1 + X2 + X3, data = test_data, method = "bart", family = 3)
  }, "family")

  # use.offset
  set.seed(123)
  expect_no_condition({
    W.off <- do.call(weightit, c(list(Ac ~ X1 + X2 + X3, data = test_data,
                                      method = "bart", use.offset = TRUE,
                                      include.obj = TRUE), fast))
  })
  expect_true(isTRUE(W.off$obj$has_offset))
  expect_true(all(is.finite(W.off$weights) & W.off$weights > 0))
})

# --- Censoring weights -------------------------------------------------------------

test_that("Censoring weights", {
  d <- test_data
  set.seed(55)
  d$C <- rbinom(nrow(d), 1, plogis(-1.5 + .5 * d$X1 + .3 * d$A))

  set.seed(123)
  expect_no_condition({
    Wc <- do.call(weightit, c(list(.cens(C) ~ X1 + X2 + X3 + X4 + A, data = d,
                                   method = "bart", include.obj = TRUE), fast))
  })
  expect_true(all(Wc$weights[d$C == 1] == 0))
  expect_true(all(Wc$weights[d$C == 0] > 0))
  expect_equal(unname(Wc$weights[d$C == 0]), unname(1 / (1 - Wc$ps[d$C == 0])))

  # Sampling weights with censoring
  set.seed(123)
  expect_no_condition({
    Wcs <- do.call(weightit, c(list(.cens(C) ~ X1 + X2 + X3 + X4 + A, data = d,
                                    method = "bart", s.weights = "SW"), fast))
  })
  expect_true(all(Wcs$weights[d$C == 1] == 0))
})

# --- Random intercepts -------------------------------------------------------------

test_that("Random intercepts are fit by bartisan directly", {
  set.seed(123)
  expect_no_condition({
    W.re <- do.call(weightit, c(list(A ~ X1 + X2 + X3 + (1 | cluster), data = test_data,
                                     method = "bart", include.obj = TRUE), fast))
  })
  expect_s3_class(W.re$obj, "bartisan_fit")
  expect_false(is.null(W.re$obj$ranef))
  expect_true(all(W.re$ps > 0 & W.re$ps < 1))

  set.seed(123)
  expect_no_condition({
    W.rec <- do.call(weightit, c(list(Ac ~ X1 + X2 + X3 + (1 | cluster), data = test_data,
                                      method = "bart", include.obj = TRUE), fast))
  })
  expect_true(all(is.finite(W.rec$weights) & W.rec$weights > 0))

  set.seed(123)
  expect_no_condition({
    W.rem <- do.call(weightit, c(list(Am ~ X1 + X2 + X3 + (1 | cluster), data = test_data,
                                      method = "bart", include.obj = TRUE), fast))
  })
  expect_true(all(is.finite(W.rem$weights) & W.rem$weights > 0))

  # A random slope is refused by bartisan, and the error is relayed
  expect_error({
    weightit(A ~ X1 + X2 + (1 + X1 | cluster), data = test_data, method = "bart",
             num_draws = 20L, num_burn = 20L)
  }, "bartisan")
})

# --- Longitudinal treatments -------------------------------------------------------

test_that("weightitMSM() works with the new method", {
  data("msmdata")
  set.seed(123)
  expect_no_condition({
    Wmsm <- do.call(weightitMSM, c(list(list(A_1 ~ X1_0 + X2_0,
                                             A_2 ~ X1_1 + X2_1 + A_1 + X1_0 + X2_0,
                                             A_3 ~ X1_2 + X2_2 + A_2 + X1_1 + X2_1 + A_1 + X1_0 + X2_0),
                                        data = msmdata, method = "bart",
                                        stabilize = TRUE), fast))
  })
  expect_true(all(is.finite(Wmsm$weights) & Wmsm$weights > 0))
  expect_lt(abs(mean(Wmsm$weights) - 1), .25)
})

# --- Comparison with the dbarts-based implementation ---------------------------------

cat("\n\n=== Comparison with the dbarts implementation on lalonde ===\n")
data("lalonde", package = "cobalt")

run_old <- function(expr) {
  # Temporarily restore the old functions
  for (fn in c("weightit2bart", "weightit2bart.multi", "weightit2bart.cont", "weightit2bart.cens")) {
    assignInNamespace(fn, get(fn, envir = old_env), ns = "WeightIt")
  }
  on.exit({
    for (fn in c("weightit2bart", "weightit2bart.multi", "weightit2bart.cont", "weightit2bart.cens")) {
      assignInNamespace(fn, get(fn, envir = new_env), ns = "WeightIt")
    }
  }, add = TRUE)
  methods_old <- ns$.weightit_methods
  methods_old$bart$packages_needed <- "dbarts"
  methods_old$bart$s.weights_ok <- FALSE
  assignInNamespace(".weightit_methods", methods_old, ns = "WeightIt")
  on.exit(assignInNamespace(".weightit_methods", methods, ns = "WeightIt"), add = TRUE)
  force(expr)
}

compare <- NULL
add_row <- function(label, old, new, t_old, t_new, bal_fun) {
  bo <- bal_fun(old); bn <- bal_fun(new)
  compare <<- rbind(compare, data.frame(
    case = label,
    max_un = round(bo[["un"]], 3),
    max_old = round(bo[["adj"]], 3), max_new = round(bn[["adj"]], 3),
    ess_old = round(sum(ESS_by_group(old))), ess_new = round(sum(ESS_by_group(new))),
    sec_old = round(t_old, 1), sec_new = round(t_new, 1)))
}
ESS_by_group <- function(W) {
  if (is.factor(W$treat) || length(unique(W$treat)) <= 5) {
    sapply(split(W$weights, W$treat), ESS)
  } else {
    ESS(W$weights)
  }
}

f_bin <- treat ~ age + educ + married + nodegree + re74 + re75 + race
f_mul <- race ~ age + educ + married + nodegree + re74 + re75 + treat
f_con <- re75 ~ age + educ + married + nodegree + re74 + race

for (case in c("binary ATE", "binary ATT", "multi ATE", "continuous")) {
  f <- switch(case, "binary ATE" = f_bin, "binary ATT" = f_bin, "multi ATE" = f_mul, f_con)
  est <- if (grepl("ATT", case)) "ATT" else "ATE"
  set.seed(1)
  t_old <- system.time(old <- run_old(weightit(f, data = lalonde, method = "bart", estimand = est,
                                               n.chains = 1, n.threads = 1)))[["elapsed"]]
  set.seed(1)
  t_new <- system.time(new <- weightit(f, data = lalonde, method = "bart", estimand = est))[["elapsed"]]
  add_row(case, old, new, t_old, t_new, if (case == "continuous") max_cor else max_smd)
}

print(compare, row.names = FALSE)
cat("\nmax_un: largest absolute SMD (or treatment-covariate correlation) before weighting;",
    "max_old/max_new: after weighting with dbarts / bartisan; ess: summed effective sample size;",
    "sec: fit time in seconds (dbarts single-threaded, bartisan one chain, both at their defaults).\n")
