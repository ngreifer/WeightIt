# Inverse probability of censoring weights (IPCW), requested by wrapping the
# censoring indicator in `.cens()` on the left side of a model formula.
#
# The target: weights on the units still under observation that make their
# covariate distribution match that of the FULL at-risk sample (censored and
# uncensored combined). Censored units get a weight of exactly 0, and the rest get
# 1/P(C = 0 | X).

skip_on_cran()

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

test_data <- readRDS(test_path("fixtures", "test_data.rds"))

# Censoring depends on X1 and X3, so the uncensored units are a biased subsample
make_cens_data <- function(seed = 123L, p = -0.9) {
  set.seed(seed)
  d <- test_data
  d$C <- rbinom(nrow(d), 1L, prob = plogis(p + 0.8 * d$X1 - 0.5 * d$X3))
  d$G <- factor(d$X5, labels = c("g0", "g1"))
  d
}

# Weighted covariate means among the uncensored vs. the s.weights-weighted means of
# the full at-risk sample: the quantity the balance-based methods drive to 0.
target_diff <- function(W, d, covs = c("X1", "X2", "X3", "X4")) {
  X <- as.matrix(d[covs])
  sw <- W$s.weights %or% rep.int(1, nrow(d))
  u <- which(.make_cens_treat(W$treat) == 0)

  tw <- W$weights[u] * sw[u]

  max(abs(colSums(X[u, , drop = FALSE] * tw) / sum(tw) -
            colSums(X * sw) / sum(sw)))
}

# ---- Formula marker and treatment type -------------------------------------

test_that(".cens() is detected, stripped, and tagged", {
  d <- make_cens_data()

  t.c <- get_covs_and_treat_from_formula2(.cens(C) ~ X1 + X2, d)

  # The marker is stripped, so treat.name is the indicator itself
  expect_identical(attr(t.c$treat, "treat.name"), "C")
  expect_identical(get_treat_type(t.c$treat), "censoring")

  # An ordinary formula on the same variable is still binary
  t.b <- get_covs_and_treat_from_formula2(C ~ X1 + X2, d)
  expect_identical(get_treat_type(as.treat(t.b$treat, process = TRUE)), "binary")

  # The zero-RHS early return also tags (a marginal censoring model)
  t.0 <- get_covs_and_treat_from_formula2(.cens(C) ~ 1, d)
  expect_identical(attr(t.0$treat, "treat.name"), "C")
  expect_identical(get_treat_type(t.0$treat), "censoring")

  # `[.treat` preserves the type through subsetting
  expect_identical(get_treat_type(t.c$treat[1:10]), "censoring")
})

test_that(".cens() validates the indicator", {
  d <- make_cens_data()

  # A single value is allowed: no unit (or every unit) censored is degenerate but
  # not malformed
  d0 <- d
  d0$C <- 0L
  expect_identical(get_treat_type(assign_treat_type(
    get_covs_and_treat_from_formula2(.cens(C) ~ X1, d0)$treat)), "censoring")

  # Values other than 0/1 are not
  d2 <- d
  d2$C[1L] <- 2L
  expect_error(assign_treat_type(
    get_covs_and_treat_from_formula2(.cens(C) ~ X1, d2)$treat),
    "only the values")

  # A factor whose levels are not 0/1 must error rather than becoming all-NA
  df <- d
  df$C <- factor(ifelse(d$C == 1L, "Yes", "No"))
  expect_error(assign_treat_type(
    get_covs_and_treat_from_formula2(.cens(C) ~ X1, df)$treat),
    "only the values")

  # Logical indicators are fine
  dl <- d
  dl$C <- as.logical(d$C)
  expect_identical(get_treat_type(assign_treat_type(
    get_covs_and_treat_from_formula2(.cens(C) ~ X1, dl)$treat)), "censoring")
})

test_that(".cens() misuse is diagnosable", {
  d <- make_cens_data()

  # A direct call is legitimate use as of the constructor change; only a bad
  # indicator errors (covered fully in the constructor tests above)
  expect_no_error(.cens(d$C))

  # More than one argument
  expect_error(get_covs_and_treat_from_formula2(.cens(C, X1) ~ X1, d),
               "exactly one censoring indicator")

  # A mistyped marker points at the right thing
  expect_error(get_covs_and_treat_from_formula2(.censor(C) ~ X1, d),
               "\\.cens")
})

test_that("methods without censoring support error informatively", {
  d <- make_cens_data()

  expect_error(weightit(.cens(C) ~ X1, data = d, method = "npcbps"),
               "cannot be used to estimate censoring weights")
})

# ---- .cens() as a constructor, and the lower-level entry points -------------

test_that(".cens() returns a tagged censoring indicator", {
  d <- make_cens_data()

  C <- .cens(d$C)

  expect_s3_class(C, "treat")
  expect_identical(get_treat_type(C), "censoring")
  expect_true(is.numeric(C))
  expect_true(all(unclass(C) %in% c(0, 1)))
  expect_identical(attr(C, "treat.name"), "d$C")

  # The tag must survive subsetting: the `by` and risk-set paths rely on it
  expect_identical(get_treat_type(C[1:10]), "censoring")

  # Logical and 0/1 factors are accepted
  expect_identical(get_treat_type(.cens(as.logical(d$C))), "censoring")
  expect_identical(get_treat_type(.cens(factor(d$C))), "censoring")

  # Anything else errors at the call site
  expect_error(.cens(c(0, 1, 2)), "only the values")
  expect_error(.cens(factor(ifelse(d$C == 1L, "Yes", "No"))), "only the values")
})

test_that("both ways of requesting censoring weights agree", {
  d <- make_cens_data()

  X <- as.matrix(d[c("X1", "X2", "X3", "X4")])

  W_formula <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "glm")
  WF_tagged <- weightit.fit(X, treat = .cens(d$C), method = "glm")

  expect_equal(unname(W_formula$weights), unname(WF_tagged$weights),
               tolerance = eps)

  # An untagged 0/1 vector is still an ordinary binary treatment
  WF_binary <- weightit.fit(X, treat = d$C, method = "glm", estimand = "ATE")
  expect_false(isTRUE(all.equal(unname(WF_binary$weights),
                                unname(WF_tagged$weights))))
})

test_that("get_w_from_ps() computes censoring weights for a tagged indicator", {
  d <- make_cens_data()

  W <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "glm")

  w <- get_w_from_ps(W$ps, treat = .cens(d$C))

  # ps is the probability of BEING CENSORED
  expect_equal(unname(w), unname((1 - d$C) / (1 - W$ps)), tolerance = eps)
  expect_equal(unname(w), unname(W$weights), tolerance = eps)
  expect_true(all(w[d$C == 1L] == 0))

  # An untagged indicator is treated as binary, not silently as censoring
  expect_false(isTRUE(all.equal(unname(get_w_from_ps(W$ps, treat = d$C)),
                                unname(w))))
})

test_that("get_w_from_ps() ignores inapplicable arguments with a warning", {
  d <- make_cens_data()

  W <- weightit(.cens(C) ~ X1, data = d, method = "glm")
  C <- .cens(d$C)

  expect_warning(get_w_from_ps(W$ps, C, estimand = "ATT"), "ignored")
  expect_warning(get_w_from_ps(W$ps, C, focal = 1), "ignored")
  expect_warning(get_w_from_ps(W$ps, C, subclass = 5), "ignored")
  expect_warning(get_w_from_ps(W$ps, C, stabilize = TRUE), "ignored")

  # A propensity score matrix makes no sense here
  expect_error(get_w_from_ps(cbind(W$ps, 1 - W$ps), C), "numeric vector")
  expect_error(get_w_from_ps(W$ps[-1L], C), "same number of units")
})

test_that("a tagged indicator routes weightit.fit() to every .cens method", {
  skip_if_not_installed("rootSolve")

  d <- make_cens_data()

  X <- as.matrix(d[c("X1", "X2", "X3")])

  for (m in c("glm", "ebal", "cbps", "ipt")) {
    WF <- weightit.fit(X, treat = .cens(d$C), method = m)
    W <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = m)

    expect_equal(unname(WF$weights), unname(W$weights), tolerance = eps,
                 label = sprintf("weightit.fit vs weightit for method = \"%s\"", m))
  }
})

# ---- weightit(): weights, output shape, and options ------------------------

test_that("glm censoring weights equal 1/P(C = 0 | X) and zero out the censored", {
  d <- make_cens_data()

  expect_no_condition({
    W <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "glm")
  })

  # The propensity score is the probability of BEING CENSORED
  ps <- fitted(glm(C ~ X1 + X2 + X3 + X4, data = d, family = quasibinomial))
  expect_equal(unname(W$ps), unname(ps), tolerance = eps)

  expect_equal(unname(W$weights), unname((1 - d$C) / (1 - ps)), tolerance = eps)

  # Exactly 0, not merely small
  expect_true(all(W$weights[d$C == 1L] == 0))
  expect_true(all(W$weights[d$C == 0L] > 0))
  expect_false(anyNA(W$weights))

  # No estimand or focal applies
  expect_null(W$estimand)
  expect_null(W$focal)

  expect_output(print(W), "censoring")
  expect_no_condition(summary(W))
})

test_that("glm censoring weights are the binary ATT weights plus one", {
  d <- make_cens_data()

  W <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "glm")
  W_att <- weightit(C ~ X1 + X2 + X3 + X4, data = d, method = "glm",
                    estimand = "ATT", focal = 1)

  # Same model, so same propensity score; only the ps -> w map differs
  expect_equal(unname(W$ps), unname(W_att$ps), tolerance = eps)
  expect_equal(unname(W$weights), unname((1 - d$C) * (W_att$weights + 1)),
               tolerance = eps)
})

test_that("method = NULL gives complete-case weights", {
  d <- make_cens_data()

  W <- weightit(.cens(C) ~ X1 + X2, data = d, method = NULL)

  expect_equal(unname(W$weights), 1 - d$C, tolerance = eps)
})

test_that("a supplied ps is used as the probability of being censored", {
  d <- make_cens_data()

  ps <- plogis(-0.9 + 0.8 * d$X1 - 0.5 * d$X3)

  W <- weightit(.cens(C) ~ X1 + X3, data = d, ps = ps)

  expect_equal(unname(W$weights), unname((1 - d$C) / (1 - ps)), tolerance = eps)
})

test_that("options that do not apply to censoring are rejected or ignored", {
  d <- make_cens_data()

  # subclass is meaningless here
  expect_error(weightit(.cens(C) ~ X1, data = d, method = "glm", subclass = 5),
               "[Ss]ubclass")

  # estimand is ignored with a warning
  expect_warning(weightit(.cens(C) ~ X1, data = d, method = "glm",
                          estimand = "ATT"),
                 "ignored for censoring")

  # NAs in the indicator are an error for point treatments (there is no risk set)
  dna <- d
  is.na(dna$C[1L]) <- TRUE
  expect_error(weightit(.cens(C) ~ X1, data = dna, method = "glm"),
               "censoring indicator")

  # `get_w_from_ps()` recognizes the stored censoring treat and returns IPCW
  # weights rather than binary ATE weights
  W <- weightit(.cens(C) ~ X1, data = d, method = "glm")
  expect_equal(unname(get_w_from_ps(W$ps, W$treat)), unname(W$weights),
               tolerance = eps)

  # Downstream functions that cannot handle censoring say so
  expect_error(sbps(W, moderator = ~X5), "cannot be used with censoring")
})

test_that("trim() leaves the censored units at zero", {
  d <- make_cens_data()

  W <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "glm")
  Wt <- trim(W, at = .9)

  expect_true(all(Wt$weights[d$C == 1L] == 0))
})

# ---- Marginal (empty) censoring models -------------------------------------

# An empty right side (`.cens(C) ~ 1`) requests a censoring model that does not
# depend on the covariates, which must still route through the whole censoring
# infrastructure. This mirrors an empty treatment formula in `weightit()`.

test_that("an empty censoring formula gives a marginal censoring model", {
  d <- make_cens_data()

  W <- weightit(.cens(C) ~ 1, data = d, method = "glm")

  expect_identical(get_treat_type(W$treat), "censoring")
  expect_true(all(W$weights[d$C == 1L] == 0))

  # 1/P(C = 0), constant across the units still under observation
  expect_equal(unname(W$weights[d$C == 0L]),
               rep.int(1 / mean(d$C == 0L), sum(d$C == 0L)),
               tolerance = eps)

  # `ps` is the (constant) marginal probability of being censored
  expect_equal(unname(W$ps), rep.int(mean(d$C == 1L), nrow(d)),
               tolerance = eps)

  # There are no covariates to report, exactly as for an empty treatment formula
  expect_length(W$covs, 0L)
  expect_length(weightit(A ~ 1, data = d, method = "glm")$covs, 0L)

  # method = NULL and a supplied `ps` also take the censoring path
  expect_equal(unname(weightit(.cens(C) ~ 1, data = d, method = NULL)$weights),
               1 - d$C, tolerance = eps)
  expect_equal(unname(weightit(.cens(C) ~ 1, data = d,
                               ps = rep.int(mean(d$C), nrow(d)))$weights),
               unname(W$weights), tolerance = eps)

  # `weightit.fit()` accepts the equivalent zero-column `covs`
  W.fit <- weightit.fit(matrix(0, nrow = nrow(d), ncol = 0L),
                        treat = .cens(d$C), method = "glm")
  expect_equal(unname(W.fit$weights), unname(W$weights), tolerance = eps)
})

test_that("bart censoring weights zero out the censored and invert P(C = 0 | X)", {
  skip_if_not_installed("dbarts")

  # `weightit2bart.cens()` was the one `.cens` method with no test at all
  d <- make_cens_data()

  set.seed(123)
  W <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "bart",
                n.trees = 20, n.samples = 50, n.burn = 50, n.chains = 1,
                n.threads = 1)

  expect_identical(get_treat_type(W$treat), "censoring")
  expect_true(all(W$weights[d$C == 1L] == 0))
  expect_true(all(W$weights[d$C == 0L] > 0))

  # `ps` is the probability of BEING censored, and the weights invert its complement
  expect_true(all(W$ps > 0 & W$ps < 1))
  expect_equal(unname(W$weights[d$C == 0L]),
               unname(1 / (1 - W$ps[d$C == 0L])),
               tolerance = eps)

  # BART supplies no M-estimation parts, for censoring as for treatments
  expect_null(attr(W, "Mparts", exact = TRUE))
})

test_that("every method gives the same marginal censoring weights", {
  d <- make_cens_data()

  # With no covariates there is nothing for any method to model or balance, so
  # every method lands on the same 1/P(C = 0). How that is computed is covered in
  # test-empty_formula.R.
  W_glm <- weightit(.cens(C) ~ 1, data = d, method = "glm")

  for (m in c("ebal", "ipt", "cbps", "energy", "optweight", "cfd")) {
    if (!all(vapply(.weightit_methods[[m]]$packages_needed,
                    rlang::is_installed, logical(1L)))) {
      next
    }

    W <- weightit(.cens(C) ~ 1, data = d, method = m)

    expect_identical(as.character(W$method), m)
    expect_equal(unname(W$weights), unname(W_glm$weights), tolerance = eps)
  }
})

test_that("a marginal censoring model composes with by, s.weights, and stabilize", {
  skip_if_not_installed("rootSolve")

  d <- make_cens_data()

  # `by`: the marginal model is fit within each stratum
  W_by <- weightit(.cens(C) ~ 1, data = d, method = "glm", by = ~G)

  for (g in levels(d$G)) {
    in.g <- which(d$G == g & d$C == 0L)
    expect_equal(unname(W_by$weights[in.g]),
                 rep.int(1 / mean(d$C[d$G == g] == 0L), length(in.g)),
                 tolerance = eps)
  }

  expect_M_parts_okay(W_by, tolerance = eps)

  # `s.weights`: the marginal rate is the sampling-weighted one
  W_sw <- weightit(.cens(C) ~ 1, data = d, method = "glm", s.weights = "SW")
  expect_equal(unname(W_sw$weights[d$C == 0L]),
               rep.int(sum(d$SW) / sum(d$SW * (d$C == 0L)), sum(d$C == 0L)),
               tolerance = eps)

  # `stabilize`: the numerator is the same marginal model, so every nonzero
  # weight is exactly 1
  W_st <- weightit(.cens(C) ~ 1, data = d, method = "glm", stabilize = TRUE)
  expect_equal(unname(W_st$weights), 1 - d$C, tolerance = eps)

  expect_M_parts_okay(W_st, tolerance = eps)
})

test_that("a marginal censoring model supports M-estimation", {
  skip_if_not_installed("rootSolve")

  d <- make_cens_data()

  for (m in c("glm", "ebal", "ipt")) {
    W <- weightit(.cens(C) ~ 1, data = d, method = m)
    expect_M_parts_okay(W, tolerance = eps)
  }

  # NAs are still tolerated in the censored units, since their weight is 0
  d$Y <- d$Y_B
  is.na(d$Y[d$C == 1L]) <- TRUE

  W <- weightit(.cens(C) ~ 1, data = d, method = "glm")
  fit <- glm_weightit(Y ~ X1, data = d, weightit = W, family = binomial)

  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
})

# ---- Degenerate risk sets --------------------------------------------------

test_that("no censored units gives all weights of 1", {
  d <- make_cens_data()
  d$C <- 0L

  for (m in c("glm", "ebal", "cbps", "ipt", "energy", "optweight")) {
    W <- expect_no_error(weightit(.cens(C) ~ X1 + X2, data = d, method = m))
    expect_true(all(W$weights == 1),
                label = sprintf("all weights 1 for method = \"%s\"", m))
  }
})

test_that("all units censored warns and gives all weights of 0", {
  d <- make_cens_data()
  d$C <- 1L

  expect_warning(W <- weightit(.cens(C) ~ X1 + X2, data = d, method = "glm"),
                 "censored")
  expect_true(all(W$weights == 0))
})

test_that("both degenerate risk sets are caught before an empty model is fit", {
  d <- make_cens_data()

  d0 <- d
  d0$C <- 0L
  expect_true(all(weightit(.cens(C) ~ 1, data = d0, method = "glm")$weights == 1))

  d1 <- d
  d1$C <- 1L
  expect_warning(W1 <- weightit(.cens(C) ~ 1, data = d1, method = "glm"),
                 "censored")
  expect_true(all(W1$weights == 0))
})

test_that("very few uncensored units still works", {
  # ~5% uncensored: the case that motivates estimating weights for one group only
  d <- make_cens_data(seed = 5L, p = 3)

  skip_if(mean(d$C == 0) > .15)

  for (m in c("glm", "ebal", "cbps", "ipt")) {
    W <- expect_no_error(weightit(.cens(C) ~ X1 + X2, data = d, method = m))
    expect_true(all(W$weights[d$C == 1L] == 0))
    expect_true(all(W$weights[d$C == 0L] > 0))
  }
})

# ---- Per-method behavior ---------------------------------------------------

test_that("balance-based methods hit the full-sample target exactly", {
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("osqp")
  skip_if_not_installed("optweight")

  d <- make_cens_data()

  # These methods solve the balance conditions directly, so the weighted
  # uncensored means equal the full at-risk sample means to solver precision.
  for (m in c("ebal", "cbps", "ipt", "optweight")) {
    W <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = m)

    expect_lt(target_diff(W, d), 1e-5)
    expect_true(all(W$weights[d$C == 1L] == 0))
  }

  # `glm` maximizes a likelihood rather than solving the balance conditions, so it
  # reduces imbalance substantially but not to 0.
  Wg <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "glm")
  W0 <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = NULL)

  expect_lt(target_diff(Wg, d), target_diff(W0, d) / 5)
})

test_that("the balance-based methods put the weights on the 1/P(C = 0 | X) scale", {
  skip_if_not_installed("rootSolve")

  d <- make_cens_data()

  # The intercept row of the moment condition forces sum(w) == n
  for (m in c("ebal", "cbps", "ipt")) {
    W <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = m)
    expect_equal(sum(W$weights), nrow(d), tolerance = eps)
  }

  # optweight normalizes to a mean of 1 among the uncensored instead
  skip_if_not_installed("optweight")
  Wo <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "optweight")
  expect_equal(mean(Wo$weights[d$C == 0L]), 1, tolerance = eps)
})

test_that("ebal censoring weights equal ipt with link = 'clog'", {
  skip_if_not_installed("rootSolve")

  # Entropy balancing solves the same problem as inverse probability tilting under
  # the complementary log link: 1/(1 - e) == exp(X'b). Two entirely different
  # algorithms, so this is a strong cross-check on both.
  d <- make_cens_data()

  W_ebal <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "ebal")
  W_ipt <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "ipt",
                    link = "clog")

  expect_equal(unname(W_ebal$weights), unname(W_ipt$weights), tolerance = 1e-6)
})

test_that("cbps censoring weights equal ipt when just-identified", {
  skip_if_not_installed("rootSolve")

  # Both reduce to sum(SW * [(1 - C)/(1 - e) - 1] * X) == 0
  d <- make_cens_data()

  W_cbps <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "cbps")
  W_ipt <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "ipt")

  expect_equal(unname(W_cbps$weights), unname(W_ipt$weights), tolerance = 1e-6)
  expect_equal(unname(W_cbps$ps), unname(W_ipt$ps), tolerance = 1e-6)
})

test_that("cfd with kernel = 'energy' equals energy", {
  skip_if_not_installed("osqp")

  d <- make_cens_data()

  W_e <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "energy")
  W_c <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "cfd",
                  kernel = "energy")

  expect_equal(unname(W_e$weights), unname(W_c$weights), tolerance = 1e-4)
})

test_that("optweight with norm = 'entropy' equals ebal up to scale", {
  skip_if_not_installed("optweight")

  d <- make_cens_data()

  W_ebal <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "ebal")
  W_ow <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "optweight",
                   norm = "entropy")

  u <- d$C == 0L

  expect_equal(W_ebal$weights[u] / mean(W_ebal$weights[u]),
               W_ow$weights[u] / mean(W_ow$weights[u]),
               tolerance = 1e-4)
})

test_that("energy and cfd reduce the distance to the full at-risk sample", {
  skip_if_not_installed("osqp")

  d <- make_cens_data()

  X <- scale(as.matrix(d[c("X1", "X2", "X3")]))
  D <- as.matrix(dist(X))
  n <- nrow(d)
  u <- which(d$C == 0L)

  edist <- function(w) {
    m <- rep.int(0, n)
    m[u] <- w[u]
    m <- m / sum(m)
    r <- rep.int(1 / n, n)

    2 * drop(t(m) %*% D %*% r) - drop(t(m) %*% D %*% m) - drop(t(r) %*% D %*% r)
  }

  W_e <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "energy")

  # The energy method minimizes exactly this criterion, so it must beat both the
  # unweighted uncensored sample and any other method's weights.
  expect_lt(edist(W_e$weights), edist(rep.int(1, n)))
  expect_lte(edist(W_e$weights),
             edist(weightit(.cens(C) ~ X1 + X2 + X3, data = d,
                            method = "glm")$weights))

  # min.w must not leak into the censored units' zeros
  expect_true(all(W_e$weights[d$C == 1L] == 0))
})

test_that("energy with moments = 1 and tols = 0 also hits the target", {
  skip_if_not_installed("osqp")

  d <- make_cens_data()

  W <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "energy",
                moments = 1)

  expect_lt(target_diff(W, d), 1e-5)
})

test_that("s.weights are respected", {
  skip_if_not_installed("rootSolve")

  d <- make_cens_data()

  for (m in c("ebal", "cbps", "ipt")) {
    W <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = m,
                  s.weights = d$SW)

    # target_diff() weights the target by s.weights too
    expect_lt(target_diff(W, d), 1e-5)
    expect_M_parts_okay(W, tolerance = eps)
  }
})

# ---- Balance-tuned methods use target balance, not a between-group proxy ----
#
# `gbm` and `super` choose tuning parameters by minimizing a balance criterion.
# For censoring the criterion must compare the weighted units still under
# observation against the FULL at-risk sample, which is what the weights target.
# cobalt expresses this as a "target" init: omit `treat` from bal.init() and pass
# the full sample as `x`, encoding the uncensored subset with zero weights.

# The target init and the augmented/stacked-dataset construction are two ways of
# writing the same quantity; this helper computes the second.
augmented_balance <- function(covs, C, w, stat = "smd.mean", s.weights = NULL) {
  n <- nrow(covs)
  u <- which(C == 0)

  if (is_null(s.weights)) s.weights <- rep.int(1, n)

  init <- cobalt::bal.init(rbind(covs[u, , drop = FALSE], covs),
                           treat = c(rep.int(0L, length(u)), rep.int(1L, n)),
                           stat = stat, estimand = "ATT", focal = 1,
                           s.weights = c(s.weights[u], s.weights))

  cobalt::bal.compute(init, weights = c(w[u], rep.int(1, n)))
}

test_that("a target init equals the augmented-dataset construction", {
  skip_if_not_installed("cobalt")

  d <- make_cens_data()

  covs <- as.matrix(d[c("X1", "X2", "X3", "X4")])
  w <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "glm")$weights

  # Exact for the stats that depend only on weighted moments/ECDFs
  for (st in c("smd.mean", "smd.max", "smd.rms", "ks.mean", "ks.max",
               "mahalanobis")) {
    init <- cobalt::bal.init(covs, stat = st)

    expect_equal(cobalt::bal.compute(init, weights = w),
                 augmented_balance(covs, d$C, w, stat = st),
                 label = sprintf("target init for stat = \"%s\"", st))
  }

  # And with sampling weights
  init_sw <- cobalt::bal.init(covs, stat = "smd.mean", s.weights = d$SW)
  expect_equal(cobalt::bal.compute(init_sw, weights = w),
               augmented_balance(covs, d$C, w, s.weights = d$SW))

  # `ovl.*` and `energy.dist` are NOT expected to agree: they depend on the row
  # multiset, and the stacked matrix has n0 + n rows where the target init has 2n.
  # Use the native target init for those, never the stacked trick.
  init_ovl <- cobalt::bal.init(covs, stat = "ovl.mean")
  expect_false(isTRUE(all.equal(cobalt::bal.compute(init_ovl, weights = w),
                                augmented_balance(covs, d$C, w, stat = "ovl.mean"))))
})

test_that("gbm reports the target criterion, not the between-group proxy", {
  skip_if_not_installed("gbm")
  skip_if_not_installed("cobalt")

  d <- make_cens_data(seed = 42L)

  covs <- as.matrix(d[c("X1", "X2", "X3", "X4")])

  W <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "gbm",
                criterion = "smd.mean", n.trees = 200)

  reported <- W$info$tree.val$smd.mean[W$info$tree.val$tree == W$info$best.tree]

  # What the criterion should be: weighted uncensored vs. the full at-risk sample
  target <- cobalt::bal.compute(cobalt::bal.init(covs, stat = "smd.mean"),
                                weights = W$weights)

  expect_equal(reported, target, tolerance = eps)
  expect_equal(reported, augmented_balance(covs, d$C, W$weights), tolerance = eps)

  # What it must NOT be: the old proxy, evaluated on ATT weights against the
  # censored group
  proxy <- cobalt::bal.compute(
    cobalt::bal.init(covs, treat = d$C, stat = "smd.mean",
                     estimand = "ATT", focal = 1),
    weights = ifelse(d$C == 1L, 1, W$weights - 1))

  expect_false(isTRUE(all.equal(reported, proxy)))
})

test_that("the target criterion selects a different tree than the old proxy", {
  skip_if_not_installed("gbm")

  # With a shrinkage large enough that the criterion has an interior minimum, the
  # two criteria order candidate trees differently. Without this the criterion is
  # still improving at `n.trees` and both pick the last tree, so the test would
  # pass even if the proxy were restored.
  set.seed(3)
  d <- test_data
  d$C <- rbinom(nrow(d), 1L,
                plogis(-1.2 + 1.2 * d$X1 - 0.8 * d$X3 + 0.6 * (d$X6 == "A")))

  args <- list(data = d, method = "gbm", criterion = "smd.mean",
               n.trees = 300, shrinkage = 0.3, interaction.depth = 4)

  W_target <- do.call(weightit, c(list(.cens(C) ~ X1 + X2 + X3 + X6), args))
  W_proxy <- do.call(weightit, c(list(C ~ X1 + X2 + X3 + X6),
                                 args, list(estimand = "ATT", focal = 1)))

  expect_false(W_target$info$best.tree == W_proxy$info$best.tree)
})

test_that("gbm censoring honors the target criterion set and cv", {
  skip_if_not_installed("gbm")

  d <- make_cens_data()

  # available.stats("target") is a strict subset of "binary": `r2` and friends are
  # not defined against a target and are now rejected
  expect_error(weightit(.cens(C) ~ X1 + X2, data = d, method = "gbm",
                        criterion = "r2", n.trees = 50),
               "criterion")

  # `cv{#}` selects trees by cross-validation rather than balance, but must still
  # return censoring weights
  W <- weightit(.cens(C) ~ X1 + X2, data = d, method = "gbm",
                criterion = "cv3", n.trees = 60)

  expect_true(all(W$weights[d$C == 1L] == 0))
  expect_true(all(W$weights[d$C == 0L] > 0))
  expect_equal(unname(W$weights), unname((1 - d$C) / (1 - W$ps)), tolerance = eps)
})

test_that("super's balance criterion is target-based for censoring", {
  skip_if_not_installed("SuperLearner")
  skip_if_not_installed("cobalt")

  d <- make_cens_data(seed = 42L)

  covs <- as.matrix(d[c("X1", "X2", "X3", "X4")])

  W <- suppressWarnings(
    weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "super",
             SL.library = c("SL.glm", "SL.mean"),
             SL.method = "method.balance", criterion = "smd.mean")
  )

  expect_true(all(W$weights[d$C == 1L] == 0))
  expect_true(all(W$weights[d$C == 0L] > 0))

  # The realized balance agrees between the two constructions
  expect_equal(cobalt::bal.compute(cobalt::bal.init(covs, stat = "smd.mean"),
                                   weights = W$weights),
               augmented_balance(covs, d$C, W$weights),
               tolerance = eps)

  # A library member that cannot model censoring at all should score worse
  expect_gt(W$info$cvRisk[["SL.mean_All"]], W$info$cvRisk[["SL.glm_All"]])

  expect_error(suppressWarnings(
    weightit(.cens(C) ~ X1, data = d, method = "super",
             SL.library = c("SL.glm", "SL.mean"),
             SL.method = "method.balance", criterion = "r2")),
    "criterion")
})

# ---- M-estimation ----------------------------------------------------------

test_that("M-estimation parts are internally consistent", {
  skip_if_not_installed("rootSolve")

  d <- make_cens_data()

  # `expect_M_parts_okay()` re-solves colSums(psi_treat) == 0 and checks that wfun
  # reproduces the weights, INCLUDING the zeros -- so it directly catches a wfun
  # that returns 1 rather than 0 for the censored units.
  for (m in c("glm", "ebal", "cbps", "ipt")) {
    W <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = m)
    expect_M_parts_okay(W, tolerance = eps)
  }

  # glm with non-default links
  for (lk in c("probit", "cloglog", "br.logit"[rlang::is_installed("brglm2")])) {
    W <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "glm", link = lk)
    expect_M_parts_okay(W, tolerance = eps)
  }
})

test_that("methods without M-estimation return no Mparts", {
  skip_if_not_installed("osqp")

  d <- make_cens_data()

  for (m in c("energy", "cfd")) {
    W <- weightit(.cens(C) ~ X1 + X2, data = d, method = m)
    expect_null(attr(W, "Mparts", exact = TRUE))
    expect_null(attr(W, "Mparts.list", exact = TRUE))
  }

  # cbps with over = TRUE is over-identified, so no Mparts (as for binary)
  skip_if_not_installed("rootSolve")
  W <- weightit(.cens(C) ~ X1 + X2, data = d, method = "cbps", over = TRUE)
  expect_null(attr(W, "Mparts", exact = TRUE))
})

test_that("censoring composes with by", {
  skip_if_not_installed("rootSolve")

  d <- make_cens_data()

  W <- weightit(.cens(C) ~ X1 + X2, data = d, method = "glm", by = ~G)

  expect_null(attr(W, "Mparts", exact = TRUE))
  expect_false(is_null(attr(W, "Mparts.list", exact = TRUE)))
  expect_length(attr(W, "Mparts.list"), nlevels(d$G))

  expect_M_parts_okay(W, tolerance = eps)

  # Estimating within strata is equivalent to interacting the stratifier with all
  # the covariates, so the asymptotic variances must agree.
  W_int <- weightit(.cens(C) ~ G * (X1 + X2), data = d, method = "glm")

  d$Y <- d$Y_B
  is.na(d$Y[d$C == 1L]) <- TRUE

  f_by <- glm_weightit(Y ~ X1, data = d, weightit = W, vcov = "asympt",
                       family = binomial)
  f_int <- glm_weightit(Y ~ X1, data = d, weightit = W_int, vcov = "asympt",
                        family = binomial)

  expect_equal(unname(coef(f_by)), unname(coef(f_int)), tolerance = eps)
  expect_equal(unname(vcov(f_by)), unname(vcov(f_int)), tolerance = 1e-5)
})

test_that("a by stratum with no uncensored units errors", {
  d <- make_cens_data()
  d$C[d$G == "g1"] <- 1L

  expect_error(weightit(.cens(C) ~ X1, data = d, method = "glm", by = ~G),
               "[Aa]ll units are censored")
})

test_that("censoring weights can be stabilized", {
  skip_if_not_installed("rootSolve")

  d <- make_cens_data()

  W <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "glm",
                stabilize = TRUE)

  # The numerator is a marginal censoring model, so w = P(C = 0)/P(C = 0 | X);
  # crucially no 0/0 for the censored units.
  expect_true(all(is.finite(W$weights)))
  expect_true(all(W$weights[d$C == 1L] == 0))

  W0 <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "glm")
  expect_equal(unname(W$weights), unname(W0$weights * mean(d$C == 0L)),
               tolerance = 1e-3)

  expect_M_parts_okay(W, tolerance = eps)

  Wn <- weightit(.cens(C) ~ X1 + X2 + X3, data = d, method = "glm",
                 stabilize = ~X1)

  # w = P(C = 0 | X1)/P(C = 0 | X) for the uncensored, exactly 0 for the censored
  den <- fitted(glm(C ~ X1 + X2 + X3, data = d, family = quasibinomial))
  num <- fitted(glm(C ~ X1, data = d, family = quasibinomial))

  expect_equal(unname(Wn$weights),
               unname(ifelse(d$C == 0L, (1 - num) / (1 - den), 0)),
               tolerance = eps)
  expect_true(all(is.finite(Wn$weights)))

  # The censoring numerator's weights are 0 for the censored units, so inverting
  # them would give `Inf` and then `0 * Inf = NaN`. They are held at a factor of 1
  # instead, which leaves the weights alone (the denominator already zeroes them)
  # and keeps `wfun` reproducing them exactly, including the zeros.
  mp <- attr(Wn, "Mparts.list")
  expect_length(mp, 2L)

  wfun_num <- mp[[2L]]$wfun(mp[[2L]]$btreat, mp[[2L]]$Xtreat, mp[[2L]]$A)
  expect_equal(unname(wfun_num[d$C == 0L]), unname((1 - num)[d$C == 0L]),
               tolerance = eps)
  expect_true(all(wfun_num[d$C == 1L] == 1))

  # (`wfun` multiplies by the reciprocal where the weights themselves divide, so
  # these agree to floating-point precision rather than bit for bit)
  wfun_prod <- Reduce("*", lapply(mp, function(m) m$wfun(m$btreat, m$Xtreat, m$A)),
                      init = 1)
  expect_equal(unname(wfun_prod), unname(Wn$weights), tolerance = eps)
  expect_identical(unname(wfun_prod == 0), unname(Wn$weights == 0))

  expect_M_parts_okay(Wn, tolerance = eps)
})

# ---- Outcome models tolerate NAs in zero-weight units ----------------------

test_that("glm_weightit tolerates NAs in censored units", {
  d <- make_cens_data()

  cens <- d$C == 1L

  # After censoring, the outcome and later covariates are unobserved
  d$Y <- d$Y_C
  is.na(d$Y[cens]) <- TRUE
  d$Xpost <- d$X2
  is.na(d$Xpost[cens]) <- TRUE

  W <- weightit(.cens(C) ~ X1 + X3, data = d, method = "glm")

  expect_no_error({
    fit <- glm_weightit(Y ~ X1 + Xpost, data = d, weightit = W, vcov = "asympt")
  })

  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))

  # Zero-weight units are excluded from n and the residual df, but no rows are
  # dropped, so the fit stays aligned with the full-length M-estimation parts.
  expect_identical(stats::nobs(fit), sum(!cens))
  expect_identical(fit$df.residual, sum(!cens) - length(coef(fit)))
  expect_null(fit$na.action)
  expect_identical(nrow(fit$model), nrow(d))

  # The returned model frame reports the data as supplied
  expect_true(anyNA(fit$model$Y))

  # Identical to fitting on the uncensored units with the same weights
  cc <- glm(Y ~ X1 + Xpost, data = d[!cens, ], weights = W$weights[!cens])
  expect_equal(unname(coef(fit)), unname(coef(cc)), tolerance = eps)

  # A missing value in a NONZERO-weight unit must still error
  dbad <- d
  is.na(dbad$Y[which(!cens)[1L]]) <- TRUE
  expect_error(glm_weightit(Y ~ X1 + Xpost, data = dbad, weightit = W,
                            vcov = "asympt"),
               "[Mm]issing")
})

test_that("other model types also tolerate NAs in censored units", {
  d <- make_cens_data()

  cens <- d$C == 1L
  d$Y <- d$Y_C
  is.na(d$Y[cens]) <- TRUE
  d$Yb <- factor(d$Y_B)
  is.na(d$Yb[cens]) <- TRUE
  d$Yo <- d$Y_O
  is.na(d$Yo[cens]) <- TRUE

  W <- weightit(.cens(C) ~ X1 + X3, data = d, method = "glm")

  expect_no_error({
    f_lm <- lm_weightit(Y ~ X1, data = d, weightit = W, vcov = "asympt")
  })
  expect_identical(stats::nobs(f_lm), sum(!cens))

  expect_no_error({
    f_mn <- multinom_weightit(Yb ~ X1, data = d, weightit = W, vcov = "asympt")
  })
  expect_identical(stats::nobs(f_mn), sum(!cens))

  expect_no_error({
    f_or <- ordinal_weightit(Yo ~ X1, data = d, weightit = W, vcov = "asympt")
  })
  expect_identical(stats::nobs(f_or), sum(!cens))
})

# `survival::coxph.fit()` refuses weights of 0, so `.coxph_weightit()` fits on the
# units with a positive weight and scatters the row-indexed components back. A unit
# with a weight of 0 contributes nothing to the weighted partial likelihood -- not
# its own term, and not to any risk set denominator -- so this is exact rather than
# an approximation.
test_that("coxph_weightit tolerates zero weights and NA event times", {
  skip_if_not_installed("survival")

  d <- make_cens_data()

  cens <- d$C == 1L

  cutoff <- quantile(d$Y_S, .8)
  d$event <- as.numeric(d$Y_S < cutoff)
  d$time <- pmin(d$Y_S, cutoff)

  # After censoring, the event time is never ascertained
  is.na(d$time[cens]) <- TRUE
  is.na(d$event[cens]) <- TRUE

  W <- weightit(.cens(C) ~ X1 + X3, data = d, method = "glm")

  expect_no_error({
    fit <- coxph_weightit(survival::Surv(time, event) ~ A + X2,
                          data = d, weightit = W, vcov = "asympt", x = TRUE)
  })

  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
  expect_null(fit$na.action)

  # No rows are dropped, so the fit stays aligned with the full-length
  # M-estimation parts; only the units that contribute are counted as events.
  expect_identical(fit$n, nrow(d))
  expect_length(fit$residuals, nrow(d))
  expect_length(fit$linear.predictors, nrow(d))
  expect_true(all(is.finite(fit$linear.predictors)))
  expect_identical(fit$nevent, sum(d$event[!cens]))

  # Exactly equal to fitting on the uncensored units with the same weights
  ref <- survival::coxph(survival::Surv(time, event) ~ A + X2,
                         data = d[!cens, ], weights = W$weights[!cens])

  expect_equal(unname(coef(fit)), unname(coef(ref)), tolerance = eps)
  expect_equal(fit$loglik, ref$loglik, tolerance = eps)

  # The estimating equation is solved, and the censored units contribute nothing
  w <- W$weights * W$s.weights
  psi <- fit$psi(coef(fit), fit$x, fit$y, w)

  expect_equal(colSums(psi), rep_with(0, coef(fit)), tolerance = eps,
               ignore_attr = TRUE)
  expect_true(all(psi[cens, ] == 0))

  # `.compute_vcov()` evaluates `psi` with the estimated weights set to 1 to build
  # the cross-derivative block. The censored units must stay out of the risk sets
  # there too, or they perturb the psi of the units that do contribute.
  psi1 <- fit$psi(coef(fit), fit$x, fit$y, W$s.weights)

  expect_true(all(psi1[cens, ] == 0))

  # Masking internally is equivalent to never having supplied the censored units:
  # a closure built on the uncensored subset alone gives the same values.
  psi1_sub <- .get_coxph_psi(list(y = fit$y[!cens]))(
    coef(fit), fit$x[!cens, , drop = FALSE], fit$y[!cens], W$s.weights[!cens])

  expect_equal(psi1[!cens, ], psi1_sub, tolerance = eps, ignore_attr = TRUE)

  # A missing value in a NONZERO-weight unit must still error
  dbad <- d
  is.na(dbad$time[which(!cens)[1L]]) <- TRUE
  expect_error(coxph_weightit(survival::Surv(time, event) ~ A + X2, data = dbad,
                              weightit = W, vcov = "asympt"),
               "[Mm]issing")

  # All weights 0 leaves nothing to fit
  W0 <- W
  W0$weights[] <- 0
  expect_error(coxph_weightit(survival::Surv(time, event) ~ A + X2, data = d,
                              weightit = W0),
               "[Aa]ll weights are 0")
})

test_that("coxph_weightit variance types work with censoring weights", {
  skip_if_not_installed("survival")

  d <- make_cens_data()

  cutoff <- quantile(d$Y_S, .8)
  d$event <- as.numeric(d$Y_S < cutoff)
  d$time <- pmin(d$Y_S, cutoff)
  is.na(d$time[d$C == 1L]) <- TRUE
  is.na(d$event[d$C == 1L]) <- TRUE

  W <- weightit(.cens(C) ~ X1 + X3, data = d, method = "glm")

  f_a <- coxph_weightit(survival::Surv(time, event) ~ A, data = d,
                        weightit = W, vcov = "asympt")
  f_hc <- coxph_weightit(survival::Surv(time, event) ~ A, data = d,
                         weightit = W, vcov = "HC0")

  expect_true(all(is.finite(sqrt(diag(vcov(f_a))))))
  expect_true(all(is.finite(sqrt(diag(vcov(f_hc))))))

  # Same point estimate regardless of how the variance is computed
  expect_identical(coef(f_hc), coef(f_a))

  # `vcov = "BS"`/`"FWB"` re-evaluate the weightit call, which contains `.cens()`
  set.seed(1)
  f_bs <- coxph_weightit(survival::Surv(time, event) ~ A, data = d,
                         weightit = W, vcov = "BS", R = 10)
  expect_true(all(is.finite(sqrt(diag(vcov(f_bs))))))

  set.seed(1)
  f_fwb <- coxph_weightit(survival::Surv(time, event) ~ A, data = d,
                          weightit = W, vcov = "FWB", R = 10)
  expect_true(all(is.finite(sqrt(diag(vcov(f_fwb))))))

  # The variance comes from `psi`, which masks the zero-weight units, rather than
  # from `residuals.coxph()`, which is undefined on a fit whose row-indexed
  # components were scattered back from a subset fit
  expect_null(f_a$gradient)
})

test_that("HC0 and bootstrap variances also work", {
  d <- make_cens_data()

  d$Y <- d$Y_C
  is.na(d$Y[d$C == 1L]) <- TRUE

  W <- weightit(.cens(C) ~ X1 + X3, data = d, method = "glm")

  f_hc <- glm_weightit(Y ~ X1, data = d, weightit = W, vcov = "HC0")
  expect_true(all(is.finite(sqrt(diag(vcov(f_hc))))))

  # `vcov = "BS"` re-evaluates the weightit call, which contains `.cens()`
  set.seed(1)
  f_bs <- glm_weightit(Y ~ X1, data = d, weightit = W, vcov = "BS", R = 10)
  expect_true(all(is.finite(sqrt(diag(vcov(f_bs))))))
})

# ---- weightitMSM() --------------------------------------------------------

make_msm_cens_data <- function(seed = 1234L) {
  set.seed(seed)
  d <- msmdata
  d$C_2 <- rbinom(nrow(d), 1L, prob = plogis(-2 + 0.2 * d$X1_1 - 0.3 * d$A_1))
  is.na(d[d$C_2 == 1L, c("X1_2", "X2_2", "A_3")]) <- TRUE
  d
}

msm_cens_formulas <- list(A_1 ~ X1_0 + X2_0,
                          A_2 ~ X1_1 + X2_1 + A_1,
                          .cens(C_2) ~ X1_1 + X2_1 + A_1,
                          A_3 ~ X1_2 + X2_2 + A_2)

test_that("weightitMSM separates censoring from treatment models", {
  d <- make_msm_cens_data()
  cens <- d$C_2 == 1L

  expect_no_condition({
    W <- weightitMSM(msm_cens_formulas, data = d, method = "glm",
                     include.obj = TRUE)
  })

  # treat.list and covs.list describe treatments only
  expect_length(W$treat.list, 3L)
  expect_named(W$treat.list, c("A_1", "A_2", "A_3"))
  expect_length(W$cens.list, 1L)
  expect_named(W$cens.list, "C_2")
  expect_equal(unname(W$cens.time), 3L)
  expect_identical(colnames(W$at.risk), c("A_1", "A_2", "C_2", "A_3"))

  # formula.list round-trips with the marker intact
  expect_length(W$formula.list, 4L)

  expect_true(all(W$weights[cens] == 0))
  expect_true(all(W$weights[!cens] > 0))
  expect_false(anyNA(W$weights))

  expect_output(print(W), "censoring")
  expect_no_condition(summary(W))

  # The censoring models get their own covariate listing, as the treatment models do
  expect_output(print(W), "censoring covariates")
  expect_output(print(W), "C_2: X1_1, X2_1, A_1")

  # ...and none appears when there is no censoring model (`d` has NAs in `A_3` for
  # the censored units, so use the untouched data here)
  W_nocens <- weightitMSM(msm_cens_formulas[-3L], data = msmdata, method = "glm")
  expect_output(print(W_nocens), "covariates")
  expect_failure(expect_output(print(W_nocens), "censoring covariates"))
})

test_that("each model is fit only among the units under observation", {
  # This is the regression test for the risk-set direction: with the comparison
  # inverted, the post-censoring model would be fit on the CENSORED units instead.
  d <- make_msm_cens_data()
  cens <- d$C_2 == 1L

  W <- weightitMSM(msm_cens_formulas, data = d, method = "glm",
                   include.obj = TRUE)

  expect_identical(stats::nobs(W$obj$A_3), sum(!cens))

  # The censoring model itself uses everyone still at risk, censored or not
  expect_identical(stats::nobs(W$obj$C_2), nrow(d))
})

test_that("weightitMSM weights equal the hand-computed product", {
  d <- make_msm_cens_data()
  cens <- d$C_2 == 1L

  W <- weightitMSM(msm_cens_formulas, data = d, method = "glm")

  w1 <- weightit(A_1 ~ X1_0 + X2_0, data = d, method = "glm")$weights
  w2 <- weightit(A_2 ~ X1_1 + X2_1 + A_1, data = d, method = "glm")$weights

  psc <- fitted(glm(C_2 ~ X1_1 + X2_1 + A_1, data = d, family = quasibinomial))
  cw <- ifelse(cens, 0, 1 / (1 - psc))

  w3 <- rep.int(1, nrow(d))
  w3[!cens] <- weightit(A_3 ~ X1_2 + X2_2 + A_2, data = d[!cens, ],
                        method = "glm")$weights

  expect_equal(unname(W$weights), unname(w1 * w2 * cw * w3), tolerance = eps)
})

test_that("weightitMSM censoring supports M-estimation", {
  skip_if_not_installed("rootSolve")

  d <- make_msm_cens_data()

  W <- weightitMSM(msm_cens_formulas, data = d, method = "glm")

  expect_length(attr(W, "Mparts.list"), 4L)
  expect_M_parts_okay(W, tolerance = eps)
})

test_that("weightitMSM censoring works with glm_weightit and NA outcomes", {
  d <- make_msm_cens_data()
  cens <- d$C_2 == 1L

  d$Y <- d$Y_B
  is.na(d$Y[cens]) <- TRUE

  W <- weightitMSM(msm_cens_formulas, data = d, method = "glm")

  expect_no_error({
    fit <- glm_weightit(Y ~ A_1 * A_2 * A_3, data = d, weightit = W,
                        vcov = "asympt", family = binomial)
  })

  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
  expect_identical(stats::nobs(fit), sum(!cens))
})

test_that("weightitMSM censoring works with coxph_weightit", {
  skip_if_not_installed("survival")

  d <- make_msm_cens_data()
  cens <- d$C_2 == 1L

  # `msmdata` has no event time, so build one
  set.seed(77)
  Y_S <- rexp(nrow(d), rate = exp(-1 + 0.3 * d$A_1 + 0.3 * d$A_2))
  cutoff <- quantile(Y_S, .8)
  d$event <- as.numeric(Y_S < cutoff)
  d$time <- pmin(Y_S, cutoff)
  is.na(d$time[cens]) <- TRUE
  is.na(d$event[cens]) <- TRUE

  W <- weightitMSM(msm_cens_formulas, data = d, method = "glm")

  expect_no_error({
    fit <- coxph_weightit(survival::Surv(time, event) ~ A_1 + A_2, data = d,
                          weightit = W, vcov = "asympt", x = TRUE)
  })

  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
  expect_identical(fit$n, nrow(d))

  w <- W$weights * W$s.weights
  psi <- fit$psi(coef(fit), fit$x, fit$y, w)

  expect_true(all(psi[cens, ] == 0))
  expect_equal(colSums(psi), rep_with(0, coef(fit)), tolerance = eps,
               ignore_attr = TRUE)
})

test_that("multiple censoring time points give nested risk sets", {
  set.seed(99)
  d <- msmdata
  n <- nrow(d)

  d$C_1 <- rbinom(n, 1L, plogis(-2.2 + 0.15 * d$X1_0))
  c1 <- d$C_1 == 1L

  d$C_2 <- rbinom(n, 1L, plogis(-2 + 0.2 * d$X1_1))
  is.na(d$C_2[c1]) <- TRUE
  c2 <- !is.na(d$C_2) & d$C_2 == 1L

  is.na(d[c1, c("X1_1", "X2_1", "A_2", "X1_2", "X2_2", "A_3")]) <- TRUE
  is.na(d[c2, c("X1_2", "X2_2", "A_3")]) <- TRUE

  ever <- c1 | c2

  W <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                        .cens(C_1) ~ X1_0 + X2_0,
                        A_2 ~ X1_1 + X2_1 + A_1,
                        .cens(C_2) ~ X1_1 + X2_1 + A_1,
                        A_3 ~ X1_2 + X2_2 + A_2),
                   data = d, method = "glm", include.obj = TRUE)

  expect_equal(unname(W$cens.time), c(2L, 4L))
  expect_length(attr(W, "Mparts.list"), 5L)

  # Once censored, always censored
  expect_true(all(W$weights[ever] == 0))
  expect_true(all(W$weights[!ever] > 0))

  # Risk sets shrink monotonically
  expect_true(all(W$at.risk[, "A_3"] <= W$at.risk[, "A_2"]))
  expect_true(all(W$at.risk[, "A_2"] <= W$at.risk[, "A_1"]))

  expect_identical(stats::nobs(W$obj$C_1), n)
  expect_identical(stats::nobs(W$obj$C_2), sum(!c1))
  expect_identical(stats::nobs(W$obj$A_3), sum(!ever))

  skip_if_not_installed("rootSolve")
  expect_M_parts_okay(W, tolerance = eps)
})

test_that(".cens() can come first or last in formula.list", {
  skip_if_not_installed("rootSolve")

  set.seed(7)
  d <- msmdata
  n <- nrow(d)

  # First: a selection/entry model
  d$S <- rbinom(n, 1L, plogis(-2.2 + 0.15 * d$X1_0))
  W_first <- weightitMSM(list(.cens(S) ~ X1_0 + X2_0,
                              A_1 ~ X1_0 + X2_0),
                         data = d, method = "glm")

  expect_true(all(W_first$weights[d$S == 1L] == 0))
  expect_M_parts_okay(W_first, tolerance = eps)

  # Last: loss to follow-up before the outcome is measured
  d$C_3 <- rbinom(n, 1L, plogis(-2 + 0.2 * d$X1_2))
  W_last <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                             A_2 ~ X1_1 + X2_1 + A_1,
                             A_3 ~ X1_2 + X2_2 + A_2,
                             .cens(C_3) ~ X1_2 + X2_2 + A_3),
                        data = d, method = "glm")

  expect_true(all(W_last$weights[d$C_3 == 1L] == 0))
  expect_M_parts_okay(W_last, tolerance = eps)
})

test_that("weightitMSM censoring composes with stabilize and by", {
  skip_if_not_installed("rootSolve")

  d <- make_msm_cens_data()
  cens <- d$C_2 == 1L

  W_stab <- weightitMSM(msm_cens_formulas, data = d, method = "glm",
                        stabilize = TRUE)

  expect_true(all(is.finite(W_stab$weights)))
  expect_length(W_stab$weights, nrow(d))
  expect_true(all(W_stab$weights[cens] == 0))
  expect_M_parts_okay(W_stab, tolerance = eps)

  # A numerator formula list has one entry per entry of `formula.list`, censoring
  # entries included
  expect_no_error({
    W_num <- weightitMSM(msm_cens_formulas, data = d, method = "glm",
                         num.formula = list(~1, ~A_1, ~A_1, ~ A_1 + A_2))
  })
  expect_true(all(W_num$weights[cens] == 0))
  expect_true(all(is.finite(W_num$weights)))

  # The censoring entry's numerator formula is honored
  expect_true(any(vapply(W_num$stabilization,
                         function(f) identical(deparse1(f), "~A_1"), logical(1L))))

  # Skipping the censoring entries is wrong
  expect_error(weightitMSM(msm_cens_formulas, data = d, method = "glm",
                           num.formula = list(~1, ~A_1, ~ A_1 + A_2)),
               "as many entries")

  d$G <- factor(d$X2_0)
  W_by <- weightitMSM(msm_cens_formulas, data = d, method = "glm", by = ~G)

  expect_length(attr(W_by, "Mparts.list"), 4L * nlevels(d$G))
  expect_M_parts_okay(W_by, tolerance = eps)
})

test_that("weightitMSM censoring handles missing values by risk set", {
  d <- make_msm_cens_data()

  # Covariates missing only for already-censored units must not trigger the
  # missingness machinery
  W <- weightitMSM(msm_cens_formulas, data = d, method = "glm")
  expect_null(W$missing)

  # A missing value in a unit still under observation is an error
  dbad <- d
  is.na(dbad$A_3[which(d$C_2 == 0L)[1L]]) <- TRUE
  expect_error(weightitMSM(msm_cens_formulas, data = dbad, method = "glm"),
               "still under observation")
})

test_that("weightitMSM accepts an empty censoring formula", {
  skip_if_not_installed("rootSolve")

  d <- make_msm_cens_data()
  cens <- d$C_2 == 1L

  f <- list(A_1 ~ X1_0 + X2_0,
            A_2 ~ X1_1 + X2_1 + A_1,
            .cens(C_2) ~ 1,
            A_3 ~ X1_2 + X2_2 + A_2)

  W <- weightitMSM(f, data = d, method = "glm", include.obj = TRUE)

  # The censoring factor is the marginal 1/P(C = 0) among those still at risk
  w1 <- weightit(A_1 ~ X1_0 + X2_0, data = d, method = "glm")$weights
  w2 <- weightit(A_2 ~ X1_1 + X2_1 + A_1, data = d, method = "glm")$weights
  cw <- ifelse(cens, 0, 1 / mean(!cens))
  w3 <- rep.int(1, nrow(d))
  w3[!cens] <- weightit(A_3 ~ X1_2 + X2_2 + A_2, data = d[!cens, ],
                        method = "glm")$weights

  expect_equal(unname(W$weights), unname(w1 * w2 * cw * w3), tolerance = eps)

  # The rest of the censoring infrastructure is untouched: the risk sets are the
  # same as with a covariate-dependent censoring model, the models after censoring
  # are fit on the at-risk units only, and NAs there are still tolerated
  W_cov <- weightitMSM(msm_cens_formulas, data = d, method = "glm")
  expect_identical(W$at.risk, W_cov$at.risk)
  expect_identical(W$cens.time, W_cov$cens.time)
  expect_identical(stats::nobs(W$obj$A_3), sum(!cens))
  expect_null(W$missing)

  expect_length(attr(W, "Mparts.list"), 4L)
  expect_M_parts_okay(W, tolerance = eps)

  # Stabilization divides by the same marginal model, so the censoring factor
  # becomes exactly 1 for the units still under observation: the stabilized
  # weights differ from the unstabilized ones only by the treatment numerators,
  # which is the factor mean(!cens) that `cw` contributes above
  W_stab <- weightitMSM(f, data = d, method = "glm", stabilize = TRUE)
  W_stab_cov <- weightitMSM(msm_cens_formulas, data = d, method = "glm",
                            stabilize = TRUE)

  expect_true(all(W_stab$weights[cens] == 0))
  expect_equal(unname(W_stab$weights / W$weights)[!cens],
               unname(W_stab_cov$weights / W_cov$weights)[!cens],
               tolerance = eps)
  expect_M_parts_okay(W_stab, tolerance = eps)

  # Empty and non-empty censoring formulas can be mixed
  d$C_3 <- 0L
  at.risk <- which(!cens)
  set.seed(8L)
  d$C_3[at.risk] <- rbinom(length(at.risk), 1L,
                           prob = plogis(-1.5 + 0.3 * d$X1_2[at.risk]))
  is.na(d$C_3[cens]) <- TRUE

  W_mix <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                            A_2 ~ X1_1 + X2_1 + A_1,
                            .cens(C_2) ~ 1,
                            A_3 ~ X1_2 + X2_2 + A_2,
                            .cens(C_3) ~ X1_2 + A_3),
                       data = d, method = "glm")

  expect_equal(unname(W_mix$cens.time), c(3L, 5L))
  expect_true(all(W_mix$weights[cens | d$C_3 == 1L] == 0))
  expect_true(all(W_mix$weights[!cens & d$C_3 == 0L] > 0))
  expect_M_parts_okay(W_mix, tolerance = eps)
})

test_that("weightitMSM accepts an empty treatment formula", {
  # An empty right side is allowed in `weightitMSM()` just as it is in
  # `weightit()`; the corresponding time point gets a marginal model
  d <- msmdata

  W <- weightitMSM(list(A_1 ~ 1, A_2 ~ X1_1 + A_1), data = d, method = "glm")

  w1 <- weightit(A_1 ~ 1, data = d, method = "glm")$weights
  w2 <- weightit(A_2 ~ X1_1 + A_1, data = d, method = "glm")$weights

  expect_equal(unname(W$weights), unname(w1 * w2), tolerance = eps)

  # Printing reports the absent covariates at every time point, not just baseline
  expect_output(print(W), "baseline: \\(none\\)")
  expect_output(print(weightitMSM(list(A_1 ~ X1_0, A_2 ~ 1), data = d,
                                  method = "glm")),
                "after time 1: \\(none\\)")
})

test_that("weightitMSM censoring error handling", {
  d <- make_msm_cens_data()

  # At least one treatment model is required, empty formula or not
  expect_error(weightitMSM(list(.cens(C_2) ~ X1_1), data = d, method = "glm"),
               "at least one treatment model")
  expect_error(weightitMSM(list(.cens(C_2) ~ 1), data = d, method = "glm"),
               "at least one treatment model")

  # `method = "cbps"` with `is.MSM.method = TRUE` supports censoring (see the
  # simultaneous CBPS tests below); it is no longer rejected
  expect_no_error(weightitMSM(msm_cens_formulas, data = d, method = "cbps",
                              is.MSM.method = TRUE))

  # Without censoring, missing treatment values remain an error
  dna <- msmdata
  is.na(dna$A_3[1L]) <- TRUE
  expect_error(weightitMSM(list(A_1 ~ X1_0, A_2 ~ X1_1 + A_1, A_3 ~ X1_2 + A_2),
                           data = dna, method = "glm"),
               "[Mm]issing values")
})

# ---- Censoring in the simultaneous CBPS MSM (is.MSM.method = TRUE) ----------
#
# `weightitMSM2cbps()` estimates one weight per unit satisfying the balance
# conditions at every time point at once. Each condition is evaluated on that time
# point's risk set: the covariates of already-censored units are missing, and
# `0 * NA` is NA, so down-weighting cannot rescue them -- their rows are zeroed.

test_that("simultaneous CBPS is unchanged when there is no censoring", {
  skip_if_not_installed("rootSolve")

  # A regression guard on the risk-set restructuring: with no censoring every unit
  # is at risk at every time point, so the computation must be identical.
  W <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                        A_2 ~ X1_1 + X2_1 + A_1,
                        A_3 ~ X1_2 + X2_2 + A_2),
                   data = msmdata, method = "cbps", is.MSM.method = TRUE,
                   include.obj = TRUE)

  # The moment conditions are solved essentially exactly
  expect_lt(W$obj$value, 1e-8)

  # Treated and control covariate means agree at every time point
  for (tp in c("A_1", "A_2", "A_3")) {
    cvs <- switch(tp, A_1 = c("X1_0", "X2_0"), A_2 = c("X1_1", "X2_1"),
                  A_3 = c("X1_2", "X2_2"))
    A <- msmdata[[tp]]
    X <- as.matrix(msmdata[cvs])

    m1 <- colSums(X[A == 1, , drop = FALSE] * W$weights[A == 1]) / sum(W$weights[A == 1])
    m0 <- colSums(X[A == 0, , drop = FALSE] * W$weights[A == 0]) / sum(W$weights[A == 0])

    expect_lt(max(abs(m1 - m0)), 1e-6)
  }
})

test_that("simultaneous CBPS solves the conditions with censoring present", {
  skip_if_not_installed("rootSolve")

  d <- make_msm_cens_data()
  cens <- d$C_2 == 1L

  expect_no_condition({
    W <- weightitMSM(msm_cens_formulas, data = d, method = "cbps",
                     is.MSM.method = TRUE, include.obj = TRUE)
  })

  expect_true(all(W$weights[cens] == 0))
  expect_true(all(W$weights[!cens] > 0))
  expect_true(all(is.finite(W$weights)))

  # All the moment conditions, treatment and censoring, are satisfied jointly
  expect_lt(W$obj$value, 1e-8)

  # Each treatment condition holds within its own risk set
  for (tp in c("A_1", "A_2", "A_3")) {
    cvs <- switch(tp, A_1 = c("X1_0", "X2_0"), A_2 = c("X1_1", "X2_1"),
                  A_3 = c("X1_2", "X2_2"))
    A <- d[[tp]]
    X <- as.matrix(d[cvs])

    ok <- W$at.risk[, tp] & !is.na(A)
    t1 <- ok & A == 1
    t0 <- ok & A == 0

    m1 <- colSums(X[t1, , drop = FALSE] * W$weights[t1]) / sum(W$weights[t1])
    m0 <- colSums(X[t0, , drop = FALSE] * W$weights[t0]) / sum(W$weights[t0])

    expect_lt(max(abs(m1 - m0)), 1e-6)
  }

  # Matching the treatment path, this method supplies no M-estimation parts
  expect_null(attr(W, "Mparts", exact = TRUE))
  expect_null(attr(W, "Mparts.list", exact = TRUE))

  # The censoring models are described separately, as on the product path
  expect_named(W$cens.list, "C_2")
  expect_length(W$treat.list, 3L)
})

test_that("simultaneous and product CBPS solve different problems", {
  skip_if_not_installed("rootSolve")

  d <- make_msm_cens_data()
  cens <- d$C_2 == 1L

  W_sim <- weightitMSM(msm_cens_formulas, data = d, method = "cbps",
                       is.MSM.method = TRUE)
  W_prod <- suppressMessages(
    weightitMSM(msm_cens_formulas, data = d, method = "cbps",
                is.MSM.method = FALSE)
  )

  # Both zero exactly the censored units, but the weights themselves differ:
  # one set satisfies all the conditions jointly, the other multiplies per-time-point
  # solutions.
  expect_identical(W_sim$weights == 0, W_prod$weights == 0)
  expect_true(all(W_sim$weights[cens] == 0))
  expect_true(all(W_prod$weights[cens] == 0))
  expect_false(isTRUE(all.equal(unname(W_sim$weights), unname(W_prod$weights))))
})

test_that("simultaneous CBPS censoring handles over-identification and degeneracy", {
  skip_if_not_installed("rootSolve")

  d <- make_msm_cens_data()
  cens <- d$C_2 == 1L

  # over = TRUE stacks the GLM scores with the balance conditions
  W_over <- weightitMSM(msm_cens_formulas, data = d, method = "cbps",
                        is.MSM.method = TRUE, over = TRUE)

  expect_true(all(W_over$weights[cens] == 0))
  expect_true(all(is.finite(W_over$weights)))

  # A censoring time point at which nobody is censored contributes weights of 1
  d0 <- d
  d0$C_2 <- 0L
  d0[c("X1_2", "X2_2", "A_3")] <- msmdata[c("X1_2", "X2_2", "A_3")]

  W0 <- weightitMSM(msm_cens_formulas, data = d0, method = "cbps",
                    is.MSM.method = TRUE)

  expect_true(all(W0$weights > 0))
  expect_true(all(is.finite(W0$weights)))
})

test_that("only methods that can receive the risk sets may fit all times at once", {
  d <- make_msm_cens_data()

  # `cbps` is currently the only method with `msm_method_available = TRUE`, and
  # user-defined methods are barred from `is.MSM.method = TRUE` independently of
  # censoring, so the censoring-specific guard in `weightitMSM()` is defensive
  # against a future method rather than reachable today. What is reachable:
  user_msm <- function(covs, treat, s.weights, subset, ...) {
    list(w = rep.int(1, length(treat[subset])))
  }

  expect_error(weightitMSM(msm_cens_formulas, data = d, method = user_msm,
                           is.MSM.method = TRUE),
               "is.MSM.method")

  # And a method with no simultaneous version simply falls back with a warning
  expect_warning(weightitMSM(msm_cens_formulas, data = d, method = "glm",
                             is.MSM.method = TRUE),
                 "cannot be used with a single model")
})

# ---- Regression tests for two latent CBPS MSM closure bugs ------------------

test_that("simultaneous CBPS handles multi-category time points correctly", {
  skip_if_not_installed("rootSolve")

  # The multi-category balance closure used to reference `treat.list[[i]]` with `i`
  # resolved lazily to the LAST time point. When that time point is binary its
  # levels are NULL and `combn(NULL, 2)` errors; when the level counts merely
  # differ, the wrong level set was used silently.
  set.seed(5)
  d <- msmdata
  d$M1 <- factor(sample(c("a", "b", "c"), nrow(d), replace = TRUE))
  d$M2 <- factor(sample(c("p", "q", "r", "s"), nrow(d), replace = TRUE))

  # A multi-category time point followed by a binary one
  expect_no_error({
    W1 <- suppressMessages(
      weightitMSM(list(M1 ~ X1_0 + X2_0, A_3 ~ X1_1 + M1),
                  data = d, method = "cbps", is.MSM.method = TRUE)
    )
  })
  expect_true(all(is.finite(W1$weights)))

  # Two multi-category time points with different numbers of levels
  expect_no_error({
    W2 <- weightitMSM(list(M1 ~ X1_0, M2 ~ X1_1 + M1),
                      data = d, method = "cbps", is.MSM.method = TRUE)
  })
  expect_true(all(is.finite(W2$weights)))
})

test_that("bal.tab's limitations with censoring are as documented", {
  skip_if_not_installed("cobalt")

  d <- make_msm_cens_data()

  # cobalt rejects the NA treatments that censoring necessarily produces; the
  # documented workaround is to subset using the `at.risk` component.
  W <- weightitMSM(msm_cens_formulas, data = d, method = "glm")
  expect_error(cobalt::bal.tab(W), "[Mm]issing")

  ar <- W$at.risk[, "A_3"]
  expect_no_error(cobalt::bal.tab(A_3 ~ X1_2 + X2_2 + A_2, data = d[ar, ],
                                  weights = W$weights[ar], s.d.denom = "pooled"))

  # For a point-treatment censoring model, every censored unit has a weight of 0,
  # which cobalt also rejects; balance is assessed on a stacked pseudo-sample.
  dp <- make_cens_data()
  Wp <- weightit(.cens(C) ~ X1 + X2 + X3, data = dp, method = "ebal")

  expect_error(cobalt::bal.tab(Wp), "zero")

  covs <- dp[c("X1", "X2", "X3")]
  u <- which(Wp$treat == 0)

  b <- cobalt::bal.tab(rbind(covs[u, ], covs),
                       treat = c(rep.int(0L, length(u)), rep.int(1L, nrow(covs))),
                       weights = c(Wp$weights[u], rep.int(1, nrow(covs))),
                       estimand = "ATT", s.d.denom = "treated")

  # ebal drives these to 0 by construction
  expect_true(all(abs(b$Balance$Diff.Adj) < 1e-6))
})

# The zero-weight path above and the handling of degenerate outcome models (no
# events, no covariates, collinear covariates) were developed independently, so
# their combinations are checked here explicitly.
test_that("degenerate Cox models still work with zero weights", {
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  d <- make_cens_data()

  cutoff <- quantile(d$Y_S, .8)
  d$event <- as.numeric(d$Y_S < cutoff)
  d$time <- pmin(d$Y_S, cutoff)
  is.na(d$time[d$C == 1L]) <- TRUE
  is.na(d$event[d$C == 1L]) <- TRUE

  #A collinear duplicate of an existing covariate
  d$X2b <- d$X2

  W <- weightit(.cens(C) ~ X1 + X3, data = d, method = "glm")

  expect_true(any(W$weights == 0))

  keep <- W$weights > 0
  d_keep <- d[keep, ]
  w_keep <- W$weights[keep]

  #No events among the units that contribute
  d_ne <- d
  d_ne$event[keep] <- 0

  expect_warning({
    fit_ne <- coxph_weightit(survival::Surv(time, event) ~ A + X2,
                             data = d_ne, weightit = W)
  }, "no events", ignore.case = TRUE)

  expect_true(all(is.na(coef(fit_ne))))
  expect_true(all(is.na(vcov(fit_ne))))
  expect_identical(fit_ne$nevent, 0)

  #No covariates
  expect_no_condition({
    fit_null <- coxph_weightit(survival::Surv(time, event) ~ 1,
                               data = d, weightit = W)
  })

  expect_length(coef(fit_null), 0L)
  expect_identical(dim(vcov(fit_null)), c(0L, 0L))

  #Collinear covariates: the aliased coefficient is NA and the estimable ones
  #match a coxph() fit to the contributing units only
  expect_no_condition({
    fit_rd <- coxph_weightit(survival::Surv(time, event) ~ A + X2 + X2b,
                             data = d, weightit = W, vcov = "HC0")
  })

  expect_identical(is.na(coef(fit_rd)), c(A = FALSE, X2 = FALSE, X2b = TRUE))
  expect_false(anyNA(vcov(fit_rd, complete = FALSE)))

  fit_rd_ref <- survival::coxph(survival::Surv(time, event) ~ A + X2 + X2b,
                                data = d_keep, weights = w_keep, robust = TRUE)

  expect_equal(coef(fit_rd), coef(fit_rd_ref), tolerance = eps)
  expect_equal(unname(vcov(fit_rd, complete = FALSE)),
               unname(vcov(fit_rd_ref)[!is.na(coef(fit_rd_ref)),
                                       !is.na(coef(fit_rd_ref))]),
               tolerance = eps)

  #vcov = "const" is the model-based variance of that same fit
  expect_warning({
    fit_const <- coxph_weightit(survival::Surv(time, event) ~ A + X2,
                                data = d, weightit = W, vcov = "const")
  }, "const.*should not be used", ignore.case = TRUE)

  fit_const_ref <- survival::coxph(survival::Surv(time, event) ~ A + X2,
                                   data = d_keep, weights = w_keep,
                                   robust = FALSE)

  expect_equal(unname(vcov(fit_const)), unname(vcov(fit_const_ref)),
               tolerance = eps)
})

# ---- Agreement with a hand-written IPCW/IPTW implementation -----------------
#
# The tests above check the censoring weights against identities and against other
# methods, which verifies internal consistency but shares WeightIt's machinery. These
# reproduce the weights from scratch instead: one `glm()` per model, fit by hand among the
# units still under observation, and the reciprocal of the predicted probability of the
# value actually observed. Nothing below calls `weightit()` to build the reference, so a
# systematic error anywhere in the censoring path -- the risk sets, which units each model
# is fit on, the direction of the censoring indicator, or how the per-time-point factors
# are combined -- shows up as a numeric disagreement.
#
# `quasibinomial()` rather than `binomial()` only to avoid the non-integer-successes
# warning when sampling weights are passed; the two give identical fitted values.

# Predicted probability of the observed value of the left-hand-side variable, from a
# logistic regression fit among `subset` only. `NA` outside `subset`.
p_observed <- function(formula, data, subset = NULL, s.weights = NULL) {
  n <- nrow(data)

  if (is_null(subset)) subset <- rep.int(TRUE, n)

  data$.sw <- s.weights %or% rep.int(1, n)

  fit <- glm(formula, data = data[subset, , drop = FALSE],
             family = quasibinomial, weights = .sw)

  p1 <- rep(NA_real_, n)
  p1[subset] <- unname(fitted(fit))

  y <- data[[all.vars(formula)[1L]]]

  ifelse(y == 1, p1, 1 - p1)
}

# Inverse probability of treatment weight for one time point: 1/P(A = a | X) among the
# units at risk, and 1 for everyone else (they contribute nothing at this time point).
iptw_factor <- function(formula, data, at.risk, s.weights = NULL) {
  f <- rep.int(1, nrow(data))
  f[at.risk] <- (1 / p_observed(formula, data, at.risk, s.weights))[at.risk]
  f
}

# Inverse probability of censoring weight for one time point: 1/P(C = 0 | X) for the
# units still under observation, exactly 0 for those censored here, and 1 for units
# outside the risk set.
ipcw_factor <- function(formula, data, at.risk, s.weights = NULL) {
  C <- data[[all.vars(formula)[1L]]]

  f <- rep.int(1, nrow(data))
  f[at.risk] <- ifelse(C[at.risk] == 1L, 0,
                       (1 / p_observed(formula, data, at.risk, s.weights))[at.risk])
  f
}

test_that("point-treatment censoring weights match a hand-written IPCW", {
  d <- make_cens_data()

  all.at.risk <- rep.int(TRUE, nrow(d))

  # Plain
  W <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "glm")

  manual <- ipcw_factor(C ~ X1 + X2 + X3 + X4, d, all.at.risk)

  expect_equal(unname(W$weights), unname(manual), tolerance = eps)
  expect_true(all(manual[d$C == 1L] == 0))
  expect_true(all(manual[d$C == 0L] > 0))

  # Sampling weights enter the censoring model as regression weights
  W_sw <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "glm",
                   s.weights = "SW")

  manual_sw <- ipcw_factor(C ~ X1 + X2 + X3 + X4, d, all.at.risk,
                           s.weights = d$SW)

  expect_equal(unname(W_sw$weights), unname(manual_sw), tolerance = eps)
  expect_not_equal(unname(manual_sw), unname(manual))

  # `by` fits a separate censoring model within each stratum
  W_by <- weightit(.cens(C) ~ X1 + X2, data = d, method = "glm", by = ~G)

  manual_by <- rep.int(NA_real_, nrow(d))
  for (g in levels(d$G)) {
    in.g <- d$G == g
    manual_by[in.g] <- ipcw_factor(C ~ X1 + X2, d, in.g)[in.g]
  }

  expect_equal(unname(W_by$weights), unname(manual_by), tolerance = eps)

  # Stabilization divides by a marginal censoring model, so the numerator is
  # P(C = 0) rather than P(C = 0 | X)
  W_st <- weightit(.cens(C) ~ X1 + X2 + X3 + X4, data = d, method = "glm",
                   stabilize = TRUE)

  manual_st <- manual * p_observed(C ~ 1, d, all.at.risk)

  expect_equal(unname(W_st$weights), unname(manual_st), tolerance = eps)
})

test_that("longitudinal censoring weights match a hand-written IPTW x IPCW product", {
  d <- make_msm_cens_data()

  # `msm_cens_formulas` is A_1, A_2, .cens(C_2), A_3: censoring occurs after the second
  # treatment, so the first three models are fit on everyone and only the A_3 model is
  # restricted to the units still under observation.
  W <- weightitMSM(msm_cens_formulas, data = d, method = "glm")

  at.risk <- rep.int(TRUE, nrow(d))

  f_a1 <- iptw_factor(A_1 ~ X1_0 + X2_0, d, at.risk)
  f_a2 <- iptw_factor(A_2 ~ X1_1 + X2_1 + A_1, d, at.risk)
  f_c2 <- ipcw_factor(C_2 ~ X1_1 + X2_1 + A_1, d, at.risk)

  at.risk <- at.risk & d$C_2 == 0L

  f_a3 <- iptw_factor(A_3 ~ X1_2 + X2_2 + A_2, d, at.risk)

  manual <- f_a1 * f_a2 * f_c2 * f_a3

  expect_equal(unname(W$weights), unname(manual), tolerance = eps)

  # The risk sets WeightIt used are the ones the manual version assumed
  expect_identical(unname(W$at.risk[, "A_3"]), unname(at.risk))
  expect_identical(unname(W$at.risk[, "C_2"]), rep.int(TRUE, nrow(d)))

  # Censored units are zeroed out by the censoring factor alone
  expect_true(all(manual[d$C_2 == 1L] == 0))
  expect_true(all(manual[d$C_2 == 0L] > 0))
  expect_identical(sum(W$weights == 0), sum(d$C_2 == 1L))

  # Sampling weights again enter every model as regression weights (`msmdata` has none
  # of its own, so one is built here)
  set.seed(99L)
  d$SW <- 1 / plogis(0.5 + 0.4 * d$X1_0 - 0.3 * d$A_1)

  W_sw <- weightitMSM(msm_cens_formulas, data = d, method = "glm",
                      s.weights = "SW")

  at.risk <- rep.int(TRUE, nrow(d))
  g_a1 <- iptw_factor(A_1 ~ X1_0 + X2_0, d, at.risk, d$SW)
  g_a2 <- iptw_factor(A_2 ~ X1_1 + X2_1 + A_1, d, at.risk, d$SW)
  g_c2 <- ipcw_factor(C_2 ~ X1_1 + X2_1 + A_1, d, at.risk, d$SW)
  at.risk <- at.risk & d$C_2 == 0L
  g_a3 <- iptw_factor(A_3 ~ X1_2 + X2_2 + A_2, d, at.risk, d$SW)

  expect_equal(unname(W_sw$weights), unname(g_a1 * g_a2 * g_c2 * g_a3),
               tolerance = eps)
  expect_not_equal(unname(W_sw$weights), unname(W$weights))
})

test_that("two censoring time points match a hand-written product with nested risk sets", {
  # The distinctive part of the longitudinal implementation is that each model is fit on
  # the units surviving every earlier censoring event, so the risk sets are nested. This
  # reproduces that by hand.
  set.seed(1234L)

  d <- msmdata
  n <- nrow(d)

  d$C_2 <- rbinom(n, 1L, prob = plogis(-2 + 0.2 * d$X1_1 - 0.3 * d$A_1))

  u <- which(d$C_2 == 0L)
  d$C_3 <- NA_integer_
  d$C_3[u] <- rbinom(length(u), 1L, prob = plogis(-1.8 + 0.25 * d$X1_2[u]))

  # Everything measured after each censoring event is unobserved for the units it removed
  is.na(d[d$C_2 == 1L, c("X1_2", "X2_2", "A_2", "A_3")]) <- TRUE
  is.na(d[which(d$C_3 == 1L), "A_3"]) <- TRUE

  formulas <- list(A_1 ~ X1_0 + X2_0,
                   .cens(C_2) ~ X1_1 + X2_1 + A_1,
                   A_2 ~ X1_1 + X2_1 + A_1,
                   .cens(C_3) ~ X1_2 + A_2,
                   A_3 ~ X1_2 + X2_2 + A_2)

  W <- weightitMSM(formulas, data = d, method = "glm")

  ar1 <- rep.int(TRUE, n)                       # A_1 and C_2: everyone
  f_a1 <- iptw_factor(A_1 ~ X1_0 + X2_0, d, ar1)
  f_c2 <- ipcw_factor(C_2 ~ X1_1 + X2_1 + A_1, d, ar1)

  ar2 <- ar1 & d$C_2 == 0L                      # A_2 and C_3: survivors of C_2
  f_a2 <- iptw_factor(A_2 ~ X1_1 + X2_1 + A_1, d, ar2)
  f_c3 <- ipcw_factor(C_3 ~ X1_2 + A_2, d, ar2)

  ar3 <- ar2 & d$C_3 == 0L                      # A_3: survivors of both
  f_a3 <- iptw_factor(A_3 ~ X1_2 + X2_2 + A_2, d, ar3)

  manual <- f_a1 * f_c2 * f_a2 * f_c3 * f_a3

  expect_equal(unname(W$weights), unname(manual), tolerance = eps)

  # The risk sets shrink as assumed, and censoring at either point zeroes the weight
  expect_identical(unname(W$at.risk[, "C_2"]), unname(ar1))
  expect_identical(unname(W$at.risk[, "A_2"]), unname(ar2))
  expect_identical(unname(W$at.risk[, "C_3"]), unname(ar2))
  expect_identical(unname(W$at.risk[, "A_3"]), unname(ar3))

  expect_true(all(manual[!ar3] == 0))
  expect_true(all(manual[ar3] > 0))
  expect_equal(unname(W$cens.time), c(2L, 4L))
})

test_that("stabilized longitudinal censoring weights match a hand-written product", {
  d <- make_msm_cens_data()
  n <- nrow(d)

  W <- weightitMSM(msm_cens_formulas, data = d, method = "glm", stabilize = TRUE)

  # With `stabilize = TRUE` and no `num.formula`, every numerator -- treatment or
  # censoring alike -- is a model saturated in the previous treatments (an intercept
  # when no treatment precedes it). Censoring indicators are not counted among the
  # previous treatments.
  at.risk <- rep.int(TRUE, n)

  a1 <- iptw_factor(A_1 ~ X1_0 + X2_0, d, at.risk) *
    p_observed(A_1 ~ 1, d, at.risk)
  a2 <- iptw_factor(A_2 ~ X1_1 + X2_1 + A_1, d, at.risk) *
    p_observed(A_2 ~ A_1, d, at.risk)
  c2 <- ipcw_factor(C_2 ~ X1_1 + X2_1 + A_1, d, at.risk) *
    p_observed(C_2 ~ A_1 * A_2, d, at.risk)

  at.risk <- at.risk & d$C_2 == 0L

  num_a3 <- rep.int(1, n)
  num_a3[at.risk] <- p_observed(A_3 ~ A_1 * A_2, d, at.risk)[at.risk]
  a3 <- iptw_factor(A_3 ~ X1_2 + X2_2 + A_2, d, at.risk) * num_a3

  manual <- a1 * a2 * c2 * a3

  expect_equal(unname(W$weights), unname(manual), tolerance = eps)

  # The censoring numerator is a marginal model of C, not a censoring model: a censoring
  # numerator would carry the same (1 - C) factor as the denominator and give 0/0 for the
  # censored units. Here they stay finite and exactly 0.
  expect_true(all(is.finite(manual)))
  expect_true(all(manual[d$C_2 == 1L] == 0))

  # Stabilizing changes the weights but not which units are zeroed
  W_un <- weightitMSM(msm_cens_formulas, data = d, method = "glm")
  expect_not_equal(unname(W$weights), unname(W_un$weights))
  expect_identical(W$weights == 0, W_un$weights == 0)
})
