# M-estimation with `by`: estimating weights separately within each `by` stratum
# should yield an asymptotic variance (via glm_weightit(..., vcov = "asympt"))
# essentially identical to estimating the weights from a single model in which
# the `by` variable is fully interacted with all the covariates. The per-stratum
# Mparts are expanded to full-sample size (.expand_Mparts_by()) and stacked
# block-diagonally by the shared M-estimation machinery, which is exactly the
# structure of the fully-interacted model.

skip_on_cran()

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

test_data <- readRDS(test_path("fixtures", "test_data.rds"))
# `by` variable must be a factor so the interacted model treats it categorically
# (rather than as a single numeric column) -- only then are the two specifications
# comparable.
test_data$G <- factor(test_data$X5, labels = c("g0", "g1"))

# Compare vcov(outcome model) using by-stratified weights vs. fully-interacted
# weights, for a given method/estimand.
expect_by_equiv <- function(method, estimand = "ATE", ..., outcome = Y_B ~ A * G, treat = "A",
                            ps_covs = "X1 + X2 + X3 + X4", tol = 1e-5) {
  ps_by  <- as.formula(sprintf("%s ~ %s", treat, ps_covs))
  ps_int <- as.formula(sprintf("%s ~ G * (%s)", treat, ps_covs))

  W_by  <- weightit(ps_by,  data = test_data, by = ~G, method = method,
                    estimand = estimand, ...)
  W_int <- weightit(ps_int, data = test_data, method = method,
                    estimand = estimand, ...)

  # The by-fit must now carry the stacked M-estimation components
  expect_null(attr(W_by, "Mparts", exact = TRUE))
  expect_false(is_null(attr(W_by, "Mparts.list", exact = TRUE)))

  expect_M_parts_okay(W_by, tolerance = eps)
  expect_M_parts_okay(W_int, tolerance = eps)

  f_by  <- glm_weightit(outcome, data = test_data, weightit = W_by,  vcov = "asympt",
                        family = binomial)
  f_int <- glm_weightit(outcome, data = test_data, weightit = W_int, vcov = "asympt",
                        family = binomial)

  expect_equal(unname(coef(f_by)), unname(coef(f_int)), tolerance = eps)
  expect_equal(unname(vcov(f_by)), unname(vcov(f_int)), tolerance = tol)
}

test_that("by-stratified M-estimation matches interacted model: glm", {
  expect_by_equiv("glm")
  expect_by_equiv("glm", link = "probit")
  expect_by_equiv("glm", link = "cloglog")
})

test_that("by-stratified M-estimation matches interacted model: glm, link = 'br.logit'", {
  skip_if_not_installed("brglm2")
  expect_by_equiv("glm", link = "br.logit")
  expect_by_equiv("glm", link = "br.probit")
  expect_by_equiv("glm", link = "br.cloglog")
})

test_that("by-stratified M-estimation matches interacted model: ebal", {
  # ebal supplies Mparts only for exact balance (tols == 0, the default).
  expect_by_equiv("ebal")
})

test_that("by-stratified M-estimation matches interacted model: cbps", {
  # just-identified (the default, over = FALSE) supplies Mparts.
  expect_by_equiv("cbps")
  expect_by_equiv("cbps", link = "probit")
  expect_by_equiv("cbps", link = "cloglog")
})

test_that("by-stratified M-estimation matches interacted model: ipt", {
  expect_by_equiv("ipt")
  expect_by_equiv("ipt", link = "probit")
  expect_by_equiv("ipt", link = "cloglog")
})

test_that("by-stratified M-estimation composes with stabilize (glm)", {
  # Within-stratum stabilization (numerator ~1 per group) is equivalent to a
  # group-specific numerator (~G) in the interacted specification.
  W_by <- weightit(A ~ X1 + X2 + X3 + X4, data = test_data, by = ~G,
                   method = "glm", estimand = "ATE", stabilize = ~1)
  W_int <- weightit(A ~ G * (X1 + X2 + X3 + X4), data = test_data,
                    method = "glm", estimand = "ATE", stabilize = ~G)

  expect_false(is_null(attr(W_by, "Mparts.list", exact = TRUE)))

  f_by  <- glm_weightit(Y_C ~ A * G, data = test_data, weightit = W_by,  vcov = "asympt")
  f_int <- glm_weightit(Y_C ~ A * G, data = test_data, weightit = W_int, vcov = "asympt")

  expect_equal(unname(coef(f_by)), unname(coef(f_int)), tolerance = eps)
  expect_equal(unname(vcov(f_by)), unname(vcov(f_int)), tolerance = 1e-7)
})

test_that("by-stratified M-estimation matches interacted model: multi-category", {
  expect_by_equiv("glm", outcome = Y_B ~ Am * G, treat = "Am")
  expect_by_equiv("ebal", outcome = Y_B ~ Am * G, treat = "Am")
  expect_by_equiv("cbps", outcome = Y_B ~ Am * G, treat = "Am")
  expect_by_equiv("ipt", outcome = Y_B ~ Am * G, treat = "Am")
})

# Note: equivalence doesn't hold for continuous treatments, so not tested


# ---- Longitudinal treatments (weightitMSM) ----
# The MSM weight is the product across time points of each time point's weight;
# with `by`, the per-(time point x stratum) Mparts are combined so the resulting
# standard errors match a specification that interacts the (baseline) `by`
# variable with all covariates at every time point.

test_that("by-stratified MSM M-estimation matches interacted model", {
  data("msmdata", package = "WeightIt", envir = environment())
  md <- msmdata
  # Baseline grouping variable; excluded from the covariate formulas since it is
  # the stratifier. X2_0 has both treatment levels in every stratum at each time.
  md$G <- factor(md$X2_0, labels = c("g0", "g1"))

  fl_by  <- list(A_1 ~ X1_0,
                 A_2 ~ X1_1 + X2_1 + A_1 + X1_0)
  fl_int <- list(A_1 ~ G * (X1_0),
                 A_2 ~ G * (X1_1 + X2_1 + A_1 + X1_0))

  expect_equiv_msm <- function(method, ...) {
    W_by  <- suppressMessages(weightitMSM(fl_by,  data = md, by = ~G,
                                          method = method, ...))
    W_int <- suppressMessages(weightitMSM(fl_int, data = md, method = method, ...))

    # by-fit carries the combined stacked Mparts: one part per (time point x stratum)
    expect_false(is_null(attr(W_by, "Mparts.list", exact = TRUE)))
    expect_identical(length(attr(W_by, "Mparts.list", exact = TRUE)),
                     length(fl_by) * nlevels(md$G))
    expect_M_parts_okay(W_by, tolerance = eps)

    f_by  <- glm_weightit(Y_B ~ A_1 * A_2 * G, data = md, weightit = W_by,
                          vcov = "asympt", family = binomial)
    f_int <- glm_weightit(Y_B ~ A_1 * A_2 * G, data = md, weightit = W_int,
                          vcov = "asympt", family = binomial)

    expect_equal(unname(coef(f_by)), unname(coef(f_int)), tolerance = eps)
    expect_equal(unname(vcov(f_by)), unname(vcov(f_int)), tolerance = 1e-4)
  }

  expect_equiv_msm("glm")
  # cbps M-estimation requires a separate model per time point (is.MSM.method =
  # FALSE); the single-model MSM CBPS does not support M-estimation.
  expect_equiv_msm("cbps", is.MSM.method = FALSE)
})

test_that("by-stratified MSM M-estimation composes with stabilization (glm)", {
  data("msmdata", package = "WeightIt", envir = environment())
  md <- msmdata
  md$G <- factor(md$X2_0, labels = c("g0", "g1"))

  fl_by  <- list(A_1 ~ X1_0,
                 A_2 ~ X1_1 + X2_1 + A_1 + X1_0)
  fl_int <- list(A_1 ~ G * (X1_0),
                 A_2 ~ G * (X1_1 + X2_1 + A_1 + X1_0))

  # by: default numerator (saturated prior treatments, fit within each stratum).
  W_by <- weightitMSM(fl_by, data = md, by = ~G, method = "glm", stabilize = TRUE)

  # Interacted analogue: the numerator must cross G with the prior treatments so
  # that it, too, is block-diagonal across strata (t = 1: ~G; t = 2: ~G * A_1).
  W_int <- weightitMSM(fl_int, data = md, method = "glm", stabilize = TRUE,
                       num.formula = list(~G, ~G * A_1))

  # Combined list has denominator + numerator parts for every (time point x stratum)
  expect_identical(length(attr(W_by, "Mparts.list", exact = TRUE)),
                   2L * length(fl_by) * nlevels(md$G))
  expect_M_parts_okay(W_by, tolerance = eps)

  f_by  <- glm_weightit(Y_B ~ A_1 * A_2 * G, data = md, weightit = W_by,
                        vcov = "asympt", family = binomial)
  f_int <- glm_weightit(Y_B ~ A_1 * A_2 * G, data = md, weightit = W_int,
                        vcov = "asympt", family = binomial)

  expect_equal(unname(coef(f_by)), unname(coef(f_int)), tolerance = eps)
  expect_equal(unname(vcov(f_by)), unname(vcov(f_int)), tolerance = 1e-4)
})