# Supplying `ps` bypasses model fitting and computes weights directly from the given
# propensity scores (`weightit2ps()` and its `.multi`/`.cont`/`.cens` variants). Before
# these tests only the binary `A ~ 1` case and the censoring case were covered, so the
# multi-category and continuous parsing were untested despite being documented.

skip_on_cran()

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

test_data <- readRDS(test_path("fixtures", "test_data.rds"))

test_that("a supplied ps gives the inverse-probability weights for a binary treatment", {
  d <- test_data

  ps <- fitted(glm(A ~ X1 + X2 + X3, data = d, family = binomial))

  W <- weightit(A ~ X1 + X2 + X3, data = d, ps = ps)

  # ATE: 1/e for the treated, 1/(1 - e) for the control
  expect_equal(unname(W$weights), unname(ifelse(d$A == 1, 1 / ps, 1 / (1 - ps))),
               tolerance = eps)
  expect_equal(unname(W$ps), unname(ps), tolerance = eps)

  # Identical to what method = "glm" produces from the same model
  W_glm <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm")
  expect_equal(unname(W$weights), unname(W_glm$weights), tolerance = eps)

  # ATT leaves the treated at 1
  W_att <- weightit(A ~ X1 + X2 + X3, data = d, ps = ps, estimand = "ATT")
  expect_ATT_weights_okay(W_att, tolerance = eps)

  # `ps` may also name a variable in `data`
  d$my_ps <- ps
  expect_equal(unname(weightit(A ~ X1 + X2 + X3, data = d, ps = "my_ps")$weights),
               unname(W$weights), tolerance = eps)

  # Supplying `ps` means no model is fit, so there are no M-estimation parts: the
  # scores are taken as fixed and known
  expect_null(attr(W, "Mparts", exact = TRUE))
})

test_that("a supplied ps works for a multi-category treatment", {
  d <- test_data

  # One binomial fit per level, normalized, then each unit's own-level probability --
  # which is the vector shape `weightit2ps.multi()` documents
  P <- vapply(levels(d$Am), function(lv) {
    fitted(glm(as.numeric(d$Am == lv) ~ X1 + X2 + X3, data = d,
               family = quasibinomial))
  }, numeric(nrow(d)))
  P <- P / rowSums(P)
  own <- P[cbind(seq_len(nrow(P)), match(as.character(d$Am), colnames(P)))]

  W <- weightit(Am ~ X1 + X2 + X3, data = d, ps = own)

  expect_identical(get_treat_type(W$treat), "multinomial")

  # ATE weights are the inverse of the probability of the treatment received
  expect_equal(unname(W$weights), unname(1 / own), tolerance = eps)
  expect_true(all(is.finite(W$weights) & W$weights > 0))

  # ATT leaves the focal group at 1
  W_att <- weightit(Am ~ X1 + X2 + X3, data = d, ps = own,
                    estimand = "ATT", focal = "T")
  expect_ATT_weights_okay(W_att, focal = "T", tolerance = eps)

  # A length mismatch is caught rather than recycled
  expect_error(weightit(Am ~ X1 + X2 + X3, data = d, ps = own[-1L]),
               "same number of units")
})

test_that("a supplied ps works for a continuous treatment", {
  d <- test_data

  # For a continuous treatment `ps` is the conditional mean, and the weights are the
  # ratio of the marginal to the conditional density evaluated at the observed treatment
  gps <- fitted(lm(Ac ~ X1 + X2 + X3, data = d))

  W <- weightit(Ac ~ X1 + X2 + X3, data = d, ps = gps)

  expect_identical(get_treat_type(W$treat), "continuous")
  expect_true(all(is.finite(W$weights) & W$weights > 0))

  # Identical to method = "glm", which fits exactly this model internally
  W_glm <- weightit(Ac ~ X1 + X2 + X3, data = d, method = "glm")
  expect_equal(unname(W$weights), unname(W_glm$weights), tolerance = eps)

  # `density` is honored on this path too, so it is not silently ignored
  W_t <- weightit(Ac ~ X1 + X2 + X3, data = d, ps = gps, density = "dt_4")
  expect_not_equal(unname(W_t$weights), unname(W$weights))
  expect_true(all(is.finite(W_t$weights) & W_t$weights > 0))
})

test_that("a supplied ps is used as the probability of being censored", {
  d <- test_data
  set.seed(123L)
  d$C <- rbinom(nrow(d), 1L, prob = plogis(-0.9 + 0.8 * d$X1))

  ps <- fitted(glm(C ~ X1 + X2, data = d, family = binomial))

  W <- weightit(.cens(C) ~ X1 + X2, data = d, ps = ps)

  expect_identical(get_treat_type(W$treat), "censoring")
  expect_equal(unname(W$weights), unname((1 - d$C) / (1 - ps)), tolerance = eps)
})

test_that("ps and method are not both used", {
  d <- test_data

  ps <- fitted(glm(A ~ X1, data = d, family = binomial))

  # A non-glm method is ignored with a warning, since the scores are already given
  expect_warning(W <- weightit(A ~ X1, data = d, ps = ps, method = "ebal"),
                 "ignored")
  expect_equal(unname(W$weights),
               unname(weightit(A ~ X1, data = d, ps = ps)$weights),
               tolerance = eps)
})
