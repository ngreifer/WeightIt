# `sbps()` implements subgroup balancing propensity score weighting: it searches over
# which covariates are balanced overall versus within moderator subgroups, and returns
# the combination with the best balance. 561 lines with one error test, so the search
# itself, the summary methods, and the M-estimation handling were all unverified -- and a
# wrongly chosen split produces plausible-looking weights with no error.

skip_on_cran()

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

test_data <- readRDS(test_path("fixtures", "test_data.rds"))

make_data <- function() {
  d <- test_data
  d$G <- factor(d$X5, labels = c("g0", "g1"))
  d
}

test_that("the three ways of specifying the moderator agree", {
  skip_if_not_installed("cobalt")

  d <- make_data()

  W1 <- weightit(A ~ X1 + X2 + X3 + G, data = d, method = "glm",
                 estimand = "ATT")
  W2 <- weightit(A ~ X1 + X2 + X3 + G, data = d, method = "glm",
                 estimand = "ATT", by = "G")

  S_obj2 <- sbps(W1, W2)
  S_str <- sbps(W1, moderator = "G")
  S_form <- sbps(W1, moderator = ~G)

  expect_s3_class(S_obj2, "weightit.sbps")
  expect_s3_class(S_obj2, "weightit")

  # Supplying the subgroup fit directly, naming the moderator, or giving it as a
  # formula all describe the same problem
  expect_equal(unname(S_str$weights), unname(S_obj2$weights), tolerance = eps)
  expect_equal(unname(S_form$weights), unname(S_obj2$weights), tolerance = eps)
})

test_that("sbps() returns weights between the overall and subgroup fits", {
  skip_if_not_installed("cobalt")

  d <- make_data()

  W1 <- weightit(A ~ X1 + X2 + X3 + G, data = d, method = "glm",
                 estimand = "ATT")
  W2 <- weightit(A ~ X1 + X2 + X3 + G, data = d, method = "glm",
                 estimand = "ATT", by = "G")

  S <- sbps(W1, W2)

  expect_length(S$weights, nrow(d))
  expect_true(all(is.finite(S$weights) & S$weights >= 0))
  expect_identical(S$estimand, "ATT")
  expect_identical(S$focal, W1$focal)

  # ATT, so the focal group is left at 1 by every candidate combination
  expect_true(all(S$weights[S$treat == S$focal] == 1))

  # The search records the moderator and, per subgroup, the share of covariates it
  # ended up balancing within subgroup rather than overall
  expect_named(S$moderator, "G")
  expect_length(S$prop.subgroup, nlevels(d$G))
  expect_true(all(S$prop.subgroup >= 0 & S$prop.subgroup <= 1))

  # And it does at least as well as the overall fit on the balance criterion it
  # optimizes, since the overall fit is one of the candidates it considers
  covs <- as.matrix(d[c("X1", "X2", "X3")])
  init <- cobalt::bal.init(covs, treat = d$A, stat = "smd.mean",
                           estimand = "ATT", focal = 1)

  expect_lte(cobalt::bal.compute(init, weights = S$weights),
             cobalt::bal.compute(init, weights = W1$weights) + eps)
})

test_that("smooth = TRUE gives a different, valid solution", {
  skip_if_not_installed("cobalt")

  d <- make_data()

  W1 <- weightit(A ~ X1 + X2 + X3 + G, data = d, method = "glm",
                 estimand = "ATT")
  W2 <- weightit(A ~ X1 + X2 + X3 + G, data = d, method = "glm",
                 estimand = "ATT", by = "G")

  S <- sbps(W1, W2)
  S_sm <- sbps(W1, W2, smooth = TRUE)

  expect_s3_class(S_sm, "weightit.sbps")
  expect_true(all(is.finite(S_sm$weights) & S_sm$weights >= 0))
  expect_not_equal(unname(S_sm$weights), unname(S$weights))

  # Smoothing interpolates the propensity scores, so it needs them
  W_ebal <- weightit(A ~ X1 + X2 + X3 + G, data = d, method = "ebal",
                     estimand = "ATT")
  expect_null(W_ebal$ps)
  expect_error(sbps(W_ebal, moderator = ~G, smooth = TRUE),
               "propensity score")
})

test_that("sbps() strips the M-estimation parts", {
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  d <- make_data()

  W1 <- weightit(A ~ X1 + X2 + X3 + G, data = d, method = "glm",
                 estimand = "ATT")

  expect_false(is_null(attr(W1, "Mparts", exact = TRUE)))

  S <- sbps(W1, moderator = ~G)

  # The subgroup search is not a smooth function of the propensity score parameters, so
  # the parts cannot be carried through; they are dropped rather than left stale
  expect_null(attr(S, "Mparts", exact = TRUE))
  expect_null(attr(S, "Mparts.list", exact = TRUE))

  # An outcome model then treats the weights as fixed rather than erroring
  fit <- glm_weightit(Y_B ~ A, data = d, weightit = S, family = binomial)
  expect_true(all(is.finite(sqrt(diag(vcov(fit))))))
  expect_error(glm_weightit(Y_B ~ A, data = d, weightit = S,
                            family = binomial, vcov = "asympt"))
})

test_that("print() and summary() work on an sbps object", {
  skip_if_not_installed("cobalt")

  d <- make_data()

  W1 <- weightit(A ~ X1 + X2 + X3 + G, data = d, method = "glm",
                 estimand = "ATT")

  S <- sbps(W1, moderator = ~G)

  expect_output(print(S), "weightit.sbps object")

  s <- summary(S)
  expect_s3_class(s, "summary.weightit.sbps")

  # One summary per subgroup
  expect_output(print(s), "g0")
  expect_output(print(s), "g1")
})

test_that("sbps() rejects what it cannot handle", {
  d <- make_data()

  W1 <- weightit(A ~ X1 + X2 + X3 + G, data = d, method = "glm",
                 estimand = "ATT")

  # One of `obj2` or `moderator` is required
  expect_error(sbps(W1), "must be specified")

  # `obj2` has to be a weightit object. The class check currently runs after `obj2` is
  # first indexed, so a non-list gives an internal subscript error rather than the
  # intended message; asserting only that it errors so this test does not pin that in
  # place.
  expect_error(sbps(W1, obj2 = "a"))

  # A moderator that cannot be found errors rather than proceeding. The message is
  # currently an internal one from the formula processing rather than a clear
  # "variable not found"; asserting only that it errors so this test does not pin the
  # wording in place.
  expect_error(sbps(W1, moderator = ~nope, data = d))

  # A continuous moderator has subgroups without both treatment levels
  expect_error(sbps(W1, moderator = ~X1, data = d), "treatment levels")

  # Censoring weights have no subgroup balancing objective
  set.seed(123L)
  d$C <- rbinom(nrow(d), 1L, prob = plogis(-0.9 + 0.8 * d$X1))
  Wc <- weightit(.cens(C) ~ X1 + X2, data = d, method = "glm")
  expect_error(sbps(Wc, moderator = ~G), "cannot be used with censoring")
})
