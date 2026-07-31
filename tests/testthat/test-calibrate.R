# `calibrate()` refits the propensity score on itself -- Platt scaling or isotonic
# regression -- and recomputes the weights. It is exported and had no tests at all, which
# matters because it produces numbers directly: a sign or weighting error in the
# calibration is invisible without assertions.

skip_on_cran()

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

test_data <- readRDS(test_path("fixtures", "test_data.rds"))

ps_fit <- function(d) {
  fitted(glm(A ~ X1 + X2 + X3, data = d, family = binomial))
}

test_that("both methods return calibrated probabilities", {
  d <- test_data
  ps <- ps_fit(d)

  for (method in c("platt", "isoreg")) {
    p <- calibrate(ps, treat = d$A, method = method)

    expect_length(p, nrow(d))
    expect_false(anyNA(p))
    expect_named(p, names(ps))
    expect_true(all(p >= 0 & p <= 1),
                label = sprintf("calibrated scores in [0, 1] for %s", method))

    # The defining property of calibration: the mean predicted probability equals the
    # observed treatment prevalence. Both methods hit this exactly, so it is a tight
    # check on the weighting inside each.
    expect_equal(mean(p), mean(d$A), tolerance = eps,
                 label = sprintf("mean calibrated score for %s", method))
  }
})

test_that("Platt scaling is a smooth monotone transformation", {
  d <- test_data
  ps <- ps_fit(d)

  p <- calibrate(ps, treat = d$A, method = "platt")

  # A single-predictor logistic fit, so strictly inside (0, 1) and order-preserving
  expect_true(all(p > 0 & p < 1))
  expect_identical(order(p), order(ps))
  expect_identical(nunique(p), nunique(ps))

  # It is a fit on the probability scale rather than the logit scale, so it is not the
  # identity even on an already-correctly-specified score
  expect_not_equal(unname(p), unname(ps))
})

test_that("isotonic calibration is a monotone step function within treatment group", {
  d <- test_data
  ps <- ps_fit(d)

  p <- calibrate(ps, treat = d$A, method = "isoreg")

  # Far fewer distinct values than the input: isotonic regression pools into blocks
  expect_lt(nunique(p), nunique(ps) / 10)

  # The treated and control scores are calibrated separately and then blended, so
  # monotonicity holds within each group rather than overall
  for (a in 0:1) {
    in_a <- d$A == a
    expect_false(is.unsorted(p[in_a][order(ps[in_a])]),
                 label = sprintf("monotone within A = %s", a))
  }

  # Isotonic calibration can reach exactly 0 or 1, but never in the group that would
  # make the weight infinite: no treated unit at 0 and no control at 1
  expect_identical(sum(p == 0 & d$A == 1), 0L)
  expect_identical(sum(p == 1 & d$A == 0), 0L)
  expect_true(all(is.finite(get_w_from_ps(p, d$A))))
})

test_that("s.weights are honored", {
  d <- test_data
  ps <- ps_fit(d)

  for (method in c("platt", "isoreg")) {
    p_unw <- calibrate(ps, treat = d$A, method = method)
    p_sw <- calibrate(ps, treat = d$A, s.weights = d$SW, method = method)

    expect_not_equal(unname(p_sw), unname(p_unw))
    expect_true(all(p_sw >= 0 & p_sw <= 1))

    # Calibration now targets the sampling-weighted prevalence
    expect_equal(weighted.mean(p_sw, d$SW), weighted.mean(d$A, d$SW),
                 tolerance = eps,
                 label = sprintf("weighted mean calibrated score for %s", method))

    # `s.weights` may also name a variable in `data`
    p_str <- calibrate(ps, treat = d$A, s.weights = "SW", data = d,
                       method = method)
    expect_equal(unname(p_str), unname(p_sw), tolerance = eps)
  }
})

test_that("calibrate.weightit() replaces the scores and the weights", {
  d <- test_data

  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATE")

  Wc <- calibrate(W)

  expect_s3_class(Wc, "weightit")
  expect_not_equal(unname(Wc$ps), unname(W$ps))
  expect_not_equal(unname(Wc$weights), unname(W$weights))
  expect_true(all(is.finite(Wc$weights) & Wc$weights > 0))

  # The new weights are exactly those implied by the new scores
  expect_equal(unname(Wc$weights),
               unname(get_w_from_ps(Wc$ps, Wc$treat, estimand = "ATE")),
               tolerance = eps)

  expect_identical(attr(Wc, "calibrate")$method, "platt")

  # Calibrating breaks the M-estimation chain, so the parts are dropped rather than
  # left stale: standard errors then treat the weights as fixed
  expect_false(is_null(attr(W, "Mparts", exact = TRUE)))
  expect_null(attr(Wc, "Mparts", exact = TRUE))

  expect_no_error(vcov(glm_weightit(Y_B ~ A, data = d, weightit = Wc,
                                    family = binomial)))

  # method = "isoreg" is carried through too
  Wi <- calibrate(W, method = "isoreg")
  expect_identical(attr(Wi, "calibrate")$method, "isoreg")
  expect_not_equal(unname(Wi$weights), unname(Wc$weights))
  expect_true(all(is.finite(Wi$weights) & Wi$weights > 0))
})

test_that("calibrate() honors the estimand of the object it is given", {
  d <- test_data

  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATT")

  Wc <- calibrate(W)

  # ATT weights stay 1 in the focal group after calibration
  expect_ATT_weights_okay(Wc, tolerance = eps)
})

test_that("calibrate() rejects what it cannot handle", {
  d <- test_data

  # Only binary treatments, on the default method
  expect_error(calibrate(runif(nrow(d)), treat = d$Am), "binary")
  expect_error(calibrate(runif(nrow(d)), treat = d$Ac), "binary")

  # `treat` is required
  expect_error(calibrate(runif(nrow(d))))

  # An unrecognized method
  expect_error(calibrate(ps_fit(d), treat = d$A, method = "nope"))

  # On a `weightit` object, estimated propensity scores are required; multi-category and
  # continuous glm fits store none, so that is the error they raise
  for (f in list(Am ~ X1 + X2, Ac ~ X1 + X2)) {
    W <- weightit(f, data = d, method = "glm")
    expect_null(W$ps)
    expect_error(calibrate(W), "propensity score")
  }

  W_ebal <- weightit(A ~ X1 + X2, data = d, method = "ebal")
  expect_null(W_ebal$ps)
  expect_error(calibrate(W_ebal), "propensity score")
})
