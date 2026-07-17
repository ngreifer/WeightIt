# msmdata (data("msmdata", package = "WeightIt")) has 7500 units, baseline
# covariates X1_0 (count)/X2_0 (binary), then at 3 time points a binary
# treatment A_t and, for t = 1, 2 only, covariates X1_t/X2_t (none after the
# last treatment), and a binary outcome Y_B. All three treatment columns are
# binary in the shipped data, so a couple of scenarios below build small
# local, in-test-only copies with a continuous or multi-category time point
# to exercise that code path; data/msmdata.rda and data-raw/ are never
# touched.

skip_on_cran()

data("msmdata", package = "WeightIt")

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

# Formula list used throughout: treatment on the LHS, cumulative covariate
# (and, from the second time point on, prior treatment) history on the RHS.
msm_formulas <- list(
  A_1 ~ X1_0 + X2_0,
  A_2 ~ X1_1 + X2_1 + A_1,
  A_3 ~ X1_2 + X2_2 + A_2
)

test_that("Baseline smoke test with method = 'glm'", {
  expect_no_error({
    W <- weightitMSM(msm_formulas, data = msmdata, method = "glm")
  })

  expect_s3_class(W, "weightitMSM")
  expect_length(W$treat.list, 3L)
  expect_length(W$weights, nrow(msmdata))
  expect_true(all(is.finite(W$weights)))
  expect_true(all(W$weights > 0))
  expect_identical(W$estimand, "ATE")
  expect_identical(names(W$treat.list), c("A_1", "A_2", "A_3"))
})

test_that("msm_valid methods fit without error: glm, gbm, cbps, ipt, super, bart", {
  skip_if_not_installed("cobalt")

  # glm, gbm, ipt, super, and bart are not `msm_method_available`, so they
  # fit a separate model at each time point and multiply the resulting
  # weights together (the `is.MSM.method = FALSE` path in weightitMSM()).
  # cbps is the only method with `msm_method_available = TRUE`, so by
  # default (is.MSM.method left unspecified) it fits all three time points
  # jointly via weightitMSM.fit()/weightitMSM2cbps() instead.

  expect_no_error({
    W_glm <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                         include.obj = TRUE)
  })
  expect_true(all(is.finite(W_glm$weights) & W_glm$weights > 0))
  # per-time-point path: one fit object per formula
  expect_length(W_glm$obj, length(msm_formulas))

  skip_if_not_installed("gbm")
  set.seed(123)
  expect_no_error({
    W_gbm <- weightitMSM(msm_formulas, data = msmdata, method = "gbm",
                         n.trees = 200, criterion = "smd.max")
  })
  expect_true(all(is.finite(W_gbm$weights) & W_gbm$weights > 0))

  expect_no_error({
    W_cbps <- weightitMSM(msm_formulas, data = msmdata, method = "cbps",
                          include.obj = TRUE)
  })
  expect_true(all(is.finite(W_cbps$weights) & W_cbps$weights > 0))

  # cbps used the joint-fit path by default (is.MSM.method resolves to TRUE
  # since msm_method_available = TRUE and is.MSM.method was left unspecified):
  # a single fit object for all 3 time points together, not a length-3 list
  # of per-time-point fits like W_glm$obj above.
  expect_false(is.list(W_cbps$obj) && length(W_cbps$obj) == length(msm_formulas))

  expect_no_error({
    W_ipt <- weightitMSM(msm_formulas, data = msmdata, method = "ipt")
  })
  expect_true(all(is.finite(W_ipt$weights) & W_ipt$weights > 0))

  skip_if_not_installed("SuperLearner")
  set.seed(123)
  expect_no_error({
    W_super <- weightitMSM(msm_formulas, data = msmdata, method = "super",
                           SL.library = c("SL.mean", "SL.glm", "SL.step.interaction"))
  })
  expect_true(all(is.finite(W_super$weights) & W_super$weights > 0))

  skip_if_not_installed("dbarts")
  set.seed(123)
  expect_no_error({
    W_bart <- weightitMSM(msm_formulas, data = msmdata, method = "bart",
                          n.trees = 20L, n.threads = 1L, seed = 123)
  })
  expect_true(all(is.finite(W_bart$weights) & W_bart$weights > 0))

  # Sanity check that these six methods don't all produce identical weights
  wts <- list(glm = W_glm$weights, gbm = W_gbm$weights, cbps = W_cbps$weights,
             ipt = W_ipt$weights, super = W_super$weights, bart = W_bart$weights)
  for (nm in setdiff(names(wts), "glm")) {
    expect_not_equal(wts[[nm]], wts[["glm"]], expected.label = "glm weights")
  }
})

test_that("msm_valid = FALSE methods error without weightit.force and succeed with it: ebal", {
  # ebal has no extra package dependency, so it can always be tested.
  expect_error({
    weightitMSM(msm_formulas, data = msmdata, method = "ebal")
  }, "has not been validated")

  expect_no_error({
    W <- weightitMSM(msm_formulas, data = msmdata, method = "ebal",
                     weightit.force = TRUE)
  })
  expect_true(all(is.finite(W$weights) & W$weights > 0))
})

test_that("msm_valid = FALSE methods error without weightit.force and succeed with it: optweight", {
  skip_if_not_installed("optweight")

  expect_error({
    weightitMSM(msm_formulas, data = msmdata, method = "optweight")
  }, "has not been validated")

  expect_no_error({
    W <- weightitMSM(msm_formulas, data = msmdata, method = "optweight",
                     weightit.force = TRUE)
  })
  expect_true(all(is.finite(W$weights) & W$weights > 0))
})

test_that("stabilize = TRUE changes weights for a stabilize_ok method (glm)", {
  W0 <- weightitMSM(msm_formulas, data = msmdata, method = "glm")

  expect_no_error({
    W1 <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                      stabilize = TRUE)
  })

  expect_not_equal(W1$weights, W0$weights)
  expect_true(all(is.finite(W1$weights) & W1$weights > 0))
  expect_false(is_null(W1$stabilization))
})

test_that("stabilize = TRUE is a no-op (with a warning) for stabilize_ok = FALSE methods: cbps", {
  # method = "cbps" is WeightIt's own implementation and does not depend on
  # the CBPS package (that's only needed for method = "npcbps"), so no
  # skip_if_not_installed("CBPS") guard is needed here.
  W0 <- weightitMSM(msm_formulas, data = msmdata, method = "cbps")

  expect_warning({
    W1 <- weightitMSM(msm_formulas, data = msmdata, method = "cbps",
                      stabilize = TRUE)
  }, "stabilize.*cannot be used")

  expect_equal(W1$weights, W0$weights, tolerance = eps)
})

test_that("stabilize = TRUE is a no-op (with a warning) for stabilize_ok = FALSE methods: ipt", {
  W0 <- weightitMSM(msm_formulas, data = msmdata, method = "ipt")

  expect_warning({
    W1 <- weightitMSM(msm_formulas, data = msmdata, method = "ipt",
                      stabilize = TRUE)
  }, "stabilize.*cannot be used")

  expect_equal(W1$weights, W0$weights, tolerance = eps)
})

test_that("is.MSM.method = FALSE changes weights for a msm_method_available method: cbps", {
  # By default (is.MSM.method left unspecified), cbps resolves to TRUE and
  # fits all 3 time points jointly (see the "msm_valid methods fit without
  # error" test above). Explicitly setting is.MSM.method = FALSE switches to
  # fitting a separate CBPS() model at each time point and multiplying the
  # weights together instead -- a message is emitted since this deviates from
  # cbps's default, and the resulting weights and fit-object structure should
  # differ from the default (TRUE) case.
  W0 <- weightitMSM(msm_formulas, data = msmdata, method = "cbps",
                    include.obj = TRUE)
  expect_false(is.list(W0$obj) && length(W0$obj) == length(msm_formulas))

  expect_message({
    W1 <- weightitMSM(msm_formulas, data = msmdata, method = "cbps",
                      is.MSM.method = FALSE, include.obj = TRUE)
  }, "single model.*model for each time point", ignore.case = TRUE)

  expect_true(all(is.finite(W1$weights) & W1$weights > 0))
  expect_not_equal(W1$weights, W0$weights)

  # is.MSM.method = FALSE uses the per-time-point path: one fit object per
  # formula, just like the non-msm_method_available methods (e.g., glm).
  expect_true(is.list(W1$obj) && length(W1$obj) == length(msm_formulas))
})

test_that("custom num.formula works with method = 'glm'", {
  W0 <- weightitMSM(msm_formulas, data = msmdata, method = "glm")

  # A single one-sided formula applied (with the appropriate saturated
  # treatment-history terms added by weightitMSM()) at every time point.
  expect_no_error({
    W1 <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                      num.formula = ~X1_0, stabilize = TRUE)
  })
  expect_true(all(is.finite(W1$weights) & W1$weights > 0))
  expect_not_equal(W1$weights, W0$weights)
  expect_false(is_null(W1$stabilization))

  # num.formula's terms must appear in the stabilization formula at every
  # time point, not just the first -- previously, a bug in the string
  # construction for time points after the first silently dropped
  # num.formula's terms entirely (see R/weightitMSM.R's handling of the
  # `else` branch under `rlang::is_formula(num.formula)`).
  for (stab_f in W1$stabilization) {
    expect_true("X1_0" %in% all.vars(stab_f))
  }

  # A list of one-sided formulas, one per time point.
  expect_no_error({
    W2 <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                      num.formula = list(~1, ~X1_0, ~X1_0 + X1_1),
                      stabilize = TRUE)
  })
  expect_true(all(is.finite(W2$weights) & W2$weights > 0))
  expect_not_equal(W2$weights, W0$weights)

  # Supplying num.formula implies stabilize = TRUE even if not set explicitly
  expect_equal(W1$weights,
              weightitMSM(msm_formulas, data = msmdata, method = "glm",
                         stabilize = TRUE, num.formula = ~X1_0)$weights,
              tolerance = eps)
})

test_that("`by` is accepted and retains combined Mparts.list across strata and time points", {
  msmdata_by <- msmdata
  msmdata_by$grp <- factor(rep(c("a", "b"), length.out = nrow(msmdata_by)))

  expect_no_error({
    W <- weightitMSM(msm_formulas, data = msmdata_by, method = "glm",
                     by = "grp")
  })

  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_false(is_null(W$by))
  expect_identical(names(W$by), "grp")

  # With `by`, the per-stratum, per-time-point Mparts are combined into a single
  # "Mparts.list": one part per (time point x by-stratum). The formula list has
  # 3 time points and `grp` has 2 levels, so 6 parts total. (Equivalence to the
  # fully-interacted specification is checked in test-by_mest.R.)
  W_nb <- weightitMSM(msm_formulas, data = msmdata, method = "glm")
  n_nb <- length(attr(W_nb, "Mparts.list", exact = TRUE))
  expect_false(is_null(attr(W, "Mparts.list", exact = TRUE)))
  expect_identical(length(attr(W, "Mparts.list", exact = TRUE)),
                   n_nb * nlevels(msmdata_by$grp))
})

test_that("s.weights works with method = 'glm'", {
  msmdata_sw <- msmdata
  set.seed(123)
  msmdata_sw$sw <- runif(nrow(msmdata_sw), 0.5, 1.5)

  W0 <- weightitMSM(msm_formulas, data = msmdata_sw, method = "glm")

  expect_no_error({
    W1 <- weightitMSM(msm_formulas, data = msmdata_sw, method = "glm",
                      s.weights = "sw")
  })

  expect_true(all(is.finite(W1$weights) & W1$weights > 0))
  expect_equal(unname(W1$s.weights), msmdata_sw$sw)
  expect_not_equal(W1$weights, W0$weights)
})

test_that("s.weights errors for a s.weights_ok = FALSE method: bart", {
  skip_if_not_installed("dbarts")

  msmdata_sw <- msmdata
  set.seed(123)
  msmdata_sw$sw <- runif(nrow(msmdata_sw), 0.5, 1.5)

  # `.weightit_methods[["bart"]]$s.weights_ok` is FALSE, and
  # `.check_method_s.weights()` is called directly inside weightitMSM(), so
  # supplying non-constant s.weights errors immediately -- it is neither
  # silently ignored nor merely warned about.
  expect_error({
    weightitMSM(msm_formulas, data = msmdata_sw, method = "bart",
               s.weights = "sw", n.trees = 20L, n.threads = 1L, seed = 123)
  }, "sampling weights cannot be used", ignore.case = TRUE)
})

test_that("Non-binary treatment time point: continuous A_2, method = 'glm'", {
  msmdata_cont <- msmdata
  set.seed(123)
  msmdata_cont$A_2 <- rnorm(nrow(msmdata_cont))

  # A_2's formula must be updated for a linear/continuous treatment model;
  # keep the same covariate/treatment-history structure.
  formulas_cont <- list(
    A_1 ~ X1_0 + X2_0,
    A_2 ~ X1_1 + X2_1 + A_1,
    A_3 ~ X1_2 + X2_2 + A_2
  )

  expect_no_error({
    W <- weightitMSM(formulas_cont, data = msmdata_cont, method = "glm")
  })

  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_identical(unname(vapply(W$treat.list, get_treat_type, character(1L))),
                   c("binary", "continuous", "binary"))
})

test_that("Non-binary treatment time point: continuous A_2, method = 'cbps'", {
  # weightit2cbps.R explicitly claims: "Any combination of treatment types is
  # supported" for longitudinal (MSM) treatments. This is a direct
  # doc-vs-behavior check using a mixed binary/continuous/binary sequence.
  # method = "cbps" is WeightIt's own implementation and does not depend on
  # the CBPS package (that's only needed for method = "npcbps"), so no
  # skip_if_not_installed("CBPS") guard is needed here.
  msmdata_cont <- msmdata
  set.seed(123)
  msmdata_cont$A_2 <- rnorm(nrow(msmdata_cont))

  formulas_cont <- list(
    A_1 ~ X1_0 + X2_0,
    A_2 ~ X1_1 + X2_1 + A_1,
    A_3 ~ X1_2 + X2_2 + A_2
  )

  expect_no_error({
    W <- weightitMSM(formulas_cont, data = msmdata_cont, method = "cbps")
  })

  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_identical(unname(vapply(W$treat.list, get_treat_type, character(1L))),
                   c("binary", "continuous", "binary"))
})
