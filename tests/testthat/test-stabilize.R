# Stabilization, i.e., replacing the numerator of the weights with a fitted model
# rather than 1. Covers `stabilize` in `weightit()` (a flag or a one-sided formula,
# possibly with random effects) and in `weightitMSM()` (where it additionally
# accepts anything `num.formula` does and supersedes it), and the `stabilization`
# component of the output, which stores the numerator formula that was used.
#
# The `by`- and censoring-specific compositions live in test-by_mest.R and
# test-censoring.R; the M-estimation checks live in test-method_glm.R.

skip_on_cran()

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

test_data <- readRDS(test_path("fixtures", "test_data.rds"))

data("msmdata", package = "WeightIt")

msm_formulas <- list(
  A_1 ~ X1_0 + X2_0,
  A_2 ~ X1_1 + X2_1 + A_1,
  A_3 ~ X1_2 + X2_2 + A_2
)

# msmdata has no grouping variable, so one is built here (in-test only) with a real
# random-intercept signal: the treatments are redrawn so that cluster membership
# actually shifts the treatment probability. Without that the variance component is
# estimated at the boundary and the multilevel fit is indistinguishable from the
# fixed-effects one, which would make the tests below vacuous.
msmdata_re <- local({
  set.seed(123)

  d <- msmdata
  nclus <- 40L
  d$clus <- factor(sample(seq_len(nclus), nrow(d), replace = TRUE))

  u <- rnorm(nclus, 0, 1.5)

  for (a in c("A_1", "A_2", "A_3")) {
    d[[a]] <- rbinom(nrow(d), 1L, plogis(-0.5 + 0.2 * d$X1_0 + u[as.integer(d$clus)]))
  }

  d
})

test_that("weightit(): stabilize = TRUE is the marginal numerator", {
  W0 <- weightit(A ~ X1 + X2 + X3, data = test_data, method = "glm",
                 estimand = "ATE")

  expect_no_condition({
    W <- weightit(A ~ X1 + X2 + X3, data = test_data, method = "glm",
                  estimand = "ATE", stabilize = TRUE)
  })

  #`stabilize = TRUE` is shorthand for an intercept-only numerator
  expect_equal(W$weights,
               weightit(A ~ X1 + X2 + X3, data = test_data, method = "glm",
                        estimand = "ATE", stabilize = ~1)$weights,
               tolerance = eps)

  #P(A = a)/P(A = a | X), computed by hand from the unstabilized weights
  p <- mean(test_data$A)
  expect_equal(W$weights,
               W0$weights * ifelse(test_data$A == 1, p, 1 - p),
               tolerance = eps)

  #Stabilization only rescales within treatment group; balance is untouched
  expect_equal(unname(W$weights / W0$weights),
               ifelse(test_data$A == 1, p, 1 - p),
               tolerance = eps)
})

test_that("weightit(): the stabilization formula is stored in `stabilization`", {
  W <- weightit(A ~ X1 + X2 + X3, data = test_data, method = "glm",
                stabilize = ~X5)

  #One-sided, and the `stabilize` component of pre-2.0.1 objects is gone
  expect_true(rlang::is_formula(W$stabilization, lhs = FALSE))
  expect_identical(deparse1(W$stabilization), "~X5")
  expect_null(W$stabilize)

  expect_identical(deparse1(weightit(A ~ X1 + X2, data = test_data,
                                     method = "glm", stabilize = TRUE)$stabilization),
                   "~1")

  #Absent when not stabilized
  expect_null(weightit(A ~ X1 + X2, data = test_data, method = "glm")$stabilization)

  #Printed, and reported by summary() as the mean of the weights
  expect_output(print(W), "stabilized; stabilization factors:")
  expect_output(print(W), "X5")

  #The stabilization block ends its line, so whatever comes after it starts on
  #its own
  out <- capture.output(print(trim(W, at = 0.9)))
  expect_identical(out[length(out)], " - weights trimmed at 90%")

  s <- summary(W)
  expect_length(s$weight.mean, 2L)
  expect_equal(unname(s$weight.mean["Treated"]),
               mean(W$weights[test_data$A == 1]),
               tolerance = eps)
})

test_that("weightit(): the stabilization formula can contain random effects", {
  skip_if_not_installed("lme4")

  W0 <- weightit(A ~ X2 + X3 + X4 + X7, data = test_data, method = "glm")

  expect_no_condition({
    W <- weightit(A ~ X2 + X3 + X4 + X7, data = test_data, method = "glm",
                  stabilize = ~ (1 | cluster))
  })

  expect_identical(deparse1(W$stabilization), "~(1 | cluster)")
  expect_true(all(is.finite(W$weights) & W$weights > 0))

  #The numerator is `lme4::glmer(A ~ 1 + (1 | cluster))`, fit the same way
  #`.make_re_data_formula()` sets it up, and the stabilized weights are the
  #unstabilized ones divided by the numerator's own weights
  num_fit <- lme4::glmer(A ~ 1 + (1 | cluster), data = test_data,
                         family = binomial("logit"))
  p_num <- as.numeric(predict(num_fit, type = "response"))

  expect_equal(W$weights,
               W0$weights / get_w_from_ps(p_num, test_data$A, estimand = "ATE"),
               tolerance = eps)

  #The random intercept is actually estimated, so this differs from `~1`
  expect_gt(as.numeric(lme4::VarCorr(num_fit)[["cluster"]]), 0.1)
  expect_not_equal(W$weights,
                   weightit(A ~ X2 + X3 + X4 + X7, data = test_data,
                            method = "glm", stabilize = ~1)$weights)

  #No M-estimation parts once a mixed model is in the stack
  expect_null(attr(W, "Mparts"))
  expect_null(attr(W, "Mparts.list"))
})

test_that("weightit(): random effects in the denominator and the numerator", {
  skip_if_not_installed("lme4")

  expect_no_condition({
    W <- weightit(A ~ X2 + X3 + X7 + (1 | cluster), data = test_data,
                  method = "glm", stabilize = ~ (1 | cluster))
  })

  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_identical(deparse1(W$stabilization), "~(1 | cluster)")

  #The denominator random effects do not leak into the numerator formula
  expect_not_equal(W$weights,
                   weightit(A ~ X2 + X3 + X7 + (1 | cluster), data = test_data,
                            method = "glm", stabilize = TRUE)$weights)
})

test_that("weightit(): stabilize is restricted to the ATE", {
  skip_if_not_installed("patrick")

  #For an estimand other than the ATE the numerator is not the marginal probability
  #of the observed treatment, so dividing by it rescales the treatment groups
  #relative to each other rather than leaving the estimand alone.
  patrick::with_parameters_test_that(
    "estimand = {estimand}, stabilize = {stab_name}",
    {
      expect_error(weightit(A ~ X1 + X2 + X3, data = test_data, method = "glm",
                            estimand = estimand, stabilize = stab),
                   'estimand = "ATE"', fixed = TRUE)
    },
    .cases = merge(patrick::cases(list(stab_name = "TRUE", stab = TRUE),
                                  list(stab_name = "~X5", stab = ~X5)),
                   data.frame(estimand = c("ATT", "ATC", "ATO", "ATM")),
                   by = NULL)
  )
})

test_that("weightit(): stabilize's ATE restriction for other treatment types", {
  #Multi-category treatments too
  expect_error(weightit(Am ~ X1 + X2, data = test_data, method = "glm",
                        estimand = "ATT", focal = "T", stabilize = TRUE),
               'estimand = "ATE"', fixed = TRUE)

  #A continuous treatment has no estimand, so it is exempt
  expect_no_condition({
    weightit(Ac ~ X1 + X2, data = test_data, method = "glm", stabilize = TRUE)
  })

  #So is a censoring model
  d_cens <- test_data
  set.seed(1234)
  d_cens$C <- rbinom(nrow(d_cens), 1L, 0.2)

  expect_no_condition({
    weightit(.cens(C) ~ X1 + X2, data = d_cens, method = "glm", stabilize = TRUE)
  })
})

test_that("weightit(): stabilize errors and warnings", {
  #Two-sided formula
  expect_error(weightit(A ~ X1, data = test_data, method = "glm",
                        stabilize = A ~ X5),
               "one-sided formula")

  #Variable not in `data`
  expect_error(weightit(A ~ X1, data = test_data, method = "glm",
                        stabilize = ~nonexistent),
               "must be variables in")

  #Not a flag or a formula
  expect_error(weightit(A ~ X1, data = test_data, method = "glm",
                        stabilize = "yes"),
               "logical value")

  #Ignored, with a warning, for methods that cannot use it
  expect_warning({
    W <- weightit(A ~ X1 + X2, data = test_data, method = "ebal",
                  stabilize = ~X5)
  }, "cannot be used with")

  expect_null(W$stabilization)
})

test_that("weightitMSM(): stabilize accepts a formula and supersedes num.formula", {
  Wf <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                    stabilize = ~X2_0)

  Wn <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                    stabilize = TRUE, num.formula = ~X2_0)

  expect_equal(Wf$weights, Wn$weights, tolerance = eps)
  expect_identical(lapply(Wf$stabilization, deparse1),
                   lapply(Wn$stabilization, deparse1))

  #The stabilization factors are added to the saturated prior-treatment model
  expect_identical(vapply(Wf$stabilization, deparse1, character(1L)),
                   c("~X2_0", "~X2_0 + A_1", "~X2_0 + A_1 + A_2 + A_1:A_2"))

  #`num.formula` is ignored when `stabilize` is a formula
  Wo <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                    stabilize = ~X2_0, num.formula = ~X1_0)
  expect_equal(Wo$weights, Wf$weights, tolerance = eps)

  #`num.formula` alone still turns stabilization on, with a message
  expect_message({
    Wm <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                      num.formula = ~X2_0)
  }, "num.formula")
  expect_equal(Wm$weights, Wf$weights, tolerance = eps)
})

test_that("weightitMSM(): stabilize accepts a list of formulas", {
  nf <- list(~1, ~X2_0, ~ X2_0 + X1_0)

  Wl <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                    stabilize = nf)

  Wn <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                    stabilize = TRUE, num.formula = nf)

  expect_equal(Wl$weights, Wn$weights, tolerance = eps)

  #A list is used as supplied: no saturated prior-treatment terms are added
  expect_identical(vapply(Wl$stabilization, deparse1, character(1L)),
                   c("~1", "~X2_0", "~X2_0 + X1_0"))

  #One entry per entry of `formula.list`
  expect_error(weightitMSM(msm_formulas, data = msmdata, method = "glm",
                           stabilize = list(~1, ~X2_0)),
               "as many entries as")

  #Entries must be one-sided formulas
  expect_error(weightitMSM(msm_formulas, data = msmdata, method = "glm",
                           stabilize = list(~1, ~X2_0, "X1_0")),
               "list thereof")

  #Neither a flag, nor a formula, nor a list
  expect_error(weightitMSM(msm_formulas, data = msmdata, method = "glm",
                           stabilize = "yes"),
               "logical value")
})

test_that("weightitMSM(): stabilize = TRUE is the saturated prior-treatment model", {
  W <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                   stabilize = TRUE)

  expect_identical(vapply(W$stabilization, deparse1, character(1L)),
                   c("~1", "~A_1", "~A_1 + A_2 + A_1:A_2"))

  #Equivalent to dividing by the weights from the numerator models fit by hand
  W0 <- weightitMSM(msm_formulas, data = msmdata, method = "glm")

  num_w <- Reduce("*", lapply(seq_along(msm_formulas), function(i) {
    treat_i <- msmdata[[all.vars(rlang::f_lhs(msm_formulas[[i]]))]]
    f_i <- update(msm_formulas[[i]], W$stabilization[[i]])
    p_i <- fitted(glm(f_i, data = msmdata, family = binomial()))
    get_w_from_ps(p_i, treat_i, estimand = "ATE")
  }))

  expect_equal(W$weights, W0$weights / num_w, tolerance = eps)
})

test_that("weightitMSM(): stabilization formulas can contain random effects", {
  skip_if_not_installed("lme4")

  W0 <- weightitMSM(msm_formulas, data = msmdata_re, method = "glm")

  #No warning: the bar term must not be evaluated as an expression when the
  #variables in `stabilize` are checked against `data`
  expect_no_condition({
    W <- weightitMSM(msm_formulas, data = msmdata_re, method = "glm",
                     stabilize = ~ (1 | clus))
  })

  expect_identical(vapply(W$stabilization, deparse1, character(1L)),
                   c("~(1 | clus)", "~(1 | clus) + A_1",
                     "~(1 | clus) + A_1 + A_2 + A_1:A_2"))

  expect_true(all(is.finite(W$weights) & W$weights > 0))
  expect_not_equal(W$weights, W0$weights)

  #Differs from the same numerators without the random intercept
  expect_not_equal(W$weights,
                   weightitMSM(msm_formulas, data = msmdata_re, method = "glm",
                               stabilize = TRUE)$weights)

  #The first time point's numerator is `glmer(A_1 ~ 1 + (1 | clus))`
  num_fit <- lme4::glmer(A_1 ~ 1 + (1 | clus), data = msmdata_re,
                         family = binomial("logit"))
  expect_gt(as.numeric(lme4::VarCorr(num_fit)[["clus"]]), 0.1)

  #A list may mix random-effects and fixed-effects entries
  expect_no_condition({
    Wl <- weightitMSM(msm_formulas, data = msmdata_re, method = "glm",
                      stabilize = list(~1, ~ (1 | clus), ~X2_0))
  })

  expect_identical(vapply(Wl$stabilization, deparse1, character(1L)),
                   c("~1", "~(1 | clus)", "~X2_0"))

  #A grouping variable that does not exist is caught by name
  expect_error(weightitMSM(msm_formulas, data = msmdata_re, method = "glm",
                           stabilize = ~ (1 | nope)),
               "nope")
})

test_that("a numerator of exactly 1 is not reported as stabilization", {
  #A continuous treatment's weights are already f_A(a)/f_{A|X}(a), so a `~1` numerator
  #models the density they divide by: it estimates nothing and comes to exactly 1. An
  #object whose numerator applied nothing is not stabilized and does not say it is.
  W_F <- weightit(Ac ~ X1 + X2, data = test_data, method = "glm")

  for (stab in list(TRUE, ~1)) {
    W <- weightit(Ac ~ X1 + X2, data = test_data, method = "glm", stabilize = stab)

    expect_null(W$stabilization, label = deparse1(stab))
    expect_identical(unname(W$weights), unname(W_F$weights), label = deparse1(stab))
  }

  #A numerator with terms in it changes the weights, so it is reported
  W_x <- weightit(Ac ~ X1 + X2, data = test_data, method = "glm", stabilize = ~X5)

  expect_identical(deparse1(W_x$stabilization), "~X5")
  expect_false(isTRUE(all.equal(unname(W_x$weights), unname(W_F$weights))))

  #The rule is about the numerator coming to 1, not about the treatment being
  #continuous: a binary treatment's `~1` numerator is the marginal probability, which
  #is not 1, and stays
  W_b <- weightit(A ~ X1 + X2, data = test_data, method = "glm", stabilize = TRUE)

  expect_identical(deparse1(W_b$stabilization), "~1")
  expect_false(isTRUE(all.equal(
    unname(W_b$weights),
    unname(weightit(A ~ X1 + X2, data = test_data, method = "glm")$weights))))
})

test_that("weightitMSM() drops the stabilization only when every numerator is 1", {
  #Continuous treatments at every time point. At the first the numerator is `~1` and so
  #is exactly 1, but at later ones it is a model for the preceding treatments and is
  #not, which is what makes `stabilize = TRUE` meaningful for an MSM where it is inert
  #for a point treatment.
  set.seed(31)
  n <- 1000

  d <- data.frame(X1_0 = rnorm(n))
  d$Ac_1 <- rnorm(n, d$X1_0)
  d$X1_1 <- rnorm(n, 0.4 * d$Ac_1)
  d$Ac_2 <- rnorm(n, d$X1_1 + 0.5 * d$Ac_1)

  fl <- list(Ac_1 ~ X1_0, Ac_2 ~ X1_1 + Ac_1)

  W2 <- weightitMSM(fl, data = d, method = "glm", stabilize = TRUE)

  expect_length(W2$stabilization, 2L)
  expect_identical(vapply(W2$stabilization, deparse1, character(1L)),
                   c("~1", "~Ac_1"))
  expect_false(isTRUE(all.equal(
    unname(W2$weights),
    unname(weightitMSM(fl, data = d, method = "glm")$weights))))

  #With one continuous time point there is no preceding treatment to condition on, so
  #the only numerator is 1 and nothing was stabilized
  W1 <- weightitMSM(fl[1L], data = d, method = "glm", stabilize = TRUE)

  expect_null(W1$stabilization)
  expect_identical(unname(W1$weights),
                   unname(weightitMSM(fl[1L], data = d, method = "glm")$weights))
  expect_failure(expect_output(print(W1), "stabilized"))

  #A binary treatment's numerators are probabilities rather than 1, so a single time
  #point keeps its stabilization
  Wb <- weightitMSM(msm_formulas[1L], data = msmdata, method = "glm", stabilize = TRUE)

  expect_length(Wb$stabilization, 1L)
  expect_output(print(Wb), "stabilized")
})

test_that("weightitMSM(): stabilization formulas are printed", {
  W <- weightitMSM(msm_formulas, data = msmdata, method = "glm",
                   stabilize = TRUE)

  out <- capture.output(print(W))

  expect_true(any(grepl("stabilized; stabilization factors:", out, fixed = TRUE)))
  expect_true(any(grepl("+ time 1 (A_1): (none)", out, fixed = TRUE)))
  expect_true(any(grepl("+ time 2 (A_2): A_1", out, fixed = TRUE)))
  expect_true(any(grepl("+ time 3 (A_3): A_1, A_2, A_1:A_2", out, fixed = TRUE)))

  #A single time point with an intercept-only numerator has no factors to name
  W1 <- weightitMSM(msm_formulas[1L], data = msmdata, method = "glm",
                    stabilize = TRUE)

  out1 <- capture.output(print(W1))
  expect_true(any(grepl(" - stabilized", out1, fixed = TRUE)))
  expect_false(any(grepl("stabilization factors", out1, fixed = TRUE)))

  #The stabilization block ends its line in both shapes, so a following line (here
  #the trimming note) is not run onto the end of it
  for (Wi in list(W, W1)) {
    out_i <- capture.output(print(trim(Wi, at = 0.9)))
    expect_identical(out_i[length(out_i)], " - weights trimmed at 90%")
  }
})
