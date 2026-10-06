# The four cases below assert that `vcov`, `R`, and `cluster` are honored identically by
# `vcov()`, `summary()`, `anova()`, and `update()` for each model class. They check
# reproducibility given a seed and well-formedness, never statistical accuracy, so `R`
# only has to exceed the number of parameters for the bootstrap covariance to be
# full-rank -- `R = 10` exercises exactly the same code paths as a realistic number of
# replicates. Every replicate re-estimates the weights as well as the outcome model, so
# this is the difference between ~560 and ~1400 model fits in this file alone.

test_that("vcov arg works in vcov(), summary(), and anova() for each model class", {
  skip_on_cran()
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  patrick::with_parameters_test_that(
    "vcov arg works in vcov(), summary(), and anova() for {fitter}",
    {
      if (fitter == "coxph_weightit") {
        skip_if_not_installed("survival")
      }

      test_data <- readRDS(test_path("fixtures", "test_data.rds"))
      set.seed(123)
      if (draw_off) {
        test_data$off <- runif(nrow(test_data))
      }
      test_data$clus <- sample(1:50, nrow(test_data), replace = TRUE)
      test_data <- prep(test_data)

      w_call <- quote(weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                               data = test_data, method = "glm", estimand = "ATE"))
      if (quick) {
        w_call$quick <- TRUE
      }

      W <- eval(w_call)

      # `update()` re-evaluates the stored call in the caller's frame, and the checks
      # below compare whole objects, `$call` included. So the call is built from the
      # fitter's name and unevaluated arguments, and evaluated here, to read exactly as
      # if it had been typed.
      fit_call <- function(rhs, ...) {
        rlang::call2(fitter, call("~", outcome, rhs),
                     data = quote(test_data), weightit = quote(W), ...)
      }

      fit_none <- eval(fit_call(quote(A * (X1)), vcov = "none"))
      fit_asympt <- eval(fit_call(quote(A * (X1)), vcov = "asympt"))
      fit_hc0 <- eval(fit_call(quote(A * (X1)), vcov = "HC0"))
      set.seed(123)
      fit_bs <- eval(fit_call(quote(A * (X1)), vcov = "BS", R = 10))

      set.seed(123)
      fit_fwb <- eval(fit_call(quote(A * (X1)), vcov = "FWB", R = 10))

      fit_asympt_clus <- eval(fit_call(quote(A * (X1)), vcov = "asympt",
                                       cluster = quote(~clus)))
      fit_hc0_clus <- eval(fit_call(quote(A * (X1)), vcov = "HC0",
                                    cluster = quote(~clus)))
      set.seed(123)
      fit_bs_clus <- eval(fit_call(quote(A * (X1)), vcov = "BS", R = 10,
                                   cluster = quote(~clus)))

      set.seed(123)
      fit_fwb_clus <- eval(fit_call(quote(A * (X1)), vcov = "FWB", R = 10,
                                    cluster = quote(~clus)))

      expect_equal(vcov(fit_none, vcov = "asympt"),
                   vcov(fit_asympt),
                   tolerance = eps)

      expect_equal(vcov(fit_none, vcov = "HC0"),
                   vcov(fit_hc0),
                   tolerance = eps)

      set.seed(123)
      expect_equal(vcov(fit_none, vcov = "BS", R = 10),
                   vcov(fit_bs),
                   tolerance = eps)

      set.seed(123)
      expect_equal(vcov(fit_none, vcov = "FWB", R = 10),
                   vcov(fit_fwb),
                   tolerance = eps)

      expect_equal(vcov(fit_none, vcov = "asympt", cluster = ~clus),
                   vcov(fit_asympt_clus),
                   tolerance = eps)

      expect_equal(vcov(fit_none, vcov = "HC0", cluster = ~clus),
                   vcov(fit_hc0_clus),
                   tolerance = eps)

      set.seed(123)
      expect_equal(vcov(fit_none, vcov = "BS", R = 10, cluster = ~clus),
                   vcov(fit_bs_clus),
                   tolerance = eps)

      set.seed(123)
      expect_equal(vcov(fit_none, vcov = "FWB", R = 10, cluster = ~clus),
                   vcov(fit_fwb_clus),
                   tolerance = eps)

      expect_equal(vcov(fit_asympt_clus, vcov = "asympt", cluster = NULL),
                   vcov(fit_asympt),
                   tolerance = eps)

      expect_equal(vcov(fit_asympt_clus, vcov = "asympt"),
                   vcov(fit_asympt_clus),
                   tolerance = eps)

      expect_equal(vcov(fit_asympt_clus, vcov = "HC0", cluster = NULL),
                   vcov(fit_hc0),
                   tolerance = eps)

      expect_equal(vcov(fit_asympt_clus, vcov = "HC0"),
                   vcov(fit_hc0_clus),
                   tolerance = eps)

      expect_equal(summary(fit_asympt, vcov = "HC0")$coef,
                   summary(fit_hc0)$coef,
                   tolerance = eps)

      expect_equal(summary(fit_asympt, vcov = "HC0", cluster = ~clus)$coef,
                   summary(fit_hc0_clus)$coef,
                   tolerance = eps)

      expect_equal(summary(fit_asympt_clus, vcov = "HC0", cluster = NULL)$coef,
                   summary(fit_hc0)$coef,
                   tolerance = eps)

      expect_equal(summary(fit_asympt_clus, vcov = "HC0")$coef,
                   summary(fit_hc0_clus)$coef,
                   tolerance = eps)

      set.seed(123)
      expect_equal(summary(fit_asympt_clus, vcov = "BS", R = 10)$coef,
                   summary(fit_bs_clus)$coef,
                   tolerance = eps)

      fit_small <- eval(fit_call(quote(A), vcov = "none"))

      expect_equal(anova(fit_asympt, fit_small),
                   anova(fit_none, fit_small, vcov = "asympt"),
                   tolerance = eps)

      expect_equal(anova(fit_hc0, fit_small),
                   anova(fit_none, fit_small, vcov = "HC0"),
                   tolerance = eps)

      expect_equal(anova(fit_asympt_clus, fit_small),
                   anova(fit_none, fit_small, vcov = "asympt", cluster = ~clus),
                   tolerance = eps)

      expect_equal(anova(fit_hc0_clus, fit_small),
                   anova(fit_none, fit_small, vcov = "HC0", cluster = ~clus),
                   tolerance = eps)

      set.seed(123)
      expect_equal(anova(fit_bs_clus, fit_small),
                   anova(fit_none, fit_small, vcov = "BS", R = 10, cluster = ~clus),
                   tolerance = eps)

      expect_error(anova(fit_asympt_clus, fit_small, vcov = "none"),
                   "no variance matrix was found", ignore.case = TRUE)

      fit_small_hc0 <- eval(fit_call(quote(A), vcov = "HC0"))

      expect_warning(anova(fit_asympt, fit_small_hc0),
                     "different `vcov` types detected", ignore.case = TRUE)

      expect_no_condition(anova(fit_hc0, fit_small_hc0))

      expect_no_condition(anova(fit_asympt, fit_small_hc0, vcov = "asympt"))

      # For ordinal and multinomial models these refit, since the fit computes the
      # Hessian only for the variance types that use it
      expect_equal(update(fit_none, vcov = "HC0"),
                   fit_hc0,
                   tolerance = eps)

      expect_equal(update(fit_none, vcov = "asympt"),
                   fit_asympt,
                   tolerance = eps)

      expect_equal(update(fit_hc0, vcov = "asympt"),
                   fit_asympt,
                   tolerance = eps)

      set.seed(123)
      expect_equal(fit_bs,
                   update(fit_none, vcov = "BS", R = 10),
                   tolerance = eps)

      set.seed(123)
      expect_equal(fit_fwb,
                   update(fit_none, vcov = "FWB", R = 10),
                   tolerance = eps)

      expect_equal(update(fit_hc0, cluster = ~clus),
                   fit_hc0_clus,
                   tolerance = eps)

      expect_identical(update(fit_asympt, cluster = ~clus),
                       fit_asympt_clus)

      expect_equal(update(fit_asympt_clus, cluster = NULL),
                   fit_asympt,
                   tolerance = eps)

      #Note: need to remove call because order of arguments is different
      .remove_call <- function(x) {
        x$call <- NULL
        x
      }

      set.seed(123)
      expect_equal(.remove_call(update(fit_bs, R = 10, cluster = ~clus)),
                   .remove_call(fit_bs_clus),
                   tolerance = eps)

      set.seed(123)
      expect_equal(.remove_call(update(fit_fwb, cluster = ~clus, R = 10)),
                   .remove_call(fit_fwb_clus),
                   tolerance = eps)
    },
    # `quick`: whether `weightit()` is passed `quick = TRUE`.
    # `draw_off`: whether an unused `off` column is drawn before `clus`, which changes
    #   `clus`; kept so each model sees the same data it always has.
    # `prep`: builds the outcome that `outcome` names.
    .cases = patrick::cases(
      list(fitter = "glm_weightit", outcome = quote(Y_C),
           quick = TRUE, draw_off = FALSE,
           prep = identity),
      list(fitter = "ordinal_weightit", outcome = quote(Y_O),
           quick = FALSE, draw_off = FALSE,
           prep = identity),
      list(fitter = "multinom_weightit", outcome = quote(Y_M),
           quick = FALSE, draw_off = TRUE,
           prep = function(d) {
             d$Y_M <- factor(d$Y_O, ordered = FALSE)
             d
           }),
      list(fitter = "coxph_weightit", outcome = quote(survival::Surv(time, event)),
           quick = FALSE, draw_off = FALSE,
           prep = function(d) {
             # Y_S is uncensored by construction; censor at the 80th percentile so
             # there's a mix of events and censoring to fit against
             d$event <- as.numeric(d$Y_S < quantile(d$Y_S, .8))
             d$time <- pmin(d$Y_S, quantile(d$Y_S, .8))
             d
           })
    )
  )
})
