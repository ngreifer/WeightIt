test_that("vcov arg works in vcov(), summary(), and anova() for glm_weightit", {
  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  set.seed(123)
  test_data$clus <- sample(1:50, nrow(test_data), replace = TRUE)

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE", quick = TRUE)

  fit_none <- glm_weightit(Y_C ~ A * (X1),
                           data = test_data, weightit = W, vcov = "none")
  fit_asympt <- glm_weightit(Y_C ~ A * (X1),
                             data = test_data, weightit = W, vcov = "asympt")
  fit_hc0 <- glm_weightit(Y_C ~ A * (X1),
                          data = test_data, weightit = W, vcov = "HC0")
  set.seed(123)
  fit_bs <- glm_weightit(Y_C ~ A * (X1),
                         data = test_data, weightit = W, vcov = "BS", R = 25)

  set.seed(123)
  fit_fwb <- glm_weightit(Y_C ~ A * (X1),
                          data = test_data, weightit = W, vcov = "FWB", R = 25)

  fit_asympt_clus <- glm_weightit(Y_C ~ A * (X1),
                                  data = test_data, weightit = W, vcov = "asympt",
                                  cluster = ~clus)
  fit_hc0_clus <- glm_weightit(Y_C ~ A * (X1),
                               data = test_data, weightit = W, vcov = "HC0",
                               cluster = ~clus)
  set.seed(123)
  fit_bs_clus <- glm_weightit(Y_C ~ A * (X1),
                              data = test_data, weightit = W, vcov = "BS", R = 25,
                              cluster = ~clus)

  set.seed(123)
  fit_fwb_clus <- glm_weightit(Y_C ~ A * (X1),
                               data = test_data, weightit = W, vcov = "FWB", R = 25,
                               cluster = ~clus)


  expect_equal(vcov(fit_none, vcov = "asympt"),
               vcov(fit_asympt),
               tolerance = eps)

  expect_equal(vcov(fit_none, vcov = "HC0"),
               vcov(fit_hc0),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "BS", R = 25),
               vcov(fit_bs),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "FWB", R = 25),
               vcov(fit_fwb),
               tolerance = eps)

  expect_equal(vcov(fit_none, vcov = "asympt", cluster = ~clus),
               vcov(fit_asympt_clus),
               tolerance = eps)

  expect_equal(vcov(fit_none, vcov = "HC0", cluster = ~clus),
               vcov(fit_hc0_clus),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "BS", R = 25, cluster = ~clus),
               vcov(fit_bs_clus),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "FWB", R = 25, cluster = ~clus),
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
  expect_equal(summary(fit_asympt_clus, vcov = "BS", R = 25)$coef,
               summary(fit_bs_clus)$coef,
               tolerance = eps)

  fit_small <- glm_weightit(Y_C ~ A,
                           data = test_data, weightit = W, vcov = "none")

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
               anova(fit_none, fit_small, vcov = "BS", R = 25, cluster = ~clus),
               tolerance = eps)

  expect_error(anova(fit_asympt_clus, fit_small, vcov = "none"),
               "No variance matrix was found")

  fit_small_hc0 <- glm_weightit(Y_C ~ A,
                            data = test_data, weightit = W, vcov = "HC0")

  expect_warning(anova(fit_asympt, fit_small_hc0),
                 "Different `vcov` types detected")

  expect_no_condition(anova(fit_hc0, fit_small_hc0))

  expect_no_condition(anova(fit_asympt, fit_small_hc0, vcov = "asympt"))

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
               update(fit_none, vcov = "BS", R = 25),
               tolerance = eps)

  set.seed(123)
  expect_equal(fit_fwb,
               update(fit_none, vcov = "FWB", R = 25),
               tolerance = eps)

  expect_equal(update(fit_hc0, cluster = ~clus),
               fit_hc0_clus,
               tolerance = eps)

  expect_equal(update(fit_asympt, cluster = ~clus),
               fit_asympt_clus,
               tolerance = eps)

  expect_equal(update(fit_asympt_clus, cluster = NULL),
               fit_asympt,
               tolerance = eps)

  #Note: need to remove call because order of arguments is different
  .remove_call <- function(x) {
    x$call <- NULL
    x
  }

  set.seed(123)
  expect_equal(.remove_call(update(fit_bs, R = 25, cluster = ~clus)),
               .remove_call(fit_bs_clus),
               tolerance = eps)

  set.seed(123)
  expect_equal(.remove_call(update(fit_fwb, cluster = ~clus, R = 25)),
               .remove_call(fit_fwb_clus),
               tolerance = eps)
})

test_that("vcov arg works in vcov(), summary(), and anova() for ordinal_weightit", {
  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  set.seed(123)
  test_data$clus <- sample(1:50, nrow(test_data), replace = TRUE)
  test_data$Y_O <- with(test_data, factor(findInterval(Y_C, quantile(Y_C, seq(0, 1, length = 5)),
                                                       all.inside = TRUE), ordered = TRUE))

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  fit_none <- ordinal_weightit(Y_O ~ A * (X1),
                           data = test_data, weightit = W, vcov = "none")
  fit_asympt <- ordinal_weightit(Y_O ~ A * (X1),
                             data = test_data, weightit = W, vcov = "asympt")
  fit_hc0 <- ordinal_weightit(Y_O ~ A * (X1),
                          data = test_data, weightit = W, vcov = "HC0")
  set.seed(123)
  fit_bs <- ordinal_weightit(Y_O ~ A * (X1),
                         data = test_data, weightit = W, vcov = "BS", R = 25)

  set.seed(123)
  fit_fwb <- ordinal_weightit(Y_O ~ A * (X1),
                          data = test_data, weightit = W, vcov = "FWB", R = 25)

  fit_asympt_clus <- ordinal_weightit(Y_O ~ A * (X1),
                                  data = test_data, weightit = W, vcov = "asympt",
                                  cluster = ~clus)
  fit_hc0_clus <- ordinal_weightit(Y_O ~ A * (X1),
                               data = test_data, weightit = W, vcov = "HC0",
                               cluster = ~clus)
  set.seed(123)
  fit_bs_clus <- ordinal_weightit(Y_O ~ A * (X1),
                              data = test_data, weightit = W, vcov = "BS", R = 25,
                              cluster = ~clus)

  set.seed(123)
  fit_fwb_clus <- ordinal_weightit(Y_O ~ A * (X1),
                               data = test_data, weightit = W, vcov = "FWB", R = 25,
                               cluster = ~clus)


  expect_equal(vcov(fit_none, vcov = "asympt"),
               vcov(fit_asympt),
               tolerance = eps)

  expect_equal(vcov(fit_none, vcov = "HC0"),
               vcov(fit_hc0),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "BS", R = 25),
               vcov(fit_bs),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "FWB", R = 25),
               vcov(fit_fwb),
               tolerance = eps)

  expect_equal(vcov(fit_none, vcov = "asympt", cluster = ~clus),
               vcov(fit_asympt_clus),
               tolerance = eps)

  expect_equal(vcov(fit_none, vcov = "HC0", cluster = ~clus),
               vcov(fit_hc0_clus),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "BS", R = 25, cluster = ~clus),
               vcov(fit_bs_clus),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "FWB", R = 25, cluster = ~clus),
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
  expect_equal(summary(fit_asympt_clus, vcov = "BS", R = 25)$coef,
               summary(fit_bs_clus)$coef,
               tolerance = eps)

  fit_small <- ordinal_weightit(Y_O ~ A,
                            data = test_data, weightit = W, vcov = "none")

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
               anova(fit_none, fit_small, vcov = "BS", R = 25, cluster = ~clus),
               tolerance = eps)

  expect_error(anova(fit_asympt_clus, fit_small, vcov = "none"),
               "No variance matrix was found")

  fit_small_hc0 <- ordinal_weightit(Y_O ~ A,
                                data = test_data, weightit = W, vcov = "HC0")

  expect_warning(anova(fit_asympt, fit_small_hc0),
                 "Different `vcov` types detected")

  expect_no_condition(anova(fit_hc0, fit_small_hc0))

  expect_no_condition(anova(fit_asympt, fit_small_hc0, vcov = "asympt"))

  expect_equal(summary(update(fit_none, vcov = "HC0")),
               summary(fit_hc0),
               tolerance = eps)

  expect_equal(summary(update(fit_none, vcov = "asympt")),
               summary(fit_asympt),
               tolerance = eps)

  expect_equal(update(fit_hc0, vcov = "asympt"),
               fit_asympt,
               tolerance = eps)

  set.seed(123)
  expect_equal(fit_bs,
               update(fit_none, vcov = "BS", R = 25),
               tolerance = eps)

  set.seed(123)
  expect_equal(fit_fwb,
               update(fit_none, vcov = "FWB", R = 25),
               tolerance = eps)

  expect_equal(update(fit_hc0, cluster = ~clus),
               fit_hc0_clus,
               tolerance = eps)

  expect_equal(update(fit_asympt, cluster = ~clus),
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
  expect_equal(.remove_call(update(fit_bs, R = 25, cluster = ~clus)),
               .remove_call(fit_bs_clus),
               tolerance = eps)

  set.seed(123)
  expect_equal(.remove_call(update(fit_fwb, cluster = ~clus, R = 25)),
               .remove_call(fit_fwb_clus),
               tolerance = eps)
})

test_that("vcov arg works in vcov(), summary(), and anova() for multinom_weightit", {
  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))
  test_data$Y_M <- with(test_data, factor(findInterval(Y_C, quantile(Y_C, seq(0, 1, length = 5)),
                                                       all.inside = TRUE)))
  set.seed(123)
  test_data$off <- runif(nrow(test_data))
  test_data$clus <- sample(1:50, nrow(test_data), replace = TRUE)

  W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                data = test_data, method = "glm", estimand = "ATE")

  fit_none <- multinom_weightit(Y_M ~ A * (X1),
                               data = test_data, weightit = W, vcov = "none")
  fit_asympt <- multinom_weightit(Y_M ~ A * (X1),
                                 data = test_data, weightit = W, vcov = "asympt")
  fit_hc0 <- multinom_weightit(Y_M ~ A * (X1),
                              data = test_data, weightit = W, vcov = "HC0")
  set.seed(123)
  fit_bs <- multinom_weightit(Y_M ~ A * (X1),
                             data = test_data, weightit = W, vcov = "BS", R = 25)

  set.seed(123)
  fit_fwb <- multinom_weightit(Y_M ~ A * (X1),
                              data = test_data, weightit = W, vcov = "FWB", R = 25)

  fit_asympt_clus <- multinom_weightit(Y_M ~ A * (X1),
                                      data = test_data, weightit = W, vcov = "asympt",
                                      cluster = ~clus)
  fit_hc0_clus <- multinom_weightit(Y_M ~ A * (X1),
                                   data = test_data, weightit = W, vcov = "HC0",
                                   cluster = ~clus)
  set.seed(123)
  fit_bs_clus <- multinom_weightit(Y_M ~ A * (X1),
                                  data = test_data, weightit = W, vcov = "BS", R = 25,
                                  cluster = ~clus)

  set.seed(123)
  fit_fwb_clus <- multinom_weightit(Y_M ~ A * (X1),
                                   data = test_data, weightit = W, vcov = "FWB", R = 25,
                                   cluster = ~clus)

  expect_equal(vcov(fit_none, vcov = "asympt"),
               vcov(fit_asympt),
               tolerance = eps)

  expect_equal(vcov(fit_none, vcov = "HC0"),
               vcov(fit_hc0),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "BS", R = 25),
               vcov(fit_bs),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "FWB", R = 25),
               vcov(fit_fwb),
               tolerance = eps)

  expect_equal(vcov(fit_none, vcov = "asympt", cluster = ~clus),
               vcov(fit_asympt_clus),
               tolerance = eps)

  expect_equal(vcov(fit_none, vcov = "HC0", cluster = ~clus),
               vcov(fit_hc0_clus),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "BS", R = 25, cluster = ~clus),
               vcov(fit_bs_clus),
               tolerance = eps)

  set.seed(123)
  expect_equal(vcov(fit_none, vcov = "FWB", R = 25, cluster = ~clus),
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
  expect_equal(summary(fit_asympt_clus, vcov = "BS", R = 25)$coef,
               summary(fit_bs_clus)$coef,
               tolerance = eps)

  fit_small <- multinom_weightit(Y_M ~ A,
                                data = test_data, weightit = W, vcov = "none")

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
               anova(fit_none, fit_small, vcov = "BS", R = 25, cluster = ~clus),
               tolerance = eps)

  expect_error(anova(fit_asympt_clus, fit_small, vcov = "none"),
               "No variance matrix was found")

  fit_small_hc0 <- multinom_weightit(Y_M ~ A,
                                    data = test_data, weightit = W, vcov = "HC0")

  expect_warning(anova(fit_asympt, fit_small_hc0),
                 "Different `vcov` types detected")

  expect_no_condition(anova(fit_hc0, fit_small_hc0))

  expect_no_condition(anova(fit_asympt, fit_small_hc0, vcov = "asympt"))

  expect_equal(summary(update(fit_none, vcov = "HC0")),
               summary(fit_hc0),
               tolerance = eps)

  expect_equal(summary(update(fit_none, vcov = "asympt")),
               summary(fit_asympt),
               tolerance = eps)

  expect_equal(update(fit_hc0, vcov = "asympt"),
               fit_asympt,
               tolerance = eps)

  set.seed(123)
  expect_equal(fit_bs,
               update(fit_none, vcov = "BS", R = 25),
               tolerance = eps)

  set.seed(123)
  expect_equal(fit_fwb,
               update(fit_none, vcov = "FWB", R = 25),
               tolerance = eps)

  expect_equal(update(fit_hc0, cluster = ~clus),
               fit_hc0_clus,
               tolerance = eps)

  expect_equal(update(fit_asympt, cluster = ~clus),
               fit_asympt_clus,
               tolerance = eps)

  expect_equal(update(fit_asympt_clus, cluster = NULL),
               fit_asympt,
               tolerance = eps)

  #Note: need to remove call because order of arguments is different
  .remove_call <- function(x) {
    x$call <- NULL
    x
  }

  set.seed(123)
  expect_equal(.remove_call(update(fit_bs, R = 25, cluster = ~clus)),
               .remove_call(fit_bs_clus),
               tolerance = eps)

  set.seed(123)
  expect_equal(.remove_call(update(fit_fwb, cluster = ~clus, R = 25)),
               .remove_call(fit_fwb_clus),
               tolerance = eps)
})