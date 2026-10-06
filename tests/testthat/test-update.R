test_that("update.glm_weightit() works", {
  skip_on_cran()
  # skip_if(!capabilities("long.double"))
  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    fit0 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial)
  })

  #Updating formula
  fit1 <- expect_no_condition({
    glm_weightit(Y_B ~ A,
                 data = test_data, family = binomial)
  })

  expect_equal(update(fit0, formula = . ~ A),
               fit1,
               tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating dataset
  sub <- test_data$X3 > 0
  test_data_s <- test_data[sub,]

  fit1 <- expect_no_condition({
    glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                 data = test_data_s, family = binomial)
  })

  expect_equal(update(fit0, data = test_data_s),
               fit1,
               tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  expect_equal(update(fit0, subset = sub)$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating vcov only
  fit1 <- expect_no_condition({
    glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                 data = test_data, family = binomial, vcov = "const")
  })

  expect_equal(update(fit0, vcov = "const"),
               fit1,
               tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  expect_not_equal(vcov(fit0), vcov(fit1))

  #Model should not be refit when only vcov is changed
  set.seed(123)
  clus <- sample(1:50, nrow(test_data), replace = TRUE)
  i <- FALSE
  expect_no_condition({
    fit2 <- glm_weightit(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                         data = test_data, family = binomial,
                         control = if (i) list(stop("bad error")) else list())
  })

  i <- TRUE
  expect_no_condition({
    update(fit2, vcov = "const")
  })

  expect_error({
    update(fit2, family = binomial("probit"))
  }, "bad error", ignore.case = TRUE)

  expect_error({
    update(fit2)
  }, "bad error", ignore.case = TRUE)

  expect_no_condition({
    fit2c <- update(fit2, cluster = ~clus)
  })

  expect_equal(vcov(fit2c),
               vcov(update(fit0, cluster = ~clus)),
               tolerance = eps)

  expect_not_equal(vcov(fit2c),
                   vcov(fit0))

  #Updating s.weights without weightit object should add one
  expect_no_condition({
    suppressWarnings({
      fitg <- glm(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                  data = test_data, family = binomial,
                  weights = SW)
    })
  })

  expect_no_condition({
    fit0s <- update(fit0, s.weights = "SW")
  })

  expect_null(fit0$weightit)
  expect_equal(fit0s$weightit,
               list(s.weights = test_data$SW,
                    weights = rep(1, nrow(test_data)),
                    method = NULL),
               ignore_attr = c("names", "class"))
  expect_equal(fit0s$weightit$s.weights,
               fit0s$prior.weights,
               ignore_attr = "names")

  expect_equal(coef(fit0s),
               coef(fitg),
               tolerance = eps)

})

test_that("update.weightit() works", {
  skip_on_cran()
  # skip_if(!capabilities("long.double"))
  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data)
  })

  sub <- test_data$X3 > 0
  test_data_s <- test_data[sub,]

  skip_if_not_installed("patrick")

  # Each case refits with the arguments in `refit` and checks that `update()` with the
  # arguments in `changes` gives the same object, `$call` included. Both calls go through
  # `do.call()` with the function's name and unevaluated arguments, so each reads exactly
  # as if it had been typed.
  patrick::with_parameters_test_that(
    "updating {change}",
    {
      W1 <- expect_no_condition({
        do.call("weightit", refit)
      })

      expect_equal(do.call("update", c(list(quote(W0)), changes)),
                   W1, tolerance = eps,
                   ignore_attr = "class",
                   ignore_formula_env = TRUE)
    },
    .cases = patrick::cases(
      #Updating formula
      list(change = "formula",
           refit = list(quote(A ~ X1 + X2 + X3 + X4),
                        data = quote(test_data)),
           changes = list(formula = quote(. ~ X1 + X2 + X3 + X4))),
      #Updating dataset
      list(change = "dataset",
           refit = list(quote(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = quote(test_data_s)),
           changes = list(data = quote(test_data_s))),
      #Updating method and estimand
      list(change = "method and estimand",
           refit = list(quote(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = quote(test_data), method = "ebal", estimand = "ATT"),
           changes = list(method = "ebal", estimand = "ATT")),
      #Updating s.weights
      list(change = "s.weights",
           refit = list(quote(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                        data = quote(test_data), s.weights = "SW"),
           changes = list(s.weights = "SW"))
    )
  )
})

test_that("update() works with weightit for each model class", {
  skip_on_cran()
  skip_if_not_installed("patrick")
  # skip_if(!capabilities("long.double"))
  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  patrick::with_parameters_test_that(
    "update.{fitter}() works with weightit",
    {
      test_data <- readRDS(test_path("fixtures", "test_data.rds"))
      test_data <- prep(test_data)

      # `update()` re-evaluates the stored call in the caller's frame, and the checks
      # below compare whole objects, `$call` included. So each call is built from the
      # fitter's name and unevaluated arguments, and evaluated here, to read exactly as
      # if it had been typed. The "should not be refit" checks depend on this: `control`
      # stays unevaluated in the stored call, so a refit re-evaluates it where `i` is
      # `TRUE`.
      fit_call <- function(data, ...) {
        cl <- rlang::call2(fitter, fml, data = data)

        if (binomial_family) {
          cl$family <- quote(binomial)
        }

        as.call(c(as.list(cl), list(...)))
      }

      expect_no_condition({
        W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                       data = test_data)
      })

      expect_no_condition({
        fit0 <- eval(fit_call(quote(test_data), weightit = quote(W0)))
      })

      #Updating weightit
      fit1 <- expect_no_condition({
        eval(fit_call(quote(test_data)))
      })

      expect_equal(update(fit0, weightit = NULL),
                   fit1, tolerance = eps,
                   ignore_attr = "class",
                   ignore_formula_env = TRUE)

      #Updating dataset
      sub <- test_data$X3 > 0
      test_data_s <- test_data[sub,]

      W1 <- expect_no_condition({
        weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                 data = test_data_s)
      })

      fit1 <- expect_no_condition({
        eval(fit_call(quote(test_data_s), weightit = quote(W1)))
      })

      expect_equal(update(fit0, data = test_data_s)$coefficients,
                   fit1$coefficients, tolerance = eps,
                   ignore_attr = "class",
                   ignore_formula_env = TRUE)

      #Updating vcov only
      fit1 <- expect_no_condition({
        eval(fit_call(quote(test_data), weightit = quote(W0),
                      vcov = "HC0"))
      })

      expect_equal(update(fit0, vcov = "HC0"),
                   fit1, tolerance = eps,
                   ignore_attr = "class",
                   ignore_formula_env = TRUE)

      expect_not_equal(vcov(fit0), vcov(fit1))

      #Model should not be refit when only vcov is changed
      set.seed(123)
      clus <- sample(1:50, nrow(test_data), replace = TRUE)
      i <- FALSE
      expect_no_condition({
        fit2 <- eval(fit_call(quote(test_data), weightit = quote(W0),
                              control = quote(if (i) list(stop("bad error")) else list())))
      })

      i <- TRUE
      expect_no_condition({
        update(fit2, vcov = "HC0")
      })

      if (binomial_family) {
        expect_error({
          update(fit2, family = binomial("probit"))
        }, "bad error", ignore.case = TRUE)
      }

      expect_error({
        update(fit2)
      }, "bad error", ignore.case = TRUE)

      expect_no_condition({
        fit2c <- update(fit2, cluster = ~clus)
      })

      expect_equal(vcov(fit2c),
                   vcov(update(fit0, cluster = ~clus)),
                   tolerance = eps)

      expect_not_equal(vcov(fit2c),
                       vcov(fit0))

      #Updating s.weights updates weightit
      W1 <- expect_no_condition({
        weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                 data = test_data, s.weights = "SW")
      })

      fit1 <- expect_no_condition({
        eval(fit_call(quote(test_data), weightit = quote(W1)))
      })

      expect_no_condition({
        fit1s <- eval(fit_call(quote(test_data),
                               weightit = quote(update(object = W0, s.weights = "SW"))))
      })

      expect_no_condition({
        fit0s <- update(fit0, s.weights = "SW")
      })

      expect_equal(fit0s$coefficients,
                   fit1$coefficients,
                   ignore_attr = c("names", "class"),
                   tolerance = eps)

      expect_equal(fit0s$coefficients,
                   fit1s$coefficients,
                   ignore_attr = c("names", "class"),
                   tolerance = eps)

      expect_equal(fit0s,
                   fit1s,
                   tolerance = eps)

      #Mimic bootstrapping
      boot_dat <- test_data[sample(nrow(test_data), replace = TRUE),]

      W1 <- expect_no_condition({
        weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                 data = boot_dat)
      })

      boot_call <- fit_call(quote(boot_dat), weightit = quote(W1))

      if (boot_vcov_none) {
        boot_call$vcov <- "none"
      }

      fit1 <- expect_no_condition({
        eval(boot_call)
      })

      expect_equal(update(fit0, data = boot_dat, vcov = "none")$coefficients,
                   fit1$coefficients, tolerance = eps,
                   ignore_attr = "class",
                   ignore_formula_env = TRUE)
    },
    # `prep`: builds the outcome that `fml` names.
    # `binomial_family`: whether each fit passes `family = binomial`, which also adds a
    #   refit check through a change of family.
    # `boot_vcov_none`: whether the bootstrap refit passes `vcov = "none"`.
    .cases = patrick::cases(
      list(fitter = "glm_weightit",
           fml = quote(Y_B ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9)),
           prep = identity,
           binomial_family = TRUE, boot_vcov_none = FALSE),
      list(fitter = "multinom_weightit",
           fml = quote(Y_M ~ A * (X1 + X2 + X3 + X4 + X5)),
           prep = function(d) {
             d$Y_M <- factor(d$Y_O, ordered = FALSE)
             d
           },
           binomial_family = FALSE, boot_vcov_none = TRUE),
      list(fitter = "ordinal_weightit",
           fml = quote(Y_O ~ A * (X1 + X2 + X3 + X4 + X5)),
           prep = identity,
           binomial_family = FALSE, boot_vcov_none = TRUE)
    )
  )
})

test_that("update.coxph_weightit() works with weightit", {
  skip_on_cran()
  # skip_if(!capabilities("long.double"))
  skip_if_not_installed("survival")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data)
  })

  expect_no_condition({
    fit0 <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                           data = test_data, weightit = W0)
  })

  #Updating weightit
  fit1 <- expect_no_condition({
    coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                   data = test_data)
  })

  expect_equal(update(fit0, weightit = NULL),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating dataset
  sub <- test_data$X3 > 0
  test_data_s <- test_data[sub,]

  W1 <- expect_no_condition({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data_s)
  })

  fit1 <- expect_no_condition({
    coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                   data = test_data_s, weightit = W1)
  })

  expect_equal(update(fit0, data = test_data_s)$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  #Updating vcov only
  set.seed(123)
  clus <- sample(1:50, nrow(test_data), replace = TRUE)

  fit1 <- expect_no_condition({
    coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                   data = test_data, weightit = W0,
                   cluster = ~clus)
  })

  expect_equal(update(fit0, cluster = ~clus),
               fit1, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)

  expect_not_equal(vcov(fit0), vcov(fit1))

  #Model should not be refit when only vcov is changed
  i <- FALSE
  expect_no_condition({
    fit2 <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                           data = test_data, weightit = W0,
                           control = if (i) list(stop("bad error")) else list())
  })

  if (exists("fit2", inherits = FALSE)) {
    i <- TRUE
    expect_error({
      update(fit2, x = TRUE)
    }, "bad error", ignore.case = TRUE)

    expect_error({
      update(fit2)
    }, "bad error", ignore.case = TRUE)

    expect_no_condition({
      fit2c <- update(fit2, cluster = ~clus)
    })

    expect_equal(vcov(fit2c),
                 vcov(update(fit0, cluster = ~clus)),
                 tolerance = eps)

    expect_not_equal(vcov(fit2c),
                     vcov(fit0))
  }

  #Updating s.weights updates weightit
  W1 <- expect_no_condition({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = test_data, s.weights = "SW")
  })

  fit1 <- expect_no_condition({
    coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                   data = test_data, weightit = W1)
  })

  expect_no_condition({
    fit1s <- coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                            data = test_data,
                            weightit = update(object = W0, s.weights = "SW"))
  })

  expect_no_condition({
    fit0s <- update(fit0, s.weights = "SW")
  })

  expect_equal(fit0s$coefficients,
               fit1$coefficients,
               ignore_attr = c("names", "class"),
               tolerance = eps)

  expect_equal(fit0s$coefficients,
               fit1s$coefficients,
               ignore_attr = c("names", "class"),
               tolerance = eps)

  expect_equal(fit0s,
               fit1s,
               tolerance = eps)

  #Mimic bootstrapping
  boot_dat <- test_data[sample(nrow(test_data), replace = TRUE),]

  W1 <- expect_no_condition({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
             data = boot_dat)
  })

  fit1 <- expect_no_condition({
    coxph_weightit(survival::Surv(Y_S) ~ A * (X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9),
                   data = boot_dat, weightit = W1)
  })

  expect_equal(update(fit0, data = boot_dat, vcov = "none")$coefficients,
               fit1$coefficients, tolerance = eps,
               ignore_attr = "class",
               ignore_formula_env = TRUE)
})
