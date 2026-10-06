test_that("Binary treatment", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "ipt", estimand = "ATE",
                   include.obj = TRUE)
  })

  expect_M_parts_okay(W0, tolerance = eps)

  # Weights from the configurations already run, so each new one can be checked
  # against all of them
  seen <- new.env()

  patrick::with_parameters_test_that(
    "IPT: sw = {sw}, estimand = {estimand}, link = {link}",
    {
      W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "ipt", estimand = estimand,
                    link = link,
                    s.weights = if (sw) "SW" else NULL,
                    include.obj = TRUE)

      expect_M_parts_okay(W, tolerance = eps)
      expect_equal(cobalt::col_w_smd(W$covs, W$treat, W$weights,
                                     s.weights = W$s.weights),
                   0 * cobalt::col_w_smd(W$covs, W$treat,
                                         s.weights = W$s.weights),
                   expected.label = "all 0s",
                   tolerance = eps)

      expect_true(is.numeric(W$ps))
      expect_false(is_null(W$obj))

      if (estimand %in% c("ATT", "ATC")) {
        expect_ATT_weights_okay(W, tolerance = eps)
      }

      for (i in 0:1) {
        e <- {
          if (estimand == "ATT" && i == 1) expect_equal
          else if (estimand == "ATC" && i == 0) expect_equal
          else expect_not_equal
        }

        e(unname(W$weights[W$treat == i]),
          rep(1, sum(W$treat == i)),
          label = sprintf("%s weights", i),
          expected.label = "all 1s",
          tolerance = eps)
      }

      for (other in ls(seen)) {
        expect_not_equal(unname(W$weights), seen[[other]],
                         expected.label = sprintf("weights for %s", other),
                         tolerance = eps)
      }

      seen[[sprintf("sw = %s, estimand = %s, link = %s", sw, estimand, link)]] <- unname(W$weights)
    },
    .cases = expand.grid(link = c("logit", "probit", "loglog", "cauchit"),
                         estimand = c("ATE", "ATT", "ATC"),
                         sw = c(FALSE, TRUE),
                         stringsAsFactors = FALSE)
  )

  # Estimands
  expect_error({
    W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                  data = test_data, method = "ipt", estimand = "ATO")
  }, "not an allowable estimand", ignore.case = TRUE)

  #Non-full rank
  W <- expect_no_condition({
    weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 +
               I(1 - X5) + I(X9 * 2),
             data = test_data, method = "ipt", estimand = "ATE",
             include.obj = TRUE)
  })

  expect_M_parts_okay(W, tolerance = eps)
  expect_equal(W$weights, W0$weights, tolerance = eps)

  #Should be equivalent to CBPS for ATT
  patrick::with_parameters_test_that(
    "IPT matches CBPS for ATT: sw = {sw}, link = {link}",
    {
      W <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                    data = test_data, method = "ipt", estimand = "ATT",
                    s.weights = if (sw) "SW" else NULL,
                    link = link,
                    include.obj = TRUE)

      Wcbps <- weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                        data = test_data, method = "cbps", estimand = "ATT",
                        s.weights = if (sw) "SW" else NULL,
                        link = link, solver = "multiroot",
                        include.obj = TRUE)

      expect_equal(ESS(W$weights[W$treat == 0] * W$s.weights[W$treat == 0]),
                   ESS(Wcbps$weights[Wcbps$treat == 0] * Wcbps$s.weights[Wcbps$treat == 0]),
                   expected.label = "ESS for CBPS",
                   tolerance = .01)
    },
    .cases = expand.grid(link = c("logit", "probit", "loglog", "cauchit"),
                         sw = c(FALSE, TRUE),
                         stringsAsFactors = FALSE)
  )
})

test_that("Multi-category treatment", {
  skip_on_cran()
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")
  skip_if_not_installed("patrick")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-3

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  expect_no_condition({
    W0 <- weightit(Am ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9,
                   data = test_data, method = "ipt", estimand = "ATE",
                   include.obj = TRUE)
  })

  # Weights from the configurations already run, so each new one can be checked
  # against all of them
  seen <- new.env()

  patrick::with_parameters_test_that(
    "IPT: sw = {sw}, estimand = {estimand}, link = {link}",
    {
      W <- weightit(Am ~ X1 + X2 + X3 + X4 + X5,
                    data = test_data, method = "ipt", estimand = estimand,
                    focal = if (estimand == "ATE") NULL else "T",
                    s.weights = if (sw) "SW" else NULL,
                    link = link,
                    include.obj = TRUE)

      expect_M_parts_okay(W, tolerance = eps)
      for (tt in combn(levels(W$treat), 2, simplify = FALSE)) {
        in_tt <- W$treat %in% tt
        expect_equal(cobalt::col_w_smd(W$covs[in_tt,], W$treat[in_tt], W$weights[in_tt],
                                       s.weights = W$s.weights[in_tt]),
                     0 * cobalt::col_w_smd(W$covs[in_tt,], W$treat[in_tt],
                                           s.weights = W$s.weights[in_tt]),
                     label = sprintf("SMDs for %s", paste(tt, collapse = " vs. ")),
                     expected.label = "all 0s",
                     tolerance = eps)
      }

      expect_true(is_null(W$ps))
      expect_false(is_null(W$obj))

      if (estimand %in% c("ATT", "ATC")) {
        expect_ATT_weights_okay(W, tolerance = eps)
      }

      for (i in levels(W$treat)) {
        e <- {
          if (estimand == "ATT" && i == W$focal) expect_equal
          else expect_not_equal
        }

        e(unname(W$weights[W$treat == i]),
          rep(1, sum(W$treat == i)),
          label = sprintf("%s weights", i),
          expected.label = "all 1s",
          tolerance = eps)
      }

      for (other in ls(seen)) {
        expect_not_equal(unname(W$weights), seen[[other]],
                         expected.label = sprintf("weights for %s", other),
                         tolerance = eps)
      }

      seen[[sprintf("sw = %s, estimand = %s, link = %s", sw, estimand, link)]] <- unname(W$weights)
    },
    .cases = expand.grid(link = c("logit", "probit", "loglog", "cauchit"),
                         estimand = c("ATE", "ATT"),
                         sw = c(FALSE, TRUE),
                         stringsAsFactors = FALSE)
  )
})
