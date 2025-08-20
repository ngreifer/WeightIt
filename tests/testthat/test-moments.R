test_that("moments works as expected, binary", {
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("optweight")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  poly_test <- cbind(test_data[c("X3", "X4")],
                     poly(test_data$X1, 2),
                     poly(test_data$X2, 3))

  init1 <- cobalt::bal.init(poly_test, treat = test_data$A, stat = "smd.max",
                            pairwise = FALSE)
  init2 <- cobalt::bal.init(poly(test_data$X1, 3), treat = test_data$A, stat = "smd.max",
                            pairwise = FALSE)

  #Moments the same as poly
  for (method in c("ipt", "ebal", "optweight")) {
    expect_no_condition({
      Wm <- weightit(A ~ X1 + X2 + X3 + X4,
                     data = test_data, method = method,
                     include.obj = TRUE,
                     moments = c(X1 = 2, X2 = 3))
    })
    nm <- sprintf("Wm_%s", method)
    assign(nm, Wm)

    expect_no_condition({
      Wp <- weightit(A ~ poly(X1, 2) + poly(X2, 3) + X3 + X4,
                     data = test_data, method = method,
                     include.obj = TRUE)
    })
    np <- sprintf("Wp_%s", method)
    assign(np, Wp)

    if (method %in% c("ipt", "ebal")) {
      expect_M_parts_okay(Wm, tolerance = eps)
      expect_M_parts_okay(Wp, tolerance = eps)
    }

    expect_equal(Wm$weights, Wp$weights,
                 label = sprintf("weights for %s", nm),
                 expected.label = sprintf("weights for %s", np),
                 tolerance = eps)


    expect_equal(cobalt::bal.compute(init1, weights = Wm$weights),
                 0,
                 label = sprintf("largest SMD for %s", nm),
                 tolerance = eps)
    expect_equal(cobalt::bal.compute(init1, weights = Wp$weights),
                 0,
                 label = sprintf("largest SMD for %s", np),
                 tolerance = eps)

    expect_not_equal(cobalt::bal.compute(init2, weights = Wm$weights),
                     0,
                     label = sprintf("largest SMD for %s with extra terms", nm),
                     tolerance = eps)
  }
})

test_that("moments works as expected, multi", {
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("optweight")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  poly_test <- cbind(test_data[c("X3", "X4")],
                     poly(test_data$X1, 2),
                     poly(test_data$X2, 3))

  init1 <- cobalt::bal.init(poly_test, treat = test_data$Am, stat = "smd.max",
                            pairwise = FALSE)
  init2 <- cobalt::bal.init(poly(test_data$X1, 3), treat = test_data$Am, stat = "smd.max",
                            pairwise = FALSE)

  #Moments the same as poly
  for (method in c("ipt", "ebal", "optweight")) {
    expect_no_condition({
      Wm <- weightit(Am ~ X1 + X2 + X3 + X4,
                     data = test_data, method = method,
                     include.obj = TRUE,
                     moments = c(X1 = 2, X2 = 3))
    })
    nm <- sprintf("Wm_%s", method)
    assign(nm, Wm)

    expect_no_condition({
      Wp <- weightit(Am ~ poly(X1, 2) + poly(X2, 3) + X3 + X4,
                     data = test_data, method = method,
                     include.obj = TRUE)
    })
    np <- sprintf("Wp_%s", method)
    assign(np, Wp)

    if (method %in% c("ipt", "ebal")) {
      expect_M_parts_okay(Wm, tolerance = eps)
      expect_M_parts_okay(Wp, tolerance = eps)
    }

    expect_equal(Wm$weights, Wp$weights,
                 label = sprintf("weights for %s", nm),
                 expected.label = sprintf("weights for %s", np),
                 tolerance = eps)


    expect_equal(cobalt::bal.compute(init1, weights = Wm$weights),
                 0,
                 label = sprintf("largest SMD for %s", nm),
                 tolerance = eps)
    expect_equal(cobalt::bal.compute(init1, weights = Wp$weights),
                 0,
                 label = sprintf("largest SMD for %s", np),
                 tolerance = eps)

    expect_not_equal(cobalt::bal.compute(init2, weights = Wm$weights),
                     0,
                     label = sprintf("largest SMD for %s with extra terms", nm),
                     tolerance = eps)
  }
})

test_that("moments works as expected, cont", {
  skip_if_not_installed("rootSolve")
  skip_if_not_installed("cobalt")

  eps <- if (capabilities("long.double")) 1e-5 else 1e-1

  test_data <- readRDS(test_path("fixtures", "test_data.rds"))

  poly_test <- cbind(test_data[c("X3", "X4")],
                     poly(test_data$X1, 2),
                     poly(test_data$X2, 3))

  init1 <- cobalt::bal.init(poly_test, treat = test_data$Ac, stat = "p.max")
  init2 <- cobalt::bal.init(poly(test_data$X1, 3), treat = test_data$Ac, stat = "p.max")

  init1t <- cobalt::bal.init(poly_test, stat = "smd.max")
  init2t <- cobalt::bal.init(poly(test_data$X1, 3), stat = "smd.max")

  #Moments the same as poly
  for (method in c("ebal", "cbps")) {
    expect_no_condition({
      Wm <- weightit(Ac ~ X1 + X2 + X3 + X4,
                     data = test_data, method = method,
                     include.obj = TRUE,
                     moments = c(X1 = 2, X2 = 3))
    })
    nm <- sprintf("Wm_%s", method)
    assign(nm, Wm)

    expect_no_condition({
      Wp <- weightit(Ac ~ poly(X1, 2) + poly(X2, 3) + X3 + X4,
                     data = test_data, method = method,
                     include.obj = TRUE)
    })
    np <- sprintf("Wp_%s", method)
    assign(np, Wp)

    if (method %in% c("ipt", "ebal")) {
      expect_M_parts_okay(Wm, tolerance = eps)
      expect_M_parts_okay(Wp, tolerance = eps)
    }

    expect_equal(Wm$weights, Wp$weights,
                 label = sprintf("weights for %s", nm),
                 expected.label = sprintf("weights for %s", np),
                 tolerance = eps)

    expect_equal(cobalt::bal.compute(init1, weights = Wm$weights),
                 0,
                 label = sprintf("largest corr for %s", nm),
                 tolerance = eps)
    expect_equal(cobalt::bal.compute(init1, weights = Wp$weights),
                 0,
                 label = sprintf("largest corr for %s", np),
                 tolerance = eps)

    expect_not_equal(cobalt::bal.compute(init2, weights = Wm$weights),
                     0,
                     label = sprintf("largest corr for %s with extra terms", nm),
                     tolerance = eps)

    if (method == "ebal") {
      expect_equal(cobalt::bal.compute(init1t, weights = Wm$weights),
                   0,
                   label = sprintf("largest target SMD for %s", nm),
                   tolerance = eps)
      expect_equal(cobalt::bal.compute(init1t, weights = Wp$weights),
                   0,
                   label = sprintf("largest target SMD for %s", np),
                   tolerance = eps)

      expect_not_equal(cobalt::bal.compute(init2t, weights = Wm$weights),
                       0,
                       label = sprintf("largest target SMD for %s with extra terms", nm),
                       tolerance = eps)
    }
  }
})