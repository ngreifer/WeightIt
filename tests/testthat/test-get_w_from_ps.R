test_that("get_w_from_ps() works for binary", {
  set.seed(1234)
  treat <- sample(0:1, 1e3, TRUE)
  ps <- runif(1e3)

  # w_ate <- rep(1, length(treat))
  # w_ate[treat == 1] <- 1 / ps[treat == 1]
  # w_ate[treat == 0] <- 1 / (1 - ps[treat == 0])

  w_ate <- treat / ps + (1 - treat) / (1 - ps)

  # w_att <- rep(1, length(treat))
  # w_att[treat == 0] <- ps[treat == 0] / (1 - ps[treat == 0])

  w_att <- ps * w_ate

  # w_atc <- rep(1, length(treat))
  # w_atc[treat == 1] <- (1 - ps[treat == 1]) / ps[treat == 1]

  w_atc <- (1 - ps) * w_ate

  # w_ato <- treat * (1 - ps) + (1 - treat) * ps

  w_ato <- w_ate * ps * (1 - ps)

  # w_atm <- rep(1, length(treat))
  # w_atm[treat == 1 & ps > .5] <- (1 - ps[treat == 1 & ps > .5]) / ps[treat == 1 & ps > .5]
  # w_atm[treat == 0 & ps < .5] <- ps[treat == 0 & ps < .5] / (1 - ps[treat == 0 & ps < .5])

  w_atm <- w_ate * pmin(ps, 1 - ps)

  expect_equal(get_w_from_ps(ps, treat, estimand = "ATE"), w_ate)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATT"), w_att)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATC"), w_atc)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATO"), w_ato)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATM"), w_atm)
})

test_that("get_w_from_ps() works for binary, PS 0/1", {
  treat <- rep(0:1, each = 2)
  ps <- rep(0:1, 2)

  w_ate <- rep(1, length(treat))
  w_ate[treat == 1] <- 1 / ps[treat == 1]
  w_ate[treat == 0] <- 1 / (1 - ps[treat == 0])

  w_att <- rep(1, length(treat))
  w_att[treat == 0] <- ps[treat == 0] / (1 - ps[treat == 0])

  w_atc <- rep(1, length(treat))
  w_atc[treat == 1] <- (1 - ps[treat == 1]) / ps[treat == 1]

  w_ato <- treat * (1 - ps) + (1 - treat) * ps

  w_atm <- rep(1, length(treat))
  w_atm[treat == 1 & ps > .5] <- (1 - ps[treat == 1 & ps > .5]) / ps[treat == 1 & ps > .5]
  w_atm[treat == 0 & ps < .5] <- ps[treat == 0 & ps < .5] / (1 - ps[treat == 0 & ps < .5])

  expect_equal(get_w_from_ps(ps, treat, estimand = "ATE"), w_ate)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATT"), w_att)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATC"), w_atc)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATO"), w_ato)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATM"), w_atm)
})

test_that("get_w_from_ps() agrees with .get_w_from_ps_internal_bin() for binary", {
  set.seed(1234)
  treat <- sample(0:1, 1e3, TRUE)
  ps <- runif(1e3)

  w_ate <- .get_w_from_ps_internal_bin(ps, treat, "ATE")
  w_att <- .get_w_from_ps_internal_bin(ps, treat, "ATT")
  w_atc <- .get_w_from_ps_internal_bin(ps, treat, "ATC")
  w_ato <- .get_w_from_ps_internal_bin(ps, treat, "ATO")
  w_atm <- .get_w_from_ps_internal_bin(ps, treat, "ATM")

  expect_equal(get_w_from_ps(ps, treat, estimand = "ATE"), w_ate)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATT"), w_att)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATC"), w_atc)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATO"), w_ato)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATM"), w_atm)
})

test_that("get_w_from_ps() agrees with .get_w_from_ps_internal_bin() for binary, PS 0/1", {
  treat <- rep(0:1, each = 2)
  ps <- rep(0:1, 2)

  w_ate <- .get_w_from_ps_internal_bin(ps, treat, "ATE")
  w_att <- .get_w_from_ps_internal_bin(ps, treat, "ATT")
  w_atc <- .get_w_from_ps_internal_bin(ps, treat, "ATC")
  w_ato <- .get_w_from_ps_internal_bin(ps, treat, "ATO")
  w_atm <- .get_w_from_ps_internal_bin(ps, treat, "ATM")

  expect_equal(get_w_from_ps(ps, treat, estimand = "ATE"), w_ate)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATT"), w_att)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATC"), w_atc)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATO"), w_ato)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATM"), w_atm)
})

test_that("get_w_from_ps() agrees with .get_w_from_ps_internal_array() for binary", {
  set.seed(1234)
  treat <- sample(0:1, 1e3, TRUE)
  ps <- matrix(runif(1e3 * 50), nrow = 1e3)

  w_ate <- .get_w_from_ps_internal_array(ps, treat, "ATE")
  w_att <- .get_w_from_ps_internal_array(ps, treat, "ATT")
  w_atc <- .get_w_from_ps_internal_array(ps, treat, "ATC")
  w_ato <- .get_w_from_ps_internal_array(ps, treat, "ATO")
  w_atm <- .get_w_from_ps_internal_array(ps, treat, "ATM")

  w_ate2 <- apply(ps, 2, get_w_from_ps, treat, estimand = "ATE")
  w_att2 <- apply(ps, 2, get_w_from_ps, treat, estimand = "ATT")
  w_atc2 <- apply(ps, 2, get_w_from_ps, treat, estimand = "ATC")
  w_ato2 <- apply(ps, 2, get_w_from_ps, treat, estimand = "ATO")
  w_atm2 <- apply(ps, 2, get_w_from_ps, treat, estimand = "ATM")

  expect_equal(w_ate, w_ate2)
  expect_equal(w_att, w_att2)
  expect_equal(w_atc, w_atc2)
  expect_equal(w_ato, w_ato2)
  expect_equal(w_atm, w_atm2)
})

test_that("get_w_from_ps() agrees with .get_w_from_ps_internal_array() for binary, PS 0/1", {
  set.seed(1234)
  treat <- sample(0:1, 1e3, TRUE)
  ps <- matrix(round(runif(1e3 * 50)), nrow = 1e3)

  w_ate <- .get_w_from_ps_internal_array(ps, treat, "ATE")
  w_att <- .get_w_from_ps_internal_array(ps, treat, "ATT")
  w_atc <- .get_w_from_ps_internal_array(ps, treat, "ATC")
  w_ato <- .get_w_from_ps_internal_array(ps, treat, "ATO")
  w_atm <- .get_w_from_ps_internal_array(ps, treat, "ATM")

  #Do same adjustment that .get_w_from_ps_internal_array() does
  ps_ <- ps
  ps_[ps_ < 1e-8] <- 1e-8
  ps_[ps_ > 1 - 1e-8] <- 1 - 1e-8

  expect_equal(apply(ps_, 2, get_w_from_ps, treat, estimand = "ATE"), w_ate)
  expect_equal(apply(ps_, 2, get_w_from_ps, treat, estimand = "ATT"), w_att)
  expect_equal(apply(ps_, 2, get_w_from_ps, treat, estimand = "ATC"), w_atc)
  expect_equal(apply(ps, 2, get_w_from_ps, treat, estimand = "ATO"), w_ato)
  expect_equal(apply(ps_, 2, get_w_from_ps, treat, estimand = "ATM"), w_atm)
})

test_that("get_w_from_ps() works for multi-category", {
  set.seed(1234)
  treat <- factor(LETTERS[sample(1:4, 1e3, TRUE)])
  ps <- matrix(runif(1e3 * nlevels(treat)), ncol = nlevels(treat),
               dimnames = list(NULL, levels(treat)))
  ps <- ps / rowSums(ps)

  w_ate <- 1 / ps[cbind(1:length(treat), match(treat, levels(treat)))]

  w_att <- rep(1, length(treat))
  for (i in levels(treat)[-1]) {
    w_att[treat == i] <- ps[treat == i, levels(treat)[1]] / ps[treat == i, i]
  }

  w_ato <- w_ate / rowSums(1 / ps)

  w_atm <- rep(1, length(treat))
  min_ind <- max.col(-ps)
  no_match <- treat != levels(treat)[min_ind]
  w_atm[no_match] <- w_ate[no_match] * ps[cbind(which(no_match), min_ind[no_match])]

  expect_equal(get_w_from_ps(ps, treat, estimand = "ATE"), w_ate)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATT", focal = levels(treat)[1]), w_att)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATO"), w_ato)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATM"), w_atm)
})

test_that("get_w_from_ps() works for multi-category, PS 0/1", {
  set.seed(1234)
  treat <- factor(LETTERS[sample(1:4, 1e3, TRUE)])
  ps <- matrix(runif(1e3 * nlevels(treat)), ncol = nlevels(treat),
               dimnames = list(NULL, levels(treat)))
  ps <- ps / rowSums(ps)
  for (i in 1:nrow(ps)) {
    ps[i,] <- as.numeric(ps[i,] == max(ps[i,]))
  }

  w_ate <- 1 / ps[cbind(1:length(treat), match(treat, levels(treat)))]

  w_att <- rep(1, length(treat))
  for (i in levels(treat)[-1]) {
    w_att[treat == i] <- ps[treat == i, levels(treat)[1]] / ps[treat == i, i]
  }

  w_ato <- w_ate / rowSums(1 / ps)

  w_atm <- rep(1, length(treat))
  min_ind <- max.col(-ps, ties.method = "first")
  no_match <- which(ps[cbind(seq_along(treat), match(treat, levels(treat)))] != ps[cbind(seq_along(treat), min_ind)])

  w_atm[no_match] <- w_ate[no_match] * ps[cbind(no_match, min_ind[no_match])]

  expect_equal(get_w_from_ps(ps, treat, estimand = "ATE"), w_ate)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATT", focal = levels(treat)[1]), w_att)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATO"), w_ato)
  expect_equal(get_w_from_ps(ps, treat, estimand = "ATM"), w_atm)
})
