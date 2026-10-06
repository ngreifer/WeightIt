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

test_that("get_w_from_ps() agrees with .get_w_from_ps_internal_bin() and .get_w_from_ps_internal_array() for binary", {
  skip_if_not_installed("patrick")

  # The random setups each reset the seed before drawing
  setups <- list(
    "internal_bin" = local({
      set.seed(1234)
      treat <- sample(0:1, 1e3, TRUE)
      ps <- runif(1e3)

      list(treat = treat, ps = ps)
    }),
    "internal_bin, PS 0/1" = local({
      treat <- rep(0:1, each = 2)
      ps <- rep(0:1, 2)

      list(treat = treat, ps = ps)
    }),
    "internal_array" = local({
      set.seed(1234)
      treat <- sample(0:1, 1e3, TRUE)
      ps <- matrix(runif(1e3 * 50), nrow = 1e3)

      list(treat = treat, ps = ps)
    }),
    "internal_array, PS 0/1" = local({
      set.seed(1234)
      treat <- sample(0:1, 1e3, TRUE)
      ps <- matrix(round(runif(1e3 * 50)), nrow = 1e3)

      #Do same adjustment that .get_w_from_ps_internal_array() does
      ps_ <- ps
      ps_[ps_ < 1e-8] <- 1e-8
      ps_[ps_ > 1 - 1e-8] <- 1 - 1e-8

      list(treat = treat, ps = ps, ps_ = ps_)
    })
  )

  # `input` is the propensity score given to get_w_from_ps(); the internal functions
  # always get the raw `ps`. A matrix of scores is passed column by column.
  patrick::with_parameters_test_that(
    "{setup}, estimand = {estimand}",
    {
      s <- setups[[setup]]

      if (is.matrix(s$ps)) {
        w_public <- apply(s[[input]], 2, get_w_from_ps, s$treat, estimand = estimand)
        w_internal <- .get_w_from_ps_internal_array(s$ps, s$treat, estimand)
      }
      else {
        w_public <- get_w_from_ps(s[[input]], s$treat, estimand = estimand)
        w_internal <- .get_w_from_ps_internal_bin(s$ps, s$treat, estimand)
      }

      expect_equal(w_public, w_internal)
    },
    .cases = rbind(expand.grid(estimand = c("ATE", "ATT", "ATC", "ATO", "ATM"),
                               setup = c("internal_bin", "internal_bin, PS 0/1",
                                         "internal_array"),
                               input = "ps",
                               stringsAsFactors = FALSE),
                   data.frame(estimand = c("ATE", "ATT", "ATC", "ATO", "ATM"),
                              setup = "internal_array, PS 0/1",
                              input = c("ps_", "ps_", "ps_", "ps", "ps_")))
  )
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

test_that("estimand = 'ATOS' is invariant to which level is treated", {
  skip_on_cran()

  # Crump et al.'s optimal subpopulation depends on the propensity score only
  # through e(1 - e), which is unchanged by swapping the treatment labels. The
  # candidate values of alpha searched are therefore the sub-.5 half of both
  # columns pooled; bounding the search by one column's count instead used to
  # truncate it, so a sample with few units below .5 got no trimming at all in
  # one orientation and substantial trimming in the other.
  set.seed(11)
  n <- 2000L
  x <- rnorm(n)
  e <- plogis(1.5 + 1.2 * x)
  A <- rbinom(n, 1L, e)

  expect_lt(mean(e < .5), .2)

  w  <- get_w_from_ps(e, A, estimand = "ATOS")
  ws <- get_w_from_ps(1 - e, 1 - A, estimand = "ATOS")

  expect_equal(attr(w, "alpha"), attr(ws, "alpha"))
  expect_identical(unname(w == 0), unname(ws == 0))

  # Something was actually dropped, so this is not vacuous
  expect_gt(sum(w == 0), 0L)

  # The retained units keep their ATE weights
  keep <- w > 0
  expect_equal(unname(w[keep]),
               unname(get_w_from_ps(e, A, estimand = "ATE")[keep]))
})
