# `summary()`, `plot()`, `nobs()`, `as.weightit()`, and `as.weightitMSM()` are all
# exported and were essentially untested: `summary.weightit()` was only ever called
# inside `expect_no_condition()`, `plot()` never at all, and `as.weightit()` never,
# despite carrying the largest block of untested input validation in the package.

skip_on_cran()

eps <- if (capabilities("long.double")) 1e-5 else 1e-3

test_data <- readRDS(test_path("fixtures", "test_data.rds"))

# ---- summary.weightit() ----------------------------------------------------

test_that("summary.weightit() reports the quantities it claims to", {
  d <- test_data

  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATE")

  s <- summary(W)

  expect_s3_class(s, "summary.weightit")
  expect_named(s, c("weight.range", "weight.top", "weight.mean", "coef.of.var",
                    "scaled.mad", "negative.entropy", "num.zeros",
                    "effective.sample.size"))

  # The ESS table matches `ESS()` computed by hand within each treatment group
  ess <- s$effective.sample.size
  expect_identical(dim(ess), c(2L, 2L))

  for (a in 0:1) {
    grp <- if (a == 1) "Treated" else "Control"
    in_a <- d$A == a
    expect_equal(ess["Unweighted", grp], sum(in_a), tolerance = eps)
    expect_equal(ess["Weighted", grp], ESS(W$weights[in_a]), tolerance = eps)
  }

  # Weighting always costs effective sample size
  expect_true(all(ess["Weighted", ] <= ess["Unweighted", ]))

  # The reported range brackets the actual weights
  expect_equal(unname(unlist(s$weight.range)),
               unname(c(range(W$weights[d$A == 0]), range(W$weights[d$A == 1]))[c(3, 4, 1, 2)]),
               tolerance = eps, ignore_attr = TRUE)

  expect_output(print(s), "Summary of weights")
  expect_identical(nobs(W), nrow(d))
})

test_that("summary.weightit() handles continuous treatments and zero weights", {
  d <- test_data

  # Continuous treatments get one "Total" column rather than one per group
  Wc <- weightit(Ac ~ X1 + X2 + X3, data = d, method = "glm")
  sc <- summary(Wc)

  expect_s3_class(sc, "summary.weightit")
  expect_identical(colnames(sc$effective.sample.size), "Total")
  expect_equal(sc$effective.sample.size["Weighted", "Total"],
               ESS(Wc$weights), tolerance = eps)
  expect_equal(sc$effective.sample.size["Unweighted", "Total"], nrow(d),
               tolerance = eps)
  expect_output(print(sc), "Summary of weights")

  # For a censoring model the weight distribution is deliberately summarized over the
  # units still under observation, so `num.zeros` is 0 even though the censored units
  # carry a weight of exactly 0. The loss they cause shows up in the ESS instead.
  set.seed(123L)
  d$C <- rbinom(nrow(d), 1L, prob = plogis(-0.9 + 0.8 * d$X1))
  Wcens <- weightit(.cens(C) ~ X1 + X2, data = d, method = "glm")

  s_cens <- summary(Wcens)
  expect_identical(sum(unlist(s_cens$num.zeros)), 0)
  expect_gt(min(unlist(s_cens$weight.range)), 0)
  expect_lt(s_cens$effective.sample.size["Weighted", "Total"], sum(d$C == 0L))
})

test_that("summary.weightitMSM() summarizes each time point", {
  d <- msmdata

  W <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                        A_2 ~ X1_1 + X2_1 + A_1),
                   data = d, method = "glm")

  s <- summary(W)

  expect_s3_class(s, "summary.weightitMSM")
  expect_length(s, 2L)
  expect_named(s, c("A_1", "A_2"))

  for (i in seq_along(s)) {
    expect_s3_class(s[[i]], "summary.weightit")
  }

  # Each heading gives the model's position in `formula.list`, whether it is a
  # treatment or a censoring model, and the variable it is about
  expect_output(print(s), "1. Treatment: A_1", fixed = TRUE)
  expect_output(print(s), "2. Treatment: A_2", fixed = TRUE)
})

test_that("summary.weightit() names its groups the way it labels them", {
  d <- test_data

  s <- summary(weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATE"))

  # Every element split by group carries the same names, which are also the row labels
  # of the printed tables. Pinned because a rename here turns `s$weight.range$<group>`
  # into `NULL` at every call site, which reads as an empty result rather than an error.
  grouped <- c("weight.range", "weight.top", "coef.of.var", "scaled.mad",
               "negative.entropy", "num.zeros")

  expect_identical(lapply(s[grouped], names),
                   setNames(rep(list(c("Treated", "Control")), length(grouped)),
                            grouped))

  out <- cli::ansi_strip(capture.output(print(s)))
  expect_match(out, "^Treated", all = FALSE)
  expect_match(out, "^Control", all = FALSE)

  # A treatment with no groups to split by is one group covering everyone
  sc <- summary(weightit(Ac ~ X1 + X2 + X3, data = d, method = "glm"))
  expect_identical(lapply(sc[grouped], names),
                   setNames(rep(list("All"), length(grouped)), grouped))

  # Stabilizing adds the mean of the weights, named the same way
  ss <- summary(weightit(A ~ X1 + X2 + X3, data = d, method = "glm", stabilize = TRUE))
  expect_named(ss$weight.mean, c("Treated", "Control"))
})

test_that("the weight-range plot marks a group whose weights are all equal", {
  d <- test_data

  # Under the ATT the focal group's weights are all exactly 1, so its range is a point.
  # A span drawn between two identical endpoints would be an empty pair of brackets, so
  # it gets a single mark instead.
  s <- summary(weightit(A ~ X1 + X2 + X3, data = d, method = "glm", estimand = "ATT"))

  expect_identical(s$weight.range$Treated, c(1, 1))

  out <- cli::ansi_strip(capture.output(print(s)))
  rows <- grep("^Treated |^Control ", out, value = TRUE)[1:2]

  # Asserted against the span character rather than the brackets, since which bracket is
  # drawn depends on whether Unicode is usable but whether a span is drawn at all does
  # not
  expect_no_match(rows[1L], cli::symbol$double_line, fixed = TRUE)
  expect_match(rows[2L], cli::symbol$double_line, fixed = TRUE)
})

test_that(".text_box_plot() follows cli's Unicode fallback", {
  #The span is a cli symbol; the brackets and the point marker are not, and used to stay
  #Unicode in output that had otherwise degraded to ASCII -- visible under `testthat`,
  #which sets `cli.unicode = FALSE`.
  cells <- function(unicode) {
    op <- options(cli.unicode = unicode)
    on.exit(options(op))

    d <- .text_box_plot(list(span = c(0, 1), point = c(0.5, 0.5)))

    setNames(d[[2L]], rownames(d))
  }

  utf8 <- cells(TRUE)

  expect_match(utf8[["span"]], "╞", fixed = TRUE)
  expect_match(utf8[["span"]], "═", fixed = TRUE)
  expect_match(utf8[["point"]], "│", fixed = TRUE)

  # Nothing outside the ASCII set survives, the brackets and marker included
  ascii <- cells(FALSE)

  expect_match(ascii[["span"]], "|", fixed = TRUE)
  expect_match(ascii[["span"]], "=", fixed = TRUE)
  expect_match(ascii[["point"]], "|", fixed = TRUE)
  expect_no_match(ascii, "[^ =|]")
})

test_that("print.summary.weightitMSM() heads the summary once and rules off each model", {
  W <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                        A_2 ~ X1_1 + X2_1 + A_1,
                        A_3 ~ X1_2 + X2_2 + A_2),
                   data = msmdata, method = "glm")

  out <- cli::ansi_strip(capture.output(print(summary(W))))

  # One header for the object as a whole, not one per model
  expect_length(grep("Summary of weights", out, fixed = TRUE), 1L)

  # A divider per model, each carrying that model's label. `cli::symbol$line` rather
  # than the character itself, since cli falls back to ASCII where Unicode is unusable.
  rule <- strrep(cli::symbol$line, 3L)
  rules <- which(startsWith(out, rule))

  expect_length(rules, 3L)
  expect_identical(sub(sprintf("^%s+ (.*?) %s+$", cli::symbol$line, cli::symbol$line),
                       "\\1", out[rules]),
                   c("1. Treatment: A_1", "2. Treatment: A_2", "3. Treatment: A_3"))

  # All drawn to one width, set by the widest line under any of them, so the models
  # read as one block rather than as a ragged set of tables
  expect_length(unique(nchar(out[rules])), 1L)

  body <- out[-seq_len(rules[1L] - 1L)]
  expect_identical(max(nchar(body)), nchar(out[rules][1L]))

  # Every model gets its own divider, including a lone one
  W1 <- weightitMSM(list(A_1 ~ X1_0 + X2_0), data = msmdata, method = "glm")
  out1 <- cli::ansi_strip(capture.output(print(summary(W1))))

  expect_length(which(startsWith(out1, rule)), 1L)
  expect_match(out1, "1. Treatment: A_1", fixed = TRUE, all = FALSE)
})

# ---- plot() ----------------------------------------------------------------

test_that("plot() works for the methods that support it and errors otherwise", {
  skip_if_not_installed("gbm")
  skip_if_not_installed("cobalt")

  d <- test_data

  # `plot.weightit_ok` is TRUE only for gbm and optweight; there it plots the
  # tuning/dual path
  Wg <- weightit(A ~ X1 + X2 + X3, data = d, method = "gbm",
                 criterion = "smd.mean", n.trees = 50)

  p <- plot(Wg)
  expect_s3_class(p, "ggplot")

  # Everything else says so rather than drawing something meaningless
  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm")
  expect_error(plot(W), "cannot be used")
})

test_that("plot.summary.weightit() returns a ggplot for both treatment types", {
  d <- test_data

  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm")
  expect_s3_class(plot(summary(W)), "ggplot")

  Wc <- weightit(Ac ~ X1 + X2 + X3, data = d, method = "glm")
  expect_s3_class(plot(summary(Wc)), "ggplot")

  Wm <- weightit(Am ~ X1 + X2 + X3, data = d, method = "glm")
  expect_s3_class(plot(summary(Wm)), "ggplot")
})

test_that("plot.summary.weightit() honors `bins` and `binwidth`", {
  d <- test_data

  Wc <- weightit(Ac ~ X1 + X2 + X3, data = d, method = "glm")
  s <- summary(Wc)

  nbins <- function(p) {
    length(unique(round(ggplot2::ggplot_build(p)$data[[1L]]$x, 8L)))
  }

  #`bins` used to be overwritten with 20 whenever `breaks` was not supplied, so
  #the user's value never reached geom_histogram()
  expect_gt(nbins(plot(s, bins = 60L)), nbins(plot(s, bins = 10L)))

  #`binwidth` is not clobbered by a default bin count
  w <- ggplot2::ggplot_build(plot(s, binwidth = 0.05))$data[[1L]]
  expect_equal(unique(round(w$xmax - w$xmin, 8L)), 0.05)
})

test_that("plot.summary.weightitMSM() returns a plot per time point", {
  W <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                        A_2 ~ X1_1 + X2_1 + A_1),
                   data = msmdata, method = "glm")

  p <- plot(summary(W))

  expect_s3_class(p, "ggplot")
})

# ---- as.weightit() --------------------------------------------------------

test_that("as.weightit() promotes a weightit.fit object", {
  d <- test_data

  covs <- d[c("X1", "X2", "X3")]

  WF <- weightit.fit(as.matrix(covs), treat = d$A, method = "glm",
                     estimand = "ATE")
  W <- as.weightit(WF, covs = covs)

  expect_s3_class(W, "weightit")
  expect_equal(unname(W$weights), unname(WF$weights), tolerance = eps)
  expect_identical(W$estimand, "ATE")
  expect_equal(unname(W$ps), unname(WF$ps), tolerance = eps)
  expect_named(W$covs, names(covs))

  # It behaves like a real weightit object downstream
  expect_no_error(summary(W))
  expect_output(print(W), "weightit object")
  expect_no_error(vcov(glm_weightit(Y_B ~ A, data = d, weightit = W,
                                    family = binomial)))

  # And matches what weightit() would have produced from the same model
  W_direct <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm",
                       estimand = "ATE")
  expect_equal(unname(W$weights), unname(W_direct$weights), tolerance = eps)
})

test_that("as.weightit.default() wraps a bare weight vector", {
  d <- test_data

  w <- runif(nrow(d), .5, 2)

  W <- as.weightit(w, treat = d$A, covs = d[c("X1", "X2")],
                   estimand = "ATE", s.weights = d$SW)

  expect_s3_class(W, "weightit")
  expect_equal(unname(W$weights), w, tolerance = eps)
  expect_identical(W$estimand, "ATE")
  expect_equal(unname(W$s.weights), unname(d$SW), tolerance = eps)
  expect_identical(get_treat_type(W$treat), "binary")

  expect_no_error(summary(W))

  # Weights supplied this way were not estimated, so there is nothing for
  # M-estimation to account for
  expect_null(attr(W, "Mparts", exact = TRUE))

  # Minimal input works too
  expect_s3_class(as.weightit(w, treat = d$Ac), "weightit")
})

test_that("as.weightit() validates its input", {
  d <- test_data

  w <- runif(nrow(d))

  expect_error(as.weightit(w))
  expect_error(as.weightit("a", treat = d$A), "numeric vector")
  expect_error(as.weightit(matrix(w, ncol = 2L), treat = d$A), "numeric vector")
  expect_error(as.weightit(w, treat = d$A[-1L]), "same length")
  expect_error(as.weightit(w, treat = d$A, covs = d[1:10, c("X1", "X2")]))
  expect_error(as.weightit(w, treat = d$A, s.weights = d$SW[-1L]))
})

test_that("as.weightitMSM() wraps weights for a longitudinal treatment", {
  d <- msmdata

  w <- runif(nrow(d), .5, 2)

  W <- as.weightitMSM(w, treat.list = list(d$A_1, d$A_2, d$A_3),
                      covs.list = list(d[c("X1_0", "X2_0")],
                                       d[c("X1_1", "X2_1")],
                                       d[c("X1_2", "X2_2")]))

  expect_s3_class(W, "weightitMSM")
  expect_equal(unname(W$weights), w, tolerance = eps)
  expect_length(W$treat.list, 3L)

  expect_no_error(summary(W))
  expect_output(print(W), "weightitMSM object")

  expect_error(as.weightitMSM(w))
  expect_error(as.weightitMSM(w, treat.list = list(d$A_1, d$A_2[-1L])))
})

test_that("as.weightitMSM() takes a censoring indicator among the treatments", {
  #Censoring is a treatment type, so an indicator goes in `treat.list` like any other
  #entry. Such an object has no `at.risk`, which `weightitMSM()` builds while fitting, so
  #printing and summarizing must read the risk set off the indicator instead: a unit
  #censored earlier has no indicator here and was never eligible.
  d <- msmdata

  C <- rbinom(nrow(d), 1L, .2)
  A_2 <- d$A_2
  is.na(A_2[C == 1L]) <- TRUE

  w <- runif(nrow(d), .5, 2)
  w[C == 1L] <- 0

  W <- as.weightitMSM(w, treat.list = list(A_1 = d$A_1, C = .cens(C), A_2 = A_2),
                      covs.list = list(d[c("X1_0", "X2_0")],
                                       d[c("X1_0", "X2_0")],
                                       d[c("X1_1", "X2_1")]))

  expect_null(W$at.risk)
  expect_identical(unname(which(is_cens_treat(W$treat.list))), 2L)

  expect_output(print(W),
                sprintf("time 2 \\(C\\): censoring \\(IPCW\\); %s of %s units censored",
                        sum(C == 1L), nrow(d)))

  s <- summary(W)
  expect_named(s, c("A_1", "C", "A_2"))
  expect_output(print(s), "2. Censoring: C")

  #The censoring entry covers the uncensored units, whose weights are all nonzero
  expect_identical(unname(s$C$num.zeros), 0)
})

# ---- tidy() and glance() ---------------------------------------------------

test_that("tidy() and glance() work for every outcome-model class", {
  skip_if_not_installed("patrick")

  d <- test_data
  d$Y_M <- factor(d$Y_O, ordered = FALSE)
  q <- quantile(d$Y_S, .8)
  d$time <- pmin(d$Y_S, q)
  d$event <- as.integer(d$Y_S <= q)

  W <- weightit(A ~ X1 + X2 + X3, data = d, method = "glm")

  patrick::with_parameters_test_that(
    "tidy() and glance() work for {nm}_weightit",
    {
      if (nm == "coxph") {
        skip_if_not_installed("survival")
      }

      fit <- switch(nm,
                    glm = glm_weightit(Y_B ~ A + X1, data = d, weightit = W,
                                       family = binomial),
                    multinom = multinom_weightit(Y_M ~ A, data = d, weightit = W),
                    ordinal = ordinal_weightit(Y_O ~ A, data = d, weightit = W),
                    coxph = coxph_weightit(survival::Surv(time, event) ~ A, data = d,
                                           weightit = W))

      tid <- generics::tidy(fit)
      expect_s3_class(tid, "tbl_df")
      expect_named(tid, c("term", "estimate", "std.error", "statistic", "p.value"))
      expect_identical(nrow(tid), length(coef(fit)))
      expect_equal(tid$estimate, unname(coef(fit)), tolerance = eps,
                   ignore_attr = TRUE)
      expect_equal(tid$std.error, unname(sqrt(diag(vcov(fit)))), tolerance = eps,
                   ignore_attr = TRUE)

      # `conf.int = TRUE` adds properly named limits that bracket the estimate
      tid_ci <- generics::tidy(fit, conf.int = TRUE)
      expect_named(tid_ci, c("term", "estimate", "std.error", "statistic",
                             "p.value", "conf.low", "conf.high"))
      expect_true(all(tid_ci$conf.low <= tid_ci$estimate),
                  label = sprintf("conf.low below estimate for %s", nm))
      expect_true(all(tid_ci$conf.high >= tid_ci$estimate),
                  label = sprintf("conf.high above estimate for %s", nm))

      # `exponentiate = TRUE` transforms the estimate and its limits together. It also
      # drops the standard error, since that is not the transform of the original, so the
      # remaining columns must still be labeled correctly rather than shifted by one.
      tid_exp <- generics::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
      expect_named(tid_exp, c("term", "estimate", "statistic", "p.value",
                              "conf.low", "conf.high"))
      expect_equal(tid_exp$estimate, exp(tid$estimate), tolerance = eps)
      expect_equal(tid_exp$statistic, tid$statistic, tolerance = eps)
      expect_equal(tid_exp$p.value, tid$p.value, tolerance = eps)
      expect_true(all(tid_exp$conf.low > 0))
      expect_equal(tid_exp$conf.low, exp(tid_ci$conf.low), tolerance = eps)

      gl <- generics::glance(fit)
      expect_s3_class(gl, "tbl_df")
      expect_identical(gl$nobs, nobs(fit))
    },
    nm = c("glm", "multinom", "ordinal", "coxph")
  )
})
