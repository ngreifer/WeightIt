expect_ATT_weights_okay <- function(W, focal = NULL, ...) {
  if (is_null(focal)) {
    focal <- W$focal
  }
  else {
    expect_equal(focal, W$focal)
  }

  expect_true(W$estimand %in% c("ATT", "ATC"))

  expect_equal(ESS(W$weights[W$treat == focal]),
               sum(W$treat == focal), ...)

  for (i in setdiff(unique(W$treat), focal)) {
    expect_not_equal(ESS(W$weights[W$treat == i]), sum(W$treat == i), ...)
  }
}

expect_not_equal <- function(object, expected, ...,
                             tolerance = if (edition_get() >= 3) testthat_tolerance(),
                             info = NULL, label = NULL, expected.label = NULL) {

  # if (!capabilities("long.double")) {
  #   return(invisible(NULL))
  # }

  act <- quasi_label(enquo(object), label, arg = "object")
  exp <- quasi_label(enquo(expected), expected.label, arg = "expected")

  if (edition_get() >= 3) {
    expect_waldo_not_equal("equal", act, exp, info, ..., tolerance = tolerance)
  }
  else {
    if (!is.null(tolerance)) {
      comp <- compare(act$val, exp$val, ..., tolerance = tolerance)
    }
    else {
      comp <- compare(act$val, exp$val, ...)
    }
    expect(!comp$equal, sprintf("%s equal to %s.\n%s", act$lab, exp$lab, comp$message), info = info)
    invisible(act$val)
  }
}

expect_waldo_not_equal <- function(type, act, exp, info, ...) {
  comp <- waldo::compare(act$val, exp$val, ..., x_arg = "actual",
                        y_arg = "expected")

  expect(length(comp) > 0, sprintf("%s (%s) is %s to %s (%s).\n\n%s",
                                    act$lab, "`actual`", type, exp$lab,
                                    "`expected`",
                                    paste0(comp, collapse = "\n\n")),
         info = info, trace_env = rlang::caller_env())
  invisible(act$val)
}

expect_balance_improved <- function(W, ...) {
  weighted <- abs(cobalt::col_w_smd(W$covs, W$treat, W$weights,
                                    s.weights = W$s.weights))
  unweighted <- abs(cobalt::col_w_smd(W$covs, W$treat,
                                      s.weights = W$s.weights))

  expect_true(max(weighted) < max(unweighted), ...)
}

inject_missingness <- function(data, cols, prop = 0.1, seed = 4321) {
  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  set.seed(seed)

  for (col in cols) {
    is.na(data[[col]]) <- sample(nrow(data), round(prop * nrow(data)))
  }

  if (is_null(old_seed)) rm(".Random.seed", envir = globalenv())
  else assign(".Random.seed", old_seed, envir = globalenv())

  data
}

# `method = "cbps"` with a continuous treatment cannot always solve its balance
# conditions -- at some samples the estimating equations have no root, and `weightit()`
# now warns rather than silently returning weights that do not balance. Whether it fires
# depends on the sample and the specification, so the affected tests cannot simply
# `expect_warning()`. This accepts that one warning and nothing else, so an unrelated
# condition still fails the test, and it keeps working unchanged if the warning stops
# firing.
expect_no_unexpected_warning <- function(expr, known = "could not be solved") {
  ws <- character()

  val <- withCallingHandlers(expr,
                             warning = function(w) {
                               ws <<- c(ws, conditionMessage(w))
                               invokeRestart("muffleWarning")
                             })

  expect_true(all(vapply(ws, grepl, logical(1L), pattern = known, fixed = TRUE)),
              label = sprintf("all warnings match \"%s\" (got: %s)", known,
                              paste(substr(ws, 1L, 40L), collapse = "; ")))

  invisible(val)
}
