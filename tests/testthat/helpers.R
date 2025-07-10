expect_ATT_weights_okay <- function(W, focal = NULL, ...) {
  if (is.null(focal)) {
    focal <- W$focal
  }
  else {
    expect_equal(focal, W$focal)
  }

  expect_true(W$estimand %in% c("ATT", "ATC"))

  expect_equal(ESS(W$weights[W$treat == focal]), sum(W$treat == focal), ...)
  for (i in setdiff(unique(W$treat), focal)) {
    expect_not_equal(ESS(W$weights[W$treat == i]), sum(W$treat == i), ...)
  }
}

expect_not_equal <- function(object, expected, ...,
                             tolerance = if (edition_get() >= 3) testthat_tolerance(),
                             info = NULL, label = NULL, expected.label = NULL) {

  if (!capabilities("long.double")) {
    return(NULL)
  }

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