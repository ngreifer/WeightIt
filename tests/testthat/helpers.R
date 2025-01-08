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
    expect_failure(expect_equal(ESS(W$weights[W$treat == i]), sum(W$treat == i), ...))
  }
}