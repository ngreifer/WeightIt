expect_M_parts_okay <- function(W, tolerance = 1e-5, ...) {
  Mparts.list <- {
    if (is_not_null(attr(W, "Mparts", exact = TRUE))) {
      Mparts.list <- list(attr(W, "Mparts"))
    }
    else if (is_not_null(attr(W, "Mparts.list", exact = TRUE))) {
      Mparts.list <- attr(W, "Mparts.list", exact = TRUE)
    }
    else {
      NULL
    }
  }

  expect_false(is_null(Mparts.list))

  psi_treat.list <- lapply(Mparts.list, `[[`, "psi_treat")
  wfun.list <- lapply(Mparts.list, `[[`, "wfun")
  Xtreat.list <- lapply(Mparts.list, `[[`, "Xtreat")
  A.list <- lapply(Mparts.list, `[[`, "A")
  btreat.list <- lapply(Mparts.list, `[[`, "btreat")
  SW <- W$s.weights

  psi_treat <- function(Btreat.list, A.list, Xtreat.list, SW) {
    do.call("cbind", lapply(seq_along(Btreat.list), function(i) {
      psi_treat.list[[i]](Btreat.list[[i]], A.list[[i]], Xtreat.list[[i]], SW)
    }))
  }

  wfun <- function(Btreat.list, A.list, Xtreat.list) {
    Reduce("*", lapply(seq_along(Btreat.list), function(i) {
      wfun.list[[i]](Btreat.list[[i]], Xtreat.list[[i]], A.list[[i]])
    }), init = 1)
  }

  psi <- function(B, A.list, Xtreat.list, SW) {
    Btreat.list <- btreat.list
    k <- 0
    for (i in seq_along(btreat.list)) {
      Btreat.list[[i]] <- B[k + seq_along(btreat.list[[i]])]
      k <- k + length(btreat.list[[i]])
    }

    psi_treat(Btreat.list, A.list, Xtreat.list, SW)
  }

  gradfun <- function(B, A, X, SW) {
    colSums(psi(B, A, X, SW))
  }

  start <- 1.01 * unlist(btreat.list)

  out <- rootSolve::multiroot(gradfun, start = start,
                              X = Xtreat.list, A = A.list, SW = SW,
                              maxiter = 1e5)

  Btreat.list <- btreat.list
  k <- 0
  for (i in seq_along(btreat.list)) {
    Btreat.list[[i]] <- out$root[k + seq_along(btreat.list[[i]])]
    k <- k + length(btreat.list[[i]])
  }

  w <- wfun(Btreat.list, A.list, Xtreat.list)

  expect_equal(unname(unlist(Btreat.list)), unname(unlist(btreat.list)),
               tolerance = tolerance, ...)

  expect_equal(unname(as.vector(w)), unname(as.vector(W$weights)),
               tolerance = tolerance, ...)

  invisible(list(solve = out,
                 b = out$root,
                 weights = w
  ))
}