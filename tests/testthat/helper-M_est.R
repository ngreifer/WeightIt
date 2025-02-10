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
  hess_treat.list <- lapply(Mparts.list, `[[`, "hess_treat")
  dw_dBtreat.list <- lapply(Mparts.list, `[[`, "dw_dBtreat")
  SW <- W$s.weights

  psi_treat <- function(Btreat.list, Xtreat.list, A.list, SW) {
    do.call("cbind", lapply(seq_along(Btreat.list), function(i) {
      psi_treat.list[[i]](Btreat.list[[i]], Xtreat.list[[i]], A.list[[i]], SW = SW)
    }))
  }

  wfun <- function(Btreat.list, Xtreat.list, A.list) {
    Reduce("*", lapply(seq_along(Btreat.list), function(i) {
      wfun.list[[i]](Btreat.list[[i]], Xtreat.list[[i]], A.list[[i]])
    }), init = 1)
  }

  psi <- function(B, Xtreat.list, A.list, SW) {
    Btreat.list <- btreat.list
    k <- 0
    for (i in seq_along(btreat.list)) {
      Btreat.list[[i]] <- B[k + seq_along(btreat.list[[i]])]
      k <- k + length(btreat.list[[i]])
    }

    psi_treat(Btreat.list, Xtreat.list, A.list, SW)
  }

  gradfun <- function(B, X, A, SW) {
    colSums(psi(B, X, A, SW))
  }

  start <- 1.01 * unlist(btreat.list)

  out <- rootSolve::multiroot(gradfun, start = start,
                              X = Xtreat.list, A = A.list, SW = SW,
                              maxiter = 1e5, rtol = 1e-8, atol = 1e-8, ctol = 1e-8)

  Btreat.list <- btreat.list
  k <- 0
  for (i in seq_along(btreat.list)) {
    Btreat.list[[i]] <- out$root[k + seq_along(btreat.list[[i]])]
    k <- k + length(btreat.list[[i]])
  }

  expect_equal(unlist(Btreat.list),
               unlist(btreat.list),
               ignore_attr = TRUE,
               tolerance = tolerance, ...)

  w <- wfun(Btreat.list, Xtreat.list, A.list)

  expect_equal(as.vector(w), as.vector(W$weights),
               ignore_attr = TRUE,
               tolerance = tolerance, ...)

  if (all(lengths(hess_treat.list) > 0)) {
    hess_numerical <- .gradient(gradfun, unlist(btreat.list),
                             X = Xtreat.list, A = A.list, SW = SW,
                             .method = "rich")

    hess_analytical <- .block_diag(lapply(seq_along(Mparts.list), function(i) {
      hess_treat.list[[i]](btreat.list[[i]], Xtreat.list[[i]], A.list[[i]], SW = SW)
    }))

    expect_equal(hess_analytical,
                 hess_numerical,
                 ignore_attr = TRUE,
                 tolerance = tolerance, ...)
  }

  if (all(lengths(dw_dBtreat.list) > 0)) {
    wfun2 <- function(B, Xtreat.list, A.list) {
      Btreat.list <- btreat.list
      k <- 0
      for (i in seq_along(btreat.list)) {
        Btreat.list[[i]] <- B[k + seq_along(btreat.list[[i]])]
        k <- k + length(btreat.list[[i]])
      }

      wfun(Btreat.list, Xtreat.list, A.list)
    }

    dwdb_numerical <- .gradient(wfun2, unlist(btreat.list),
                                Xtreat.list = Xtreat.list,
                                A.list = A.list)

    w.list <- c(lapply(seq_along(btreat.list), function(i) {
      wfun.list[[i]](btreat.list[[i]], Xtreat.list[[i]], A.list[[i]])
    }), list(rep(1, length(A.list[[1]]))))

    dwdb_analytical <- do.call("cbind", lapply(seq_along(btreat.list), function(i) {
      dw_dBtreat.list[[i]](btreat.list[[i]], X = Xtreat.list[[i]], A = A.list[[i]], SW = SW) *
        Reduce("*", w.list[-i])
        # wfun(btreat.list[-i], Xtreat.list[-i], A.list[-i])
    }))

    expect_equal(dwdb_analytical,
                 dwdb_numerical,
                 ignore_attr = TRUE,
                 tolerance = tolerance, ...)
  }

  invisible(list(solve = out,
                 b = out$root,
                 weights = w
  ))
}