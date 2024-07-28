#Original version of .multinom_weightit.fit() that uses rootSolve
..multinom_weightit.fit <- function(x, y, weights = NULL, offset = NULL, start = NULL, hess = TRUE, ...) {
  rlang::check_installed("rootSolve")

  chk::chk_atomic(y)
  y <- as.factor(y)
  chk::chk_numeric(x)
  chk::chk_matrix(x)

  nobs <- length(y)

  if (is_null(weights)) weights <- rep.int(1, nobs)
  else chk::chk_numeric(weights)

  if (is.null(offset)) offset <- rep.int(0, nobs)
  else chk::chk_numeric(offset)

  chk::chk_all_equal(c(length(y), nrow(x), length(weights)))

  QR <- qr(x)

  aliased_X <- !colnames(x) %in% colnames(make_full_rank(x, with.intercept = FALSE))
  aliased_B <- rep(aliased_X, nlevels(y) - 1)

  n <- length(y)
  k0 <- (nlevels(y) - 1) * ncol(x)
  k <- sum(!aliased_B)

  if (is_null(start)) {
    start <- rep(0, k0)
  }
  else {
    chk::chk_numeric(start)
    chk::chk_length(start, k0)
  }

  nm <- unlist(lapply(levels(y)[-1], function(i) paste(i, colnames(x), sep = "~")))
  names(start) <- nm

  x_ <- x[, !aliased_X, drop = FALSE]

  coef_ind <- setNames(lapply(seq_len(nlevels(y) - 1), function(i) {
    (i - 1) * sum(!aliased_X) + seq_len(sum(!aliased_X))
  }), levels(y)[-1])

  get_pp <- function(B, X, offset = NULL) {
    if (length(offset) == 0) offset <- 0
    qq <- lapply(levels(y)[-1], function(i) {
      exp(offset + drop(X %*% B[coef_ind[[i]]]))
    })

    pden <- 1 + rowSums(do.call("cbind", qq))

    pp <- do.call("cbind", c(list(1), qq)) / pden

    colnames(pp) <- levels(y)

    pp
  }

  #Multinomial logistic regression score
  psi <- function(B, X, y, weights, offset = NULL) {
    pp <- get_pp(B, X, offset)

    out <- do.call("cbind", lapply(levels(y)[-1], function(i) {
      weights * ((y == i) - pp[,i]) * X
    }))

    colnames(out) <- names(B)
    out
  }

  f <- function(B, X, y, weights, offset) {
    .colSums(psi(B, X, y, weights, offset), n, k)
  }

  out <- rootSolve::multiroot(f,
                              start = start[!aliased_B],
                              X = x_,
                              y = y,
                              weights = weights,
                              offset = offset,
                              ...)

  grad <- psi(out$root, X = x_, y = y,
              weights = weights, offset = offset)

  hessian <- NULL
  if (hess) {
    hessian <- gradient(f,
                        .x = out$root,
                        X = x_,
                        y = y, weights = weights, offset = offset)
    colnames(hessian) <- rownames(hessian) <- nm[!aliased_B]
  }

  pp <- get_pp(out$root, x_, offset)
  res <- rep(0, length(y))
  for (i in seq_len(nlevels(y))) {
    res[y == levels(y)[i]] <- 1 - pp[y == levels(y)[i], i]
  }

  coefs <- setNames(rep(NA_real_, length(start)), names(start))
  coefs[!aliased_B] <- out$root

  aliased <- rep(TRUE, ncol(x))
  aliased[QR$pivot[seq_len(QR$rank)]] <- FALSE
  attr(QR$qr, "aliased") <- aliased

  list(coefficients = coefs,
       residuals = res,
       fitted.values = pp,
       rank = QR$rank,
       qr = QR,
       solve = out,
       psi = psi,
       f = f,
       df.residual = length(res) - QR$rank,
       x = x,
       y = y,
       weights = weights,
       gradient = grad,
       hessian = hessian)
}
