# Multinomial logistic regression
.multinom_weightit.fit <- function(x, y, weights = NULL, offset = NULL, start = NULL,
                                   hess = TRUE, control = list(), ...) {
  chk::chk_atomic(y)
  y <- as.factor(y)
  chk::chk_numeric(x)
  chk::chk_matrix(x)

  if (is_null(colnames(x))) {
    colnames(x) <- paste0("x", seq_col(x))
  }

  N <- length(y)

  if (is_null(weights)) weights <- rep.int(1, N)
  else chk::chk_numeric(weights)

  if (is.null(offset)) offset <- rep.int(0, N)
  else chk::chk_numeric(offset)

  chk::chk_all_equal(c(length(y), nrow(x), length(weights), length(offset)))

  K <- nlevels(y) - 1L

  aliased_X <- colnames(x) %nin% colnames(make_full_rank(x, with.intercept = FALSE))
  aliased_B <- rep.int(aliased_X, K)

  k0 <- K * ncol(x)

  if (is_null(start)) {
    start <- rep.int(0, k0)
  }
  else {
    chk::chk_numeric(start)
    chk::chk_length(start, k0)
  }

  nm <- unlist(lapply(levels(y)[-1L], function(i) paste(i, colnames(x), sep = "~")))
  names(start) <- nm

  x_ <- x[, !aliased_X, drop = FALSE]

  get_pp <- function(B, X, offset = NULL) {
    if (length(offset) == 0L) offset <- 0

    qq <- exp(offset + X %*% matrix(B, nrow = ncol(X)))

    pp <- cbind(1, qq) / (1 + rowSums(qq))

    colnames(pp) <- levels(y)
    rownames(pp) <- rownames(X)

    pp
  }

  #Multinomial logistic regression score
  psi <- function(B, X, y, weights, offset = NULL) {
    pp <- get_pp(B, X, offset)

    out <- do.call("cbind", lapply(levels(y)[-1L], function(i) {
      weights * ((y == i) - pp[, i]) * X
    }))

    if (is_not_null(names(B)))
      colnames(out) <- names(B)

    out
  }

  gr <- function(B, X, y, weights, offset) {
    colSums(psi(B, X, y, weights, offset))
  }

  ind_mat <- cbind(seq_along(y), as.integer(y))

  ll <- function(B, X, y, weights, offset) {
    p <- get_pp(B, X, offset)[ind_mat]

    sum(weights * log(p))
  }

  m_control <- list(fnscale = -1, #maximize likelihood; optim() minimizes by default
                  trace = 0,
                  maxit = 1e3L,
                  reltol = 1e-12)

  control <- utils::modifyList(m_control, control)

  out <- optim(par = start[!aliased_B],
               ll,
               X = x_,
               y = y,
               weights = weights,
               offset = offset,
               gr = gr,
               method = "BFGS",
               control = control)

  grad <- psi(out$par, X = x_, y = y,
              weights = weights, offset = offset)

  pp <- get_pp(out$par, x_, offset)

  res <- setNames(1 - pp[ind_mat], rownames(x))

  coefs <- rep_with(NA_real_, start)
  coefs[!aliased_B] <- out$par

  fit <- list(coefficients = coefs,
              residuals = res,
              fitted.values = pp,
              solve = out,
              psi = psi,
              f = gr,
              get_p = get_pp,
              df.residual = length(res) - sum(!is.na(coefs)),
              x = x,
              y = y,
              weights = weights,
              gradient = grad)

  if (hess) {
    hessian <- matrix(NA_real_, nrow = sum(!aliased_B), ncol = sum(!aliased_B))

    for (i in seq_len(K)) {
      i_ind <- (i - 1L) * ncol(x_) + seq_len(ncol(x_))
      for (j in seq_len(i)) {
        if (i == j) {
          hessian[i_ind, i_ind] <- -crossprod(x_ * ((1 - pp[, i + 1L]) * pp[, i + 1L] * weights), x_)
        }
        else {
          j_ind <- (j - 1L) * ncol(x_) + seq_len(ncol(x_))

          hessian[i_ind, j_ind] <- -crossprod(x_ * (-pp[, i + 1L] * pp[, j + 1L] * weights), x_)
          hessian[j_ind, i_ind] <- t(hessian[i_ind, j_ind])
        }
      }
    }

    colnames(hessian) <- rownames(hessian) <- nm[!aliased_B]

    fit$hessian <- hessian
  }

  fit
}

.multinom_weightit <- function(formula, data, weights, subset, start = NULL, na.action,
                               hess = TRUE, control = list(), model = TRUE,
                               x = FALSE, y = TRUE, contrasts = NULL, ...) {
  cal <- match.call()

  chk::chk_flag(hess)
  chk::chk_flag(model)
  chk::chk_flag(x)
  chk::chk_flag(y)

  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- {
    if (is.empty.model(mt)) matrix(NA_real_, NROW(Y), 0L)
    else model.matrix(mt, mf, contrasts)
  }

  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    .err("`weights` must be a numeric vector")

  if (!is.null(weights) && any(weights < 0))
    .err("negative weights not allowed")

  offset <- as.vector(model.offset(mf))
  if (is_not_null(offset) && length(offset) != NROW(Y)) {
      .err(gettextf("number of offsets is %d; should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }

  fit <- eval(call(".multinom_weightit.fit",
                   x = X, y = Y, weights = weights,
                   offset = offset, start = start,
                   hess = hess, control = control))

  if (model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (!x) fit$x <- NULL
  if (!y) fit$y <- NULL

  c(fit,
    list(call = cal, formula = formula, terms = mt,
         data = data, offset = offset,
         contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)))
}

.get_hess_multinom <- function(fit) {
  x <- if_null_then(fit[["x"]], model.matrix(fit))
  y <- if_null_then(fit[["y"]], model.response(model.frame(fit)))
  weights <- weights(fit)
  coefs <- coef(fit)

  y <- as.factor(y)

  if (is_null(colnames(x))) {
    colnames(x) <- paste0("x", seq_col(x))
  }

  if (is_null(weights)) {
    weights <- rep.int(1, length(y))
  }

  K <- nlevels(y) - 1L

  aliased_X <- colnames(x) %nin% colnames(make_full_rank(x, with.intercept = FALSE))

  x_ <- x[, !aliased_X, drop = FALSE]

  theta0 <- na.rem(coefs)

  pp <- fit$fitted.values

  hessian <- matrix(NA_real_, nrow = length(theta0), ncol = length(theta0))

  for (i in seq_len(K)) {
    i_ind <- (i - 1L) * ncol(x_) + seq_len(ncol(x_))
    for (j in seq_len(i)) {
      if (i == j) {
        hessian[i_ind, i_ind] <- -crossprod(x_ * ((1 - pp[, i + 1L]) * pp[, i + 1L] * weights), x_)
      }
      else {
        j_ind <- (j - 1L) * ncol(x_) + seq_len(ncol(x_))

        hessian[i_ind, j_ind] <- -crossprod(x_ * (-pp[, i + 1L] * pp[, j + 1L] * weights), x_)
        hessian[j_ind, i_ind] <- t(hessian[i_ind, j_ind])
      }
    }
  }

  colnames(hessian) <- rownames(hessian) <- names(theta0)

  hessian
}
