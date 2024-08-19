# Multinomial logistic regression
.multinom_weightit.fit <- function(x, y, weights = NULL, offset = NULL, start = NULL, hess = TRUE, ...) {
  chk::chk_atomic(y)
  y <- as.factor(y)
  chk::chk_numeric(x)
  chk::chk_matrix(x)

  if (is_null(colnames(x))) {
    colnames(x) <- paste0("x", seq_col(x))
  }

  nobs <- length(y)

  if (is_null(weights)) weights <- rep.int(1, nobs)
  else chk::chk_numeric(weights)

  if (is.null(offset)) offset <- rep.int(0, nobs)
  else chk::chk_numeric(offset)

  chk::chk_all_equal(c(length(y), nrow(x), length(weights), length(offset)))

  aliased_X <- colnames(x) %nin% colnames(make_full_rank(x, with.intercept = FALSE))
  aliased_B <- rep.int(aliased_X, nlevels(y) - 1)

  k0 <- (nlevels(y) - 1) * ncol(x)

  if (is_null(start)) {
    start <- rep.int(0, k0)
  }
  else {
    chk::chk_numeric(start)
    chk::chk_length(start, k0)
  }

  nm <- unlist(lapply(levels(y)[-1], function(i) paste(i, colnames(x), sep = "~")))
  names(start) <- nm

  x_ <- x[, !aliased_X, drop = FALSE]

  get_pp <- function(B, X, offset = NULL) {
    if (length(offset) == 0) offset <- 0

    qq <- exp(offset + X %*% matrix(B, nrow = ncol(X)))

    pp <- cbind(1, qq) / (1 + rowSums(qq))

    colnames(pp) <- levels(y)
    rownames(pp) <- rownames(X)

    pp
  }

  #Multinomial logistic regression score
  psi <- function(B, X, y, weights, offset = NULL) {
    pp <- get_pp(B, X, offset)

    out <- do.call("cbind", lapply(levels(y)[-1], function(i) {
      weights * ((y == i) - pp[,i]) * X
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

  out <- optim(par = start[!aliased_B],
               ll,
               X = x_,
               y = y,
               weights = weights,
               offset = offset,
               gr = gr,
               method = "BFGS",
               hessian = hess,
               control = list(fnscale = -1, #maximize likelihood; optim() minimizes by default
                              maxit = 1e3,
                              reltol = 1e-14))

  grad <- psi(out$par, X = x_, y = y,
              weights = weights, offset = offset)

  hessian <- NULL
  if (hess) {
    hessian <- out$hessian
    colnames(hessian) <- rownames(hessian) <- nm[!aliased_B]
  }

  pp <- get_pp(out$par, x_, offset)

  res <- setNames(1 - pp[ind_mat], rownames(x))

  coefs <- setNames(rep.int(NA_real_, length(start)), names(start))
  coefs[!aliased_B] <- out$par

  list(coefficients = coefs,
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
       gradient = grad,
       hessian = hessian)
}

.multinom_weightit <- function(formula, data, weights, subset, start = NULL, na.action,
                             hess = TRUE, model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...) {
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
    if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
    else matrix(NA_real_, NROW(Y), 0L)
  }

  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights))
    .err("`weights` must be a numeric vector")

  if (!is.null(weights) && any(weights < 0))
    .err("negative weights not allowed")

  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y))
      .err(gettextf("number of offsets is %d; should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }

  fit <- eval(call(".multinom_weightit.fit",
                   x = X, y = Y, weights = weights,
                   offset = offset, start = start,
                   hess = hess))

  if (model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (!x) fit$x <- NULL
  if (!y) fit$y <- NULL

  c(fit,
    list(call = cal, formula = formula, terms = mt,
         data = data, offset = offset,
         contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)))
}

