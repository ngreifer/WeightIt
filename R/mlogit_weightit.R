# Multinomial logistic regression by finding roots of score equations
.mlogit_weightit.fit <- function(x, y, weights = NULL, offset = NULL, start = NULL, hess = TRUE, ...) {
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

  gr <- function(B, X, y, weights, offset) {
    .colSums(psi(B, X, y, weights, offset), n, k)
  }

  w_y_mat <- weights * do.call("cbind", lapply(levels(y), function(i) y == i))

  ll <- function(B, X, y, weights, offset) {
    pp <- get_pp(B, X, offset)

    sum(w_y_mat * log(pp))
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
  res <- rep(0, length(y))
  for (i in seq_len(nlevels(y))) {
    res[y == levels(y)[i]] <- 1 - pp[y == levels(y)[i], i]
  }

  coefs <- setNames(rep(NA_real_, length(start)), names(start))
  coefs[!aliased_B] <- out$par

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
       f = gr,
       get_p = get_pp,
       df.residual = length(res) - QR$rank,
       x = x,
       y = y,
       weights = weights,
       gradient = grad,
       hessian = hessian)
}

.mlogit_weightit <- function(formula, data, weights, subset, start = NULL, na.action,
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
      .err(gettextf("number of offsets is %d should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }

  fit <- eval(call(".mlogit_weightit.fit",
                   x = X, y = Y, weights = weights,
                   offset = offset, start = start,
                   hess = hess))

  if (model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (!x) fit$x <- NULL
  if (!y) fit$y <- NULL

  structure(c(fit, list(call = cal, formula = formula, terms = mt,
                        data = data, offset = offset,
                        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf))),
            class = c("mlogit_weightit"))
}

#' @exportS3Method vcov mlogit_weightit
vcov.mlogit_weightit <- function(object, ...) {
  if (is_null(object$hessian)) {
    .err("`hess = TRUE` must be specified in the original fitting call to use `vcov()`")
  }

  solve(-object$hessian)
}

#' @exportS3Method predict mlogit_weightit
predict.mlogit_weightit <- function(object, newdata = NULL, na.action = na.pass, ...) {

  na.act <- object$na.action
  object$na.action <- NULL

  if (is_null(newdata)) {
    pred <- object$fitted.values
    if (!is.null(na.act))
      pred <- napredict(na.act, pred)
    return(pred)
  }

  tt <- terms(object)

  Terms <- delete.response(tt)
  m <- model.frame(tt, newdata, na.action = na.action,
                   xlev = object$xlevels)
  y <- as.factor(model.response(m))
  if (!is.null(cl <- attr(Terms, "dataClasses")))
    .checkMFClasses(cl, m)
  x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)

  offset <- model.offset(m)
  if (!is.null(addO <- object$call$offset)) {
    addO <- eval(addO, newdata, environment(tt))
    offset <- {
      if (length(offset) > 0) offset + addO
      else addO
    }
  }

  beta <- object$coefficients

  coef_ind <- setNames(lapply(seq_len(nlevels(y) - 1), function(i) {
    (i - 1) * ncol(x) + seq_len(ncol(x))
  }), levels(y)[-1])

  qq <- lapply(levels(y)[-1], function(i) {
    exp(offset + drop(x %*% beta[coef_ind[[i]]]))
  })

  pden <- 1 + rowSums(do.call("cbind", qq))

  pp <- do.call("cbind", c(list(1), qq)) / pden

  colnames(pp) <- levels(y)

  pp
}

#' @exportS3Method sandwich::estfun mlogit_weightit
estfun.mlogit_weightit <- function(x, ...) {
  x$gradient
}

#' @exportS3Method model.matrix mlogit_weightit
model.matrix.mlogit_weightit <- function(object, ...) {
  class(object) <- "lm"
  model.matrix(object, ...)
}