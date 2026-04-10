star_fit <- function(x, y, weights = NULL, start = NULL, offset = NULL,
                     transformation = "identity", hess = TRUE, control = list(),
                     bin = "round", ...) {
  chk::chk_atomic(y)
  chk::chk_numeric(x)
  chk::chk_matrix(x)

  n <- length(y)

  if (is_null(weights)) weights <- rep.int(1, n)
  else chk::chk_numeric(weights)

  if (is.null(offset)) offset <- rep.int(0, n)
  else chk::chk_numeric(offset)

  chk::chk_all_equal(c(length(y), nrow(x), length(weights), length(offset)))

  k0 <- ncol(x) + 1

  aliased_X <- !colnames(x) %in% colnames(make_full_rank(x, with.intercept = FALSE))
  aliased_B <- c(aliased_X, FALSE)

  x_ <- x[, !aliased_X, drop = FALSE]

  no_x <- is_null(x_)

  nm <- c(colnames(x), "log(SD)")

  start <- c(rep(0, ncol(x_)), log(sd(y)))

  names(start) <- nm[!aliased_B]

  lower_bin <- switch(bin,
                      "round" = function(z) z - .5,
                      "floor" = function(z) z,
                      "ceiling" = function(z) z - 1)
  upper_bin <- function(z) lower_bin(z) + 1

  trans <- switch(transformation,
                  "identity" = identity,
                  "log" = log,
                  "sqrt" = sqrt)

  transinv <- switch(transformation,
                     "identity" = identity,
                     "log" = exp,
                     "sqrt" = function(z) z^2)

  lli <- function(B, X, y, weights, offset) {

    b <- B[seq_col(X)]
    sd <- exp(B[length(B)])
    xb <- offset + drop(X %*% b)

    # eta <- transinv(xb)
    eta <- xb

    # Probability of observed outcome
    uy <- trans(upper_bin(y))
    ly <- trans(lower_bin(y))

    uz <- (uy - eta) / sd
    lz <- (ly - eta) / sd

    p <- pnorm(uz) - pnorm(lz)

    weights * log(p)
  }

  ll <- function(B, X, y, weights, offset) {
    sum(lli(B, X, y, weights, offset))
  }

  m_control <- list(fnscale = -1, #maximize likelihood; optim() minimizes by default
                    trace = 0,
                    maxit = 1e3,
                    reltol = 1e-12)

  control <- utils::modifyList(m_control, control)

  # Psi function and gradient using natural parameterization
  psi <- function(B, X, y, weights, offset = NULL) {
    if (is_null(offset)) offset <- rep.int(0, length(y))

    b <- B[seq_col(X)]
    sd <- exp(B[length(B)])
    xb <- offset + drop(X %*% b)

    # eta <- transinv(xb)
    eta <- xb

    uy <- trans(upper_bin(y))
    ly <- trans(lower_bin(y))

    uz <- (uy - eta) / sd
    lz <- (ly - eta) / sd

    # Probability of observed outcome
    p <- pnorm(uz) - pnorm(lz)

    gju <- dnorm(uz)
    gjl <- dnorm(lz)

    out <- cbind(X * (gjl - gju) / sd, (gjl * lz - gju * uz)) * (weights / p)

    -((dnorm(uz) * uz - dnorm(lz) * lz) * x)/(p)

    colnames(out) <- names(B)
    out
  }

  gr <- function(B, X, y, weights, offset = NULL) {
    colSums(psi(B, X, y, weights, offset))
  }

  out0 <- optim(par = start,
                ll,
                gr = gr,
                X = x_,
                y = y,
                weights = weights,
                offset = offset,
                method = "BFGS",
                hessian = FALSE,
                control = control)

  theta0 <- out0$par

  hessian <- NULL
  if (hess) {
    # Estimate using natural parameterization to get hessian
    hessian <- try(optimHess(par = theta0,
                             ll,
                             X = x_,
                             y = y,
                             weights = weights,
                             offset = offset,
                             gr = gr,
                             control = list(fnscale = -1)),
                   silent = TRUE)

    # If optimization fails, use numeric differentiation of gradient to get hessian
    if (null_or_error(hessian)) {
      hessian <- .gradient(gr, theta0,
                           X = x_,
                           y = y,
                           weights = weights,
                           offset = offset)
    }

    colnames(hessian) <- rownames(hessian) <- names(theta0)
  }

  theta <- theta0

  grad <- psi(theta, X = x_, y = y,
              weights = weights, offset = offset)

  xb <- offset + drop(x_ %*% theta[seq_col(x_)])
  sd <- exp(theta[length(theta)])
  k <- qnorm(1e-6)

  # Fitted values
  fitted <- vapply((xb), function(e) {
    y_range <- seq(floor(e - k * sd),
                   ceiling(e + k * sd))

    uy <- trans(upper_bin(y_range))
    ly <- trans(lower_bin(y_range))

    uz <- (uy - e) / sd
    lz <- (ly - e) / sd

    sum((pnorm(uz) - pnorm(lz)) * y_range)
  }, numeric(1L))

  # Residuals
  res <- y - fitted

  coefs <- setNames(rep.int(NA_real_, ncol(x) + 1L), nm)
  coefs[!aliased_B] <- theta

  list(coefficients = coefs,
       residuals = res,
       fitted.values = fitted,
       linear.predictors = xb,
       solve = out0,
       psi = psi,
       f = gr,
       df.residual = length(res) - ncol(x_) - 1L,
       x = x,
       y = y,
       weights = weights,
       gradient = grad,
       hessian = hessian,
       transformation = transformation,
       bin = bin)
}

star_lm <- function(formula, data, weights, subset, start = NULL, na.action,
                    hess = TRUE, control = list(), model = TRUE,
                    x = FALSE, y = TRUE, contrasts = NULL, transformation = "identity",
                    bin = "round", ...) {
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
  if (!is.null(weights)) {
    chk::chk_numeric(weights)
    chk::chk_gte(weights)
  }

  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    chk::chk_numeric(offset)
    if (length(offset) != NROW(Y))
      .err(gettextf("number of offsets is %d; should equal %d (number of observations)",
                    length(offset), NROW(Y)), domain = NA)
  }

  fit <- eval(call("star_fit",
                   x = X, y = Y, transformation = transformation, weights = weights,
                   offset = offset, start = start,
                   hess = hess, control = control, bin = bin))

  if (model) fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (!x) fit$x <- NULL
  if (!y) fit$y <- NULL

  out <- c(fit,
           list(call = cal, formula = formula, terms = mt,
                data = data, offset = offset,
                contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)))

  class(out) <- c("star_lm", "glm")
  out
}

predict.star_lm <- function(object, newdata = NULL, type = "mean",
                            na.action = na.pass, ...) {
  chk::chk_string(type)
  type <- switch(type, "probs" = "response", "lp" = "link", type)
  type <- match_arg(type, c("mean", "response", "link"))

  na.act <- object$na.action
  object$na.action <- NULL

  if (is_null(newdata)) {
    if (type == "link") {
      out <- object$linear.predictors
    }
    else if (type == "response") {
      bin <- match.fun(object$bin)
      out <- bin(object$linear.predictors)
    }
    else if (type == "mean") {
      out <- object$fitted.values
    }

    if (is_not_null(na.act)) {
      out <- napredict(na.act, out)
    }

    return(out)
  }

  tt <- terms(object)

  Terms <- delete.response(tt)
  m <- model.frame(Terms, newdata, na.action = na.action,
                   xlev = object$xlevels)

  cl <- attr(Terms, "dataClasses")
  if (is_not_null(cl)) {
    .checkMFClasses(cl, m)
  }
  x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)

  offset <- model.offset(m)
  addO <- object$call$offset

  if (is_not_null(addO)) {
    addO <- eval(addO, newdata, environment(tt))
    offset <- {
      if (is_null(offset)) addO
      else offset + addO
    }
  }
  if (is_null(offset)) offset <- rep.int(0, nrow(x))

  xb <- offset + drop(x %*% object$coefficients[seq_col(x)])

  if (type == "link") {
    return(xb)
  }

  if (type == "response") {
    bin <- match.fun(object$bin)
    return(setNames(bin(xb), rownames(x)))
  }

  sd <- exp(object$coefficients[length(object$coefficients)])

  trans <- switch(object$transformation,
                  "identity" = identity,
                  "log" = exp,
                  "sqrt" = function(z) z^2)

  eta <- xb

  lower_bin <- switch(object$bin,
                      "round" = function(z) z - .5,
                      "floor" = function(z) z,
                      "ceiling" = function(z) z - 1)
  upper_bin <- function(z) lower_bin(z) + 1

  y_range <- seq(floor(min(eta) - 4 * sd),
                 ceiling(max(eta) + 4 * sd))

  p <- do.call("cbind", lapply(y_range, function(y_) {
    uy <- trans(upper_bin(y_))
    ly <- trans(lower_bin(y_))

    uz <- (uy - eta) / sd
    lz <- (ly - eta) / sd

    pnorm(uz) - pnorm(lz)
  }))

  no_dens <- colSums(p) < 1e-8

  p <- p[,!no_dens, drop = FALSE]
  y_range <- y_range[!no_dens]

  setNames(drop(p %*% y_range), rownames(x))
}

vcov.star_lm <- function(object, complete = TRUE, ...) {
  .vcov.aliased(is.na(object$coefficients),
                -solve(object[["hessian"]]),
                complete = complete)
}

