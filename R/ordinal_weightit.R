# Ordinal regression
.ordinal_weightit.fit <- function(x, y, weights = NULL, start = NULL, offset = NULL,
                                  link = "logit", hess = TRUE, ...) {
  chk::chk_atomic(y)
  chk::chk_numeric(x)
  chk::chk_matrix(x)

  chk::chk_string(link)

  family <- binomial(link)

  if (!is.function(family$linkinv)) {
    .err("the supplied link seems not to create a valid binomial family object")
  }

  .linkfun <- family$linkfun
  .linkinv <- family$linkinv
  .mu.eta <- family$mu.eta

  y <- droplevels(as.factor(y))
  n <- length(y)

  if (is_null(weights)) weights <- rep.int(1, n)
  else chk::chk_numeric(weights)

  if (is.null(offset)) offset <- rep.int(0, n)
  else chk::chk_numeric(offset)

  chk::chk_all_equal(c(length(y), nrow(x), length(weights), length(offset)))

  x <- x[,colnames(x) != "(Intercept)", drop = FALSE]

  m <- nlevels(y) #num. thresholds
  k0 <- ncol(x) + m - 1 #num. params

  nm <- c(colnames(x), paste(levels(y)[-m], levels(y)[-1], sep = "|"))

  aliased_X <- !colnames(x) %in% colnames(make_full_rank(x, with.intercept = TRUE))
  aliased_B <- c(aliased_X, rep.int(FALSE, m - 1L))

  x_ <- x[, !aliased_X, drop = FALSE]
  y_ <- as.integer(y)

  no_x <- is_null(x_)

  if (no_x) {
    start <- .linkfun(cumsum(tabulate(y_)[-m]/n))
  }
  else if (is_null(start)) {
    q1 <- floor(median(y_))
    y1 <- as.numeric(y_ > q1)
    X <- cbind(1, x_)
    fit <- suppressWarnings(glm.fit(X, y1, weights, family = family, offset = offset))

    coefs <- {
      if (!fit$converged) c(.linkfun(weighted.mean(y1, weights)), rep.int(0, ncol(x_)))
      else fit$coefficients
    }

    if (m > 2L) {
      spacing <- .linkfun(cumsum(tabulate(y_)[-m]/n))
      start <- c(coefs[-1L], -coefs[1L] + spacing - spacing[q1])
    }
    else {
      start <- c(coefs[-1L], -coefs[1L])
    }
  }
  else {
    chk::chk_numeric(start)
    chk::chk_length(start, k0)

    start <- start[!aliased_B]

    if (any(diff(start[-seq_col(x)]) <= 0)) {
      .err("starting values for the thresholds must be in ascending order")
    }
  }

  # Adjust start to use cumsum parameterization
  if (m > 2) {
    if (no_x) {
      start[-1L] <- log(diff(start))
    }
    else {
      start[-seq_len(ncol(x_) + 1L)] <- log(diff(start[-seq_len(ncol(x_))]))
    }
  }

  names(start) <- nm[!aliased_B]

  ind_mat <- cbind(seq_along(y), y_)

  y_mat <- 1 * vapply(seq_len(m), function(j) y_ == j, logical(length(y_)))
  y_mat0 <- y_mat[, -m, drop = FALSE]
  y_mat1 <- y_mat[, -1L, drop = FALSE]

  #Get predictors on a smaller scale
  if (!no_x) {
    sds <- apply(x_, 2L, sd)
    x_ <- sweep(x_, 2L, sds, "/")
    start <- start * c(sds, rep.int(1, m - 1L))
  }

  # Get predicted probabilities for all units for category y
  get_p <- function(y, xb, a) {
    pmax(.linkinv(c(a, Inf)[y] - xb) - .linkinv(c(-Inf, a)[y] - xb),
         1e-16)
  }

  #Ordinal regression LL (cumsum paramaterization)
  ll <- function(B, X, y, weights, offset, .cumsum_param = TRUE) {
    if (no_x) {
      a <- B
      Xb <- offset
    }
    else {
      a <- B[-seq_col(X)]
      b <- B[seq_col(X)]
      Xb <- offset + drop(X %*% b)
    }

    if (.cumsum_param && m > 2L) {
      a <- c(a[1L], a[1L] + cumsum(exp(a[-1L])))
    }

    # Probability of observed outcome
    p <- get_p(y, Xb, a)

    sum(weights * log(p))
  }

  dots <- list(...)

  control <- list(fnscale = -1, #maximize likelihood; optim() minimizes by default
                  trace = 0,
                  maxit = 1e3,
                  reltol = 1e-12)

  control <- utils::modifyList(control,
                               dots[intersect(names(dots), c("trace", "maxit", "reltol", "ndeps", "REPORT"))])

  # Estimate using cumsum parameterization to get estimates
  out0 <- optim(par = start,
                ll,
                X = x_,
                y = y_,
                weights = weights,
                offset = offset,
                .cumsum_param = TRUE,
                method = "BFGS",
                hessian = FALSE,
                control = control)

  theta0 <- out0$par

  # Convert to natural param
  if (m > 2L) {
    if (no_x) {
      theta0[-1L] <- theta0[1L] + cumsum(exp(theta0[-1L]))
    }
    else {
      a1 <- theta0[ncol(x_) + 1L]
      theta0[-seq_len(ncol(x_))] <- c(a1, a1 + cumsum(exp(theta0[-seq_len(ncol(x_) + 1L)])))
    }
  }

  # Psi function and gradient using natural parameterization
  psi <- function(B, X, y, weights, offset = NULL) {
    if (is_null(offset)) offset <- rep.int(0, length(y))

    if (no_x) {
      a <- B
      Xb <- offset
    }
    else {
      a <- B[-seq_col(X)]
      b <- B[seq_col(X)]
      Xb <- offset + drop(X %*% b)
    }

    pp <- get_p(y, Xb, a)

    gj <- .mu.eta(c(a, Inf)[y] - Xb)
    gj1 <- .mu.eta(c(-Inf, a)[y] - Xb)

    .psi_a <- gj * y_mat0 - gj1 * y_mat1

    if (no_x) {
      out <- as.matrix(.psi_a) * (weights / pp)
    }
    else {
      .psi_b <- X * (gj1 - gj)
      out <- cbind(.psi_b, .psi_a) * (weights / pp)
    }

    colnames(out) <- names(B)
    out
  }

  gr <- function(B, X, y, weights, offset = NULL) {
    colSums(psi(B, X, y, weights, offset))
  }

  hessian <- NULL
  if (hess) {
    # Estimate using natural parameterization to get hessian
    hessian <- try(optimHess(par = theta0,
                             function(...) ll(..., .cumsum_param = FALSE),
                             X = x_,
                             y = y_,
                             weights = weights,
                             offset = offset,
                             gr = gr,
                             control = list(fnscale = -1)),
                   silent = TRUE)

    # If optimization fails, use numeric differentiation of gradient to get hessian
    if (null_or_error(hessian)) {
      hessian <- .gradient(gr, theta0,
                           X = x_,
                           y = y_,
                           weights = weights,
                           offset = offset)
    }

    if (!no_x) {
      hessian <- hessian * tcrossprod(c(sds, rep.int(1, m - 1L)))
    }

    colnames(hessian) <- rownames(hessian) <- names(theta0)
  }

  theta <- theta0

  grad <- psi(theta, X = x_, y = y_,
              weights = weights, offset = offset)

  # Get predicted probabilities for all units for all categories,
  # natural parameterization of `a`
  get_pp <- function(B, X, offset = NULL) {
    if (length(offset) == 0L) offset <- rep.int(0, n)

    if (ncol(X) == 0L) {
      a <- B
      Xb <- offset
    }
    else {
      a <- B[-seq_col(X)]
      b <- B[seq_col(X)]
      Xb <- offset + drop(X %*% b)
    }

    GG <- vapply(a, function(a_) .linkinv(a_ - Xb),
                 numeric(length(Xb)))

    pp <- cbind(GG, 1) - cbind(0, GG)

    dimnames(pp) <- list(rownames(X), levels(y))

    pp
  }

  # Fitted values
  pp <- get_pp(theta, x_, offset)

  # Residuals
  res <- setNames(1 - pp[ind_mat], rownames(x))

  # Adjust estimates and gradient to be put on original scale
  if (!no_x) {
    theta <- theta / c(sds, rep.int(1, m - 1L))
    grad <- sweep(grad, 2, c(sds, rep.int(1, m - 1L)), "*")
  }

  coefs <- setNames(rep.int(NA_real_, ncol(x) + m - 1L), nm)
  coefs[!aliased_B] <- theta

  list(coefficients = coefs,
       residuals = res,
       fitted.values = pp,
       family = family,
       linear.predictors = offset + drop(x_ %*% theta[seq_col(x_)]),
       solve = out0,
       psi = psi,
       f = gr,
       get_p = get_pp,
       df.residual = length(res) - ncol(x_) - (m - 1),
       x = x,
       y = y,
       weights = weights,
       gradient = grad,
       hessian = hessian)
}

.ordinal_weightit <- function(formula, data, link = "logit", weights, subset, start = NULL, na.action,
                              hess = TRUE, model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...) {
  cal <- match.call()

  chk::chk_flag(hess)
  chk::chk_flag(model)
  chk::chk_flag(x)
  chk::chk_flag(y)
  chk::chk_string(link)

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

  fit <- eval(call(".ordinal_weightit.fit",
                   x = X, y = Y, link = link, weights = weights,
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

.get_hess_ordinal <- function(fit) {
  x <- if_null_then(fit[["x"]], model.matrix(fit))
  y <- if_null_then(fit[["y"]], model.response(model.frame(fit)))
  family <- fit[["family"]]
  weights <- weights(fit)
  offset <- fit$offset
  coefs <- coef(fit)

  .linkinv <- family$linkinv
  .mu.eta <- family$mu.eta

  y <- droplevels(as.factor(y))
  n <- length(y)

  if (is_null(weights)) weights <- rep.int(1, n)

  if (is.null(offset)) offset <- rep.int(0, n)

  m <- nlevels(y) #num. thresholds

  aliased_X <- !colnames(x) %in% colnames(make_full_rank(x, with.intercept = TRUE))

  x_ <- x[, !aliased_X, drop = FALSE]
  y_ <- as.integer(y)

  no_x <- is_null(x_)

  y_mat <- 1 * vapply(seq_len(m), function(j) y_ == j, logical(length(y_)))
  y_mat0 <- y_mat[, -m, drop = FALSE]
  y_mat1 <- y_mat[, -1L, drop = FALSE]

  #Get predictors on a smaller scale
  if (!no_x) {
    sds <- apply(x_, 2L, sd)
    x_ <- sweep(x_, 2L, sds, "/")
  }

  # Get predicted probabilities for all units for category y
  get_p <- function(y, xb, a) {
    pmax(.linkinv(c(a, Inf)[y] - xb) - .linkinv(c(-Inf, a)[y] - xb),
         1e-16)
  }

  #Ordinal regression LL (cumsum paramaterization)
  ll <- function(B, X, y, weights, offset, .cumsum_param = TRUE) {
    if (no_x) {
      a <- B
      Xb <- offset
    }
    else {
      a <- B[-seq_col(X)]
      b <- B[seq_col(X)]
      Xb <- offset + drop(X %*% b)
    }

    if (.cumsum_param && m > 2L) {
      a <- c(a[1L], a[1L] + cumsum(exp(a[-1L])))
    }

    # Probability of observed outcome
    p <- get_p(y, Xb, a)

    sum(weights * log(p))
  }

  theta0 <- na.rem(coefs) * c(sds, rep.int(1, m - 1L))

  # Psi function and gradient using natural parameterization
  psi <- function(B, X, y, weights, offset = NULL) {
    if (is_null(offset)) offset <- rep.int(0, length(y))

    if (no_x) {
      a <- B
      Xb <- offset
    }
    else {
      a <- B[-seq_col(X)]
      b <- B[seq_col(X)]
      Xb <- offset + drop(X %*% b)
    }

    pp <- get_p(y, Xb, a)

    gj <- .mu.eta(c(a, Inf)[y] - Xb)
    gj1 <- .mu.eta(c(-Inf, a)[y] - Xb)

    .psi_a <- gj * y_mat0 - gj1 * y_mat1

    if (no_x) {
      out <- as.matrix(.psi_a) * (weights / pp)
    }
    else {
      .psi_b <- X * (gj1 - gj)
      out <- cbind(.psi_b, .psi_a) * (weights / pp)
    }

    colnames(out) <- names(B)
    out
  }

  gr <- function(B, X, y, weights, offset = NULL) {
    colSums(psi(B, X, y, weights, offset))
  }

  # Estimate using natural parameterization to get hessian
  hessian <- try(optimHess(par = theta0,
                           function(...) ll(..., .cumsum_param = FALSE),
                           X = x_,
                           y = y_,
                           weights = weights,
                           offset = offset,
                           gr = gr,
                           control = list(fnscale = -1)),
                 silent = TRUE)

  # If optimization fails, use numeric differentiation of gradient to get hessian
  if (null_or_error(hessian)) {
    hessian <- .gradient(gr, theta0,
                         X = x_,
                         y = y_,
                         weights = weights,
                         offset = offset)
  }

  if (!no_x) {
    hessian <- hessian * tcrossprod(c(sds, rep.int(1, m - 1L)))
  }

  colnames(hessian) <- rownames(hessian) <- names(theta0)

  hessian
}