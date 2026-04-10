#NPCBPS
f <- function(A, X) {
  X <- scale(X[,-1])
  X_ <- cbind(A * X, (1 - A) * X)
  n <- length(A)

  log0 <- function(x, e = 1e-5) {
    z <- log(e) + x/e - 1 - (x - e)^2/(2*e^2)
    z[x > e] <- log(x[x > e])
    z
  }

  obj <- function(gamma, .X, e = 1e-5) {
    XG <- drop(1 - .X %*% gamma)
    sum(-log0(XG, e))
  }

  opt <- optim(rep(c(-.1, .1), each = ncol(X)),
               obj,
               .X = X_,
               e = 1/n,
               method = "BFGS",
               hessian = TRUE,
               control = list(reltol = 1e-10,
                              maxit = 1e4))

  c(list(weights = drop(1/(1 - X_ %*% opt$par))),
    opt)
}

data("lalonde")
A <- lalonde$treat
X <- model.matrix(~age + educ + married + race + re74 + re75, data = lalonde)

res <- f(A, X)

bal.tab(X, treat = A, weights = res$weights)

weightit2npcbps2 <- function(treat, covs, ...) {
  n <- length(treat)

  f <- function(A, X, e = 1/n, ...) {
    X <- scale(X[,-1])
    X_ <- cbind(A * X, (1 - A) * X)

    log0 <- function(x, e = 1e-5) {
      z <- log(e) + x/e - 1 - (x - e)^2/(2*e^2)
      z[x > e] <- log(x[x > e])
      z
    }

    obj <- function(gamma, .X, e = 1e-5) {
      XG <- drop(1 - .X %*% gamma)
      sum(-log0(XG, e))
    }

    opt <- optim(rep(c(-.1, .1), each = ncol(X)),
                 obj,
                 .X = X_,
                 e = e,
                 method = "BFGS",
                 hessian = TRUE,
                 control = list(reltol = 1e-10,
                                maxit = 1e4))

    c(list(weights = drop(1/(1 - X_ %*% opt$par))),
      opt)
  }

  out <- f(A, X, ...)

  list(w = out$weights)
}