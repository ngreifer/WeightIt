#NPCBPS

f <- function(A, X) {

  # targets <- colMeans(X)

  X <- scale(X)

  X_ <- cbind(A * X, (1 - A) * X)

  n <- length(A)

  obj <- function(gamma, .X, e) {
    eta <- drop(1 - .X %*% gamma)

    z <- log(e) + (eta - e)/e - (eta - e)^2/(2 * e^2)
    z[eta > e] <- log(eta[eta > e])

    -sum(z)
  }

  opt <- optim(rep(c(-.1, .1), each = ncol(X)),
               obj,
               .X = X_,
               e = 1/n,
               method = "BFGS",
               control = list(reltol = 1e-10,
                              maxit = 1e4))

  c(list(weights = drop(1/(1 - X_ %*% opt$par))),
    opt)
}

data("lalonde")
A <- lalonde$treat
X <- model.matrix(~0 + age + educ + married + race + re74 + re75, data = lalonde)

res <- f(A, X)

bal.tab(X, treat = A, weights = res$weights)

f2 <- function(A, X) {

  X <- scale(X[,-1])

  X_ <- cbind(A * X, (1 - A) * X)

  n <- length(A)

  log0 <- function(x, e = 1e-5) {
    z <- log(e) + x/e - 1 - (x - e)^2/(2*e^2)
    z[x > e] <- log(x[x > e])
    z
  }

  obj <- function(gamma, .X) {
    XG <- drop(1 - .X %*% gamma)

    sum(-log0(XG))
  }

  opt <- optim(rep(c(-.1, .1), each = ncol(X)),
               obj,
               .X = X_,
               method = "BFGS",
               control = list(reltol = 1e-10,
                              maxit = 1e4))

  c(list(weights = drop(1/(1 - X_ %*% opt$par))),
    opt)
}

data("lalonde")
A <- lalonde$treat
X <- model.matrix(~1 + age + educ + married + race + re74 + re75, data = lalonde)

res <- f2(A, X)

bal.tab(X, treat = A, weights = res$weights)
