#Cham et al. (2016)

#p(v); p'(v); D(v)

ebcw <- function(A, X, p, pp) {
  X[,-1] <- scale(X[,-1])

  p_ <- function(v, a, e = 1/length(A)) {
    suppressWarnings(o <- p(v, a))
    if (any(bad <- !is.finite(o))) {
      o[bad] <- p(e, a) + pp(e, a) * (v[bad] - e)
    }

    o
  }

  f <- function(beta, a, X) {
    eta <- drop(X %*% beta)
    mean(a * p_(eta, a) - eta)
  }

  opt1 <- optim(rep(.4, ncol(X)),
                f,
                a = A,
                X = X,
                method = "BFGS",
                control = list(fnscale = -1,
                               maxit = 1e3,
                               reltol = 1e-12))

  opt0 <- optim(rep(.4, ncol(X)),
                f,
                a = 1 - A,
                X = X,
                method = "BFGS",
                control = list(fnscale = -1,
                               maxit = 1e3,
                               reltol = 1e-12))

  weights <- numeric(length(A))

  weights[A == 1] <- pp(X[A==1,,drop = FALSE] %*% opt1$par, A)
  weights[A == 0] <- pp(X[A==0,,drop = FALSE] %*% opt0$par, 1 - A)

  list(weights = weights,
       opt1 = opt1, opt0 = opt0)
}

p <- function(v, a) {
  log(1/v) + v - 1
}

pp <- function(v, a) {
  1 - 1/v
}

A <- lalonde$treat
X <- model.matrix(~1 + age + educ + married + race + re74 + re75, data = lalonde)

res <- ebcw(A, X, p, pp)

bal.tab(scale(X[,-1]), treat = A, weights = list(ebcw = res$weights, weightit = W$weights))

W <- weightit(treat ~ age + educ + married + race + re74 + re75, data = lalonde,
              method = "ebal", link = "probit")
W <- lmw::lmw(~ age + educ + married + race + re74 + re75, data = lalonde,
              method = "MRI", estimand = "ATE", treat = "treat")

#Ebal
p <- function(v, a) {
  -exp(-v)
}

pp <- function(v, a) {
  exp(-v)
}

D <- function(v) {
  sum(x * log(x))
}

#Logistic
p <- function(v, a) {
  v - exp(-v)
}

pp <- function(v, a) {
  1 + exp(-v)
}

#SBW
p <- function(v, a) {
  -(v^2)/4
}

pp <- function(v, a) {
  -v/2
}