#M-estimation when only the likelihood is known
#What I learned: it is possible to estimate the M-estimation variance,
#but to estimate the gradient and hessian requires a numerical approximation
#that is liable to be inaccurate. Using optimHess() seems to be
#the most accurate for calculating the hessian.
lli <- function(b, y, x, w = NULL) {
  if (is.null(w)) w <- 1

  p <- plogis(drop(x %*% b))

  w * dbinom(y, 1, p, log = TRUE)
}

ll <- function(b, y, x, w = NULL) {
  sum(lli(b, y, x, w))
}

getw <- function(b, x, a) {
  p <- plogis(drop(x %*% b))

  w <- rep(1, length(p))
  w[a == 0] <- p[a == 0] / (1 - p[a == 0])

  w
}

X <- cbind(1, lalonde$age, lalonde$educ)
A <- lalonde$treat
Y <- as.integer(lalonde$re78 > 0)

#PS model for weights
opt.ps <- optim(rep(0, ncol(X)),
                ll,
                y = A,
                x = X,
                method = "BFGS",
                control = list(fnscale = -1, reltol = 1e-12))
b.ps <- opt.ps$par

psw <- getw(b.ps, X, A)

#Outcome model for weights
opt.out <- optim(rep(0, ncol(X) + 1),
                 ll,
                 y = Y,
                 x = cbind(X, A),
                 w = psw,
                 method = "BFGS",
                 control = list(fnscale = -1, maxit = 1e3, reltol = 1e-16))
b.out <- opt.out$par
b.out <- fit$coefficients
psi_b <- .gradient(function(b, a, xw, y, xy) {
  lli(b[-(1:ncol(xw))], y = Y, x = xy, w = getw(b[1:ncol(xw)], xw, a))
},
c(b.ps, b.out), a = A, xw = X, y = Y, xy = cbind(X, A))

hess <- .gradient(function(b, a, xw, y, xy) {
    .gradient(function(.b, a, xw, y, xy) {
      ll(.b[-(1:ncol(xw))], y = Y, x = xy, w = getw(.b[1:ncol(xw)], xw, a))
    }, b, a = a, xw = xw, y = y, xy = xy)
},
c(b.ps, b.out), a = A, xw = X, y = Y, xy = cbind(X, A))

hess <- rootSolve::hessian(function(b, a, xw, y, xy) {
  ll(b[-(1:ncol(xw))], y = Y, x = xy, w = getw(b[1:ncol(xw)], xw, a))
},
c(b.ps, b.out), a = A, xw = X, y = Y, xy = cbind(X, A),
centered = T)


hess <- optimHess(c(b.ps, b.out),
                  function(b, a, xw, y, xy) {
                    ll(b[-(1:ncol(xw))], y = Y, x = xy, w = getw(b[1:ncol(xw)], xw, a))
                  },
                  a = A, xw = X, y = Y, xy = cbind(X, A),
                  control = list(fnscale = -1))

A1 <- solve(hess)
B <- crossprod(psi_b)

V <- A1 %*% tcrossprod(B, A1)
V[-c(1:3), -c(1:3)] |> diag() |> sqrt()

W <- weightit(treat ~ age + educ, data = lalonde, estimand = "ATT")
fit <- glm_weightit(re78 > 0 ~ age + educ + treat, data = lalonde,
                    weightit = W, family = binomial, vcov = "HC0")

vcov(fit) |> diag() |> sqrt()
coef(fit)
