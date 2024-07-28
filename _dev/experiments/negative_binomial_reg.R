#Negative binomial regression using MLE and M-estimation
#Could go into glm_weightit()

a <- 1

x <- rnorm(1e4)
y <- rnbinom(1e4, size = a, mu = exp(.5 + .5 * x))

fit <- MASS::glm.nb(y ~ x)
c((coef(fit)), fit$theta)
c(sqrt(diag(vcov(fit))), fit$SE.theta)

ll <- function(b, x, y) {
  lp <- drop(x %*% b[-length(b)])
  sum(dnbinom(y, exp(b[length(b)]), mu = exp(lp), log = TRUE))
}

opt <- optim(ll, par = c(0, 0, 1),
             x = cbind(1,x), y = y, hessian = TRUE,
             method = "BFGS",
             control = list(fnscale = -1, reltol = 1e-14))

c(opt$par[1:2], exp(opt$par[3]))
sqrt(diag(solve(-opt$hessian)))



psi <- function(b, x, y) {
  lp <- drop(x %*% b[-length(b)])
  th <- b[length(b)]

  .fam <- MASS::negative.binomial(th)
  mu <- .fam$linkinv(lp)

  psi_b <-  x * (.fam$mu.eta(lp) * (y - mu) / .fam$variance(mu))

  # psi_b <- (y - (y + th) * mu / (th + mu)) * x
  psi_th <- (digamma(th + y) - digamma(th) + log(th) + 1 - log(th + mu) - (y + th)/(mu + th))

  cbind(psi_b, psi_th)
}

grad <- function(b, x, y) {
  colSums(psi(b, x, y))
}

r <- rootSolve::multiroot(grad, c(0, 0, .5),
                          x = cbind(1, x), y = y)
r$root
