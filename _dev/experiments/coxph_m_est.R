#Cox PH model using M-estimation
#Can estimate parameters, but can't get sandwich vcov, so not going anywhere

#Generating data similar to Austin (2009) for demonstrating treatment effect estimation
gen_X <- function(n) {
  X <- matrix(rnorm(9 * n), nrow = n, ncol = 9)
  X[,5] <- as.numeric(X[,5] < .5)
  X
}

gen_Ac <- function(X) {
  LP_A <- -1.2 + log(2)*X[,1] - log(1.5)*X[,2] + log(2)*X[,4] - log(2.4)*X[,5] + log(2)*X[,7] - log(1.5)*X[,8]
  LP_A + rlogis(nrow(X))
}

#~20% treated
gen_A <- function(Ac) {
  1 * (Ac > 0)
}

gen_Am <- function(A) {
  factor(ifelse(A == 1, "T", sample(c("C1", "C2"), length(A), TRUE)))
}

# Continuous outcome
gen_Y_C <- function(A, X) {
  2*A + 2*X[,1] + 2*X[,2] + 2*X[,3] + 1*X[,4] + 2*X[,5] + 1*X[,6] + rnorm(length(A), 0, 5)
}
#Conditional:
#  MD: 2
#Marginal:
#  MD: 2

# Binary outcome
gen_Y_B <- function(A, X) {
  LP_B <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
  P_B <- plogis(LP_B)
  rbinom(length(A), 1, P_B)
}
#Conditional:
#  OR:   2.4
#  logOR: .875
#Marginal:
#  RD:    .144
#  RR:   1.54
#  logRR: .433
#  OR:   1.92
#  logOR  .655

# Survival outcome
gen_Y_S <- function(A, X) {
  LP_S <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
  sqrt(-log(runif(length(A)))*2e4*exp(-LP_S))
}
#Conditional:
#  HR:   2.4
#  logHR: .875
#Marginal:
#  HR:   1.57
#  logHR: .452

n <- 200
X <- gen_X(n)
Ac <- gen_Ac(X)
A <- gen_A(Ac)
Am <- gen_Am(A)

Y_C <- gen_Y_C(A, X)
Y_B <- gen_Y_B(A, X)
Y_S <- gen_Y_S(A, X)

d <- data.frame(A, Am, Ac, X, Y_C, Y_B, Y_S)
d$sw <- 1#runif(nrow(d))

library(survival)
fit <- coxph(Surv(Y_S) ~ X1 + X2, data = d, weights = d$sw, robust = T)

X <- model.matrix(fit)
Y <- d$Y_S
ord <- order(Y, decreasing = T)

psi <- function(B, X, SW) {
  psi_ <- 0 * X
  X <- X[ord,, drop = FALSE]
  SW <- SW[ord]

  swexlp <- SW * exp(drop(X %*% B))

  out <- SW * (X - apply(X, 2, function(x) cumsum(x * swexlp)) /
                 cumsum(swexlp))
  psi_[ord,] <- out
  psi_
}

gradfun <- function(B, X, SW) {
  colSums(psi(B, X, SW))
}

out <- rootSolve::multiroot(gradfun,
                            rep(0, ncol(X)),
                            X = X,
                            SW = d$sw)
out$root
unname(fit$coefficients)


psi_b <- psi(out$root, X = X, SW = d$sw)
B <- crossprod(psi_b)

hess <- gradient(gradfun,
                 .x = out$root,
                 X = X,
                 SW = d$sw)

A1 <- solve(-hess)
V <- A1 %*% tcrossprod(B, A1)

V
fit$var
crossprod(residuals(fit, type = "dfbeta",
                    weighted = TRUE))

A1
fit$naive.var

head(psi_b)
head(residuals(fit, type = "score"))

#Lin & Wei 1989
getW <- function(B, X, SW) {
  X <- X[ord,, drop = FALSE]
  SW <- SW[ord]

  swexlp <- SW * exp(drop(X %*% B))

  S1 <- apply(X, 2, function(x) cumsum(x * swexlp))
  S0 <- cumsum(swexlp)

  (X - S1/S0) - cumsum(1 - )

  out <- SW * (X - apply(X, 2, function(x) cumsum(x * swexlp)) /
                 cumsum(swexlp))

  (X - apply(X, 2, function(x) cumsum(x * swexlp)) /
      cumsum(swexlp))
}