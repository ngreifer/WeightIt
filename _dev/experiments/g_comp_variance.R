#Comparing M-estimation, bootstrap, and Hansen & Overgaard (2024) for g-comutation variance
library(sandwich)
n <- 200

p_A <- .5
m_B <- 0
m_C <- 0
g_C <- 3
alpha <- c(-1, -5, 2)

B <- rnorm(n, m_B)
C <- rnorm(n, m_C)

p_A <- plogis(g_C * C)

A <- rbinom(n, 1, p_A)

Y <- alpha[1] * A + alpha[2] * A * B + alpha[3] * B * C + rnorm(n)

d <- data.frame(A, B, C, Y)

library(marginaleffects)
#Standard g-computation
fit <- lm(Y ~ A * B + C, data = d)
fit <- WeightIt::lm_weightit(Y ~ A * B + C, data = d)
# est <- avg_comparisons(fit, variables = list(A = 0:1))

# se1 <- est$std.error
#Robust g-computation
avg_predictions(fit, variables = list(A = 0:1),
                vcov = "HC0") |> vcov()

#Bootstrap
bootfun <- function(i) {
  boot_data <- d[sample.int(n, replace = TRUE),]

  fit_boot <- update(fit, data = boot_data)

  p <- predict(fit_boot, newdata = rbind(transform(boot_data, A = 0),
                                    transform(boot_data, A = 1)))

  c(mean(p[1:n]),  mean(p[-(1:n)]))
}

cov(do.call("rbind", purrr::map(1:1000, bootfun)))

#M-estimation
est <- avg_predictions(fit, variables = list(A = 0:1))

psi <- function(beta, X, Y, X0, X1) {
  yhat <- drop(X %*% beta[1:ncol(X)])

  p0 <- drop(X0 %*% beta[1:ncol(X)])
  p1 <- drop(X1 %*% beta[1:ncol(X)])
  EP0 <- p0
  EP1 <- p1

  cbind(X * (Y - yhat),
        EP0 - beta[ncol(X) + 1],
        EP1 - beta[ncol(X) + 2])
}

X <- model.matrix(fit$terms, data = d)
X0 <- model.matrix(fit$terms, data = transform(d, A = 0))
X1 <- model.matrix(fit$terms, data = transform(d, A = 1))

beta <- c(coef(fit), est$estimate)

psi_b <- psi(beta, X, d$Y, X0, X1)
B_ <- crossprod(psi_b)
M <- WeightIt:::.gradient(function(b) colSums(psi(b, X, d$Y, X0, X1)),
                          .x = beta)

A1 <- WeightIt:::.solve_hessian(M)

V <- A1 %*% tcrossprod(B_, A1)

V[ncol(X) + 1:2,ncol(X) + 1:2]

r <- rootSolve::multiroot(function(b) colSums(psi(b, X, d$Y, X0, X1)),
                          rep(0, ncol(X) + 2))

#using formula in paper
beta_dot <- -solve(M[1:ncol(X), 1:ncol(X)], t(psi_b[,1:ncol(X)])) * n
p0 <- drop(X0 %*% beta[1:ncol(X)])
p1 <- drop(X1 %*% beta[1:ncol(X)])
dp1dbeta <- WeightIt:::.gradient(function(b) drop(X1 %*% b),
                                 .x = beta[1:ncol(X)])
dp0dbeta <- WeightIt:::.gradient(function(b) drop(X0 %*% b),
                                 .x = beta[1:ncol(X)])

#Equal to M-estimation
mean((p0 - mean(p0) + drop(colMeans(dp0dbeta) %*% beta_dot))^2) / n
mean((p1 - mean(p1) + drop(colMeans(dp1dbeta) %*% beta_dot))^2) / n

#Equal to sandwich + marginaleffects
# sqrt(mean((drop(S1 %*% beta_dot))^2) / n)

S <- do.call("cbind", lapply(0:1, function(t) {
  Xt <- model.matrix(fit$terms, data = transform(d, A = t))
  pt <- drop(Xt %*% beta[1:ncol(X)])
  dptdbeta <- WeightIt:::.gradient(function(b) drop(Xt %*% b), .x = beta[1:ncol(X)])

  (pt - mean(pt) + drop((colMeans(dptdbeta) %*% beta_dot))) / n
}))

crossprod(S)
sqrt(c(-1, 1) %*% (crossprod(S)) %*% c(-1, 1))

br <- bread(fit) #-solve(M)
me <- meat(fit)
Vf <- br %*% me %*% br / n #sandwich(fit), tcrossprod(beta_dot)

gr <- WeightIt:::.gradient(function(b) mean(get_predict(set_coef(fit, b), newdata = transform(d, A = 1), type = "response")$est),
                           .x = beta[1:ncol(X)])
#gr == colMeans(dp1dbeta)
gr %*% Vf %*% t(gr)

mean((drop(colMeans(dp1dbeta) %*% beta_dot))^2) / n

#psi_b[,1:ncol(X)] == estfun(fit)
#beta_dot = tcrossprod(bread(fit), estfun(fit))

mean((drop(gr %*% tcrossprod(bread(fit), estfun(fit))))^2) / n
mean((p1 - mean(p1) + drop(gr %*% tcrossprod(bread(fit), estfun(fit))))^2) / n

newdata <- datagrid(model = fit, A = 0:1, grid_type = "counterfactual")
variable <- "A"
byfun <- mean

p <- predictions(fit, newdata = newdata, numderiv=list("richardson"))

destdbeta <- attr(p, "jacobian")

beta_dot <- tcrossprod(sandwich::bread(fit), sandwich::estfun(fit))

S1 <- do.call("cbind", lapply(unique(newdata[[variable]]), function(v) {
  n <- ncol(beta_dot)
  in_v <- newdata[[variable]] == v
  (1 / n) * ((p$estimate[in_v] - byfun(p$estimate[in_v])) +
               drop((apply(destdbeta[in_v,,drop = FALSE], 2, byfun) %*% beta_dot)))
}))

V <- crossprod(S1)
V

p <- avg_predictions(fit, variables = variable, vcov  = "HC0", numderiv=list("richardson"))
destdbeta <- attr(p, "jacobian")
S1 <- do.call("cbind", lapply(1:2, function(v) {
  n <- ncol(beta_dot)
  # in_v <- newdata[[variable]] == v
  (1 / n) * ((p$estimate[v] - byfun(p$estimate[v])) +
               drop((destdbeta[v,] %*% beta_dot)))
}))

crossprod(S)

p0 <- p$estimate[newdata[[variable]] == 0]
dp0dbeta <- destdbeta[newdata[[variable]] == 0,]


sum((p0 - mean(p0))^2 / n^2) + sum((dp0dbeta %*% beta_dot)^2/n^2)

sum((((p0 - mean(p0)) + colMeans(dp0dbeta) %*% beta_dot) / n)^2)
mean((p0 - mean(p0))^2)/n +
  mean(p0 * colMeans(dp0dbeta) %*% beta_dot) * 2 / n +
  mean((colMeans(dp0dbeta) %*% beta_dot)^2) / n

#ATT
#Note: for M-estimation, need to include P(A=1) as estimated parameter
est <- avg_predictions(fit, variables = list(A = 0:1), newdata = subset(A == 1))

psi <- function(beta, X, Y, X0, X1) {
  yhat <- drop(X %*% beta[1:ncol(X)])

  p0 <- drop(X0 %*% beta[1:ncol(X)])
  p1 <- drop(X1 %*% beta[1:ncol(X)])

  PA1 <- as.numeric(A == 1)

  EP0 <- p0 * PA1 / beta[ncol(X) + 3]
  EP1 <- p1 * PA1 / beta[ncol(X) + 3]

  cbind(X * (Y - yhat),
        EP0 - beta[ncol(X) + 1],
        EP1 - beta[ncol(X) + 2],
        PA1 - beta[ncol(X) + 3])
}

X <- model.matrix(fit$terms, data = d)
X0 <- model.matrix(fit$terms, data = transform(d, A = 0))
X1 <- model.matrix(fit$terms, data = transform(d, A = 1))

beta <- c(coef(fit), est$estimate, mean(A == 1))

psi_b <- psi(beta, X, d$Y, X0, X1)
B_ <- crossprod(psi_b)
M <- WeightIt:::.gradient(function(b) colSums(psi(b, X, d$Y, X0, X1)),
                          .x = beta)
A1 <- solve(M, t(WeightIt:::.chol2(B_)))

tcrossprod(A1)[-(1:ncol(X)), -(1:ncol(X))]

r <- rootSolve::multiroot(function(b) colSums(psi(b, X, d$Y, X0, X1)),
                          rep(0, ncol(X) + 2))

beta_dot <- -solve(M[1:ncol(X), 1:ncol(X)], t(psi_b[,1:ncol(X)])) * n

S <- do.call("cbind", lapply(0:1, function(t) {
  Xt <- model.matrix(fit$terms, data = transform(d, A = t))
  pt <- drop(Xt %*% beta[1:ncol(X)])
  dptdbeta <- WeightIt:::.gradient(function(b) drop(Xt %*% b), .x = beta[1:ncol(X)])

  ((pt - mean(pt[A == 1])) * (A == 1) / (mean(A == 1)) + drop((colMeans(dptdbeta[A == 1,]) %*% beta_dot))) / n
}))

(crossprod(S))

z0 <- numeric(n)
z1 <- numeric(n)
p1 <- predict(fit, newdata = transform(d, A = 1))
p0 <- predict(fit, newdata = transform(d, A = 0))
m <- sum(A == 1)
for (i in 1:n) {
  z0[i] <- ((p0[i] - mean(p0[d$A == 1])) * (A[i] == 1) / (m / n) + (1/m * colSums(WeightIt:::.gradient(function(b) (X0[d$A == 1,] %*% b), .x = beta[1:ncol(X)]))) %*% beta_dot[,i])^2

  z1[i] <- ((p1[i] - mean(p1[d$A == 1])) * (A[i] == 1) / (m / n) + (1/m * colSums(WeightIt:::.gradient(function(b) (X1[d$A == 1,] %*% b), .x = beta[1:ncol(X)]))) %*% beta_dot[,i])^2
}
mean(z0)/n
mean(z1)/n
