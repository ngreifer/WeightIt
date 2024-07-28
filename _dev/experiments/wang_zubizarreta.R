#SBW using algorithms in Wang & Zubizarreta (2020) and Kallberg & Waernbaum (2023)

#Fast optweight with exact balance

xl <- colnames(data)[substr(colnames(data),1,1) == "X"]
tr <- data$tr
N  <- length(tr)

X <- cbind(1,as.matrix(data[xl]))

# sum of complete sample
b <- colSums(X)

X1 <- X*tr # treated
X0 <- X*(1-tr) # control

A1 <- t(X1) %*% X1
A0 <- t(X0) %*% X0

# dual parameters

k1 <- c(mean(tr),0*is.na(xl))
k0 <- c(mean(1-tr),0*is.na(xl))

conv <- F

ch.solve <- function(A,b) {
  chA <- chol(A)
  backsolve(chA, forwardsolve(t(chA), b, transpose = F))
}

try({k1 <- ch.solve(A1,b); k0 <- ch.solve(A0,b); conv <- T})

w1 <- X1 %*% k1
w0 <- X0 %*% k0

w <- w1 + w0

covs <- names(lalonde)[-c(1, 4, 9)]
X <- cbind(1, as.matrix(lalonde[covs]))

W <- weightit(reformulate(covs, "treat"), lalonde, method = "optweight")

######
#Wang & Zubzizarreta (2018)
covs <- names(lalonde)[-c(1, 4, 9)]
X <- cbind(1, scale(as.matrix(lalonde[covs]), center = T, scale = MatchIt:::pooled_sd(as.matrix(lalonde[covs]),
                                                                                      lalonde$treat, cont = "equal")))
A <- lalonde$treat

p <- function(x, A) {
  # -exp(-x-1)
  -x^2/4 + x/sum(A)
}
p_ <- function(x, A) {
  # exp(-x-1)
  -x/2 + 1/sum(A)
}

fun <- function(g, X, A, delta = 0) {
  Xg <- drop(X %*% g)
  delta <- c(0, rep(delta, length(g) - 1))
  sum(-A * p(Xg, A) + Xg) + sum(abs(g) * delta)
}
f <- function(g, X, A, delta = rep(0, ncol(X))) {
  fun(g[1:ncol(X)], X, A == 1, delta) +
    fun(g[-(1:ncol(X))], X, A == 0, delta)
}

wfun <- function(g, X, A) {
  p_(drop(X %*% g), A)
}

opt <- optim(rep(0, 2 * ncol(X)),
             f, X = X, A = A, delta = .1,
             method = "BFGS",
             control = list(maxit = 1e3, reltol = 1e-10))

w <- 0 * A
w[A == 1] <- wfun(opt$par[1:ncol(X)], X[A == 1,], A == 1)
w[A == 0] <- wfun(opt$par[-(1:ncol(X))], X[A == 0,], A == 0)

opt1 <- optim(rep(0, ncol(X)),
              fun, X = X, A = A == 1, delta = .2,
              method = "BFGS",
              control = list(maxit = 1e3, reltol = 1e-10))
opt0 <- optim(rep(0, ncol(X)),
              fun, X = X, A = A == 0, delta = .2,
              method = "BFGS",
              control = list(maxit = 1e3, reltol = 1e-10))

w <- 0 * A
w[A == 1] <- wfun(opt1$par, X[A == 1,], A == 1)
w[A == 0] <- wfun(opt0$par, X[A == 0,], A == 0)

# w <- 0 * A
# w[A == 1] <- wfun(attr(W, "Mparts")$b[-(1:ncol(X))], X[A == 1,])
# w[A == 0] <- wfun(attr(W, "Mparts")$b[1:ncol(X)], X[A == 0,])

W <- weightit(reformulate(covs, "treat"), lalonde, method = "optweight", tols = 0)
bal.tab(W, weights = w, cont  = "std")
