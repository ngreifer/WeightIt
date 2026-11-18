# Checks for _dev/plan-modified-treatment-policies.md (2026-11-18). Needs adrftools (for its data file), cobalt, and osqp.
# Part 1: density-ratio vs. stacked weights, weighted g-computation, ADE identity.
suppressPackageStartupMessages({library(WeightIt)})
set.seed(1)
n <- 4000
X1 <- rnorm(n); X2 <- rbinom(n, 1, .4)
A <- 1 + .8*X1 - .5*X2 + rnorm(n)        # normal linear GPS, sd 1
m <- function(a, x1, x2) 2 + a + .5*a^2 + x1 + a*x1 - x2
Y <- m(A, X1, X2) + rnorm(n)
d <- data.frame(A, X1, X2, Y)
delta <- .5
truth <- mean(m(A + delta, X1, X2))

## (1) density-ratio weights from a normal linear GPS
fit <- lm(A ~ X1 + X2, d); mu <- fitted(fit); s <- sqrt(mean(resid(fit)^2))
w_dr <- dnorm(A - delta, mu, s) / dnorm(A, mu, s)
## log-linear in (A, X)?
print(coef(lm(log(w_dr) ~ A + X1 + X2)))
cat("delta/s^2 =", delta/s^2, "\n")

## (2) stacked classification with WeightIt binary methods, ATT toward shifted copy
st <- rbind(transform(d, S = 0), transform(d, A = A + delta, S = 1))
W_glm  <- weightit(S ~ A + X1 + X2, data = st, method = "glm",  estimand = "ATT")
W_ebal <- weightit(S ~ A + X1 + X2 + A:X1 + A:X2 + I(A^2), data = st, method = "ebal", estimand = "ATT")
w_glm <- W_glm$weights[st$S == 0]; w_ebal <- W_ebal$weights[st$S == 0]
cat("cor(density-ratio, stacked glm) =", cor(w_dr, w_glm), "\n")

est <- function(w) sum(w * Y) / sum(w)
print(round(c(truth = truth, naive_gcomp_unadj = mean(Y) + delta * 0,
              ipw_density = est(w_dr), ipw_stack_glm = est(w_glm), ipw_stack_ebal = est(w_ebal)), 3))

## (3) DR via weighted regression + g-computation at shifted values (misspecified outcome model)
o <- lm(Y ~ A + X1 + X2, data = d, weights = w_dr)
o0 <- lm(Y ~ A + X1 + X2, data = d)
gc <- function(fit) mean(predict(fit, transform(d, A = A + delta)))
print(round(c(gcomp_unweighted_misspec = gc(o0), gcomp_mtpweighted_misspec = gc(o)), 3))

## (4) ADE: normal-GPS score weights reproduce the OLS coefficient
omega <- (A - mu) / s^2
cat("ADE via score weights:", mean(omega * Y), " OLS coef on A:", coef(o0)["A"],
    " true ADE:", mean(1 + A + X1), "\n")
## derivative of shift weights at 0 equals omega
eps <- 1e-5
cat("max |d/ddelta w at 0 - omega| =", max(abs((dnorm(A - eps, mu, s)/dnorm(A, mu, s) - 1)/eps - omega)), "\n")

# Part 2: stacked energy weights and balance against the shifted sample in cobalt.
rm(list = ls())
suppressPackageStartupMessages({library(WeightIt); library(cobalt)})
set.seed(2)
n <- 1000
X1 <- rnorm(n); X2 <- rbinom(n, 1, .4)
A <- 1 + .8*X1 - .5*X2 + rnorm(n)
m <- function(a, x1, x2) 2 + a + .5*a^2 + x1 + a*x1 - x2
Y <- m(A, X1, X2) + rnorm(n)
d <- data.frame(A, X1, X2, Y); delta <- .5
truth <- mean(m(A + delta, X1, X2))
st <- rbind(transform(d, S = 0), transform(d, A = A + delta, S = 1))

W_en <- weightit(S ~ A + X1 + X2, data = st, method = "energy", estimand = "ATT")
w_en <- W_en$weights[st$S == 0]
fit <- lm(A ~ X1 + X2, d); s <- sqrt(mean(resid(fit)^2))
w_dr <- dnorm(A - delta, fitted(fit), s) / dnorm(A, fitted(fit), s)
est <- function(w) sum(w * Y) / sum(w)
print(round(c(truth = truth, unweighted = mean(Y), density = est(w_dr), energy_stack = est(w_en)), 3))
cat("ESS/n density:", ESS(w_dr)/n, " theory exp(-(delta/s)^2):", exp(-(delta/s)^2), " energy:", ESS(w_en)/n, "\n")

## balance against the shifted target via stacking, in cobalt as it is
b <- bal.tab(S ~ A + X1 + X2 + A:X1 + A:X2, data = st,
             weights = list(density = c(w_dr, rep(1, n)), energy = c(w_en, rep(1, n))),
             estimand = "ATT", s.d.denom = "treated", stats = c("m", "ks"), un = TRUE)
print(b)
for (w in list(rep(1, n), w_dr, w_en))
  print(bal.compute(st[c("A", "X1", "X2")], treat = st$S, weights = c(w, rep(1, n)), stat = "energy.dist"))

# Part 3: ADRF weights used for a shift (unadjusted weighted g-computation).
rm(list = ls())
## Do correctly specified ADRF (GPS) weights work for a shift intervention?
set.seed(3)
n <- 2e5                                  # large n so differences are bias, not noise
X <- rnorm(n)
A <- .8 * X + rnorm(n)                    # A | X ~ N(.8 X, 1); A ~ N(0, 1.64)
delta <- .5
w_adrf <- dnorm(A, 0, sqrt(1.64)) / dnorm(A, .8 * X, 1)   # true ADRF weights
w_mtp  <- dnorm(A - delta, .8 * X, 1) / dnorm(A, .8 * X, 1) # true MTP weights

run <- function(m, label) {
  Y <- m(A, X) + rnorm(n)
  truth_level <- mean(m(A + delta, X)); truth_contrast <- truth_level - mean(m(A, X))
  ## weighted g-computation with no covariates: flexible marginal model in A
  f <- lm(Y ~ poly(A, 4, raw = TRUE), weights = w_adrf)
  th <- function(a) predict(f, data.frame(A = a))
  c(truth_level = truth_level,
    adrf_level = mean(th(A + delta)),
    mtp_level = sum(w_mtp * Y) / sum(w_mtp),
    truth_contrast = truth_contrast,
    adrf_contrast_vs_thetaA = mean(th(A + delta)) - mean(th(A)),
    adrf_contrast_vs_meanY = mean(th(A + delta)) - mean(Y),
    mtp_contrast = sum(w_mtp * Y) / sum(w_mtp) - mean(Y)) |> round(3) |> print()
}
cat("No effect modification: m = A + .5 A^2 + X\n");      run(\(a, x) a + .5 * a^2 + x)
cat("Modification linear in A: m = A + A*X + X\n");        run(\(a, x) a + a * x + x)
cat("Modification nonlinear in A: m = A + .5 A^2 X + X\n"); run(\(a, x) a + .5 * a^2 * x + x)

# Part 4: ADE weights on adrftools' nhanes3lead (constraints, two-group form).
rm(list = ls())
## ADE weights on nhanes3lead: A = logBLL, Y = Math
suppressPackageStartupMessages(library(cobalt))
load("/Users/NoahGreifer/Dropbox/Research/R/adrftools/data/nhanes3lead.rda")
d <- nhanes3lead; n <- nrow(d)
X <- model.matrix(~ Age + Male + Race + PIR + Enough_Food + Smoke_in_Home + Smoke_Pregnant + NICU, d)[, -1]
A <- d$logBLL; Y <- d$Math

## basis of the outcome model and the derivative of each column in A
G  <- cbind(1, X, A = A, A * X, A2 = A^2)
dG <- cbind(0, X * 0, A = 1, X, A2 = 2 * A)
nm <- c("1", colnames(X), "A", paste0("A:", colnames(X)), "A2"); colnames(G) <- colnames(dG) <- nm
tgt <- colMeans(dG)                            # what (1/n) sum(w * g) must equal

## (1) GPS score weights, normal linear model for A | X
fA <- lm(A ~ X); s2 <- mean(resid(fA)^2)
w_gps <- resid(fA) / s2
## (2) minimum-norm balancing weights: min sum w^2 s.t. (1/n) G'w = tgt
P <- solve(crossprod(G), n * tgt)
w_bal <- drop(G %*% P)
## (3) GPS weights calibrated to the same constraints (smallest change)
w_cal <- w_gps + drop(G %*% solve(crossprod(G), n * tgt - crossprod(G, w_gps)))

cons <- function(w) drop(crossprod(G, w)) / n
tab <- cbind(target = tgt, gps = cons(w_gps), balancing = cons(w_bal), calibrated = cons(w_cal))
print(round(tab[c("1", "Age", "PIR", "Male", "A", "A:Age", "A:PIR", "A2"), ], 3))

## regression ADE from the outcome model with the same basis
fY <- lm(Y ~ G - 1); ade_reg <- mean(dG %*% coef(fY))
ade_dr <- ade_reg + mean(w_gps * (Y - fitted(fY)))
print(round(c(gps = mean(w_gps * Y), balancing = mean(w_bal * Y), regression = ade_reg,
              calibrated = mean(w_cal * Y), augmented = ade_dr), 4))

## the weights as a balanced comparison of two groups
grp <- function(w) {
  pos <- w > 0; p <- abs(w)
  dA <- weighted.mean(A[pos], p[pos]) - weighted.mean(A[!pos], p[!pos])
  dY <- weighted.mean(Y[pos], p[pos]) - weighted.mean(Y[!pos], p[!pos])
  c(share_pos = mean(pos), mass = sum(p[pos]) / n, mean_A_pos = weighted.mean(A[pos], p[pos]),
    mean_A_neg = weighted.mean(A[!pos], p[!pos]), dA = dA, dY = dY, ratio = dY / dA)
}
print(round(rbind(gps = grp(w_gps), balancing = grp(w_bal), calibrated = grp(w_cal)), 3))
print(summary(w_bal)); cat("mean w^2:", mean(w_gps^2), mean(w_bal^2), mean(w_cal^2), "\n")

b <- bal.tab(as.data.frame(X), treat = as.integer(w_bal > 0), weights = abs(w_bal),
             s.d.denom = "pooled", binary = "std", un = TRUE, stats = c("m", "ks"))
print(b)

# Part 5: kernel and moment balancing ADE weights in simulation.
rm(list = ls())
## Kernel balancing for the ADE: Gaussian vs. energy kernel, against moment balancing
gen <- function(n) {
  X <- matrix(rnorm(n * 3), n)
  A <- .6 * X[, 1] + .4 * X[, 2]^2 - .4 + rnorm(n)
  m  <- function(a) a + .3 * a^2 * X[, 1] + sin(2 * X[, 2]) + .5 * a * X[, 3] + X[, 3]
  dm <- function(a) 1 + .6 * a * X[, 1] + .5 * X[, 3]
  list(X = X, A = A, Y = m(A) + rnorm(n), ade = mean(dm(A)))
}
## worst-case imbalance over the unit ball of kernel k:
##   (1/n^2) w'Kw - (2/n^2) w'b, b_i = sum_j d/da_j k(Z_i, Z_j); add (lambda/n^2) w'w; sum(w) = 0
kern_w <- function(Z, sdA, type, lambda = 1) {
  n <- nrow(Z); D <- as.matrix(dist(Z)); dA <- outer(Z[, 1], Z[, 1], "-")  # a_i - a_j
  if (type == "energy") { K <- -D; Db <- dA / D; Db[!is.finite(Db)] <- 0 }
  else { h <- median(D); K <- exp(-D^2 / (2 * h^2)); Db <- K * dA / h^2 }
  b <- rowSums(Db) / sdA                       # derivative w.r.t. the unstandardized A
  M <- rbind(cbind(K + lambda * diag(n), 1), c(rep(1, n), 0))
  solve(M, c(b, 0))[1:n]
}
one <- function(n = 500) {
  s <- gen(n); X <- s$X; A <- s$A; Y <- s$Y
  Z <- scale(cbind(A, X)); sdA <- sd(A)
  G <- cbind(1, X, A, A * X, A^2); dG <- cbind(0, X * 0, 1, X, 2 * A)
  w_bal <- drop(G %*% solve(crossprod(G), colSums(dG)))
  fA <- lm(A ~ X); w_gps <- resid(fA) / mean(resid(fA)^2)
  fA2 <- lm(A ~ X + I(X[, 2]^2)); w_gps2 <- resid(fA2) / mean(resid(fA2)^2)
  est <- c(gps_misspec = mean(w_gps * Y), gps_correct = mean(w_gps2 * Y),
           moment_balance = mean(w_bal * Y),
           gaussian = mean(kern_w(Z, sdA, "gaussian") * Y),
           energy = mean(kern_w(Z, sdA, "energy") * Y))
  est - s$ade
}
set.seed(4); t0 <- Sys.time()
R <- t(replicate(100, one()))
print(round(rbind(bias = colMeans(R), sd = apply(R, 2, sd), rmse = sqrt(colMeans(R^2))), 3))
cat("time:", format(Sys.time() - t0), "\n")

# Part 6: energy shift weights as the shift shrinks (uses gen() and kern_w() from Part 5).
## Energy shift weights (Jiang & Huling) as the shift shrinks: does (mean_w Y - mean Y)/delta approach the ADE?
suppressPackageStartupMessages(library(osqp))
energy_shift_w <- function(A, X, delta, lambda = 1) {
  n <- length(A); mA <- mean(A); sA <- sd(A); sX <- apply(X, 2, sd)
  Z  <- cbind((A - mA) / sA, sweep(X, 2, sX, "/"))
  Zq <- cbind((A + delta - mA) / sA, sweep(X, 2, sX, "/"))
  D  <- as.matrix(dist(Z))
  B  <- sqrt(outer(rowSums(Z^2), rowSums(Zq^2), "+") - 2 * tcrossprod(Z, Zq)); B[is.nan(B)] <- 0
  ## n^2 E = -w'Dw + 2 w'b + const; on sum(w) = n, -w'Dw = -(w-1)'D(w-1) - 2 (D1)'w + const, and -JDJ is PSD
  J <- diag(n) - 1 / n; P <- (-J %*% D %*% J + lambda * diag(n)) / n^2; q <- (2 * rowSums(B) - 2 * rowSums(D)) / n^2
  Am <- rbind(rep(1, n), diag(n))
  s <- solve_osqp(2 * P, q, Am, c(n, rep(0, n)), c(n, rep(Inf, n)),
                  osqpSettings(verbose = FALSE, eps_abs = 1e-9, eps_rel = 1e-9, max_iter = 2e5))
  pmax(s$x, 0)
}
set.seed(5); t0 <- Sys.time()
cs <- c(1, .5, .25, .1, .05)            # shift in conditional-SD units (sigma = 1 here)
one <- function() {
  s <- gen(500); A <- s$A; X <- s$X; Y <- s$Y
  r <- sapply(cs, function(cc) { w <- energy_shift_w(A, X, cc); c(fd = (sum(w * Y) / sum(w) - mean(Y)) / cc, ess = sum(w)^2 / sum(w^2)) })
  wk <- kern_w(scale(cbind(A, X)), sd(A), "energy")
  c(r["fd", ] - s$ade, mean(wk * Y) - s$ade, r["ess", ])
}
R <- t(replicate(12, one()))
colnames(R) <- c(paste0("fd_", cs), "energy_ADE_weights", paste0("ess_", cs))
print(round(rbind(bias = colMeans(R), sd = apply(R, 2, sd)), 3))
cat("time:", format(Sys.time() - t0), "\n")
