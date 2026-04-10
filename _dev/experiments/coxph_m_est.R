#Cox PH model using M-estimation

# Shu et al (2021)
# Y_S = T (event time)
# Y_event = δ (event occurred = 1)

test_data <- readRDS(testthat::test_path("fixtures", "test_data.rds"))
# cens_time <- runif(nrow(test_data), min(test_data$Y_S), max(test_data$Y_S))
cens_time <- 400
test_data$Y_death <- as.numeric(test_data$Y_S <= cens_time)
test_data$Y_S <- pmin(test_data$Y_S, cens_time)

X <- as.matrix(test_data[c("A", "X5")])

m_estimate <- function(start, psi, ...) {
  gr <- function(B, ...) {
    colSums(psi(B = B, ...))
  }

  out <- rootSolve::multiroot(gr, start = 1.1 * start,
                              ...,
                              maxiter = 1e5, rtol = 1e-9, atol = 1e-9, ctol = 1e-9)

  hess <- WeightIt:::.gradient(gr, .x = out$root, .method = "rich", ...)
  meat <- crossprod(psi(B = out$root, ...))

  V <- tryCatch({solve(hess, t(solve(hess, meat)))},
                error = function(e) {
                  matrix(NA_real_, length(out$root), length(out$root))
                })

  list(est = out$root,
       vcov = V)
}

gr <- function(B, ..., .psi) {
  colSums(.psi(B = B, ...))
}



get_V <- function(hess, meat) {
  tryCatch({solve(hess, t(solve(hess, meat)))},
           error = function(e) {
             matrix(NA_real_, length(out$root), length(out$root))
           })
}

vcov(fit)

psi_1 <- function(B, Y_S, Y_death, X, W) {
  p <- exp(drop(X %*% B))

  S0 <- vapply(Y_S, function(y) sum((Y_S >= y) * W * p), numeric(1L))
  S1 <- do.call("cbind", lapply(1:ncol(X), function(k) {
    vapply(Y_S, function(y) sum((Y_S >= y) * W * p * X[,k]), numeric(1L))
  }))

  M <- Y_death * (X - S1 / S0)

  for (i in 1:nrow(X)) {
    z <- (Y_S <= Y_S[i]) & Y_death

    M[i, ] <- M[i, , drop = FALSE] -
      p[i] * X[i, , drop = FALSE] * sum(W[z] / S0[z]) +
      p[i] * colSums(W[z] * S1[z, , drop = FALSE] / S0[z]^2)
  }

  W * M
}

w <- test_data$SW

out_1 <- m_estimate(start = rep(0, ncol(X)),
                    psi = psi_1,
                    Y_S = test_data$Y_S,
                    Y_death = test_data$Y_death,
                    X = X, W = w)
out_1

fit <- coxph(Surv(Y_S, Y_death) ~ A + X5, data = test_data,
             weights = SW, ties = "bres",
             robust = TRUE, control = coxph.control(eps = 1e-12, iter.max = 500, toler.chol = 1e-13))
coef(fit)
vcov(fit)

opt_1 <- rootSolve::multiroot(gr, start = rep(0, ncol(X)), .psi = psi_1,
                            Y_S = test_data$Y_S,
                            Y_death = test_data$Y_death,
                            X = X, W = w)

hess_1 <- WeightIt:::.gradient(gr, .x = opt_1$root, .psi = psi_1,
                               Y_S = test_data$Y_S,
                               Y_death = test_data$Y_death,
                               X = X, W = w,
                               .method = "rich")

meat_1 <- crossprod(psi_1(B = opt_1$root, Y_S = test_data$Y_S,
                            Y_death = test_data$Y_death,
                            X = X, W = w))

get_V(hess_1, meat_1)

all.equal(residuals(fit, type = "score", weighted = T),
          psi_1(B = opt_1$root, Y_S = test_data$Y_S, Y_death = test_data$Y_death, X = X, W = w),
          check.attributes = F)

psi_1b <- function(B, Y_S, Y_death, X, W) {
  n  <- length(Y_S)
  p  <- exp(drop(X %*% B))
  Wp <- W * p

  # Map each unique Y_S value to an integer rank (ties get same rank)
  # This lets us treat the risk set condition as a comparison of integer ranks
  ranks <- rank(Y_S, ties.method = "min")

  # Aggregate Wp and Wp*X by rank, then compute cumulative sums from the top
  # S0[i] = sum of Wp over all j where Y_S[j] >= Y_S[i]
  #       = sum of rank_Wp[r] for r >= ranks[i]  (cumsum from top)
  rank_Wp   <- rowsum(Wp,     ranks, reorder = TRUE)        # (max(ranks) x 1)
  rank_WpX  <- rowsum(Wp * X, ranks, reorder = TRUE)        # (max(ranks) x ncol(X))

  cum_Wp  <- rev(cumsum(rev(rank_Wp)))                      # cumsum from top rank down
  cum_WpX <- apply(rank_WpX, 2L, function(col) rev(cumsum(rev(col))))

  S0 <- cum_Wp[ranks]                                       # (n x 1)
  S1 <- cum_WpX[ranks, , drop = FALSE]                      # (n x ncol(X))

  M <- Y_death * (X - S1 / S0)

  # For the loop correction terms, we need for each i:
  #   term1[i] = sum_{j: Y_S[j] <= Y_S[i], Y_death[j]} W[j] / S0[j]
  #   term2[i] = colSums of W[j] * S1[j,] / S0[j]^2 over same j

  WdS0     <- W * Y_death / S0
  WS1dS0sq <- W * Y_death * S1 / S0^2

  rank_WdS0     <- rowsum(WdS0,     ranks, reorder = TRUE)
  rank_WS1dS0sq <- rowsum(WS1dS0sq, ranks, reorder = TRUE)

  cum_WdS0     <- cumsum(rank_WdS0)                          # cumsum from bottom rank up
  cum_WS1dS0sq <- apply(rank_WS1dS0sq, 2L, cumsum)

  term1 <- cum_WdS0[ranks]
  term2 <- cum_WS1dS0sq[ranks, , drop = FALSE]

  M <- M - p * term1 * X + p * term2

  W * M
}
