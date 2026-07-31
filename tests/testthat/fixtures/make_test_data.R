.gen_data <- function(n) {
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
    2*A + 2*X[,1] + 2*X[,2] + 2*X[,3] + 1*X[,4] + 2*X[,5] + 1*X[,6] +
      .5*A*X[,1] + .5*A*X[,2] - .25*A*X[,3] + A*(X[,5] - .5) +
      rnorm(length(A), 0, 5)
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

  # Ordered outcome
  gen_Y_O <- function(Y_C) {
    factor(findInterval(Y_C, quantile(Y_C, seq(0, 1, length = 5)),
                        all.inside = TRUE),
           ordered = TRUE)
  }

  # Survival outcome
  gen_Y_S <- function(A, X) {
    LP_S <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]

    sqrt(-log(runif(length(A)))*2e4*exp(-LP_S))
  }

  # Sampling weight
  gen_SW <- function(A, X) {
    LP_SP <- .5 + log(.4) * X[,1] - log(.7) * X[,2] + log(.5) * X[,3] - log(.1) * A
    P_SW <- plogis(LP_SP)

    1 / P_SW
  }

  X <- gen_X(n)
  Ac <- gen_Ac(X)
  A <- gen_A(Ac)
  Am <- gen_Am(A)

  Y_C <- gen_Y_C(A, X)
  Y_B <- gen_Y_B(A, X)
  Y_O <- gen_Y_O(Y_C)
  Y_S <- gen_Y_S(A, X)

  SW <- gen_SW(A, X)

  d <- data.frame(A, Am, Ac, X, Y_C, Y_B, Y_O, Y_S, SW)

  d$X6 <- factor(cut(d$X6, 4), labels = LETTERS[1:4])

  #Grouping factor for multilevel/random-effects models. Defined as a fine
  #quantile-binning of X1 (a confounder for treatment and all outcomes), so that
  #cluster membership is highly correlated with X1 and can substitute for it when
  #adjusting for confounding. This is a deterministic function of X1 (no RNG
  #draw), so every other variable above is unchanged from the fixed-effects
  #fixture.
  #
  #The bin count scales with n to keep ~30 units per cluster. Below roughly that,
  #`lmer()` cannot pin down the cluster variance well enough for the density-ratio
  #weights of a continuous treatment to improve balance, so the multilevel tests
  #start failing for reasons that have nothing to do with the code under test.
  nclust <- max(5L, round(n / 30))

  d$cluster <- factor(cut(d$X1,
                          breaks = quantile(d$X1, seq(0, 1, length.out = nclust + 1L)),
                          include.lowest = TRUE, labels = FALSE))

  d
}

set.seed(1234)

#n is a test-suite runtime lever, not a statistical choice. The optimization-based
#methods (`energy`, `cfd`) build a dense n x n kernel matrix and hand it to OSQP, and
#`expect_M_parts_okay()` builds n x p numerical Jacobians, so their cost grows faster
#than n. Going from 2000 to 750 cut those files by more than 10x while leaving ~130
#units in the smallest `Am` level -- enough for every method to converge across all
#estimands and for the balance assertions to stay meaningful. Regenerate rather than
#subset if changing it: `Y_O` is defined by quantiles of `Y_C` and `cluster` by
#quantiles of `X1`, so taking the first k rows distorts both.
test_data <- .gen_data(750)

saveRDS(test_data, testthat::test_path("fixtures", "test_data.rds"))