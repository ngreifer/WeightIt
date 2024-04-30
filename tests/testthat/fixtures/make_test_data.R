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

  # Survival outcome
  gen_Y_S <- function(A, X) {
    LP_S <- -2 + log(2.4)*A + log(2)*X[,1] + log(2)*X[,2] + log(2)*X[,3] + log(1.5)*X[,4] + log(2.4)*X[,5] + log(1.5)*X[,6]
    sqrt(-log(runif(length(A)))*2e4*exp(-LP_S))
  }

  gen_SW <- function(A, X) {
    LP_SP <- .5 + log(.4) * X[,1] - log(.7) * X[,2] + log(.5) * X[,3] - log(.1) * A
    P_SW <- plogis(LP_SP)

    1/P_SW
  }

  X <- gen_X(n)
  Ac <- gen_Ac(X)
  A <- gen_A(Ac)
  Am <- gen_Am(A)

  Y_C <- gen_Y_C(A, X)
  Y_B <- gen_Y_B(A, X)
  Y_S <- gen_Y_S(A, X)

  SW <- gen_SW(A, X)

  d <- data.frame(A, Am, Ac, X, Y_C, Y_B, Y_S, SW)

  d$X6 <- factor(cut(d$X6, 4), labels = LETTERS[1:4])

  d
}

set.seed(1234)

test_data <- .gen_data(2000)

saveRDS(test_data, testthat::test_path("fixtures", "test_data.rds"))