n = 1e4
gen_X1_0 <- function(n) {rpois(n, 4)}
gen_X2_0 <- function(n, X1_0)
  rbinom(n, 1, .3 + .2 * (X1_0 < 3))

gen_A_1 <- function(n, X1_0, X2_0)
  rbinom(n, 1, plogis(-.4 + .3 * X1_0 - 1.9 * X2_0 + .1 * X1_0 * X2_0))

gen_X1_1 <- function(n, X1_0, X2_0, A_1)
  rpois(n, exp(.5 + .1 * A_1 + .7 * log(X1_0 + .2) - .1 * X2_0))
gen_X2_1 <- function(n, X1_0, X2_0, A_1, X1_1)
  rbinom(n, 1, plogis(.35 + .2 * (X1_0 < 3) + .2 * (X1_1 < 3) - .2 * A_1))

gen_A_2 <- function(n, X1_1, X2_1)
  rbinom(n, 1, plogis(-.4 + .3 * X1_1 - 1.9 * X2_1 + .1 * X1_1 * X2_1))

gen_X1_2 <- function(n, A_2, X1_1, X2_1)
  rpois(n, exp(.5 + .1 * A_2 + .7 * log(X1_1 + .2) - .1 * X2_1))
gen_X2_2 <- function(n, X1_1, X1_2, A_2)
  rbinom(n, 1, plogis(.35 + .2 * (X1_2 < 3) + .2 * (X1_1 < 3) - .2 * A_2))

gen_A_3 <- function(n, X1_2, X2_2, X1_1, X2_1)
  rbinom(n, 1, plogis(-.4 + .2 * X1_2 - 1.9 * X2_2 - .25 * X1_2 * X2_2 +
                        .1 * X1_1 - 1.3 * X2_1 + .25 * X1_1 * X2_1))

gen_Y_B <- function(n, A_1, A_2, A_3, X1_0, X2_0, X1_1, X2_1, X1_2, X2_2)
  rbinom(n, 1, plogis(1 - .1 * A_1 - .5 * A_2 - .8 * A_3 -
                        .4 * A_2 * A_3 - .2 * A_1 * A_3 +
                        .1 * X1_0 + .1 * X2_0 +
                        .2 * X1_1 + .2 * X2_1 +
                        .4 * X1_2 + .4 * X2_2 +
                        .25 * A_3 * X1_0 + .05 * A_3 * X1_2))

gen_dat <- function(n) {
  X1_0 <- gen_X1_0(n)
  X2_0 <- gen_X2_0(n, X1_0)

  A_1 <- gen_A_1(n, X1_0, X2_0)

  X1_1 <- gen_X1_1(n, X1_0, X2_0, A_1)
  X2_1 <- gen_X2_1(n, X1_0, X2_0, A_1, X1_1)

  A_2 <- gen_A_2(n, X1_1, X2_1)

  X1_2 <- gen_X1_2(n, A_2, X1_1, X2_1)
  X2_2 <- gen_X2_2(n, X1_1, X1_2, A_2)

  A_3 <- gen_A_3(n, X1_2, X2_2, X1_1, X2_1)

  Y_B <- gen_Y_B(n, A_1, A_2, A_3, X1_0, X2_0, X1_1, X2_1, X1_2, X2_2)

  data.frame(X1_0, X2_0, A_1, X1_1, X2_1, A_2, X1_2, X2_2, A_3, Y_B)
}

d <- gen_dat(1e4)

W <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                      A_2 ~ X1_1 + X2_1 +
                        A_1 + X1_0 + X2_0,
                      A_3 ~ X1_2 + X2_2 +
                        A_2 + X1_1 + X2_1 +
                        A_1 + X1_0 + X2_0),
                 data = d, method = "glm")
