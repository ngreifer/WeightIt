## code to prepare `msmdata` dataset goes here
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
  rbinom(n, 1, plogis(-2 - .8 * A_1 - 1.5 * A_2 - 2.2 * A_3 +
                        .7 * A_2 * A_3 + .5 * A_1 * A_3 +
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

set.seed(1234)
msmdata <- gen_dat(7500)

usethis::use_data(msmdata, overwrite = TRUE)

W <- weightitMSM(list(A_1 ~ X1_0 + X2_0,
                      A_2 ~ X1_1 + X2_1 +
                        A_1 + X1_0 + X2_0,
                      A_3 ~ X1_2 + X2_2 +
                        A_2 + X1_1 + X2_1 +
                        A_1 + X1_0 + X2_0),
                 data = msmdata, method = "glm")

#Generate counterfactual outcomes
gen_dat_cf <- function(n, A_1 = 0, A_2 = 0, A_3 = 0) {
  X1_0 <- gen_X1_0(n)
  X2_0 <- gen_X2_0(n, X1_0)

  A_1 <- rep(A_1, n)

  X1_1 <- gen_X1_1(n, X1_0, X2_0, A_1)
  X2_1 <- gen_X2_1(n, X1_0, X2_0, A_1, X1_1)

  A_2 <- rep(A_2, n)

  X1_2 <- gen_X1_2(n, A_2, X1_1, X2_1)
  X2_2 <- gen_X2_2(n, X1_1, X1_2, A_2)

  A_3 <- rep(A_3, n)

  Y_B <- gen_Y_B(n, A_1, A_2, A_3, X1_0, X2_0, X1_1, X2_1, X1_2, X2_2)

  data.frame(X1_0, X2_0, A_1, X1_1, X2_1, A_2, X1_2, X2_2, A_3, Y_B)
}

grid <- expand.grid(A_1 = 0:1, A_2 = 0:1, A_3 = 0:1)
apply(grid, 1, function(i) do.call(gen_dat_cf, c(list(n = 1e6), i))$Y_B |> mean()) |>
  matrix(ncol = 1, dimnames = list(do.call(paste, c(as.list(grid), list(sep = "|"))), "E[Y_B(A)]"))

marginaleffects::avg_predictions(fit, newdata = datagridcf(A_1 = 0:1, A_2 = 0:1, A_3 = 0:1),
                                 vcov = "HC3", wts = "w", by = c("A_1", "A_2", "A_3"))