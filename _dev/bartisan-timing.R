# Timing of the bartisan-based method: where the time goes in the continuous case,
# and how the method scales with the sample size. Run from the package root.
source("_dev/bartisan-inject.R")
suppressPackageStartupMessages(library("cobalt"))
test_data <- readRDS("tests/testthat/fixtures/test_data.rds")
f_c <- Ac ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9

for (reps in c(1L, 4L)) {
  d <- test_data[rep(seq_len(nrow(test_data)), reps), ]
  cat(sprintf("--- n = %d ---\n", nrow(d)))
  set.seed(1)
  t_fit <- system.time(fit <- bartisan::bartisan(f_c, data = d, family = bartisan::dpm()))[["elapsed"]]
  set.seed(1)
  t_total <- system.time(W <- weightit(f_c, data = d, method = "bart"))[["elapsed"]]
  r <- c(un = max(abs(col_w_cov(W$covs, W$treat, std = TRUE))),
         adj = max(abs(col_w_cov(W$covs, W$treat, weights = W$weights, std = TRUE))))
  cat(sprintf("continuous: dpm fit alone %.1f s | weightit() total %.1f s | max |cor| %.3f -> %.3f | ESS %.0f\n",
              t_fit, t_total, r[1], r[2], ESS(W$weights)))
  set.seed(1); t_b <- system.time(weightit(A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9, data = d, method = "bart"))[["elapsed"]]
  set.seed(1); t_m <- system.time(weightit(Am ~ X1 + X2 + X3 + X4 + X5, data = d, method = "bart"))[["elapsed"]]
  cat(sprintf("binary %.1f s | multi-category %.1f s\n", t_b, t_m))
}
