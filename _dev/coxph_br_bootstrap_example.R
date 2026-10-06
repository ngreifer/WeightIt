# Bootstrapping coxph_weightit() with a rare outcome: the traditional bootstrap
# fails without bias reduction and works with it.
#
# Only 2 of the 123 control units have an event. The full-sample maximum partial
# likelihood estimate of the log hazard ratio is finite, but about 1 in 7 bootstrap
# resamples (the chance that neither event is drawn, (1 - 2/200)^200 = .13) contains
# no control events. In those resamples the likelihood is monotone, the estimate
# diverges (coxph() stops it at around 20 with a warning), and those few values
# dominate the bootstrap variance. With `br = TRUE`, every resample has a finite
# estimate, so the bootstrap variance is usable.

library("WeightIt")

# Simulated observational data with a rare outcome among the control units
set.seed(6)
n <- 200
X1 <- rnorm(n)
X2 <- rnorm(n)
A <- rbinom(n, 1, plogis(-.5 + .6 * X1 + .4 * X2))
time <- rexp(n, exp(-6.5 + 3 * A + .4 * X1))

d <- data.frame(A, X1, X2,
                time = pmin(time, 12),
                event = as.numeric(time <= 12))

with(d, table(A, event))

# ATE weights from a logistic regression propensity score model
W <- weightit(A ~ X1 + X2, data = d, method = "glm", estimand = "ATE")

# Traditional bootstrap without bias reduction
set.seed(123)
fit_ml <- coxph_weightit(survival::Surv(time, event) ~ A, data = d,
                         weightit = W, vcov = "BS", R = 500)

summary(fit_ml, ci = TRUE)

# Traditional bootstrap with bias reduction
set.seed(123)
fit_br <- coxph_weightit(survival::Surv(time, event) ~ A, data = d,
                         weightit = W, vcov = "BS", R = 500, br = TRUE)

summary(fit_br, ci = TRUE)

# Where the failure comes from: refitting by hand in the same resamples (the
# weights are re-estimated in each, as `vcov = "BS"` does) shows the replicate
# estimates
set.seed(123)
reps <- t(replicate(500, {
  ind <- sample.int(n, replace = TRUE)
  db <- d[ind, ]
  Wb <- weightit(A ~ X1 + X2, data = db, method = "glm", estimand = "ATE")

  c(control_events = sum(db$event[db$A == 0]),
    ml = suppressWarnings(coef(coxph_weightit(survival::Surv(time, event) ~ A, data = db,
                                              weightit = Wb, vcov = "none"))),
    br = coef(coxph_weightit(survival::Surv(time, event) ~ A, data = db,
                             weightit = Wb, vcov = "none", br = TRUE)))
}))

colnames(reps) <- c("control_events", "ml", "br")

# Resamples with no control events, and the estimates in them
table(no_control_events = reps[, "control_events"] == 0)
summary(reps[reps[, "control_events"] == 0, c("ml", "br")])

# The rest of the replicates agree closely
summary(reps[reps[, "control_events"] > 0, c("ml", "br")])

# Bootstrap SDs with and without the divergent resamples
apply(reps[, c("ml", "br")], 2L, sd)
apply(reps[reps[, "control_events"] > 0, c("ml", "br")], 2L, sd)
