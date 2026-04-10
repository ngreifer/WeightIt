#Hypothesis testing for functions
d <- readRDS(test_path("fixtures", "test_data.rds"))

fit <- glm(Y_B ~ poly(Ac, 4) * (X1 + X2 + X3 + X4), data = d)

library(marginaleffects)
plot_predictions(fit,
                 by = "Ac",
                 newdata = datagrid(Ac = v, grid_type = "counterfactual"))

v <- seq(-9, 6, length.out = 101)

s <- avg_slopes(fit, variables = "Ac",
                by = "Ac", vcov = "HC3",
                newdata = datagrid(Ac = v, grid_type = "counterfactual"))

ggplot(s, aes(x = Ac)) + geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .3)


f <- s$estimate
g <- 0 * f
V <- vcov(s)

de <- v[2] - v[1]
T_ <- sum(.5 * (f[-1]^2 + f[-length(f)]^2) * diff(v))
ev <- eigen(V * de, T, T)$values

(f %*% MASS::ginv(V) %*% f) |> pchisq(df = length(f), lower.tail = F)

CompQuadForm::davies(T_, ev)

########
n <- 1e3
A <- runif(n, -2, 2)
X1 <- rbinom(n, 1, .5)
X2 <- rnorm(n)

Y <- (X1 == 1) * ((A + .75)^2 - 1) + (X1 == 0) * (-A^2 + 1) +
  .5 * X2 * A + X2 + X1 + 7 * rnorm(n)

d <- data.frame(A, X1, X2, Y)

fit <- lm(Y ~ poly(A, 2) * (X1 + X2), data = d)

v <- seq(-2, 2, by = .05)
p1 <- avg_predictions(fit, variables = list(A = v, X1 = 0:1),
                      vcov = "HC3")

ggplot(p1, aes(x = A, fill = factor(X1))) +
  geom_line(aes(y = estimate, color = factor(X1))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .3)

p0 <- avg_predictions(fit, variables = list(A = v),
                      vcov = "HC3")
ggplot(p0, aes(x = A)) +
  geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .3)

fitr <- lm(Y ~ (X1 + X2), data = d)
anova(fit, fitr)

s1 <- avg_slopes(fit, variables = "A",
                 by = c("A", "X1"), vcov = "HC3",
                 newdata = datagrid(A = v, X1 = 0:1, grid_type = "counterfactual"))

ggplot(s1, aes(x = A, fill = factor(X1))) +
  geom_line(aes(y = estimate, color = factor(X1))) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .3)


s0 <- avg_slopes(fit, variables = "A",
                 by = "A",
                 newdata = datagrid(model = fit, A = v, grid_type = "counterfactual"))

ggplot(s0, aes(x = A)) + geom_line(aes(y = estimate)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .3)

f <- s0$estimate^2
V <- vcov(s0)

de <- mean(diff(v))
T_ <- sum(.5 * (f[-1] + f[-length(f)]) * diff(v))
ev <- eigen(V * de, T, T)$values

CompQuadForm::imhof(T_, ev)

est <- function(a, fit, d) {
  n <- nrow(d)
  ii <- seq_len(n)
  new <- d[rep(ii, 2 * length(a)),]
  iii <- seq_len(n * length(a))
  new$A <- c(rep(a + 1e-5, each = n), rep(a - 1e-5, each = n))
  xx <- model.matrix(fit$terms, data = new)

  p <- drop(xx %*% coef(fit))

  (tapply(p[iii], rep(seq_along(a), each = n), mean) -
      tapply(p[-iii], rep(seq_along(a), each = n), mean)) / 2e-5
}

covf <- function(a, fit, d) {
  V <- vcov(fit)
  n <- nrow(d)
  ii <- seq_len(n)
  new <- d[rep(ii, 2 * length(a)),]
  iii <- seq_len(n * length(a))
  new$A <- c(rep(a + 1e-5, each = n), rep(a - 1e-5, each = n))
  xx <- model.matrix(fit$terms, data = new)

  gr <- (rowsum(xx[iii,], rep(seq_along(a), each = n), FALSE) -
           rowsum(xx[-iii,], rep(seq_along(a), each = n), FALSE)) / (n * 2e-5)

  gr %*% V %*% t(gr)
}

v <- seq(-2, 2, length.out = 1000)
f <- est(v, fit = fit, d= d)^2
V <- covf(v, fit = fit, d= d)
de <- (v[length(v)] - v[1]) / length(v)
# T_ <- sum(.5 * (f[-1] + f[-length(f)]) * diff(v))
T_ <- integrate(\(...) est(..., fit = fit, d= d)^2, v[1], v[length(v)])$value

ev <- (eigen(V * de, T, T)$values)
CompQuadForm::imhof(T_, ev[ev > 1e-10])

#Satterthwaite/Patnaik approx
df <- sum(ev)^2/sum(ev^2)
sc <- sum(ev^2)/sum(ev)

df <- sum(diag(V))^2/sum(V^2)
sc <- sum(diag(V)) * de/df

pchisq(T_ / sc, df = df, lower.tail = F)

get_p_vals <- function() {
  n <- 1e3
  A <- runif(n, -2, 2)
  X1 <- rbinom(n, 1, .5)
  X2 <- rnorm(n)

  Y <- (X1 == 1) * (A^2 - 1) + (X1 == 0) * (-A^2 + 1) +
    .5 * X2 * A + X2 + X1 + 7 * rnorm(n)

  d <- data.frame(A, X1, X2, Y)

  fit <- lm(Y ~ poly(A, 2) * (X1 + X2), data = d)

  v <- seq(-2, 2, length.out = 200)
  f <- est(v, fit = fit, d = d)^2
  V <- covf(v, fit = fit, d = d)
  de <- (v[length(v)] - v[1]) / length(v)

  T_ <- integrate(\(...) est(..., fit = fit, d = d)^2, v[1], v[length(v)])$value

  ev <- eigen(V * de, T, T)$values
  pq <- CompQuadForm::imhof(T_, ev[ev > 1e-10])$Qq

  #Satterthwaite approx
  df <- sum(diag(V))^2/sum(V^2)
  sc <- sum(diag(V)) * de/df

  pc <- pchisq(T_ / sc, df = df, lower.tail = F)

  c(pq = pq, pc = pc)
}

res <- do.call("rbind", pbapply::pblapply(1:5000, function(i) get_p_vals(),
                                          cl = 3))
