#CBPS as penalized likelihood

##CBPS objective for logistic link, ATE; Zhao (2019)
L <- function(theta, X, A) {
  p <- plogis(drop(X %*% theta))

  sum(A * log(p) + (1 - A) * log(1 - p) - A * log((1-p)*exp(1/p)) - (1 - A) * log(p*exp(1/(1-p))))
}

L <- function(theta, X, A) {
  p <- plogis(drop(X %*% theta))

  sum(A * log(p/((1-p)*exp(1/p))) + (1 - A) * log(1 - p) - A * log((1-p)*exp(1/p)) - (1 - A) * log(p*exp(1/(1-p))))
}

##Logistic regression
L <- function(theta, X, A) {
  p <- plogis(drop(X %*% theta))

  sum(A * log(p) + (1 - A) * log(1 - p))
}

fit <- glm(treat ~ age + educ + married, data = lalonde, family = binomial)

opt <- optim(fn = L, par = fit$coefficients,
             X = model.matrix(fit), A = fit$y,
             control = list(fnscale = -1, reltol = 1e-16,
                            maxit = 1e5))

L(opt$par, X = model.matrix(fit), A = fit$y) -
L(coef(fit), X = model.matrix(fit), A = fit$y)

w <- get_w_from_ps(plogis(drop(model.matrix(fit) %*% opt$par)),
                   treat = fit$y, estimand = "ATE")

bal.tab(treat ~ age + educ + married, data = lalonde, weights = w)

W <- weightit(treat ~ age + educ + married, data = lalonde, method = "cbps")

bal.tab(treat ~ age + educ + married, data = lalonde, weights = list(man = w, cbps = W), pairwise =F)

#Balance representation
W <- function(p, t, family, estimand = "ATE") {
  linkp <- function(p) family$mu.eta(family$linkfun(p))
  switch(estimand,
         "ATE" = (t * (1 - p) / linkp(p)) +
           (1 - t) * p / linkp(p),
         "ATT" = p * ((t * (1 - p) / linkp(p)) +
                        (1 - t) * p / linkp(p)),
         "ATC" = (1 - p) * ((t * (1 - p) / linkp(p)) +
                              (1 - t) * p / linkp(p)),
         "ATO" = p * (1 - p) * ((t * (1 - p) / linkp(p)) +
                                  (1 - t) * p / linkp(p)))

}

S <- function(theta, X, A, family, estimand) {
  p <- family$linkinv(drop(X %*% theta))
  w <- W(p, A, family, estimand)

  colMeans((A - (1 - A)) * w * X)
}

estimand <- "ATE"
link <- "cloglog"

library(rootSolve)
out <- multiroot(f = S,
                 start = fit$coefficients,
                 X = model.matrix(fit),
                 A = fit$y,
                 family = binomial(link),
                 estimand = estimand,
                 rtol = 1e-10)

w <- do.call(W, list(p = binomial(link)$linkinv(drop(model.matrix(fit) %*% out$root)),
                     t = fit$y, family = binomial(link),
                     estimand = estimand))
w <- get_w_from_ps(binomial(link)$linkinv(drop(model.matrix(fit) %*% out$root)),
                   treat = fit$y, estimand = estimand)

bal.tab(treat ~ age + educ + married, data = lalonde, weights = w, estimand = estimand)
