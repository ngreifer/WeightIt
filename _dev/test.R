#Test
# for (i in dir("R/")) source(paste0("R/", i))
# library(ggplot2); library(crayon)
# stop("Done sourcing.", call. = FALSE)

devtools::load_all(".")

#Tests things quickly
library("cobalt")
data("lalonde", package = "cobalt")
data("lalonde_mis", package = "cobalt")
covs <- subset(lalonde, select = -c(re78, treat))
lalonde$treat3 <- factor(ifelse(lalonde$treat == 1, "A", sample(c("B", "C"), nrow(lalonde), T)), ordered = F)
lalonde$treat5 <- factor(sample(c(LETTERS[1:5]), nrow(lalonde), T))

s <- optweight::optweight.svy(~age + educ + married + re74, data = lalonde,
                              targets = c(age = 26, educ = 12, married = .5, re74 = 4000),
                              min.w = 0)$weights

W <- weightit(treat ~ covs, data = lalonde, method = "bart", estimand = "ATT")
for (i in 1:30) {cat(i, "| "); W <- weightit(treat ~ covs, data = lalonde, method = "bart", estimand = "ATT")}

#method = "glm"
W <- weightit(treat ~ covs, data = lalonde, method = "glm", estimand = "ATE")
W <- weightit(lalonde$treat ~ covs, method = "glm", estimand = "ATT", link = "probit")
W <- weightit(f.build("treat", covs), data = lalonde, method = "glm", estimand = "ATC", s.weights = s)
W <- weightit(f.build("treat", covs), data = lalonde, method = "ps", estimand = "ATO", stabilize = T,
              link = "br.logit")
W <- weightit(f.build("treat", covs), data = lalonde, method = "glm", estimand = "ATM")

W <- weightit(f.build("treat3", covs), data = lalonde, method = "ps", estimand = "ATE",
              link = "bayes.probit")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "ps", estimand = "ATT",
              focal = "A", s.weights = s)

W <- weightit(f.build("treat5", covs), data = lalonde, method = "ps", estimand = "ATE",
              link = "logit", focal = "A")

W <- weightit(f.build("re78", covs), data = lalonde, method = "glm")
W <- weightit(f.build("re78", covs), data = lalonde, method = "ps", use.kernel = T, plot = T)

#method = "gbm"
W <- weightit(f.build("treat", covs), data = lalonde, method = "gbm",
              criterion = c("es.mean", "ks.mean"), estimand = "ATE")
W <- weightit(f.build("treat", covs), data = lalonde, method = "gbm",
              stop.method = c("es.mean"), estimand = "ATE")
W <- weightit(f.build("treat", covs), data = lalonde, method = "gbm",
              criterion = c("smd.mean"), estimand = "ATT",
              subclass = 40)
W <- weightit(f.build("treat", covs), data = lalonde_mis, method = "gbm", estimand = "ATT", missing = "surr",
              criterion = "r2.2")
W <- weightit(f.build("treat", covs), data = lalonde, method = "gbr", estimand = "ATC", s.weights = s,
              criterion = "ks.rms")

W <- weightit(f.build("treat3", covs), data = lalonde, method = "gbm", estimand = "ATE",
              criterion = "energy.dist")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "gbm", estimand = "ATT",
              focal = "A", s.weights = s, criterion = "ks.max")

W <- weightit(f.build("re78", covs), data = lalonde, method = "gbm", criterion = "p.max")
W <- weightit(f.build("re78", covs), data = lalonde, method = "gbm", criterion = "distance.cov")
W <- weightit(f.build("re78", covs), data = lalonde_mis, method = "gbm", criterion = "p.max", use.kernel = TRUE)


#method = "cbps"
W <- weightit(f.build("treat", covs), data = lalonde, method = "cbps", estimand = "ATE", over = FALSE)
W <- weightit(f.build("treat", covs), data = lalonde, method = "cbps", estimand = "ATT")
W <- weightit(f.build("treat", covs), data = lalonde_mis, method = "cbps", estimand = "ATC", s.weights = s)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "cbps", estimand = "ATE", s.weights = s)
W <- weightit(f.build("treat3", covs), data = lalonde, method = "cbps", estimand = "ATT")
W <- weightit(f.build("treat5", covs), data = lalonde, method = "cbps", estimand = "ATE")

W <- weightit(f.build("re78", covs), data = lalonde, method = "cbps", over = FALSE, s.weights = s)

#method = "npcbps"
W <- weightit(f.build("treat", covs), data = lalonde, method = "npcbps", estimand = "ATE")
W <- weightit(f.build("treat", covs), data = lalonde, method = "npcbps", estimand = "ATT")
W <- weightit(f.build("treat", covs), data = lalonde, method = "npcbps", estimand = "ATE", s.weights = s)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "npcbps", estimand = "ATE")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "npcbps", estimand = "ATT")

W <- weightit(f.build("re78", covs), data = lalonde, method = "npcbps")

#method = "ebal"
W <- weightit(f.build("treat", covs), data = lalonde, method = "ebal", estimand = "ATE", moments = 2)
W <- weightit(f.build("treat", covs), data = lalonde, method = "ebal", estimand = "ATC",
              stabilize = TRUE)
W <- weightit(f.build("treat", covs), data = lalonde, method = "entropy", estimand = "ATT",
              s.weights = s, int = T)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "ebal", estimand = "ATE")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "ebal", estimand = "ATT",
              focal = "A", s.weights = s)

W <- weightit(f.build("re78", covs), data = lalonde, method = "ebal")
W <- weightit(f.build("re78", covs), data = lalonde, method = "ebal", moments = 2, d.moments = 3)


# #method = "ebcw"
# W <- weightit(f.build("treat", covs), data = lalonde, method = "ebcw", estimand = "ATE")
# W <- weightit(f.build("treat", covs), data = lalonde, method = "ebcw", estimand = "ATC")
# W <- weightit(f.build("treat", covs), data = lalonde, method = "ebcw", estimand = "ATT",
#               s.weights = s)
#
# W <- weightit(f.build("treat3", covs), data = lalonde, method = "ebcw", estimand = "ATE")
# W <- weightit(f.build("treat3", covs), data = lalonde, method = "ebcw", estimand = "ATT",
#               focal = "A", s.weights = s)

#method = "optweight"
W <- weightit(f.build("treat", covs), data = lalonde, method = "optweight", estimand = "ATE")
W <- weightit(f.build("treat", covs), data = lalonde, method = "optweight", estimand = "ATC", tols = .1)
W <- weightit(f.build("treat", covs), data = lalonde, method = "optweight", estimand = "ATT",
              s.weights = s, tols = .1)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "optweight", estimand = "ATE")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "optweight", estimand = "ATT",
              focal = "A", s.weights = s)

W <- weightit(f.build("re78", covs), data = lalonde, method = "optweight")

#method = "super"
W <- weightit(f.build("treat", covs), data = lalonde, method = "super", estimand = "ATE",
              SL.library = c("SL.glm", "SL.lda"))
W <- weightit(f.build("treat", covs), data = lalonde, method = "super", estimand = "ATC",
              SL.library = c("SL.glm", "SL.lda"), SL.method = "method.balance",
              criterion = "ks.max")
W <- weightit(f.build("treat", covs), data = lalonde, method = "super", estimand = "ATT",
              SL.library = c("SL.glm", "SL.lda", "SL.gbm"), SL.method = "method.balance",
              stop.method = "es.max")
W <- weightit(f.build("treat", covs), data = lalonde, method = "super", estimand = "ATO",
              s.weights = s, SL.library = c("SL.glm", "SL.lda"))

W <- weightit(f.build("treat3", covs), data = lalonde, method = "super", estimand = "ATE",
              SL.library = c("SL.glm", "SL.lda"))
W <- weightit(f.build("treat3", covs), data = lalonde, method = "super", estimand = "ATT",
              focal = "A", s.weights = s, SL.library = c("SL.glm", "SL.lda"))

W <- weightit(f.build("re78", covs), data = lalonde, method = "super",
              SL.library = c("SL.glm", "SL.stepAIC"), use.kernel = TRUE)
W <- weightit(f.build("re78", covs), data = lalonde, method = "super",
              SL.library = c("SL.glm", "SL.stepAIC"), SL.method = "method.balance",
              criterion = "p.max", density = "dt_3", s.weights = s)

#method = "bart"
W <- weightit(treat ~ covs, data = lalonde, method = "bart", estimand = "ATE")
W <- weightit(lalonde$treat ~ covs, method = "bart", estimand = "ATT")
all.equal(weightit(lalonde$treat ~ covs, method = "bart", estimand = "ATT", rngSeed = 100, n.threads = 1, n.trees = 20),
          weightit(lalonde$treat ~ covs, method = "bart", estimand = "ATT", rngSeed = 100, n.threads = 1, n.trees = 20))
W <- weightit(f.build("treat", covs), data = lalonde, method = "bart", estimand = "ATC", s.weights = s)
W <- weightit(f.build("treat", covs), data = lalonde, method = "bart", estimand = "ATO")

W <- weightit(f.build("treat3", covs), data = lalonde, method = "bart", estimand = "ATE")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "bart", estimand = "ATT",
              focal = "A")

W <- weightit(f.build("treat5", covs), data = lalonde, method = "bart", estimand = "ATE",
              focal = "A")

W <- weightit(f.build("re78", covs), data = lalonde, method = "bart", density = "dt_3")
W <- weightit(f.build("re78", covs), data = lalonde, method = "bart", use.kernel = T, plot = T)
all.equal(weightit(f.build("re78", covs), data = lalonde, method = "bart", rngSeed = 100, n.threads = 4, n.trees = 20),
          weightit(f.build("re78", covs), data = lalonde, method = "bart", rngSeed = 100, n.threads = 4, n.trees = 20))

#method = "energy"
W <- weightit(f.build("treat", covs), data = lalonde, method = "energy", estimand = "ATE")
W <- weightit(f.build("treat", covs), data = lalonde, method = "energy", estimand = "ATE", dist.mat = "mahal")
W <- weightit(f.build("treat", covs), data = lalonde, method = "energy", estimand = "ATE", improved = FALSE)
W <- weightit(f.build("treat", covs), data = lalonde, method = "energy", estimand = "ATC", moments = 1)
W <- weightit(f.build("treat", covs), data = lalonde, method = "energy", estimand = "ATT",
              s.weights = s, moments = 1)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "energy", estimand = "ATE")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "energy", estimand = "ATT",
              focal = "A", s.weights = s)

W <- weightit(f.build("re78", covs), data = lalonde, method = "energy")
W <- weightit(f.build("re78", covs), data = lalonde, method = "energy", moments = 1)

#user defined

myfun <- function(treat, covs, estimand, r = FALSE) {
  d <- data.frame(treat, covs)
  f <- formula(d)
  ps <- glm(f, data = d, family = "binomial")$fitted
  if (estimand == "ATT") w <- ifelse(treat == 1, 1, ps/(1-ps))
  else if (estimand == "ATE") w <- ifelse(treat == 1, 1/ps, 1/(1-ps))

  if (r) w[1] <- 100
  return(list(w = w, ps = ps))
}
W <- weightit(f.build("treat", covs), data = lalonde, method = myfun, estimand = "ATE")

#as.weightit
W <- weightit(treat ~ covs, data = lalonde, method = "ps", estimand = "ATE")
W_ <- as.weightit(W$weights, lalonde$treat, covs = covs, estimand = "ATE", s.weights = NULL, ps = NULL)
summary(W_)
bal.tab(W_)

#rms in formula
library(rms)
W <- weightit(treat ~ rcs(age) + poly(educ,2) + catg(race) + married + nodegree + re74 + re75,
              data = lalonde)

library("twang")
data(iptwExWide, package = "twang")
iptwExWide$tx2_ <- factor(ifelse(iptwExWide$tx2 == 0, 0, sample(1:2, 1000, replace = T)))
psmsm <- iptw(list(tx1 ~ use0 + gender + age,
                   tx2 ~ use1 + use0 + tx1 + gender + age,
                   tx3 ~ use2 + use1 + use0 + tx2 + tx1 + gender + age),
              timeInvariant = ~ gender + age,
              data = iptwExWide,
              cumulative = FALSE,
              priorTreatment = FALSE,
              verbose = FALSE,
              criterion = "es.max",
              n.trees = 200)

Wmsm <- weightitMSM(list(tx1 ~ use0 + gender + age,
                         tx2_ ~ use1 + use0 + tx1 + gender + age,
                         tx3 ~ use2 + use1 + use0 + tx2 + tx1 + gender + age),
                    data = iptwExWide,
                    verbose = FALSE,
                    method = "ps", stabilize = TRUE)
data(mnIptwExWide, package = "twang")
Wmsm <- weightitMSM(list(tx1 ~ use0 + gender + age,
                         tx2 ~ use1 + use0 + tx1 + gender + age,
                         tx3 ~ use2 + use1 + use0 + tx2 + tx1 + gender + age),
                    data = mnIptwExWide,
                    verbose = FALSE,
                    method = "ps", link = "logit")

data("iptwExLong")
pre <- iptwExLong$covariates
wide <- reshape(pre, timevar = "time",
                v.names = c("use", "tx"),
                idvar = c("ID"),
                direction = "wide",
                sep = "")
long <- reshape(iptwExWide, timevar = "time",
                        v.names = c("use", "tx"),
                        idvar = c("ID"),
                        direction = "long",
                        sep = "",
                varying = list(use = c("use0", "use1", "use2"),
                               tx = c("tx1", "tx2", "tx3")))

library("cobalt")
data("lalonde_mis")
covs_mis <- subset(lalonde_mis, select = -c(re78, treat))

W <- weightit(f.build("treat", covs_mis), data = lalonde_mis, estimand = "ATE",
              method = "ps")

k <- 0;
at <- .95
while(k<.1) {
  at <- at-.01
  k <- max(abs(bal.tab(trim(W, at))$Balance$Corr.Adj[-1]))
}

library("KBAL")
covs_ <- splitfactor(covs, drop.first = F)
invisible(capture.output(k <- kbal(covs_,
          lalonde$treat,
          whiten = F,
          method = "ebal")))
bal.tab(lalonde[-c(1,9)],
        lalonde$treat,
        weights = k$w)

kbal.fun <- function(treat, covs, estimand, focal, ...) {
  args <- list(...)

  if (is_not_null(focal)) treat <- as.numeric(treat == focal)
  else stop("estimand must be 'ATT' or 'ATC'.", call. = FALSE)

  if ("kbal.method" %in% names(args)) {
    names(args)[names(args) == "kbal.method"] <- "method"
  }

  args[names(args) %nin% setdiff(names(formals(KBAL::kbal)), c("X", "D"))] <- NULL
  k.out <- do.call(KBAL::kbal, c(list(X = covs, D = treat), args))

  w <- k.out$w
  return(list(w = w))
}
W <- weightit(f.build("treat", covs), data = lalonde, method = kbal.fun, estimand = "ATT")

data("lalonde_mis")
library(mice); library(survey)
imp <- mice(lalonde_mis, 4)

imp.data <- complete(imp, "long")
w <- weightit(treat ~ age + educ + race + married + nodegree + re74 + re75,
              data = imp.data, method = "ps", estimand = "ATE",
              by = ".imp")

fit.list <- lapply(unique(imp.data$.imp), function(i) {
  svyglm(re78 ~ treat, design = svydesign(~1, weights = w$weights[imp.data$.imp == i],
                                          data = imp.data[imp.data$.imp == i,]))
})

summary(pool(as.mira(fit.list)))

library(MatchThem)
w.imp <- weightthem(treat ~ age + educ + race + married + nodegree + re74 + re75,
                   imp, method = "ps", estimand = "ATE")
summary(pool(with(w.imp, glm(re78 ~ treat))))

#testing get_w_from_ps
data("lalonde", package = "cobalt")
t01 <- lalonde$treat
tf01 <- factor(lalonde$treat)
tfab <- factor(ifelse(lalonde$treat == 1, "a", "b"))

p <- glm(treat ~ age + educ + married + race, data = lalonde, family = binomial)$fitted

wate <- t01/p + (1-t01)/(1-p)
watt <- wate*p
watc <- wate*(1-p)

all.equal(wate, get_w_from_ps(p, t01, "ATE"))
all.equal(wate, get_w_from_ps(p, tf01, "ATE"))
all.equal(wate, get_w_from_ps(p, tfab, "ATE", treated = "a"))

all.equal(watt, get_w_from_ps(p, t01, "ATT"))
all.equal(watt, get_w_from_ps(p, tf01, "ATT", focal = 1))
all.equal(watt, get_w_from_ps(p, tfab, "ATT", focal = "a", treated = "a"))

all.equal(watc, get_w_from_ps(p, t01, "ATC"))
all.equal(watc, get_w_from_ps(p, tf01, "ATC"))
all.equal(watc, get_w_from_ps(p, tfab, "ATC"))