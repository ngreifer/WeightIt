#Test
source('~/Dropbox (Personal)/Research/R/WeightIt/R/weightit2method.R')
source('~/Dropbox (Personal)/Research/R/WeightIt/R/functions_for_processing.R')
source('~/Dropbox (Personal)/Research/R/WeightIt/R/weightit.R')
#source('~/Dropbox (Personal)/Research/R/WeightIt/R/weightitMSM.R')
source('~/Dropbox (Personal)/Research/R/WeightIt/R/weightitMSM2.R')
source('~/Dropbox (Personal)/Research/R/WeightIt/R/weightit.fit.R')

#Tests things quickly
library("cobalt")
data("lalonde", package = "cobalt")
covs <- subset(lalonde, select = -c(re78, treat))
lalonde$treat3 <- factor(ifelse(lalonde$treat == 1, "A", sample(c("B", "C"), nrow(lalonde), T)))
lalonde$treat5 <- factor(sample(c(LETTERS[1:5]), nrow(lalonde), T))

s <- runif(nrow(lalonde), 0, 2)

#method = "ps"
W <- weightit(treat ~ covs, data = lalonde, method = "ps", estimand = "ATE")
W <- weightit(lalonde$treat ~ covs, method = "ps", estimand = "ATT", link = "probit")
W <- weightit(f.build("treat", covs), data = lalonde, method = "ps", estimand = "ATC", s.weights = s)
W <- weightit(f.build("treat", covs), data = lalonde, method = "ps", estimand = "ATO", stabilize = T)
W <- weightit(f.build("treat", covs), data = lalonde, method = "ps", estimand = "ATM")

W <- weightit(f.build("treat3", covs), data = lalonde, method = "ps", estimand = "ATE",
              link = "bayes.probit")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "ps", estimand = "ATT",
              focal = "A", s.weights = s)

W <- weightit(f.build("re78", covs), data = lalonde, method = "ps")
W <- weightit(f.build("re78", covs), data = lalonde, method = "ps",
              num.formula = ~ age + educ + race + married)

#method = "gbm"
W <- weightit(f.build("treat", covs), data = lalonde, method = "gbm",
              stop.method = c("es.mean", "ks.mean"), estimand = "ATE")
W <- weightit(f.build("treat", covs), data = lalonde, method = "twang", estimand = "ATT")
W <- weightit(f.build("treat", covs), data = lalonde, method = "gbr", estimand = "ATC", s.weights = s)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "gbm", estimand = "ATE")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "gbm", estimand = "ATT",
              focal = "A", s.weights = s)

#method = "cbps"
W <- weightit(f.build("treat", covs), data = lalonde, method = "cbps", estimand = "ATE", over = FALSE)
W <- weightit(f.build("treat", covs), data = lalonde, method = "cbps", estimand = "ATT")
W <- weightit(f.build("treat", covs), data = lalonde, method = "cbps", estimand = "ATC", s.weights = s)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "cbps", estimand = "ATE", s.weights = s)
W <- weightit(f.build("treat3", covs), data = lalonde, method = "cbps", estimand = "ATT")
W <- weightit(f.build("treat5", covs), data = lalonde, method = "cbps", estimand = "ATE")

W <- weightit(f.build("re78", covs), data = lalonde, method = "cbps", over = FALSE)

#method = "npcbps"
W <- weightit(f.build("treat", covs), data = lalonde, method = "npcbps", estimand = "ATE")
W <- weightit(f.build("treat", covs), data = lalonde, method = "npcbps", estimand = "ATT")
W <- weightit(f.build("treat", covs), data = lalonde, method = "npcbps", estimand = "ATC", s.weights = s)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "npcbps", estimand = "ATE")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "npcbps", estimand = "ATT")

W <- weightit(f.build("re78", covs), data = lalonde, method = "npcbps")

#method = "ebal"
W <- weightit(f.build("treat", covs), data = lalonde, method = "ebal", estimand = "ATE")
W <- weightit(f.build("treat", covs), data = lalonde, method = "ebal", estimand = "ATC",
              stabilize = TRUE)
W <- weightit(f.build("treat", covs), data = lalonde, method = "entropy", estimand = "ATT",
              s.weights = s)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "ebal", estimand = "ATE")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "ebal", estimand = "ATT",
              focal = "A", s.weights = s)

#method = "ebcw"
W <- weightit(f.build("treat", covs), data = lalonde, method = "ebcw", estimand = "ATE")
W <- weightit(f.build("treat", covs), data = lalonde, method = "ebcw", estimand = "ATC")
W <- weightit(f.build("treat", covs), data = lalonde, method = "ebcw", estimand = "ATT",
              s.weights = s)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "ebcw", estimand = "ATE")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "ebcw", estimand = "ATT",
              focal = "A", s.weights = s)

m <- model.matrix(formula, data=lalonde,
                  contrasts.arg = lapply(lalonde[sapply(lalonde, is.factor)],
                                         contrasts, contrasts=FALSE))

est.fun <- function(data, which.sample, stabilize) {
  W <- weightit(f.build("treat", covs), data = data[which.sample,], method = "ps", stabilize = stabilize)
  est <- coef(lm(re78 ~ treat, data = data[which.sample,], weights = get.w(W)))["treat"]
  return(est)
}

library("twang")
data(iptwExWide)
iptwExWide$tx2_ <- factor(ifelse(iptwExWide$tx2 == 0, 0, sample(1:2, 1000, replace = T)))
psmsm <- iptw(list(tx1 ~ use0 + gender + age,
                   tx2 ~ use1 + use0 + tx1 + gender + age,
                   tx3 ~ use2 + use1 + use0 + tx2 + tx1 + gender + age),
              timeInvariant = ~ gender + age,
              data = iptwExWide,
              cumulative = FALSE,
              priorTreatment = FALSE,
              verbose = FALSE,
              stop.method = "es.max",
              n.trees = 200)

Wmsm <- weightitMSM(list(tx1 ~ use0 + gender + age,
                         tx2 ~ use1 + use0 + tx1 + gender + age,
                         tx3 ~ use2 + use1 + use0 + tx2 + tx1 + gender + age),
                    data = iptwExWide,
                    verbose = FALSE,
                    method = "ebal")

data("iptwExLong")
pre <- iptwExLong$covariates
wide <- reshape(pre, timevar = "time",
                v.names = c("use", "tx"),
                idvar = c("ID"),
                direction = "wide",
                sep = "")

library("cobalt")
data("lalonde_mis")
covs_mis <- subset(lalonde_mis, select = -c(re78, treat))

W <- weightit(f.build("treat", covs_mis), data = lalonde_mis, estimand = "ATE",
              method = "ps")
