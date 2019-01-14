#Test
for (i in dir("R/")) source(paste0("R/", i))


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
W <- weightit(f.build("re78", covs), data = lalonde, method = "ps", use.kernel = T, plot = T)

#method = "gbm"
W <- weightit(f.build("treat", covs), data = lalonde, method = "gbm",
              stop.method = c("es.mean", "ks.mean"), estimand = "ATE")
W <- weightit(f.build("treat", covs), data = lalonde, method = "twang", estimand = "ATT")
W <- weightit(f.build("treat", covs), data = lalonde, method = "gbr", estimand = "ATC", s.weights = s)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "gbm", estimand = "ATE")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "gbm", estimand = "ATT",
              focal = "A", s.weights = s, stop.method = "es.blarg")

W <- weightit(f.build("re78", covs), data = lalonde, method = "gbm", stop.method = "p.max", use.optimize = 2)
W <- weightit(f.build("re78", covs), data = lalonde, method = "gbm", stop.method = "p.max", use.kernel = TRUE)


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
W <- weightit(f.build("treat", covs), data = lalonde, method = "ebal", estimand = "ATE", moments = 2)
W <- weightit(f.build("treat", covs), data = lalonde, method = "ebal", estimand = "ATC",
              stabilize = TRUE)
W <- weightit(f.build("treat", covs), data = lalonde, method = "entropy", estimand = "ATT",
              s.weights = s, int = T)

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

#method = "optweight"
W <- weightit(f.build("treat", covs), data = lalonde, method = "optweight", estimand = "ATE")
W <- weightit(f.build("treat", covs), data = lalonde, method = "optweight", estimand = "ATC", tols = .1)
W <- weightit(f.build("treat", covs), data = lalonde, method = "optweight", estimand = "ATT",
              s.weights = s, tols = .1)

W <- weightit(f.build("treat3", covs), data = lalonde, method = "optweight", estimand = "ATE")
W <- weightit(f.build("treat3", covs), data = lalonde, method = "optweight", estimand = "ATT",
              focal = "A", s.weights = s)

W <- weightit(f.build("re78", covs), data = lalonde, method = "optweight")

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
                    method = "ps")
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

k <- kbal(lalonde[-c(1,4,9)],
          lalonde$treat,
          whiten = F,
          method = "ebal")
bal.tab(lalonde[-c(1,4,9)],
        lalonde$treat,
        weights = k$w)
