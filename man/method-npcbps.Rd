\name{method_npcbps}
\alias{method_npcbps}

\title{Nonparametric Covariate Balancing Propensity Score Weighting}

\description{
This page explains the details of estimating weights from nonparametric covariate balancing propensity scores by setting \code{method = "npcbps"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating weights by maximizing the empirical likelihood of the data subject to balance constraints. This method relies on \code{\link[CBPS]{npCBPS}} from the \pkg{CBPS} package.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the weights using \code{\link[CBPS]{npCBPS}}. The ATE is the only estimand allowed. The weights are taken from the output of the \code{npCBPS} fit object. When \code{include.obj = TRUE}, the returned object is the \code{npCBPS} fit.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, this method estimates the weights using \code{\link[CBPS]{npCBPS}}. The ATE is the only estimand allowed. The weights are taken from the output of the \code{npCBPS} fit object. When \code{include.obj = TRUE}, the returned object is the \code{npCBPS} fit.
}
\subsection{Continuous Treatments}{
For continuous treatments, this method estimates the weights using \code{\link[CBPS]{npCBPS}}. The weights are taken from the output of the \code{npCBPS} fit object. When \code{include.obj = TRUE}, the returned object is the \code{npCBPS} fit.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point. This is not how \code{\link[CBPS]{CBMSM}} in the \pkg{CBPS} package estimates weights for longitudinal treatments.
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios.
}
\subsection{Missing Data}{
Missing data is not compatible with the nonparametric CBPS algorithm, so a few extra things happen when \code{NA}s are present in the covariates. First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is \code{NA} and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with 0s (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting \code{weightit} object will be the original covariates with the \code{NA}s.
}
}
\section{Additional Arguments}{
All argument to \code{npCBPS} can be passed through \code{weightit} or \code{weightitMSM}.

All arguments take on the defaults of those in \code{npCBPS}.

}
\references{

Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity score for a continuous treatment: Application to the efficacy of political advertisements. The Annals of Applied Statistics, 12(1), 156–177. \doi{10.1214/17-AOAS1101}

}
\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}
}

\examples{
# Examples take a long time to run
\dontrun{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "npcbps", estimand = "ATE"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "npcbps", estimand = "ATE"))
summary(W2)
bal.tab(W2)
}
}