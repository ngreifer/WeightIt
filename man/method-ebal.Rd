\name{method_ebal}
\alias{method_ebal}

\title{Entropy Balancing}

\description{
This page explains the details of estimating weights using entropy balancing by setting \code{method = "ebal"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary and multinomial treatments.

In general, this method relies on estimating weights by minimizing the entropy of the weights subject to exact moment balancing constraints. This method relies on \code{\link[ebal]{ebalance}} from the \pkg{ebal} package.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the weights using \code{\link[ebal]{ebalance}}. The following estimands are allowed: ATE, ATT, and ATC. The weights are taken from the output of the \code{ebalance} fit object. When the ATE is requested, \code{ebalance} is run twice, once for each treatment group. When \code{include.obj = TRUE}, the returned object is the \code{ebal} fit (or a list of the two fits when the estimand is the ATE).
}
\subsection{Multinomial Treatments}{
For multinomial treatments, this method estimates the weights using \code{\link[ebal]{ebalance}}. The following estimands are allowed: ATE and ATT. The weights are taken from the output of the \code{ebalance} fit objects. When the ATE is requested, \code{ebalance} is run once for each treatment group. When the ATT is requested, \code{ebalance} is run once for each non-focal (i.e., control) group. When \code{include.obj = TRUE}, the returned object is the list of \code{ebalance} fits.
}
\subsection{Continuous Treatments}{
Continuous treatments are not supported.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point. This method is not guaranteed to yield exact balance at each time point. NOTE: the use of entropy balancing with longitudinal treatments has not been validated!
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios.
}
\subsection{Missing Data}{
Missing data is not compatible with the entropy balancing algorithm, so a few extra things happen when \code{NA}s are present in the covariates. First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is \code{NA} and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with 0s (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting \code{weightit} object will be the original covariates with the \code{NA}s.
}
}
\section{Additional Arguments}{
All argument to \code{ebalance} can be passed through \code{weightit} or \code{weightitMSM}, with the following exceptions:

\code{base.weight} is ignored because sampling weights are passed using \code{s.weights}.

All arguments take on the defaults of those in \code{ebalance}.

When \code{standardize = TRUE} in the call to \code{weightit}, \code{\link[ebal]{ebalance.trim}} is run on the resulting \code{ebalance} fit objects. Doing so can reduce the variability of the weights while maintaining covariate balance.

}
\references{

Hainmueller, J. (2012). Entropy Balancing for Causal Effects: A Multivariate Reweighting Method to Produce Balanced Samples in Observational Studies. Political Analysis, 20(1), 25–46. \doi{10.1093/pan/mpr025}

}
\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}
}

\examples{\dontrun{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATT"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATE",
                standardize = TRUE))
summary(W2)
bal.tab(W2)

}}