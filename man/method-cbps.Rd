\name{method_cbps}
\alias{method_cbps}

\title{Covariate Balancing Propensity Score Weighting}

\description{
This page explains the details of estimating weights from covariate balancing propensity scores by setting \code{method = "cbps"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating propensity scores using generalized method of moments and then converting those propensity scores into weights using a formula that depends on the desired estimand. This method relies on \code{\link[CBPS]{CBPS}} from the \pkg{CBPS} package.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the propensity scores and weights using \code{\link[CBPS]{CBPS}}. The following estimands are allowed: ATE, ATT, and ATC. The weights are taken from the output of the \code{CBPS} fit object. When the estimand is the ATE, the return propensity score is the probability of being in the "second" treatment group, i.e., \code{levels(factor(treat))[2]}; when the estimand is the ATC, the returned propensity score is the probability of being in the control (i.e., non-focal) group. When \code{include.obj = TRUE}, the returned object is the \code{CBPS} fit.
}
\subsection{Multinomial Treatments}{
For multinomial treatments with three or four categories and when the estimand is the ATE, this method estimates the propensity scores and weights using one call to \code{\link[CBPS]{CBPS}}. For multinomial treatments with three or four categories or when the estimand is the ATT, this method estimates the propensity scores and weights using multiple calls to \code{\link[CBPS]{CBPS}}. The following estimands are allowed: ATE and ATT. The weights are taken from the output of the \code{CBPS} fit objects. When \code{include.obj = TRUE}, the returned object is either a single \code{CBPS} fit object if only one is used or a list of the \code{CBPS} fits if multiple are used.
}
\subsection{Continuous Treatments}{
For continuous treatments, the generalized propensity score and weights are estimated using \code{\link[CBPS]{CBPS}}. When \code{include.obj = TRUE}, the returned object is the \code{CBPS} fit.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point. This is not how \code{\link[CBPS]{CBMSM}} in the \pkg{CBPS} package estimates weights for longitudinal treatments.
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios. See Note about sampling weights.
}
\subsection{Missing Data}{
Missing data is not compatible with the CBPS algorithm, so a few extra things happen when \code{NA}s are present in the covariates. First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is \code{NA} and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with 0s (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting \code{weightit} object will be the original covariates with the \code{NA}s.
}
}
\section{Additional Arguments}{
All arguments to \code{CBPS} can be passed through \code{weightit} or \code{weightitMSM}, with the following exceptions:

\code{method} in \code{CBPS} is replaced with the argument \code{over} in \code{weightit}. Setting \code{over = FALSE} in \code{weightit} is the equivalent of setting \code{method = "exact"} in \code{CBPS}.

\code{sample.weights} is ignored because sampling weights are passed using \code{s.weights}.

\code{standardize} is ignored.

All arguments take on the defaults of those in \code{CBPS}. It may be useful in many cases to set \code{over = FALSE}, especially with continuous treatments.

}
\references{
\strong{Binary treatments}

Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(1), 243–263.

\strong{Multinomial Treatments}

Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(1), 243–263.


\strong{Continuous treatments}

Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity score for a continuous treatment: Application to the efficacy of political advertisements. The Annals of Applied Statistics, 12(1), 156–177. \doi{10.1214/17-AOAS1101}


}
\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}
}
\note{
When sampling weights are used with \code{CBPS}, the estimated weights already incorporate the sampling weights. When \code{weightit} is used with \code{method = "cbps"}, the estimated weights are separated from the sampling weights, as they are with all other methods.
}
\examples{\dontrun{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cbps", estimand = "ATT"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cbps", estimand = "ATE"))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cbps", over = FALSE))
summary(W3)
bal.tab(W3)
}}