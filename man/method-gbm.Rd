\name{method_gbm}
\alias{method_gbm}

\title{Propensity Score Weighting Using Generalized Boosted Models}

\description{
This page explains the details of estimating weights from generalized boosted model-based propensity scores by setting \code{method = "gbm"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating propensity scores using generalized boosted modeling and then converting those propensity scores into weights using a formula that depends on the desired estimand. For binary and multinomial treatments, this method relies on \code{\link[twang]{ps}} from the \pkg{twang} package. For continuous treatments, this method relies on \code{\link{ps.cont}} from \pkg{WeightIt}. The algorithm involves choosing a balance criterion to optimize so that balance, rather than prediction, is prioritized.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the propensity scores using \code{\link[twang]{ps}}. The following estimands are allowed: ATE, ATT, and ATC. The weights for the ATE, ATT, and ATC are computed from the estimated propensity scores using the standard formulas. When the estimand is the ATE, the return propensity score is the probability of being in the "second" treatment group, i.e., \code{levels(factor(treat))[2]}; when the estimand is the ATC, the returned propensity score is the probability of being in the control (i.e., non-focal) group. When \code{include.obj = TRUE}, the returned object is the \code{ps} fit.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, this method estimates the propensity scores using multiple calls to \code{\link[twang]{ps}} (which is the same thing that \code{\link[twang]{mnps}} in \pkg{twang} does behind the scenes). The following estimands are allowed: ATE and ATT. The weights are computed from the estimated propensity scores using the standard formulas. When \code{include.obj = TRUE}, the returned object is a list of the \code{ps} fits.
}
\subsection{Continuous Treatments}{
For continuous treatments, the generalized propensity score is estimated using \code{\link{ps.cont}}. When \code{include.obj = TRUE}, the returned object is the \code{ps.cont} fit.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point, which is what \code{\link[twang]{iptw}} in \pkg{twang} does behind the scenes.
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios.
}
}
\section{Additional Arguments}{
Any argument to \code{ps} or \code{ps.cont} can be passed through \code{weightit} or \code{weightitMSM}. The argument \code{stop.method} is required, and, in contrast to \code{ps} and \code{ps.cont}, can only be of length 1. If not provided, a default will be used (\code{"es.mean"} for binary and multinomial treatments and \code{"s.mean.z"} for continuous treatments).

\code{sampw} is ignored because sampling weights are passed using \code{s.weights}.

All other arguments take on the defaults of those in \code{ps} and \code{ps.cont}. The main arguments of interest are \code{n.trees} and \code{interaction.depth}, and increasing them improves performance.

}
\references{
\strong{Binary treatments}

McCaffrey, D. F., Ridgeway, G., & Morral, A. R. (2004). Propensity Score Estimation With Boosted Regression for Evaluating Causal Effects in Observational Studies. Psychological Methods, 9(4), 403–425. \doi{10.1037/1082-989X.9.4.403}

\strong{Multinomial Treatments}

McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. Statistics in Medicine, 32(19), 3388–3414. \doi{10.1002/sim.5753}

\strong{Continuous treatments}

Zhu, Y., Coffman, D. L., & Ghosh, D. (2015). A Boosting Algorithm for Estimating Generalized Propensity Scores with Continuous Treatments. Journal of Causal Inference, 3(1).

}
\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}
}
\note{
The weights may be very slightly different from those estimated by a call to \code{ps}, but this difference is inconsequential.
}
\examples{\donttest{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "gbm", estimand = "ATT",
                stop.method = "es.max"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "gbm", estimand = "ATE",
                stop.method = "ks.mean"))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "gbm", use.kernel = TRUE,
                stop.method = "p.max"))
summary(W3)
bal.tab(W3)
}}