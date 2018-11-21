\name{method_ebcw}
\alias{method_ebcw}

\title{Empirical Balancing Calibration Weighting}

\description{
This page explains the details of estimating weights using entropy balancing by setting \code{method = "ebal"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary and multinomial treatments.

In general, this method relies on estimating weights by minimizing a function of the weights subject to exact moment balancing constraints. This method relies on \code{\link[ATE]{ATE}} from the \pkg{ATE} package.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the weights using \code{\link[ATE]{ATE}} with \code{ATT = TRUE}. The following estimands are allowed: ATE, ATT, and ATC. The weights are taken from the output of the \code{ATE} fit object. When the ATE is requested, \code{ATE} is run twice, once for each treatment group. When \code{include.obj = TRUE}, the returned object is the \code{ebal} fit (or a list of the two fits when the estimand is the ATE).
}
\subsection{Multinomial Treatments}{
For multinomial treatments, this method estimates the weights using \code{\link[ATE]{ATE}} with \code{ATT = TRUE}. The following estimands are allowed: ATE and ATT. The weights are taken from the output of the \code{ATE} fit objects. When the ATE is requested, \code{ATE} is run once for each treatment group. When the ATT is requested, \code{ATE} is run once for each non-focal (i.e., control) group. When \code{include.obj = TRUE}, the returned object is the list of \code{ATE} fits.
}
\subsection{Continuous Treatments}{
Continuous treatments are not supported.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point. This method is not guaranteed to yield exact balance at each time point.
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios.
}
}
\section{Additional Arguments}{
All argument to \code{ATE} can be passed through \code{weightit} or \code{weightitMSM}, with the following exceptions:

\code{ATT} is ignored because the estimand is passed using \code{estimand}.

All arguments take on the defaults of those in \code{ATE}.

}
\references{

Chan, K. C. G., Yam, S. C. P., & Zhang, Z. (2016). Globally efficient non-parametric inference of average treatment effects by empirical balancing calibration weighting. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(3), 673–700. \doi{10.1111/rssb.12129}

}
\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}
}

\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "ebcw", estimand = "ATT"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "ebcw", estimand = "ATE"))
summary(W2)
bal.tab(W2)

}