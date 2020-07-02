\name{method_ebal}
\alias{method_ebal}

\title{Entropy Balancing}

\description{
This page explains the details of estimating weights using entropy balancing by setting \code{method = "ebal"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating weights by minimizing the entropy of the weights subject to exact moment balancing constraints. This method relies on \code{\link[ebal]{ebalance}} from the \pkg{ebal} package for binary and multinomial treatments and on code written for \pkg{WeightIt} using \code{\link{optim}} for continuous treatments.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the weights using \code{\link[ebal]{ebalance}}. The following estimands are allowed: ATE, ATT, and ATC. The weights are taken from the output of the \code{ebalance} fit object. When the ATE is requested, \code{ebalance} is run twice, once for each treatment group.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, this method estimates the weights using \code{\link[ebal]{ebalance}}. The following estimands are allowed: ATE and ATT. The weights are taken from the output of the \code{ebalance} fit objects. When the ATE is requested, \code{ebalance} is run once for each treatment group. When the ATT is requested, \code{ebalance} is run once for each non-focal (i.e., control) group.
}
\subsection{Continuous Treatments}{
For continuous treatments, this method estimates the weights using \code{optim} using formulas described by Hainmueller (2012) and Tübbicke (2020).
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point. This method is not guaranteed to yield exact balance at each time point. NOTE: the use of entropy balancing with longitudinal treatments has not been validated!
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios.
}
\subsection{Missing Data}{
In the presence of missing data, the following value(s) for \code{missing} are allowed:
\describe{
\item{\code{"ind"} (default)}{
First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is \code{NA} and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with 0s (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting \code{weightit} object will be the original covariates with the \code{NA}s.
}
}
}
}
\section{Additional Arguments}{
All arguments to \code{ebalance} can be passed through \code{weightit} or \code{weightitMSM} and take on the defaults of those in \code{ebalance}.

For \code{base.weights}, a vector with length equal to the total number of units can be supplied, in contrast to \code{ebalance()}, which requires a vector with length equal to the number of controls. This works for continuous treatments as well. These correspond to the base weights \emph{q} in Hainmueller (2012). The estimated weights minimize the Kullback entropy divergence from the base weights, defined as \eqn{w*log(w/q)}, subject to exact balance constraints. These can be used to supply previously estimated weights so that the newly estimated weights retain the some of the proprties of the original weights while ensuring the balance constraints are met.

When \code{standardize = TRUE} in the call to \code{weightit} for binary and multinomial treatments, \code{\link[ebal]{ebalance.trim}} is run on the resulting \code{ebalance} fit objects. Doing so can reduce the variability of the weights while maintaining covariate balance.

}
\section{Additional Outputs}{
\describe{
\item{\code{obj}}{
When \code{include.obj = TRUE}, the entropy balancing model fit. For binary treatments with \code{estimand = "ATT"}, the output of the call to \code{\link[ebal:ebalance]{ebal::ebalance}} or \code{\link[ebal:ebalance.trim]{ebal::ebalance.trim}} when \code{stabilize = TRUE}. For binary treatments with \code{estimand = "ATE"} and multinomial treatments, a list of outputs of calls to \code{\link[ebal:ebalance]{ebal::ebalance}} or \code{\link[ebal:ebalance.trim]{ebal::ebalance.trim}}. For continuous treatments, the output of the call to \code{optim}, which contains the dual variables and convergence information.
}
}
}
\details{
Entropy balancing involves the specification of an optimization problem, the solution to which is then used to compute the weights. The constraints of the primal optimization problem correspond to covariate balance on the means (for binary and multinomial tretaments) or treatment-covariate covariances (for continuous treatments), positivity of the weights, and that the weights sum to a certain value. It turns out that the dual optimization problem is much easier to solve because it is over only as many variables as there are balance constraints rather than over the weights for each unit and it is unconstrained. Zhao and Percival (2017) found that entropy balancing for the ATT of a binary treatment actually involves the estimation of the coeffcients of a logistic regression propensity score model but using a specialized loss function different from that optimized with maximum likelihood. Entropy balancing is doubly robust (for the ATT) in the sense that it is consistent either when the true propensity score model is a logistic regression of the treatment on the covariates or when the trueoutcome model for the control units is a linear regression of the outcome on the covariates, and it attains a semi-parametric efficiency bound when both are true. Entropy balancing will always yield exact balance on the included terms.
}
\author{
The code for continuous treatments was kindly contributed by Stefan Tübbicke (\email{tuebbicke@empwifo.uni-potsdam.de}).
}
\references{

- Binary Treatments

Hainmueller, J. (2012). Entropy Balancing for Causal Effects: A Multivariate Reweighting Method to Produce Balanced Samples in Observational Studies. Political Analysis, 20(1), 25–46. \doi{10.1093/pan/mpr025}

Zhao, Q., & Percival, D. (2017). Entropy balancing is doubly robust. Journal of Causal Inference, 5(1). \doi{10.1515/jci-2016-0010}

- Continuous Treatments

Tübbicke, S. (2020). Entropy Balancing for Continuous Treatments. ArXiv:2001.06281 [Econ, Stat]. \url{http://arxiv.org/abs/2001.06281}

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

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal"))
summary(W3)
bal.tab(W3)

}