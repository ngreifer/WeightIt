\name{method_ebal}
\alias{method_ebal}

\title{Entropy Balancing}

\description{
This page explains the details of estimating weights using entropy balancing by setting \code{method = "ebal"} in the call to \fun{weightit} or \fun{weightitMSM}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating weights by minimizing the negative entropy of the weights subject to exact moment balancing constraints. This method relies on code written for \pkg{WeightIt} using \fun{optim}.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the weights using \code{optim()} using formulas described by Hainmueller (2012). The following estimands are allowed: ATE, ATT, and ATC. When the ATE is requested, the optimization is run twice, once for each treatment group.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, this method estimates the weights using \code{optim()}. The following estimands are allowed: ATE and ATT. When the ATE is requested, \code{optim()} is run once for each treatment group. When the ATT is requested, \code{optim()} is run once for each non-focal (i.e., control) group.
}
\subsection{Continuous Treatments}{
For continuous treatments, this method estimates the weights using \code{optim()} using formulas described by Tübbicke (2020).
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
\code{moments} and \code{int} are accepted. See \fun{weightit} for details.
\describe{
\item{\code{base.weights}}{
A vector of base weights, one for each unit. This works for continuous treatments as well. These correspond to the base weights \eqn{q} in Hainmueller (2012). The estimated weights minimize the Kullback entropy divergence from the base weights, defined as \eqn{\sum w \log(w/q)}, subject to exact balance constraints. These can be used to supply previously estimated weights so that the newly estimated weights retain the some of the properties of the original weights while ensuring the balance constraints are met. Sampling weights should not be passed to \code{base.weights} but can be included in a \code{weightit()} call that includes \code{s.weights}.
}
\item{\code{d.moments}}{
With continuous treatments, the number of moments of the treatment and covariate distributions that are constrained to be the same in the weighted sample as in the original sample. For example, setting \code{d.moments = 3} ensures that the mean, variance, and skew of the treatment and covariates are the same in the weighted sample as in the unweighted sample. \code{d.moments} should be greater than or equal to \code{moments} and will be automatically set accordingly if not (or if not specified). Vegetabile et al. (2020) recommend setting \code{d.moments = 3}, even if \code{moments} is less than 3. This argument corresponds to the tuning parameters \eqn{r} and \eqn{s} in Vegetabile et al. (2020) (which here must be equal). Ignored for binary and multi-category treatments.
}
}

The arguments \code{maxit} and \code{reltol} can be supplied and are passed to the \code{control} argument of \fun{optim}. The \code{"BFGS"} method is used, so the defaults correspond to this.

The \code{stabilize} argument is ignored; in the past it would reduce the variability of the weights through an iterative process. If you want to minimize the variance of the weights subject to balance constraints, use \code{method = "optweight"}.

}
\section{Additional Outputs}{
\describe{
\item{\code{obj}}{
When \code{include.obj = TRUE}, the output of the call to \fun{optim}, which contains the dual variables and convergence information. For ATE fits or with multinomial treatments, a list of \code{optim()} outputs, one for each weighted group.
}
}
}
\details{
Entropy balancing involves the specification of an optimization problem, the solution to which is then used to compute the weights. The constraints of the primal optimization problem correspond to covariate balance on the means (for binary and multinomial treatments) or treatment-covariate covariances (for continuous treatments), positivity of the weights, and that the weights sum to a certain value. It turns out that the dual optimization problem is much easier to solve because it is over only as many variables as there are balance constraints rather than over the weights for each unit and it is unconstrained. Zhao and Percival (2017) found that entropy balancing for the ATT of a binary treatment actually involves the estimation of the coefficients of a logistic regression propensity score model but using a specialized loss function different from that optimized with maximum likelihood. Entropy balancing is doubly robust (for the ATT) in the sense that it is consistent either when the true propensity score model is a logistic regression of the treatment on the covariates or when the true outcome model for the control units is a linear regression of the outcome on the covariates, and it attains a semi-parametric efficiency bound when both are true. Entropy balancing will always yield exact mean balance on the included terms.
}

\references{

\strong{Binary Treatments}

Hainmueller, J. (2012). Entropy Balancing for Causal Effects: A Multivariate Reweighting Method to Produce Balanced Samples in Observational Studies. Political Analysis, 20(1), 25–46. \doi{10.1093/pan/mpr025}

Källberg, D., & Waernbaum, I. (2022). Large Sample Properties of Entropy Balancing Estimators of Average Causal Effects. ArXiv:2204.10623 [Stat]. \url{https://arxiv.org/abs/2204.10623}

Zhao, Q., & Percival, D. (2017). Entropy balancing is doubly robust. Journal of Causal Inference, 5(1). \doi{10.1515/jci-2016-0010}

\strong{Continuous Treatments}

Tübbicke, S. (2022). Entropy Balancing for Continuous Treatments. Journal of Econometric Methods, 11(1), 71–89. \doi{10.1515/jem-2021-0002}

Vegetabile, B. G., Griffin, B. A., Coffman, D. L., Cefalu, M., Robbins, M. W., & McCaffrey, D. F. (2021). Nonparametric estimation of population average dose-response curves using entropy balancing weights for continuous exposures. Health Services and Outcomes Research Methodology, 21(1), 69–110. \doi{10.1007/s10742-020-00236-2}

}
\seealso{
\fun{weightit}, \fun{weightitMSM}

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
                method = "ebal", estimand = "ATE"))
summary(W2)
bal.tab(W2)

#Balancing covariates and squares with respect to
#re75 (continuous), maintaining 3 moments of the
#covariate and treatment distributions
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", moments = 2,
                d.moments = 3))
summary(W3)
bal.tab(W3)

}