% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weightit2ebal.R
\name{method_ebal}
\alias{method_ebal}
\title{Entropy Balancing}
\description{
This page explains the details of estimating weights using
entropy balancing by setting \code{method = "ebal"} in the call to \code{\link[=weightit]{weightit()}} or
\code{\link[=weightitMSM]{weightitMSM()}}. This method can be used with binary, multi-category, and
continuous treatments.

In general, this method relies on estimating weights by minimizing the
negative entropy of the weights subject to exact moment balancing
constraints. This method relies on code written for \pkg{WeightIt} using
\code{\link[=optim]{optim()}}.
\subsection{Binary Treatments}{

For binary treatments, this method estimates the weights using \code{optim()}
using formulas described by Hainmueller (2012). The following estimands are
allowed: ATE, ATT, and ATC. When the ATE is requested, the optimization is
run twice, once for each treatment group.
}

\subsection{Multi-Category Treatments}{

For multi-category treatments, this method estimates the weights using
\code{optim()}. The following estimands are allowed: ATE and ATT. When the ATE is
requested, \code{optim()} is run once for each treatment group. When the ATT is
requested, \code{optim()} is run once for each non-focal (i.e., control) group.
}

\subsection{Continuous Treatments}{

For continuous treatments, this method estimates the weights using \code{optim()}
using formulas described by Tübbicke (2022) and Vegetabile et al. (2021).
}

\subsection{Longitudinal Treatments}{

For longitudinal treatments, the weights are the product of the weights
estimated at each time point. This method is not guaranteed to yield exact
balance at each time point. NOTE: the use of entropy balancing with
longitudinal treatments has not been validated!
}

\subsection{Sampling Weights}{

Sampling weights are supported through \code{s.weights} in all scenarios.
}

\subsection{Missing Data}{

In the presence of missing data, the following value(s) for \code{missing} are
allowed:
\describe{
\item{\code{"ind"} (default)}{
First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is \code{NA} and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with the covariate medians (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting \code{weightit} object will be the original covariates with the \code{NA}s.
}
}
}

\subsection{M-estimation}{

M-estimation is supported for all scenarios. See \code{\link[=glm_weightit]{glm_weightit()}} and
\code{vignette("estimating-effects")} for details.
}
}
\details{
Entropy balancing involves the specification of an optimization
problem, the solution to which is then used to compute the weights. The
constraints of the primal optimization problem correspond to covariate
balance on the means (for binary and multi-category treatments) or
treatment-covariate covariances (for continuous treatments), positivity of
the weights, and that the weights sum to a certain value. It turns out that
the dual optimization problem is much easier to solve because it is over only
as many variables as there are balance constraints rather than over the
weights for each unit and it is unconstrained. Zhao and Percival (2017) found
that entropy balancing for the ATT of a binary treatment actually involves
the estimation of the coefficients of a logistic regression propensity score
model but using a specialized loss function different from that optimized
with maximum likelihood. Entropy balancing is doubly robust (for the ATT) in
the sense that it is consistent either when the true propensity score model
is a logistic regression of the treatment on the covariates or when the true
outcome model for the control units is a linear regression of the outcome on
the covariates, and it attains a semi-parametric efficiency bound when both
are true. Entropy balancing will always yield exact mean balance on the
included terms.
}
\section{Additional Arguments}{
 \code{moments} and \code{int} are accepted. See
\code{\link[=weightit]{weightit()}} for details.

\describe{
\item{\code{base.weights}}{
A vector of base weights, one for each unit. These correspond to the base weights $q$ in Hainmueller (2012). The estimated weights minimize the Kullback entropy divergence from the base weights, defined as \eqn{\sum w \log(w/q)}, subject to exact balance constraints. These can be used to supply previously estimated weights so that the newly estimated weights retain the some of the properties of the original weights while ensuring the balance constraints are met. Sampling weights should not be passed to \code{base.weights} but can be included in a \code{weightit()} call that includes \code{s.weights}.
}
\item{\code{reltol}}{the relative tolerance for convergence of the optimization. Passed to the \code{control} argument of \code{optim()}. Default is \code{1e-10}.
}
\item{\code{maxit}}{the maximum number of iterations for convergence of the optimization. Passed to the \code{control} argument of \code{optim()}. Default is 1000 for binary and multi-category treatments and 10000 for continuous and longitudinal treatments.
}
\item{\code{solver}}{the solver to use to estimate the parameters of the just-identified CBPS. Allowable options include \code{"multiroot"} to use \pkgfun{rootSolve}{multiroot} and \code{"optim"} to use \code{\link[stats:optim]{stats::optim()}}. \code{"multiroot"} is the default when \pkg{rootSolve} is installed, as it tends to be much faster and more accurate; otherwise, \code{"optim"} is the default and requires no dependencies. Regardless of \code{solver}, the output of \code{optim()} is returned when \code{include.obj = TRUE} (see below). When \code{over = TRUE}, the parameter estimates of the just-identified CBPS are used as starting values for the over-identified CBPS.
}
\item{\code{quantile}}{
A named list of quantiles (values between 0 and 1) for each continuous covariate, which are used to create additional variables that when balanced ensure balance on the corresponding quantile of the variable. For example, setting \verb{quantile = list(x1 = c(.25, .5. , .75))} ensures the 25th, 50th, and 75th percentiles of \code{x1} in each treatment group will be balanced in the weighted sample. Can also be a single number (e.g., \code{.5}) or an unnamed list of length 1 (e.g., \code{list(c(.25, .5, .75))}) to request the same quantile(s) for all continuous covariates, or a named vector (e.g., \code{c(x1 = .5, x2 = .75)} to request one quantile for each covariate. Only allowed with binary and multi-category treatments.
}
\item{\code{d.moments}}{
With continuous treatments, the number of moments of the treatment and covariate distributions that are constrained to be the same in the weighted sample as in the original sample. For example, setting \code{d.moments = 3} ensures that the mean, variance, and skew of the treatment and covariates are the same in the weighted sample as in the unweighted sample. \code{d.moments} should be greater than or equal to \code{moments} and will be automatically set accordingly if not (or if not specified). Vegetabile et al. (2021) recommend setting \code{d.moments = 3}, even if \code{moments} is less than 3. This argument corresponds to the tuning parameters $r$ and $s$ in Vegetabile et al. (2021) (which here are set to be equal). Ignored for binary and multi-category treatments.
}
}

The \code{stabilize} argument is ignored; in the past it would reduce the
variability of the weights through an iterative process. If you want to
minimize the variance of the weights subject to balance constraints, use
\code{method = "optweight"}.
}

\section{Additional Outputs}{

\describe{
\item{\code{obj}}{When \code{include.obj = TRUE}, the output of the call to \code{\link[=optim]{optim()}}, which contains the dual variables and convergence information. For ATE fits or with multi-category treatments, a list of \code{optim()} outputs, one for each weighted group.
}
}
}

\examples{
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATT"))
summary(W1)
cobalt::bal.tab(W1)

#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATE"))
summary(W2)
cobalt::bal.tab(W2)

#Balancing covariates and squares with respect to
#re75 (continuous), maintaining 3 moments of the
#covariate and treatment distributions
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", moments = 2,
                d.moments = 3))
summary(W3)
cobalt::bal.tab(W3)
}
\references{
\subsection{Binary Treatments}{
\subsection{\code{estimand = "ATT"} Hainmueller, J. (2012). Entropy Balancing for Causal}{

Effects: A Multivariate Reweighting Method to Produce Balanced Samples in
Observational Studies. \emph{Political Analysis}, 20(1), 25–46.
\doi{10.1093/pan/mpr025}

Zhao, Q., & Percival, D. (2017). Entropy balancing is doubly robust. \emph{Journal
of Causal Inference}, 5(1). \doi{10.1515/jci-2016-0010}
}

\subsection{\code{estimand = "ATE"}}{

Källberg, D., & Waernbaum, I. (2023). Large Sample Properties of Entropy
Balancing Estimators of Average Causal Effects. \emph{Econometrics and
Statistics}. \doi{10.1016/j.ecosta.2023.11.004}
}

}

\subsection{Continuous Treatments}{

Tübbicke, S. (2022). Entropy Balancing for Continuous Treatments. \emph{Journal of
Econometric Methods}, 11(1), 71–89. \doi{10.1515/jem-2021-0002}

Vegetabile, B. G., Griffin, B. A., Coffman, D. L., Cefalu, M., Robbins, M.
W., & McCaffrey, D. F. (2021). Nonparametric estimation of population average
dose-response curves using entropy balancing weights for continuous
exposures. \emph{Health Services and Outcomes Research Methodology}, 21(1),
69–110. \doi{10.1007/s10742-020-00236-2}
}
}
\seealso{
\code{\link[=weightit]{weightit()}}, \code{\link[=weightitMSM]{weightitMSM()}}

\link{method_ipt} and \link{method_cbps} for inverse probability tilting and CBPS,
which work similarly.
}
