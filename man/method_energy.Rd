\name{method_energy}
\alias{method_energy}

\title{Energy Balancing}

\description{
This page explains the details of estimating weights using energy balancing by setting \code{method = "energy"} in the call to \fun{weightit} or \fun{weightitMSM}. This method can be used with binary and multinomial treatments. %This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating weights by minimizing an energy statistic related to covariate balance. For binary and multinomial treatments, this is the energy distance, which is a multivariate distance between distributions, between treatment groups.
%For continuous treatments, this is the energy covariance between the treatment variable and the covariates.
This method relies on code written for \pkg{WeightIt} using \pkgfun{osqp}{osqp} from the \CRANpkg{osqp} package to perform the optimization. This method may be slow or memory-intensive for large datasets.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the weights using \code{osqp()} using formulas described by Huling and Mak (2020). The following estimands are allowed: ATE, ATT, and ATC.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, this method estimates the weights using \code{osqp()} using formulas described by Huling and Mak (2020). The following estimands are allowed: ATE and ATT.
}
\subsection{Continuous Treatments}{
Continuous treatments are not currently supported.
%For continuous treatments, this method estimates the weights using \code{osqp} using (unpublished) formulas based on the original formulas for energy covariance. The means of the included covariates are constrained to be equal to their unweighted means. \strong{This method has not been studied and should be used with caution.}
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point. This method is not guaranteed to yield optimal balance at each time point. NOTE: the use of energy balancing with longitudinal treatments has not been validated!
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios. In some cases, sampling weights will cause the optimization to fail due to lack of convexity or infeasible constraints.
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
For binary and multinomial treatments, the following additional arguments can be specified:
\describe{
\item{\code{improved}}{\code{logical}; whether to use the improved energy balancing weights as described by Huling and Mak (2020) when \code{estimand = "ATE"}. This involves optimizing balance not only between each treatment group and the overall sample, but also between each pair of treatment groups. Huling and Mak (2020) found that the improved energy balancing weights generally outperformed standard energy balancing. Default is \code{TRUE}; set to \code{FALSE} to use the standard energy balancing weights instead (not recommended).
}
\item{\code{dist.mat}}{the name of the method used to compute the distance matrix or the numeric distance matrix itself. Allowable options include \code{"scaled_euclidean"} for the Euclidean (L2) distance on the scaled covariates (the default), \code{"mahalanobis"} for the Mahalanobis distance, and \code{"euclidean"} for the raw Euclidean distance. Abbreviations allowed. Note that some user-supplied distance matrices can cause the R session to abort due to a bug within \pkg{osqp}, so this argument should be used with caution. A distance matrix must be a square, symmetric, numeric matrix with zeros along the diagonal and a row and column for each unit. Can also be supplied as the output of a call to \fun{dist}.
}
\item{\code{lambda}}{a positive numeric scalar used to penalize the square of the weights. This value divided by the square of the total sample size is added to the diagonal of the quadratic part of the loss function. Higher values favor weights with less variability. Note this is distinct from the lambda value described in Huling and Mak (2020), which penalizes the complexity of individual treatment rules rather than the weights.
}
}
The \code{moments} argument functions differently for \code{method = "energy"} from how it does with other methods. When unspecified or set to zero, energy balancing weights are estimated as described by Huling and Mak (2020). When \code{moments} is set to an integer larger than 0, additional balance constraints on the requested moments of the covariates are also included, guaranteeing exact moment balance on these covariates while minimizing the energy distance of the weighted sample. For binary and multinomial treatments, this involves exact balance on the means of the entered covariates.
%For binary and multinomial treatments, this involves exact balance on the means of the entered covariates; for continuous treatments, this involves exact balance on the treatment-covariate correlations of the entered covariates.

}
\section{Additional Outputs}{
\describe{
\item{\code{obj}}{
When \code{include.obj = TRUE}, the output of the call to \pkgfun{osqp}{solve_osqp}, which contains the dual variables and convergence information.
}
}
}
\note{
Sometimes the optimization can fail to converge because the problem is not convex. A warning will be displayed if so. In these cases, try simply re-fitting the weights without changing anything. If the method repeatedly fails, you should try another method or change the supplied parameters (though this is uncommon).
}
\author{
Noah Greifer
}
\references{

%\strong{Binary and Multinomial treatments}

Huling, J. D., & Mak, S. (2020). Energy Balancing of Covariate Distributions. ArXiv:2004.13962 [Stat]. \url{https://arxiv.org/abs/2004.13962}

%\strong{Continuous treatments}

%This method has not been described in the literature yet. Zhu, Coffman, and Ghosh (2015) considered the distance covariance as an objective function and omnibus measure of balance but used a weighted boostrap to calculate it rather than using a formula for the weighted distance covariance. Andrei and McCarthy (2019) proposed the (unweighted) distance covariance as a way to measure balance in matched samples for binary treatments, but didn't propose a weighted version or consider its use with continuous treatments.

%Andrei, A.-C., & McCarthy, P. M. (2019). An omnibus approach to assess covariate balance in observational studies using the distance covariance. Statistical Methods in Medical Research, 096228021987821. \doi{10.1177/0962280219878215}

%Zhu, Y., Coffman, D. L., & Ghosh, D. (2015). A Boosting Algorithm for Estimating Generalized Propensity Scores with Continuous Treatments. Journal of Causal Inference, 3(1). \doi{10.1515/jci-2014-0022}

}
\seealso{
\fun{weightit}, \fun{weightitMSM}
}

\examples{\dontshow{if (requireNamespace("osqp", quietly = TRUE)) force(\{ # examplesIf}
library("cobalt")
data("lalonde", package = "cobalt")

#Examples may not converge, but may after several runs
\dontrun{
#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "energy", estimand = "ATE"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "energy", estimand = "ATT",
                focal = "black"))
summary(W2)
bal.tab(W2)
}
\dontshow{\}) # examplesIf}}

%#Balancing covariates with respect to re75 (continuous)
%(W3 <- weightit(re75 ~ age + educ + married +
%                  nodegree + re74, data = lalonde,
%                method = "energy", moments = 1))
%summary(W3)
%bal.tab(W3, poly = 2)