\name{method_optweight}
\alias{method_optweight}

\title{Optimization-Based Weighting}

\description{
This page explains the details of estimating optimization-based weights by setting \code{method = "optweight"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating weights by solving a quadratic programming problem subject to approximate or exact balance constraints. This method relies on \code{\link[optweight]{optweight}} from the \pkg{optweight} package.

Because \code{optweight} offers finer control and uses the same syntax as \code{weightit}, it is recommended that \code{\link[optweight]{optweight}} be used instead of \code{weightit} with \code{method = "optweight"}.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the weights using \code{\link[optweight]{optweight}}. The following estimands are allowed: ATE, ATT, and ATC. The weights are taken from the output of the \code{optweight} fit object.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, this method estimates the weights using \code{\link[optweight]{optweight}}. The following estimands are allowed: ATE and ATT. The weights are taken from the output of the \code{optweight} fit object.
}
\subsection{Continuous Treatments}{
For binary treatments, this method estimates the weights using \code{\link[optweight]{optweight}}. The weights are taken from the output of the \code{optweight} fit object.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, \code{optweight} estimates weights that simultaneously satisfy balance constraints at all time points, so only one model is fit to obtain the weights. Using \code{method = "optweight"} in \code{weightitMSM} cause \code{is.MSM.method} to be set to \code{TRUE} by default. Setting it to \code{FALSE} will run one model for each time point and multiply the weights together, a method that is not recommended. NOTE: neither use of optimization-based weights with longitudinal treatments has been validated!
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
All arguments to \code{optweight} can be passed through \code{weightit} or \code{weightitMSM}, with the following exception:

\code{targets} cannot be used and is ignored.

All arguments take on the defaults of those in \code{optweight}.

}
\section{Additional Outputs}{
\describe{
\item{\code{info}}{
A list with one entry:
\describe{
\item{\code{duals}}{
A data frame of dual variables for each balance constraint.
}

}
}
\item{\code{obj}}{
When \code{include.obj = TRUE}, the output of the call to \code{\link[optweight:optweight]{optweight::optweight}}.
}
}
}
\note{
The specification of \code{tols} differs between \code{weightit} and \code{optweight}. In \code{weightit}, one tolerance value should be included per level of each factor variable, whereas in \code{optweight}, all levels of a factor are given the same tolerance, and only one value needs to be supplied for a factor variable. Because of the potential for confusion and ambiguity, it is recommended to only supply one value for \code{tols} in \code{weightit} that applies to all variables. For finer control, use \code{optweight} directly.

Seriously, just use \code{\link[optweight]{optweight}}. The syntax is almost identical and it's compatible with \pkg{cobalt}, too.
}

\references{

Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. Journal of the American Statistical Association, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}

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
                method = "optweight", estimand = "ATT",
                tols = 0))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "optweight", estimand = "ATE",
                tols = .01))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "optweight", tols = .05))
summary(W3)
bal.tab(W3)
}