\name{weightitMSM}
\alias{weightitMSM}
\alias{print.weightitMSM}
\title{Generate Balancing Weights}
\usage{
weightitMSM(formula.list,
            data, method,
            stabilize = FALSE,
            exact = NULL,
            s.weights = NULL,
            verbose = FALSE,
            ...)

\method{print}{weightitMSM}(x, ...)

}
\description{
\code{weightitMSM()} allows for the easy generation of balancing weights for marginal structural models for time-varying treatments using a variety of available methods for binary, continuous, and multinomial treatments. Many of these methods exist in other packages, which \code{weightit()} calls; these packages must be installed to use the desired method. Also included are \code{print} and \code{summary} methods for examining the output.

Currently only "wide" data sets, where each row corresponds to a unit's entire variable history, are supported. You can use \code{\link{reshape}} or other functions to transform your data into this format.
}
\arguments{
  \item{formula.list}{
a list of formulas corresponding to each time point with the time-specific treatment variable on the left hand side and pre-treatment covariates to be balanced on the right hand side. The formulas must be in temporal order, and must contain all covariates to be balanced at that time point (i.e., treatments and covariates featured in early formulas should appear in later ones). Interactions and functions of covariates are allowed.
}
  \item{data}{
a data set in the form of a data frame that contains the variables in the formulas in \code{formula.list}. This must be a wide data set with exactly one row per unit.
}
  \item{method}{
a string of length 1 containing the name of the method that will be used to estimate weights. See \code{\link{weightit}} for allowable options. The default is \code{"ps"}, which estimates the weights using binary (e.g., logistic) regression.
}
  \item{stabilize}{
\code{logical}; whether or not to stabilize the weights. Stabilizing the weights involves fitting a model predicting treatment at each time point from treatment status at prior time points. If \code{TRUE}, a saturated model will be fit, essentially using the observed treatment probabilities in the numerator (for binary and multinomial treatments). This may yield an error if some combinations are not observed. Default is \code{FALSE}. To manually specify stabilization model formulas, use \code{num.formula}.
}
\item{num.formula} {
optional; a list of formulas used in the numerator of the weights if stabilization is to occur. These formulas should specify the relationship between the treatment and the stabilization factors (including previous treatments) at each time point. They should not contain time-varying covariates except those measured prior to the first treatment. Can also be specified as a single one-sided formula with the stabilization factors (other than the previous treatments) on the right hand side, which adds, for each time point, the stabilization factors to a model saturated with previous treatments. See Cole & Hernan (2008) for a discussion of how to specifiy this model; including stabilization factors can change the estimand without proper adjustment, and should be done with caution. Unless you know what you are doing, we recommend setting \code{stabilize = TRUE} and ignoring \code{num.formula}.
}
  \item{exact}{
a vector or the names of variables in data for which weighting is to be done within catgories. For example, if \code{exact = "gender"}, weights will be generated seperately within each level of the variable \code{"gender"}.
}
  \item{s.weights}{
a vector of sampling weights or the name of a variable in \code{data} that contains sampling weights. These are ignored for some methods.
}
  \item{verbose}{
whether to print additional information output by the fitting function.
}
  \item{...}{
other arguments for functions called by \code{weightit} that control aspects of fitting that are not covered by the above arguments. See Details at \code{\link{weightit}}.
}
  \item{x}{
a \code{weightitMSM} object; the output of a call to \code{weightitMSM()}.
}

}
\value{
  A \code{weightit} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{treat.list}{A list of the values of the time-varying treatment variables.}
  \item{covs.list}{A list of the covariates used in the fitting at each time point. Only includes the raw covariates, which may have been altered in the fitting process.}
  \item{data}{The data.frame originally entered to \code{weightitMSM()}.}
  \item{estimand}{"ATE", currently the only estimand for MSMs.}
  \item{method}{The weight estimation method specified.}
  \item{ps.list}{A list of the estimated propensity scores (if any) at each time point.}
  \item{s.weights}{The provided sampling weights.}
  \item{treat.type}{The type of treatment: binary, continuous, or multinomial (\code{"multi"}).}
}
\details{
\code{weightitMSM()} works by seperating the estimation of weights into seperate procedures for each time period based on the formulas provided. For each formula, \code{weightitMSM()} simple calls \code{weightit()} to that formula, collects the weights for each time period, and multiplies them together to arrive at longitudinal balancing  weights.

Each formula should contain all the covariates to be balanced on. For example, the formula corresponding to the second time period should contain all the baseline covariates, the treatment variable at the first time period, and the time-varying covariates that took on values after the first treatment and before the second. Currently, only wide data sets are supported, where each unit is represented by exactly one row that contains the covariate and treatment history encoded in seperate variables.

The \code{"gbm"} method, which calls \code{ps()} in \pkg{twang}, yields the same results to a call to \code{iptw()} in \pkg{twang}. However, the \code{cbps} method, which calls \code{CBPS()} in \pkg{CBPS}, will yield different results from \code{CPMSM()} in \pkg{CBPS} because \code{CPMSM()} takes a different approach to generating weights than simply estimating several time-specific models.

}

\author{
Noah Greifer \email{noah@unc.edu}
}

\references{

Cole, S. R., & Hernan, M. A. (2008). Constructing Inverse Probability Weights for Marginal Structural Models. American Journal of Epidemiology, 168(6), 656–664. <doi:10.1093/aje/kwn164>




}

\note{

}
\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups
(W1 <- weightit(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "ps", estimand = "ATT"))
summary(W1)
bal.tab(W1)

#Balancing covariates among races (3-level)
(W2 <- weightit(race ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "cbps", estimand = "ATE",
                over = FALSE))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re78 (continuous)
(W3 <- weightit(re78 ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "ps", num.formula = ~ age +
                educ + married))
summary(W3)
bal.tab(W3)

}