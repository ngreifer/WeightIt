\name{method_super}
\alias{method_super}

\title{Propensity Score Weighting Using SuperLearner}

\description{
This page explains the details of estimating weights from SuperLearner-based propensity scores by setting \code{method = "super"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating propensity scores using the SuperLearner algorithm for stacking predictions and then converting those propensity scores into weights using a formula that depends on the desired estimand. For binary and multinomial treatments, one or more binary classification algorithms are used to estimate the propensity scores as the predicted probability of being in each treatment given the covariates. For continuous treatments, a regression algorithm is used to estimate generalized propensity scores as the conditional density of treatment given the covariates.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the propensity scores using \code{\link[SuperLearner]{SuperLearner}} in the \pkg{SuperLearner} package. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for the ATE, ATT, and ATC are computed from the estimated propensity scores using the standard formulas, the weights for the ATO are computed as in Li & Li (2018), and the weights for the ATM (i.e., average treatment effect in the equivalent sample "pair-matched" with calipers) are computed as in Yoshida et al (2017). When \code{include.obj = TRUE}, the returned object is the \code{SuperLearner} fit.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, the propensity scores are estimated using several calls to \code{\link[SuperLearner]{SuperLearner}}, one for each treatment group, and the treatment probabilities are normalized to sum to 1. The following estimands are allowed: ATE, ATT, ATO, and ATM. The weights for each estimand are computed using the standard formulas or those mentioned above. When \code{include.obj = TRUE}, the returned object is the list of fit object from the \code{SuperLearner} calls.
}
\subsection{Continuous Treatments}{
For continuous treatments, the generalized propensity score is estimated using \code{\link[SuperLearner]{SuperLearner}}. In addition, kernel density estimation can be used instead of assuming a normal density for the numerator and denominator of the generalized propensity score by setting \code{use.kernel = TRUE}. Other arguments to \code{\link{density}} can be specified to refine the density estimation parameters. \code{plot = TRUE} can be specified to plot the density for the numerator and denominator, which can be helpful in diagnosing extreme weights. When \code{include.obj = TRUE}, the returned object is the \code{\link[SuperLearner]{SuperLearner}} fit from denominator model.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point.
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios.
}
\subsection{Missing Data}{
Missing data is not compatible with SuperLearner, so a few extra things happen when \code{NA}s are present in the covariates. First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is \code{NA} and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with 0s (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting \code{weightit} object will be the original covariates with the \code{NA}s.
}
}
\section{Additional Arguments}{
An argument to \code{SL.library} must be supplied. To see a list of available entries, use \code{\link[SuperLearner]{listWrappers}}.

All arguments to \code{SuperLearner} can be passed through \code{weightit} or \code{weightitMSM}, with the following exceptions:

\code{method} in \code{SuperLearner} is replaced with the argument \code{SL.method} in \code{weightit}.

\code{obsWeights} is ignored because sampling weights are passed using \code{s.weights}.

The following additional arguments can be specified:
\describe{
\item{\code{use.kernel}}{If \code{TRUE}, uses kernel density estimation through \code{\link{density}} to estimate the numerator and denominator densities for the weights with continuous treatments. If \code{FALSE}, assumes a normal distribution.
}
\item{\code{bw}, \code{adjust}, \code{kernel}, \code{n}}{If \code{use.kernel = TRUE} with continuous treatments, the arguments to \code{\link{density}}. The defaults are the same as those in \code{density} except that \code{n} is 10 times the number of units in the sample.
}
\item{\code{plot}}{If \code{use.kernel = TRUE} with continuous treatments, whether to plot the estimated density.
}
}
}
\references{
Pirracchio, R., Petersen, M. L., & van der Laan, M. (2015). Improving Propensity Score Estimators’ Robustness to Model Misspecification Using Super Learner. American Journal of Epidemiology, 181(2), 108–119. \doi{10.1093/aje/kwu253}

}
\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}
}
\note{
Some methods formerly available in \pkg{SuperLearner} are now in SuperLearnerExtra, which can be found on GitHub at \url{https://github.com/ecpolley/SuperLearnerExtra}.
}
\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", estimand = "ATT",
                SL.library = c("SL.glm", "SL.gam",
                               "SL.knn")))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", estimand = "ATE",
                SL.library = c("SL.glm", "SL.gam",
                               "SL.knn")))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", use.kernel = TRUE,
                SL.library = c("SL.glm", "SL.gam",
                               "SL.ridge")))
summary(W3)
bal.tab(W3)
}