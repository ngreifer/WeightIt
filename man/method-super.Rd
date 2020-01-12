\name{method_super}
\alias{method_super}

\title{Propensity Score Weighting Using SuperLearner}

\description{
This page explains the details of estimating weights from SuperLearner-based propensity scores by setting \code{method = "super"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating propensity scores using the SuperLearner algorithm for stacking predictions and then converting those propensity scores into weights using a formula that depends on the desired estimand. For binary and multinomial treatments, one or more binary classification algorithms are used to estimate the propensity scores as the predicted probability of being in each treatment given the covariates. For continuous treatments, a regression algorithm is used to estimate generalized propensity scores as the conditional density of treatment given the covariates.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the propensity scores using \code{\link[SuperLearner]{SuperLearner}} in the \pkg{SuperLearner} package. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for the ATE, ATT, and ATC are computed from the estimated propensity scores using the standard formulas, the weights for the ATO are computed as in Li & Li (2018), and the weights for the ATM (i.e., average treatment effect in the equivalent sample "pair-matched" with calipers) are computed as in Yoshida et al (2017). Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See \code{\link{get_w_from_ps}} for details.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, the propensity scores are estimated using several calls to \code{\link[SuperLearner]{SuperLearner}}, one for each treatment group, and the treatment probabilities are normalized to sum to 1. The following estimands are allowed: ATE, ATT, ATO, and ATM. The weights for each estimand are computed using the standard formulas or those mentioned above. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See \code{\link{get_w_from_ps}} for details.
}
\subsection{Continuous Treatments}{
For continuous treatments, the generalized propensity score is estimated using \code{\link[SuperLearner]{SuperLearner}}. In addition, kernel density estimation can be used instead of assuming a normal density for the numerator and denominator of the generalized propensity score by setting \code{use.kernel = TRUE}. Other arguments to \code{\link{density}} can be specified to refine the density estimation parameters. \code{plot = TRUE} can be specified to plot the density for the numerator and denominator, which can be helpful in diagnosing extreme weights.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point.
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios.
}
\subsection{Missing Data}{
In the presence of missing data, the following value(s) for \code{missing} are allowed:
\describe{
\item{\code{"ind"} (default)}{
First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is \code{NA} and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with 0s. The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting \code{weightit} object will be the original covariates with the \code{NA}s.
}
}
}
}
\section{Additional Arguments}{
\describe{
\item{\code{discrete}}{if \code{TRUE}, uses discrete SuperLearner, which simply selects the best performing method. Default \code{FALSE}, which finds the optimal combination of predictions for the libraries using \code{SL.method}.}
}

An argument to \code{SL.library} \strong{must} be supplied. To see a list of available entries, use \code{\link[SuperLearner]{listWrappers}}.

All arguments to \code{\link[SuperLearner]{SuperLearner}} can be passed through \code{weightit} or \code{weightitMSM}, with the following exceptions:

\code{method} in \code{SuperLearner} is replaced with the argument \code{SL.method} in \code{weightit}.

\code{obsWeights} is ignored because sampling weights are passed using \code{s.weights}.

For continuous treatments only, the following arguments may be supplied:
\describe{
\item{\code{density}}{A function corresponding the conditional density of the treatment. The standardized residuals of the treatment model will be fed through this function to produce the numerator and denominator of the generalized propensity score weights. If blank, \code{\link{dnorm}} is used as recommended by Robins et al. (2000). This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if \code{density = "dt_2"} is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme outcome values are observed (Naimi et al., 2014). Ignored if \code{use.kernel = TRUE} (described below).
}
\item{\code{use.kernel}}{If \code{TRUE}, uses kernel density estimation through the \code{\link{density}} function to estimate the numerator and denominator densities for the weights. If \code{FALSE}, the argument to the \code{density} parameter is used instead.
}
\item{\code{bw}, \code{adjust}, \code{kernel}, \code{n}}{If \code{use.kernel = TRUE}, the arguments to the \code{\link{density}} function. The defaults are the same as those in \code{density} except that \code{n} is 10 times the number of units in the sample.
}
\item{\code{plot}}{If \code{use.kernel = TRUE}, whether to plot the estimated density.
}
}
}
\section{Additional Outputs}{
\describe{
\item{\code{info}}{
For binary treatments, a list with two entries, \code{coef} and \code{cvRisk}. For multinomial treatments, a list of lists with these two entries, one for each treatment level:
\describe{
\item{\code{coef}}{
The coefficients in the linear combination of the predictions from each method in \code{SL.library}. Higher values indicate that the corresponding method plays a larger role in determining the resulting predicted value, and values close to zero indicate that the method plays little role in determining the predicted value. When \code{discrete = TRUE}, these correspond to the coefficients that would have been estimated had \code{discrete} been \code{TRUE}.
}
\item{\code{cvRisk}}{
The cross-validation risk for each method in \code{SL.library}. Higher values indicate that the method has worse cross-validation accuracy.
}
}
}
\item{\code{obj}}{
When \code{include.obj = TRUE}, the SuperLearner fit(s) used to generate the predicted values. For binary and continuous treatments, the output of the call to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}. For multinomial treatments, a list of outputs to calls to \code{SuperLearner::SuperLearner}.
}
}
}
\references{
- Binary treatments

Pirracchio, R., Petersen, M. L., & van der Laan, M. (2015). Improving Propensity Score Estimators’ Robustness to Model Misspecification Using Super Learner. American Journal of Epidemiology, 181(2), 108–119. \doi{10.1093/aje/kwu253}

- Continuous treatments

Kreif, N., Grieve, R., Díaz, I., & Harrison, D. (2015). Evaluation of the Effect of a Continuous Treatment: A Machine Learning Approach with an Application to Treatment for Traumatic Brain Injury. Health Economics, 24(9), 1213–1228. \doi{10.1002/hec.3189}

See \code{\link{method_ps}} for additional references.
}
\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}, \code{\link{get_w_from_ps}}
}
\note{
Some methods formerly available in \pkg{SuperLearner} are now in \pkg{SuperLearnerExtra}, which can be found on GitHub at \url{https://github.com/ecpolley/SuperLearnerExtra}.
}
\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", estimand = "ATT",
                SL.library = c("SL.glm", "SL.stepAIC",
                               "SL.glm.interaction")))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", estimand = "ATE",
                SL.library = c("SL.glm", "SL.stepAIC",
                               "SL.glm.interaction")))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
#assuming t(8) conditional density for treatment
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "super", density = "dt_8",
                SL.library = c("SL.glm", "SL.ridge",
                               "SL.glm.interaction")))
summary(W3)
bal.tab(W3)
}