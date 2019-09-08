\name{method_ps}
\alias{method_ps}

\title{Propensity Score Weighting Using Generalized Linear Models}

\description{
This page explains the details of estimating weights from generalized linear model-based propensity scores by setting \code{method = "ps"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating propensity scores with a parametric generalized linear model and then converting those propensity scores into weights using a formula that depends on the desired estimand. For binary and multinomial treatments, a binomial or multinomial regression model is used to estimate the propensity scores as the predicted probability of being in each treatment given the covariates. For continuous treatments, a generalized linear model is used to estimate generalized propensity scores as the conditional density of treatment given the covariates.


\subsection{Binary Treatments}{
For binary treatments, this method estimates the propensity scores using \code{\link{glm}}. An additional argument is \code{link}, which uses the same options as \code{link} in \code{\link{family}}. The default link is "logit", but others, including "probit", are allowed. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for the ATE, ATT, and ATC are computed from the estimated propensity scores using the standard formulas, the weights for the ATO are computed as in Li & Li (2018), and the weights for the ATM (i.e., average treatment effect in the equivalent sample "pair-matched" with calipers) are computed as in Yoshida et al (2017). When \code{include.obj = TRUE}, the returned object is the \code{glm} fit.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, the propensity scores are estimated using multinomial regression from one of two functions depending on the requested link: for logit (\code{"logit"}) and probit (\code{"probit"}) links, \code{mlogit} from the \pkg{mlogit} package is used, and for the Bayesian probit  (\code{"bayes.probit"}) link, \code{mnp} from the \pkg{MNP} package is used. If \code{mlogit} in not installed, a series of binomial regressions using \code{glm} will be run instead, with estimated propensities normalized to sum to 1. These are the only three links allowed for multinomial treatments at this time. (These methods can fail to converge, yielding errors that may seem foreign.) The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for each estimand are computed using the standard formulas or those mentioned above. When \code{include.obj = TRUE}, the returned object is the fit object from \code{mlogit} or \code{mnp} or the list of \code{glm} fit objects if binomial models are used.
}
\subsection{Continuous Treatments}{
For continuous treatments, the generalized propensity score is estimated using linear regression. In addition, kernel density estimation can be used instead of assuming a normal density for the numerator and denominator of the generalized propensity score by setting \code{use.kernel = TRUE}. Other arguments to \code{\link{density}} can be specified to refine the density estimation parameters. \code{plot = TRUE} can be specified to plot the density for the numerator and denominator, which can be helpful in diagnosing extreme weights. When \code{include.obj = TRUE}, the returned object is the \code{glm} fit from denominator model.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point.
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios except for multinomial treatments with \code{link = "bayes.probit"}. Warning messages may appear otherwise about non-integer successes, and these can be ignored.
}
\subsection{Missing Data}{
Missing data is not directly compatible with estimating propensity scores, so a few extra things happen when \code{NA}s are present in the covariates. First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is \code{NA} and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with 0s (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting \code{weightit} object will be the original covariates with the \code{NA}s.
}
}
\section{Additional Arguments}{
The following additional arguments can be specified:
\describe{
\item{\code{link}}{The link used in the generalized linear model for the propensity scores. For binary treatments, \code{link} can be any of those allowed by \code{binomial}. For multinomial treatments, \code{link} can be "logit", "probit", or "bayes.probit". For continuous treatments, \code{link} can be any of those allowed by \code{gaussian}.
}
\item{\code{use.kernel}}{If \code{TRUE}, uses kernel density estimation through \code{\link{density}} to estimate the numerator and denominator densities for the weights with continuous treatments. If \code{FALSE}, assumes a normal distribution.
}
\item{\code{bw}, \code{adjust}, \code{kernel}, \code{n}}{If \code{use.kernel = TRUE} with continuous treatments, the arguments to \code{\link{density}}. The defaults are the same as those in \code{density} except that \code{n} is 10 times the number of units in the sample.
}
\item{\code{plot}}{If \code{use.kernel = TRUE} with continuous treatments, whether to plot the estimated density.
}
}
}
\references{
\strong{Binary treatments}

 - \code{estimand = "ATO"}

Li, F., Morgan, K. L., & Zaslavsky, A. M. (2016). Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association, 0(ja), 0–0.

 - \code{estimand = "ATM"}

Li, L., & Greene, T. (2013). A Weighting Analogue to Pair Matching in Propensity Score Analysis. The International Journal of Biostatistics, 9(2). \doi{10.1515/ijb-2012-0030}

 - Other estimands

Austin, P. C. (2011). An Introduction to Propensity Score Methods for Reducing the Effects of Confounding in Observational Studies. Multivariate Behavioral Research, 46(3), 399–424. \doi{10.1080/00273171.2011.568786}

\strong{Multinomial Treatments}

 - \code{estimand = "ATO"}

Li, F., & Li, F. (2018). Propensity Score Weighting for Causal Inference with Multi-valued Treatments. ArXiv:1808.05339 [Stat]. Retrieved from http://arxiv.org/abs/1808.05339


 - \code{estimand = "ATM"}

Yoshida, K., Hernández-Díaz, S., Solomon, D. H., Jackson, J. W., Gagne, J. J., Glynn, R. J., & Franklin, J. M. (2017). Matching weights to simultaneously compare three treatment groups: Comparison to three-way matching. Epidemiology (Cambridge, Mass.), 28(3), 387–395. \doi{10.1097/EDE.0000000000000627}


 - Other estimands

McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. Statistics in Medicine, 32(19), 3388–3414. \doi{10.1002/sim.5753}

\strong{Continuous treatments}

\code{method = "ps"}

Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural Models and Causal Inference in Epidemiology. Epidemiology, 11(5), 550–560.

}
\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}
}
\examples{\dontrun{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ps", estimand = "ATT",
                link = "probit"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ps", estimand = "ATE"))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ps", use.kernel = TRUE))
summary(W3)
bal.tab(W3)
}}