\name{method_ps}
\alias{method_ps}

\title{Propensity Score Weighting Using Generalized Linear Models}

\description{
This page explains the details of estimating weights from generalized linear model-based propensity scores by setting \code{method = "ps"} in the call to \fun{weightit} or \fun{weightitMSM}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating propensity scores with a parametric generalized linear model and then converting those propensity scores into weights using a formula that depends on the desired estimand. For binary and multinomial treatments, a binomial or multinomial regression model is used to estimate the propensity scores as the predicted probability of being in each treatment given the covariates. For ordinal treatments, an ordinal regression model is used to estimate generalized propensity scores. For continuous treatments, a generalized linear model is used to estimate generalized propensity scores as the conditional density of treatment given the covariates.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the propensity scores using \fun{glm}. An additional argument is \code{link}, which uses the same options as \code{link} in \fun{family}. The default link is "logit", but others, including "probit", are allowed. The following estimands are allowed: ATE, ATT, ATC, ATO, ATM, and ATOS. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See \fun{get_w_from_ps} for details.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, the propensity scores are estimated using multinomial regression from one of a few functions depending on the requested link: for logit (\code{"logit"}) and probit (\code{"probit"}) links, \pkgfun{mlogit}{mlogit} from the \pkg{mlogit} package is used; for the Bayesian probit (\code{"bayes.probit"}) link, \pkgfun{MNP}{mnp} from the \pkg{MNP} package is used; and for the biased-reduced multinomial logistic regression (\code{"br.logit"}), \pkgfun{brglm2}{brmultinom} from the \pkg{brglm2} package is used. If the treatment variable is an ordered factor, \pkgfun{MASS}{polr} from the \pkg{MASS} package is used to fit ordinal regression. Any of the methods allowed in the \code{method} argument of \code{polr()} can be supplied to \code{link}. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for each estimand are computed using the standard formulas or those mentioned above. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See \fun{get_w_from_ps} for details.
}
\subsection{Continuous Treatments}{
For continuous treatments, the generalized propensity score is estimated using linear regression. The conditional density can be specified as normal or another distribution. In addition, kernel density estimation can be used instead of assuming a specific density for the numerator and denominator of the generalized propensity score by setting \code{use.kernel = TRUE}. Other arguments to \fun{density} can be specified to refine the density estimation parameters. \code{plot = TRUE} can be specified to plot the density for the numerator and denominator, which can be helpful in diagnosing extreme weights.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point.
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios except for multinomial treatments with \code{link = "bayes.probit"}. Warning messages may appear otherwise about non-integer successes, and these can be ignored.
}
\subsection{Missing Data}{
In the presence of missing data, the following value(s) for \code{missing} are allowed:
\describe{
\item{\code{"ind"} (default)}{
First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is \code{NA} and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with 0s (this value is arbitrary and does not affect estimation). The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting \code{weightit} object will be the original covariates with the \code{NA}s.
}
}
%\item{\code{"saem"}}{
%For binary treatments with \code{link = "logit"} or continuous treatments, a stochastic approximation version of the EM algorithm (SAEM) is used via the \CRANpkg{misaem} package. No additional covariates are created. See Jiang et al. (2019) for information on this method. In some cases, this is a suitable alternative to multiple imputation. NOTE: because the \pkg{misaem} package has been archived, this method is not available.
%}
}
}
\section{Additional Arguments}{
The following additional arguments can be specified:
\describe{
\item{\code{link}}{The link used in the generalized linear model for the propensity scores. For binary treatments, \code{link} can be any of those allowed by \fun{binomial}. A \code{br.} prefix can be added (e.g., \code{"br.logit"}); this changes the fitting method to the \link[brglm2:brglmFit]{bias-corrected generalized linear models} implemented in the \CRANpkg{brglm2} package. For multinomial treatments, \code{link} can be "logit", "probit", "bayes.probit", or "br.logit". For ordered treatments, \code{link} can be any of those allowed by the \code{method} argument of \pkgfun{MASS}{polr}. For continuous treatments, \code{link} can be any of those allowed by \fun{gaussian}.
}
}
For continuous treatments only, the following arguments may be supplied:
\describe{
\item{\code{density}}{A function corresponding the conditional density of the treatment. The standardized residuals of the treatment model will be fed through this function to produce the numerator and denominator of the generalized propensity score weights. If blank, \fun{dnorm} is used as recommended by Robins et al. (2000). This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if \code{density = "dt_2"} is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme outcome values are observed (Naimi et al., 2014). Ignored if \code{use.kernel = TRUE} (described below).
}
\item{\code{use.kernel}}{If \code{TRUE}, uses kernel density estimation through the \fun{density} function to estimate the numerator and denominator densities for the weights. If \code{FALSE}, the argument to the \code{density} parameter is used instead.
}
\item{\code{bw}, \code{adjust}, \code{kernel}, \code{n}}{If \code{use.kernel = TRUE}, the arguments to the \fun{density} function. The defaults are the same as those in \code{density} except that \code{n} is 10 times the number of units in the sample.
}
\item{\code{plot}}{If \code{use.kernel = TRUE} with continuous treatments, whether to plot the estimated density.
}
}

For binary treatments, additional arguments to \code{glm()} can be specified as well. The \code{method} argument in \code{glm()} is renamed to \code{glm.method}. This can be used to supply alternative fitting functions, such as those implemented in the \CRANpkg{glm2} package. Other arguments to \code{weightit()} are passed to \code{...} in \code{glm()}. %In the presence of missing data with \code{link = "logit"} and \code{missing = "saem"}, additional arguments are passed to \code{\link[misaem]{miss.glm}} and \code{\link[misaem]{predict.miss.glm}}, except the \code{method} argument in \code{\link[misaem]{predict.miss.glm}} is replaced with \code{saem.method}.

For multi-category treatments with \code{link = "logit"} or \code{"probit"}, the default is to use multinomial logistic or probit regression using the \CRANpkg{mlogit} package. To request that separate binary logistic or probit regressions are run instead, set \code{use.mlogit = FALSE}. This can be helpful when \code{mlogit} is slow or fails to converge. With \code{link = "logit"}, the option \code{use.mclogit = TRUE} can be specified to request that \pkgfun{mclogit}{mblogit} from the \CRANpkg{mclogit} package is used instead, which can be faster and is recommended.

%For continuous treatments in the presence of missing data with \code{missing = "saem"}, additional arguments are passed to \code{\link[misaem]{miss.lm}} and \code{\link[misaem]{predict.miss.lm}}.

}
\section{Additional Outputs}{
\describe{
\item{\code{obj}}{
When \code{include.obj = TRUE}, the (generalized) propensity score model fit. For binary treatments, the output of the call to \fun{glm}. For ordinal treatments, the output of the call to \pkgfun{MASS}{polr}. For multinomial treatments with \code{link = "logit"} or \code{"probit"} and \code{use.mlogit = TRUE}, the output of the call to \pkgfun2{mlogit}{mlogit}{mlogit::mlogit}. For multinomial treatments with \code{use.mlogit = FALSE}, a list of the \code{glm()} fits. For multinomial treatments with \code{link = "br.logit"}, the output of the call to \pkgfun2{brglm2}{brmultinom}{brglm2::brmultinom}. For multinomial treatments with \code{link = "bayes.probit"}, the output of the call to \pkgfun2{MNP}{mnp}{MNP::mnp}. For continuous treatments, the output of the call to \code{glm()} for the predicted values in the denominator density.
}
}
}
\references{
\strong{Binary treatments}

 - \code{estimand = "ATO"}

Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates via propensity score weighting. Journal of the American Statistical Association, 113(521), 390–400. \doi{10.1080/01621459.2016.1260466}

 - \code{estimand = "ATM"}

Li, L., & Greene, T. (2013). A Weighting Analogue to Pair Matching in Propensity Score Analysis. The International Journal of Biostatistics, 9(2). \doi{10.1515/ijb-2012-0030}

 - \code{estimand = "ATOS"}

Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009). Dealing with limited overlap in estimation of average treatment effects. Biometrika, 96(1), 187–199. \doi{10.1093/biomet/asn055}

 - Other estimands

Austin, P. C. (2011). An Introduction to Propensity Score Methods for Reducing the Effects of Confounding in Observational Studies. Multivariate Behavioral Research, 46(3), 399–424. \doi{10.1080/00273171.2011.568786}

 - Marginal mean weighting through stratification

Hong, G. (2010). Marginal mean weighting through stratification: Adjustment for selection bias in multilevel data. Journal of Educational and Behavioral Statistics, 35(5), 499–531. \doi{10.3102/1076998609359785}

- Bias-reduced logistic regression

See references for the \pkg{brglm2} \link[brglm2:brglm2]{package}.

%- SAEM logistic regression for missing data

%Jiang, W., Josse, J., & Lavielle, M. (2019). Logistic regression with missing covariates — Parameter estimation, model selection and prediction within a joint-modeling framework. Computational Statistics & Data Analysis, 106907. \doi{10.1016/j.csda.2019.106907}

\strong{Multinomial Treatments}

 - \code{estimand = "ATO"}

Li, F., & Li, F. (2019). Propensity score weighting for causal inference with multiple treatments. The Annals of Applied Statistics, 13(4), 2389–2415. \doi{10.1214/19-AOAS1282}

 - \code{estimand = "ATM"}

Yoshida, K., Hernández-Díaz, S., Solomon, D. H., Jackson, J. W., Gagne, J. J., Glynn, R. J., & Franklin, J. M. (2017). Matching weights to simultaneously compare three treatment groups: Comparison to three-way matching. Epidemiology (Cambridge, Mass.), 28(3), 387–395. \doi{10.1097/EDE.0000000000000627}

 - Other estimands

McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. Statistics in Medicine, 32(19), 3388–3414. \doi{10.1002/sim.5753}

 - Marginal mean weighting through stratification

Hong, G. (2012). Marginal mean weighting through stratification: A generalized method for evaluating multivalued and multiple treatments with nonexperimental data. Psychological Methods, 17(1), 44–60. \doi{10.1037/a0024918}

\strong{Continuous treatments}

Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural Models and Causal Inference in Epidemiology. Epidemiology, 11(5), 550–560.

- Using non-normal conditional densities

Naimi, A. I., Moodie, E. E. M., Auger, N., & Kaufman, J. S. (2014). Constructing Inverse Probability Weights for Continuous Exposures: A Comparison of Methods. Epidemiology, 25(2), 292–299. \doi{10.1097/EDE.0000000000000053}

}
\seealso{
\fun{weightit}, \fun{weightitMSM}, \fun{get_w_from_ps}
}
\examples{
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
}
