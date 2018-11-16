\name{weightit}
\alias{weightit}
\alias{WeightIt}
\alias{print.weightit}
\title{Generate Balancing Weights}
\usage{
weightit(formula,
         data = NULL,
         method = "ps",
         estimand = "ATE",
         stabilize = FALSE,
         focal = NULL,
         exact = NULL,
         s.weights = NULL,
         ps = NULL,
         moments = 1,
         int = FALSE,
         verbose = FALSE,
         ...)

\method{print}{weightit}(x, ...)

}
\description{
\code{weightit()} allows for the easy generation of balancing weights using a variety of available methods for binary, continuous, and multinomial treatments. Many of these methods exist in other packages, which \code{weightit()} calls; these packages must be installed to use the desired method. Also included are \code{print} and \code{summary} methods for examining the output.
}
\arguments{
  \item{formula}{
a formula with a treatment variable on the left hand side and the covariates to be balanced on the right hand side. See \code{\link{glm}} for more details. Interactions and functions of covariates are allowed.
}
  \item{data}{
an optional data set in the form of a data frame that contains the variables in \code{formula}.
}
  \item{method}{
a string of length 1 containing the name of the method that will be used to estimate weights. See Details below for allowable options. The default is \code{"ps"}.
}
  \item{estimand}{
the desired estimand. For binary treatments, can be "ATE", "ATT", "ATC", and, for some methods, "ATO" or "ATM". For multinomial treatments, can be "ATE" or "ATT". The default for both is "ATE". This argument is ignored for continuous treatments. See Details for more information.
}
  \item{stabilize}{
\code{logical}; whether or not to stabilize the weights. For the methods that involve estimating propensity scores, this involves multiplying each unit's weight by the sum of the weights in that unit's treatment group. For the \code{"ebal"} method, this involves using \code{ebalance.trim()} to reduce the variance of the weights. Default is \code{FALSE}.
}
  \item{focal}{
when multinomial treatments are used and the "ATT" is requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group.
}
  \item{exact}{
a vector or the names of variables in data for which weighting is to be done within catgories. For example, if \code{exact = "gender"}, weights will be generated seperately within each level of the variable \code{"gender"}.
}
  \item{s.weights}{
A vector of sampling weights or the name of a variable in \code{data} that contains sampling weights. These can also be matching weights if weighting is to be used on matched data.
}
  \item{ps}{
A vector of propensity scores or the name of a variable in \code{data} containing propensity scores. If not \code{NULL}, \code{method} is ignored, and the propensity scores will be used to create weights. \code{formula} must include the treatment variable in \code{data}, but the listed covariates will play no role in the weight estimation.
}
  \item{moments}{
\code{numeric}; for entropy balancing, empirical balancing callibration weights, optimization-based weights, the greatest moment of the covariate distribution to be balanced. For example, if \code{moments = 3}, for all non-categorical covariates, the mean, second moment (variance), and third moments (skew) of the covariates will be balanced. This argument is ignored for other methods; to balance powers of the covariates, appropriate functions must be entered in \code{formula}.
}
  \item{int}{
\code{logical}; for entropy balancing, empirical balancing callibration weights, and soptimization-based weights, whether first-order interactons of the covariates are to be balanced (essentially balancing the covariances between covariates). This argument is ignored for other methods; to balance interactions between the variables, appropriate functions must be entered in \code{formula}.
}
  \item{verbose}{
whether to print additional information output by the fitting function.
}
  \item{...}{
other arguments for functions called by \code{weightit} that control aspects of fitting that are not covered by the above arguments. See Details.
}
  \item{x}{
a \code{weightit} object; the output of a call to \code{weightit()}.
}

}
\value{
  A \code{weightit} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{treat}{The values of the treatment variable.}
  \item{covs}{The covariates used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
  \item{estimand}{The estimand requested.}
  \item{method}{The weight estimation method specified.}
  \item{ps}{The estimated or provided propensity scores.}
  \item{s.weights}{The provided sampling weights.}
  \item{focal}{The focal variable if the ATT was requested with a multinomial treatment.}
}
\details{
The primary purpose of \code{weightit()} is as a dispatcher to other functions in other packages that perform the estimation of balancing weights. These functions are identified by a name, which is used in \code{method} to request them. Each method has some slight distinctions in how it is called, but in general, simply entering the method will cause \code{weightit()} to generate the weights correctly using the function. To use each method, the package containing the function must be installed, or else an error will appear. Below are the methods allowed and some information about each.
\describe{
\item{\code{"ps"}}{\strong{Propensity score weighting using GLM.} For binary treatments, this method estimates the propensity scores using \code{glm()}. An additional argument is \code{link}, which uses the same options as \code{link} in \code{\link{family}}. The default link is "logit", but others, including "probit", are allowed. The weights for the ATE, ATT, and ATC are computed from the estimated propensity scores using the standard formulas, the weights for the ATO are computed as in Li & Li (2018), and the weights for the ATM (i.e., average treatment effect in the equivalent sample "pair-matched" with calipers) are computed as in Yoshida et al (2017).

For multinomial treatments, the propensity scores are estimated using multinomial regression from one of two functions depending on the requested link: for logit (\code{"logit"}) and probit (\code{"probit"}) links, \code{mlogit()} from the \pkg{mlogit} package is used, and for the Bayesian probit  (\code{"bayes.probit"}) link, \code{MNP()} from the \pkg{MNP} package is used. If \code{mlogit} in not installed, a series of binary regressions using \code{glm()} will be run instead, with estimated propensities normalized to sum to 1. These are the only three links allowed for multinomial treatments at this time. (These methods can fail to converge, yielding errors that may seem foreign.)

For continuous treatments, the generalized propensity score is estimated using linear regression with a normal density, but other families and links are allowed, such as poisson for count treatments, using the \code{family} and \code{link} arguments. In addition, kernel density estimation can be used instead of assuming a normal density for the numerator and denominator of the generalized propensity score by setting \code{use.kernel = TRUE}. Other arguments to \code{\link{density}} can be specified to refine the density estimation parameters. \code{plot = TRUE} can be specified to plot the density for the numerator and denominator, which can be helpful in diagnosing extreme weights.

For all treatment types except multinomial treatments with a Bayesian probit link, sampling weights are supported, but a warning message from \code{glm()} may appear.}

\item{\code{"gbm"}}{\strong{Propensity score weighting using generalized boosted modeling.} This method, which can also be requested as \code{"gbr"} or \code{"twang"}, uses functions from the \pkg{twang} package to perform generalized boosted modeling to estimate propensity scores that yield balance on the requested covariates. For binary treatments, \code{ps()} is used, and the ATE, ATT, and ATC can be requested. For multinomial treatments, \code{mnps()} is used, and the ATE or ATT can be requested. For both, the \code{weightit()} argument \code{s.weights} corresponds to the \code{ps()} and \code{mnps()} argument \code{sampw}. The \code{weightit()} argument \code{focal} corresponds to the \code{mnps()} argument \code{treatATT}. For both, a single stop method must be supplied to \code{stop.method}; only one can be entered at a time. The other arguments to \code{ps()} and \code{mnps} can be specified in the call to \code{weightit()}. See \code{\link[twang]{ps}} and \code{\link[twang]{mnps}} for details.

For continuous treatments, the generalized propensity score is estimated using \code{ps.cont()}, which exists in \code{WeightIt}. Balance is optimized by a function of the correlation between the treatment and covariates in the weighted sample. See \code{\link{ps.cont}} for details.}

\item{\code{"cbps"}}{\strong{Covariate Balancing Propensity Score weighting.} This method uses the \code{CBPS()} function from the \pkg{CBPS} package to estimate propensity scores and weights. It works with binary, multinomial, and continuous treatments. For binary treatments, the ATE, ATT, and ATC can be requested. For multinomial treatments, only the ATE can be requested. The \code{weightit()} argument \code{s.weights} corresponds to the \code{CBPS()} argument \code{sampling.weights}. \code{CBPS()} can fit either an over-identified model or a model that only contains covariate balancing conditions; this option is typically specified with the \code{method} argument to \code{CBPS()}, but because this argument is already used in \code{weightit()}, a new argument, \code{over}, can be specified. \code{over = FALSE} in \code{weightit()} is equivalent to \code{method = "exact"} in \code{CBPS()}. The other arguments to \code{CBPS()} can be specified in the call to \code{weightit()}. See \code{\link[CBPS]{CBPS}} for details.}

\item{\code{"npcbps"}}{\strong{Non-parametric Covariate Balancing Propensity Score weighting.} This method uses the \code{npCBPS()} function from the \pkg{CBPS} package to estimate weights. It works with binary, multinomial, and continuous treatments. For binary and multinomial treatments, only the ATE can be requested. Sampling weights are not supported. The other arguments to \code{npCBPS()} can be specified in the call to \code{weightit()}. See \code{\link[CBPS]{npCBPS}} for details.}

\item{\code{"ebal"}}{\strong{Entropy balancing.} This method uses the \code{ebalance()} function from the \pkg{ebal} package to estimate weights. It works with binary and multinomial treatments. For binary treatments, the ATE, ATT, and ATC can be requested. For multinomial treatments, the ATE and ATT can be requested. If the ATT is requetsed with a mutlinomial treatment, one treatment level must be entered to \code{focal} to serve as the "treated". Sampling weights are supported and are automatically entered into \code{base.weight} in \code{ebal()}. When \code{stabilize = TRUE}, \code{ebalance.trim()} is used to trim and reduce the variance of the weights. The other arguments to \code{ebalance()} can be specified in the call to \code{weightit()}. See \code{\link[ebal]{ebalance}} for details.}

%\item{\code{"sbw"}}{\strong{Stable balancing weights.} This method uses the \code{sbw()} function from the \pkg{sbw} package to estimate weights. For binary treatments, the ATE, ATT, and ATC can be requested. For multinomial treatments, the ATE and ATT can be requested. If the ATT is requested with a multinomial treatment, one treatment level must be entered to \code{focal} to serve as the "treated". Sampling weights are supported for all estimands. The other arguments to \code{sbw()} can be specified in the call to \code{weightit()}. See \code{sbw} for details. There are some slight difference between the default options in \code{weightit()} and \code{sbw()}; importantly, in \code{weightit()} when \code{bal_tols_sd} is \code{TRUE} (the default), the standardized mean difference is not used for categorical variables, and the demoninator of the standardized mean difference corresponds to the standard deviation of the target group (e.g., for the ATT, the denominator is the standard deviation of the treated group).}

\item{\code{"ebcw"}}{\strong{Empirical balancing calibration weighting.} This method uses the \code{ATE()} function from the \pkg{ATE} package to estimate weights. It works with binary and multinomial treatments. For binary treatments, the ATE, ATT, and ATC can be requested. For multinomial treatments, the ATE and ATT can be requested, and an argument to \code{focal} must be specified for the ATT. Sampling weights are supported for all estimands. The other arguments to \code{ATE()} can be specified in the call to \code{weightit()}. See \code{\link[ATE]{ATE}} for details.}

\item{\code{"optweight"}}{\strong{Optimization-based weights.} This method uses the \code{optweight()} function from the \pkg{optweight} package to estimate weights. It works with binary, multinomial, continuous, and longitudinal treatments. For binary treatments, the ATE, ATT, and ATC can be requested. For multinomial treatments, the ATE and ATT can be requested, and an argument to \code{focal} must be specified for the ATT. Sampling weights are supported for all estimands. The other arguments to \code{optweight()} can be specified in the call to \code{optweight()}. Note that \code{weightit()} processes the formula containing the covariates differently from how \code{optweight()} processes it; importantly, when using \code{tols} to specify balance thresholds in \code{weightit()}, a threshold must be set for each level of factor variables. This problem is aleviated if only one tolerance value is set for all covariates. When using \code{moments} or \code{int}, it is recommended to supply only one tolerance value. See \code{\link[optweight]{optweight}} for details on other options and arguments. Because \code{optweight()} offers very similar functionality to \code{weightit()}, it is recommended to use that function directly. The syntax and output are almost identical.}

\item{\code{"super"}}{\strong{Propensity score weighting using SuperLearner.} This method uses the \code{SuperLearner()} function from the \pkg{SuperLearner} package to estimate propensity scores, which are used to compute weights. It works with binary, multinomial, and continuous treatments. For binary treatments, the ATE, ATT, ATC, ATO and ATM can be requested. For multinomial treatments, the ATE, ATT, ATO, and ATM can be requested. An argument to \code{focal} must be specified for the ATT. Sampling weights are supported for all estimands using \code{s.weights}, which replaces the \code{obsWeights} argument in \code{SuperLearner()}. The other arguments to \code{SuperLearner()} can be specified in the call to \code{weightit()}. An argument to \code{SL.library} must be supplied; use \code{SuperLearner::listWrappers()} to see all options (\pkg{SuperLearner} may to be loaded for this to work). To pass an argument to \code{method} in \code{SuperLearner()}, an argument to \code{SL.method} must be used instead (because \code{method} is already used in \code{weightit()}). See \code{\link[SuperLearner]{SuperLearner}} for details. All arguments and options that work for \code{method = "ps"} work here as well.}
}
}

\author{
Noah Greifer

The code used for using generalized boosted modeling with continuous treatments was adapted from that which appeared in Zhu, Coffman, & Ghosh (2015).

}

\references{

\strong{Binary treatments}

\code{method = "ps"}

 - \code{estimand = "ATO"}

Li, F., Morgan, K. L., & Zaslavsky, A. M. (2016). Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association, 0(ja), 0–0.

 - \code{estimand = "ATM"}

Li, L., & Greene, T. (2013). A Weighting Analogue to Pair Matching in Propensity Score Analysis. The International Journal of Biostatistics, 9(2). \doi{10.1515/ijb-2012-0030}

 - Other estimands

Austin, P. C. (2011). An Introduction to Propensity Score Methods for Reducing the Effects of Confounding in Observational Studies. Multivariate Behavioral Research, 46(3), 399–424. \doi{10.1080/00273171.2011.568786}

\code{method = "gbm"}

McCaffrey, D. F., Ridgeway, G., & Morral, A. R. (2004). Propensity Score Estimation With Boosted Regression for Evaluating Causal Effects in Observational Studies. Psychological Methods, 9(4), 403–425. \doi{10.1037/1082-989X.9.4.403}

\code{method = "cbps"}

Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 76(1), 243–263.

\code{method = "npcbps"}

Fong, Christian, Chad Hazlett, and Kosuke Imai (2015). Parametric and Nonparametric Covariate Balancing Propensity Score for General Treatment Regimes. Unpublished Manuscript. <\url{http://imai.princeton.edu/research/files/CBGPS.pdf}>

\code{method = "ebal"}

Hainmueller, J. (2012). Entropy Balancing for Causal Effects: A Multivariate Reweighting Method to Produce Balanced Samples in Observational Studies. Political Analysis, 20(1), 25–46. \doi{10.1093/pan/mpr025}

\code{method = "ebcw"}

Chan, K. C. G., Yam, S. C. P., & Zhang, Z. (2016). Globally efficient non-parametric inference of average treatment effects by empirical balancing calibration weighting. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 78(3), 673–700. \doi{10.1111/rssb.12129}

\code{method = "optweight"}

Zubizarreta, J. R. (2015). Stable Weights that Balance Covariates for Estimation With Incomplete Outcome Data. Journal of the American Statistical Association, 110(511), 910–922. \doi{10.1080/01621459.2015.1023805}

\code{method = "super"}

Pirracchio, R., Petersen, M. L., & van der Laan, M. (2015). Improving Propensity Score Estimators’ Robustness to Model Misspecification Using Super Learner. American Journal of Epidemiology, 181(2), 108–119. \doi{10.1093/aje/kwu253}

\strong{Multinomial Treatments}

\code{method = "ps"}

 - \code{estimand = "ATO"}

Li, F., & Li, F. (2018). Propensity Score Weighting for Causal Inference with Multi-valued Treatments. ArXiv:1808.05339 [Stat]. Retrieved from http://arxiv.org/abs/1808.05339


 - \code{estimand = "ATM"}

Yoshida, K., Hernández-Díaz, S., Solomon, D. H., Jackson, J. W., Gagne, J. J., Glynn, R. J., & Franklin, J. M. (2017). Matching weights to simultaneously compare three treatment groups: Comparison to three-way matching. Epidemiology (Cambridge, Mass.), 28(3), 387–395. \doi{10.1097/EDE.0000000000000627}


 - Other estimands

McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. Statistics in Medicine, 32(19), 3388–3414. \doi{10.1002/sim.5753}

\code{method = "gbm"}

McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. Statistics in Medicine, 32(19), 3388–3414. \doi{10.1002/sim.5753}

\strong{Continuous treatments}

\code{method = "ps"}

Robins, J. M., Hernán, M. Á., & Brumback, B. (2000). Marginal Structural Models and Causal Inference in Epidemiology. Epidemiology, 11(5), 550–560.

\code{method = "gbm"}

Zhu, Y., Coffman, D. L., & Ghosh, D. (2015). A Boosting Algorithm for Estimating Generalized Propensity Scores with Continuous Treatments. Journal of Causal Inference, 3(1).

\code{method = "cbps"}

Fong, C., Hazlett, C. , and Imai, K. (2015). Parametric and Nonparametric Covariate Balancing Propensity Score for General Treatment Regimes. Unpublished Manuscript. <\url{http://imai.princeton.edu/research/files/CBGPS.pdf}>

}

\note{
\code{method = "sbw"} has been deprecated in favor of \code{optweight}, which does the same thing.
}
\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "ps", estimand = "ATT"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATE"))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = "cbps", over = FALSE))
summary(W3)
bal.tab(W3)
}
