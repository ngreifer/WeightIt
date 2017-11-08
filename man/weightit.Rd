\name{weightit}
\alias{weightit}
\alias{print.weightit}
\title{Generate Balancing Weights}
\usage{
weightit(formula, data, method, estimand = "ATE",
         stabilize = FALSE, focal = NULL, exact = NULL,
         s.weights = NULL, ps = NULL, verbose = FALSE,
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
a data set in the form of a data frame that contains the variables in \code{formula}
}
  \item{method}{
a string of length 1 containing the name of the method that will be used to estimate weights. See Details below for allowable options. The default is \code{"ps"}.
}
  \item{estimand}{
the desired estimand. For binary treatments, can be "ATE", "ATT", "ATC", and, for some methods, "ATO". For multinomial treatments, can be "ATE" or "ATT". The default for both is "ATE". This argument is ignored for continuous treatments.
}
  \item{stabilize}{
\code{logical}; whether or not to stabilize the weights. For the \code{"ps"} method, this involves multiplying each unit's weight by the marginal probability of being in its treatment group. For the \code{"ebal"} method, this involves using \code{ebalance.trim()} to reduce the variance of the weights. Default is \code{FALSE}.
}
  \item{focal}{
when multinomial treatments are used and the "ATT" is requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group.
}
  \item{exact}{
a vector or the names of variables in data for which weighting is to be done within catgories. For example, if \code{exact = "gender"}, weights will be generated seperately within each level of the variable \code{"gender"}.
}
  \item{s.weights}{
A vector of sampling weights or the name of a variable in \code{data} that contains sampling weights. These are ignored for some methods. These can also be matching weights if weighting is to be used on matched data.
}
  \item{ps}{
A vector of propensity scores or the name of a variable in \code{data} containing propensity scores. If not \code{NULL}, \code{method} is ignored, and the propensity scores will be used to create weights. \code{formula} must include the treatment variable in \code{data}, but the listed covariates will play no role in the weight estimation.
}
  \item{verbose}{
whether to print additional information output but the fitting function.
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
  \item{data}{The data.frame originally entered to \code{weightit()}.}
  \item{estimand}{The estimand requested.}
  \item{method}{The weight estimation method specified.}
  \item{ps}{The estimated or provided propensity scores.}
  \item{s.weights}{The provided sampling weights.}
  \item{discarded}{Unused.}
  \item{treat.type}{The type of treatment: binary, continuous, or multinomial (\code{"multi"}).}
  \item{focal}{The focal variable if the ATT was requested with a multinomial treatment.}
}
\details{
The primary purpose of \code{weightit()} is as a dispatcher to other functions in other packages that perform the estimation of balancing weights. These functions are identified by a name, which is used in \code{method} to request them. Each method has some slight distinctions in how it is called, but in general, simply entering the method will cause \code{weightit()} to generate the weights correctly using the function. To use each method, the package containing the function must be installed, or else an error will appear. Below are the methods allowed and some information about each.
\describe{
\item{\code{"ps"}}{\strong{Propensity score weighting using GLM.} For binary treatments, this method estimates the propensity scores using \code{glm()}. An additional argument is \code{link}, which uses the same options as \code{link} in \code{\link{family}}. The default link is "logit", but others, including "probit", are allowed. The weights for the ATE, ATT, and ATC are computed from the estimated propensity scores using the standard formulas, and the weights for the ATO are computed as in Li, Morgan, & Zaslavsky (2016). For multinomial treatments, the propensity scores are estimated using multinomial regression from one of two functions depending on the requested link: for logit (\code{"logit"}) and probit (\code{"probit"}) links, \code{mlogit()} from the \pkg{mlogit} package is used, and for the Bayesian probit  (\code{"bayes.probit"}) link, \code{MNP()} from the \pkg{MNP} package is used. These are the only three links allowed for multinomial treatmets at this time. (These methods can fail to converge, yielding errors that may seem foreign.) For continuous treatments, the generalized propensity score is estimated using linear regression with a normal density, but other families and links are allowed, such as poisson for count data, using the \code{family} and \code{link} arguments. An additional argument, \code{num.formula}, may be specified, containing the stabilization variables on the right hand side. For all treatment types except multinomial treatments with a Bayesian probit link, sampling weights are supported, but a warning message from \code{glm()} may appear.}

\item{\code{"gbm"}}{\strong{Propensity score weighting using generalized boosted modeling.} This method, which can also be requested as \code{"gbr"} or \code{"twang"}, uses functions from the \pkg{twang} package to perform generalized boosted modeling to estimate propensity scores that yield balance on the requested covariates. For binary treatments, \code{ps()} is used, and the ATE, ATT, and ATC can be requested. For multinomial treatments, \code{mnps()} is used, and the ATE or ATT can be requested. For both, the \code{weightit()} argument \code{s.weights} corresponds to the \code{ps()} and \code{mnps()} argument \code{sampw}. The \code{weightit()} argument \code{focal} corresponds to the \code{mnps()} argument \code{treatATT}. For both, a single stop method must be supplied to \code{stop.method}; only one can be entered at a time. The other arguments to \code{ps()} and \code{mnps} can be specified in the call to \code{weightit()}. See \code{\link[twang]{ps}} and \code{\link[twang]{mnps}} for details.}

\item{\code{"cbps"}}{\strong{Covariate Balancing Propensity Score weighting.} This method uses the \code{CBPS()} function from the \pkg{CBPS} package to estimate propensity scores and weights. It works with binary, multinomial, and continuous treatments. For binary treatments, the ATE, ATT, and ATC can be requested. For multinomial treatments, only the ATE can be requested. The \code{weightit()} argument \code{s.weights} corresponds to the \code{CBPS()} argument \code{sampling.weights}. \code{CBPS()} can fit either an over-identified model or a model that only contains covariate balancing conditions; this option is typically specified with the \code{method} argument to \code{CBPS()}, but because this argument is already used in \code{weightit()}, a new argument, \code{over}, can be specified. \code{over = FALSE} in \code{weightit()} is equivalent to \code{method = "exact"} in \code{CBPS()}. The other arguments to \code{CBPS()} can be specified in the call to \code{weightit()}. See \code{\link[CBPS]{CBPS}} for details.}

\item{\code{"npcbps"}}{\strong{Non-parametric Covariate Balancing Propensity Score weighting.} This method uses the \code{npCBPS()} function from the \pkg{CBPS} package to estimate weights. It works with binary, multinomial, and continuous treatments. For binary and multinomial treatments, only the ATE can be requested. Sampling weights are not supported. The other arguments to \code{npCBPS()} can be specified in the call to \code{weightit()}. See \code{\link[CBPS]{npCBPS}} for details.}

\item{\code{"ebal"}}{\strong{Entropy balancing.} This method uses the \code{ebalance()} function from the \pkg{ebal} package to estimate weights. It works with binary and multinomial treatments. For binary treatments, the ATE, ATT, and ATC can be requested. For multinomial treatments, the ATE and ATT can be requested. If the ATT is requetsed with a mutlinomial treatment, one treatment level must be entered to \code{focal} to serve as the "treated". Sampling weights are supported and are automatically entered into \code{base.weight} in \code{ebal()}. When \code{stabilize = TRUE}, \code{ebalance.trim()} is used to trim and reduce the variance of the weights. The other arguments to \code{ebalance()} can be specified in the call to \code{weightit()}. See \code{\link[ebal]{ebalance}} for details.}

\item{\code{"ebcw"}}{\strong{Empirical balancing calibration weighting.} This method uses the \code{ATE()} function from the \pkg{ATE} package to estimate weights. It works with binary and multinomial treatments. For binary treatments, the ATE, ATT, and ATC can be requested. For multinomial treatments, only the ATE can be requested. Sampling weights are not supported. The other arguments to \code{ATE()} can be specified in the call to \code{weightit()}. See \code{\link[ATE]{ATE}} for details.}
}
}

\author{
Noah Greifer \email{noah@unc.edu}
}

\references{
Li, F., Morgan, K. L., & Zaslavsky, A. M. (2016). Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association, 0(ja), 0–0.

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
