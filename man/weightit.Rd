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
         by = NULL,
         s.weights = NULL,
         ps = NULL,
         moments = 1,
         int = FALSE,
         verbose = FALSE,
         include.obj = FALSE,
         ...)

\method{print}{weightit}(x, ...)

}
\description{
\code{weightit} allows for the easy generation of balancing weights using a variety of available methods for binary, continuous, and multinomial treatments. Many of these methods exist in other packages, which \code{weightit} calls; these packages must be installed to use the desired method. Also included are \code{print} and \code{summary} methods for examining the output.
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
\code{logical}; whether or not to stabilize the weights. For the methods that involve estimating propensity scores, this involves multiplying each unit's weight by the sum of the weights in that unit's treatment group. For the \code{"ebal"} method, this involves using \code{ebalance.trim} to reduce the variance of the weights. Default is \code{FALSE}.
}
  \item{focal}{
when multinomial treatments are used and the "ATT" is requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group.
}
  \item{by}{
a vector or the names of variables in data for which weighting is to be done within categories. For example, if \code{by = "gender"}, weights will be generated separately within each level of the variable \code{"gender"}. The argument used to be called \code{exact}, which will still work but with a message.
}
  \item{s.weights}{
A vector of sampling weights or the name of a variable in \code{data} that contains sampling weights. These can also be matching weights if weighting is to be used on matched data.
}
  \item{ps}{
A vector of propensity scores or the name of a variable in \code{data} containing propensity scores. If not \code{NULL}, \code{method} is ignored, and the propensity scores will be used to create weights. \code{formula} must include the treatment variable in \code{data}, but the listed covariates will play no role in the weight estimation.
}
  \item{moments}{
\code{numeric}; for entropy balancing, empirical balancing calibration weights, optimization-based weights, the greatest moment of the covariate distribution to be balanced. For example, if \code{moments = 3}, for all non-categorical covariates, the mean, second moment (variance), and third moments (skew) of the covariates will be balanced. This argument is ignored for other methods; to balance powers of the covariates, appropriate functions must be entered in \code{formula}.
}
  \item{int}{
\code{logical}; for entropy balancing, empirical balancing calibration weights, and optimization-based weights, whether first-order interactions of the covariates are to be balanced (essentially balancing the covariances between covariates). This argument is ignored for other methods; to balance interactions between the variables, appropriate functions must be entered in \code{formula}.
}
  \item{verbose}{
whether to print additional information output by the fitting function.
}
  \item{include.obj}{
whether to include in the output any fit objects created in the process of estimating the weights. For example, with \code{method = "ps"}, the \code{glm} objects containing the propensity score model will be included. See Details for information on what object will be included if \code{TRUE}.
}
  \item{...}{
other arguments for functions called by \code{weightit} that control aspects of fitting that are not covered by the above arguments. See Details.
}
  \item{x}{
a \code{weightit} object; the output of a call to \code{weightit}.
}

}
\value{
  A \code{weightit} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{treat}{The values of the treatment variable.}
  \item{covs}{The covariates used in the fitting. Only includes the raw covariates, which may have been altered in the fitting process.}
  \item{estimand}{The estimand requested.}
  \item{method}{The weight estimation method specified.}
  \item{ps}{The estimated or provided propensity scores. Estimated propensity scores are returned for binary treatments and only when \code{method} is \code{"ps"}, \code{"gbm"}, \code{"cbps"}, or \code{"super"}.}
  \item{s.weights}{The provided sampling weights.}
  \item{focal}{The focal variable if the ATT was requested with a multinomial treatment.}
  \item{obj}{When \code{include.obj = TRUE}, the fit object.}
}
\details{
The primary purpose of \code{weightit} is as a dispatcher to other functions in other packages that perform the estimation of balancing weights. These functions are identified by a name, which is used in \code{method} to request them. Each method has some slight distinctions in how it is called, but in general, simply entering the method will cause \code{weightit} to generate the weights correctly using the function. To use each method, the package containing the function must be installed, or else an error will appear. Below are the methods allowed and links to pages containing more information about them, including additional arguments and outputs (e.g., when \code{include.obj = TRUE}).
\itemize{
\item{\code{\link[=method_ps]{"ps"}} - Propensity score weighting using generalized linear models.}
\item{\code{\link[=method_gbm]{"gbm"}} - Propensity score weighting using generalized boosted modeling.}
\item{\code{\link[=method_cbps]{"cbps"}} - Covariate Balancing Propensity Score weighting.}
\item{\code{\link[=method_npcbps]{"npcbps"}} - Non-parametric Covariate Balancing Propensity Score weighting.}
\item{\code{\link[=method_ebal]{"ebal"}} - Entropy balancing.}
\item{\code{\link[=method_ebcw]{"ebcw"}} - Empirical balancing calibration weighting.}
\item{\code{\link[=method_optweight]{"optweight"}} - Optimization-based weighting.}
\item{\code{\link[=method_super]{"super"}} - Propensity score weighting using SuperLearner.}
}
}

\author{
Noah Greifer
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
