\name{weightit}
\alias{weightit}
\alias{WeightIt}
\alias{print.weightit}
\title{Generate Balancing Weights}
\description{
\code{weightit()} allows for the easy generation of balancing weights using a variety of available methods for binary, continuous, and multi-category treatments. Many of these methods exist in other packages, which \code{weightit()} calls; these packages must be installed to use the desired method. Also included are \code{print()} and \code{summary()} methods for examining the output.
}
\usage{
weightit(formula,
         data = NULL,
         method = "glm",
         estimand = "ATE",
         stabilize = FALSE,
         focal = NULL,
         by = NULL,
         s.weights = NULL,
         ps = NULL,
         moments = NULL,
         int = FALSE,
         subclass = NULL,
         missing = NULL,
         verbose = FALSE,
         include.obj = FALSE,
         ...)

\method{print}{weightit}(x, ...)
}
\arguments{
  \item{formula}{
a formula with a treatment variable on the left hand side and the covariates to be balanced on the right hand side. See \fun{glm} for more details. Interactions and functions of covariates are allowed.
}
  \item{data}{
an optional data set in the form of a data frame that contains the variables in \code{formula}.
}
  \item{method}{
a string of length 1 containing the name of the method that will be used to estimate weights. See Details below for allowable options. The default is \code{"glm"} for propensity score weighting using a generalized linear model to estimate the propensity score.
}
  \item{estimand}{
the desired estimand. For binary and multi-category treatments, can be \code{"ATE"}, \code{"ATT"}, \code{"ATC"}, and, for some methods, \code{"ATO"}, \code{"ATM"}, or \code{"ATOS"}. The default for both is \code{"ATE"}. This argument is ignored for continuous treatments. See the individual pages for each method for more information on which estimands are allowed with each method and what literature to read to interpret these estimands.
}
  \item{stabilize}{
\code{logical}; whether or not to stabilize the weights. For the methods that involve estimating propensity scores, this involves multiplying each unit's weight by the proportion of units in their treatment group. Default is \code{FALSE}.
}
  \item{focal}{
when multi-category treatments are used and ATT weights are requested, which group to consider the "treated" or focal group. This group will not be weighted, and the other groups will be weighted to be more like the focal group. If specified, \code{estimand} will automatically be set to \code{"ATT"}.
}
  \item{by}{
a string containing the name of the variable in \code{data} for which weighting is to be done within categories or a one-sided formula with the stratifying variable on the right-hand side. For example, if \code{by = "gender"} or \code{by = ~gender}, weights will be generated separately within each level of the variable \code{"gender"}. (The argument used to be called \code{exact}, which will still work but with a message.) Only one \code{by} variable is allowed; to stratify by multiply variables simultaneously, create a new variable that is a full cross of those variables using \fun{interaction}.
}
  \item{s.weights}{
A vector of sampling weights or the name of a variable in \code{data} that contains sampling weights. These can also be matching weights if weighting is to be used on matched data. See the individual pages for each method for information on whether sampling weights can be supplied.
}
  \item{ps}{
A vector of propensity scores or the name of a variable in \code{data} containing propensity scores. If not \code{NULL}, \code{method} is ignored, and the propensity scores will be used to create weights. \code{formula} must include the treatment variable in \code{data}, but the listed covariates will play no role in the weight estimation. Using \code{ps} is similar to calling \fun{get_w_from_ps} directly, but produces a full \code{weightit} object rather than just producing weights.
}
  \item{moments}{
\code{numeric}; for some methods, the greatest power of each covariate to be balanced. For example, if \code{moments = 3}, for each non-categorical covariate, the covariate, its square, and its cube will be balanced. This argument is ignored for other methods; to balance powers of the covariates, appropriate functions must be entered in \code{formula}. See the individual pages for each method for information on whether they accept \code{moments}.
}
  \item{int}{
\code{logical}; for some methods, whether first-order interactions of the covariates are to be balanced. This argument is ignored for other methods; to balance interactions between the variables, appropriate functions must be entered in \code{formula}. See the individual pages for each method for information on whether they accept \code{int}.
}
  \item{subclass}{
\code{numeric}; the number of subclasses to use for computing weights using marginal mean weighting with subclasses (MMWS). If \code{NULL}, standard inverse probability weights (and their extensions) will be computed; if a number greater than 1, subclasses will be formed and weights will be computed based on subclass membership. Attempting to set a non-\code{NULL} value for methods that don't compute a propensity score will result in an error; see each method's help page for information on whether MMWS weights are compatible with the method. See \fun{get_w_from_ps} for details and references.
}
  \item{missing}{
\code{character}; how missing data should be handled. The options and defaults depend on the \code{method} used. Ignored if no missing data is present. It should be noted that multiple imputation outperforms all available missingness methods available in \code{weightit()} and should probably be used instead. Consider the \CRANpkg{MatchThem} package for the use of \code{weightit()} with multiply imputed data.
}
  \item{verbose}{
\code{logical}; whether to print additional information output by the fitting function.
}
  \item{include.obj}{
\code{logical}; whether to include in the output any fit objects created in the process of estimating the weights. For example, with \code{method = "glm"}, the \code{glm} objects containing the propensity score model will be included. See the individual pages for each method for information on what object will be included if \code{TRUE}.
}
  \item{...}{
other arguments for functions called by \code{weightit()} that control aspects of fitting that are not covered by the above arguments. See Details.
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
  \item{ps}{The estimated or provided propensity scores. Estimated propensity scores are returned for binary treatments and only when \code{method} is \code{"glm"}, \code{"gbm"}, \code{"cbps"}, \code{"super"}, or \code{"bart"}.}
  \item{s.weights}{The provided sampling weights.}
  \item{focal}{The focal variable if the ATT was requested with a multi-category treatment.}
  \item{by}{A data.frame containing the \code{by} variable when specified.}
  \item{obj}{When \code{include.obj = TRUE}, the fit object.}
  \item{info}{Additional information about the fitting. See the individual methods pages for what is included.}
}
\details{
The primary purpose of \code{weightit()} is as a dispatcher to functions that perform the estimation of balancing weights using the requested \code{method}. Below are the methods allowed and links to pages containing more information about them, including additional arguments and outputs (e.g., when \code{include.obj = TRUE}), how missing values are treated, which estimands are allowed, and whether sampling weights are allowed.
\itemize{
\item{\code{\link[=method_glm]{"glm"}} - Propensity score weighting using generalized linear models.}
\item{\code{\link[=method_gbm]{"gbm"}} - Propensity score weighting using generalized boosted modeling.}
\item{\code{\link[=method_cbps]{"cbps"}} - Covariate Balancing Propensity Score weighting.}
\item{\code{\link[=method_npcbps]{"npcbps"}} - Non-parametric Covariate Balancing Propensity Score weighting.}
\item{\code{\link[=method_ebal]{"ebal"}} - Entropy balancing.}
%\item{\code{\link[=method_ebcw]{"ebcw"}} - Empirical balancing calibration weighting.}
\item{\code{\link[=method_optweight]{"optweight"}} - Optimization-based weighting.}
\item{\code{\link[=method_super]{"super"}} - Propensity score weighting using SuperLearner.}
\item{\code{\link[=method_bart]{"bart"}} - Propensity score weighting using Bayesian additive regression trees (BART).}
\item{\code{\link[=method_energy]{"energy"}} - Energy balancing.}
}

\code{method} can also be supplied as a user-defined function; see \code{\link{method_user}} for instructions and examples.

When using \code{weightit()}, please cite both the \pkg{WeightIt} package (using \code{citation("WeightIt")}) and the paper(s) in the references section of the method used.

}
\seealso{
\fun{weightitMSM} for estimating weights with sequential (i.e., longitudinal) treatments for use in estimating marginal structural models (MSMs).

\fun{weightit.fit}, which is a lower-level dispatcher function that accepts a matrix of covariates and a vector of treatment statuses rather than a formula and data frame and performs minimal argument checking and processing. It may be useful for speeding up simulation studies for which the correct arguments are known. In general \code{weightit()} should be used.
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
                method = "glm", estimand = "ATT"))
summary(W1)
bal.tab(W1)
\dontshow{if (requireNamespace("osqp", quietly = TRUE)) force(\{ # examplesIf}
#Balancing covariates with respect to race (multi-category)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "ebal", estimand = "ATE"))
summary(W2)
bal.tab(W2)
\dontshow{\}) # examplesIf}
\dontshow{if (requireNamespace("CBPS", quietly = TRUE)) force(\{ # examplesIf}
#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "cbps", over = FALSE))
summary(W3)
bal.tab(W3)
\dontshow{\}) # examplesIf}}
