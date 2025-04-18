% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weightit.fit.R
\name{weightit.fit}
\alias{weightit.fit}
\title{Generate Balancing Weights with Minimal Input Processing}
\usage{
weightit.fit(
  covs,
  treat,
  method = "glm",
  s.weights = NULL,
  by.factor = NULL,
  estimand = "ATE",
  focal = NULL,
  stabilize = FALSE,
  ps = NULL,
  moments = NULL,
  int = FALSE,
  subclass = NULL,
  missing = NULL,
  verbose = FALSE,
  include.obj = FALSE,
  ...
)
}
\arguments{
\item{covs}{a numeric matrix of covariates.}

\item{treat}{a vector of treatment statuses.}

\item{method}{a string containing the name of the method that will be used to
estimate weights. See \code{\link[=weightit]{weightit()}} for allowable options. The default is
\code{"glm"} for propensity score weighting using a generalized linear model to
estimate the propensity score.}

\item{s.weights}{a numeric vector of sampling weights. See the individual
pages for each method for information on whether sampling weights can be
supplied.}

\item{by.factor}{a factor variable for which weighting is to be done within
levels. Corresponds to the \code{by} argument in \code{\link[=weightit]{weightit()}}.}

\item{estimand}{the desired estimand. For binary and multi-category
treatments, can be \code{"ATE"}, \code{"ATT"}, \code{"ATC"}, and, for some methods,
\code{"ATO"}, \code{"ATM"}, or \code{"ATOS"}. The default for both is \code{"ATE"}. This
argument is ignored for continuous treatments. See the individual pages for
each method for more information on which estimands are allowed with each
method and what literature to read to interpret these estimands.}

\item{focal}{when \code{estimand} is set to \code{"ATT"} or \code{"ATC"}, which group to
consider the "treated" or "control" group. This group will not be weighted,
and the other groups will be weighted to resemble the focal group. If
specified, \code{estimand} will automatically be set to \code{"ATT"} (with a warning
if \code{estimand} is not \code{"ATT"} or \code{"ATC"}). See section \emph{\code{estimand} and
\code{focal}} in Details at \code{\link[=weightit]{weightit()}}.}

\item{stabilize}{\code{logical}; whether or not to stabilize the weights. For the
methods that involve estimating propensity scores, this involves
multiplying each unit's weight by the proportion of units in their
treatment group. Default is \code{FALSE}. Note this differs from its use with
\code{\link[=weightit]{weightit()}}.}

\item{ps}{a vector of propensity scores. If specified, \code{method} will be
ignored and set to \code{"glm"}.}

\item{moments, int, subclass}{arguments to customize the weight estimation. See
\code{\link[=weightit]{weightit()}} for details.}

\item{missing}{\code{character}; how missing data should be handled. The options
depend on the \code{method} used. If \code{NULL}, \code{covs} will be checked for \code{NA}
values, and if present, \code{missing} will be set to \code{"ind"}. If \code{""}, \code{covs}
will not be checked for \code{NA} values; this can be faster when it is known
there are none.}

\item{verbose}{\code{logical}; whether to print additional information output by
the fitting function.}

\item{include.obj}{\code{logical}; whether to include in the output any fit
objects created in the process of estimating the weights. For example, with
\code{method = "glm"}, the \code{glm} objects containing the propensity score model
will be included. See the individual pages for each method for information
on what object will be included if \code{TRUE}.}

\item{...}{other arguments for functions called by \code{weightit.fit()} that
control aspects of fitting that are not covered by the above arguments.}
}
\value{
A \code{weightit.fit} object with the following elements:
\item{weights}{The estimated weights, one for each unit.} \item{treat}{The
values of the treatment variable.} \item{estimand}{The estimand requested.}
\item{method}{The weight estimation method specified.}
\item{ps}{The estimated or provided propensity scores. Estimated propensity scores are
returned for binary treatments and only when \code{method} is \code{"glm"}, \code{"gbm"}, \code{"cbps"}, \code{"ipt"}, \code{"super"}, or \code{"bart"}. The propensity score corresponds to the predicted probability of being treated; see section \emph{\code{estimand} and \code{focal}} in Details at \code{\link[=weightit]{weightit()}} for how the treated group is determined.}
\item{s.weights}{The provided sampling weights.} \item{focal}{The focal
treatment level if the ATT or ATC was requested.} \item{fit.obj}{When
\code{include.obj = TRUE}, the fit object.} \item{info}{Additional information
about the fitting. See the individual methods pages for what is included.}

The \code{weightit.fit} object does not have specialized \code{print()}, \code{summary()},
or \code{plot()} methods. It is simply a list containing the above components. Use
\code{\link[=as.weightit]{as.weightit()}} to convert it to a \code{weightit} object, which does have these
methods. See Examples.
}
\description{
\code{weightit.fit()} dispatches one of the weight estimation methods
determined by \code{method}. It is an internal function called by \code{\link[=weightit]{weightit()}} and
should probably not be used except in special cases. Unlike \code{weightit()},
\code{weightit.fit()} does not accept a formula and data frame interface and
instead requires the covariates and treatment to be supplied as a numeric
matrix and atomic vector, respectively. In this way, \code{weightit.fit()} is to
\code{weightit()} what \code{\link[=lm.fit]{lm.fit()}} is to \code{\link[=lm]{lm()}} - a thinner, slightly faster
interface that performs minimal argument checking.
}
\details{
\code{weightit.fit()} is called by \code{\link[=weightit]{weightit()}} after the arguments to
\code{weightit()} have been checked and processed. \code{weightit.fit()} dispatches the
function used to actually estimate the weights, passing on the supplied
arguments directly. \code{weightit.fit()} is not meant to be used by anyone other
than experienced users who have a specific use case in mind. The returned
object contains limited information about the supplied arguments or details
of the estimation method; all that is processed by \code{weightit()}.

Less argument checking or processing occurs in \code{weightit.fit()} than does in
\code{weightit()}, which means supplying incorrect arguments can result in errors,
crashes, and invalid weights, and error and warning messages may not be
helpful in diagnosing the problem. \code{weightit.fit()} does check to make sure
weights were actually estimated, though.

\code{weightit.fit()} may be most useful in speeding up simulation simulation
studies that use \code{weightit()} because the covariates can be supplied as a
numeric matrix, which is often how they are generated in simulations, without
having to go through the potentially slow process of extracting the
covariates and treatment from a formula and data frame. If the user is
certain the arguments are valid (e.g., by ensuring the estimated weights are
consistent with those estimated from \code{weightit()} with the same arguments),
less time needs to be spent on processing the arguments. Also, the returned
object is much smaller than a \code{weightit} object because the covariates are
not returned alongside the weights.
}
\examples{

library("cobalt")
data("lalonde", package = "cobalt")

# Balancing covariates between treatment groups (binary)
covs <- lalonde[c("age", "educ", "race", "married",
                  "nodegree", "re74", "re75")]
## Create covs matrix, splitting any factors using
## cobalt::splitfactor()
covs_mat <- as.matrix(splitfactor(covs))

WF1 <- weightit.fit(covs_mat, treat = lalonde$treat,
                    method = "glm", estimand = "ATT")
str(WF1)

# Converting to a weightit object for use with
# summary() and bal.tab()
W1 <- as.weightit(WF1, covs = covs)
W1
summary(W1)
bal.tab(W1)

}
\seealso{
\code{\link[=weightit]{weightit()}}, which you should use for estimating weights unless you
know better.

\code{\link[=as.weightit]{as.weightit()}} for converting a \code{weightit.fit} object to a \code{weightit}
object.
}
