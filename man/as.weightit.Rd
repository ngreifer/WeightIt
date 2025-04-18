% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.weightit.R
\name{as.weightit}
\alias{as.weightit}
\alias{as.weightit.weightit.fit}
\alias{as.weightit.default}
\alias{as.weightitMSM}
\alias{as.weightitMSM.default}
\title{Create a \code{weightit} object manually}
\usage{
as.weightit(x, ...)

\method{as.weightit}{weightit.fit}(x, covs = NULL, ...)

\method{as.weightit}{default}(
  x,
  treat,
  covs = NULL,
  estimand = NULL,
  s.weights = NULL,
  ps = NULL,
  ...
)

as.weightitMSM(x, ...)

\method{as.weightitMSM}{default}(
  x,
  treat.list,
  covs.list = NULL,
  estimand = NULL,
  s.weights = NULL,
  ps.list = NULL,
  ...
)
}
\arguments{
\item{x}{required; a \code{numeric} vector of weights, one for each unit, or a
\code{weightit.fit} object from \code{\link[=weightit.fit]{weightit.fit()}}.}

\item{...}{additional arguments. These must be named. They will be included
in the output object.}

\item{covs}{an optional \code{data.frame} of covariates. For using \pkg{WeightIt}
functions, this is not necessary, but for use with \pkg{cobalt} it is. Note
that when using with a \code{weightit.fit} object, this should not be the matrix
supplied to the \code{covs} argument of \code{weightit.fit()} unless there are no
factor/character variables in it. Ideally this is the original, unprocessed
covariate data frame with factor variables included.}

\item{treat}{a vector of treatment statuses, one for each unit. Required when
\code{x} is a vector of weights.}

\item{estimand}{an optional \code{character} of length 1 giving the estimand. The
text is not checked.}

\item{s.weights}{an optional \code{numeric} vector of sampling weights, one for
each unit.}

\item{ps}{an optional \code{numeric} vector of propensity scores, one for each
unit.}

\item{treat.list}{a list of treatment statuses at each time point.}

\item{covs.list}{an optional list of \code{data.frame}s of covariates of
covariates at each time point. For using \pkg{WeightIt} functions, this is
not necessary, but for use with \pkg{cobalt} it is.}

\item{ps.list}{an optional list of \code{numeric} vectors of propensity scores at
each time point.}
}
\value{
An object of class \code{weightit} (for \code{as.weightit()}) or \code{weightitMSM}
(for \code{as.weightitMSM()}).
}
\description{
This function allows users to get the benefits of a \code{weightit}
object when using weights not estimated with \code{\link[=weightit]{weightit()}} or \code{\link[=weightitMSM]{weightitMSM()}}.
These benefits include diagnostics, plots, and direct compatibility with
\pkg{cobalt} for assessing balance.
}
\examples{

treat <- rbinom(500, 1, .3)
weights <- rchisq(500, df = 2)

W <- as.weightit(weights, treat = treat, estimand = "ATE")
summary(W)

# See ?weightit.fit for using as.weightit() with a
# weightit.fit object.

}
