\name{as.weightit}
\alias{as.weightit}
\alias{as.weightitMSM}
\alias{as.weightit.default}
\alias{as.weightitMSM.default}

\title{
Create a \code{weightit} object manually
}
\description{
This function allows users to get the benefits of a \code{weightit} object when using weights not estimated with \code{weightit()} or \code{weightitMSM()}. These benefits include diagnostics, plots, and direct compatibility with \pkg{cobalt} for assessing balance.
}
\usage{
as.weightit(...)

\method{as.weightit}{default}(weights,
            treat,
            covs = NULL,
            estimand = NULL,
            s.weights = NULL,
            ps = NULL,
            ...)

as.weightitMSM(...)

\method{as.weightitMSM}{default}(weights,
               treat.list,
               covs.list = NULL,
               estimand = NULL,
               s.weights = NULL,
               ps.list = NULL,
               ...)
}

\arguments{
  \item{weights}{
required; a \code{numeric} vector of weights, one for each unit.
}
  \item{treat}{
required; a vector of treatment statuses, one for each unit.
}
  \item{covs}{
an optional \code{data.frame} of covariates. For using \pkg{WeightIt} functions, this is not necessary, but for use with \pkg{cobalt} it is.
}
  \item{estimand}{
an optional \code{character} of length 1 giving the estimand. The text is not checked.
}
  \item{s.weights}{
an optional \code{numeric} vector of sampling weights, one for each unit.
}
  \item{ps}{
an optional \code{numeric} vector of propensity scores, one for each unit.
}
  \item{treat.list}{
a list of treatment statuses at each time point..
}
  \item{covs.list}{
an optional list of \code{data.frame}s of covariates of covariates at each time point. For using \pkg{WeightIt} functions, this is not necessary, but for use with \pkg{cobalt} it is.
}
  \item{ps.list}{
an optional list of \code{numeric} vectors of propensity scores at each time point.
}
  \item{...}{
additional arguments. These must be named. They will be included in the output object.
}
}

\value{
An object of class \code{weightit} (for \code{as.weightit}) or \code{weightitMSM} (for \code{as.weightitMSM}).
}

\author{
Noah Greifer
}
\examples{
treat <- rbinom(500, 1, .3)
weights <- rchisq(500, df = 2)
W <- as.weightit(weights= weights, treat = treat,
                 estimand = "ATE")
summary(W)
}
