\name{summary.weightit}
\alias{print.summary.weightit}
\alias{summary.weightit}
\alias{print.summary.weightitMSM}
\alias{summary.weightitMSM}
\title{Print and Summarize Output}
\usage{

\method{summary}{weightit}(object, top = 5,
        ignore.s.weights = FALSE, ...)

\method{print}{summary.weightit}(x, ...)

\method{summary}{weightitMSM}(object, top = 5,
        ignore.s.weights = FALSE, ...)

\method{print}{summary.weightitMSM}(x, ...)
}
\description{
\code{summary()} generates a summary of the \code{weightit} or \code{weightitMSM} object to evaluate the properties of the estimated weights.
}
\arguments{
  \item{object}{
a \code{weightit} or \code{weightitMSM} object; the output of a call to \code{weightit()} or \code{weightitMSM()}.
}
  \item{top}{
how many of the largest and smallest weights to display. Default is 5.
}
  \item{ignore.s.weights}{
whether or not to ignore sampling weights when computing the weight summary. If \code{FALSE}, the default, the estimated weights will be multiplied by the sampling weights (if any) before values are computed.
}
  \item{x}{
a \code{summary.weightit} or \code{summary.weightitMSM} object; the output of a call to \code{summary.weightit()} or \code{summary.weightitMSM()}.
}
  \item{...}{
  arguments passed to \code{\link{print}}.
  }
}

\value{
For point treatments (i.e., \code{weightit} objects), a \code{summary.weightit} object with the following elements:
\item{weight.range}{The range (minimum and maximum) weight for each treatment group.}
\item{weight.top}{The units with the greatest weights in each treatment group; how many are included is determined by \code{top}.}
\item{weight.ratio}{The ratio of the largest weight to the smallest weight in each treatment group and overall.}
\item{coef.of.var}{The coefficient of variation (standard deviation divided by mean) of the weights in each treatment group and overall.}
\item{effective.sample.size}{The effective sample size for each treatment group before and after weighting.}

For longitudinal treatments (i.e., \code{weightitMSM} objects), a list of the above elements for each treatment period.
}

\author{
Noah Greifer
}

\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}, \code{\link[base]{summary}}
}

\examples{
# See example at ?weightit or ?weightitMSM
}
