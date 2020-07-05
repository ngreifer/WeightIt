\name{summary.weightit}
\alias{print.summary.weightit}
\alias{plot.summary.weightit}
\alias{summary.weightit}
\alias{print.summary.weightitMSM}
\alias{summary.weightitMSM}
\title{Print and Summarize Output}
\usage{

\method{summary}{weightit}(object, top = 5,
        ignore.s.weights = FALSE, ...)

\method{print}{summary.weightit}(x, ...)

\method{plot}{summary.weightit}(x, ...)

\method{summary}{weightitMSM}(object, top = 5,
        ignore.s.weights = FALSE, ...)

\method{print}{summary.weightitMSM}(x, ...)
}
\description{
\code{summary} generates a summary of the \code{weightit} or \code{weightitMSM} object to evaluate the properties of the estimated weights. \code{plot} plots the distribution of the weights.
}
\arguments{
  \item{object}{
a \code{weightit} or \code{weightitMSM} object; the output of a call to \code{weightit} or \code{weightitMSM}.
}
  \item{top}{
how many of the largest and smallest weights to display. Default is 5.
}
  \item{ignore.s.weights}{
whether or not to ignore sampling weights when computing the weight summary. If \code{FALSE}, the default, the estimated weights will be multiplied by the sampling weights (if any) before values are computed.
}
  \item{x}{
a \code{summary.weightit} or \code{summary.weightitMSM} object; the output of a call to \code{summary.weightit} or \code{summary.weightitMSM}.
}
  \item{...}{
for \code{print}, arguments passed to \code{\link{print}}. For \code{plot}, additional arguments passed to \code{\link[graphics]{hist}} to determine the number of bins, though \code{\link[ggplot2]{geom_histogram}} from \pkg{ggplot2} is actually used to create the plot.
  }
}

\value{
For point treatments (i.e., \code{weightit} objects), a \code{summary.weightit} object with the following elements:
\item{weight.range}{The range (minimum and maximum) weight for each treatment group.}
\item{weight.top}{The units with the greatest weights in each treatment group; how many are included is determined by \code{top}.}
\item{coef.of.var (Coef of Var)}{The coefficient of variation (standard deviation divided by mean) of the weights in each treatment group and overall.}
\item{scaled.mad (MAD)}{The mean absolute deviation of the weights in each treatment group and overall divided by the mean of the weights in the corresponding group.}
\item{negative entropy (Entropy)}{The negative entropy (\eqn{\sum w log(w)}) of the weights in each treatment group and overall divided by the mean of the weights in the corresponding group.}
\item{num.zeros}{The number of weights equal to zero.}
\item{effective.sample.size}{The effective sample size for each treatment group before and after weighting. See \code{\link{ESS}}.}

For longitudinal treatments (i.e., \code{weightitMSM} objects), a list of the above elements for each treatment period.

\code{plot} returns a \code{ggplot} object with a histogram displaying the distribution of the estimated weights. If the estimand is the ATT or ATC, only the weights for the non-focal group(s) will be displayed (since the weights for the focal group are all 1). A dotted line is displayed at the mean of the weights.
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
