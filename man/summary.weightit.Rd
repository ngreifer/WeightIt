\name{summary.weightit}
\alias{print.summary.weightit}
\alias{summary.weightit}
\title{Print and Summarize Output}
\usage{

\method{summary}{weightit}(object, top = 5,
        ignore.s.weights = FALSE, ...)

\method{print}{summary.weightit}(x, ...)
}
\description{
\code{summary.weightit()} generates a summary of the \code{weightit} object to evaluate the properties of the estimated weights.
}
\arguments{
  \item{object}{
a \code{weightit} object; the output of a call to \code{weightit()}.
}
  \item{top}{
how many of the largest and smallest weights to display. Default is 5.
}
  \item{ignore.s.weights}{
whether or not to ignore sampling weights when computing the weight summary. If \code{FALSE}, the default, the estimated weights will be multiplied by the sampling weights (if any) before values are computed.
}
  \item{x}{
a \code{summary.weightit} object; the output of a call to \code{summary.weightit()}.
}
  \item{...}{
  aruments passed to \code{\link{print}}.
  }
}

\value{
A \code{summary.weightit} object with the following elements:
\item{weight.range}{The range (minimum and maximum) weight for each treatment group.}
\item{weight.top}{The units with the greatest weights in each treatment group; how many are included is determined by \code{top}.}
\item{weight.ratio}{The ratio of the largest weight to the smallest weight in each treatment group and overall.}
\item{coef.of.var}{The coefficient of variation (standard deviation divided by mean) of the weights in each treatment group and overall.}
\item{effective.sample.size}{The effective sample size for each treatment group before and after weighting.}
}

\author{
Noah Greifer \email{noah@unc.edu}
}

\seealso{
\code{\link{weightit}}, \code{\link[base]{summary}}
}
\references{
Li, F., Morgan, K. L., & Zaslavsky, A. M. (2016). Balancing Covariates via Propensity Score Weighting. Journal of the American Statistical Association, 0(ja), 0–0.

}
\examples{
# See example at ?weightit.
}
