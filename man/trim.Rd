\name{trim}
\alias{trim}
\alias{trim.weightit}
\alias{trim.numeric}

\title{
Trim Large Weights
}
\description{
Trims (i.e., truncates) large weights by setting all weights higher than that at a given quantile to the weight at the quantile. This can be useful in controlling extreme weights, which can reduce effective sample size by enlarging the variability of the weights.
}
\usage{
\method{trim}{weightit}(w, at = 0, lower = FALSE, ...)

\method{trim}{numeric}(w, at = 0, lower = FALSE, treat = NULL, ...)
}
\arguments{
  \item{w}{
A \code{weightit} object or a vector of weights.
}
  \item{at}{
\code{numeric}; either the quantile of the weights above which weights are to be trimmed. A single number between .5 and 1, or the number of weights to be trimmed (e.g., \code{at = 3} for the top 3 weights to be set to the 4th largest weight).
}
  \item{lower}{
\code{logical}; whether also to trim at the lower quantile (e.g., for \code{at = .9}, trimming at both .1 and .9, or for \code{at = 3}, trimming the top and bottom 3 weights).
}
  \item{treat}{
A vector of treatment status for each unit. This should always be included when \code{w} is numeric, but you can get away with leaving it out if the treatment is continuous or the estimand is the ATE for binary or multinomial treatments.
}
  \item{\dots}{
Not used.
}
}
\details{
\code{trim} takes in a \code{weightit} object (the output of a call to \code{weightit} or \code{weightitMSM}) or a numeric vector of weights and trims them to the specified quantile. All weights above that quantile are set to the weight at that quantile. If \code{lower = TRUE}, all weights below 1 minus the quantile are to set the weight at 1 minus the quantile. In general, trimming weights decreases balance but also decreases the variability of the weights, improving precision at the potential expense of unbiasedness (Cole & Hernán, 2008). See Lee, Lessler, and Stuart (2011) and Thoemmes and Ong (2015) for discussions and simulation results of trimming weights at various quantiles. Note that trimming weights can also change the target population and therefore the estimand.

When using \code{trim} on a numeric vector of weights, it is helpful to include the treatment vector as well. The helps determine the type of treatment and estimand, which are used to specify how trimming is performed. In particular, if the estimand is determined to be the ATT or ATC, the weights of the target (i.e., focal) group are ignored, since they should all be equal to 1. Otherwise, if the estimand is the ATE or the treatment is continuous, all weights are considered for trimming. In general, weights for any group for which all the weights are the same will not be considered in the trimming.
}
\value{
If the input is a \code{weightit} object, the output will be a \code{weightit} object with the weights replaced by the trimmed weights, while will have an additional attribute, \code{"trim"}, equal to the quantile of trimming.

If the input is a numeric vector of weights, the output will be a numeric vector of the trimmed weights, again with the aforementioned attribute.
}
\references{
Cole, S. R., & Hernán, M. Á. (2008). Constructing Inverse Probability Weights for Marginal Structural Models. American Journal of Epidemiology, 168(6), 656–664.

Lee, B. K., Lessler, J., & Stuart, E. A. (2011). Weight Trimming and Propensity Score Weighting. PLoS ONE, 6(3), e18174.

Thoemmes, F., & Ong, A. D. (2016). A Primer on Inverse Probability of Treatment Weighting and Marginal Structural Models. Emerging Adulthood, 4(1), 40–59.

}
\author{
Noah Greifer
}

\seealso{
\code{link{weightit}}, \code{link{weightitMSM}}
}
\examples{
library("cobalt")
data("lalonde", package = "cobalt")

(W <- weightit(treat ~ age + educ + married +
                 nodegree + re74, data = lalonde,
               method = "ps", estimand = "ATT"))
summary(W)

#Trimming the top and bottom 5 weights
trim(W, at = 5, lower = TRUE)

#Trimming at 90th percentile
(W.trim <- trim(W, at = .9))

summary(W.trim)
#Note that only the control weights were trimmed

#Trimming a numeric vector of weights
all.equal(trim(W$weights, at = .9, treat = lalonde$treat),
          W.trim$weights)

#Using made up data and as.weightit()
treat <- rbinom(500, 1, .3)
weights <- rchisq(500, df = 2)
W <- as.weightit(weights = weights, treat = treat,
                 estimand = "ATE")
summary(W)
summary(trim(W, at = .95))
}
