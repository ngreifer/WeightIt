\name{sbps}
\alias{sbps}
\alias{SBPS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subgroup Balancing Propensity Score
}
\description{
Implements the subgroup balancing propensity score (SBPS), which is an algorithm that attempts to achieve balance in subgroups by sharing information from the overall sample (Dong, Zhang, Zeng, & Li, 2019; DZZL). Each subgroup can use either weights estimated using the whole sample, weights estimated using just that subgroup, or a combination of the two. The optimal combination is chosen as that which minimizes an imbalance criterion that includes subgroup as well as overall balance.
}
\usage{
sbps(obj, obj2 = NULL,
    moderator = NULL,
    formula = NULL,
    data = NULL,
    smooth = FALSE,
    full.search)

\method{print}{weightit.sbps}(x, ...)

\method{summary}{weightit.sbps}(object, top = 5,
        ignore.s.weights = FALSE, ...)

\method{print}{summary.weightit.sbps}(x, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
a \code{weightit} object used to estimate weights in the overall sample.
}
  \item{obj2}{
a \code{weightit} object used to estimate weights in the subgroups. Typically this has been estimated by include \code{by} in the call to \code{weightit()}. Either \code{obj2} or \code{moderator} must be specified.
}
  \item{moderator}{
optional; a vector or the names of variables for which weighting is to be done within subgroups. For example, if \code{by = "gender"}, weights will be generated separately within each level of the variable \code{"gender"}. This argument is analogous to the \code{by} argument in \code{weightit()}, and in fact it is passed on to \code{by}. Either \code{obj2} or \code{moderator} must be specified.
}
  \item{formula}{
an optional formula with the covariates whose balance is to be optimized. If not specified, the formula in \code{obj} will be used.
}
  \item{data}{
an optional data set in the form of a data frame that contains the variables in \code{formula} or \code{moderator}.
}
  \item{smooth}{
\code{logical}; whether the smooth version of the SBPS should be used. This is only compatible with \code{weightit} methods that estimate a propensity score.
}
  \item{full.search}{
\code{logical}; when \code{smooth = FALSE}, whether every combination of subgroup and overall weights should be evaluated. If \code{FALSE}, a stochastic search as described in DZZL will be used instead. If \code{TRUE}, all 2^R combinations will be checked, where R is the number of subgroups. If unspecified, will default to \code{TRUE} if R <= 8 and \code{FALSE} otherwise.
}
  \item{x}{
a \code{weightit.sbps}, \code{summary.weightit} or \code{summary.weightitMSM} object; the output of a call to \code{sbps}, \code{summary.weightit} or \code{summary.weightitMSM}.
}
  \item{object}{
a \code{weightit.sbps} object; the output of a call to \code{sbps}.
}
  \item{top}{
how many of the largest and smallest weights to display. Default is 5.
}
  \item{ignore.s.weights}{
whether or not to ignore sampling weights when computing the weight summary. If \code{FALSE}, the default, the estimated weights will be multiplied by the sampling weights (if any) before values are computed.
}

  \item{...}{
for \code{print}, arguments passed to \code{\link{print}}.
}
}
\details{
The SBPS relies on two sets of weights: one estimated in the overall sample and one estimated within each subgroup. The algorithm decides whether each subgroup should use the weights estimated in the overall sample or those estimated in the subgroup. Where R is the number of subgroups, there are 2^R permutations of overall and subgroup weights. The optimal permutation is chosen as that which minimizes a balcne criterion as described in DZZL. The balance criterion used here is the sum of the squared standardized mean differences within subgroups and overall, which are computed using \code{\link[cobalt]{col_w_smd}} in \pkg{cobalt}.

The smooth version estimates weights that determine the relative contribution of the overall and subgroup propensity scores to a weighted average propensity score for each subgroup. If P_O are the propensity scores estimated in the overall sample and P_S are the propensity scores estimated in each subgroups, the smooth SBPS finds coefficients C so that for each subgroup, the ultimate propensity score is \eqn{C*P_S + (1-C)*P_O}, and weights are computed from this propensity score. When C is estimated to be 1 or 0, the smooth SBPS coincides with the standard SBPS.

If \code{obj2} is not specified and \code{moderator} is, \code{sbps()} will attempt to refit the model specified in \code{obj} with the \code{moderator} in the \code{by} argument. This relies on the environment in which \code{obj} was created to be intact and can take some time if \code{obj} was hard to fit. It's safer to estimate \code{obj} and \code{obj} (the latter simply by including the moderator in the \code{by} argument) and supply these to \code{SBPS}.

}
\value{
A \code{weightit.sbps} obejct, which inherits from \code{weightit}. This contains all the information in \code{obj} with the weights, propensity scores, call, and possibly covariates updated from \code{sbps()}. In addition, the \code{prop.subgroup} contains the values of the coefficients C for the subgroups (which are either 0 or 1 for the standard SBPS).

This object has its own summary methods and is compatible with \pkg{cobalt} functions. The \code{cluster} argument should be used with \pkg{cobalt} functions to accurately reflect the performance of the weights in balancing the subgroups.
}
\references{
Dong, J., Zhang, J. L., Zeng, S., & Li, F. (2019). Subgroup Balancing Propensity Score.
}
\author{
Noah Greifer
}

\seealso{
\code{\link{weightit}}, \code{link{summary.weightit}}
}
\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups
(W1 <- weightit(treat ~ age + educ + married +
                nodegree + race + re74, data = lalonde,
                method = "ps", estimand = "ATT"))

(W2 <- weightit(treat ~ age + educ + married +
                nodegree + race + re74, data = lalonde,
                method = "ps", estimand = "ATT",
                by = "race"))
S <- sbps(W1, W2)
print(S)
summary(S)
bal.tab(S, cluster = "race")

S_ <- sbps(W1, W2, smooth = TRUE)
print(S_)
summary(S_)
bal.tab(S_, cluster = "race")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
