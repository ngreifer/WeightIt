% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sbps.R
\name{sbps}
\alias{sbps}
\title{Subgroup Balancing Propensity Score}
\usage{
sbps(
  obj,
  obj2 = NULL,
  moderator = NULL,
  formula = NULL,
  data = NULL,
  smooth = FALSE,
  full.search
)
}
\arguments{
\item{obj}{a \code{weightit} object containing weights estimated in the overall
sample.}

\item{obj2}{a \code{weightit} object containing weights estimated in the
subgroups. Typically this has been estimated by including \code{by} in the call
to \code{\link[=weightit]{weightit()}}. Either \code{obj2} or \code{moderator} must be specified.}

\item{moderator}{optional; a string containing the name of the variable in
\code{data} for which weighting is to be done within subgroups or a one-sided
formula with the subgrouping variable on the right-hand side. This argument
is analogous to the \code{by} argument in \code{weightit()}, and in fact it is passed
on to \code{by}. Either \code{obj2} or \code{moderator} must be specified.}

\item{formula}{an optional formula with the covariates for which balance is
to be optimized. If not specified, the formula in \code{obj$call} will be used.}

\item{data}{an optional data set in the form of a data frame that contains
the variables in \code{formula} or \code{moderator}.}

\item{smooth}{\code{logical}; whether the smooth version of the SBPS should be
used. This is only compatible with \code{weightit} methods that return a
propensity score.}

\item{full.search}{\code{logical}; when \code{smooth = FALSE}, whether every
combination of subgroup and overall weights should be evaluated. If
\code{FALSE}, a stochastic search as described in DZZL will be used instead. If
\code{TRUE}, all \eqn{2^R} combinations will be checked, where \eqn{R} is the
number of subgroups, which can take a long time with many subgroups. If
unspecified, will default to \code{TRUE} if \eqn{R <= 8} and \code{FALSE} otherwise.}
}
\value{
A \code{weightit.sbps} object, which inherits from \code{weightit}. This
contains all the information in \code{obj} with the weights, propensity scores,
call, and possibly covariates updated from \code{sbps()}. In addition, the
\code{prop.subgroup} component contains the values of the coefficients C for the
subgroups (which are either 0 or 1 for the standard SBPS), and the
\code{moderator} component contains a data.frame with the moderator.

This object has its own summary method and is compatible with \pkg{cobalt}
functions. The \code{cluster} argument should be used with \pkg{cobalt} functions
to accurately reflect the performance of the weights in balancing the
subgroups.
}
\description{
Implements the subgroup balancing propensity score (SBPS), which
is an algorithm that attempts to achieve balance in subgroups by sharing
information from the overall sample and subgroups (Dong, Zhang, Zeng, & Li,
2020; DZZL). Each subgroup can use either weights estimated using the whole
sample, weights estimated using just that subgroup, or a combination of the
two. The optimal combination is chosen as that which minimizes an imbalance
criterion that includes subgroup as well as overall balance.
}
\details{
The SBPS relies on two sets of weights: one estimated in the overall
sample and one estimated within each subgroup. The algorithm decides whether
each subgroup should use the weights estimated in the overall sample or those
estimated in the subgroup. There are 2^R permutations of overall and subgroup
weights, where R is the number of subgroups. The optimal permutation is
chosen as that which minimizes a balance criterion as described in DZZL. The
balance criterion used here is, for binary and multi-category treatments, the
sum of the squared standardized mean differences within subgroups and
overall, which are computed using \code{\link[cobalt:balance-summary]{cobalt::col_w_smd()}}, and for continuous
treatments, the sum of the squared correlations between each covariate and
treatment within subgroups and overall, which are computed using
\code{\link[cobalt:balance-summary]{cobalt::col_w_corr()}}.

The smooth version estimates weights that determine the relative contribution
of the overall and subgroup propensity scores to a weighted average
propensity score for each subgroup. If P_O are the propensity scores
estimated in the overall sample and P_S are the propensity scores estimated
in each subgroup, the smooth SBPS finds R coefficients C so that for each
subgroup, the ultimate propensity score is \eqn{C*P_S + (1-C)*P_O}, and
weights are computed from this propensity score. The coefficients are
estimated using \code{\link[=optim]{optim()}} with \code{method = "L-BFGS-B"}. When C is estimated to
be 1 or 0 for each subgroup, the smooth SBPS coincides with the standard
SBPS.

If \code{obj2} is not specified and \code{moderator} is, \code{sbps()} will attempt to refit
the model specified in \code{obj} with the \code{moderator} in the \code{by} argument. This
relies on the environment in which \code{obj} was created to be intact and can
take some time if \code{obj} was hard to fit. It's safer to estimate \code{obj} and
\code{obj2} (the latter simply by including the moderator in the \code{by} argument)
and supply these to \code{sbps()}.
}
\examples{

library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups within races
(W1 <- weightit(treat ~ age + educ + married +
                nodegree + race + re74, data = lalonde,
                method = "glm", estimand = "ATT"))

(W2 <- weightit(treat ~ age + educ + married +
                nodegree + race + re74, data = lalonde,
                method = "glm", estimand = "ATT",
                by = "race"))
S <- sbps(W1, W2)
print(S)
summary(S)
bal.tab(S, cluster = "race")

#Could also have run
#  sbps(W1, moderator = "race")

S_ <- sbps(W1, W2, smooth = TRUE)
print(S_)
summary(S_)
bal.tab(S_, cluster = "race")
}
\references{
Dong, J., Zhang, J. L., Zeng, S., & Li, F. (2020). Subgroup
balancing propensity score. Statistical Methods in Medical Research, 29(3),
659–676. \doi{10.1177/0962280219870836}
}
\seealso{
\code{\link[=weightit]{weightit()}}, \code{\link[=summary.weightit]{summary.weightit()}}
}
