\name{get_w_from_ps}
\alias{get_w_from_ps}
\title{
Get weights from propensity scores corresponding to different estimands
}
\description{
Given a vector or matrix of propensity scores, output a vector of weights that target the provided estimand.
}
\usage{
get_w_from_ps(ps,
              treat,
              estimand = "ATE",
              focal = NULL,
              treated = NULL)
}
\arguments{
  \item{ps}{
A vector, matrix, or data frame of propensity scores. See Details.
}
  \item{treat}{
A vector of treatment statuses for each individual. See Details.
}
  \item{estimand}{
The desired estimand that the weights should target. Current options include "ATE" (average treatment effect), "ATT" (average treatment effect on the treated), "ATC" (average treatment effect on the control), "ATO" (average treatment effect in the overlap), and "ATM" (average treatment effect in the matched sample). See \code{\link{method_ps}} for additional information and references.
}
  \item{focal}{
When the estimand is the ATT or ATC, which group should be consider the (focal) "treated" or "control" group, respectively. If not \code{NULL} and \code{estimand} is not "ATT" or "ATC", \code{estimand} will automatically be set to "ATT".
}
  \item{treated}{
When treatment is binary, the value of \code{treat} that is considered the "treated" group (i.e., the group for which the propensity scores are the probability of being in). If \code{NULL}, \code{get_w_from_ps} will attempt to figure it out on its own using some heuristics. This really only matters when \code{treat} has values other than 0 and 1 and when \code{ps} is given as a vector or an unnamed single-column matrix or data frame.
}
}
\details{
\code{get_w_from_ps} applies the formula for computing weights from propensity scores for the desired estimand. See the References section at \code{\link{method_ps}} for information on these estimands and the formulas.

\code{ps} can be entered a variety of ways. For binary treatments, when \code{ps} is entered as a vector or unnamed single-column matrix or data frame, \code{get_w_from_ps} has to know which value of \code{treat} corresponds to the "treated" group. For 0/1 variables, 1 will be considered treated. For other types of variables, \code{get_w_from_ps} will try to figure it out using heuristics, but it's safer to supply an argument to \code{treated}. When \code{estimand} is "ATT" or "ATC", supplying a value to \code{focal} is sufficient (for ATT, \code{focal} is the treated group, and for ATC, \code{focal} is the control group). When entered as a matrix or data frame, the columns must be named with the levels of the treatment, and it is assumed that each column corresponds to the probability of being in that treatment group. This is the safest way to supply \code{ps} unless \code{treat} is a 0/1 variable.

For multinomial treatments, \code{ps} can be entered as a vector or a matrix or data frame. When entered as a vector, it is assumed the value corresponds to the probability of being in the treatment actually received; this is only possible when the estimand is "ATE". Otherwise, \code{ps} must be entered as a named matrix or data frame as described above for binary treatments. When the estimand is "ATT" or "ATC", a value for \code{focal} must be specified.

\code{get_w_from_ps} is not compatible with continuous treatments.
}
\value{
A vector of weights.
}

\author{
Noah Greifer
}

\seealso{
\code{\link{method_ps}}
}
\examples{
library("cobalt")
data("lalonde", package = "cobalt")

ps.fit <- glm(treat ~ age + educ + race + married +
                nodegree + re74 + re75, data = lalonde,
              family = binomial)
ps <- ps.fit$fitted.values

w1 <- get_w_from_ps(ps, treat = lalonde$treat,
                    estimand = "ATT")

treatAB <- factor(ifelse(lalonde$treat == 1, "A", "B"))
w2 <- get_w_from_ps(ps, treat = treatAB,
                    estimand = "ATT", focal = "A")
all.equal(w1, w2)
w3 <- get_w_from_ps(ps, treat = treatAB,
                    estimand = "ATT", treated = "A")
all.equal(w1, w3)

#A multinomial example using GBM predicted probabilities
library(gbm)
T3 <- factor(sample(c("A", "B", "C"), nrow(lalonde), replace = TRUE))

gbm.fit <- gbm(T3 ~ age + educ + race + married +
                 nodegree + re74 + re75, data = lalonde,
               distribution = "multinomial", n.trees = 200,
               interaction.depth = 3)
ps.multi <- drop(predict(gbm.fit, type = "response",
                         n.trees = 200))
w <- get_w_from_ps(ps.multi, T3, estimand = "ATE")

}