\name{get_w_from_ps}
\alias{get_w_from_ps}
\title{
Compute weights from propensity scores
}
\description{
Given a vector or matrix of propensity scores, outputs a vector of weights that target the provided estimand.
}
\usage{
get_w_from_ps(ps,
              treat,
              estimand = "ATE",
              focal = NULL,
              treated = NULL,
              subclass = NULL,
              stabilize = FALSE)
}
\arguments{
  \item{ps}{
A vector, matrix, or data frame of propensity scores. See Details.
}
  \item{treat}{
A vector of treatment status for each individual. See Details.
}
  \item{estimand}{
The desired estimand that the weights should target. Current options include "ATE" (average treatment effect), "ATT" (average treatment effect on the treated), "ATC" (average treatment effect on the control), "ATO" (average treatment effect in the overlap), "ATM" (average treatment effect in the matched sample), and "ATOS" (average treatment effect in the optimal subset).
}
  \item{focal}{
When the estimand is the ATT or ATC, which group should be consider the (focal) "treated" or "control" group, respectively. If not \code{NULL} and \code{estimand} is not "ATT" or "ATC", \code{estimand} will automatically be set to "ATT".
}
  \item{treated}{
When treatment is binary, the value of \code{treat} that is considered the "treated" group (i.e., the group for which the propensity scores are the probability of being in). If \code{NULL}, \code{get_w_from_ps} will attempt to figure it out on its own using some heuristics. This really only matters when \code{treat} has values other than 0 and 1 and when \code{ps} is given as a vector or an unnamed single-column matrix or data frame.
}
  \item{subclass}{
\code{numeric}; the number of subclasses to use when computing weights using marginal mean weighting through stratification (also known as fine stratification). If \code{NULL}, standard inverse probability weights (and their extensions) will be computed; if a number greater than 1, subclasses will be formed and weights will be computed based on subclass membership. \code{estimand} must be ATE, ATT, or ATC if \code{subclass} is non-\code{NULL}. See Details.
}
  \item{stabilize}{
\code{logical}; whether to compute stabilized weights or not. This simply involves multiplying each unit's weight by the proportion of units in their treatment group. For saturated outcome models and in balance checking, this won't make a difference; otherwise, this can improve performance.
}
}
\details{
\code{get_w_from_ps} applies the formula for computing weights from propensity scores for the desired estimand. See the References section for information on these estimands and the formulas.

\code{ps} can be entered a variety of ways. For binary treatments, when \code{ps} is entered as a vector or unnamed single-column matrix or data frame, \code{get_w_from_ps} has to know which value of \code{treat} corresponds to the "treated" group. For 0/1 variables, 1 will be considered treated. For other types of variables, \code{get_w_from_ps} will try to figure it out using heuristics, but it's safer to supply an argument to \code{treated}. When \code{estimand} is "ATT" or "ATC", supplying a value to \code{focal} is sufficient (for ATT, \code{focal} is the treated group, and for ATC, \code{focal} is the control group). When entered as a matrix or data frame, the columns must be named with the levels of the treatment, and it is assumed that each column corresponds to the probability of being in that treatment group. This is the safest way to supply \code{ps} unless \code{treat} is a 0/1 variable.

For multi-category treatments, \code{ps} can be entered as a vector or a matrix or data frame. When entered as a vector, it is assumed the value corresponds to the probability of being in the treatment actually received; this is only possible when the estimand is "ATE". Otherwise, \code{ps} must be entered as a named matrix or data frame as described above for binary treatments. When the estimand is "ATT" or "ATC", a value for \code{focal} must be specified.

When \code{subclass} is not \code{NULL}, marginal mean weighting through stratification (MMWS) weights are computed. The implementation differs slightly from that described in Hong (2010, 2012). First, subclasses are formed by finding the quantiles of the propensity scores in the target group (for the ATE, all units; for the ATT or ATC, just the units in the focal group). Any subclasses lacking members of a treatment group will be filled in with them from neighboring subclasses so each subclass will always have at least one member of each treatment group. A new subclass-propensity score matrix is formed, where each unit's subclass-propensity score for each treatment value is computed as the proportion of units with that treatment value in the unit's subclass. For example, if a subclass had 10 treated units and 90 control units in it, the subclass-propensity score for being treated would be .1 and the subclass-propensity score for being control would be .9 for all units in the subclass. For multi-category treatments, the propensity scores for each treatment are stratified separately as described in Hong (2012); for binary treatments, only one set of propensity scores are stratified and the subclass-propensity scores for the other treatment are computed as the complement of the propensity scores for the stratified treatment. After the subclass-propensity scores have been computed, the standard propensity score weighting formulas are used to compute the unstabilized MMWS weights. To estimate MMWS weights equivalent to those described in Hong (2010, 2012), \code{stabilize} must be set to \code{TRUE}, but, as with standard propensity score weights, this is optional. Note that MMWS weights are also known as fine stratification weights and described by Desai et al. (2017).

\code{get_w_from_ps} is not compatible with continuous treatments.
}
\value{
A vector of weights. When \code{subclass} is not \code{NULL}, the subclasses are returned as the \code{"subclass"} attribute. When \code{estimand = "ATOS"}, the chosen value of \code{alpha} (the smallest propensity score allowed to remain in the sample) is returned in the \code{"alpha"} attribute.
}

\author{
Noah Greifer
}
\references{
\strong{Binary treatments}

 - \code{estimand = "ATO"}

Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates via propensity score weighting. Journal of the American Statistical Association, 113(521), 390–400. \doi{10.1080/01621459.2016.1260466}

 - \code{estimand = "ATM"}

Li, L., & Greene, T. (2013). A Weighting Analogue to Pair Matching in Propensity Score Analysis. The International Journal of Biostatistics, 9(2). \doi{10.1515/ijb-2012-0030}

 - \code{estimand = "ATOS"}

Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009). Dealing with limited overlap in estimation of average treatment effects. Biometrika, 96(1), 187–199. \doi{10.1093/biomet/asn055}

 - Other estimands

Austin, P. C. (2011). An Introduction to Propensity Score Methods for Reducing the Effects of Confounding in Observational Studies. Multivariate Behavioral Research, 46(3), 399–424. \doi{10.1080/00273171.2011.568786}

 - Marginal mean weighting through stratification (MMWS)

Hong, G. (2010). Marginal mean weighting through stratification: Adjustment for selection bias in multilevel data. Journal of Educational and Behavioral Statistics, 35(5), 499–531. \doi{10.3102/1076998609359785}

Desai, R. J., Rothman, K. J., Bateman, B. . T., Hernandez-Diaz, S., & Huybrechts, K. F. (2017). A Propensity-score-based Fine Stratification Approach for Confounding Adjustment When Exposure Is Infrequent: Epidemiology, 28(2), 249–257. \doi{10.1097/EDE.0000000000000595}

\strong{Multinomial Treatments}

 - \code{estimand = "ATO"}

Li, F., & Li, F. (2019). Propensity score weighting for causal inference with multiple treatments. The Annals of Applied Statistics, 13(4), 2389–2415. \doi{10.1214/19-AOAS1282}

 - \code{estimand = "ATM"}

Yoshida, K., Hernández-Díaz, S., Solomon, D. H., Jackson, J. W., Gagne, J. J., Glynn, R. J., & Franklin, J. M. (2017). Matching weights to simultaneously compare three treatment groups: Comparison to three-way matching. Epidemiology (Cambridge, Mass.), 28(3), 387–395. \doi{10.1097/EDE.0000000000000627}

 - Other estimands

McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. Statistics in Medicine, 32(19), 3388–3414. \doi{10.1002/sim.5753}

 - Marginal mean weighting through stratification

Hong, G. (2012). Marginal mean weighting through stratification: A generalized method for evaluating multivalued and multiple treatments with nonexperimental data. Psychological Methods, 17(1), 44–60. \doi{10.1037/a0024918}
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

#Using MMWS
w4 <- get_w_from_ps(ps, treat = lalonde$treat,
                    estimand = "ATE", subclass = 20,
                    stabilize = TRUE)

#A multi-category example using GBM predicted probabilities
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