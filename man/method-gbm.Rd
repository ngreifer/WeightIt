\name{method_gbm}
\alias{method_gbm}

\title{Propensity Score Weighting Using Generalized Boosted Models}

\description{
This page explains the details of estimating weights from generalized boosted model-based propensity scores by setting \code{method = "gbm"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating propensity scores using generalized boosted modeling and then converting those propensity scores into weights using a formula that depends on the desired estimand. The algorithm involves choosing a balance criterion to optimize so that balance, rather than prediction, is prioritized.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the propensity scores using \code{\link[gbm]{gbm.fit}} and then optimizes balance using \code{\link[cobalt]{col_w_smd}} for standardized mean differences and \code{\link[cobalt]{col_w_ks}} for Kolmogorov-Smirnoff statistics, both from \pkg{cobalt}. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights are computed from the estimated propensity scores using \code{\link{get_w_from_ps}}, which implements the standard formulas.  When \code{include.obj = TRUE}, the returned object is the \code{gbm} fit.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, this method estimates the propensity scores using \code{\link[gbm]{gbm.fit}} with \code{distribution = "multinomial"} and then optimizes balance using \code{\link[cobalt]{col_w_smd}} for standardized mean differences and \code{\link[cobalt]{col_w_ks}} for Kolmogorov-Smirnoff statistics, both from \pkg{cobalt}. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights are computed from the estimated propensity scores using \code{\link{get_w_from_ps}}, which implements the standard formulas. The balance that is optimized is that between each non-focal treatment and the focal treatment for the ATT and ATC and between each treatment and the overall unweighted sample for other estimands. When \code{include.obj = TRUE}, the returned object is the \code{gbm} fit.
}
\subsection{Continuous Treatments}{
For continuous treatments, the generalized propensity score is estimated using \code{\link[gbm]{gbm.fit}}. When \code{include.obj = TRUE}, the returned object is the \code{gbm} fit.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point.
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios.
}
\subsection{Missing Data}{
Missing data is compatible with generalized boosted modeling. \code{NA} values are simply considered their own plausible value for the covariate. When balance statistics are computed to choose the optimal tree, they are computed using only the non-missing values for the variable in question.
}
}
\section{Additional Arguments}{

The following additional arguments can be specified:
\describe{
\item{\code{stop.method}}{A string describing the criterion used to select the best weights. This has two parts: a statistics and a summarizing function, which should be combined as "\{stat\},\{summary\}". For binary and multinomial treatments, the available stats are \code{"es"} for absolute standardized mean differences and \code{"ks"} for Kolmogorov-Smirnoff statistics; for continuous treatments, the available stats are \code{"p"} for Pearson correlations between each covariate and the treatment and \code{"s"} for Spearman correlations. For all treatment types, the available summaries are \code{"mean"} for the mean of the statistics, \code{"max"} for the maximum of the statistics, and \code{"rms"} for the root mean square of the statistics (i.e., the L2 norm). The default for binary and multinomial treatments is \code{"es.mean"} and the default for continuous treatments is \code{"p.mean"}.
}
\item{\code{trim.at}}{A number supplied to \code{at} in \code{\link{trim}} which trims the weights from all the trees before choosing the best tree. This can be valuable when some weights are extreme, which occur especially with continuous treatment. The default is 0 (i.e., no trimming).
}
\item{\code{distribution}}{A string with the distribution used in the loss function of the boosted model. This is supplied to the \code{distribution} argument in \code{\link[gbm]{gbm.fit}}. For binary treatments, \code{"bernoulli"} and \code{"adaboost"} are available, with \code{"bernoulli"} the default. For multinomial treatments, only \code{"multinomial"} is allowed. For continuous treatments \code{"gaussian"}, \code{"laplace"}, and \code{"tdist"} are available, with \code{"gaussian"} the default.
}
\item{\code{n.trees}}{The maximum number of trees used. This is passed onto the \code{n.trees} argument in \code{gbm.fit}. The default is 10000 for binary and multinomial treatments and 20000 for continuous treatments.
}
\item{\code{interaction.depth}}{The depth of the trees. This is passed onto the \code{interaction.depth} argument in \code{gbm.fit}. Higher values indicate better ability to capture nonlinear and nonadditive relationships. The default is 3 for binary and multinomial treatments and 4 for continuous treatments.
}
\item{\code{shrinkage}}{The shrinkage parameter applied to the trees. This is passed onto the \code{shrinkage} argument in \code{gbm.fit}. The default is .01 for binary and multinomial treatments and .0005 for continuous treatments.
}
\item{\code{bag.fraction}}{The fraction of the units randomly selected to propose the next tree in the expansion. This is passed onto the \code{bag.fraction} argument in \code{gbm.fit}. The default is 1, but smaller values should be tried.
}
}

All other arguments take on the defaults of those in \code{\link[gbm]{gbm.fit}}, and some are not used at all.

The \code{w} argument in \code{gbm.fit} is ignored because sampling weights are passed using \code{s.weights}.

}
\references{
\strong{Binary treatments}

McCaffrey, D. F., Ridgeway, G., & Morral, A. R. (2004). Propensity Score Estimation With Boosted Regression for Evaluating Causal Effects in Observational Studies. Psychological Methods, 9(4), 403–425. \doi{10.1037/1082-989X.9.4.403}

\strong{Multinomial Treatments}

McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. Statistics in Medicine, 32(19), 3388–3414. \doi{10.1002/sim.5753}

\strong{Continuous treatments}

Zhu, Y., Coffman, D. L., & Ghosh, D. (2015). A Boosting Algorithm for Estimating Generalized Propensity Scores with Continuous Treatments. Journal of Causal Inference, 3(1).

}
\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}, \code{\link{method_twang}}
}
\note{
This method used to implement what is now implemented by \code{method = "twang"}. Using \code{method = "gbm"} is \emph{much} faster than using \code{method = "twang"} or \code{twang::ps} because of efficiency improvements that are not present in \pkg{twang}, but the results will differ slightly between the two methods. Notably, for multinomial treatments, the methods use different strategies: \code{method = "gbm"} estimates multinomial propensity scores (one for each individual for each treatment level) and chooses the tree that optimizes balance between each treatment and the target population, while \code{method = "twang"} estimates a set of binary propensity scores for each treatment, each optimizing balance for the treatment and the target population, and then combines the weights. For continuous treatments, the results will be identical.
}
\examples{
# Examples take a long time to run
\dontrun{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "gbm", estimand = "ATE",
                stop.method = "es.max"))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "gbm", estimand = "ATT",
                focal = "hispan", stop.method = "ks.mean"))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "gbm", use.kernel = TRUE,
                stop.method = "p.rms", trim.at = .97))
summary(W3)
bal.tab(W3)
}
}