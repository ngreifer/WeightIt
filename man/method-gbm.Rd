\name{method_gbm}
\alias{method_gbm}

\title{Propensity Score Weighting Using Generalized Boosted Models}

\description{
This page explains the details of estimating weights from generalized boosted model-based propensity scores by setting \code{method = "gbm"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating propensity scores using generalized boosted modeling and then converting those propensity scores into weights using a formula that depends on the desired estimand. The algorithm involves using a balance-based or prediction-based criterion to optimize in choosing the value of a tuning parameter (the number of trees). This method mimics the functionality of functions in the \pkg{twang} package, but has improved performance and more flexible options. See Note for more details.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the propensity scores using \code{\link[gbm]{gbm.fit}} and then optimizes balance using \code{\link[cobalt]{col_w_smd}} for standardized mean differences and \code{\link[cobalt]{col_w_ks}} for Kolmogorov-Smirnov statistics, both from \pkg{cobalt}. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights are computed from the estimated propensity scores using \code{\link{get_w_from_ps}}, which implements the standard formulas. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See \code{\link{get_w_from_ps}} for details.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, this method estimates the propensity scores using \code{\link[gbm]{gbm.fit}} with \code{distribution = "multinomial"} and then optimizes balance using \code{\link[cobalt]{col_w_smd}} for standardized mean differences and \code{\link[cobalt]{col_w_ks}} for Kolmogorov-Smirnov statistics, both from \pkg{cobalt}. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights are computed from the estimated propensity scores using \code{\link{get_w_from_ps}}, which implements the standard formulas. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See \code{\link{get_w_from_ps}} for details. The balance that is optimized is that between each non-focal treatment and the focal treatment for the ATT and ATC, between each treatment and the overall unweighted sample for the ATE, and between each treatment and the overall weighted sample for other estimands.
}
\subsection{Continuous Treatments}{
For continuous treatments, this method estimates the generalized propensity score using \code{\link[gbm]{gbm.fit}} and then optimizes balance using \code{\link[cobalt]{col_w_corr}} for treatment-covariate correlations from \pkg{cobalt}.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point.
}
\subsection{Sampling Weights}{
Sampling weights are supported through \code{s.weights} in all scenarios.
}
\subsection{Missing Data}{
In the presence of missing data, the following value(s) for \code{missing} are allowed:
\describe{
\item{\code{"ind"} (default)}{
First, for each variable with missingness, a new missingness indicator variable is created that takes the value 1 if the original covariate is \code{NA} and 0 otherwise. The missingness indicators are added to the model formula as main effects. The weight estimation then proceeds with this new formula and set of covariates using surrogate splitting as described below. The covariates output in the resulting \code{weightit} object will be the original covariates with the \code{NA}s.
}
\item{\code{"surr"}}{
Surrogate splitting is used to process \code{NA}s. No missingness indicators are created. Nodes are split using only the non-missing values of each variable. To generate predicted values for each unit, a non-missing variable that operates similarly to the variable with missingness is used as a surrogate. Missing values are ignored when calculating balance statistics to choose the optimal tree.
}
}
}
}
\section{Additional Arguments}{

The following additional arguments can be specified:
\describe{
\item{\code{stop.method}}{A string describing the criterion used to select the best weights. When optimizing for balance, this has two parts, a statistic to be computed and a summarizing function, which should be combined as \code{"{stat}.{summary}"}. For binary and multinomial treatments, the available stats are \code{"es"} for absolute standardized mean differences and \code{"ks"} for Kolmogorov-Smirnov statistics; for continuous treatments, the available stats are \code{"p"} for Pearson correlations between each covariate and the treatment and \code{"s"} for Spearman correlations. For all treatment types, the available summaries are \code{"mean"} for the mean of the statistics, \code{"max"} for the maximum of the statistics, and \code{"rms"} for the root mean square of the statistics (i.e., the L2 norm). The default for binary and multinomial treatments is \code{"es.mean"} and the default for continuous treatments is \code{"p.mean"}.

In addition, to optimize the cross-validation error, \code{stop.method} can be set as \code{"cv{#}"}, where \code{{#}} is replaced by a number representing the number of cross-validation folds used (e.g., \code{"cv5"} for 5-fold cross-validation).
}
\item{\code{trim.at}}{A number supplied to \code{at} in \code{\link{trim}} which trims the weights from all the trees before choosing the best tree. This can be valuable when some weights are extreme, which occurs especially with continuous treatments. The default is 0 (i.e., no trimming).
}
\item{\code{distribution}}{A string with the distribution used in the loss function of the boosted model. This is supplied to the \code{distribution} argument in \code{\link[gbm]{gbm.fit}}. For binary treatments, \code{"bernoulli"} and \code{"adaboost"} are available, with \code{"bernoulli"} the default. For multinomial treatments, only \code{"multinomial"} is allowed. For continuous treatments \code{"gaussian"}, \code{"laplace"}, and \code{"tdist"} are available, with \code{"gaussian"} the default.
}
\item{\code{n.trees}}{The maximum number of trees used. This is passed onto the \code{n.trees} argument in \code{gbm.fit}. The default is 10000 for binary and multinomial treatments and 20000 for continuous treatments.
}
\item{\code{start.tree}}{The tree at which to start balance checking. If you know the best balance isn't in the first 100 trees, for example, you can set \code{start.tree = 101} so that balance statistics are not computed on the first 100 trees. This can save some time since balance checking takes up the bulk of the run time for balance-based stopping methods, and is especially useful when running the same model adding more and more trees. The default is 1, i.e., to start from the very first tree in assessing balance.
}
\item{\code{interaction.depth}}{The depth of the trees. This is passed onto the \code{interaction.depth} argument in \code{gbm.fit}. Higher values indicate better ability to capture nonlinear and nonadditive relationships. The default is 3 for binary and multinomial treatments and 4 for continuous treatments.
}
\item{\code{shrinkage}}{The shrinkage parameter applied to the trees. This is passed onto the \code{shrinkage} argument in \code{gbm.fit}. The default is .01 for binary and multinomial treatments and .0005 for continuous treatments. The lower this value is, the more trees one may have to include to reach the optimum.
}
\item{\code{bag.fraction}}{The fraction of the units randomly selected to propose the next tree in the expansion. This is passed onto the \code{bag.fraction} argument in \code{gbm.fit}. The default is 1, but smaller values should be tried.
}
}

All other arguments take on the defaults of those in \code{\link[gbm]{gbm.fit}}, and some are not used at all.

The \code{w} argument in \code{gbm.fit} is ignored because sampling weights are passed using \code{s.weights}.

For continuous treatments only, the following arguments may be supplied:
\describe{
\item{\code{density}}{A function corresponding the conditional density of the treatment. The standardized residuals of the treatment model will be fed through this function to produce the numerator and denominator of the generalized propensity score weights. If blank, \code{\link{dnorm}} is used as recommended by Robins et al. (2000). This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if \code{density = "dt_2"} is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme outcome values are observed (Naimi et al., 2014). Ignored if \code{use.kernel = TRUE} (described below).
}
\item{\code{use.kernel}}{If \code{TRUE}, uses kernel density estimation through the \code{\link{density}} function to estimate the numerator and denominator densities for the weights. If \code{FALSE} (the default), the argument to the \code{density} parameter is used instead.
}
\item{\code{bw}, \code{adjust}, \code{kernel}, \code{n}}{If \code{use.kernel = TRUE}, the arguments to the \code{\link{density}} function. The defaults are the same as those in \code{density} except that \code{n} is 10 times the number of units in the sample.
}
}
}
\section{Additional Outputs}{
\describe{
\item{\code{info}}{
A list with two entries:
\describe{
\item{\code{best.tree}}{
The number of trees at the optimum. If this is close to \code{n.trees}, \code{weightit} should be rerun with a larger value for \code{n.trees}, and \code{start.tree} can be set to just below \code{best.tree}. See example.
}
\item{\code{tree.val}}{
A data.frame with two columns: the first is the number of trees and the second is the value of the criterion corresponding to that tree. Running \code{\link{plot}} on this object will plot the criterion by the number of trees and is a good way to see patterns in the relationship between them and to determine if more trees are needed. See example.
}
}
}
\item{\code{obj}}{
When \code{include.obj = TRUE}, the \code{gbm} fit used to generate the predicted values.
}
}
}
\references{
\strong{Binary treatments}

McCaffrey, D. F., Ridgeway, G., & Morral, A. R. (2004). Propensity Score Estimation With Boosted Regression for Evaluating Causal Effects in Observational Studies. Psychological Methods, 9(4), 403–425. \doi{10.1037/1082-989X.9.4.403}

\strong{Multinomial Treatments}

McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand, R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for Multiple Treatments Using Generalized Boosted Models. Statistics in Medicine, 32(19), 3388–3414. \doi{10.1002/sim.5753}

\strong{Continuous treatments}

Zhu, Y., Coffman, D. L., & Ghosh, D. (2015). A Boosting Algorithm for Estimating Generalized Propensity Scores with Continuous Treatments. Journal of Causal Inference, 3(1). \doi{10.1515/jci-2014-0022}

}

\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}, \code{\link{method_twang}}
}
\note{
This method used to implement what is now implemented by \code{method = "twang"}. Using \code{method = "gbm"} is \emph{much} faster than using \code{method = "twang"} or \code{twang::ps} because of efficiency improvements that are not present in \pkg{twang}, but the results will differ slightly between the two methods. Notably, for multinomial treatments, the methods use different strategies: \code{method = "gbm"} estimates multinomial propensity scores (one for each individual for each treatment level) and chooses the tree that optimizes balance between each treatment and the target population, while \code{method = "twang"} estimates a set of binary propensity scores for each treatment, each optimizing balance for the treatment and the target population, and then combines the weights. For continuous treatments, the results will be identical because both use the same architecture found in \code{\link{ps.cont}}.

When standardized mean differences are used for the stopping method, standardized mean differences are computed for binary variables as well as for continuous variables, unlike in \pkg{cobalt}, which computes raw differences in proportion for binary variables. In addition, when \code{estimand = "ATE"}, \code{s.d.denom} is set to "all" in the call to \code{\link[cobalt]{col_w_smd}}, which mirrors the method in \pkg{twang} but differs from the default in balance checking with \code{\link[cobalt]{bal.tab}} in \pkg{cobalt}, which typically uses \code{s.d.denom = "pooled"} in this case.
}
\examples{

library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "gbm", estimand = "ATE",
                stop.method = "es.max"))
summary(W1)
bal.tab(W1)

\dontrun{
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

#Using a t(3) density and illustrating the search for
#more trees.
W4a <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "gbm", density = "dt_3",
                stop.method = "p.max",
                n.trees = 10000)

W4a$info$best.tree #10000; optimum hasn't been found
plot(W4a$info$tree.val) #decreasing at right edge

W4b <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "gbm", density = "dt_3",
                stop.method = "p.max",
                start.tree = 10000,
                n.trees = 20000)

W4b$info$best.tree #13417; optimum has been found
plot(W4b$info$tree.val) #increasing at right edge

bal.tab(W4b)

}
}