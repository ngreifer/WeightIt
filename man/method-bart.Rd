\name{method_bart}
\alias{method_bart}

\title{Propensity Score Weighting Using BART}

\description{
This page explains the details of estimating weights from Bayesian additive regression trees (BART)-based propensity scores by setting \code{method = "bart"} in the call to \code{\link{weightit}} or \code{\link{weightitMSM}}. This method can be used with binary, multinomial, and continuous treatments.

In general, this method relies on estimating propensity scores using BART and then converting those propensity scores into weights using a formula that depends on the desired estimand. This method mainly relies on \code{\link[BART]{gbart}} from the \CRANpkg{BART} package.

\subsection{Binary Treatments}{
For binary treatments, this method estimates the propensity scores using \code{\link[BART]{gbart}}. The following estimands are allowed: ATE, ATT, ATC, ATO, and ATM. The weights for the ATE, ATT, and ATC are computed from the estimated propensity scores using the standard formulas, the weights for the ATO are computed as in Li & Li (2018), and the weights for the ATM (i.e., average treatment effect in the equivalent sample "pair-matched" with calipers) are computed as in Yoshida et al. (2017). Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See \code{\link{get_w_from_ps}} for details.
}
\subsection{Multinomial Treatments}{
For multinomial treatments, this method estimates the propensity scores using \code{\link[BART]{mbart}}. The following estimands are allowed: ATE, ATT, ATO, and ATM. The weights for each estimand are computed using the standard formulas or those mentioned above. Weights can also be computed using marginal mean weighting through stratification for the ATE, ATT, and ATC. See \code{\link{get_w_from_ps}} for details.
}
\subsection{Continuous Treatments}{
For continuous treatments, the generalized propensity score is estimated using \code{\link[BART]{gbart}}. In addition, kernel density estimation can be used instead of assuming a normal density for the numerator and denominator of the generalized propensity score by setting \code{use.kernel = TRUE}. Other arguments to \code{\link{density}} can be specified to refine the density estimation parameters. \code{plot = TRUE} can be specified to plot the density for the numerator and denominator, which can be helpful in diagnosing extreme weights.
}
\subsection{Longitudinal Treatments}{
For longitudinal treatments, the weights are the product of the weights estimated at each time point.
}
\subsection{Sampling Weights}{
Sampling weights are not supported.
}
\subsection{Missing Data}{
In the presence of missing data, the following value(s) for \code{missing} are allowed:
\describe{
\item{\code{"ind"} (default)}{
First, for each variable with missingness, a new missingness indicator variable is created which takes the value 1 if the original covariate is \code{NA} and 0 otherwise. The missingness indicators are added to the model formula as main effects. The missing values in the covariates are then replaced with 0s. The weight estimation then proceeds with this new formula and set of covariates. The covariates output in the resulting \code{weightit} object will be the original covariates with the \code{NA}s.
}
\item{\code{"hot"}}{
The function that fits the BART model (i.e., \code{gbart} or \code{mbart}) performs hot-deck imputation, which involves randomly selecting a value to impute the missing covariate with from another unit in the dataset with a nonmissing value for that variable. The \pkg{BART} authors recommend only doing this for small amounts of missingness. Doing so will slightly weaken the estimated relationship between the covariate with missingness and the outcome (what the authors describe as being "biased toward the null").
}
}
}
}
\section{Additional Arguments}{
\describe{
\item{\code{link}}{The link used in the BART model for the propensity scores with binary or multinomial treatments. The allowable options are \code{"logit"} and \code{"probit"}. \code{"probit"} is much faster and is the default; there is generally little reason to choose \code{"logit"}. This option is used to supply the argument to \code{type} in \code{gbart} and \code{mbart}.}
\item{\code{use.mbart2}}{For multinomial treatments, whether to use \code{\link[BART]{mbart2}} or \code{\link[BART]{mbart}} to fit the BART model for the propensity score. \code{mbart} uses a conditional probability specification to fit the models; first the probability of one category versus all the others is modeled, then the probability of the second category versus all the other is model using only those who are not in the first category, and so on. \code{mbart2} uses a series of binary BART models for each category versus the others. The package authors recommend using \code{mbart} for "few" categories, so this is the default (i.e., \code{use.mbart2 = FALSE}). \code{mbart} tends to be faster.}
\item{\code{mc.cores}}{The number of cores to enage when performing parallelization (not available on Windows). If greater than 1, the multicore version of the corresponding \pkg{BART} function will be used (e.g., \code{\link[BART]{mc.gbart}} instead of \code{gbart} and \code{\link[BART]{mc.mbart}} instead of \code{mbart}). Parallelization requires no extra effort or knowledge to perform and can dramatically speed up model fitting, so it is highly recommended. It's as simple as setting \code{mc.cores = 4} in the call to \code{weightit}. This value is passed on to the \code{mc.cores} argument in the fitting function. The default is 1 for no parallelization. When using parallelization, an argument should be passed to \code{seed} to ensure the estimation is replicable; setting the seed outside the function will not always work. \code{seed} should probably be set differently each time \code{weightit} is called if replicability is not important because otherwise it will remain at its default value of 99 every time, which guarantees replicability but also dependence on that particular (arbitrary) value.}
}

All arguments to the \pkg{BART} fitting functions can be passed through \code{weightit} or \code{weightitMSM}, with the following exceptions:
\itemize{
\item{\code{type} and \code{ntype} are ignored because the type is determined by the form of the treatment variable and the argument supplied to \code{link}.}
}

The \code{ntree} argument can optionally be supplied as \code{n.trees} to be consistent with \code{method = "gbm"}.

For continuous treatments only, the following arguments may be supplied:
\describe{
\item{\code{density}}{A function corresponding to the conditional density of the treatment. The standardized residuals of the treatment model will be fed through this function to produce the numerator and denominator of the generalized propensity score weights. If blank, \code{\link{dnorm}} is used as recommended by Robins et al. (2000). This can also be supplied as a string containing the name of the function to be called. If the string contains underscores, the call will be split by the underscores and the latter splits will be supplied as arguments to the second argument and beyond. For example, if \code{density = "dt_2"} is specified, the density used will be that of a t-distribution with 2 degrees of freedom. Using a t-distribution can be useful when extreme outcome values are observed (Naimi et al., 2014). Ignored if \code{use.kernel = TRUE} (described below).
}
\item{\code{use.kernel}}{If \code{TRUE}, uses kernel density estimation through the \code{\link{density}} function to estimate the numerator and denominator densities for the weights. If \code{FALSE}, the argument to the \code{density} parameter is used instead.
}
\item{\code{bw}, \code{adjust}, \code{kernel}, \code{n}}{If \code{use.kernel = TRUE}, the arguments to the \code{\link{density}} function. The defaults are the same as those in \code{density} except that \code{n} is 10 times the number of units in the sample.
}
\item{\code{plot}}{If \code{use.kernel = TRUE}, whether to plot the estimated density.
}
}

}
\section{Additional Outputs}{
\describe{
\item{\code{obj}}{
When \code{include.obj = TRUE}, the BART fit(s) used to generate the predicted values. The specific function used to fit the model depends on the inputs.
}
}
}
\details{
BART works by fitting a sum-of-trees model for the treatment or probability of treatment. The number of trees is determiend by the \code{ntree} argument. Bayesian priors are used for the hyperparameters, so the result is a posterior distribution of predicted values for each unit. The mean of these for each unit is taken as the (generalized) propensity score. Although the hyperparameters governing the priors can be modified by supplying arguments to \code{weightit} that are passed to the BART fitting function, the default values tend to work well and require little modification (though the defaults differ for continuous and categorical treatments; see the \code{\link[BART]{gbart}} documentation for details). Unlike many other machine learning methods, no loss function is optimized and the hyperparameters do not need to be tuned (e.g., using cross-validation). BART tends to balance sparseness with flexibility by using very weak learners as the trees, which makes it suitable for capturing complex functions without specifying a particular functional form beyond the link and without overfitting.
}
\references{

Hill, J., Weiss, C., & Zhai, F. (2011). Challenges With Propensity Score Strategies in a High-Dimensional Setting and a Potential Alternative. Multivariate Behavioral Research, 46(3), 477–513. \doi{10.1080/00273171.2011.570161}

Chipman, H. A., George, E. I., & McCulloch, R. E. (2010). BART: Bayesian additive regression trees. The Annals of Applied Statistics, 4(1), 266–298. \doi{10.1214/09-AOAS285}

Note that many references that deal with BART for causal inference focus on estimating potential outcomes with BART, not the propensity scores, and so are not directly relevant when using BART to estimate propensity scores for weights.

See \code{\link{method_ps}} for additional references.
}
\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}, \code{\link{get_w_from_ps}}

\code{\link{method_super}} for stacking predictions from several machine learning methods, including BART.
}
\examples{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "bart", estimand = "ATT",
                mc.cores = 4))
summary(W1)
bal.tab(W1)

#Balancing covariates with respect to race (multinomial)
(W2 <- weightit(race ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "bart", estimand = "ATE",
                mc.cores = 4))
summary(W2)
bal.tab(W2)

#Balancing covariates with respect to re75 (continuous)
#assuming t(3) conditional density for treatment
(W3 <- weightit(re75 ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = "bart", density = "dt_3",
                mc.cores = 4))
summary(W3)
bal.tab(W3)
}