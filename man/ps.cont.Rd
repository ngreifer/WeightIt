\name{ps.cont}
\alias{ps.cont}
\alias{summary.ps.cont}
\alias{plot.ps.cont}
\alias{boxplot.ps.cont}
\title{
Generalized Propensity Score Estimation using GBM
}
\description{
\code{ps.cont} calculates generalized propensity scores and corresponding weights using boosted linear regression as implemented in \code{\link[gbm]{gbm}}. This function extends \code{\link[twang]{ps}} in \pkg{twang} to continuous treatments. The syntax and output are largely the same. The GBM parameter defaults are those found in Zhu, Coffman, & Ghosh (2015).

Note: \code{ps.cont} will phased out when \pkg{twang} adds functionality for continuous treatments. All functionality and more is already present in \code{\link{weightit}} with \code{\link{method_gbm}[method = "gbm"]}.
}
\usage{
ps.cont(formula, data,
        n.trees = 20000,
        interaction.depth = 4,
        shrinkage = 0.0005,
        bag.fraction = 1,
        print.level = 0,
        verbose = FALSE,
        stop.method,
        sampw = NULL,
        optimize = 1,
        use.kernel = FALSE,
        ...)
\method{summary}{ps.cont}(object, ...)
\method{plot}{ps.cont}(x, ...)
\method{boxplot}{ps.cont}(x, ...)

}
\arguments{
  \item{formula}{
A formula for the propensity score model with the treatment indicator on the left side of the formula and the potential confounding variables on the right side.
}
  \item{data}{
The dataset in the form of a data frame, which should include treatment assignment as well as the covariates specified in \code{formula}.
}
  \item{n.trees}{
The number of GBM iterations passed on to \code{\link[gbm]{gbm}}. The more, the better the final solution will be, but the more time it will take.
}
  \item{interaction.depth}{
The \code{interaction.depth} passed on to \code{\link[gbm]{gbm}}.
}
\item{shrinkage}{
The \code{shrinkage} passed on to \code{\link[gbm]{gbm}}.
}
\item{bag.fraction}{
The \code{bag.fraction} passed on to \code{\link[gbm]{gbm}}.
}
  \item{print.level}{
Currently ignored.
}
  \item{verbose}{
If \code{TRUE}, information will be printed to monitor the the progress of the fitting.
}
\item{stop.method}{A method or methods of measuring and summarizing balance across pretreatment variables.  Current options are \code{p.max}, \code{p.mean}, \code{p.rms}, \code{s.max}, \code{s.mean}, and \code{s.rms}. \code{p} refers to the Pearson correlation and \code{s} refers to the Spearman correlation. These are summarized across the pretreatment variables by the maximum (\code{max}), the mean (\code{mean}), or the square root of the mean of the squares (\code{rms}).
}
\item{sampw}{Optional sampling weights.
}
\item{optimize}{
A numeric value, either \code{0}, \code{1}, or \code{2}. If \code{0}, balance will be checked for every tree, and the tree with the best balance will be the one used to generate the final weights. If \code{1}, the default, balance will be checked for a subset of trees, and then \code{\link{optimize}} will be used to find the tree with the best balance within the tree interval chosen. If \code{2}, \code{optimize} will be used to find the tree that yields the best balance. \code{0} takes the longest but is guaranteed to find the best balance among the trees. \code{2} is the quickest but will often choose a tree that that suboptimal balance, though not by much. \code{1} is a compromise between speed and comprehensiveness and is the algorithm implemented in \pkg{twang}.
}
\item{use.kernel}{
Whether to use kernel density estimation as implemented in \code{\link{density}} to estimate the numerator of the weights. If \code{TRUE}, \code{density} will be used. If \code{FALSE}, the default, a normal density will be assumed and will be estimated using \fun{dnorm}.
}
\item{object, x}{
A \code{ps.cont} object.
}
\item{\dots}{
For \code{ps.cont}, if \code{use.density = TRUE}, additional arguments to \code{\link{density}}, which is used to produce the density for the numerator of the weights. These include \code{bw}, \code{adjust}, \code{kernel}, and \code{n}. The default values are the defaults for \code{density}, except \code{n}, which is 10 times the number of units.

For \code{summary.ps.cont}, additional arguments affecting the summary produced.}
}
\details{
\code{ps.cont} extends \code{ps} in \pkg{twang} to continuous treatments. It estimates weights from a series of trees and then outputs the weights that optimize a user-set criterion. The criterion employed involves the correlation between the treatment and each covariate. In a fully balanced sample, the treatment will have a correlation of 0 with covariates sufficient for removing confounding. Zhu, Coffman, & Ghosh (2015), who were the first to describe GBM for propensity score weighting with continuous treatments, recommend this procedure and provided R code to implement the methods they describe. \code{ps.cont} adapts their syntax to make it consistent with that of \code{ps} in \pkg{twang}. As in Zhu et al. (2015), when the Pearson correlation is requested, weighted biserial correlations will be computed for binary covariates.

The weights are estimated as the marginal density of the treatment divided by the conditional density of the treatment on the covariates for each unit. For the marginal density, a kernel density estimator can be implemented using the \code{\link{density}} function. For the conditional density, a Gaussian density is assumed. Note that with treatment with outlying values, extreme weights can be produced, so it is important to examine the weights and trim them if necessary.

It is recommended to use as many trees as possible, though this requires more computation time, especially with \code{use.optimize} set to \code{0}. There is little difference between using Pearson and Spearman correlations or between using the raw correlations and the Z-transformed correlations. Typically the only \code{gbm}-related options that should be changed are the interaction depth and number of trees.

Missing data is not allowed in the covariates because of the ambiguity in computing correlations with missing values.

\code{summary.ps.cont} compresses the information in the \code{desc} component of the \code{ps.cont} object into a short summary table describing the size of the dataset and the quality of the generalized propensity score weights, in a similar way to \code{\link[twang]{summary.ps}}.

\code{plot.ps.cont} and \code{boxplot.ps.cont} function almost identically to \code{\link[twang]{plot.ps}} and \code{\link[twang]{boxplot.ps}}. See the help pages there for more information. Note that for \code{plot.ps}, only options 1, 2, and 6 are available for the \code{plots} argument. When \code{use.optimize = 2}, option 1 is not available.

}
\value{
Returns an object of class \code{ps} and \code{ps.cont}, a list containing
  \item{gbm.obj}{The returned \code{\link[gbm]{gbm}} object.}
  \item{treat}{The treatment variable.}
   \item{desc}{a list containing balance tables for each method selected in
              \code{stop.method}. Includes a component for the unweighted
              analysis names \dQuote{unw}. Each \code{desc} component includes
              a list with the following components:
  \describe{
     \item{ess}{The effective sample size}
     \item{n}{The number of subjects}
     \item{max.p.cor}{The largest absolute Pearson correlation across the covariates}
     \item{mean.p.cor}{The mean absolute Pearson correlation of the covariates}
     \item{rmse.p.cor}{The root mean squared Pearson correlation across the covariates}
     \item{max.s.cor}{The largest absolute Spearman correlation across the covariates}
     \item{mean.s.cor}{The mean absolute Spearman correlation of the covariates}
     \item{rmse.s.cor}{The root mean squared Spearman correlation across the covariates}
     \item{bal.tab}{a table summarizing the quality of the weights for yielding low treatment-covariate correlations. This table is best extracted using \code{\link[twang]{bal.table}}.}
     \item{n.trees}{The estimated optimal number of \code{\link[gbm]{gbm}} iterations to optimize the loss function for the associated \code{stop.method}s}
   }
}
  \item{ps}{a data frame containing the estimated generalized propensity scores. Each column is associated with one of the methods selected in \code{stop.methods}.}
  \item{w}{a data frame containing the propensity score weights. Each column is associated with one of the methods selected in \code{stop.methods}. If sampling weights are given then these are incorporated into the weights.}
\item{estimand}{\code{NULL}}
  \item{datestamp}{Records the date of the analysis.}
  \item{parameters}{Saves the \code{ps.cont} call.}
  \item{alerts}{\code{NULL}}
\item{iters}{A sequence of iterations used in the GBM fits used by \code{plot.ps.cont}.}
\item{balance}{The balance summary for each tree examined, with a column for each stop.method. If \code{optimize = 0}, this will contain balance summaries for all trees. If \code{optimize = 1}, this will contain balance summaries for the subset of trees corresponding to \code{iters}. If \code{optimize = 2}, this will be NULL.}
\item{n.trees}{Maximum number of trees considered in GBM fit.}
\item{data}{Data as specified in the \code{data} argument.}

The \code{NULL} entries exist so the output object is similar to that of \code{ps} in \pkg{twang}.
}

\references{
Zhu, Y., Coffman, D. L., & Ghosh, D. (2015). A Boosting Algorithm for Estimating Generalized Propensity Scores with Continuous Treatments. Journal of Causal Inference, 3(1). \doi{10.1515/jci-2014-0022}

}
\author{
Noah Greifer

\code{ps.cont} is heavily adapted from the R code in Zhu, Coffman, & Ghosh (2015). In contrast with their code, \code{ps.cont} uses weighted Pearson and Spearman correlations rather than probability weighted bootstrapped correlations, allows for different degrees of optimization in searching for the best solution, and allows for the use of kernel density estimation for the generalized propensity score. \code{ps.cont} also takes inspiration from \code{ps} in \pkg{twang}.
}

\seealso{
\code{\link{weightit}} and \code{\link{method_gbm}} for its implementation using \code{weightit} syntax.

\code{\link[twang]{ps}} and \code{\link[twang]{mnps}} for GBM with binary and multinomial treatments.

\code{\link[gbm]{gbm}} for the underlying machinery and explanation of the parameters.
}
\examples{
# Examples take a long time
\dontrun{
library("cobalt")
data("lalonde", package = "cobalt")

#Balancing covariates with respect to re75
psc.out <- ps.cont(re75 ~ age + educ + married +
                nodegree + race + re74, data = lalonde,
                stop.method = c("p.mean", "p.max"),
                use.optimize = 2)
summary(psc.out)
twang::bal.table(psc.out) #twang's bal.table
}
}
