\name{weightitMSM}
\alias{weightitMSM}
\alias{print.weightitMSM}
\title{Generate Balancing Weights}
\usage{
weightitMSM(formula.list,
            data = NULL,
            method = "ps",
            stabilize = FALSE,
            by = NULL,
            s.weights = NULL,
            num.formula = NULL,
            moments = 1,
            int = FALSE,
            verbose = FALSE,
            include.obj = FALSE,
            is.MSM.method,
            weightit.force = FALSE,
            ...)

\method{print}{weightitMSM}(x, ...)

}
\description{
\code{weightitMSM} allows for the easy generation of balancing weights for marginal structural models for time-varying treatments using a variety of available methods for binary, continuous, and multinomial treatments. Many of these methods exist in other packages, which \code{\link{weightit}} calls; these packages must be installed to use the desired method. Also included are \code{print} and \code{summary} methods for examining the output.

Currently only "wide" data sets, where each row corresponds to a unit's entire variable history, are supported. You can use \code{\link{reshape}} or other functions to transform your data into this format; see example below.
}
\arguments{
  \item{formula.list}{
a list of formulas corresponding to each time point with the time-specific treatment variable on the left hand side and pre-treatment covariates to be balanced on the right hand side. The formulas must be in temporal order, and must contain all covariates to be balanced at that time point (i.e., treatments and covariates featured in early formulas should appear in later ones). Interactions and functions of covariates are allowed.
}
  \item{data}{
an optional data set in the form of a data frame that contains the variables in the formulas in \code{formula.list}. This must be a wide data set with exactly one row per unit.
}
  \item{method}{
a string of length 1 containing the name of the method that will be used to estimate weights. See \code{\link{weightit}} for allowable options. The default is \code{"ps"}, which estimates the weights using generalized linear models.
}
  \item{stabilize}{
\code{logical}; whether or not to stabilize the weights. Stabilizing the weights involves fitting a model predicting treatment at each time point from treatment status at prior time points. If \code{TRUE}, a saturated model will be fit, essentially using the observed treatment probabilities in the numerator (for binary and multinomial treatments). This may yield an error if some combinations are not observed. Default is \code{FALSE}. To manually specify stabilization model formulas, use \code{num.formula}.
}
\item{num.formula}{
optional; a one-sided formula with the stabilization factors (other than the previous treatments) on the right hand side, which adds, for each time point, the stabilization factors to a model saturated with previous treatments. See Cole & Hernán (2008) for a discussion of how to specify this model; including stabilization factors can change the estimand without proper adjustment, and should be done with caution. Unless you know what you are doing, we recommend setting \code{stabilize = TRUE} and ignoring \code{num.formula}.
}
  \item{by}{
a vector or the name of the variable in data for which weighting is to be done within categories. For example, if \code{by = "gender"}, weights will be generated separately within each level of the variable \code{"gender"}. The argument used to be called \code{exact}, which will still work but with a message. Only one \code{by} variable is allowed.
}
  \item{s.weights}{
a vector of sampling weights or the name of a variable in \code{data} that contains sampling weights. These are ignored for some methods.
}
  \item{moments}{
\code{numeric}; for some methods, the greatest moment of the covariate distribution to be balanced. For example, if \code{moments = 3}, for all non-categorical covariates, the mean, second moment (variance), and third moments (skew) of the covariates will be balanced. This argument is ignored for other methods; to balance powers of the covariates, appropriate functions must be entered in \code{formula}. See the specific methods help pages for information on whether they accept \code{moments}.
}
  \item{int}{
\code{logical}; for some methods, whether first-order interactions of the covariates are to be balanced (essentially balancing the covariances between covariates). This argument is ignored for other methods; to balance interactions between the variables, appropriate functions must be entered in \code{formula}. See the specific methods help pages for information on whether they accept \code{int}.
}
  \item{verbose}{
whether to print additional information output by the fitting function.
}
  \item{include.obj}{
whether to include in the output a list of the fit objects created in the process of estimating the weights at each time point. For example, with \code{method = "ps"}, a list of the \code{glm} objects containing the propensity score models at each time point will be included. See the help pages for each method for information on what object will be included if \code{TRUE}.
}
  \item{is.MSM.method}{
whether the method estimates weights for multiple time points all at once (\code{TRUE}) or by estimating weights at each time point and then multiplying them together (\code{FALSE}). This is only relevant for \code{method = "optweight")}, which estimates weights for longitudinal treatments all at once, and for user-specified functions.
}
  \item{weightit.force}{
several methods are not valid for estimating weights with longitudinal treatments, and will produce an error message if attempted. Set to \code{TRUE} to bypass this error message.
}
  \item{...}{
other arguments for functions called by \code{weightit} that control aspects of fitting that are not covered by the above arguments. See Details at \code{\link{weightit}}.
}
  \item{x}{
a \code{weightitMSM} object; the output of a call to \code{weightitMSM}.
}

}
\value{
  A \code{weightitMSM} object with the following elements:
  \item{weights}{The estimated weights, one for each unit.}
  \item{treat.list}{A list of the values of the time-varying treatment variables.}
  \item{covs.list}{A list of the covariates used in the fitting at each time point. Only includes the raw covariates, which may have been altered in the fitting process.}
  \item{data}{The data.frame originally entered to \code{weightitMSM}.}
  \item{estimand}{"ATE", currently the only estimand for MSMs with binary or multinomial treatments.}
  \item{method}{The weight estimation method specified.}
  \item{ps.list}{A list of the estimated propensity scores (if any) at each time point.}
  \item{s.weights}{The provided sampling weights.}
  \item{by}{A data.frame containing the \code{by} variable when specified.}
  \item{stabilization}{The stabilization factors, if any.}
}
\details{
In general, \code{weightitMSM} works by separating the estimation of weights into separate procedures for each time period based on the formulas provided. For each formula, \code{weightitMSM} simply calls \code{weightit} to that formula, collects the weights for each time period, and multiplies them together to arrive at longitudinal balancing weights.

Each formula should contain all the covariates to be balanced on. For example, the formula corresponding to the second time period should contain all the baseline covariates, the treatment variable at the first time period, and the time-varying covariates that took on values after the first treatment and before the second. Currently, only wide data sets are supported, where each unit is represented by exactly one row that contains the covariate and treatment history encoded in separate variables.

The \code{"gbm"} method, which calls \code{ps} in \pkg{twang}, yields the same results to a call to \code{iptw} in \pkg{twang}. However, the \code{cbps} method, which calls \code{CBPS} in \pkg{CBPS}, will yield different results from \code{CBMSM} in \pkg{CBPS} because \code{CBMSM} takes a different approach to generating weights than simply estimating several time-specific models.

}

\author{
Noah Greifer
}
\seealso{
\code{\link{weightit}} for information on the allowable methods
}
\references{
Cole, S. R., & Hernán, M. A. (2008). Constructing Inverse Probability Weights for Marginal Structural Models. American Journal of Epidemiology, 168(6), 656–664. \doi{10.1093/aje/kwn164}
}

\examples{
library("twang")
# library("cobalt")

data("iptwExWide", package = "twang")
(W <- weightitMSM(list(tx1 ~ age + gender + use0,
                      tx2 ~ tx1 + use1 + age + gender + use0,
                      tx3 ~ tx2 + use2 + tx1 + use1 + age + gender + use0),
                 data = iptwExWide,
                 method = "ps"))
summary(W)
# bal.tab(W)

##Going from long to wide data
data("iptwExLong", package = "twang")
wide_data <- reshape(iptwExLong$covariates,    #long data
                     timevar = "time",         #time variable
                     v.names = c("use", "tx"), #time-varying
                     idvar = "ID",             #time-stable
                     direction = "wide",
                     sep = "")

(W2 <- weightitMSM(list(tx1 ~ age + gender + use1,
                      tx2 ~ tx1 + use2 + age + gender + use1,
                      tx3 ~ tx2 + use3 + tx1 + use2 + age +
                          gender + use1),
                 data = wide_data,
                 method = "ps"))
summary(W2)

all.equal(W$weights, W2$weights)
}
