\name{method_user}
\alias{method_user}

\title{User-Defined Functions for Estimating Weights}

\description{
This page explains the details of estimating weights using a user-defined function. The function must take in arguments that are passed to it by \code{weightit} or \code{weightitMSM} and return a vector of weights or a list containing the weights.

To supply  user-defined function, the function object should be entered directly to \code{method}; for example, for a function \code{fun}, \code{method = fun}.

\subsection{Point Treatments}{
The following arguments are automatically passed to the user-defined function, which should have named parameters corresponding to them:
\itemize{
\item{\code{treat}: a vector of treatment status for each unit. This comes directly from the left hand side of the formula passed to \code{weightit} and so will have it's type (e.g., numeric, factor, etc.), which may need to be converted.}
\item{\code{covs}: a data frame of covariate values for each unit. This comes directly from the right hand side of the formula passed to \code{weightit}. The covariates are processed so that all columns are numeric; all factor variables are split into dummies and all interactions are evaluated. All levels of factor variables are given dummies, so the matrix of the covariates is not full rank. Users can use \code{\link{make_full_rank}}, which accepts a numeric matrix or data frame and removes columns to make it full rank, if a full rank covariate matrix is desired.}
\item{\code{s.weights}: a numeric vector of sampling weights, one for each unit.}
\item{\code{ps}: a numeric vector of propensity scores.}
\item{\code{subset}: a logical vector the same length as \code{treat} that is \code{TRUE} for units to be included in the estimation and \code{FALSE} otherwise. This is used to subset the input objects when \code{exact} is used. \code{treat}, \code{covs}, \code{s.weights}, and \code{ps}, if supplied, will already have been subsetted by \code{subset}.}
\item{\code{estimand}: a character vector of length 1 containing the desired estimand. The characters will have been converted to uppercase. If "ATC" was supplied to estimand, \code{weightit} sets \code{focal} to the control level (usually 0 or the lowest level of \code{treat}) and sets \code{estimand} to "ATT".}
\item{\code{focal}: a character vector of length 1 containing the focal level of the treatment when the estimand is the ATT (or the ATC as detailed above). \code{weightit} ensures the value of focal is a level of \code{treat}.}
\item{\code{stabilize}: a logical vector of length 1. It is not processed by \code{weightit} before it reaches the fitting function.}
\item{\code{moments}: a numeric vector of length 1. It is not processed by \code{weightit} before it reaches the fitting function except that \code{as.integer} is applied to it. This is used in other methods to determine whether polynomials of the entered covariates are to be used in the weight estimation.}
\item{\code{int}: a logical vector of length 1. It is not processed by \code{weightit} before it reaches the fitting function. This is used in other methods to determine whether interactions of the entered covariates are to be used in the weight estimation.}
}
None of these parameters are required to be in the fitting function. These are simply those that are automatically available.

In addition, any additional arguments supplied to \code{weightit} will be passed on to the fitting function. \code{weightit} ensures the arguments correspond to the parameters of the fitting function and throws an error if an incorrectly named argument is supplied and the fitting function doesn't include \code{\dots} as a parameter.

The fitting function must output either a numeric vector of weights or a list (or list-like object) with an entry named wither "w" or "weights". If a list, the list can contain other named entries, but only entries named "w", "weights", "ps", and "fit.obj" will be processed. "ps" is a vector of propensity scores and "fit.obj" should be an object used in the fitting process that a user may want to examine and that is included in the \code{weightit} output object as "obj" when \code{include.obj = TRUE}. The "ps" and "fit.obj" components are optional, but "weights" or "w" is required.
}
\subsection{Longitudinal Treatments}{
Longitudinal treatments can be handled either by running the fitting function for point treatments for each time point and multiplying the resulting weights together or by running a method that accommodates multiple time points and outputs a single set of weights. For the former, \code{weightitMSM} can be used with the user-defined function just as it is with \code{weightit}. The latter method is not yet accommodated by \code{weightitMSM}, but will be someday, maybe.
}
}

\seealso{
\code{\link{weightit}}, \code{\link{weightitMSM}}
}

\examples{

library("cobalt")
data("lalonde", package = "cobalt")

#A user-defined version of method = "ps"
my.ps <- function(treat, covs, estimand, focal = NULL) {
  covs <- make_full_rank(covs)
  d <- data.frame(treat, covs)
  f <- formula(d)
  ps <- glm(f, data = d, family = "binomial")$fitted
  w <- get_w_from_ps(ps, treat = treat, estimand = estimand,
                     focal = focal)

  return(list(w = w, ps = ps))
}

#Balancing covariates between treatment groups (binary)
(W1 <- weightit(treat ~ age + educ + married +
                  nodegree + re74, data = lalonde,
                method = my.ps, estimand = "ATT"))
summary(W1)
bal.tab(W1)

#Balancing covariates for longitudinal treatments
# my.ps is used at each time point.
library("twang")
data("iptwExWide", package = "twang")
(W2 <- weightitMSM(list(tx1 ~ age + gender + use0,
                        tx2 ~ tx1 + use1 + age + gender + use0,
                        tx3 ~ tx2 + use2 + tx1 + use1 + age + gender + use0),
                   data = iptwExWide,
                   method = my.ps))
summary(W2)
bal.tab(W2)

# Kernel balancing using the KBAL package, available
# using devtools::install_github("chadhazlett/KBAL").
# Only the ATT and ATC are available. Use 'kbal.method'
# instead of 'method' in weightit() to choose between
# "ebal" and "el".

\dontrun{
kbal.fun <- function(treat, covs, estimand, focal, ...) {
    args <- list(...)
    if (is_not_null(focal))
        treat <- as.numeric(treat == focal)
    else if (estimand != "ATT")
        stop("estimand must be 'ATT' or 'ATC'.", call. = FALSE)
    if ("kbal.method" \%in\% names(args)) {
        names(args)[names(args) == "kbal.method"] <- "method"
    }
    args[names(args) \%nin\% setdiff(names(formals(KBAL::kbal)),
        c("X", "D"))] <- NULL
    k.out <- do.call(KBAL::kbal, c(list(X = covs, D = treat),
        args))
    w <- k.out$w
    return(list(w = w))
}

(Wk <- weightit(treat ~ age + educ + married +
                nodegree + re74, data = lalonde,
                method = kbal.fun, estimand = "ATT",
                kbal.method = "ebal"))
summary(Wk)
bal.tab(Wk, disp.ks = TRUE)
}
}