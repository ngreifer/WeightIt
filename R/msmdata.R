#' Simulated data for a 3 time point sequential study
#'
#' @description
#' This is a simulated dataset of 7500 units with covariates and treatment
#' measured three times and the outcome measured at the end from a
#' hypothetical observational study examining the effect of treatment
#' delivered at each time point on an adverse event.
#'
#' The data were generated using a simple simulation mechanism.
#' For further details on how the dataset was built, see the code
#' at \href{https://github.com/ngreifer/Weightit/blob/master/data-raw/msmdata.R}{data-raw/msmdata.R}.
#'
#' The dataset is provided to illustrate the features of
#' `weightitMSM()` and is not based on a realistic data-generating
#' process, so it should not be used as a benchmark.
#'
#' For simulating realistic data with a
#' known data-generating mechanism, consider using the \pkg{simcausal} package.
#'
#' @name msmdata
#' @docType data
#' @format A data frame with 7500 observations on the following 10 variables.
#' \describe{
#' \item{`X1_0`}{a count covariate measured at baseline}
#' \item{`X2_0`}{a binary covariate measured at baseline}
#' \item{`A_1`}{a binary indicator of treatment status at the first time point}
#' \item{`X1_1`}{a count covariate measured at the first time point (after the first treatment)}
#' \item{`X2_1`}{a binary covariate measured at the first time point (after the first treatment)}
#' \item{`A_2`}{a binary indicator of treatment status at the second time point}
#' \item{`X1_2`}{a count covariate measured at the second time point (after the second treatment)}
#' \item{`X2_2`}{a binary covariate measured at the first time point (after the first treatment)}
#' \item{`A_3`}{a binary indicator of treatment status at the third time point}
#' \item{`Y_B`}{a binary indicator of the outcome event (e.g., death)}
#' }
#' @keywords datasets
#' @examples
#'
#' data("msmdata")
#'
"msmdata"
