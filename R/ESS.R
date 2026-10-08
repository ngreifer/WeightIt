#' Compute effective sample size of weighted sample
#'
#' @description
#' Computes the effective sample size (ESS) of a weighted sample,
#' which represents the size of an unweighted sample with approximately the same
#' amount of precision as the weighted sample under consideration.
#'
#' @param w a vector of weights.
#'
#' @details
#' The ESS is calculated as \eqn{(\sum w)^2/\sum w^2}. It is invariant to multiplicative scaling of the weights (i.e., multiplying all weights by a nonzero scalar).
#'
#' @returns
#' A single number, the effective sample size. For non-negative weights it lies between 1 and `length(w)`, and equals `length(w)` only when all the weights are equal. It is `NA` if any weight is missing.
#'
#' @seealso [summary.weightit()]
#'
#' `vignette("weighting-methods")` for the role of the effective sample size in choosing a weighting specification.
#'
#' @references
#' McCaffrey, D. F., Ridgeway, G., & Morral, A. R. (2004).
#' Propensity Score Estimation With Boosted Regression for Evaluating Causal
#' Effects in Observational Studies. *Psychological Methods*, 9(4), 403–425.
#' \doi{10.1037/1082-989X.9.4.403}
#'
#' Shook-Sa, B. E., & Hudgens, M. G. (2020). Power and sample size for
#' observational studies of point exposure effects. *Biometrics*, biom.13405.
#' \doi{10.1111/biom.13405}
#'
#' @examples
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' #Balancing covariates between treatment groups (binary)
#' (W1 <- weightit(treat ~ age + educ + married +
#'                   nodegree + re74, data = lalonde,
#'                 method = "glm", estimand = "ATE"))
#'
#' summary(W1)
#'
#' ESS(W1$weights[W1$treat == 0])
#' ESS(W1$weights[W1$treat == 1])

#' @export
ESS <- function(w) {
  arg::arg_supplied(w)
  arg::arg_numeric(w)

  sum(w)^2 / sum(w^2)
}
