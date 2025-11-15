#' Weighting methods
#'
#' @description
#' `.weightit_methods` is a list containing the allowable weighting methods that can be supplied by name to the `method` argument of [weightit()], [weightitMSM()], and [weightit.fit()]. Each entry corresponds to an allowed method and contains information about what options are and are not allowed for each method. While this list is primarily for internal use by checking functions in \pkg{WeightIt}, it might be of use for package authors that want to support different weighting methods.
#'
#' @details
#' Each component is itself a list containing the following components:
#' \describe{
#' \item{`treat_type`}{at least one of `"binary"`, `"multinomial"`, or `"continuous"` indicating which treatment types are available for this method.}
#' \item{`estimand`}{which estimands are available for this method. All methods that support binary and multi-category treatments accept `"ATE"`, `"ATT"`, and `"ATC"`, as well as some other estimands depending on the method. See [get_w_from_ps()] for more details about what each estimand means.}
#' \item{`alias`}{a character vector of aliases for the method. When an alias is supplied, the corresponding method will still be dispatched. For example, the canonical method to request entropy balancing is `"ebal"`, but `"ebalance"` and `"entropy"` also work. The first value is the canonical name.}
#' \item{`description`}{a string containing the description of the name in English.}
#' \item{`ps`}{a logical for whether propensity scores are returned by the method for binary treatments. Propensity scores are never returned for multi-category or continuous treatments.}
#' \item{`msm_valid`}{a logical for whether the method can be validly used with longitudinal treatments.}
#' \item{`msm_method_available`}{a logical for whether a version of the method can be used that estimates weights using a single model rather than multiplying the weights across time points. This is related to the `is.MSM.method` argument of `weightitMSM()`.}
#' \item{`subclass_ok`}{a logical for whether `subclass` can be supplied to compute subclassification weights from the propensity scores.}
#' \item{`packages_needed`}{a character vector of the minimal packages required to use the method. Some methods may require additional packages for certain options.}
#' \item{`package_versions_needed`}{a named character vector of the minimal package versions required to use the method. Only present when `packages_needed` is not empty.}
#' \item{`s.weights_ok`}{a logical for whether sampling weights can be used with the method.}
#' \item{`missing`}{a character vector of the allowed options that can be supplied to `missing` when missing data is present. All methods accept `"ind"` for the missingness indicator approach; some other methods accept additional values.}
#' \item{`moments_int_ok`}{a logical for whether `moments`, `int`, and `quantile` can be used with the method.}
#' \item{`moments_default`}{when `moments_int_ok` is `TRUE`, the default value of `moments` used with the method. For most methods, this is 1.}
#' \item{`density_ok`}{a logical for whether arguments that control the density can be used with the method when used with a continuous treatment.}
#' \item{`stabilize_ok`}{a logical for whether the `stabilize` argument (and `num.formula` for longitudinal treatments) can be used with the method.}
#' \item{`plot.weightit_ok`}{a logical for whether `plot()` can be used on the `weightit` output with the method.}
#' }
#'
#' @seealso
#' [weightit()] and [weightitMSM()] for how the methods are used. Also see the individual methods pages for information on whether and how each option can be used.
#'
#' @examples
#' # Get all acceptable names
#' names(.weightit_methods)
#'
#' # Get all acceptable names and aliases
#' lapply(.weightit_methods, `[[`, "alias")
#'
#' # Which estimands are allowed with `method = "bart"`
#' .weightit_methods[["bart"]]$estimand
#'
#' # All methods that support continuous treatments
#' supp <- sapply(.weightit_methods, function(x) {
#'   "continuous" %in% x$treat_type
#' })
#' names(.weightit_methods)[supp]
#'
#' # All methods that return propensity scores (for
#' # binary treatments only)
#' supp <- sapply(.weightit_methods, `[[`, "ps")
#' names(.weightit_methods)[supp]

#' @export
.weightit_methods <- {list(
  "glm" = list(
    treat_type = c("binary", "multinomial", "continuous"),
    estimand = c("ATE", "ATT", "ATC", "ATO", "ATM", "ATOS"),
    alias = c("glm", "ps"),
    description = "propensity score weighting with GLM",
    ps = TRUE,
    msm_valid = TRUE,
    msm_method_available = FALSE,
    subclass_ok = TRUE,
    packages_needed = character(),
    s.weights_ok = TRUE,
    missing = c("ind", "saem"),
    moments_int_ok = FALSE,
    moments_default = NULL,
    density_ok = TRUE,
    stabilize_ok = TRUE,
    plot.weightit_ok = FALSE
  ),
  "bart" = list(
    treat_type = c("binary", "multinomial", "continuous"),
    estimand = c("ATE", "ATT", "ATC", "ATO", "ATM", "ATOS"),
    alias = "bart",
    description = "propensity score weighting with BART",
    ps = TRUE,
    msm_valid = TRUE,
    msm_method_available = FALSE,
    subclass_ok = TRUE,
    packages_needed = "dbarts",
    package_versions_needed = c(dbarts = "0.9-29"),
    s.weights_ok = FALSE,
    missing = "ind",
    moments_int_ok = FALSE,
    moments_default = NULL,
    density_ok = TRUE,
    stabilize_ok = TRUE,
    plot.weightit_ok = FALSE
  ),
  "cbps" = list(
    treat_type = c("binary", "multinomial", "continuous"),
    estimand = c("ATE", "ATT", "ATC", "ATO"),
    alias = c("cbps", "cbgps"),
    description = "covariate balancing propensity score weighting",
    ps = TRUE,
    msm_valid = TRUE,
    msm_method_available = TRUE,
    subclass_ok = FALSE,
    packages_needed = character(),
    s.weights_ok = TRUE,
    missing = "ind",
    moments_int_ok = TRUE,
    moments_default = 1,
    density_ok = FALSE,
    stabilize_ok = FALSE,
    plot.weightit_ok = FALSE
  ),
  "ebal" = list(
    treat_type = c("binary", "multinomial", "continuous"),
    estimand = c("ATE", "ATT", "ATC"),
    alias = c("ebal", "ebalance", "entropy"),
    description = "entropy balancing",
    ps = FALSE,
    msm_valid = FALSE,
    msm_method_available = FALSE,
    subclass_ok = FALSE,
    packages_needed = character(),
    s.weights_ok = TRUE,
    missing = "ind",
    moments_int_ok = TRUE,
    moments_default = 1,
    density_ok = FALSE,
    stabilize_ok = FALSE,
    plot.weightit_ok = FALSE
  ),
  "energy" = list(
    treat_type = c("binary", "multinomial", "continuous"),
    estimand = c("ATE", "ATT", "ATC"),
    alias = c("energy", "dcows"),
    description = "energy balancing",
    ps = FALSE,
    msm_valid = FALSE,
    msm_method_available = FALSE,
    subclass_ok = FALSE,
    packages_needed = "osqp",
    s.weights_ok = TRUE,
    missing = "ind",
    moments_int_ok = TRUE,
    moments_default = 0,
    density_ok = FALSE,
    stabilize_ok = FALSE,
    plot.weightit_ok = FALSE
  ),
  "gbm" = list(
    treat_type = c("binary", "multinomial", "continuous"),
    estimand = c("ATE", "ATT", "ATC", "ATO", "ATM", "ATOS"),
    alias = c("gbm", "gbr"),
    description = "propensity score weighting with GBM",
    ps = TRUE,
    msm_valid = TRUE,
    msm_method_available = FALSE,
    subclass_ok = TRUE,
    packages_needed = "gbm",
    s.weights_ok = TRUE,
    missing = c("ind", "surr"),
    moments_int_ok = FALSE,
    moments_default = NULL,
    density_ok = TRUE,
    stabilize_ok = TRUE,
    plot.weightit_ok = TRUE
  ),
  "ipt" = list(
    treat_type = c("binary", "multinomial"),
    estimand = c("ATE", "ATT", "ATC"),
    alias = "ipt",
    description = "inverse probability tilting",
    ps = TRUE,
    msm_valid = TRUE,
    msm_method_available = FALSE,
    subclass_ok = FALSE,
    packages_needed = "rootSolve",
    s.weights_ok = TRUE,
    missing = "ind",
    moments_int_ok = TRUE,
    moments_default = 1,
    density_ok = FALSE,
    stabilize_ok = FALSE,
    plot.weightit_ok = FALSE
  ),
  "npcbps" = list(
    treat_type = c("binary", "multinomial", "continuous"),
    estimand = "ATE",
    alias = c("npcbps", "npcbgps"),
    description = "non-parametric covariate balancing propensity score weighting",
    ps = FALSE,
    msm_valid = FALSE,
    msm_method_available = FALSE,
    subclass_ok = FALSE,
    packages_needed = "CBPS",
    s.weights_ok = FALSE,
    missing = "ind",
    moments_int_ok = TRUE,
    moments_default = 1,
    density_ok = FALSE,
    stabilize_ok = FALSE,
    plot.weightit_ok = FALSE
  ),
  "optweight" = list(
    treat_type = c("binary", "multinomial", "continuous"),
    estimand = c("ATE", "ATT", "ATC"),
    alias = c("optweight", "sbw"),
    description = "stable balancing weights",
    ps = FALSE,
    msm_valid = FALSE,
    msm_method_available = FALSE,
    subclass_ok = FALSE,
    packages_needed = "optweight",
    package_versions_needed = c(optweight = "1.0.0"),
    s.weights_ok = TRUE,
    missing = "ind",
    moments_int_ok = TRUE,
    moments_default = 1,
    density_ok = FALSE,
    stabilize_ok = FALSE,
    plot.weightit_ok = TRUE
  ),
  "super" = list(
    treat_type = c("binary", "multinomial", "continuous"),
    estimand = c("ATE", "ATT", "ATC", "ATO", "ATM", "ATOS"),
    alias = c("super", "superlearner"),
    description = "propensity score weighting with SuperLearner",
    ps = TRUE,
    msm_valid = TRUE,
    msm_method_available = FALSE,
    subclass_ok = TRUE,
    packages_needed = "SuperLearner",
    s.weights_ok = TRUE,
    missing = "ind",
    moments_int_ok = FALSE,
    moments_default = NULL,
    density_ok = TRUE,
    stabilize_ok = TRUE,
    plot.weightit_ok = FALSE
  )
)}
