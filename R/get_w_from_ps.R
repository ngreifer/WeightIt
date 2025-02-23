#' Compute weights from propensity scores
#'
#' @description Given a vector or matrix of propensity scores, outputs a vector
#' of weights that target the provided estimand.
#'
#' @param ps a vector, matrix, or data frame of propensity scores. See Details.
#' @param treat a vector of treatment status for each individual. See Details.
#' @param estimand the desired estimand that the weights should target. Current
#'   options include `"ATE"` (average treatment effect), `"ATT"` (average
#'   treatment effect on the treated), `"ATC"` (average treatment effect on the
#'   control), `"ATO"` (average treatment effect in the overlap), `"ATM"`
#'   (average treatment effect in the matched sample), and `"ATOS"` (average
#'   treatment effect in the optimal subset). See Details.
#' @param focal when `estimand` is `"ATT"` or `"ATC"`, which group should be
#'   consider the (focal) "treated" or "control" group, respectively. If not
#'   `NULL` and `estimand` is not `"ATT"` or `"ATC"`, `estimand` will
#'   automatically be set to `"ATT"`.
#' @param treated when treatment is binary, the value of `treat` that is
#'   considered the "treated" group (i.e., the group for which the propensity
#'   scores are the probability of being in). If `NULL`, `get_w_from_ps()` will
#'   attempt to figure it out on its own using some heuristics. This really only
#'   matters when `treat` has values other than 0 and 1 and when `ps` is given
#'   as a vector or an unnamed single-column matrix or data frame.
#' @param subclass `numeric`; the number of subclasses to use when computing
#'   weights using marginal mean weighting through stratification (also known as
#'   fine stratification). If `NULL`, standard inverse probability weights (and
#'   their extensions) will be computed; if a number greater than 1, subclasses
#'   will be formed and weights will be computed based on subclass membership.
#'   `estimand` must be `"ATE"`, `"ATT"`, or `"ATC"` if `subclass` is
#'   non-`NULL`. See Details.
#' @param stabilize `logical`; whether to compute stabilized weights or not.
#'   This simply involves multiplying each unit's weight by the proportion of
#'   units in their treatment group. For saturated outcome models and in balance
#'   checking, this won't make a difference; otherwise, this can improve
#'   performance.
#'
#' @returns A vector of weights. When `subclass` is not `NULL`, the subclasses
#' are returned as the `"subclass"` attribute. When `estimand = "ATOS"`, the
#' chosen value of `alpha` (the smallest propensity score allowed to remain in
#' the sample) is returned in the `"alpha"` attribute.
#'
#' @details `get_w_from_ps()` applies the formula for computing weights from
#' propensity scores for the desired estimand. The formula for each estimand is
#' below, with \eqn{A_i} the treatment value for unit \eqn{i} taking on values
#' \eqn{\mathcal{A} = (1, \ldots, g)}, \eqn{p_{a, i}} the probability of
#' receiving treatment level \eqn{a} for unit \eqn{i}, and \eqn{f} is the focal
#' group (the treated group for the ATT and the control group for the ATC):
#'
#' \deqn{
#' \begin{aligned}
#' w^{ATE}_i &= 1 / p_{A_i, i} \\
#' w^{ATT}_i &= w^{ATE}_i \times p_{f, i} \\
#' w^{ATO}_i &= w^{ATE}_i / \sum_{a \in \mathcal{A}}{1/p_{a, i}} \\
#' w^{ATM}_i &= w^{ATE}_i \times \min(p_{1, i}, \ldots, p_{g, i}) \\
#' w^{ATOS}_i &= w^{ATE}_i \times \mathbb{1}\left(\alpha < p_{2, i} < 1 - \alpha\right)
#' \end{aligned}
#' }
#'
#' `get_w_from_ps()` can only be used with binary and multi-category treatments.
#'
#' ## Supplying the `ps` argument
#'
#' The `ps` argument can be entered in two ways:
#' * A numeric matrix with a row for each unit and a (named) column for each treatment level, with each cell corresponding to the probability of receiving the corresponding treatment level
#' * A numeric vector with a value for each unit corresponding to the probability of being "treated" (only allowed for binary treatments)
#'
#' When supplied as a vector, `get_w_from_ps()` has to know which value of
#' `treat` corresponds to the "treated" group. For 0/1 variables, 1 will be
#' considered treated. For other types of variables, `get_w_from_ps()` will try
#' to figure it out using heuristics, but it's safer to supply an argument to
#' `treated`. When `estimand` is `"ATT"` or `"ATC"`, supplying a value to
#' `focal` is sufficient (for ATT, `focal` is the treated group, and for ATC,
#' `focal` is the control group).
#'
#' When supplied as a matrix, the columns must be named with the levels of the
#' treatment, and it is assumed that each column corresponds to the probability
#' of being in that treatment group. This is the safest way to supply `ps`
#' unless `treat` is a 0/1 variable. When `estimand` is `"ATT"` or `"ATC"`, a
#' value for `focal` must be specified.
#'
#' ## Marginal mean weighting through stratification (MMWS)
#'
#' When `subclass` is not `NULL`, MMWS weights are computed. The implementation
#' differs slightly from that described in Hong (2010, 2012). First, subclasses
#' are formed by finding the quantiles of the propensity scores in the target
#' group (for the ATE, all units; for the ATT or ATC, just the units in the
#' focal group). Any subclasses lacking members of a treatment group will be
#' filled in with them from neighboring subclasses so each subclass will always
#' have at least one member of each treatment group. A new subclass-propensity
#' score matrix is formed, where each unit's subclass-propensity score for each
#' treatment value is computed as the proportion of units with that treatment
#' value in the unit's subclass. For example, if a subclass had 10 treated units
#' and 90 control units in it, the subclass-propensity score for being treated
#' would be .1 and the subclass-propensity score for being control would be .9
#' for all units in the subclass.
#'
#' For multi-category treatments, the propensity scores for each treatment are
#' stratified separately as described in Hong (2012); for binary treatments,
#' only one set of propensity scores are stratified and the subclass-propensity
#' scores for the other treatment are computed as the complement of the
#' propensity scores for the stratified treatment.
#'
#' After the subclass-propensity scores have been computed, the standard
#' propensity score weighting formulas are used to compute the unstabilized MMWS
#' weights. To estimate MMWS weights equivalent to those described in Hong
#' (2010, 2012), `stabilize` must be set to `TRUE`, but, as with standard
#' propensity score weights, this is optional. Note that MMWS weights are also
#' known as fine stratification weights and described by Desai et al. (2017).
#'
#' @seealso [`method_glm`]
#'
#' @references ## Binary treatments
#'
#' - `estimand = "ATO"`
#'
#' Li, F., Morgan, K. L., & Zaslavsky, A. M. (2018). Balancing covariates via
#' propensity score weighting. Journal of the American Statistical Association,
#' 113(521), 390–400. \doi{10.1080/01621459.2016.1260466}
#'
#' - `estimand = "ATM"`
#'
#' Li, L., & Greene, T. (2013). A Weighting Analogue to Pair Matching in
#' Propensity Score Analysis. The International Journal of Biostatistics, 9(2).
#' \doi{10.1515/ijb-2012-0030}
#'
#' - `estimand = "ATOS"`
#'
#' Crump, R. K., Hotz, V. J., Imbens, G. W., & Mitnik, O. A. (2009). Dealing
#' with limited overlap in estimation of average treatment effects. Biometrika,
#' 96(1), 187–199. \doi{10.1093/biomet/asn055}
#'
#' - Other estimands
#'
#' Austin, P. C. (2011). An Introduction to Propensity Score Methods for
#' Reducing the Effects of Confounding in Observational Studies. Multivariate
#' Behavioral Research, 46(3), 399–424. \doi{10.1080/00273171.2011.568786}
#'
#' - Marginal mean weighting through stratification (MMWS)
#'
#' Hong, G. (2010). Marginal mean weighting through stratification: Adjustment
#' for selection bias in multilevel data. Journal of Educational and Behavioral
#' Statistics, 35(5), 499–531. \doi{10.3102/1076998609359785}
#'
#' Desai, R. J., Rothman, K. J., Bateman, B. . T., Hernandez-Diaz, S., &
#' Huybrechts, K. F. (2017). A Propensity-score-based Fine Stratification
#' Approach for Confounding Adjustment When Exposure Is Infrequent:
#' Epidemiology, 28(2), 249–257. \doi{10.1097/EDE.0000000000000595}
#'
#' ## Multi-Category Treatments
#'
#' - `estimand = "ATO"`
#'
#' Li, F., & Li, F. (2019). Propensity score weighting for causal inference with
#' multiple treatments. The Annals of Applied Statistics, 13(4), 2389–2415.
#' \doi{10.1214/19-AOAS1282}
#'
#' - `estimand = "ATM"`
#'
#' Yoshida, K., Hernández-Díaz, S., Solomon, D. H., Jackson, J. W., Gagne, J.
#' J., Glynn, R. J., & Franklin, J. M. (2017). Matching weights to
#' simultaneously compare three treatment groups: Comparison to three-way
#' matching. Epidemiology (Cambridge, Mass.), 28(3), 387–395.
#' \doi{10.1097/EDE.0000000000000627}
#'
#' - Other estimands
#'
#' McCaffrey, D. F., Griffin, B. A., Almirall, D., Slaughter, M. E., Ramchand,
#' R., & Burgette, L. F. (2013). A Tutorial on Propensity Score Estimation for
#' Multiple Treatments Using Generalized Boosted Models. Statistics in Medicine,
#' 32(19), 3388–3414. \doi{10.1002/sim.5753}
#'
#' - Marginal mean weighting through stratification
#'
#' Hong, G. (2012). Marginal mean weighting through stratification: A
#' generalized method for evaluating multivalued and multiple treatments with
#' nonexperimental data. *Psychological Methods*, 17(1), 44–60.
#' \doi{10.1037/a0024918}
#'
#' @examples
#'
#' library("cobalt")
#' data("lalonde", package = "cobalt")
#'
#' ps.fit <- glm(treat ~ age + educ + race + married +
#'                 nodegree + re74 + re75, data = lalonde,
#'               family = binomial)
#' ps <- ps.fit$fitted.values
#'
#' w1 <- get_w_from_ps(ps, treat = lalonde$treat,
#'                     estimand = "ATT")
#'
#' treatAB <- factor(ifelse(lalonde$treat == 1, "A", "B"))
#' w2 <- get_w_from_ps(ps, treat = treatAB,
#'                     estimand = "ATT", focal = "A")
#' all.equal(w1, w2)
#' w3 <- get_w_from_ps(ps, treat = treatAB,
#'                     estimand = "ATT", treated = "A")
#' all.equal(w1, w3)
#'
#' # Using MMWS
#' w4 <- get_w_from_ps(ps, treat = lalonde$treat,
#'                     estimand = "ATE", subclass = 20,
#'                     stabilize = TRUE)
#'
#' # A multi-category example using predicted probabilities
#' # from multinomial logistic regression
#' T3 <- factor(sample(c("A", "B", "C"), nrow(lalonde),
#'                     replace = TRUE))
#'
#' multi.fit <- multinom_weightit(
#'   T3 ~ age + educ + race + married +
#'     nodegree + re74 + re75, data = lalonde,
#'   vcov = "none"
#' )
#'
#' ps.multi <- fitted(multi.fit)
#' head(ps.multi)
#'
#' w5 <- get_w_from_ps(ps.multi, treat = T3,
#'                     estimand = "ATE")

#' @export
get_w_from_ps <- function(ps, treat, estimand = "ATE", focal = NULL, treated = NULL,
                          subclass = NULL, stabilize = FALSE) {

  if (!has_treat_type(treat)) treat <- assign_treat_type(treat)
  treat.type <- get_treat_type(treat)

  if (treat.type == "continuous") {
    .err("`get_w_from_ps()` can only be used with binary or multi-category treatments")
  }

  estimand <- .process_estimand(estimand, method = "glm", treat.type = treat.type)

  processed.estimand <- .process_focal_and_estimand(focal, estimand, treat, treated)
  estimand <- toupper(processed.estimand$estimand)
  focal <- processed.estimand$focal
  assumed.treated <- processed.estimand$treated

  ps_mat <- .ps_to_ps_mat(ps, treat, assumed.treated, treat.type, treated, estimand)

  if (nrow(ps_mat) != length(treat)) {
    .err("`ps` and `treat` must have the same number of units")
  }

  if (is_not_null(subclass)) {
    #Get MMW subclass propensity scores
    ps_mat <- .subclass_ps_multi(ps_mat, treat, estimand, focal, subclass)
  }

  n <- length(treat)
  w <- rep.int(1, n)

  if (is.factor(treat)) {
    treat <- as.integer(factor(treat, levels = colnames(ps_mat)))
  }
  else {
    treat <- match(as.character(treat), colnames(ps_mat))
  }

  focal <- match(as.character(focal), colnames(ps_mat))

  if (estimand == "ATE") {
    w[] <- 1 / ps_mat[cbind(seq_len(n), treat)]
  }
  else if (estimand %in% c("ATT", "ATC")) {
    not_focal <- which(treat != focal)
    w[not_focal] <- ps_mat[not_focal, focal] / ps_mat[cbind(not_focal, treat[not_focal])]
  }
  else if (estimand == "ATO") {
    if (treat.type == "binary") {
      w[] <- ps_mat[cbind(seq_len(n), 3L - treat)]
    }
    else {
      #Li & Li (2019)
      w[] <- 1 / (ps_mat[cbind(seq_len(n), treat)] * rowSums(1 / ps_mat))
    }
  }
  else if (estimand == "ATM") {
    min_ind <- max.col(-ps_mat, ties.method = "first")
    no_match <- which(ps_mat[cbind(seq_len(n), treat)] != ps_mat[cbind(seq_len(n), min_ind)])

    if (length(no_match) > 0L) {
      w[no_match] <- ps_mat[cbind(no_match, min_ind[no_match])] /
        ps_mat[cbind(no_match, treat[no_match])]
    }
  }
  else if (estimand == "ATOS") {
    #Crump et al. (2009)
    ps.sorted <- sort(ps_mat)
    z <- ps_mat[, 1L] * ps_mat[, 2L]
    alpha.opt <- 0
    for (i in seq_len(sum(ps_mat[, 2L] < .5))) {
      if (i == 1L || !check_if_zero(ps.sorted[i] - ps.sorted[i - 1L])) {
        alpha <- ps.sorted[i]
        a <- alpha * (1 - alpha)
        if (2 * a * sum(1 / z[z >= a]) / sum(z >= a) >= 1) {
          alpha.opt <- alpha
          break
        }
      }
    }

    w[] <- 1 / ps_mat[cbind(seq_len(n), treat)]
    w[!between(ps_mat[, 2L], c(alpha.opt, 1 - alpha.opt))] <- 0
  }

  if (stabilize) w <- stabilize_w(w, treat)

  names(w) <- if_null_then(rownames(ps_mat), names(treat), NULL)

  attr(w, "subclass") <- attr(ps_mat, "sub_mat")
  if (estimand == "ATOS") attr(w, "alpha") <- alpha.opt

  w
}

.ps_to_ps_mat <- function(ps, treat, assumed.treated = NULL, treat.type = NULL,
                          treated = NULL, estimand = NULL) {

  t.levels <- {
    if (is.factor(treat)) levels(treat)
    else unique(treat, nmax = switch(treat.type, binary = 2L, length(treat) / 4))
  }

  if (treat.type == "binary") {
    if (is.matrix(ps)) {
      if (!is.numeric(ps)) {
        .err("`ps` must be numeric when supplied as a matrix")
      }
      ps.names <- rownames(ps)
    }
    else if (is.data.frame(ps)) {
      if (!all_apply(ps, is.numeric)) {
        .err("all columns of `ps` must be numeric when supplied as a data.frame")
      }
      ps.names <- rownames(ps)
      ps <- as.matrix(ps)
    }
    else if (is.numeric(ps) && is_null(dim(ps))) {
      ps.names <- names(ps)
      ps <- matrix(ps, ncol = 1L)
    }
    else {
      .err("`ps` must be a matrix, data frame, or vector of propensity scores")
    }

    if (ncol(ps) == 1L) {
      if (is_not_null(treated)) {
        if (treated %nin% t.levels) {
          .err("the argument to `treated` must be a value in `treat`")
        }
        treated.level <- treated
      }
      else if (is_not_null(assumed.treated)) {
        treated.level <- assumed.treated
      }
      else if (can_str2num(treat) &&
               all(check_if_zero(binarize(treat) - str2num(treat)))) {
        treated.level <- 1
      }
      else if (is_not_null(colnames(ps)) && colnames(ps) %in% as.character(t.levels)) {
        treated.level <- colnames(ps)
      }
      else {
        .err("if the treatment has two non-0/1 levels and `ps` is a vector or has only one column, an argument to `treated` must be supplied")
      }

      t.levels <- c(setdiff(t.levels, treated.level), treated.level)
      ps <- matrix(c(1 - ps[, 1L], ps[, 1L]), ncol = 2L,
                   dimnames = list(ps.names, as.character(t.levels)))
    }
    else if (ncol(ps) == 2L) {
      if (!all(as.character(t.levels) %in% colnames(ps))) {
        .err("if `ps` has two columns, they must be named with the treatment levels")
      }
    }
    else {
      .err("`ps` cannot have more than two columns if the treatment is binary")
    }

  }
  else if (treat.type == "multi-category") {
    if (is.matrix(ps)) {
      if (!is.numeric(ps)) {
        .err("`ps` must be numeric when supplied as a matrix")
      }
      ps.names <- rownames(ps)
    }
    else if (is.data.frame(ps)) {
      if (!all_apply(ps, is.numeric)) {
        .err("all columns of `ps` must be numeric when supplied as a data.frame")
      }
      ps.names <- rownames(ps)
      ps <- as.matrix(ps)
    }
    else {
      .err("`ps` must be a matrix or data frame of propensity scores")
    }

    if (ncol(ps) != nunique(treat)) {
      .err("`ps` must have as many columns as there are treatment levels")
    }

    if (!all(t.levels %in% colnames(ps))) {
      .err("the columns of `ps` must be named with the treatment levels")
    }

  }

  ps
}
