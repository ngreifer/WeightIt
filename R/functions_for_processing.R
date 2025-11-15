.method_to_proper_method <- function(method) {
  if (is_null(method)) {
    return(NULL)
  }

  if (!is.character(method)) {
    return(method)
  }

  method <- tolower(method)

  if (method %nin% unlist(grab(.weightit_methods, "alias"))) {
    return(method)
  }

  .allowable.methods <- unlist(lapply(names(.weightit_methods), function(m) {
    aliases <- .weightit_methods[[m]]$alias
    rep_with(m, aliases) |>
      setNames(aliases)
  }))

  unname(.allowable.methods[method])
}

.check_acceptable_method <- function(method, msm = FALSE, force = FALSE) {

  if (missing(method)) {
    method <- "glm"
  }
  else if (is_null(method)) {
    return(invisible(NULL))
  }

  if (identical(method, "twang")) {
    .err('"twang" is no longer an acceptable argument to `method`. Please use "gbm" for generalized boosted modeling')
  }

  if ((!is.character(method) && !is.function(method)) ||
      (is.character(method) && (length(method) > 1L ||
                                !utils::hasName(.weightit_methods,
                                                .method_to_proper_method(method))))) {
    .err(sprintf("`method` must be a string of length 1 containing the name of an acceptable weighting method or a function that produces weights. Allowable methods:\n%s",
                 word_list(names(.weightit_methods), and.or = FALSE, quotes = 2L)),
         tidy = FALSE)
  }

  if (msm && !force && is.character(method)) {
    m <- .method_to_proper_method(method)
    if (!.weightit_methods[[m]]$msm_valid) {
      .err(sprintf("the use of %s with longitudinal treatments has not been validated. Set `weightit.force = TRUE` to bypass this error",
                   .method_to_phrase(m)))
    }
  }
}

.check_method_treat.type <- function(method, treat.type) {
  if (is_not_null(method) && is.character(method) &&
      utils::hasName(.weightit_methods, method) &&
      (treat.type %nin% .weightit_methods[[method]]$treat_type)) {
    .err(sprintf("%s can only be used with a %s treatment",
                 .method_to_phrase(method),
                 word_list(.weightit_methods[[method]]$treat_type, and.or = "or")))
  }
}

.check_required_packages <- function(method) {
  if (!chk::vld_string(method) ||
      !utils::hasName(.weightit_methods, method)) {
    return(invisible(NULL))
  }

  pkgs <- .weightit_methods[[method]]$packages_needed

  if (is_not_null(pkgs)) {
    versions_needed <- .weightit_methods[[method]]$package_versions_needed

    versions <- NULL
    if (is_not_null(versions_needed)) {
      versions <- rep_with(NA_character_, pkgs) |>
        setNames(pkgs)
      versions[names(versions_needed)] <- versions_needed
    }

    rlang::check_installed(pkgs, version = versions)
  }

  invisible(NULL)
}

.process.s.weights <- function(s.weights, data = NULL) {
  #Process s.weights
  if (is_null(s.weights)) {
    return(NULL)
  }

  if (is.numeric(s.weights)) {
    return(s.weights)
  }

  if (!chk::vld_string(s.weights)) {
    .err("the argument to `s.weights` must be a vector or data frame of sampling weights or the (quoted) names of the variable in `data` that contains sampling weights")
  }

  if (is_null(data)) {
    .err("`s.weights` was specified as a string but there was no argument to `data`")
  }

  if (!utils::hasName(data, s.weights)) {
    .err("the name supplied to `s.weights` is not the name of a variable in `data`")
  }

  data[[s.weights]]
}

.check_method_s.weights <- function(method, s.weights) {
  if (is_not_null(method) &&
      !is.function(method) &&
      !.weightit_methods[[method]]$s.weights_ok &&
      !all_the_same(s.weights)) {
    .err(sprintf("sampling weights cannot be used with %s",
                 .method_to_phrase(method)))
  }
}

.method_to_phrase <- function(method) {

  if (is_null(method)) {
    return("no weighting")
  }

  if (is.function(method)) {
    return("a user-defined method")
  }

  method <- .method_to_proper_method(method)

  if (!utils::hasName(.weightit_methods, method)) {
    return("the chosen method of weighting")
  }

  .weightit_methods[[method]]$description
}

.process_estimand <- function(estimand, method, treat.type) {

  if (is.function(method)) {
    chk::chk_null_or(estimand, vld = chk::vld_string)
    return(toupper(estimand))
  }

  if (treat.type == "continuous") {
    if (is_not_null(estimand) && (!chk::vld_string(estimand) || !identical(toupper(estimand), "ATE"))) {
      .wrn("`estimand` is ignored for continuous treatments")
    }

    return("ATE")
  }

  chk::chk_string(estimand)

  allowable_estimands <- {
    if (is_null(method)) unique(unlist(grab(.weightit_methods, "estimand")))
    else .weightit_methods[[method]]$estimand
  }

  if (treat.type %in% c("multinomial", "multi-category")) {
    allowable_estimands <- setdiff(allowable_estimands, "ATOS")
  }

  if (toupper(estimand) %nin% allowable_estimands) {
    .err(sprintf("%s is not an allowable estimand for %s with a %s treatment. Only %s allowed",
                 add_quotes(estimand), .method_to_phrase(method), treat.type,
                 word_list(allowable_estimands, quotes = TRUE, and.or = "and", is.are = TRUE)))
  }

  toupper(estimand)
}

.check_subclass <- function(subclass = NULL, method, treat.type) {
  if (is_not_null(subclass) && is_not_null(method) && !is.function(method)) {

    subclass_ok <- .weightit_methods[[method]]$subclass_ok

    if (treat.type == "continuous" || !subclass_ok) {
      .err(sprintf("subclasses are not compatible with %s with a %s treatment",
                   .method_to_phrase(method), treat.type))
    }
  }
}

.process_moments_int_quantile <- function(moments = NULL, int = FALSE, quantile = NULL, method = NULL, ...) {
  if (is.function(method)) {
    return(list(moments = moments,
                int = int,
                quantile = quantile))
  }

  if (is_null(method) || !.weightit_methods[[method]]$moments_int_ok) {
    if (is_not_null(method)) {
      mi0 <- c(is_not_null(moments),
               is_not_null(int) && !isFALSE(int),
               is_not_null(quantile))

      if (any(mi0)) {
        .wrn(sprintf("%s not compatible with %s. Ignoring %s",
                     word_list(c("moments", "int", "quantile")[mi0], and.or = "and", is.are = TRUE, quotes = "`"),
                     .method_to_phrase(method),
                     word_list(c("moments", "int", "quantile")[mi0], and.or = "and", quotes = "`")))
      }
    }

    return(list(moments = integer(),
                int = FALSE,
                quantile = list()))
  }

  if (is_null(int)) {
    int <- FALSE
  }
  else {
    chk::chk_flag(int)
  }

  if (is_not_null(quantile)) {
    .vld_qu <- function(x) {
      is.numeric(x) && all(x >= 0) && all(x <= 1)
    }

    bad.q <- FALSE
    if (.vld_qu(quantile)) {
      if (length(quantile) == 1L || (is_not_null(names(quantile)) && all(nzchar(names(quantile))))) {
        quantile <- as.list(quantile)
      }
      else if (is_null(names(quantile)) && anyDuplicated(quantile) == 0L) {
        quantile <- list(quantile)
      }
      else {
        bad.q <- TRUE
      }
    }
    else if (is.list(quantile)) {
      if ((length(quantile) > 1L && (is_null(names(quantile)) || !all(nzchar(names(quantile))))) ||
          !all_apply(quantile, .vld_qu)) {
        bad.q <- TRUE
      }
    }
    else {
      bad.q <- TRUE
    }

    if (bad.q) {
      .err("`quantile` must be a number between 0 and 1, a named list or vector of such values, or a named list of vectors of such values")
    }
  }

  if (is_not_null(moments)) {
    chk::chk_whole_numeric(moments)

    if (length(moments) > 1L && (is_null(names(moments)) || !all(nzchar(names(moments))))) {
      .err("`moments` must be an integer or a named vector of integers")
    }

    chk::chk_gte(moments,
                 if (is_null(quantile)) .weightit_methods[[method]]$moments_default
                 else 0)

    if (int && any(moments < 1)) {
      .err("when `int = TRUE`, `moments` must be greater than or equal to 1")
    }

    moments[] <- as.integer(moments)
  }
  else {
    moments <- {
      if (int) 1L
      else .weightit_methods[[method]]$moments_default
    }
  }

  attr(moments, "moments_default") <- .weightit_methods[[method]]$moments_default

  list(moments = moments, int = int, quantile = quantile)
}

.process_MSM_method <- function(is.MSM.method, method) {
  if (is_null(method)) {
    return(FALSE)
  }

  if (is.function(method)) {
    if (isTRUE(is.MSM.method)) {
      .err("currently, only user-defined methods that work with `is.MSM.method = FALSE` are allowed")
    }

    return(FALSE)
  }

  if (.weightit_methods[[method]]$msm_method_available) {
    if (is_null(is.MSM.method)) {
      return(TRUE)
    }

    chk::chk_flag(is.MSM.method)

    if (!is.MSM.method) {
      .msg(sprintf("%s can be used with a single model when multiple time points are present. Using a seperate model for each time point. To use a single model, set `is.MSM.method` to `TRUE`",
                   .method_to_phrase(method)))
    }

    return(is.MSM.method)
  }

  if (is_not_null(is.MSM.method)) {

    chk::chk_flag(is.MSM.method)

    if (is.MSM.method) {
      .wrn(sprintf("%s cannot be used with a single model when multiple time points are present. Using a seperate model for each time point",
                   .method_to_phrase(method)))
    }
  }

  FALSE
}

.process_missing <- function(missing, method) {
  if (is_null(method)) {
    return("")
  }

  allowable.missings <- .weightit_methods[[method]]$missing

  if (is_null(missing)) {
    .wrn(sprintf("missing values are present in the covariates. See `?WeightIt::method_%s` for information on how these are handled",
                 method))
    return(allowable.missings[1L])
  }

  chk::chk_string(missing)

  if (missing %nin% allowable.missings) {
    .err(sprintf("only %s allowed for `missing` with %s",
                 word_list(allowable.missings, quotes = 2L, is.are = TRUE),
                 .method_to_phrase(method)))
  }

  missing
}

.missing_to_phrase <- function(missing) {
  switch(missing,
         ind = "missingness indicators",
         saem = "SAEM",
         surr = "surrogate splitting",
         missing)
}

.process_missing2 <- function(missing, covs) {
  if (is_null(missing) || identical(missing, "") || !anyNA(covs)) {
    return("")
  }

  missing
}

.check_user_method <- function(method) {
  #Check to make sure it accepts treat and covs
  if (all(c("covs", "treat") %in% names(formals(method)))) {
  }
  # else if (all(c("covs.list", "treat.list") %in% names(formals(method)))) {
  # }
  else {
    .err("the user-provided function to `method` must contain `covs` and `treat` as named parameters")
  }
}

.process_ps <- function(ps, data = NULL, treat = NULL) {
  if (is_null(ps)) {
    return(NULL)
  }

  if (chk::vld_string(ps)) {
    if (is_null(data)) {
      .err("`ps` was specified as a string but there was no argument to `data`")
    }

    if (!utils::hasName(data, ps)) {
      .err("the name supplied to `ps` is not the name of a variable in `data`")
    }

    ps <- data[[ps]]

    if (!is.numeric(ps)) {
      .err("the name supplied to `ps` must correspond to a numeric variable in `data`")
    }
  }
  else if (is.numeric(ps)) {
    if (length(ps) != length(treat)) {
      .err("`ps` must have the same number of units as the treatment")
    }
  }
  else {
    .err("the argument to `ps` must be a vector of propensity scores or the (quoted) name of a numeric variable in `data` that contains propensity scores")
  }

  ps
}

.process_focal_and_estimand <- function(focal, estimand, treat, treated = NULL) {
  reported.estimand <- estimand

  if (!has_treat_type(treat)) treat <- assign_treat_type(treat)
  treat.type <- get_treat_type(treat)

  if (treat.type == "continuous") {
    return(list(focal = NULL,
                estimand = "ATE",
                reported.estimand = "ATE"))
  }

  if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
    .wrn(sprintf('`estimand = %s` is not compatible with `focal`. Setting `estimand` to "ATT"',
                 add_quotes(estimand)))
    reported.estimand <- estimand <- "ATT"
  }

  if (treat.type == "binary") {
    ct <- .get_control_and_treated_levels(treat, estimand, focal, treated)

    focal <- switch(estimand,
                    ATT = ct["treated"],
                    ATC = ct["control"],
                    NULL)

    treated <- ct["treated"]
  }
  else {
    unique.vals <- {
      if (chk::vld_character_or_factor(treat))
        levels(factor(treat, nmax = ceiling(length(treat) / 4)))
      else
        sort(unique(treat, nmax = ceiling(length(treat) / 4)))
    }

    #Check focal
    if (is_not_null(focal) && (length(focal) > 1L || focal %nin% unique.vals)) {
      .err("the argument supplied to `focal` must be the name of a level of treatment")
    }

    if (estimand == "ATT") {
      if (is_null(focal)) {
        if (is_null(treated) || treated %nin% unique.vals) {
          .err('when `estimand = "ATT"` for multi-category treatments, an argument must be supplied to `focal`')
        }

        focal <- treated
      }
    }
    else if (estimand == "ATC") {
      if (is_null(focal)) {
        .err('when `estimand = "ATC"` for multi-category treatments, an argument must be supplied to `focal`')
      }

      estimand <- "ATT"
    }
  }

  list(focal = unname(focal),
       estimand = estimand,
       reported.estimand = reported.estimand,
       treated = switch(treat.type, binary = unname(treated), NULL))
}

.get_control_and_treated_levels <- function(treat, estimand, focal = NULL, treated = NULL) {

  if (is_not_null(.attr(treat, "control")) &&
      is_not_null(.attr(treat, "treated"))) {
    return(setNames(c(.attr(treat, "control"), .attr(treat, "treated")),
                    c("control", "treated")))
  }

  control <- NULL
  throw_message <- FALSE

  unique.vals <- {
    if (chk::vld_character_or_factor(treat))
      levels(factor(treat, nmax = 2L))
    else
      sort(unique(treat, nmax = 2L))
  }

  if (is_not_null(focal)) {
    if (length(focal) > 1L || focal %nin% unique.vals) {
      .err("the argument supplied to `focal` must be the name of a level of treatment")
    }

    if (estimand == "ATC") {
      control <- focal
      treated <- NULL
    }
    else {
      treated <- focal
    }
  }
  else if (is_not_null(.attr(treat, "treated"))) {
    treated <- .attr(treat, "treated")
  }
  else if (is_not_null(.attr(treat, "control"))) {
    control <- .attr(treat, "control")
  }
  else if (is_not_null(treated)) {
    if (length(treated) > 1L || treated %nin% unique.vals) {
      .err("the argument supplied to `treated` must be the name of a level of treatment")
    }
  }
  else if (is.logical(treat)) {
    treated <- TRUE
    control <- FALSE
  }
  else if (is.numeric(unique.vals)) {
    control <- unique.vals[unique.vals == 0]

    if (is_null(control)) {
      control <- unique.vals[1L]

      throw_message <- TRUE
    }
  }
  else if (can_str2num(unique.vals)) {
    unique.vals.numeric <- str2num(unique.vals)

    control <- unique.vals[unique.vals.numeric == 0]

    if (is_null(control)) {
      control <- unique.vals[which.min(unique.vals.numeric)]

      throw_message <- TRUE
    }
  }
  else {
    treated_options <- c("t", "tr", "treat", "treated", "exposed")
    control_options <- c("c", "co", "ctrl", "control", "unexposed")

    t_match <- which(unique.vals %in% treated_options)
    c_match <- which(unique.vals %in% control_options)

    if (length(t_match) == 1L) {
      treated <- unique.vals[t_match]
    }
    else if (length(c_match) == 1L) {
      control <- unique.vals[c_match]
    }
  }

  if (is_null(control) && is_null(treated)) {
    control <- unique.vals[1L]
    treated <- unique.vals[2L]

    throw_message <- TRUE
  }
  else if (is_null(control)) {
    control <- setdiff(unique.vals, treated)
  }
  else if (is_null(treated)) {
    treated <- setdiff(unique.vals, control)
  }

  if (throw_message) {
    if (estimand == "ATT") {
      .msg(sprintf("assuming %s is the treated level. If not, supply an argument to `focal`",
                   add_quotes(treated, !is.numeric(unique.vals))))

    }
    else if (estimand == "ATC") {
      .msg(sprintf("assuming %s is the control level. If not, supply an argument to `focal`",
                   add_quotes(control, !is.numeric(unique.vals))))
    }
    else {
      .msg(sprintf("assuming %s is the treated level. If not, recode the treatment so that 1 is treated and 0 is control",
                   add_quotes(treated, !is.numeric(unique.vals))))
    }
  }

  setNames(c(control, treated),
           c("control", "treated"))
}

get_treated_level <- function(treat, estimand, focal = NULL) {
  ct <- .get_control_and_treated_levels(treat, estimand, focal)

  unname(ct["treated"])
}

.process_by <- function(by, data, treat, treat.name = NULL, by.arg = "by") {

  ##Process by
  bad.by <- FALSE
  n <- length(treat)

  if (!has_treat_type(treat)) treat <- assign_treat_type(treat)
  treat.type <- get_treat_type(treat)

  if (missing(by)) {
    bad.by <- TRUE
  }
  else if (is_null(by)) {
    by <- NULL
    by.name <- NULL
  }
  else if (chk::vld_string(by) && utils::hasName(data, by)) {
    by.name <- by
    by <- data[[by]]
  }
  else if (length(dim(by)) == 2L && len(by) == n) {
    by.name <- colnames(by)[1L]
    by <- drop(by[, 1L])
  }
  else if (rlang::is_formula(by, lhs = FALSE)) {
    t.c <- get_covs_and_treat_from_formula2(by, data)
    by <- t.c[["reported.covs"]]
    if (NCOL(by) != 1L) {
      .err(sprintf("only one variable can be on the right-hand side of the formula for `%s`",
                   by.arg))
    }
    by.name <- colnames(by)
  }
  else {
    bad.by <- TRUE
  }

  if (bad.by) {
    .err(sprintf("`%s` must be a string containing the name of the variable in data for which weighting is to occur within strata or a one-sided formula with the stratifying variable on the right-hand side",
                 by.arg))
  }

  if (anyNA(by)) {
    .err(sprintf("the variable supplied to `%s` cannot contain any missing (NA) values",
                 by.arg))
  }

  by.components <- data.frame(by)

  names(by.components) <- colnames(by) %or% by.name

  by.factor <- {
    if (is_null(by)) gl(1, n)
    else factor(by.components[[1L]], levels = sort(unique(by.components[[1L]])),
                labels = paste(names(by.components), "=", sort(unique(by.components[[1L]]))))
  }

  if (treat.type != "continuous" &&
      any_apply(levels(by.factor), function(x) nunique(treat) != nunique(treat[by.factor == x]))) {
    .err(sprintf("not all the groups formed by `%s` contain all treatment levels%s. Consider coarsening `%s`",
                 by.arg,
                 if (is_not_null(treat.name)) sprintf(" in %s", treat.name) else "",
                 by.arg))
  }

  attr(by.components, "by.factor") <- by.factor

  by.components
}

.check_num.formula <- function(num.formula, data, env, formula.list) {
  if (rlang::is_formula(num.formula)) {
    if (!rlang::is_formula(num.formula, lhs = FALSE)) {
      .err("the argument to `num.formula` must have right hand side variables but not a response variable (e.g., ~ V1 + V2)")
    }

    rhs.vars.mentioned.lang <- .attr(terms(num.formula), "variables")[-1L]
    rhs.vars.mentioned <- vapply(rhs.vars.mentioned.lang, deparse1, character(1L))
    rhs.vars.failed <- vapply(rhs.vars.mentioned.lang, function(v) {
      null_or_error(try(eval(v, c(data, env)), silent = TRUE))
    }, logical(1L))

    if (any(rhs.vars.failed)) {
      .err(sprintf("All variables in `num.formula` must be variables in `data` or objects in the global environment.\nMissing variables: %s",
                   word_list(rhs.vars.mentioned[rhs.vars.failed], and.or = FALSE)),
           tidy = FALSE)
    }
  }
  else if (is.list(num.formula)) {
    if (length(num.formula) != length(formula.list)) {
      .err("when supplied as a list, `num.formula` must have as many entries as `formula.list`")
    }

    if (!all_apply(num.formula, rlang::is_formula, lhs = FALSE)) {
      .err("`num.formula` must be a single formula with no response variable and with the stabilization factors on the right hand side or a list thereof")
    }

    rhs.vars.mentioned.lang.list <- lapply(num.formula, function(nf) .attr(terms(nf), "variables")[-1L])
    rhs.vars.mentioned <- unique(unlist(lapply(rhs.vars.mentioned.lang.list,
                                               function(r) vapply(r, deparse1, character(1L)))))
    rhs.vars.failed <- vapply(rhs.vars.mentioned, function(v) {
      null_or_error(try(eval(parse(text = v), c(data, env)), silent = TRUE))
    }, logical(1L))

    if (any(rhs.vars.failed)) {
      .err(sprintf("All variables in `num.formula` must be variables in `data` or objects in the global environment.\nMissing variables: %s",
                   word_list(rhs.vars.mentioned[rhs.vars.failed], and.or = FALSE)),
           tidy = FALSE)
    }
  }
  else {
    .err("`num.formula` must be a single formula with no response variable and with the stabilization factors on the right hand side or a list thereof")
  }
}

.make_covs_full_rank <- function(covs) {
  if (ncol(covs) <= 1L) {
    return(covs)
  }

  colinear.covs.to.remove <- setdiff(colnames(covs), colnames(make_full_rank(covs)))

  covs[, colnames(covs) %nin% colinear.covs.to.remove, drop = FALSE]
}

.make_closer_to_1 <- function(x) {
  if (chk::vld_character_or_factor(x) || all_the_same(x)) {
    return(x)
  }

  if (is_binary(x)) {
    return(as.numeric(x == max(x, na.rm = TRUE)))
  }

  (x - mean_fast(x, TRUE)) / sd(x, na.rm = TRUE)
}

.make_covs_closer_to_1 <- function(covs) {
  for (i in seq_col(covs)) {
    covs[, i] <- .make_closer_to_1(covs[, i])
  }

  covs
}

.apply_moments_int_quantile <- function(d, moments = integer(), int = FALSE, quantile = list(),
                                        center = TRUE, orthogonal_poly = TRUE, s.weights = NULL,
                                        focal = NULL, treat = NULL, const = 2000) {

  if (is_null(d)) {
    return(d)
  }

  if (!is.numeric(d)) {
    .err("an error occurred, probably a bug")
  }

  binary.vars <- is_binary_col(d)

  if (center && (int || !orthogonal_poly) && !all(binary.vars)) {
    d[, !binary.vars] <- center_w(d[, !binary.vars, drop = FALSE],
                                  s.weights)
  }

  nd <- NCOL(d)

  # moments
  default_moments <- .attr(moments, "moments_default") %or% 1L
  poly <- setNames(rep.int(default_moments, nd),
                   colnames(d))

  if (length(moments) == 1L && (is_null(names(moments)) || !nzchar(names(moments)))) {
    poly[] <- moments
  }
  else {
    in_poly_and_moments <- intersect(names(poly), names(moments))
    poly[in_poly_and_moments] <- moments[in_poly_and_moments]
  }

  # Binary variables get no poly
  poly[binary.vars] <- pmin(poly[binary.vars], 1L)

  poly_terms <- poly_co.names <- make_list(colnames(d)[poly > 0L])

  for (i in names(poly_terms)) {
    if (poly[i] == 1L) {
      poly_terms[[i]] <- d[, i]
      poly_co.names[[i]] <- i
    }
    else {
      poly_terms[[i]] <- stats::poly(d[, i], degree = poly[i],
                                     raw = !orthogonal_poly, simple = TRUE)
      poly_co.names[[i]] <- sprintf("%s%s%s",
                                    if (orthogonal_poly) "orth_" else "",
                                    i,
                                    num_to_superscript(seq_len(poly[i])))
    }
  }

  # int
  if (int && nd > 1L) {
    ints_to_make <- utils::combn(colnames(d), 2L, simplify = FALSE)

    int_terms <- int_co.names <- make_list(length(ints_to_make))

    for (i in seq_along(ints_to_make)) {
      int_i <- d[, ints_to_make[[i]][1L]] * d[, ints_to_make[[i]][2L]]

      if (!all_the_same(int_i)) {
        int_terms[[i]] <- int_i
        int_co.names[[i]] <- sprintf("%s * %s",
                                     ints_to_make[[i]][1L],
                                     ints_to_make[[i]][2L])
      }
    }
  }
  else {
    int_terms <- int_co.names <- list()
  }

  # quantile
  if (all(binary.vars)) {
    qu <- list()
  }
  else {
    qu <- make_list(colnames(d)[!binary.vars])

    if (length(quantile) == 1L && (is_null(names(quantile)) || !nzchar(names(quantile)))) {
      qu[] <- quantile[rep.int(1L, sum(!binary.vars))]
    }
    else {
      if (any(names(quantile) %in% colnames(d)[binary.vars])) {
        .wrn("ignoring `quantile` constraints on binary and categorical variables")
      }

      in_qu_and_quantile <- intersect(names(qu), names(quantile))
      qu[in_qu_and_quantile] <- quantile[in_qu_and_quantile]
    }

    qu <- clear_null(qu)
  }

  if (is_not_null(qu)) {
    if (is_not_null(focal) && is_not_null(s.weights)) {
      s.weights <- s.weights[treat == focal]
    }

    quantile_terms <- quantile_co.names <- make_list(names(qu))

    for (i in names(qu)) {
      target <- if (is_null(focal)) d[, i] else d[treat == focal, i]

      quantile_terms[[i]] <- do.call("cbind", lapply(qu[[i]], function(q) {
        plogis(const * (d[, i] - w.quantile(target, q, s.weights)))
      }))

      quantile_co.names[[i]] <- paste0(i, "_", qu[[i]])
    }
  }
  else {
    quantile_terms <- quantile_co.names <- list()
  }

  if (is_null(poly_terms) && is_null(int_terms) && is_null(quantile_terms)) {
    return(matrix(ncol = 0L, nrow = nrow(d), dimnames = list(rownames(d), NULL)))
  }

  out <- do.call("cbind", c(poly_terms, int_terms, quantile_terms))
  colnames(out) <- c(unlist(poly_co.names), unlist(int_co.names), unlist(quantile_co.names))

  # Remove single values
  single_value <- apply(out, 2L, all_the_same)

  if (!any(single_value)) {
    return(out)
  }

  out[, !single_value, drop = FALSE]
}

get.s.d.denom.weightit <- function(s.d.denom = NULL, estimand = NULL, weights = NULL, treat = NULL, focal = NULL) {
  s.d.denom.specified <- is_not_null(s.d.denom)
  estimand.specified <- is_not_null(estimand)
  if (!is.factor(treat)) treat <- factor(treat)

  if (s.d.denom.specified) {
    allowable.s.d.denoms <- c("treated", "control", "pooled", "all", "weighted", "hedges")

    try.s.d.denom <- try(match_arg(s.d.denom, allowable.s.d.denoms), silent = TRUE)

    if (!null_or_error(try.s.d.denom)) {
      return(try.s.d.denom)
    }
  }

  if (estimand.specified) {
    allowable.estimands <- c("ATT", "ATC", "ATE", "ATO", "ATM")

    try.estimand <- try(match_arg(toupper(estimand), allowable.estimands), silent = TRUE)

    if (!null_or_error(try.estimand) && try.estimand %nin% c("ATC", "ATT")) {
      s.d.denom <- switch(try.estimand,
                          ATO = "weighted",
                          ATM = "weighted",
                          "pooled")
      return(s.d.denom)
    }
  }

  if (is_not_null(focal)) {
    return(focal)
  }

  if (is_null(weights) || all_the_same(weights)) {
    return("pooled")
  }

  for (tv in levels(treat)) {
    if (all_the_same(weights[treat == tv]) &&
        !all_the_same(weights[treat != tv])) {
      return(tv)
    }
  }

  "pooled"
}

.check_estimated_weights <- function(w, treat, treat.type, s.weights) {

  tw <- w * s.weights

  extreme.warn <- FALSE
  if (all_the_same(w)) {
    .wrn(sprintf("all weights are %s, possibly indicating an estimation failure", w[1L]))
  }
  else if (treat.type == "continuous") {
    w.cv <- sd(tw, na.rm = TRUE) / mean(tw, na.rm = TRUE)
    if (!is.finite(w.cv) || w.cv > 4) extreme.warn <- TRUE
  }
  else {
    t.levels <- unique(treat)
    bad.treat.groups <- setNames(rep.int(FALSE, length(t.levels)), t.levels)
    for (i in t.levels) {
      ti <- which(treat == i)
      if (all(is.na(w[ti])) || all(check_if_zero(w[ti]))) {
        bad.treat.groups[as.character(i)] <- TRUE
      }
      else if (!extreme.warn && sum(is.finite(tw[ti])) > 1L) {
        w.cv <- sd(tw[ti], na.rm = TRUE) / mean(tw[ti], na.rm = TRUE)
        if (!is.finite(w.cv) || w.cv > 4) {
          extreme.warn <- TRUE
        }
      }
    }

    if (any(bad.treat.groups)) {
      n <- sum(bad.treat.groups)
      .wrn(sprintf("all weights are `NA` or 0 in treatment %s %s",
                   ngettext(n, "group", "groups"),
                   word_list(t.levels[bad.treat.groups], quotes = TRUE)))
    }
  }

  if (extreme.warn) {
    .wrn("some extreme weights were generated. Examine them with `summary()` and maybe trim them with `trim()`")
  }

  if (any(tw < 0)) {
    .wrn("some weights are negative; these cannot be used in most model fitting functions")
  }
}

.subclass_ps_multi <- function(ps_mat, treat, estimand = "ATE", focal = NULL, subclass) {
  chk::chk_count(subclass)
  subclass <- round(subclass)

  estimand <- toupper(estimand)

  if (estimand %nin% c("ATE", "ATT")) {
    .err("only the ATE, ATT, and ATC are compatible with stratification weights")
  }

  if (is_not_null(focal)) {
    ps_mat <- ps_mat[, c(focal, setdiff(colnames(ps_mat), focal))]
  }

  ps_sub <- sub_mat <- ps_mat * 0

  for (i in colnames(ps_mat)) {
    if (estimand == "ATE") {
      sub <- findInterval(ps_mat[, as.character(i)],
                          quantile(ps_mat[, as.character(i)],
                                   seq(0, 1, length.out = subclass + 1L)),
                          all.inside = TRUE) |>
        as.integer()
    }
    else if (estimand == "ATT") {
      if (i != focal) {
        ps_mat[, as.character(i)] <- 1 - ps_mat[, as.character(i)]
      }

      sub <- findInterval(ps_mat[, as.character(i)],
                          quantile(ps_mat[treat == focal, as.character(i)],
                                   seq(0, 1, length.out = subclass + 1L)),
                          all.inside = TRUE) |>
        as.integer()
    }

    sub_tab <- table(treat, sub)

    if (any(sub_tab == 0L)) {
      sub <- .subclass_scoot(sub, treat, ps_mat[, i])
      sub_tab <- table(treat, sub)
    }

    sub <- as.character(sub)

    sub_totals <- colSums(sub_tab)
    sub_ps <- setNames(sub_tab[as.character(i), ] / sub_totals,
                       colnames(sub_tab))

    ps_sub[, i] <- sub_ps[sub]
    sub_mat[, i] <- sub

    if (ncol(ps_sub) == 2L) {
      ps_sub[, colnames(ps_sub) != i] <- 1 - ps_sub[, i]
      sub_mat[, colnames(sub_mat) != i] <- sub
      break
    }
  }

  attr(ps_sub, "sub_mat") <- sub_mat

  ps_sub
}

.subclass_ps_bin <- function(ps, treat, estimand = "ATE", subclass) {
  chk::chk_count(subclass)
  subclass <- round(subclass)

  estimand <- toupper(estimand)

  if (estimand %nin% c("ATE", "ATT", "ATC")) {
    .err("only the ATE, ATT, and ATC are compatible with stratification weights")
  }

  sub <- findInterval(ps,
                      quantile(switch(estimand,
                                      ATE = ps,
                                      ATT = ps[treat == 1],
                                      ATC = ps[treat == 0]),
                               seq(0, 1, length.out = subclass + 1L)),
                      all.inside = TRUE) |>
    as.integer()

  max_sub <- max(sub)
  sub_tab1 <- tabulate(sub[treat == 1], max_sub)
  sub_tab0 <- tabulate(sub[treat == 0], max_sub)

  if (any(sub_tab1 == 0L) || any(sub_tab0 == 0L)) {
    sub <- .subclass_scoot(sub, treat, ps)
    sub_tab1 <- tabulate(sub[treat == 1], max_sub)
    sub_tab0 <- tabulate(sub[treat == 0], max_sub)
  }

  sub_totals <- sub_tab1 + sub_tab0
  sub1_prop <- sub_tab1 / sub_totals

  sub_ps <- sub1_prop[sub]

  attr(sub_ps, "sub") <- sub

  sub_ps
}

.subclass_scoot <- function(sub, treat, x, min.n = 1L) {
  #Reassigns subclasses so there are no empty subclasses
  #for each treatment group. min.n is the smallest a
  #subclass is allowed to be.
  treat <- as.character(treat)
  unique.treat <- unique(treat, nmax = 2L)

  names(x) <- seq_along(x)
  names(sub) <- seq_along(sub)
  original.order <- names(x)

  nsub <- nunique(sub)

  #Turn subs into a contiguous sequence
  sub <- setNames(seq_len(nsub), sort(unique(sub)))[as.character(sub)] |>
    setNames(original.order)

  if (any(table(treat) < nsub * min.n)) {
    .err("too many subclasses were requested")
  }

  for (t in unique.treat) {
    if (sum(treat == t) == nsub) {
      sub[treat == t] <- seq_len(nsub)
    }
  }

  sub_tab <- table(treat, sub)

  if (all(sub_tab > 0L)) {
    return(sub)
  }

  .soft_thresh <- function(x, minus = 1) {
    x <- x - minus
    x[x < 0] <- 0
    x
  }

  for (t in unique.treat) {
    for (n in seq_len(min.n)) {
      while (any(sub_tab[t, ] == 0L)) {
        first_0 <- which(sub_tab[t, ] == 0L)[1L]

        if (first_0 == nsub ||
            (first_0 != 1L &&
             sum(.soft_thresh(sub_tab[t, seq_len(first_0 - 1L)]) / abs(first_0 - seq_len(first_0 - 1L))) >=
             sum(.soft_thresh(sub_tab[t, seq(first_0 + 1L, nsub)]) / abs(first_0 - seq(first_0 + 1L, nsub))))) {
          #If there are more and closer nonzero subs to the left...
          first_non0_to_left <- max(seq_len(first_0 - 1L)[sub_tab[t, seq_len(first_0 - 1L)] > 0L])

          name_to_move <- names(sub)[which(x == max(x[treat == t & sub == first_non0_to_left]) &
                                             treat == t & sub == first_non0_to_left)[1L]]

          sub[name_to_move] <- first_0
          sub_tab[t, first_0] <- 1L
          sub_tab[t, first_non0_to_left] <- sub_tab[t, first_non0_to_left] - 1L

        }
        else {
          #If there are more and closer nonzero subs to the right...
          first_non0_to_right <- min(seq(first_0 + 1L, nsub)[sub_tab[t, seq(first_0 + 1L, nsub)] > 0L])
          name_to_move <- names(sub)[which(x == min(x[treat == t & sub == first_non0_to_right]) &
                                             treat == t & sub == first_non0_to_right)[1L]]
          sub[name_to_move] <- first_0
          sub_tab[t, first_0] <- 1L
          sub_tab[t, first_non0_to_right] <- sub_tab[t, first_non0_to_right] - 1L
        }
      }

      sub_tab[t, ] <- sub_tab[t, ] - 1L
    }
  }

  #Unsort
  sub[names(sub)]
}

stabilize_w <- function(weights, treat) {
  t.levels <- {
    if (is.factor(treat)) levels(treat)
    else unique(treat)
  }

  w.names <- names(weights)

  tab <- vapply(t.levels, function(x) mean_fast(treat == x), numeric(1L)) |>
    setNames(t.levels)

  setNames(weights * tab[as.character(treat)], w.names)
}

.get_dens_fun <- function(use.kernel = FALSE, bw = NULL, adjust = NULL, kernel = NULL,
                          n = NULL, treat = NULL, density = NULL, weights = NULL) {
  if (is_null(n)) n <- 10L * length(treat)
  if (is_null(adjust)) adjust <- 1

  if (!isFALSE(use.kernel)) {
    if (isTRUE(use.kernel)) {
      .wrn('`use.kernel` is deprecated; use `density = "kernel"` instead. Setting `density = "kernel"`')
      density <- "kernel"
    }
    else {
      .wrn("`use.kernel` is deprecated")
    }
  }

  if (identical(density, "kernel")) {
    if (is_null(bw)) bw <- "nrd0"
    if (is_null(kernel)) kernel <- "gaussian"

    densfun <- function(p, log = FALSE) {
      d <- stats::density(p, n = n,
                          weights = weights / sum(weights),
                          give.Rkern = FALSE,
                          bw = bw,
                          adjust = adjust,
                          kernel = kernel)

      out <- with(d, approxfun(x = x, y = y))(p)

      if (log) out <- log(out)

      attr(out, "density") <- d
      out
    }
  }
  else {
    if (is_null(density)) .density <- function(x, log = FALSE) dnorm(x, log = log)
    else if (is.function(density)) .density <- function(x, log = FALSE) {
      if (utils::hasName(formals(density), "log")) density(x, log = log)
      else if (log) log(density(x))
      else density(x)
    }
    else if (identical(density, "dlaplace")) .density <- function(x, log = FALSE) {
      mu <- 0
      b <- 1
      if (log)
        -abs(x - mu) / b - log(2 * b)
      else
        exp(-abs(x - mu) / b) / (2 * b)
    }
    else if (chk::vld_string(density)) {
      splitdens <- strsplit(density, "_", fixed = TRUE)[[1L]]

      splitdens1 <- get0(splitdens[1L], mode = "function", envir = parent.frame())

      if (is_null(splitdens1)) {
        .err(sprintf("%s is not an appropriate argument to `density` because %s is not an available function",
                     density, splitdens[1L]))
      }

      if (length(splitdens) > 1L && !can_str2num(splitdens[-1L])) {
        .err(sprintf("%s is not an appropriate argument to `density` because %s cannot be coerced to numeric",
                     density, word_list(splitdens[-1L], and.or = "or", quotes = TRUE)))
      }

      .density <- function(x, log = FALSE) {
        if (utils::hasName(formals(splitdens1), "log")) {
          out <- tryCatch(do.call(splitdens1, c(list(x, log = log), as.list(str2num(splitdens[-1L])))),
                          error = function(e) {
                            .err(sprintf("Error in applying density:\n  %s",
                                         conditionMessage(e)),
                                 tidy = FALSE)
                          })
        }
        else {
          out <- tryCatch(do.call(splitdens1, c(list(x), as.list(str2num(splitdens[-1L])))),
                          error = function(e) {
                            .err(sprintf("Error in applying density:\n  %s",
                                         conditionMessage(e)),
                                 tidy = FALSE)
                          })

          if (log) out <- log(out)
        }

        out
      }
    }
    else {
      .err("the argument to `density` cannot be evaluated as a density function")
    }

    densfun <- function(p, log = FALSE) {
      # sd <- sd(p)
      # sd <- sqrt(col.w.v(p, s.weights))
      dens <- .density(p, log = log)
      if (is_null(dens) || !is.numeric(dens) || anyNA(dens)) {
        .err("there was a problem with the output of `density`. Try another density function or leave it blank to use the Gaussian density")
      }

      if ((log && !all(is.finite(dens))) ||
          (!log && any(dens <= 0))) {
        .err("the input to density may not accept the full range of standardized treatment values or residuals")
      }

      x <- seq.int(min(p) - 3 * adjust * bw.nrd0(p),
                   max(p) + 3 * adjust * bw.nrd0(p),
                   length.out = n)
      attr(dens, "density") <- data.frame(x = x,
                                          y = .density(x, log = log))
      dens
    }
  }

  densfun
}

.get_w_from_ps_internal_bin <- function(ps, treat, estimand = "ATE",
                                        subclass = NULL, stabilize = FALSE) {

  estimand <- toupper(estimand)
  w <- rep_with(1, treat)

  #Assume treat is binary
  if (is_not_null(subclass)) {
    #Get MMW subclass propensity scores
    ps <- .subclass_ps_bin(ps, treat, estimand, subclass)
  }

  i0 <- which(treat == 0)

  if (estimand == "ATE") {
    w[i0] <- 1 / (1 - ps[i0])
    w[-i0] <- 1 / ps[-i0]
  }
  else if (estimand == "ATT") {
    w[i0] <- .p2o(ps[i0])
  }
  else if (estimand == "ATC") {
    w[-i0] <- .p2o(1 - ps[-i0])
  }
  else if (estimand == "ATO") {
    w[i0] <- ps[i0]
    w[-i0] <- 1 - ps[-i0]
  }
  else if (estimand == "ATM") {
    w[i0][ps[i0] < .5] <- .p2o(ps[i0][ps[i0] < .5])
    w[-i0][ps[-i0] > .5] <- .p2o(1 - ps[-i0][ps[-i0] > .5])
  }
  else if (estimand == "ATOS") {
    w[i0] <- 1 / (1 - ps[i0])
    w[-i0] <- 1 / ps[-i0]

    ps.sorted <- sort(c(ps, 1 - ps))
    z <- ps * (1 - ps)
    alpha.opt <- 0
    for (i in seq_len(sum(ps < .5))) {
      if (i == 1L || !check_if_zero(ps.sorted[i] - ps.sorted[i - 1L])) {
        alpha <- ps.sorted[i]
        a <- alpha * (1 - alpha)
        if (2 * a * sum(1 / z[z >= a]) / sum(z >= a) >= 1) {
          alpha.opt <- alpha
          break
        }
      }
    }
    w[!between(ps, c(alpha.opt, 1 - alpha.opt))] <- 0
  }

  names(w) <- names(treat) %or% NULL

  if (stabilize) {
    w <- stabilize_w(w, treat)
  }

  w
}

.get_w_from_ps_internal_multi <- function(ps, treat, estimand = "ATE", focal = NULL,
                                          subclass = NULL, stabilize = FALSE) {

  estimand <- toupper(estimand)
  w <- rep_with(0, treat)

  ps_mat <- ps

  if (is_not_null(subclass)) {
    #Get MMW subclass propensity scores
    ps_mat <- .subclass_ps_multi(ps_mat, treat, estimand, focal, subclass)
  }

  for (i in colnames(ps_mat)) {
    w[treat == i] <- 1 / ps_mat[treat == i, i]
  }

  if (estimand == "ATE") {
    # w <- w
  }
  else if (estimand %in% c("ATT", "ATC")) {
    in_f <- which(treat == focal)
    w[in_f] <- 1
    w[-in_f] <- w[-in_f] * ps_mat[-in_f, as.character(focal)]
  }
  else if (estimand == "ATO") {
    w <- w / rowSums(1 / ps_mat) #Li & Li (2019)
  }
  else if (estimand == "ATM") {
    w <- w * do.call("pmin", lapply(seq_col(ps_mat), function(x) ps_mat[, x]), quote = TRUE)
  }
  else if (estimand == "ATOS") {
    #Crump et al. (2009)
    ps.sorted <- sort(c(ps_mat[, 2L], 1 - ps_mat[, 2L]))
    z <- ps_mat[, 2L] * (1 - ps_mat[, 2L])
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
    w[!between(ps_mat[, 2L], c(alpha.opt, 1 - alpha.opt))] <- 0
  }
  else {
    return(numeric(0L))
  }

  if (stabilize) {
    w <- stabilize_w(w, treat)
  }

  names(w) <- rownames(ps_mat) %or% names(treat) %or% NULL

  w
}

.get_w_from_ps_internal_array <- function(ps, treat, estimand = "ATE", focal = NULL,
                                          subclass = NULL, stabilize = FALSE) {
  #Batch turn PS into weights; primarily for output of predict.gbm
  # Assumes a (0,1) treatment if binary
  if (length(dim(ps)) <= 1L) {
    ps <- matrix(ps, ncol = 1L)
  }

  eps <- 1e-8

  if (length(dim(ps)) == 2L) {
    #Binary treatment, vector ps

    w <- ps
    w[] <- 0

    if (is_not_null(subclass)) {
      #Get MMW subclass propensity scores
      for (p in seq_col(ps)) {
        ps[, p] <- .subclass_ps_bin(ps[, p], treat, estimand, subclass)
      }
    }

    t1 <- which(treat == 1)
    t0 <- which(treat == 0)

    if (estimand == "ATE") {
      ps[t1, ][ps[t1, ] < eps] <- eps
      ps[t0, ][ps[t0, ] > 1 - eps] <- 1 - eps

      w[t1, ] <- 1 / ps[t1, ]
      w[t0, ] <- 1 / (1 - ps[t0, ])
    }
    else if (estimand == "ATT") {
      ps[t0, ][ps[t0, ] > 1 - eps] <- 1 - eps

      w[t1, ] <- 1
      w[t0, ] <- .p2o(ps[t0, ])
    }
    else if (estimand == "ATC") {
      ps[t1, ][ps[t1, ] < eps] <- eps

      w[t1, ] <- .p2o(1 - ps[t1, ])
      w[t0, ] <- 1
    }
    else if (estimand == "ATO") {
      w[t1, ] <- 1 - ps[t1, ]
      w[t0, ] <- ps[t0, ]
    }
    else if (estimand == "ATM") {
      pslt.5 <- ps < .5
      w[t1, ][pslt.5[t1, ]] <- 1
      w[t1, ][!pslt.5[t1, ]] <- .p2o(1 - ps[t1, ][!pslt.5[t1, ]])

      w[t0, ][pslt.5[t0, ]] <- .p2o(ps[t0, ][pslt.5[t0, ]])
      w[t0, ][!pslt.5[t0, ]] <- 1
    }

    if (stabilize) {
      w[t1] <- w[t1] * length(t1) / length(treat)
      w[t0] <- w[t0] * length(t0) / length(treat)
    }
  }
  else if (length(dim(ps)) == 3L) {
    #Multi-category treatment, matrix PS

    if (is_not_null(subclass)) {
      #Get MMW subclass propensity scores
      for (p in seq_len(last(dim(ps))))
        ps[, , p] <- .subclass_ps_multi(ps[, , p], treat, estimand, focal, subclass)
    }

    ps <- squish(ps, eps)

    w <- matrix(0.0, ncol = dim(ps)[3L], nrow = dim(ps)[1L])
    t.levs <- unique(treat)

    for (i in t.levs) {
      w[treat == i, ] <- 1 / ps[treat == i, as.character(i), ]
    }

    if (estimand == "ATE") {
      #Do nothing
    }
    else if (estimand %in% c("ATT", "ATC")) {
      not_focal <- which(treat != focal)
      w[-not_focal, ] <- 1
      w[not_focal, ] <- w[not_focal, ] * ps[not_focal, as.character(focal), ]
    }
    else if (estimand == "ATO") {
      w <- w / colSums(aperm(1 / ps, c(2L, 1L, 3L)))
    }
    else if (estimand == "ATM") {
      treat <- as.integer(treat)

      for (p in seq_len(dim(ps)[3L])) {
        ps_p <- ps[, , p]
        min_ind <- max.col(-ps_p, ties.method = "first")
        no_match <- which(ps_p[cbind(seq_along(treat), treat)] != ps_p[cbind(seq_along(treat), min_ind)])

        if (length(no_match) < length(treat)) {
          w[-no_match, p] <- 1
        }

        if (is_not_null(no_match)) {
          w[no_match, p] <- w[no_match, p] * ps_p[cbind(no_match, min_ind[no_match])]
        }
      }
    }

    if (stabilize) {
      for (i in t.levs) {
        w[treat == i, ] <- mean_fast(treat == i) * w[treat == i, ]
      }
    }
  }
  else {
    .err("don't know how to process more than 3 dims (likely a bug)")
  }

  w
}

#Derivative of weights wrt ps for different estimands
.dw_dp_bin <- function(p, treat, estimand = "ATE") {
  estimand <- toupper(estimand)

  dw <- rep_with(0, treat)

  i0 <- which(treat == 0)

  if (estimand == "ATE") {
    dw[i0] <- (1 - p[i0])^(-2)
    dw[-i0] <- -p[-i0]^(-2)
  }
  else if (estimand == "ATT") {
    dw[i0] <- (1 - p[i0])^(-2)
  }
  else if (estimand == "ATC") {
    dw[-i0] <- -p[-i0]^(-2)
  }
  else if (estimand == "ATO") {
    dw[i0] <- 1
    dw[-i0] <- -1
  }
  else if (estimand == "ATM") {
    dw[i0][p[i0] < .5] <- (1 - p[i0][p[i0] < .5])^(-2)
    dw[-i0][p[-i0] > .5] <- -p[-i0][p[-i0] > .5]^(-2)
  }

  dw
}

#Derivative of weights wrt ps for different estimands
.dw_dp_multi <- function(p, treat, estimand = "ATE", focal = NULL) {
  estimand <- toupper(estimand)

  dw <- array(0, dim = dim(p), dimnames = dimnames(p))

  pA <- numeric(nrow(p))

  for (k in levels(treat)) {
    pA[treat == k] <- p[treat == k, k]
  }

  if (is_not_null(focal)) {
    pF <- p[, focal]

    for (i in setdiff(levels(treat), focal)) {
      dw[treat == i, focal] <- 1 / pA[treat == i]
      dw[treat == i, i] <- -pF[treat == i] / pA[treat == i]^2
    }
  }
  else if (estimand == "ATE") {
    for (k in levels(treat)) {
      dw[treat == k, k] <- -1 / pA[treat == k]^2
    }
  }
  else if (estimand == "ATO") {
    S <- 1 / rowSums(1 / p)
    for (k in levels(treat)) {
      dw[treat == k, k] <- (S[treat == k] / pA[treat == k]^2) * (S[treat == k] / pA[treat == k] - 1)

      for (j in setdiff(levels(treat), k)) {
        dw[treat == k, j] <- (1 / pA[treat == k]) * (S[treat == k] / p[treat == k, j])^2
      }
    }
  }
  else if (estimand == "ATM") {
    M <- do.call("pmin", as.data.frame(p))

    for (k in levels(treat)) {
      m1 <- p[, k] == M & treat != k
      dw[m1, k] <- 1 / pA[m1]

      m2 <- p[, k] != M & treat == k
      dw[m2, k] <- -M[m2] / pA[m2]^2
    }
  }

  dw
}

plot_density <- function(d.n, d.d, log = FALSE) {
  d.d <- cbind(as.data.frame(d.d[c("x", "y")]), dens = "Denominator Density", stringsAsfactors = FALSE)
  d.n <- cbind(as.data.frame(d.n[c("x", "y")]), dens = "Numerator Density", stringsAsfactors = FALSE)
  d.all <- rbind(d.d, d.n)
  d.all$dens <- factor(d.all$dens, levels = c("Numerator Density", "Denominator Density"))

  if (log) {
    d.all$x <- exp(d.all$x)
  }

  pl <- ggplot(d.all, aes(x = .data$x, y = .data$y)) +
    geom_line() +
    labs(title = "Weight Component Densities", x = "E[Treat|X]", y = "Density") +
    facet_grid(rows = vars(.data$dens)) +
    theme(panel.background = element_rect(fill = "white"),
          panel.border = element_rect(fill = NA, color = "black"),
          axis.text.x = element_text(color = "black"),
          axis.text.y = element_text(color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  print(pl)
}

neg_ent <- function(w) {
  w <- w[w > 0]
  w <- w / mean_fast(w)
  mean(w * log(w))
}

replace_na_with <- function(covs, with = "median") {
  if (is.na(with) || !anyNA(covs)) {
    return(covs)
  }

  if (is.character(with)) {
    .with <- match.fun(with)
    for (i in colnames(covs)[anyNA_col(covs)]) {
      if (all(is.na(covs[, i]))) covs <- covs[, colnames(covs) != i, drop = FALSE]
      else covs[is.na(covs[, i]), i] <- .with(covs[, i], na.rm = TRUE)
    }
  }
  else {
    covs[is.na(covs)] <- with
  }

  covs
}

add_missing_indicators <- function(covs, replace_with = "median") {
  covs_w_missing <- which(anyNA_col(covs))
  if (is_null(covs_w_missing)) {
    return(covs)
  }

  missing_ind <- apply(covs[, covs_w_missing, drop = FALSE], 2L, function(x) as.numeric(is.na(x)))

  colnames(missing_ind) <- paste0(colnames(missing_ind), ":<NA>")
  covs <- cbind(covs, missing_ind)

  if (is_null(replace_with) || is.na(replace_with)) {
    return(covs)
  }

  replace_na_with(covs, replace_with)
}

verbosely <- function(expr, verbose = TRUE) {
  if (verbose) {
    return(expr)
  }

  invisible(utils::capture.output({
    out <- invisible(expr)
  }))

  out
}

#Generalized matrix inverse (port of MASS::ginv)
generalized_inverse <- function(sigma, .try = TRUE) {

  if (!.try) {
    sigmasvd <- svd(sigma)
    pos <- sigmasvd$d > max(1e-9 * sigmasvd$d[1L], 0)
    sigma_inv <- sigmasvd$v[, pos, drop = FALSE] %*% (sigmasvd$d[pos]^(-1) * t(sigmasvd$u[, pos, drop = FALSE]))

    return(sigma_inv)
  }

  tryCatch(solve(sigma),
           error = function(e) {
             generalized_inverse(sigma, .try = FALSE)
           })
}

#Compute gradient numerically using centered difference
.gradient <- function(.f, .x, .eps = 1e-8, .parm = NULL, .direction = "center", .method = "fd", ...) {
  .method <- match_arg(.method, c("fd", "richardson"))

  if (.method == "fd") {
    .gradientFD(.f = .f, .x = .x, .eps = .eps, .parm = .parm, .direction = .direction, ...)
  }
  else if (.method == "richardson") {
    .gradientRich(.f = .f, .x = .x, .eps = .eps, .parm = .parm, .direction = .direction, ...)
  }
}

#Finite difference gradient
.gradientFD <- function(.f, .x, .eps = 1e-8, .parm = NULL, .direction = "center", ...) {

  .direction <- match_arg(.direction, c("center", "left", "right"))

  if (is_null(.parm)) {
    .parm <- seq_along(.x)
  }

  .x0 <- .x

  .eps <- squish(abs(.x) * .eps, lo = .eps, hi = Inf)

  if (.direction != "center") {
    .f0 <- .f(.x0, ...)
  }

  for (jj in seq_along(.parm)) {
    j <- .parm[jj]

    if (.direction == "center") {
      .x[j] <- .x0[j] + .eps[j] / 2

      f_new_r <- .f(.x, ...)
    }
    else if (.direction == "left") {
      f_new_r <- .f0
    }
    else if (.direction == "right") {
      .x[j] <- .x0[j] + .eps[j]

      f_new_r <- .f(.x, ...)
    }

    if (j == 1L) {
      jacob <- matrix(0, nrow = length(f_new_r), ncol = length(.parm),
                      dimnames = list(names(f_new_r), names(.x)[.parm]))
    }

    if (.direction == "center") {
      .x[j] <- .x0[j] - .eps[j] / 2

      f_new_l <- .f(.x, ...)
    }
    else if (.direction == "left") {
      .x[j] <- .x0[j] - .eps[j]

      f_new_l <- .f(.x, ...)
    }
    else if (.direction == "right") {
      f_new_l <- .f0
    }

    jacob[, jj] <- (f_new_r - f_new_l) / .eps[j]

    .x[j] <- .x0[j]
  }

  jacob
}

#Using Richardson extrapolation
.gradientRich <- function(.f, .x, .eps = 1e-8, .parm = NULL, .direction = "center", ...) {

  .direction <- match_arg(.direction, c("center", "left", "right"))

  if (is_null(.parm)) {
    .parm <- seq_along(.x)
  }

  if (.direction != "center") {
    .f0 <- .f(.x, ...)
  }

  n <- length(.x)

  d <- 1e-4
  r <- 4
  v <- 2
  a <- NULL
  h <- abs(d * .x) + .eps * (abs(.x) < 1e-5)

  for (k in seq_len(r)) {
    eps_i <- rep_with(0, .parm)

    for (ii in seq_along(.parm)) {
      i <- .parm[ii]

      eps_i[i] <- h[i]

      a_k_ii <- switch(.direction,
                       center = (.f(.x + eps_i, ...) - .f(.x - eps_i, ...)) / (2 * h[i]),
                       right = (.f(.x + 2 * eps_i, ...) - .f0) / (2 * h[i]),
                       left = (.f0 - .f(.x - 2 * eps_i, ...)) / (2 * h[i]))

      if (is_null(a)) {
        a <- array(NA_real_, dim = c(length(a_k_ii), r, n))
      }

      a[, k, ii] <- a_k_ii

      eps_i[i] <- 0
    }

    h <- h / v
  }

  for (m in seq_len(r - 1L)) {
    a <- (a[, 1L + seq_len(r - m), , drop = FALSE] * (4^m) - a[, seq_len(r - m), , drop = FALSE]) / (4^m - 1)
  }

  array(a, dim = dim(a)[c(1L, 3L)],
        dimnames = list(names(a_k_ii), names(.x)[.parm]))
}

#Convert probability to odds
.p2o <- function(p) {
  p / (1 - p)
}

#Get psi function (individual contributions to gradient) from glm fit
.get_glm_psi <- function(fit) {
  fam <- fit$family

  if (is_null(fam) || identical(fam, gaussian(), ignore.environment = TRUE)) {
    psi <- function(B, X, y, weights, offset = 0) {
      p <- drop(X %*% B) + offset

      X * (weights * (y - p))
    }
  }
  else if (inherits(fit, "brglmFit") &&
           !identical(fit$type, "ML") &&
           !identical(fit$type, "correction")) {
    br_type <- fit$type

    if (is_null(fit$control[["a"]])) {
      rlang::check_installed("brglm2")
      fit$control[["a"]] <- eval(formals(brglm2::brglmControl)[["a"]])
    }

    br_psi <- function(X, W, d, p, XB, V) {
      DD <- fam$d2mu.deta(XB)
      Wt <- W * d^2 / V #"working weight"

      ## Compute hat values
      XWt <- sqrt(Wt) * X
      qrXWT <- qr(XWt)
      Q <- qr.Q(qrXWT)
      H <- rowSums(Q * Q)

      if (br_type %in% c("AS_mixed", "AS_mean")) {
        AA <- .5 * X * H * DD / d
        return(AA)
      }

      V1 <- fam$d1variance(p)

      if (br_type == "MPL_Jeffreys") {
        return(fit$control[["a"]] * X * H * (2 * DD / d - V1 * d / V))
      }

      #br_type == "AS_median"
      R_matrix <- qr.R(qrXWT)
      info_unscaled <- crossprod(R_matrix)
      inverse_info_unscaled <- chol2inv(R_matrix)

      b_vector <- vapply(seq_col(X), function(j) {
        inverse_info_unscaled_j <- inverse_info_unscaled[j, ]
        vcov_j <- tcrossprod(inverse_info_unscaled_j) / inverse_info_unscaled_j[j]
        hats_j <- rowSums((X %*% vcov_j) * X) * Wt

        inverse_info_unscaled_j %*% colSums(X * (hats_j * (d * V1 / (6 * V) - DD / (2 * d))))
      }, numeric(1L))

      AA <- .5 * X * H * DD / d
      sweep(AA, 2L, info_unscaled %*% b_vector / nrow(X), "+")
    }

    psi <-  function(B, X, y, weights, offset = 0) {
      XB <- drop(X %*% B) + offset
      p <- fam$linkinv(XB)
      d <- fam$mu.eta(XB)
      V <- fam$variance(p)

      .psi <- X * (weights * d * (y - p) / V)

      .psi + br_psi(X, weights, d, p, XB, V)
    }
  }
  else if (identical(fam, binomial(), ignore.environment = TRUE) ||
           identical(fam, quasibinomial(), ignore.environment = TRUE) ||
           identical(fam, poisson(), ignore.environment = TRUE) ||
           identical(fam, quasipoisson(), ignore.environment = TRUE)) {
    psi <- function(B, X, y, weights, offset = 0) {
      XB <- drop(X %*% B) + offset
      p <- fam$linkinv(XB)

      X * (weights * (y - p))
    }
  }
  else {
    psi <- function(B, X, y, weights, offset = 0) {
      XB <- drop(X %*% B) + offset
      p <- fam$linkinv(XB)
      d <- fam$mu.eta(XB)
      V <- fam$variance(p)

      X * (weights * d * (y - p) / V)
    }
  }

  psi
}

.make_link <- function(link) {
  link0 <- try(make.link(link), silent = TRUE)

  if (!null_or_error(link0)) {
    return(link)
  }

  if (!chk::vld_string(link) || !link %in% c("clog", "loglog")) {
    .err("link function not recognized")
  }

  if (link == "clog") {
    linkfun <- function(mu) -log(1 - mu)
    linkinv <- function(eta) (1 - exp(-eta)) |> squish(-Inf, 1 - .Machine$double.eps)
    mu.eta <- function(eta) exp(-eta) |> squish(.Machine$double.eps, Inf)
    valideta <- function(eta) TRUE
    name <- "clog"
  }
  else if (link == "loglog") {
    linkfun <- function(mu) -log(-log(mu))
    linkinv <- function(eta) exp(-exp(-eta)) |> squish(.Machine$double.eps)
    mu.eta <- function(eta) {
      eta <- squish(eta, -Inf, 700)
      exp(-eta - exp(-eta)) |> squish(.Machine$double.eps, Inf)
    }
    valideta <- function(eta) TRUE
    name <- "loglog"
  }

  out <- list(linkfun = linkfun,
              linkinv = linkinv,
              mu.eta = mu.eta,
              valideta = valideta,
              name = name)
  class(out) <- "link-glm"
  out
}

.get_glm_starting_values <- function(X, Y, w, family, offset = NULL) {

  if (is_null(w)) {
    w <- rep_with(1, Y)
  }

  if (is_null(offset)) {
    offset <- rep_with(0, Y)
  }

  mustart <- .25 + .5 * Y

  suppressWarnings({
    fit <- try(glm.fit(X, Y, weights = w, offset = offset, family = family,
                       mustart = mustart, control = list(maxit = 1e4L)), silent = TRUE)
  })

  if (!null_or_error(fit) && isTRUE(fit$converged)) {
    return(fit$coefficients)
  }

  coef_start <- c(family$linkfun(w.m(Y, w)), rep.int(0, ncol(X) - 1L))

  suppressWarnings({
    fit <- try(glm.fit(X, Y, weights = w, offset = offset, family = family,
                       start = coef_start, control = list(maxit = 1e4L)), silent = TRUE)
  })

  if (null_or_error(fit) || !isTRUE(fit$converged)) {
    return(coef_start)
  }

  fit$coefficients
}

came_from_weightit <- function(obj) {
  cl <- getCall(obj)

  if (is_null(cl) || !is.call(cl)) {
    return(FALSE)
  }

  env <- obj[["env"]] %or% globalenv()

  fn <- try(eval(cl[[1L]], envir = env), silent = TRUE)

  if (!is.function(fn)) {
    return(FALSE)
  }

  identical(fn, weightit, ignore.environment = TRUE) ||
    identical(fn, weightitMSM, ignore.environment = TRUE)
}