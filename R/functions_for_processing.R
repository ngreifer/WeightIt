.method_to_proper_method <- function(method) {
  if (is_null(method)) return(NULL)

  if (!is.character(method)) {
    return(method)
  }

  method <- tolower(method)

  if (method %nin% unlist(grab(.weightit_methods, "alias"))) {
    return(method)
  }

  .allowable.methods <- unlist(lapply(names(.weightit_methods), function(m) {
    alias <- .weightit_methods[[m]]$alias
    setNames(rep(m, length(alias)), alias)
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
                                .method_to_proper_method(method) %nin% names(.weightit_methods)))) {
    .err(sprintf("`method` must be a string of length 1 containing the name of an acceptable weighting method or a function that produces weights. Allowable methods:\n%s",
                 word_list(names(.weightit_methods), and.or = FALSE, quotes = 2)),
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
      (method %in% names(.weightit_methods)) &&
      (treat.type %nin% .weightit_methods[[method]]$treat_type)) {
    .err(sprintf("%s can only be used with a %s treatment",
                 .method_to_phrase(method),
                 word_list(.weightit_methods[[method]]$treat_type, and.or = "or")))
  }
}

.process.s.weights <- function(s.weights, data = NULL) {
  #Process s.weights
  if (is_null(s.weights)) return(NULL)

  if (is.numeric(s.weights)) return(s.weights)

  if (!is.character(s.weights) || length(s.weights) != 1) {
    .err("the argument to `s.weights` must be a vector or data frame of sampling weights or the (quoted) names of the variable in `data` that contains sampling weights")
  }

  if (is_null(data)) {
    .err("`s.weights` was specified as a string but there was no argument to `data`")
  }

  if (s.weights %nin% names(data)) {
    .err("the name supplied to `s.weights` is not the name of a variable in `data`")
  }

  data[[s.weights]]
}

.check_method_s.weights <- function(method, s.weights) {
  if (is_not_null(method) &&
      !is.function(method) &&
      !.weightit_methods[[method]]$s.weights_ok &&
      !all_the_same(s.weights)) {
    .err(sprintf("sampling weights cannot be used with %s", .method_to_phrase(method)))
  }
}

.method_to_phrase <- function(method) {

  if (is_null(method))
    return("no weighting")

  if (is.function(method))
    return("a user-defined method")

  method <- .method_to_proper_method(method)

  if (method %nin% names(.weightit_methods))
    return("the chosen method of weighting")

  .weightit_methods[[method]]$description
}

.process_estimand <- function(estimand, method, treat.type) {

  if (is.function(method)) {
    chk::chk_null_or(estimand, vld = chk::vld_string)
    return(toupper(estimand))
  }

  if (treat.type == "continuous") {
    if (is_not_null(estimand) && !identical(toupper(estimand), "ATE")) {
      .wrn("`estimand` is ignored for continuous treatments")
    }

    return("ATE")
  }

  chk::chk_string(estimand)
  estimand <- toupper(estimand)

  allowable_estimands <- {
    if (is_null(method)) unique(unlist(grab(.weightit_methods, "estimand")))
    else .weightit_methods[[method]]$estimand
  }

  if (treat.type == "multi-category") {
    allowable_estimands <- setdiff(allowable_estimands, "ATOS")
  }

  if (estimand %nin% allowable_estimands) {
    .err(sprintf("%s is not an allowable estimand for %s with a %s treatment. Only %s allowed",
                 add_quotes(estimand), .method_to_phrase(method), treat.type,
                 word_list(allowable_estimands, quotes = TRUE, and.or = "and", is.are = TRUE)))
  }

  estimand
}

.check_subclass <- function(method, treat.type) {
  if (is_not_null(method) && !is.function(method)) {

    subclass_ok <- .weightit_methods[[method]]$subclass_ok

    if (treat.type == "continuous" || !subclass_ok) {
      .err(sprintf("subclasses are not compatible with %s with a %s treatment",
                   .method_to_phrase(method), treat.type))
    }
  }
}

.process_moments_int_quantile <- function(moments, int, quantile = NULL, method) {
  if (is.function(method)) {
    return(list(moments = moments, int = int, quantile = quantile))
  }

  if (is_null(method) || !.weightit_methods[[method]]$moments_int_ok) {
    if (is_not_null(method) &&
      any(mi0 <- c(is_not_null(moments), is_not_null(int) && !isFALSE(int), is_not_null(quantile)))) {
      .wrn(sprintf("%s not compatible with %s. Ignoring %s",
                   word_list(c("moments", "int", "quantile")[mi0], and.or = "and", is.are = TRUE, quotes = "`"),
                   .method_to_phrase(method),
                   word_list(c("moments", "int", "quantile")[mi0], and.or = "and", quotes = "`")))
    }

    return(list(moments = integer(), int = FALSE, quantile = list()))
  }

  chk::chk_flag(int)

  if (is_not_null(quantile)) {
    .vld_qu <- function(x) {
      is.numeric(x) && all(x >= 0) && all(x <= 1)
    }

    bad.q <- FALSE
    if (is.numeric(quantile) && .vld_qu(quantile)) {
      if (length(quantile) == 1L || (is_not_null(names(quantile)) && !any(names(quantile) == ""))) {
        quantile <- as.list(quantile)
      }
      else {
        bad.q <- TRUE
      }
    }
    else if (is.list(quantile)) {
      if ((length(quantile) > 1L && (is_null(names(quantile)) || any(names(quantile) == ""))) ||
          !all(vapply(quantile, .vld_qu, logical(1L)))) {
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
    chk::chk_whole_number(moments)

    chk::chk_gte(moments,
                 if (is_null(quantile)) .weightit_methods[[method]]$moments_default
                 else 0)

    if (int && moments < 1) {
      .wrn("when `int = TRUE`, `moments` must be greater than or equal to 1. Setting `moments = 1`")
      moments <- 1L
    }
    else {
      moments <- as.integer(moments)
    }
  }
  else {
    moments <- {
      if (int) 1L
      else .weightit_methods[[method]]$moments_default
    }
  }

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
    if (is_null(is.MSM.method)) return(TRUE)

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
    return(allowable.missings[1])
  }

  chk::chk_string(missing)

  if (missing %nin% allowable.missings) {
    .err(sprintf("only %s allowed for `missing` with %s",
                 word_list(allowable.missings, quotes = 2, is.are = TRUE),
                 .method_to_phrase(method)))
  }

  missing
}

.missing_to_phrase <- function(missing) {
  switch(missing,
         "ind" = "missingness indicators",
         "saem" = "SAEM",
         "surr" = "surrogate splitting",
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

.process_ps <- function(ps, data = NULL, treat) {
  if (is_null(ps)) return(NULL)

  if (chk::vld_string(ps)) {
    if (is_null(data)) {
      .err("`ps` was specified as a string but there was no argument to `data`")
    }

    if (ps %nin% names(data)) {
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

  unique.treat <- unique(treat, nmax = switch(treat.type, "binary" = 2, "multi-category" = length(treat)/4))

  #Check focal
  if (is_not_null(focal) && (length(focal) > 1L || focal %nin% unique.treat)) {
    .err("the argument supplied to `focal` must be the name of a level of treatment")
  }

  if (treat.type == "multi-category") {

    if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
      .wrn(sprintf("`estimand = %s` is not compatible with `focal`. Setting `estimand` to \"ATT\"",
                   add_quotes(estimand)))
      reported.estimand <- estimand <- "ATT"
    }

    if (estimand == "ATT") {
      if (is_null(focal)) {
        if (is_null(treated) || treated %nin% unique.treat) {
          .err("when `estimand = \"ATT\"` for multi-category treatments, an argument must be supplied to `focal`")
        }
        focal <- treated
      }
    }
    else if (estimand == "ATC") {
      if (is_null(focal)) {
        .err("when `estimand = \"ATC\"` for multi-category treatments, an argument must be supplied to `focal`")
      }
      estimand <- "ATT"
    }
  }
  else if (treat.type == "binary") {
    unique.treat.bin <- unique(binarize(treat), nmax = 2)

    if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
      .wrn(sprintf("`estimand = %s` is not compatible with `focal`. Setting `estimand` to \"ATT\"",
                   add_quotes(estimand)))
      reported.estimand <- estimand <- "ATT"
    }

    if (is_null(treated) || treated %nin% unique.treat) {
      if (is_null(focal)) {
        if (all(as.character(unique.treat.bin) == as.character(unique.treat))) {
          treated <- unique.treat[unique.treat.bin == 1]
        }
        else {
          treated <- {
            if (is.factor(treat)) levels(treat)[2]
            else unique.treat[unique.treat.bin == 1]
          }

          if (estimand == "ATT") {
            .msg(sprintf("assuming %s the treated level. If not, supply an argument to `focal`",
                         word_list(treated, quotes = !is.numeric(treat), is.are = TRUE)))

          }
          else if (estimand == "ATC") {
            .msg(sprintf("assuming %s the control level. If not, supply an argument to `focal`",
                         word_list(setdiff(unique.treat, treated), quotes = !is.numeric(treat), is.are = TRUE)))
          }

        }

        focal <- switch(estimand,
                        "ATT" = treated,
                        "ATC" = setdiff(unique.treat, treated))
      }
      else {
        treated <- switch(estimand,
                          "ATT" = focal,
                          "ATC" = setdiff(unique.treat, focal))
      }

    }
    else if (is_null(focal)) {
      focal <- switch(estimand,
                      "ATT" = treated,
                      "ATC" = setdiff(unique.treat, treated))
    }
  }

  list(focal = focal,
       estimand = estimand,
       reported.estimand = reported.estimand,
       treated = if (is.factor(treated)) as.character(treated) else treated)
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
  else if (chk::vld_string(by) && by %in% names(data)) {
    by.name <- by
    by <- data[[by]]
  }
  else if (length(dim(by)) == 2L && len(by) == n) {
    by.name <- colnames(by)[1]
    by <- drop(by[, 1])
  }
  else if (rlang::is_formula(by, lhs = FALSE)) {
    t.c <- get_covs_and_treat_from_formula(by, data)
    by <- t.c[["reported.covs"]]
    if (NCOL(by) != 1) {
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

  names(by.components) <- {
    if (is_not_null(colnames(by))) colnames(by)
    else by.name
  }

  by.factor <- {
    if (is_null(by)) factor(rep.int(1L, n), levels = 1L)
    else factor(by.components[[1]], levels = sort(unique(by.components[[1]])),
                labels = paste(names(by.components), "=", sort(unique(by.components[[1]]))))
  }

  if (treat.type != "continuous" &&
      any(vapply(levels(by.factor), function(x) nunique(treat) != nunique(treat[by.factor == x]), logical(1L)))) {
    .err(sprintf("not all the groups formed by `%s` contain all treatment levels%s. Consider coarsening `%s`",
                 by.arg,
                 if (is_not_null(treat.name)) sprintf(" in %s", treat.name) else "",
                 by.arg))
  }

  attr(by.components, "by.factor") <- by.factor

  by.components
}

.make_closer_to_1 <- function(x) {
  if (chk::vld_character_or_factor(x) || all_the_same(x)) {
    return(x)
  }

  if (is_binary(x)) {
    return(as.numeric(x == max(x, na.rm = TRUE)))
  }

  (x - mean_fast(x, TRUE))/sd(x, na.rm = TRUE)
}

.int_poly_f <- function(d, ex = NULL, int = FALSE, poly = 1, center = TRUE, orthogonal_poly = TRUE) {
  #Adds to data frame interactions and polynomial terms
  #d=matrix input
  #ex=names of variables to exclude in interactions and polynomials; a subset of df
  #int=whether to include interactions or not; currently only 2-way are supported
  #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included

  if (!is.matrix(d)) {
    if (!is.numeric(d))
      .err("an error occurred, probably a bug")

    matrix(d, ncol = 1, dimnames = list(NULL, "x"))
  }

  if (is_null(ex)) ex <- rep.int(FALSE, ncol(d))

  binary.vars <- is_binary_col(d)

  if (center && (int || !orthogonal_poly)) {
    d[,!binary.vars] <- center(d[, !binary.vars, drop = FALSE])
  }

  nd <- NCOL(d)

  if (poly == 0 || nd == 0L) {
    poly_terms <- poly_co.names <- list()
  }
  else if (poly == 1) {
    poly_terms <- list(d)
    poly_co.names <- list(colnames(d))
  }
  else {
    poly_terms <- poly_co.names <- make_list(nd)

    for (i in seq_col(d)) {
      if (ex[i] || binary.vars[i]) {
        poly_terms[[i]] <- d[, i]
        poly_co.names[[i]] <- colnames(d)[i]
      }
      else {
        poly_terms[[i]] <- poly(d[, i], degree = poly, raw = !orthogonal_poly, simple = TRUE)
        poly_co.names[[i]] <- sprintf("%s%s%s",
                                      if (orthogonal_poly) "orth_" else "",
                                      colnames(d)[i],
                                      num_to_superscript(seq_len(poly)))
      }
    }
  }

  if (int && nd > 1L) {
    int_terms <- int_co.names <- make_list(1)
    ints_to_make <- utils::combn(colnames(d)[!ex], 2, simplify = FALSE)

    if (is_not_null(ints_to_make)) {
      int_terms[[1]] <- do.call("cbind", lapply(ints_to_make, function(i) d[,i[1]] * d[,i[2]]))

      int_co.names[[1]] <- vapply(ints_to_make, paste, character(1L), collapse = " * ")
    }
  }
  else {
    int_terms <- int_co.names <- list()
  }

  if (is_null(poly_terms) && is_null(int_terms)) {
    return(matrix(ncol = 0L, nrow = nrow(d), dimnames = list(rownames(d), NULL)))
  }

  out <- do.call("cbind", c(poly_terms, int_terms))
  out_co.names <- c(unlist(poly_co.names), unlist(int_co.names))

  colnames(out) <- out_co.names

  #Remove single values
  if (is_not_null(out)) {
    single_value <- apply(out, 2, all_the_same)
    out <- out[, !single_value, drop = FALSE]
  }

  out
}

.quantile_f <- function(d, qu = NULL, s.weights = NULL, focal = NULL, treat = NULL, const = 2000) {
  # Creates new variables for use in balance quantiles. `qu` is a list of quantiles for each
  # continuous variable in `d`, and returns a matrix with a column for each requested quantile
  # of each variable, taking on 0 for values less than the quantile, .5 for values at the quantile,
  # and 1 for values greater than the quantile. The mean of each variable is equal to the quantile.

  if (is_null(qu)) {
    return(matrix(ncol = 0L, nrow = nrow(d), dimnames = list(rownames(d), NULL)))
  }

  vld_qu <- function(x) {
    is.numeric(x) && all(x >= 0) && all(x <= 1)
  }

  binary.vars <- is_binary_col(d)

  if (length(qu) == 1L && is_null(names(qu))) {
    qu <- setNames(qu[rep.int(1L, sum(!binary.vars))],
                   colnames(d)[!binary.vars])
  }

  if (!all(names(qu) %in% colnames(d)[!binary.vars])) {
    .err("all names of `quantile` must refer to continuous covariates")
  }

  for (i in qu) {
    if (!vld_qu(i)) {
      .err("`quantile` must be a number between 0 and 1 or a named list thereof")
    }
  }

  if (is_not_null(focal) && is_not_null(s.weights)) {
    s.weights <- s.weights[treat == focal]
  }

  do.call("cbind", lapply(names(qu), function(i) {
    target <- if (is_null(focal)) d[,i] else d[treat == focal, i]
    out <- do.call("cbind", lapply(qu[[i]], function(q) {
      plogis(const * (d[,i] - w.quantile(target, q, s.weights)))
    }))

    colnames(out) <- paste0(i, "_", qu[[i]])
    out
  }))
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

get.s.d.denom.cont.weightit <- function(s.d.denom = NULL) {
  s.d.denom.specified <- is_not_null(s.d.denom)

  if (!s.d.denom.specified) {
    return("all")
  }

  allowable.s.d.denoms <- c("all", "weighted")

  try.s.d.denom <- try(match_arg(s.d.denom, allowable.s.d.denoms), silent = TRUE)

  if (!null_or_error(try.s.d.denom)) {
    return(try.s.d.denom)
  }

  "all"
}

.check_estimated_weights <- function(w, treat, treat.type, s.weights) {

  tw <- w * s.weights

  extreme.warn <- FALSE
  if (all_the_same(w)) {
    .wrn(sprintf("all weights are %s, possibly indicating an estimation failure", w[1]))
  }
  else if (treat.type == "continuous") {
    w.cv <- sd(tw, na.rm = TRUE)/mean(tw, na.rm = TRUE)
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
      else if (!extreme.warn && sum(is.finite(tw[ti])) > 1) {
        w.cv <- sd(tw[ti], na.rm = TRUE)/mean(tw[ti], na.rm = TRUE)
        if (!is.finite(w.cv) || w.cv > 4) extreme.warn <- TRUE
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
    .wrn("some weights are negative; these cannot be used in most model fitting functions, including `(g)lm_weightit()`")
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
    ps_mat <- ps_mat[,c(focal, setdiff(colnames(ps_mat), focal))]
  }

  ps_sub <- sub_mat <- ps_mat * 0

  for (i in colnames(ps_mat)) {
    if (estimand == "ATE") {
      sub <- as.integer(findInterval(ps_mat[, as.character(i)],
                                     quantile(ps_mat[, as.character(i)],
                                              seq(0, 1, length.out = subclass + 1)),
                                     all.inside = TRUE))
    }
    else if (estimand == "ATT") {
      if (i != focal) ps_mat[, as.character(i)] <- 1 - ps_mat[, as.character(i)]
      sub <- as.integer(findInterval(ps_mat[, as.character(i)],
                                     quantile(ps_mat[treat == focal, as.character(i)],
                                              seq(0, 1, length.out = subclass + 1)),
                                     all.inside = TRUE))
    }

    sub_tab <- table(treat, sub)

    if (any(sub_tab == 0)) {
      # .err("Too many subclasses were requested")
      sub <- .subclass_scoot(sub, treat, ps_mat[,i])
      sub_tab <- table(treat, sub)
    }

    sub <- as.character(sub)

    sub_totals <- colSums(sub_tab)
    sub_ps <- setNames(sub_tab[as.character(i), ] / sub_totals,
                       colnames(sub_tab))

    ps_sub[,i] <- sub_ps[sub]
    sub_mat[,i] <- sub

    if (ncol(ps_sub) == 2L) {
      ps_sub[,colnames(ps_sub) != i] <- 1 - ps_sub[,i]
      sub_mat[,colnames(sub_mat) != i] <- sub
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

  sub <- as.integer(findInterval(ps,
                                 quantile(switch(estimand,
                                                 "ATE" = ps,
                                                 "ATT" = ps[treat == 1],
                                                 "ATC" = ps[treat == 0]),
                                          seq(0, 1, length.out = subclass + 1)),
                                 all.inside = TRUE))

  max_sub <- max(sub)
  sub_tab1 <- tabulate(sub[treat == 1], max_sub)
  sub_tab0 <- tabulate(sub[treat == 0], max_sub)

  if (any(sub_tab1 == 0) || any(sub_tab0 == 0)) {
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

.subclass_scoot <- function(sub, treat, x, min.n = 1) {
  #Reassigns subclasses so there are no empty subclasses
  #for each treatment group. min.n is the smallest a
  #subclass is allowed to be.
  treat <- as.character(treat)
  unique.treat <- unique(treat, nmax = 2)

  names(x) <- seq_along(x)
  names(sub) <- seq_along(sub)
  original.order <- names(x)

  nsub <- nunique(sub)

  #Turn subs into a contiguous sequence
  sub <- setNames(setNames(seq_len(nsub), sort(unique(sub)))[as.character(sub)],
                  original.order)

  if (any(table(treat) < nsub * min.n)) {
    .err("too many subclasses were requested")
  }

  for (t in unique.treat) {
    if (length(x[treat == t]) == nsub) {
      sub[treat == t] <- seq_len(nsub)
    }
  }

  sub_tab <- table(treat, sub)

  if (!any(sub_tab == 0)) {
    return(sub)
  }

  .soft_thresh <- function(x, minus = 1) {
    x <- x - minus
    x[x < 0] <- 0
    x
  }

  for (t in unique.treat) {
    for (n in seq_len(min.n)) {
      while (any(sub_tab[t,] == 0)) {
        first_0 <- which(sub_tab[t,] == 0)[1]

        if (first_0 == nsub ||
            (first_0 != 1 &&
             sum(.soft_thresh(sub_tab[t, seq(1, first_0 - 1)]) / abs(first_0 - seq(1, first_0 - 1))) >=
             sum(.soft_thresh(sub_tab[t, seq(first_0 + 1, nsub)]) / abs(first_0 - seq(first_0 + 1, nsub))))) {
          #If there are more and closer nonzero subs to the left...
          first_non0_to_left <- max(seq(1, first_0 - 1)[sub_tab[t, seq(1, first_0 - 1)] > 0])

          name_to_move <- names(sub)[which(x == max(x[treat == t & sub == first_non0_to_left]) &
                                             treat == t & sub == first_non0_to_left)[1]]

          sub[name_to_move] <- first_0
          sub_tab[t, first_0] <- 1L
          sub_tab[t, first_non0_to_left] <- sub_tab[t, first_non0_to_left] - 1L

        }
        else {
          #If there are more and closer nonzero subs to the right...
          first_non0_to_right <- min(seq(first_0 + 1, nsub)[sub_tab[t, seq(first_0 + 1, nsub)] > 0])
          name_to_move <- names(sub)[which(x == min(x[treat == t & sub == first_non0_to_right]) &
                                             treat == t & sub == first_non0_to_right)[1]]
          sub[name_to_move] <- first_0
          sub_tab[t, first_0] <- 1L
          sub_tab[t, first_non0_to_right] <- sub_tab[t, first_non0_to_right] - 1L
        }
      }

      sub_tab[t,] <- sub_tab[t,] - 1
    }
  }

  #Unsort
  sub[names(sub)]
}

stabilize_w <- function(weights, treat) {
  if (is.factor(treat)) t.levels <- levels(treat)
  else t.levels <- unique(treat)

  w.names <- names(weights)
  tab <- setNames(vapply(t.levels, function(x) mean_fast(treat == x), numeric(1L)), t.levels)

  setNames(weights * tab[as.character(treat)], w.names)
}

.get_dens_fun <- function(use.kernel = FALSE, bw = NULL, adjust = NULL, kernel = NULL,
                          n = NULL, treat = NULL, density = NULL, weights = NULL) {
  if (is_null(n)) n <- 10 * length(treat)
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

    densfun <- function(p) {
      d <- stats::density(p, n = n,
                          weights = weights/sum(weights), give.Rkern = FALSE,
                          bw = bw, adjust = adjust, kernel = kernel)
      out <- with(d, approxfun(x = x, y = y))(p)
      attr(out, "density") <- d
      out
    }
  }
  else {
    if (is_null(density)) .density <- function(x) dnorm(x)
    else if (is.function(density)) .density <- function(x) density(x)
    else if (identical(density, "dlaplace")) {
      .density <- function(x) {
        mu <- 0
        b <- 1
        exp(-abs(x - mu)/b)/(2 * b)
      }
    }
    else if (is.character(density) && length(density) == 1L) {
      splitdens <- strsplit(density, "_", fixed = TRUE)[[1]]

      if (is_null(splitdens1 <- get0(splitdens[1], mode = "function", envir = parent.frame()))) {
        .err(sprintf("%s is not an appropriate argument to `density` because %s is not an available function",
                     density, splitdens[1]))
      }

      if (length(splitdens) > 1L && !can_str2num(splitdens[-1])) {
        .err(sprintf("%s is not an appropriate argument to `density` because %s cannot be coerced to numeric",
                     density, word_list(splitdens[-1], and.or = "or", quotes = TRUE)))
      }

      .density <- function(x) {
        tryCatch(do.call(splitdens1, c(list(x), as.list(str2num(splitdens[-1])))),
                 error = function(e) {
                   .err(sprintf("Error in applying density:\n  %s", conditionMessage(e)), tidy = FALSE)
                 })
      }

    }
    else {
      .err("the argument to `density` cannot be evaluated as a density function")
    }

    densfun <- function(p) {
      # sd <- sd(p)
      # sd <- sqrt(col.w.v(p, s.weights))
      dens <- .density(p)
      if (is_null(dens) || !is.numeric(dens) || anyNA(dens)) {
        .err("there was a problem with the output of `density`. Try another density function or leave it blank to use the Gaussian density")
      }
      if (any(dens <= 0)) {
        .err("the input to density may not accept the full range of standardized treatment values or residuals")
      }

      x <- seq.int(min(p) - 3 * adjust * bw.nrd0(p),
                   max(p) + 3 * adjust * bw.nrd0(p),
                   length.out = n)
      attr(dens, "density") <- data.frame(x = x,
                                          y = .density(x))
      dens
    }
  }

  densfun
}

.get_w_from_ps_internal_bin <- function(ps, treat, estimand = "ATE",
                                        subclass = NULL, stabilize = FALSE) {

  estimand <- toupper(estimand)
  w <- rep.int(1, length(treat))

  #Assume treat is binary
  if (is_not_null(subclass)) {
    #Get MMW subclass propensity scores
    ps <- .subclass_ps_bin(ps, treat, estimand, subclass)
  }

  i0 <- which(treat == 0)
  i1 <- which(treat == 1)

  if (estimand == "ATE") {
    w[i0] <- 1 / (1 - ps[i0])
    w[i1] <- 1 / ps[i1]
  }
  else if (estimand == "ATT") {
    w[i0] <- .p2o(ps[i0])
  }
  else if (estimand == "ATC") {
    w[i1] <- .p2o(1 - ps[i1])
  }
  else if (estimand == "ATO") {
    w[i0] <- ps[i0]
    w[i1] <- (1 - ps[i1])
  }
  else if (estimand == "ATM") {
    w[i0][ps[i0] < .5] <- .p2o(ps[i0][ps[i0] < .5])
    w[i1][ps[i1] > .5] <- .p2o(1 - ps[i1][ps[i1] > .5])
  }
  else if (estimand == "ATOS") {
    w[i0] <- 1 / (1 - ps[i0])
    w[i1] <- 1 / ps[i1]

    ps.sorted <- sort(c(ps, 1 - ps))
    q <- ps * (1 - ps)
    alpha.opt <- 0
    for (i in seq_len(sum(ps < .5))) {
      if (i == 1 || !check_if_zero(ps.sorted[i] - ps.sorted[i-1])) {
        alpha <- ps.sorted[i]
        a <- alpha * (1 - alpha)
        if (1/a <= 2*sum(1/q[q >= a])/sum(q >= a)) {
          alpha.opt <- alpha
          break
        }
      }
    }
    w[!between(ps, c(alpha.opt, 1 - alpha.opt))] <- 0
  }

  names(w) <- if_null_then(names(treat), NULL)

  if (stabilize) w <- stabilize_w(w, treat)

  w
}

.get_w_from_ps_internal_multi <- function(ps, treat, estimand = "ATE", focal = NULL,
                                          subclass = NULL, stabilize = FALSE) {

  estimand <- toupper(estimand)
  w <- rep.int(0, length(treat))

  ps_mat <- ps

  if (is_not_null(subclass)) {
    #Get MMW subclass propensity scores
    ps_mat <- .subclass_ps_multi(ps_mat, treat, estimand, focal, subclass)
  }

  for (i in colnames(ps_mat)) {
    w[treat == i] <- 1/ps_mat[treat == i, i]
  }

  if (estimand == "ATE") {
    # w <- w
  }
  else if (estimand %in% c("ATT", "ATC")) {
    w <- w * ps_mat[, as.character(focal)]
  }
  else if (estimand == "ATO") {
    w <- w / rowSums(1 / ps_mat) #Li & Li (2019)
  }
  else if (estimand == "ATM") {
    w <- w * do.call("pmin", lapply(seq_col(ps_mat), function(x) ps_mat[,x]), quote = TRUE)
  }
  else if (estimand == "ATOS") {
    #Crump et al. (2009)
    ps.sorted <- sort(c(ps_mat[,2], 1 - ps_mat[,2]))
    q <- ps_mat[,2]*(1-ps_mat[,2])
    alpha.opt <- 0
    for (i in 1:sum(ps_mat[,2] < .5)) {
      if (i == 1 || !check_if_zero(ps.sorted[i] - ps.sorted[i-1])) {
        alpha <- ps.sorted[i]
        a <- alpha*(1-alpha)
        if (1/a <= 2*sum(1/q[q >= a])/sum(q >= a)) {
          alpha.opt <- alpha
          break
        }
      }
    }
    w[!between(ps_mat[,2], c(alpha.opt, 1 - alpha.opt))] <- 0
  }
  else return(numeric(0))

  if (stabilize) w <- stabilize_w(w, treat)

  names(w) <- if_null_then(rownames(ps_mat), names(treat), NULL)

  w
}

.get_w_from_ps_internal_array <- function(ps, treat, estimand = "ATE", focal = NULL,
                                          subclass = NULL, stabilize = FALSE) {
  #Batch turn PS into weights; primarily for output of predict.gbm
  # Assumes a (0,1) treatment if binary
  if (is_null(dim(ps))) {
    ps <- matrix(ps, ncol = 1)
  }

  eps <- 1e-8

  if (length(dim(ps)) == 2) {
    #Binary treatment, vector ps

    w <- ps
    w[] <- 0

    if (is_not_null(subclass)) {
      #Get MMW subclass propensity scores
      for (p in seq_col(ps)) {
        ps[,p] <- .subclass_ps_bin(ps[,p], treat, estimand, subclass)
      }
    }

    t1 <- which(treat == 1)
    t0 <- which(treat == 0)

    if (estimand == "ATE") {
      ps[t1,][ps[t1,] < eps] <- eps
      ps[t0,][ps[t0,] > 1 - eps] <- 1 - eps

      w[t1,] <- 1 / ps[t1,]
      w[t0,] <- 1 / (1 - ps[t0,])
    }
    else if (estimand == "ATT") {
      ps[t0,][ps[t0,] > 1 - eps] <- 1 - eps

      w[t1,] <- 1
      w[t0,] <- .p2o(ps[t0,])
    }
    else if (estimand == "ATC") {
      ps[t1,][ps[t1,] < eps] <- eps

      w[t1,] <- .p2o(1 - ps[t1,])
      w[t0,] <- 1
    }
    else if (estimand == "ATO") {
      w[t1,] <- 1 - ps[t1,]
      w[t0,] <- ps[t0,]
    }
    else if (estimand == "ATM") {
      pslt.5 <- ps < .5
      w[t1,][pslt.5[t1,]] <- 1
      w[t1,][!pslt.5[t1,]] <- .p2o(1 - ps[t1,][!pslt.5[t1,]])

      w[t0,][pslt.5[t0,]] <- .p2o(ps[t0,][pslt.5[t0,]])
      w[t0,][!pslt.5[t0,]] <- 1
    }

    if (stabilize) {
      w[t1] <- w[t1] * length(t1) / length(treat)
      w[t0] <- w[t0] * length(t0) / length(treat)
    }
  }
  else if (length(dim(ps)) == 3) {
    #Multi-category treatment, matrix PS

    if (is_not_null(subclass)) {
      #Get MMW subclass propensity scores
      for (p in seq_len(last(dim(ps))))
        ps[,,p] <- .subclass_ps_multi(ps[,,p], treat, estimand, focal, subclass)
    }

    ps <- squish(ps, eps)

    w <- matrix(0, ncol = dim(ps)[3], nrow = dim(ps)[1])
    t.levs <- unique(treat)

    for (i in t.levs) {
      w[treat == i,] <- 1 / ps[treat == i, as.character(i),]
    }

    if (estimand == "ATE") {
      #Do nothing
    }
    else if (estimand %in% c("ATT", "ATC")) {
      not_focal <- which(treat != focal)
      w[-not_focal,] <- 1
      w[not_focal,] <- w[not_focal,] * ps[not_focal, as.character(focal),]
    }
    else if (estimand == "ATO") {
      w <- w / colSums(aperm(1/ps, c(2, 1, 3)))
    }
    else if (estimand == "ATM") {
      treat <- as.integer(treat)

      for (p in seq_len(dim(ps)[3])) {
        # w[,p] <- w[,p] * do.call("pmin", lapply(seq_len(dim(ps)[2]), function(i) ps[,i,p]))

        ps_p <- ps[,,p]
        min_ind <- max.col(-ps_p, ties.method = "first")
        no_match <- which(ps_p[cbind(seq_along(treat), treat)] != ps_p[cbind(seq_along(treat), min_ind)])

        if (length(no_match) < length(treat)) {
          w[-no_match, p] <- 1
        }

        if (length(no_match) > 0) {
          w[no_match, p] <- w[no_match, p] * ps_p[cbind(no_match, min_ind[no_match])]
        }
      }
    }

    if (stabilize) {
      for (i in t.levs) {
        w[treat == i,] <- mean_fast(treat == i) * w[treat == i,]
      }
    }
  }
  else {
    .err("don't know how to process more than 3 dims (likely a bug)")
  }

  w
}

plot_density <- function(d.n, d.d) {
  d.d <- cbind(as.data.frame(d.d[c("x", "y")]), dens = "Denominator Density", stringsAsfactors = FALSE)
  d.n <- cbind(as.data.frame(d.n[c("x", "y")]), dens = "Numerator Density", stringsAsfactors = FALSE)
  d.all <- rbind(d.d, d.n)
  d.all$dens <- factor(d.all$dens, levels = c("Numerator Density", "Denominator Density"))
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
  w <- w/mean_fast(w)
  mean(w*log(w))
}

replace_na_with <- function(covs, with = "median") {
  if (is.na(with) || !anyNA(covs)) return(covs)

  if (is.character(with)) {
    .with <- match.fun(with)
    for (i in colnames(covs)[anyNA_col(covs)]) {
      if (all(is.na(covs[,i]))) covs <- covs[, colnames(covs) != i, drop = FALSE]
      else covs[is.na(covs[,i]), i] <- .with(covs[, i], na.rm = TRUE)
    }
  }
  else {
    covs[is.na(covs)] <- with
  }

  covs
}

add_missing_indicators <- function(covs, replace_with = "median") {
  covs_w_missing <- which(anyNA_col(covs))
  if (is_null(covs_w_missing)) return(covs)

  missing_ind <- apply(covs[, covs_w_missing, drop = FALSE], 2, function(x) as.numeric(is.na(x)))

  colnames(missing_ind) <- paste0(colnames(missing_ind), ":<NA>")
  covs <- cbind(covs, missing_ind)

  if (is_null(replace_with) || is.na(replace_with)) return(covs)

  replace_na_with(covs, replace_with)
}

verbosely <- function(expr, verbose = TRUE) {
  if (verbose) return(expr)

  void <- utils::capture.output({
    out <- invisible(expr)
  })

  out
}

#Generalized matrix inverse (port of MASS::ginv)
generalized_inverse <- function(sigma) {
  sigmasvd <- svd(sigma)
  pos <- sigmasvd$d > max(1e-9 * sigmasvd$d[1L], 0)
  sigma_inv <- sigmasvd$v[, pos, drop = FALSE] %*% (sigmasvd$d[pos]^-1 * t(sigmasvd$u[, pos, drop = FALSE]))

  sigma_inv
}

#Compute gradient numerically using centered difference
.gradient <- function(.f, .x, .eps = 1e-8, parm = NULL, ...) {
  if (is_null(parm)) {
    parm <- seq_along(.x)
  }

  .x0 <- .x

  .eps <- pmax(abs(.x) * .eps, .eps)

  for (j in parm) {
    # forward
    .x[j] <- .x0[j] + .eps[j]/2

    # recalculate model function value
    f_new_forward <- .f(.x, ...)

    if (j == 1L) {
      jacob <- matrix(0, nrow = length(f_new_forward), ncol = length(parm),
                      dimnames = list(names(f_new_forward), names(.x)[parm]))
    }

    # backward
    .x[j] <- .x0[j] - .eps[j]/2

    # recalculate model function value
    f_new_backward  <- .f(.x, ...)

    jacob[,j] <- (f_new_forward - f_new_backward) / .eps[j]

    .x[j] <- .x0[j]
  }

  jacob
}

#Convert probability to odds
.p2o <- function(p) {
  p / (1 - p)
}

#Get psi function (individual contributions to gradient) from glm fit
.get_glm_psi <- function(fit) {
  family <- fit$family

  if (!identical(fit$class, "brglmFit") ||
      identical(fit$type, "ML") || identical(fit$type, "correction")) {
    psi <- function(B, X, y, weights, offset = 0) {
      XB <- drop(X %*% B) + offset
      p <- family$linkinv(XB)
      D <- family$mu.eta(XB)
      V <- family$variance(p)

      X * (weights * D * (y - p) / V)
    }
  }
  else {
    br_type <- fit$type

    if (is_null(fit$control[["a"]])) {
      rlang::check_installed("brglm2")
      fit$control[["a"]] <- formals(brglm2::brglmControl)[["a"]]
    }

    br_psi <- function(X, W, D, p, XB, V) {
      DD <- family$d2mu.deta(XB)
      Wt <- W * D^2 / V #"working weight"

      ## Compute hat values
      XWt <- sqrt(Wt) * X
      q <- qr(XWt)
      Q <- qr.Q(q)
      H <- rowSums(Q * Q)

      if (br_type %in% c("AS_mixed", "AS_mean")) {
        AA <- .5 * X * H * DD / D
        return(AA)
      }

      V1 <- family$d1variance(p)

      if (br_type == "MPL_Jeffreys") {
        return(fit$control[["a"]] * X * H * (2 * DD / D - V1 * D / V))
      }

      #br_type == "AS_median"
      R_matrix <- qr.R(q)
      info_unscaled <- crossprod(R_matrix)
      inverse_info_unscaled <- chol2inv(R_matrix)

      b_vector <- vapply(seq_col(X), function(j) {
        inverse_info_unscaled_j <- inverse_info_unscaled[j, ]
        vcov_j <- tcrossprod(inverse_info_unscaled_j) / inverse_info_unscaled_j[j]
        hats_j <- rowSums((X %*% vcov_j) * X) * Wt

        inverse_info_unscaled_j %*% colSums(X * (hats_j * (D * V1 / (6 * V) - DD / (2 * D))))
      }, numeric(1L))

      AA <- .5 * X * H * DD / D
      sweep(AA, 2, info_unscaled %*% b_vector / nrow(X), "+")
    }

    psi <-  function(B, X, y, weights, offset = 0) {
      XB <- drop(X %*% B) + offset
      p <- family$linkinv(XB)
      D <- family$mu.eta(XB)
      V <- family$variance(p)

      .psi <- X * (weights * D * (y - p) / V)

      .psi + br_psi(X, weights, D, p, XB, V)
    }
  }

  psi
}
