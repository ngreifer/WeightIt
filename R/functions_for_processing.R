allowable.methods <- function() {
  c("glm" = "glm", "ps" = "glm",
    "gbm" = "gbm", "gbr" = "gbm",
    "cbps" = "cbps",
    "npcbps" = "npcbps",
    "ebal" = "ebal", "entropy" = "ebal", "ebalance" = "ebal",
    # "ebcw" = "ebcw", "ate" = "ebcw",
    "optweight" = "optweight", "sbw" = "optweight",
    "super" = "super", "superlearner" = "super",
    "bart" = "bart",
    "energy" = "energy")
}
method.to.proper.method <- function(method) {
  method <- tolower(method)
  allowable.methods()[method]
}
check.acceptable.method <- function(method, msm = FALSE, force = FALSE) {
  bad.method <- FALSE

  if (missing(method)) method <- "glm"
  else if (is_null(method) || length(method) > 1) bad.method <- TRUE
  else if (is.character(method)) {
    if (tolower(method) %nin% names(allowable.methods())) bad.method <- TRUE
  }
  else if (!is.function(method)) bad.method <- TRUE

  if (bad.method) {
    if (identical(method, "twang")) {
      .err('"twang" is no longer an acceptable argument to `method`. Please use "gbm" for generalized boosted modeling')
    }

    .err(paste0("`method` must be a string of length 1 containing the name of an acceptable weighting\n\tmethod or a function that produces weights. Allowable methods:\n", paste(add_quotes(unique(allowable.methods())), collapse = ", ")), tidy = FALSE)
  }

  if (msm && !force && is.character(method)) {
    m <- method.to.proper.method(method)
    if (m %in% c("nbcbps", "ebal", "ebcw", "optweight", "energy", "kbal")) {
      .err(sprintf("the use of %s with longitudinal treatments has not been validated. Set `weightit.force = TRUE` to bypass this error",
                   method.to.phrase(m)))
    }
  }
}
check.user.method <- function(method) {
  #Check to make sure it accepts treat and covs
  if (all(c("covs", "treat") %in% names(formals(method)))) {
  }
  # else if (all(c("covs.list", "treat.list") %in% names(formals(method)))) {
  # }
  else {
    .err("the user-provided function to `method` must contain `covs` and `treat` as named parameters")
  }
}
method.to.phrase <- function(method) {

  if (is.function(method)) return("a user-defined method")

  method <- method.to.proper.method(method)
  if (method %in% c("glm")) return("propensity score weighting with GLM")
  if (method %in% c("gbm")) return("propensity score weighting with GBM")
  if (method %in% c("cbps")) return("covariate balancing propensity score weighting")
  if (method %in% c("npcbps")) return("non-parametric covariate balancing propensity score weighting")
  if (method %in% c("ebal")) return("entropy balancing")
  # if (method %in% c("ebcw")) return("empirical balancing calibration weighting")
  if (method %in% c("optweight")) return("targeted stable balancing weights")
  if (method %in% c("super")) return("propensity score weighting with SuperLearner")
  if (method %in% c("bart")) return("propensity score weighting with BART")
  if (method %in% c("energy")) return("energy balancing")
  # if (method %in% c("kbal")) return("kernel balancing")

  "the chosen method of weighting"
}
process.estimand <- function(estimand, method, treat.type) {
  #Allowable estimands
  AE <- list(
    binary = list(  glm = c("ATT", "ATC", "ATE", "ATO", "ATM", "ATOS")
                    , gbm = c("ATT", "ATC", "ATE", "ATO", "ATM")
                    , cbps = c("ATT", "ATC", "ATE")
                    , npcbps = c("ATE")
                    , ebal = c("ATT", "ATC", "ATE")
                    # , ebcw = c("ATT", "ATC", "ATE")
                    , optweight = c("ATT", "ATC", "ATE")
                    , super = c("ATT", "ATC", "ATE", "ATO", "ATM")
                    , energy = c("ATT", "ATC", "ATE")
                    , bart = c("ATT", "ATC", "ATE", "ATO", "ATM")
                    # , kbal = c("ATT", "ATC", "ATE")
    ),
    multinomial = list(  glm = c("ATT", "ATC", "ATE", "ATO", "ATM")
                         , gbm = c("ATT", "ATC", "ATE", "ATO", "ATM")
                         , cbps = c("ATT", "ATC", "ATE")
                         , npcbps = c("ATE")
                         , ebal = c("ATT", "ATC", "ATE")
                         # , ebcw = c("ATT", "ATC", "ATE")
                         , optweight = c("ATT", "ATC", "ATE")
                         , super = c("ATT", "ATC", "ATE", "ATO", "ATM")
                         , energy = c("ATT", "ATC", "ATE")
                         , bart = c("ATT", "ATC", "ATE", "ATO", "ATM")
                         # , kbal = c("ATT", "ATE")
    ))


  if (treat.type == "continuous" || is.function(method)) {
    .chk_null_or(estimand, chk::chk_string)
    return(toupper(estimand))
  }

  chk::chk_string(estimand)
  estimand <- toupper(estimand)

  if (estimand %nin% AE[[treat.type]][[method]]) {
    .err(sprintf("%s is not an allowable estimand for %s with %s treatments. Only %s allowed",
                 add_quotes(estimand), method.to.phrase(method), treat.type,
                 word_list(AE[[treat.type]][[method]], quotes = TRUE, and.or = "and", is.are = TRUE)))
  }

  estimand
}
check.subclass <- function(method, treat.type) {
  #Allowable estimands
  AE <- list(
    binary = list(  glm = TRUE
                    , gbm = TRUE
                    , cbps = TRUE
                    , npcbps = FALSE
                    , ebal = FALSE
                    # , ebcw = FALSE
                    , optweight = FALSE
                    , super = TRUE
                    , energy = FALSE
                    , bart = TRUE
                    # , kbal = FALSE
    ),
    multinomial = list(  glm = TRUE
                         , gbm = TRUE
                         , cbps = FALSE
                         , npcbps = FALSE
                         , ebal = FALSE
                         # , ebcw = FALSE
                         , optweight = FALSE
                         , super = TRUE
                         , energy = FALSE
                         , bart = TRUE
    ))

  if (treat.type != "continuous" && !is.function(method) &&
      !AE[[treat.type]][[method]]) {
    .err(sprintf("subclasses are not compatible with %s with %s treatments",
                 method.to.phrase(method), treat.type))
  }
}
process.ps <- function(ps, data = NULL, treat) {
  if (is_null(ps)) return(NULL)

  if (is.character(ps) && length(ps) == 1L) {
    if (is_null(data)) {
      .err("`ps` was specified as a string but there was no argument to `data`")
    }
    else if (ps %in% names(data)) {
      ps <- data[[ps]]
    }
    else {
      .err("the name supplied to `ps` is not the name of a variable in `data`")
    }
  }
  else if (is.numeric(ps)) {
    if (length(ps) != length(treat)) {
      .err("`ps` must have the same number of units as the treatment")
    }
  }
  else {
    .err("the argument to `ps` must be a vector of propensity scores or the (quoted) name of the variable in `data` that contains propensity scores")
  }

  ps
}
process.focal.and.estimand <- function(focal, estimand, treat, treated = NULL) {
  reported.estimand <- estimand

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  unique.treat <- unique(treat, nmax = switch(treat.type, "binary" = 2, "multinomial" = length(treat)/4))

  #Check focal
  if (is_not_null(focal) && (length(focal) > 1L || focal %nin% unique.treat)) {
    .err("the argument supplied to `focal` must be the name of a level of treatment")
  }

  if (treat.type == "multinomial") {

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
          if (is.factor(treat)) treated <- levels(treat)[2]
          else treated <- unique.treat[unique.treat.bin == 1]

          if (estimand == "ATT") {
            .msg(sprintf("assuming %s the treated level. If not, supply an argument to `focal`",
                         word_list(treated, quotes = !is.numeric(treat), is.are = TRUE)))

          }
          else if (estimand == "ATC") {
            .msg(sprintf("assuming %s the control level. If not, supply an argument to `focal`",
                         word_list(unique.treat[unique.treat %nin% treated], quotes = !is.numeric(treat), is.are = TRUE)))
          }

        }

        if (estimand == "ATT")
          focal <- treated
        else if (estimand == "ATC")
          focal <- unique.treat[unique.treat %nin% treated]
      }
      else {
        if (estimand == "ATT")
          treated <- focal
        else if (estimand == "ATC")
          treated <- unique.treat[unique.treat %nin% focal]
      }
      if (estimand == "ATC") estimand <- "ATT"
    }
    else {
      if (is_null(focal)) {
        if (estimand == "ATT")
          focal <- treated
        else if (estimand == "ATC")
          focal <- unique.treat[unique.treat %nin% treated]
      }
      if (estimand == "ATC") estimand <- "ATT"
    }
  }

  list(focal = as.character(focal),
       estimand = estimand,
       reported.estimand = reported.estimand,
       treated = if (is.factor(treated)) as.character(treated) else treated)
}
process.by <- function(by, data, treat, treat.name = NULL, by.arg = "by") {

  ##Process by
  bad.by <- FALSE
  n <- length(treat)

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  if (missing(by)) {
    bad.by <- TRUE
  }
  else if (is_null(by)) {
    by <- NULL
    by.name <- NULL
  }
  else if (is.character(by) && length(by) == 1 && by %in% names(data)) {
    by.name <- by
    by <- data[[by]]
  }
  else if (length(dim(by)) == 2L && len(by) == n) {
    by.name <- colnames(by)[1]
    by <- drop(by[, 1])
  }
  else if (rlang::is_formula(by, lhs = FALSE)) {
    t.c <- get.covs.and.treat.from.formula(by, data)
    by <- t.c[["reported.covs"]]
    if (NCOL(by) != 1) {
      .err(sprintf("only one variable can be on the right-hand side of the formula for `%s`",
                   by.arg))
    }
    by.name <- colnames(by)
  }
  else bad.by <- TRUE

  if (!bad.by) {
    by.components <- data.frame(by)

    if (is_not_null(colnames(by))) names(by.components) <- colnames(by)
    else names(by.components) <- by.name

    if (is_null(by)) by.factor <- factor(rep(1L, n), levels = 1L)
    else by.factor <- factor(by.components[[1]], levels = unique(by.components[[1]]),
                             labels = paste(names(by.components), "=", unique(by.components[[1]])))
    # by.vars <- acceptable.bys[vapply(acceptable.bys, function(x) equivalent.factors(by, data[[x]]), logical(1L))]
  }
  else {
    .err(sprintf("`%s` must be a string containing the name of the variable in data for which weighting is to occur within strata or a one-sided formula with the stratifying variable on the right-hand side",
                 by.arg))
  }

  if (treat.type != "continuous" && any(vapply(levels(by.factor), function(x) nunique(treat) != nunique(treat[by.factor == x]), logical(1L)))) {
    .err(sprintf("Not all the groups formed by `%s` contain all treatment levels%s. Consider coarsening `%s`",
                 by.arg, if (is_not_null(treat.name)) paste(" in", treat.name) else "", by.arg))
  }
  attr(by.components, "by.factor") <- by.factor

  by.components
}
process.moments.int <- function(moments, int, method) {

  if (is.function(method) || method %in% c("npcbps", "ebal", "ebcw", "optweight", "energy")) {
    chk::chk_flag(int)

    if (is_not_null(moments)) {
      if (length(moments) != 1 || !is.numeric(moments) ||
          !check_if_zero(moments - round(moments))) {
        chk::chk_whole_number(moments)
        if (method == "energy") {
          chk::chk_gte(moments, 0)
        }
        else if (method %in% c("npcbps", "ebal", "ebcw", "optweight")) {
          chk::chk_gt(moments, 0)
        }
        moments <- as.integer(moments)
      }
    }
    else {
      moments <- {
        if (!is.function(method) && method == "energy" && !int) 0L
        else 1L
      }
    }
  }
  else if (is_not_null(moments) && any(mi0 <- c(as.integer(moments) != 1L, int))) {
    .wrn(sprintf("%s not compatible with %s. Ignoring %s",
                 word_list(c("moments", "int")[mi0], and.or = "and", is.are = TRUE, quotes = "`"),
                 method.to.phrase(method),
                 word_list(c("moments", "int")[mi0], and.or = "and", quotes = "`")))
    moments <- NULL
    int <- FALSE
  }
  moments <- as.integer(moments)

  list(moments = moments, int = int)
}
process.MSM.method <- function(is.MSM.method, method) {
  methods.with.MSM <- c("optweight")
  if (is.function(method)) {
    if (isTRUE(is.MSM.method)) .err("currently, only user-defined methods that work with `is.MSM.method = FALSE` are allowed")
    is.MSM.method <- FALSE
  }
  else if (method %in% methods.with.MSM) {
    if (is_null(is.MSM.method)) is.MSM.method <- TRUE
    else if (!isTRUE(is.MSM.method)) {
      .msg(paste0("%s can be used with a single model when multiple time points are present.\nUsing a seperate model for each time point. To use a single model, set `is.MSM.method` to `TRUE`",
                  method.to.phrase(method)))
    }
  }
  else {
    if (isTRUE(is.MSM.method)) {
      .wrn(sprintf("%s cannot be used with a single model when multiple time points are present.\nUsing a seperate model for each time point",
                   method.to.phrase(method)))
    }
    is.MSM.method <- FALSE
  }

  is.MSM.method
}
process.missing <- function(missing, method, treat.type) {
  #Allowable estimands
  AE <- list(binary = list(glm = c("ind", "saem")
                           , gbm = c("ind", "surr")
                           , cbps = c("ind")
                           , npcbps = c("ind")
                           , ebal = c("ind")
                           # , ebcw = c("ind")
                           , optweight = c("ind")
                           , super = c("ind")
                           , bart = c("ind")
                           , energy = c("ind")
                           # , kbal = c("ind")
  ),
  multinomial = list(glm = c("ind")
                     , gbm = c("ind", "surr")
                     , cbps = c("ind")
                     , npcbps = c("ind")
                     , ebal = c("ind")
                     # , ebcw = c("ind")
                     , optweight = c("ind")
                     , super = c("ind")
                     , bart = c("ind")
                     , energy = c("ind")
                     # , kbal = c("ind")
  ),
  continuous = list(glm = c("ind", "saem")
                    , gbm = c("ind", "surr")
                    , cbps = c("ind")
                    , npcbps = c("ind")
                    , ebal = c("ind")
                    # , ebcw = c("ind")
                    , optweight = c("ind")
                    , super = c("ind")
                    , bart = c("ind")
                    , energy = c("ind")
                    # , kbal = c("ind")
  ))

  allowable.missings <- AE[[treat.type]][[method]]

  if (is_null(missing)) {
    .wrn(sprintf("missing values are present in the covariates. See `?WeightIt::method_%s` for information on how these are handled",
                 method))
    return(allowable.missings[1])
  }

  chk::chk_string(missing)
  if (!missing %pin% allowable.missings) {
    .err(sprintf("only %s allowed for `missing` with `method = %s` for %s treatments",
                 word_list(allowable.missings, quotes = 2, is.are = TRUE),
                 add_quotes(method),
                 treat.type))
    return(allowable.missings[1])
  }

  allowable.missings[pmatch(missing, allowable.missings)]
}
make.closer.to.1 <- function(x) {
  if (chk::vld_character_or_factor(x) || all_the_same(x)) {
    return(x)
  }

  if (is_binary(x)) {
    return(as.numeric(x == x[!is.na(x)][1]))
  }

  (x - mean_fast(x, TRUE))/sd(x, na.rm = TRUE)
}
int.poly.f <- function(d, ex = NULL, int = FALSE, poly = 1, center = TRUE, orthogonal_poly = TRUE) {
  #Adds to data frame interactions and polynomial terms
  #d=matrix input
  #ex=names of variables to exclude in interactions and polynomials; a subset of df
  #int=whether to include interactions or not; currently only 2-way are supported
  #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included

  if (is_null(ex)) ex <- rep(FALSE, ncol(d))

  binary.vars <- is_binary_col(d)

  if (center && (int || !orthogonal_poly)) {
    d[,!binary.vars] <- center(d[, !binary.vars, drop = FALSE])
  }
  nd <- NCOL(d)

  if (poly > 1) {
    make.poly <- which(!binary.vars & !ex)
    npol <- length(make.poly)
    poly_terms <- poly_co.names <- make_list(npol)
    if (npol > 0) {
      for (i in seq_along(make.poly)) {
        poly_terms[[i]] <- poly(d[, make.poly[i]], degree = poly, raw = !orthogonal_poly, simple = TRUE)[,-1, drop = FALSE]
        poly_co.names[[i]] <- paste0(if (orthogonal_poly) "orth_", colnames(d)[make.poly[i]], num_to_superscript(2:poly))
      }
    }
  }
  else poly_terms <- poly_co.names <- list()

  if (int && nd > 1) {
    int_terms <- int_co.names <- make_list(1)
    ints_to_make <- combn(colnames(d)[!ex], 2, simplify = FALSE)

    if (is_not_null(ints_to_make)) {
      int_terms[[1]] <- do.call("cbind", lapply(ints_to_make, function(i) d[,i[1]]*d[,i[2]]))

      int_co.names[[1]] <- vapply(ints_to_make, paste, character(1L), collapse = " * ")
    }
  }
  else int_terms <- int_co.names <- list()

  if (is_not_null(poly_terms) || is_not_null(int_terms)) {
    out <- do.call("cbind", c(poly_terms, int_terms))
    out_co.names <- c(do.call("c", poly_co.names), do.call("c", int_co.names))

    colnames(out) <- unlist(out_co.names)

    #Remove single values
    if (is_not_null(out)) {
      single_value <- apply(out, 2, all_the_same)
      out <- out[, !single_value, drop = FALSE]
    }
  }
  else {
    out <- matrix(ncol = 0, nrow = nrow(d), dimnames = list(rownames(d), NULL))
  }

  return(out)
}
get.s.d.denom.weightit <- function(s.d.denom = NULL, estimand = NULL, weights = NULL, treat = NULL, focal = NULL) {
  check.estimand <- check.weights <- check.focal <- FALSE
  s.d.denom.specified <- is_not_null(s.d.denom)
  estimand.specified <- is_not_null(estimand)
  if (!is.factor(treat)) treat <- factor(treat)

  if (s.d.denom.specified) {
    allowable.s.d.denoms <- c("treated", "control", "pooled", "all", "weighted", "hedges")
    try.s.d.denom <- tryCatch(match_arg(s.d.denom, allowable.s.d.denoms),
                              error = function(cond) NA_character_)
    if (anyNA(try.s.d.denom)) {
      check.estimand <- TRUE
    }
    else {
      s.d.denom <- try.s.d.denom
    }
  }
  else {
    check.estimand <- TRUE
  }

  if (check.estimand) {
    if (estimand.specified) {
      allowable.estimands <- c("ATT", "ATC", "ATE", "ATO", "ATM")
      try.estimand <- tryCatch(match_arg(toupper(estimand), allowable.estimands),
                               error = function(cond) NA_character_)
      if (anyNA(try.estimand) || try.estimand %in% c("ATC", "ATT")) {
        check.focal <- TRUE
      }
      else {
        s.d.denom <- vapply(try.estimand, switch, FUN.VALUE = character(1L),
                            ATO = "weighted", ATM = "weighted", "pooled")
      }
    }
    else {
      check.focal <- TRUE
    }
  }
  if (check.focal) {
    if (is_not_null(focal)) {
      s.d.denom <- focal
    }
    else check.weights <- TRUE
  }
  if (check.weights) {
    if (is_null(weights)) {
      s.d.denom <- "pooled"
    }
    else {
      for (tv in levels(treat)) {
        if (all_the_same(weights[treat == tv]) &&
            !all_the_same(weights[treat != tv])) {
          s.d.denom <- tv
        }
        else if (tv == last(levels(treat))) {
          s.d.denom <- "pooled"
        }
      }
    }
  }

  s.d.denom
}
get.s.d.denom.cont.weightit <- function(s.d.denom = NULL) {
  s.d.denom.specified <- is_not_null(s.d.denom)

  if (!s.d.denom.specified) {
    return("all")
  }

  allowable.s.d.denoms <- c("all", "weighted")

  try.s.d.denom <- tryCatch(match_arg(s.d.denom, allowable.s.d.denoms),
                            error = function(cond) NA_character_)
  if (anyNA(try.s.d.denom)) {
    return("all")
  }

  try.s.d.denom
}
check_estimated_weights <- function(w, treat, treat.type, s.weights) {

  tw <- w * s.weights

  extreme.warn <- FALSE
  if (treat.type == "continuous") {
    if (all_the_same(w)) {
      .wrn(sprintf("all weights are %s, possibly indicating an estimation failure", w[1]))
    }
    else {
      w.cv <- sd(tw, na.rm = TRUE)/mean(tw, na.rm = TRUE)
      if (!is.finite(w.cv) || w.cv > 4) extreme.warn <- TRUE
    }
  }
  else {
    if (all_the_same(w)) {
      .wrn(sprintf("all weights are %s, possibly indicating an estimation failure", w[1]))
    }
    else {
      t.levels <- unique(treat)
      bad.treat.groups <- setNames(rep(FALSE, length(t.levels)), t.levels)
      for (i in t.levels) {
        ti <- which(treat == i)
        if (all(is.na(w[ti])) || all(check_if_zero(w[ti]))) bad.treat.groups[as.character(i)] <- TRUE
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
  }

  if (extreme.warn) {
    .wrn("some extreme weights were generated. Examine them with `summary()` and maybe trim them with `trim()`")
  }
}
ps_to_ps_mat <- function(ps, treat, assumed.treated = NULL, treat.type = NULL, treated = NULL, estimand = NULL) {
  if (is_(ps, c("matrix", "data.frame"))) {
    ps.names <- rownames(ps)
  }
  else if (is_(ps, "numeric")) {
    ps.names <- names(ps)
    ps <- matrix(ps, ncol = 1)
  }

  if (is.factor(treat)) t.levels <- levels(treat)
  else t.levels <- unique(treat, nmax = length(treat)/4)

  if (treat.type == "binary") {
    if ((!is.matrix(ps) || !is.numeric(ps)) &&
        (!is.data.frame(ps) || !all(vapply(ps, is.numeric, logical(1L))))) {
      .err("`ps` must be a matrix, data frame, or vector of propensity scores")
    }

    if (ncol(ps) == 1) {
      if (can_str2num(treat) &&
          all(check_if_zero(binarize(treat) - str2num(treat)))) {
        treated.level <- 1
      }
      else if (is_not_null(treated)) {
        if (!treated %in% treat) {
          .err("the argument to `treated` must be a value in `treat`")
        }
        treated.level <- treated
      }
      else if (is_not_null(assumed.treated)) {
        treated.level <- assumed.treated
      }
      else if (is_not_null(colnames(ps)) && colnames(ps) %in% as.character(t.levels)) {
        treated.level <- colnames(ps)
      }
      else {
        .err("if the treatment has two non-0/1 levels and `ps` is a vector or has only one column, an argument to `treated` must be supplied")
      }

      t.levels <- c(treated.level, t.levels[t.levels != treated.level])
      ps <- matrix(c(ps[, 1], 1 - ps[, 1]), ncol = 2, dimnames = list(ps.names, as.character(t.levels)))
    }
    else if (ncol(ps) == 2) {
      if (!all(as.character(t.levels) %in% colnames(ps))) {
        .err("if `ps` has two columns, they must be named with the treatment levels")
      }
    }
    else {
      .err("`ps` cannot have more than two columns if the treatment is binary")
    }

  }
  else if (treat.type == "multinomial") {
    if ((!is.matrix(ps) || !is.numeric(ps)) &&
        (!is.data.frame(ps) || !all(vapply(ps, is.numeric, logical(1L))))) {
      .err("`ps` must be a matrix or data frame of propensity scores")
    }

    if (ncol(ps) == 1) {
      if (toupper(estimand) != "ATE") {
        .err("with multinomial treatments, `ps` can be a vector or have only one column only if the estimand is the ATE")
      }

      ps <- matrix(rep(ps, nunique(treat)), nrow = length(treat), dimnames = list(ps.names, t.levels))
    }
    else if (ncol(ps) == nunique(treat)) {
      if (!all(t.levels %in% colnames(ps))) {
        .err("the columns of `ps` must be named with the treatment levels")
      }
    }
    else {
      .err("`ps` must have as many columns as there are treatment levels")
    }

  }

  ps
}
subclass_ps <- function(ps_mat, treat, estimand = "ATE", focal = NULL, subclass) {
  chk::chk_count(subclass)
  subclass <- round(subclass)

  if (!toupper(estimand) %in% c("ATE", "ATT")) {
    .err("only the ATE, ATT, and ATC are compatible with stratification weights")
  }

  if (is_not_null(focal)) {
    ps_mat <- ps_mat[,c(focal, setdiff(colnames(ps_mat), focal))]
  }

  ps_sub <- sub_mat <- ps_mat * 0

  for (i in colnames(ps_mat)) {
    if (toupper(estimand) == "ATE") {
      sub <- as.integer(findInterval(ps_mat[, as.character(i)],
                                     quantile(ps_mat[, as.character(i)],
                                              seq(0, 1, length.out = subclass + 1)),
                                     all.inside = TRUE))
    }
    else if (toupper(estimand) == "ATT") {
      if (i != focal) ps_mat[, as.character(i)] <- 1 - ps_mat[, as.character(i)]
      sub <- as.integer(findInterval(ps_mat[, as.character(i)],
                                     quantile(ps_mat[treat == focal, as.character(i)],
                                              seq(0, 1, length.out = subclass + 1)),
                                     all.inside = TRUE))
    }

    sub.tab <- table(treat, sub)

    if (any(sub.tab == 0)) {
      # .err("Too many subclasses were requested")
      sub <- subclass_scoot(sub, treat, ps_mat[,i])
      sub.tab <- table(treat, sub)
    }

    sub <- as.character(sub)

    sub.totals <- colSums(sub.tab)
    sub.ps <- setNames(sub.tab[as.character(i), ] / sub.totals,
                       colnames(sub.tab))

    ps_sub[,i] <- sub.ps[sub]
    sub_mat[,i] <- sub

    if (ncol(ps_sub) == 2) {
      ps_sub[,colnames(ps_sub) != i] <- 1 - ps_sub[,i]
      sub_mat[,colnames(sub_mat) != i] <- sub
      break
    }
  }

  attr(ps_sub, "sub_mat") <- sub_mat

  ps_sub
}
subclass_scoot <- function(sub, treat, x, min.n = 1) {
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

  if (any({sub_tab <- table(treat, sub)} == 0)) {

    soft_thresh <- function(x, minus = 1) {
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
               sum(soft_thresh(sub_tab[t, seq(1, first_0 - 1)]) / abs(first_0 - seq(1, first_0 - 1))) >=
               sum(soft_thresh(sub_tab[t, seq(first_0 + 1, nsub)]) / abs(first_0 - seq(first_0 + 1, nsub))))) {
            #If there are more and closer nonzero subs to the left...
            first_non0_to_left <- max(seq(1, first_0 - 1)[sub_tab[t, seq(1, first_0 - 1)] > 0])

            name_to_move <- names(sub)[which(x == max(x[treat == t & sub == first_non0_to_left]) & treat == t & sub == first_non0_to_left)[1]]

            sub[name_to_move] <- first_0
            sub_tab[t, first_0] <- 1L
            sub_tab[t, first_non0_to_left] <- sub_tab[t, first_non0_to_left] - 1L

          }
          else {
            #If there are more and closer nonzero subs to the right...
            first_non0_to_right <- min(seq(first_0 + 1, nsub)[sub_tab[t, seq(first_0 + 1, nsub)] > 0])
            name_to_move <- names(sub)[which(x == min(x[treat == t & sub == first_non0_to_right]) & treat == t & sub == first_non0_to_right)[1]]
            sub[name_to_move] <- first_0
            sub_tab[t, first_0] <- 1L
            sub_tab[t, first_non0_to_right] <- sub_tab[t, first_non0_to_right] - 1L
          }
        }

        sub_tab[t,] <- sub_tab[t,] - 1
      }
    }

    #Unsort
    sub <- sub[names(sub)]
  }

  sub
}
stabilize_w <- function(weights, treat) {
  if (is.factor(treat)) t.levels <- levels(treat)
  else t.levels <- unique(treat)

  w.names <- names(weights)
  tab <- setNames(vapply(t.levels, function(x) mean_fast(treat == x), numeric(1L)), t.levels)

  setNames(weights * tab[as.character(treat)], w.names)
}
`%+%` <- function(...) {
  if (is.atomic(..1) && is.atomic(..2)) crayon::`%+%`(as.character(..1), as.character(..2))
  else ggplot2::`%+%`(...)
}

get_cont_weights <- function(ps, treat, s.weights, dens.num, densfun = dnorm, use.kernel = FALSE,
                             densControl = list(bw = "nrd0", n = 10 * length(treat),
                                                adjust = 1, kernel = "gaussian")) {

  if (!is.matrix(ps)) ps <- matrix(ps, ncol = 1)

  p.denom <- treat - ps

  if (use.kernel) {
    s.weights <- s.weights/sum(s.weights)
    dens.denom <- apply(p.denom, 2, function(p) {
      d.d <- density(p, n = densControl[["n"]],
                     weights = s.weights, give.Rkern = FALSE,
                     bw = densControl[["bw"]], adjust = densControl[["adjust"]],
                     kernel = densControl[["kernel"]])
      with(d.d, approxfun(x = x, y = y))(p)
    })
  }
  else {
    dens.denom <- densfun(mat_div(p.denom, sqrt(col.w.v(p.denom))))
    if (is_null(dens.denom) || !is.numeric(dens.denom) || anyNA(dens.denom)) {
      .err("there was a problem with the output of `density`. Try another density function or leave it blank to use the normal density")
    }
    if (any(dens.denom <= 0)) {
      .err("the input to `density` may not accept the full range of standardized residuals of the treatment model")
    }

  }

  drop(dens.num/dens.denom)
}

get.w.from.ps <- function(ps, treat, estimand = "ATE", focal = NULL, subclass = NULL, stabilize = FALSE) {
  #Batch turn PS into weights; primarily for output of predict.gbm
  # Assumes a (0,1) treatment if binary, with ATT already processed
  if (is_null(dim(ps))) {
    ps <- matrix(ps, ncol = 1)
  }

  if (length(dim(ps)) == 2) {
    #Binary treatment, vector ps
    if (is_not_null(focal)) focal <- "1"
    if (is_not_null(subclass)) {
      #Get MMW subclass propensity scores
      for (p in seq_col(ps)) {
        ps_mat <- matrix(c(1 - ps[,p], ps[,p]), ncol = 2, dimnames = list(rownames(ps), c("0", "1")))
        ps[,p] <- subclass_ps(ps_mat, treat, estimand, focal, subclass)[, 2]
      }
    }

    if (estimand == "ATE") {
      w <- treat/ps + (1-treat)/(1-ps)
    }
    else if (estimand == "ATT") {
      w <- treat + (1-treat)*ps/(1-ps)
    }
    else if (estimand == "ATO") {
      w <- ps * (1-ps)
    }
    else if (estimand == "ATM") {
      w <- (treat/ps + (1-treat)/(1-ps))
      w <- w * pmin(ps, 1 - ps)
    }

    if (stabilize) {
      for (i in 0:1) {
        w[treat == i] <- mean_fast(treat == i) * w[treat == i]
      }
    }
  }
  else if (length(dim(ps)) == 3) {
    #Multi-category treatment, matrix PS

    if (is_not_null(subclass)) {
      #Get MMW subclass propensity scores
      for (p in seq_len(last(dim(ps))))
        ps[,,p] <- subclass_ps(ps[,,p], treat, estimand, focal, subclass)
    }

    w <- matrix(0, ncol = dim(ps)[3], nrow = dim(ps)[1])
    t.levs <- unique(treat)
    for (i in t.levs) w[treat == i,] <- 1/ps[treat == i, as.character(i),]

    if (estimand == "ATE") {
    }
    else if (estimand == "ATT") {
      w <- w * ps[, as.character(focal),]
    }
    else if (estimand == "ATO") {
      w <- w / colSums(aperm(1/ps, c(2,1,3)))
    }
    else if (estimand == "ATM") {
      for (p in seq_len(dim(ps)[3])) {
        w[,p] <- w[,p] * do.call("pmin", lapply(seq_len(dim(ps)[2]), function(i) ps[,i,p]))
      }
    }

    if (stabilize) {
      for (i in t.levs) {
        w[treat == i,] <- mean_fast(treat == i)*w[treat == i,]
      }
    }
  }
  else .err("don't know how to process more than 3 dims (likely a bug)")

  w
}

plot_density <- function(d.n, d.d) {
  d.d_ <- cbind(as.data.frame(d.d[c("x", "y")]), dens = "Denominator Density", stringsAsfactors = FALSE)
  d.n_ <- cbind(as.data.frame(d.n[c("x", "y")]), dens = "Numerator Density", stringsAsfactors = FALSE)
  d.all <- rbind(d.d_, d.n_)
  d.all$dens <- factor(d.all$dens, levels = c("Numerator Density", "Denominator Density"))
  pl <- ggplot(d.all, aes(x = x, y = y)) + geom_line() +
    labs(title = "Weight Component Densities", x = "E[Treat|X]", y = "Density") +
    facet_grid(rows = vars(dens)) + theme(panel.background = element_rect(fill = "white"),
                                          panel.border = element_rect(fill = NA, color = "black"),
                                          axis.text.x = element_text(color = "black"),
                                          axis.text.y = element_text(color = "black"),
                                          panel.grid.major = element_blank(),
                                          panel.grid.minor = element_blank()
    )
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
    for (i in colnames(covs)[anyNA_col(covs)]) {
      if (all(is.na(covs[,i]))) covs <- covs[, colnames(covs) != i, drop = FALSE]
      else covs[is.na(covs[,i]),i] <- match.fun(with)(covs[,i], na.rm = TRUE)
    }
    return(covs)
  }

  covs[is.na(covs)] <- with

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

#Choleski decomp for non-negative definite matrices
chol2 <- function(Sinv) {
  ch <- suppressWarnings(chol(Sinv, pivot = TRUE))
  p <- order(attr(ch, "pivot"))
  ch[, p, drop = FALSE]
}

#For balance SuperLearner
method.balance <- function() {

  out <- list(
    # require allows you to pass a character vector with required packages
    # use NULL if no required packages
    require = "cobalt",

    # computeCoef is a function that returns a list with two elements:
    # 1) coef: the weights (coefficients) for each algorithm
    # 2) cvRisk: the V-fold CV risk for each algorithm
    computeCoef = function(Z, Y, libraryNames, obsWeights, control, verbose, ...) {
      estimand <- attr(control$trimLogit, "vals")$estimand
      init <- attr(control$trimLogit, "vals")$init

      tol <- .001
      for (i in seq_col(Z)) {
        Z[Z[,i] < tol, i] <- tol
        Z[Z[,i] > 1-tol, i] <- 1-tol
      }
      w_mat <- get.w.from.ps(Z, treat = Y, estimand = estimand)
      cvRisk <- apply(w_mat, 2, cobalt::bal.compute, init = init)

      names(cvRisk) <- libraryNames

      loss <- function(coefs) {
        ps <- crossprod(t(Z), coefs/sum(coefs))
        w <- get_w_from_ps(ps, Y, estimand)
        cobalt::bal.compute(init, weights = w)
      }
      fit <- optim(rep(1/ncol(Z), ncol(Z)), loss, method = "L-BFGS-B", lower = 0, upper = 1)
      coef <- fit$par

      list(cvRisk = cvRisk, coef = coef/sum(coef))
    },

    # computePred is a function that takes the weights and the predicted values
    # from each algorithm in the library and combines them based on the model to
    # output the super learner predicted values
    computePred = function(predY, coef, control, ...) {
      crossprod(t(predY), coef/sum(coef))
    }
  )

  out
}

method.balance.cont <- function() {

  out <- list(
    # require allows you to pass a character vector with required packages
    # use NULL if no required packages
    require = "cobalt",

    # computeCoef is a function that returns a list with two elements:
    # 1) coef: the weights (coefficients) for each algorithm
    # 2) cvRisk: the V-fold CV risk for each algorithm
    computeCoef = function(Z, Y, libraryNames, obsWeights, control, verbose, ...) {
      dens.num <- attr(control$trimLogit, "vals")$dens.num
      densfun <- attr(control$trimLogit, "vals")$densfun
      use.kernel <- attr(control$trimLogit, "vals")$use.kernel
      densControl <- attr(control$trimLogit, "vals")$densControl
      init <- attr(control$trimLogit, "vals")$init

      w_mat<- get_cont_weights(Z, treat = Y, s.weights = obsWeights,
                               dens.num = dens.num, densfun = densfun, use.kernel = use.kernel,
                               densControl = densControl)
      cvRisk <- apply(w_mat, 2, cobalt::bal.compute, init = init)
      names(cvRisk) <- libraryNames

      loss <- function(coefs) {
        ps <- crossprod(t(Z), coefs/sum(coefs))
        w <- get_cont_weights(ps, treat = Y, s.weights = obsWeights,
                              dens.num = dens.num, densfun = densfun,
                              use.kernel = use.kernel,
                              densControl = densControl)
        cobalt::bal.compute(init, weights = w)
      }
      fit <- optim(rep(1/ncol(Z), ncol(Z)), loss, method = "L-BFGS-B", lower = 0, upper = 1)
      coef <- fit$par

      list(cvRisk = cvRisk, coef = coef/sum(coef))
    },

    # computePred is a function that takes the weights and the predicted values
    # from each algorithm in the library and combines them based on the model to
    # output the super learner predicted values
    computePred = function(predY, coef, control, ...) {
      crossprod(t(predY), coef/sum(coef))
    }
  )

  out
}

.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}

#To pass CRAN checks:
utils::globalVariables(c(".s.weights", "dens", "x", "y"))
