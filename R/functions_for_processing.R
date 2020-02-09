method.to.proper.method <- function(method) {
  method <- tolower(method)
  if      (method %in% c("ps")) return("ps")
  else if (method %in% c("gbm", "gbr")) return("gbm")
  else if (method %in% c("twang")) return("twang")
  else if (method %in% c("cbps", "cbgps")) return("cbps")
  else if (method %in% c("npcbps", "npcbgps")) return("npcbps")
  else if (method %in% c("entropy", "ebal", "ebalance")) return("ebal")
  else if (method %in% c("ebcw", "ate")) return("ebcw")
  else if (method %in% c("optweight", "opt", "sbw")) return("optweight")
  else if (method %in% c("super", "superlearner")) return("super")
  # else if (method %in% c("kbal")) return("kbal")
  else return(method)
}
check.acceptable.method <- function(method, msm = FALSE, force = FALSE) {
  bad.method <- FALSE
  acceptable.methods <- c("ps"
                          , "gbm", "gbr"
                          , "twang"
                          , "cbps", "cbgps"
                          , "npcbps", "npcbgps"
                          , "ebal", "entropy", "ebalance"
                          , "sbw"
                          , "ebcw", "ate"
                          , "optweight", "opt", "sbw"
                          , "super", "superlearner"
                          # "kbal",
  )

  if (missing(method)) method <- "ps"
  else if (is_null(method) || length(method) > 1) bad.method <- TRUE
  else if (is.character(method)) {
    if (tolower(method) %nin% acceptable.methods) bad.method <- TRUE
  }
  else if (!is.function(method)) bad.method <- TRUE

  if (bad.method) stop("method must be a string of length 1 containing the name of an acceptable weighting method or a function that produces weights.", call. = FALSE)

  if (msm && !force && is.character(method)) {
    m <- method.to.proper.method(method)
    if (m %in% c("nbcbps", "ebal", "ebcw", "optweight", "kbal")) {
      stop(paste0("The use of ", method.to.phrase(m), " with longitudinal treatments has not been validated. Set weightit.force = TRUE to bypass this error."), call. = FALSE)
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
    stop("The user-provided function to method must contain \"covs\" and \"treat\" as named parameters.", call. = FALSE)
  }
}
method.to.phrase <- function(method) {

  if (is.function(method)) return("a user-defined method")
  else {
    method <- method.to.proper.method(method)
    if (method %in% c("ps")) return("propensity score weighting")
    else if (method %in% c("gbm")) return("generalized boosted modeling")
    else if (method %in% c("twang")) return("generalized boosted modeling with TWANG")
    else if (method %in% c("cbps")) return("covariate balancing propensity score weighting")
    else if (method %in% c("npcbps")) return("non-parametric covariate balancing propensity score weighting")
    else if (method %in% c("ebal")) return("entropy balancing")
    else if (method %in% c("ebcw")) return("empirical balancing calibration weighting")
    else if (method %in% c("optweight")) return("targeted stable balancing weights")
    else if (method %in% c("super")) return("propensity score weighting with SuperLearner")
    # else if (method %in% c("kbal")) return("kernel balancing")
    else return("the chosen method of weighting")
  }
}
process.estimand <- function(estimand, method, treat.type) {
  #Allowable estimands
  AE <- list(
    binary = list(ps = c("ATT", "ATC", "ATE", "ATO", "ATM")
                  , gbm = c("ATT", "ATC", "ATE", "ATO", "ATM")
                  , twang = c("ATT", "ATC", "ATE")
                  , cbps = c("ATT", "ATC", "ATE")
                  , npcbps = c("ATE")
                  , ebal = c("ATT", "ATC", "ATE")
                  , ebcw = c("ATT", "ATC", "ATE")
                  , optweight = c("ATT", "ATC", "ATE")
                  , super = c("ATT", "ATC", "ATE", "ATO", "ATM")
                  # , kbal = c("ATT", "ATC", "ATE")
    ),
    multinomial = list(ps = c("ATT", "ATC", "ATE", "ATO", "ATM")
                       , gbm = c("ATT", "ATC", "ATE", "ATO", "ATM")
                       , twang = c("ATT", "ATC", "ATE")
                       , cbps = c("ATT", "ATC", "ATE")
                       , npcbps = c("ATE")
                       , ebal = c("ATT", "ATC", "ATE")
                       , ebcw = c("ATT", "ATC", "ATE")
                       , optweight = c("ATT", "ATC", "ATE")
                       , super = c("ATT", "ATC", "ATE", "ATO", "ATM")
                       # , kbal = c("ATT", "ATE")
    ))

  if (treat.type != "continuous" && !is.function(method)) {
    if (is_null(estimand)) stop(paste0("estimand must be one of ", word_list(AE[[treat.type]][[method]], quotes = TRUE, and.or = "or"), "."), call. = FALSE)
    else if (toupper(estimand) %nin% AE[[treat.type]][[method]]) {
      stop(paste0("\"", estimand, "\" is not an allowable estimand for ", method.to.phrase(method),
                  " with ", treat.type, " treatments. Only ", word_list(AE[[treat.type]][[method]], quotes = TRUE, and.or = "and", is.are = TRUE),
                  " allowed."), call. = FALSE)
    }
  }
  return(toupper(estimand))
}
check.subclass <- function(method, treat.type) {
  #Allowable estimands
  AE <- list(
    binary = list(ps = TRUE
                  , gbm = TRUE
                  , twang = FALSE
                  , cbps = TRUE
                  , npcbps = FALSE
                  , ebal = FALSE
                  , ebcw = FALSE
                  , optweight = FALSE
                  , super = TRUE
                  # , kbal = FALSE
    ),
    multinomial = list(ps = TRUE
                       , gbm = TRUE
                       , twang = FALSE
                       , cbps = FALSE
                       , npcbps = FALSE
                       , ebal = FALSE
                       , ebcw = FALSE
                       , optweight = FALSE
                       , super = TRUE
    ))

  if (treat.type != "continuous" && !is.function(method) &&
      !AE[[treat.type]][[method]]) {
    stop(paste0("subclasses are not compatible with ", method.to.phrase(method),
                " with ", treat.type, " treatments."), call. = FALSE)
  }
}
process.focal.and.estimand <- function(focal, estimand, treat, treat.type, treated = NULL) {
  reported.estimand <- estimand

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  #Check focal
  if (treat.type == "multinomial") {
    unique.treat <- unique(treat, nmax = length(treat)/4)
    if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
      warning(paste(estimand, "is not compatible with focal. Setting estimand to \"ATT\"."), call. = FALSE)
      reported.estimand <- estimand <- "ATT"
    }

    if (estimand == "ATT") {
      if (is_null(focal)) {
        if (is_null(treated) || treated %nin% unique.treat) {
          stop("When estimand = \"ATT\" for multinomial treatments, an argument must be supplied to focal.", call. = FALSE)
        }
        focal <- treated
      }
    }
    else if (estimand == "ATC") {
      if (is_null(focal)) {
        stop("When estimand = \"ATC\" for multinomial treatments, an argument must be supplied to focal.", call. = FALSE)
      }
    }
  }
  else if (treat.type == "binary") {
    unique.treat <- unique(treat, nmax = 2)
    unique.treat.bin <- unique(binarize(treat), nmax = 2)
    if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
      warning(paste(estimand, "is not compatible with focal. Setting estimand to \"ATT\"."), call. = FALSE)
      reported.estimand <- estimand <- "ATT"
    }

    if (estimand == "ATT") {
      if (is_null(focal)) {
        if (is_null(treated) || treated %nin% unique.treat) {
          if (all(as.character(unique.treat.bin) == as.character(unique.treat))) {
            #If is 0/1
            treated <- unique.treat[unique.treat.bin == 1]
          }
          else {
            treated <- names(which.min(table(treat))) #Smaller group is treated
            message(paste0("Assuming ", word_list(treated, quotes = !is.numeric(treat), is.are = TRUE),
                           " the treated level. If not, supply an argument to focal."))
          }
        }
        focal <- treated
      }
      else {
        if (is_null(treated) || treated %nin% unique.treat) {
          treated <- focal
        }
      }
    }
    else if (estimand == "ATC") {
      if (is_null(focal)) {
        if (is_null(treated) || treated %nin% unique.treat) {

          if (all(as.character(unique.treat.bin) == as.character(unique.treat))) {
            treated <- unique.treat[unique.treat.bin == 1]
          }
          else {
            treated <- names(which.min(table(treat))) #Smaller group is treated
            message(paste0("Assuming ", word_list(unique.treat[unique.treat %nin% treated], quotes = !is.numeric(treat), is.are = TRUE),
                           " the control level. If not, supply an argument to focal."))
          }
        }
        focal <- unique.treat[unique.treat %nin% treated]
      }
      else {
        if (is_null(treated) || treated %nin% unique.treat) {
          treated <- unique.treat[unique.treat %nin% focal]
        }
      }
      estimand <- "ATT"
    }
  }

  if (is_not_null(focal) && (length(focal) > 1L || focal %nin% unique.treat)) {
    stop("The argument supplied to focal must be the name of a level of treat.", call. = FALSE)
  }

  return(list(focal = as.character(focal),
              estimand = estimand,
              reported.estimand = reported.estimand,
              treated = if (is.factor(treated)) as.character(treated) else treated))
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
  else if (is.formula(by, 1)) {
    t.c <- get.covs.and.treat.from.formula(by, data)
    by <- t.c[["reported.covs"]]
    if (NCOL(by) != 1) stop(paste0("Only one variable can be on the right-hand side of the formula for ", by.arg, "."), call. = FALSE)
    else by.name <- colnames(by)
  }
  else bad.by <- TRUE

  if (!bad.by) {
    by.components <- data.frame(by)

    if (is_not_null(colnames(by))) names(by.components) <- colnames(by)
    else names(by.components) <- by.name

    if (is_null(by)) by.factor <- factor(rep(1, n))
    else by.factor <- factor(by.components[[1]], levels = unique(by.components[[1]]),
                             labels = paste(names(by.components), "=", unique(by.components[[1]])))
    # by.vars <- acceptable.bys[vapply(acceptable.bys, function(x) equivalent.factors(by, data[[x]]), logical(1L))]
  }
  else {
    stop(paste(by.arg, "must be a string containing the name of the variable in data for which weighting is to occur within strata or a one-sided formula with the stratifying variable on the right-hand side."), call. = FALSE)
  }

  if (treat.type != "continuous" && any(vapply(levels(by.factor), function(x) nunique(treat) != nunique(treat[by.factor == x]), logical(1L)))) {
    stop(paste0("Not all the groups formed by ", by.arg, " contain all treatment levels", if (is_not_null(treat.name)) paste("in", treat.name) else "", ". Consider coarsening ", by.arg, "."), call. = FALSE)
  }
  attr(by.components, "by.factor") <- by.factor
  return(by.components)
}
process.moments.int <- function(moments, int, method) {
  if (!is.function(method)) {
    if (method %in% c("npcbps", "ebal", "ebcw", "optweight")) {
      if (length(int) != 1 || !is.logical(int)) {
        stop("int must be a logical (TRUE/FALSE) of length 1.", call. = FALSE)
      }
      if (length(moments) != 1 || !is.numeric(moments) ||
          !check_if_zero(moments - round(moments)) ||
          moments < 1) {
        stop("moments must be a positive integer of length 1.", call. = FALSE)
      }
    }
    else if (any(mi0 <- c(as.integer(moments) != 1L, int))) {
      warning(paste0(word_list(c("moments", "int")[mi0], and.or = "and", is.are = TRUE),
                     " not compatible with ", method.to.phrase(method), ". Ignoring ", word_list(c("moments", "int")[mi0], and.or = "and"), "."), call. = FALSE)
      moments <- 1
      int <- FALSE
    }
  }
  return(c(moments = as.integer(moments), int = int))
}
process.MSM.method <- function(is.MSM.method, method) {
  methods.with.MSM <- c("optweight")
  if (is.function(method)) {
    if (isTRUE(is.MSM.method)) stop("Currently, only user-defined methods that work with is.MSM.method = FALSE are allowed.", call. = FALSE)
    is.MSM.method <- FALSE
  }
  else if (method %in% methods.with.MSM) {
    if (is_null(is.MSM.method)) is.MSM.method <- TRUE
    else if (!isTRUE(is.MSM.method)) {
      message(paste0(method.to.phrase(method), " can be used with a single model when multiple time points are present.\nUsing a seperate model for each time point. To use a single model, set is.MSM.method to TRUE."))
    }
  }
  else {
    if (isTRUE(is.MSM.method)) warning(paste0(method.to.phrase(method), " cannot be used with a single model when multiple time points are present.\nUsing a seperate model for each time point."),
                                       call. = FALSE, immediate. = TRUE)
    is.MSM.method <- FALSE
  }

  return(is.MSM.method)

}
process.missing <- function(missing, method, treat.type) {
  #Allowable estimands
  AE <- list(binary = list(ps = c("ind", "saem")
                           , gbm = c("ind", "surr")
                           , twang = c("ind", "surr")
                           , cbps = c("ind")
                           , npcbps = c("ind")
                           , ebal = c("ind")
                           , ebcw = c("ind")
                           , optweight = c("ind")
                           , super = c("ind")
                           # , kbal = c("ind")
  ),
  multinomial = list(ps = c("ind")
                     , gbm = c("ind", "surr")
                     , twang = c("ind", "surr")
                     , cbps = c("ind")
                     , npcbps = c("ind")
                     , ebal = c("ind")
                     , ebcw = c("ind")
                     , optweight = c("ind")
                     , super = c("ind")
                     # , kbal = c("ind")
  ),
  continuous = list(ps = c("ind")
                    , gbm = c("ind", "surr")
                    , twang = c("ind", "surr")
                    , cbps = c("ind")
                    , npcbps = c("ind")
                    , ebal = c("ind")
                    , ebcw = c("ind")
                    , optweight = c("ind")
                    , super = c("ind")
                    # , kbal = c("ind")
  ))

  if (!is.character(missing) || length(missing) != 1) stop("missing must be a string of length 1.", call. = FALSE)

  allowable.missings <- AE[[treat.type]][[method]]
  if (is_null(missing)) {
    missing <- allowable.missings[1]
    warning(paste0("Missing values are present in the covariates. See ?WeightIt::method_", method, " for information on how these are handled."), call. = FALSE)
  }
  else if (missing %pin% allowable.missings) {
    missing <- allowable.missings[pmatch(missing, allowable.missings)]
  }
  else {
    missing <- allowable.missings[1]
    warning(paste0("Only ", word_list(allowable.missings, quotes = TRUE, is.are = TRUE), " allowed for missing with ",
                   treat.type,
                   " treatments. Using link = ", word_list(allowable.missings[1], quotes = TRUE), "."),
            call. = FALSE, immediate. = TRUE)
  }
  return(missing)
}
check.package <- function(package.name, alternative = FALSE) {
  packages.not.installed <- package.name[package.name %nin% .packages(all.available = TRUE)]
  if (is_not_null(packages.not.installed)) {
    if (alternative) return(FALSE)
    else {
      plural <- length(packages.not.installed) > 1
      stop(paste0("Package", if (plural) "s " else " ",
                  word_list(packages.not.installed, quotes = TRUE, is.are = TRUE),
                  " needed for this function to work. Please install ",
                  if (plural) "them" else "it","."),
           call. = FALSE)
    }
  }
  else return(invisible(TRUE))
}
make.closer.to.1 <- function(x) {
  if (is.factor(x) || is.character(x) || all_the_same(x[!is.na(x)])) return(x)
  else if (is_binary(x)) {
    return(as.numeric(x == x[!is.na(x)][1]))
  }
  else {
    (x - mean_fast(x, TRUE))/sd(x, na.rm = TRUE)
  }
}
int.poly.f <- function(mat, ex = NULL, int = FALSE, poly = 1, center = FALSE, sep = " * ") {
  #Adds to data frame interactions and polynomial terms; interaction terms will be named "v1 * v2" and polynomials will be named "v1_2"
  #mat=matrix input
  #ex=matrix of variables to exclude in interactions and polynomials; a subset of df
  #int=whether to include interactions or not; currently only 2-way are supported
  #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included
  #nunder=number of underscores between variables

  if (is_not_null(ex)) d <- mat[, colnames(mat) %nin% colnames(ex), drop = FALSE]
  else d <- mat
  binary.vars <- apply(d, 2, is_binary)
  if (center) {
    d[,!binary.vars] <- center(d[, !binary.vars, drop = FALSE])
  }
  nd <- ncol(d)
  nrd <- nrow(d)
  no.poly <- binary.vars
  npol <- nd - sum(no.poly)
  new <- matrix(0, ncol = (poly-1)*npol + int*(.5*(nd)*(nd-1)), nrow = nrd)
  nc <- ncol(new)
  new.names <- character(nc)
  if (poly > 1 && npol != 0) {
    for (i in 2:poly) {
      new[, (1 + npol*(i - 2)):(npol*(i - 1))] <- apply(d[, !no.poly, drop = FALSE], 2, function(x) x^i)
      new.names[(1 + npol*(i - 2)):(npol*(i - 1))] <- paste0(colnames(d)[!no.poly], num_to_superscript(i))
    }
  }
  if (int && nd > 1) {
    new[,(nc - .5*nd*(nd-1) + 1):nc] <- matrix(t(apply(d, 1, combn, 2, prod)), nrow = nrd)
    new.names[(nc - .5*nd*(nd-1) + 1):nc] <- combn(colnames(d), 2, paste, collapse = sep)
  }

  single.value <- apply(new, 2, all_the_same)
  colnames(new) <- new.names

  return(new[, !single.value, drop = FALSE])
}
get.s.d.denom.weightit <- function(s.d.denom = NULL, estimand = NULL, weights = NULL, treat = NULL, focal = NULL) {
  check.estimand <- check.weights <- check.focal <- bad.s.d.denom <- bad.estimand <- FALSE
  s.d.denom.specified <- is_not_null(s.d.denom)
  estimand.specified <- is_not_null(estimand)

  if (is_not_null(weights) && !is.data.frame(weights)) weights <- data.frame(weights)

  if (s.d.denom.specified) {
    try.s.d.denom <- tryCatch(match_arg(s.d.denom, c("treated", "control", "pooled", "all"), several.ok = TRUE),
                              error = function(cond) FALSE)
    if (any(try.s.d.denom == FALSE)) {
      check.estimand <- TRUE
      bad.s.d.denom <- TRUE
    }
    else {
      if (length(try.s.d.denom) > 1 && length(try.s.d.denom) != ncol(weights)) {
        stop("s.d.denom must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
      }
      else s.d.denom <- try.s.d.denom
    }
  }
  else {
    check.estimand <- TRUE
  }

  if (check.estimand == TRUE) {
    if (estimand.specified) {
      try.estimand <- tryCatch(match_arg(tolower(estimand), c("att", "atc", "ate"), several.ok = TRUE),
                               error = function(cond) FALSE)
      if (any(try.estimand == FALSE)) {
        check.focal <- TRUE
        bad.estimand <- TRUE
      }
      else {
        if (length(try.estimand) > 1 && length(try.estimand) != ncol(weights)) {
          stop("estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
        }
        else s.d.denom <- vapply(try.estimand, switch, character(1L), att = "treated", atc = "control", ate = "pooled")
      }
    }
    else {
      check.focal <- TRUE
    }
  }
  if (check.focal == TRUE) {
    if (is_not_null(focal)) {
      s.d.denom <- "treated"
      estimand <- "att"
    }
    else check.weights <- TRUE
  }
  if (check.weights == TRUE) {
    if (is_null(weights)) {
      s.d.denom <- "pooled"
      estimand <- "ate"
    }
    else {
      s.d.denom <- estimand <- character(ncol(weights))
      for (i in seq_len(ncol(weights))) {
        if (is_binary(treat)) {
          if (all_the_same(weights[[i]][treat==1 & !check_if_zero(weights[[i]])]) &&
              !all_the_same(weights[[i]][treat==0 & !check_if_zero(weights[[i]])])
          ) { #if treated weights are the same and control weights differ; ATT
            estimand[i] <- "att"
            s.d.denom[i] <- "treated"
          }
          else if (all_the_same(weights[[i]][treat==0 & !check_if_zero(weights[[i]])]) &&
                   !all_the_same(weights[[i]][treat==1 & !check_if_zero(weights[[i]])])
          ) { #if control weights are the same and treated weights differ; ATC
            estimand[i] <- "atc"
            s.d.denom[i] <- "control"
          }
          else {
            estimand[i] <- "ate"
            s.d.denom[i] <- "pooled"
          }
        }
        else {
          if (length(focal) == 1) {
            estimand[i] <- "att"
            s.d.denom[i] <- "treated"
          }
          else {
            estimand[i] <- "ate"
            s.d.denom[i] <- "pooled"
          }
        }
      }
    }
  }
  if (is_not_null(weights) && length(s.d.denom) == 1) s.d.denom <- rep(s.d.denom, ncol(weights))

  if (s.d.denom.specified && bad.s.d.denom && (!estimand.specified || bad.estimand)) {
    # message("Warning: s.d.denom should be one of \"treated\", \"control\", \"pooled\", or \"all\".\n         Using \"", word_list(s.d.denom), "\" instead.")
  }
  else if (estimand.specified && bad.estimand) {
    # message("Warning: the supplied estimand is not allowed. Using \"", ifelse(all_the_same(estimand), toupper(estimand)[1], word_list(toupper(estimand))), "\" instead.")
  }
  else if (check.focal || check.weights) {
    # message("Note: estimand not specified; assuming ", ifelse(all_the_same(toupper(estimand)), toupper(unique(estimand)), word_list(toupper(estimand))), ".")
  }

  if (is_not_null(weights) && length(s.d.denom) != ncol(weights)) {
    # stop("Valid inputs to s.d.denom or estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
  }

  return(s.d.denom)
}
ps_to_ps_mat <- function(ps, treat, assumed.treated = NULL, treat.type = NULL, treated = NULL, estimand = NULL) {
  if (is.vector(ps, mode = "numeric")) {
    ps.names <- names(ps)
    ps <- matrix(ps, ncol = 1)
  }
  else if (is.matrix(ps) || is.data.frame(ps)) {
    ps.names <- rownames(ps)
  }

  if (is.factor(treat)) t.levels <- levels(treat)
  else t.levels <- unique(treat, nmax = length(treat)/4)

  if (treat.type == "binary") {
    if ((is.matrix(ps) && all(is.numeric(ps))) ||
        (is.data.frame(ps) && all(vapply(ps, is.numeric, logical(1L))))) {
      if (ncol(ps) == 1) {
        if (can_str2num(treat) &&
            all(check_if_zero(binarize(treat) - str2num(treat)))) treated.level <- 1
        else if (is_not_null(treated)) {
          if (treated %in% treat) treated.level <- treated
          else stop("The argument to treated must be a value in treat.", call. = FALSE)
        }
        else if (is_not_null(assumed.treated)) {
          treated.level <- assumed.treated
        }
        else if (is_not_null(colnames(ps)) && colnames(ps) %in% as.character(t.levels)) {
          treated.level <- colnames(ps)
        }
        else {
          stop("If the treatment has two non-0/1 levels and ps is a vector or has only one column, an argument to treated must be supplied.", call. = FALSE)
        }

        t.levels <- c(treated.level, t.levels[t.levels != treated.level])
        ps <- matrix(c(ps[,1], 1-ps[,1]), ncol = 2, dimnames = list(ps.names, as.character(t.levels)))
      }
      else if (ncol(ps) == 2) {
        if (!all(as.character(t.levels) %in% colnames(ps))) {
          stop("If ps has two columns, they must be named with the treatment levels.", call. = FALSE)
        }
        else ps <- ps
      }
      else {
        stop("ps cannot have more than two columns if the treatment is binary.", call. = FALSE)
      }
    }
    else {
      stop("ps must be a matrix, data.frame, or vector of propensity scores.", call. = FALSE)
    }
  }
  else if (treat.type == "multinomial") {
    if ((is.matrix(ps) && all(is.numeric(ps))) ||
        (is.data.frame(ps) && all(vapply(ps, is.numeric, logical(1L))))) {
      if (ncol(ps) == 1) {
        if (toupper(estimand) == "ATE") {
          ps <- matrix(rep(ps, nunique(treat)), ncol = nunique(treat), dimnames = list(ps.names, t.levels))
        }
        else {
          stop("With multinomial treatments, ps can be a vector or have only one column only if the estimand is the ATE.", call. = FALSE)
        }
      }
      else if (ncol(ps) == nunique(treat)) {
        if (!all(t.levels %in% colnames(ps))) {
          stop("The columns of ps must be named with the treatment levels.", call. = FALSE)
        }
        else ps <- ps
      }
      else {
        stop("ps must have as many columns as there are treatment levels.", call. = FALSE)
      }
    }
    else {
      stop("ps must be a matrix or data.frame of propensity scores.", call. = FALSE)
    }
  }
  return(ps)
}
subclass_ps <- function(ps_mat, treat, estimand = "ATE", focal = NULL, subclass) {
  if (!length(subclass) == 1 || !is.numeric(subclass)) {
    stop("subclass must be a single number.", call. = FALSE)
  }
  else if (round(subclass) < 1) {
    stop("subclass must be greater than 1.", call. = FALSE)
  }
  subclass <- round(subclass)

  if (is_not_null(focal)) ps_mat <- ps_mat[,c(focal, setdiff(colnames(ps_mat), focal))]

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
    else {
      stop("Only the ATE, ATT, and ATC are compatible with stratification weights.")
    }

    sub.tab <- table(treat, sub)

    if (any(sub.tab == 0)) {
      # stop("Too many subclasses were requested.", call. = FALSE)
      sub <- subclass_scoot(sub, treat, ps_mat[,i])
      sub.tab <- table(treat, sub)
    }

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
  return(ps_sub)
}
subclass_scoot <- function(sub, treat, x) {

  treat <- as.character(treat)
  unique.treat <- unique(treat, nmax = 2)

  names(x) <- seq_along(x)
  names(sub) <- seq_along(sub)
  original.order <- names(x)

  nsub <- nunique(sub)

  #Turn subs into a contiguous sequence
  sub <- setNames(setNames(seq_len(nsub), sort(unique(sub)))[as.character(sub)],
                  original.order)

  if (any(table(treat) < nsub)) {
    stop("Too many subclasses were requested.", call. = FALSE)
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
    }

    #Unsort
    sub <- sub[names(sub)]
  }

  return(sub)
}
stabilize_w <- function(weights, treat) {
  if (is.factor(treat)) t.levels <- levels(treat)
  else t.levels <- unique(treat)
  w.names <- names(weights)
  tab <- vapply(t.levels, function(x) mean_fast(treat == x), numeric(1L))
  return(setNames(weights * tab[as.character(treat)], w.names))
}
`%+%` <- function(rhs, lhs) {
  if (is_(rhs, "gg") && is_(lhs, "gg")) ggplot2::`%+%`(rhs, lhs)
  else crayon::`%+%`(as.character(rhs), as.character(lhs))
}

get_cont_weights <- function(ps, treat, s.weights, dens.num, densfun = dnorm, use.kernel = FALSE,
                             densControl = list(bw = "nrd0", n = 10*length(treat),
                                                adjust = 1, kernel = "gaussian")) {
  p.denom <- treat - ps

  if (isTRUE(densControl[["use.kernel"]])) {
    d.d <- density(p.denom, n = densControl[["n"]],
                   weights = s.weights/sum(s.weights), give.Rkern = FALSE,
                   bw = densControl[["bw"]], adjust = densControl[["adjust"]],
                   kernel = densControl[["kernel"]])
    dens.denom <- with(d.d, approxfun(x = x, y = y))(p.denom)
  }
  else {
    dens.denom <- densfun(p.denom/sd(p.denom))
    if (is_null(dens.denom) || !is.atomic(dens.denom) || anyNA(dens.denom)) {
      stop("There was a problem with the output of density. Try another density function or leave it blank to use the normal density.", call. = FALSE)
    }
    else if (any(dens.denom <= 0)) {
      stop("The input to density may not accept the full range of standardized residuals of the treatment model.", call. = FALSE)
    }

  }

  w <- dens.num/dens.denom

  return(w)
}

method.balance <- function(stop.method) {

  if (startsWith(stop.method, "es.")) {
    stop.fun <- function(mat, treat, weights, s.d.denom, s.weights = NULL, bin.vars, subset = NULL) {
      cobalt::col_w_smd(mat, treat, weights, std = rep(TRUE, ncol(mat)), s.d.denom = s.d.denom, abs = TRUE,
                        s.weights = s.weights, subset = subset, bin.vars = bin.vars)
    }
  }
  else if (startsWith(stop.method, "ks.")) {
    stop.fun <- function(mat, treat, weights, s.d.denom, s.weights = NULL, bin.vars, subset = NULL) {
      cobalt::col_w_ks(mat, treat, weights, s.weights = s.weights, subset = subset, bin.vars = bin.vars)
    }
  }

  if (endsWith(stop.method, ".mean")) stop.sum <- mean
  else if (endsWith(stop.method, ".max")) stop.sum <- max
  else if (endsWith(stop.method, ".rms")) stop.sum <- function(x, ...) sqrt(mean(x^2, ...))

  out <- list(
  # require allows you to pass a character vector with required packages
  # use NULL if no required packages
  require = "cobalt",

  # computeCoef is a function that returns a list with two elements:
  # 1) coef: the weights (coefficients) for each algorithm
  # 2) cvRisk: the V-fold CV risk for each algorithm
  computeCoef = function(Z, Y, libraryNames, obsWeights, control, verbose, ...) {
    covs <- attr(control$trimLogit, "vals")$covs
    estimand <- attr(control$trimLogit, "vals")$estimand
    s.d.denom <- get.s.d.denom.weightit(estimand = estimand, treat = Y)
    bin.vars <- apply(covs, 2, is_binary)

    for (i in 1:ncol(Z)) {
      Z[Z[,i]<.001,i] <- .001
      Z[Z[,i]>1-.001,i] <- 1-.001
    }
    w_mat<- apply(Z, 2, get_w_from_ps, Y, estimand)
    cvRisk <- apply(w_mat, 2, function(w) stop.sum(stop.fun(covs, treat = Y,
                                                            weights = w,
                                                            s.weights = obsWeights,
                                                            s.d.denom = s.d.denom,
                                                            bin.vars = bin.vars)))
    names(cvRisk) <- libraryNames

    loss <- function(coefs) {
      ps <- crossprod(t(Z), coefs/sum(coefs))
      w <- get_w_from_ps(ps, Y, estimand)
      out <- stop.sum(stop.fun(covs, treat = Y,
                               weights = w,
                               s.weights = obsWeights,
                               s.d.denom = s.d.denom,
                               bin.vars = bin.vars))
      out
    }
    fit <- optim(rep(1, ncol(Z)), loss, method = "L-BFGS-B", lower = 0)
    coef <- fit$par
    out <- list(cvRisk = cvRisk, coef = coef/sum(coef))
    return(out)
  },

  # computePred is a function that takes the weights and the predicted values
  # from each algorithm in the library and combines them based on the model to
  # output the super learner predicted values
  computePred = function(predY, coef, control, ...) {
    out <- crossprod(t(predY), coef/sum(coef))
    return(out)
  }
  )
  # attr(out, "stop.method") <- stop.method
  # class(out) <- "method.balance"
  return(out)
}

method.balance.cont <- function(stop.method) {

  if (startsWith(stop.method, "s.")) {
    stop.fun <- function(mat, treat, weights, s.weights, bin.vars, subset = NULL) {
      cobalt::col_w_corr(mat, treat, weights, type = "spearman", abs = TRUE,
                         s.weights = s.weights, bin.vars = bin.vars, subset = subset)
    }
  }
  else if (startsWith(stop.method, "p.")) {
    stop.fun <- function(mat, treat, weights, s.weights, bin.vars, subset = NULL) {
      cobalt::col_w_corr(mat, treat, weights, type = "pearson", abs = TRUE,
                         s.weights = s.weights, bin.vars = bin.vars, subset = subset)
    }
  }

  if (endsWith(stop.method, ".mean")) stop.sum <- mean
  else if (endsWith(stop.method, ".max")) stop.sum <- max
  else if (endsWith(stop.method, ".rms")) stop.sum <- function(x, ...) sqrt(mean(x^2, ...))

  out <- list(
    # require allows you to pass a character vector with required packages
    # use NULL if no required packages
    require = "cobalt",

    # computeCoef is a function that returns a list with two elements:
    # 1) coef: the weights (coefficients) for each algorithm
    # 2) cvRisk: the V-fold CV risk for each algorithm
    computeCoef = function(Z, Y, libraryNames, obsWeights, control, verbose, ...) {
      covs <- attr(control$trimLogit, "vals")$covs
      dens.num <- attr(control$trimLogit, "vals")$dens.num
      densfun <- attr(control$trimLogit, "vals")$densfun
      use.kernel <- attr(control$trimLogit, "vals")$use.kernel
      densControl <- attr(control$trimLogit, "vals")$densControl

      bin.vars <- apply(covs, 2, is_binary)

      w_mat<- apply(Z, 2, get_cont_weights, treat = Y, s.weights = obsWeights,
                    dens.num = dens.num, densfun = densfun, use.kernel = use.kernel,
                    densControl = densControl)
      cvRisk <- apply(w_mat, 2, function(w) stop.sum(stop.fun(covs, treat = Y,
                                                              weights = w,
                                                              s.weights = obsWeights,
                                                              bin.vars = bin.vars)))
      names(cvRisk) <- libraryNames

      loss <- function(coefs) {
        ps <- crossprod(t(Z), coefs/sum(coefs))
        w <- get_cont_weights(ps, treat = Y, s.weights = obsWeights,
                              dens.num = dens.num, densfun = densfun,
                              use.kernel = use.kernel,
                              densControl = densControl)
        out <- stop.sum(stop.fun(covs, treat = Y,
                                 weights = w,
                                 s.weights = obsWeights,
                                 bin.vars = bin.vars))
        out
      }
      fit <- optim(rep(1, ncol(Z)), loss, method = "L-BFGS-B", lower = 0)
      coef <- fit$par
      out <- list(cvRisk = cvRisk, coef = coef/sum(coef))
      return(out)
    },

    # computePred is a function that takes the weights and the predicted values
    # from each algorithm in the library and combines them based on the model to
    # output the super learner predicted values
    computePred = function(predY, coef, control, ...) {
      out <- crossprod(t(predY), coef/sum(coef))
      return(out)
    }
  )
  # attr(out, "stop.method") <- stop.method
  # class(out) <- "method.balance"
  return(out)
}

.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}

#To pass CRAN checks:
utils::globalVariables(c(".s.weights", "dens", "x", "y"))
