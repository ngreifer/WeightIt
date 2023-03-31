#Turn a vector into a 0/1 vector. 'zero' and 'one' can be supplied to make it clear which is
#which; otherwise, a guess is used.
binarize <- function(variable, zero = NULL, one = NULL) {
  var.name <- deparse1(substitute(variable))
  if (length(unique(variable)) > 2) {
    stop(sprintf("Cannot binarize %s: more than two levels.", var.name), call. = FALSE)
  }
  if (is.character(variable) || is.factor(variable)) {
    variable <- factor(variable, nmax = 2)
    unique.vals <- levels(variable)
  }
  else {
    unique.vals <- unique(variable, nmax = 2)
  }

  if (is.null(zero)) {
    if (is.null(one)) {
      if (can_str2num(unique.vals)) {
        variable.numeric <- str2num(variable)
      }
      else {
        variable.numeric <- as.numeric(variable)
      }

      if (0 %in% variable.numeric) zero <- 0
      else zero <- min(variable.numeric, na.rm = TRUE)

      out <- setNames(as.integer(variable.numeric != zero), names(variable))
    }
    else {
      if (one %in% unique.vals) out <- setNames(as.integer(variable == one), names(variable))
      else stop("The argument to 'one' is not the name of a level of variable.", call. = FALSE)
    }
  }
  else {
    if (zero %in% unique.vals) out <- setNames(as.integer(variable != zero), names(variable))
    else stop("The argument to 'zero' is not the name of a level of variable.", call. = FALSE)
  }

  return(out)
}

#Get covariates (RHS) vars from formula
get.covs.matrix <- function(formula = NULL, data = NULL) {

  if (is.null(formula)) {
    fnames <- colnames(data)
    fnames[!startsWith(fnames, "`")] <- paste0("`", fnames[!startsWith(fnames, "`")], "`")
    formula <- reformulate(fnames)
  }
  else formula <- update(terms(formula, data = data), NULL ~ . + 1)

  mf <- model.frame(terms(formula, data = data), data,
                    na.action = na.pass)

  chars.in.mf <- vapply(mf, is.character, logical(1L))
  mf[chars.in.mf] <- lapply(mf[chars.in.mf], factor)

  X <- model.matrix(formula, data = mf,
                    contrasts.arg = lapply(Filter(is.factor, mf),
                                           contrasts, contrasts = FALSE))
  assign <- attr(X, "assign")[-1]
  X <- X[,-1, drop = FALSE]
  attr(X, "assign") <- assign

  return(X)
}

#Add missing indicators (is.na(x)) to formula
add_miss_ind_to_formula <- function(formula, data = NULL) {
  mf <- model.frame(delete.response(terms(formula, data = data)),
                    data = data, na.action = "na.pass")
  missing_vars <- names(mf)[vapply(mf, anyNA, logical(1L))]

  s <- do.call("paste",
               c(list(". ~ ."),
                 lapply(missing_vars,
                        sprintf,
                        fmt = "+ is.na(%s)")))

  update(formula, s)
}

#treat.type processing
assign.treat.type <- function(treat, use.multi = FALSE) {
  #Returns treat with treat.type attribute
  nunique.treat <- nunique(treat)

  if (nunique.treat < 2) {
    stop("The treatment must have at least two unique values.", call. = FALSE)
  }
  else if (!use.multi && nunique.treat == 2) {
    treat.type <- "binary"
  }
  else if (use.multi || is_(treat, c("factor", "character"))) {
    treat.type <- "multinomial"
    if (!is_(treat, "processed.treat")) treat <- factor(treat)
  }
  else {
    treat.type <- "continuous"
  }
  attr(treat, "treat.type") <- treat.type
  return(treat)
}
get.treat.type <- function(treat) {
  return(attr(treat, "treat.type"))
}
has.treat.type <- function(treat) {
  is_not_null(get.treat.type(treat))
}
get.treated.level <- function(treat) {
  if (!is_binary(treat)) stop("'treat' must be a binary variable.")
  if (is.character(treat) || is.factor(treat)) {
    treat <- factor(treat, nmax = 2)
    unique.vals <- levels(treat)
  }
  else {
    unique.vals <- unique(treat, nmax = 2)
  }

  if (can_str2num(unique.vals)) {
    unique.vals.numeric <- str2num(unique.vals)
  }
  else {
    unique.vals.numeric <- seq_along(unique.vals)
  }

  if (0 %in% unique.vals.numeric) treated <- unique.vals[unique.vals.numeric != 0]
  else treated <- unique.vals[which.max(unique.vals.numeric)]

  return(treated)
}

#Converting string to numeric
can_str2num <- function(x) {
  if (is.numeric(x) || is.logical(x)) return(TRUE)
  nas <- is.na(x)
  suppressWarnings(x_num <- as.numeric(as.character(x[!nas])))
  return(!anyNA(x_num))
}
str2num <- function(x) {
  nas <- is.na(x)
  if (!is_(x, c("numeric", "logical"))) x <- as.character(x)
  suppressWarnings(x_num <- as.numeric(x))
  is.na(x_num[nas]) <- TRUE
  return(x_num)
}

#Weighted mean faster than weighted.mean()
w.m <- function(x, w = NULL, na.rm = TRUE) {
  if (is.null(w)) {
    if (anyNA(x)) {
      if (!na.rm) return(NA_real_)
      nas <- which(is.na(x))
      x <- x[-nas]
    }
    return(sum(x)/length(x))
  }
  else {
    if (anyNA(x) || anyNA(w)) {
      if (!na.rm) return(NA_real_)
      nas <- which(is.na(x) | is.na(w))
      x <- x[-nas]
      w <- w[-nas]
    }
    return(sum(x*w)/sum(w))
  }
}

#Input processing
method.to.proper.method <- function(method) {
  method <- tolower(method)
  if      (method %in% c("ps")) return("ps")
  else if (method %in% c("gbm", "gbr")) return("gbm")
  else if (method %in% c("cbps", "cbgps")) return("cbps")
  else if (method %in% c("npcbps", "npcbgps")) return("npcbps")
  else if (method %in% c("entropy", "ebal", "ebalance")) return("ebal")
  else if (method %in% c("ebcw", "ate")) return("ebcw")
  else if (method %in% c("optweight", "opt", "sbw")) return("optweight")
  else if (method %in% c("super", "superlearner")) return("super")
  else if (method %in% c("bart")) return("bart")
  else if (method %in% c("energy")) return("energy")
  # else if (method %in% c("kbal")) return("kbal")
  else return(method)
}
check.acceptable.method <- function(method, msm = FALSE, force = FALSE) {
  bad.method <- FALSE
  acceptable.methods <- c("ps"
                          , "gbm", "gbr"
                          , "cbps", "cbgps"
                          , "npcbps", "npcbgps"
                          , "ebal", "entropy", "ebalance"
                          , "sbw"
                          , "ebcw", "ate"
                          , "optweight", "opt", "sbw"
                          , "super", "superlearner"
                          , "energy"
                          , "bart"
                          # "kbal",
  )

  if (missing(method)) method <- "ps"
  else if (is_null(method) || length(method) > 1) bad.method <- TRUE
  else if (is.character(method)) {
    if (tolower(method) %nin% acceptable.methods) bad.method <- TRUE
  }
  else if (!is.function(method)) bad.method <- TRUE

  if (bad.method) {
    if (identical(method, "twang")) stop('"twang" is no longer an acceptable argument to \'method\'. Please use "gmb" for generalized boosted modeling.', call. = FALSE)
    stop("'method' must be a string of length 1 containing the name of an acceptable weighting method or a function that produces weights.", call. = FALSE)
  }

  if (msm && !force && is.character(method)) {
    m <- method.to.proper.method(method)
    if (m %in% c("nbcbps", "ebal", "ebcw", "optweight", "energy", "kbal")) {
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
    stop("The user-provided function to 'method' must contain \"covs\" and \"treat\" as named parameters.", call. = FALSE)
  }
}
method.to.phrase <- function(method) {

  if (is.function(method)) return("a user-defined method")
  else {
    method <- method.to.proper.method(method)
    if (method %in% c("ps")) return("propensity score weighting")
    else if (method %in% c("gbm")) return("propensity score weighting with GBM")
    else if (method %in% c("cbps")) return("covariate balancing propensity score weighting")
    else if (method %in% c("npcbps")) return("non-parametric covariate balancing propensity score weighting")
    else if (method %in% c("ebal")) return("entropy balancing")
    else if (method %in% c("ebcw")) return("empirical balancing calibration weighting")
    else if (method %in% c("optweight")) return("targeted stable balancing weights")
    else if (method %in% c("super")) return("propensity score weighting with SuperLearner")
    else if (method %in% c("bart")) return("propensity score weighting with BART")
    else if (method %in% c("energy")) return("energy balancing")
    # else if (method %in% c("kbal")) return("kernel balancing")
    else return("the chosen method of weighting")
  }
}
process.estimand <- function(estimand, method, treat.type) {
  #Allowable estimands
  AE <- list(
    binary = list(ps = c("ATT", "ATC", "ATE", "ATO", "ATM", "ATOS")
                  , gbm = c("ATT", "ATC", "ATE", "ATO", "ATM")
                  , cbps = c("ATT", "ATC", "ATE")
                  , npcbps = c("ATE")
                  , ebal = c("ATT", "ATC", "ATE")
                  , ebcw = c("ATT", "ATC", "ATE")
                  , optweight = c("ATT", "ATC", "ATE")
                  , super = c("ATT", "ATC", "ATE", "ATO", "ATM")
                  , energy = c("ATT", "ATC", "ATE")
                  , bart = c("ATT", "ATC", "ATE", "ATO", "ATM")
                  # , kbal = c("ATT", "ATC", "ATE")
    ),
    multinomial = list(ps = c("ATT", "ATC", "ATE", "ATO", "ATM")
                       , gbm = c("ATT", "ATC", "ATE", "ATO", "ATM")
                       , cbps = c("ATT", "ATC", "ATE")
                       , npcbps = c("ATE")
                       , ebal = c("ATT", "ATC", "ATE")
                       , ebcw = c("ATT", "ATC", "ATE")
                       , optweight = c("ATT", "ATC", "ATE")
                       , super = c("ATT", "ATC", "ATE", "ATO", "ATM")
                       , energy = c("ATT", "ATC", "ATE")
                       , bart = c("ATT", "ATC", "ATE", "ATO", "ATM")
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
                  , cbps = TRUE
                  , npcbps = FALSE
                  , ebal = FALSE
                  , ebcw = FALSE
                  , optweight = FALSE
                  , super = TRUE
                  , energy = FALSE
                  , bart = TRUE
                  # , kbal = FALSE
    ),
    multinomial = list(ps = TRUE
                       , gbm = TRUE
                       , cbps = FALSE
                       , npcbps = FALSE
                       , ebal = FALSE
                       , ebcw = FALSE
                       , optweight = FALSE
                       , super = TRUE
                       , energy = FALSE
                       , bart = TRUE
    ))

  if (treat.type != "continuous" && !is.function(method) &&
      !AE[[treat.type]][[method]]) {
    stop(paste0("subclasses are not compatible with ", method.to.phrase(method),
                " with ", treat.type, " treatments."), call. = FALSE)
  }
}
process.ps <- function(ps, data = NULL, treat) {
  if (is_not_null(ps)) {
    if (is.character(ps) && length(ps)==1) {
      if (is_null(data)) {
        stop("'ps' was specified as a string but there was no argument to 'data'.", call. = FALSE)
      }
      else if (ps %in% names(data)) {
        ps <- data[[ps]]
      }
      else stop("The name supplied to 'ps' is not the name of a variable in 'ps'.", call. = FALSE)
    }
    else if (is.numeric(ps)) {
      if (length(ps) != length(treat)) {
        stop("'ps' must have the same number of units as the treatment.", call. = FALSE)
      }
    }
    else {
      stop("The argument to 'ps' must be a vector of propensity scores or the (quoted) names of the variable in 'data' that contains sampling weights.", call. = FALSE)
    }
  }
  else ps <- NULL
  return(ps)
}
process.focal.and.estimand <- function(focal, estimand, treat, treated = NULL) {
  reported.estimand <- estimand

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
  treat.type <- get.treat.type(treat)

  unique.treat <- unique(treat, nmax = switch(treat.type, "binary" = 2, "multinomial" = length(treat)/4))

  #Check focal
  if (is_not_null(focal) && (length(focal) > 1L || focal %nin% unique.treat)) {
    stop("The argument supplied to 'focal' must be the name of a level of treatment.", call. = FALSE)
  }

  if (treat.type == "multinomial") {

    if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
      warning(paste(estimand, "is not compatible with 'focal'. Setting 'estimand' to \"ATT\"."), call. = FALSE, immediate. = TRUE)
      reported.estimand <- estimand <- "ATT"
    }

    if (estimand == "ATT") {
      if (is_null(focal)) {
        if (is_null(treated) || treated %nin% unique.treat) {
          stop("When estimand = \"ATT\" for multinomial treatments, an argument must be supplied to 'focal'.", call. = FALSE)
        }
        focal <- treated
      }
    }
    else if (estimand == "ATC") {
      if (is_null(focal)) {
        stop("When estimand = \"ATC\" for multinomial treatments, an argument must be supplied to 'focal'.", call. = FALSE)
      }
    }
  }
  else if (treat.type == "binary") {
    unique.treat.bin <- unique(binarize(treat), nmax = 2)
    if (estimand %nin% c("ATT", "ATC") && is_not_null(focal)) {
      warning(paste(estimand, "is not compatible with 'focal'. Setting 'estimand' to \"ATT\"."), call. = FALSE, immediate. = TRUE)
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
            message(paste0("Assuming ", word_list(treated, quotes = !is.numeric(treat), is.are = TRUE),
                           " the treated level. If not, supply an argument to 'focal'."))

          }
          else if (estimand == "ATC") {
            message(paste0("Assuming ", word_list(unique.treat[unique.treat %nin% treated], quotes = !is.numeric(treat), is.are = TRUE),
                           " the control level. If not, supply an argument to 'focal'."))
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
    if (NCOL(by) != 1) stop(paste0("Only one variable can be on the right-hand side of the formula for '", by.arg, "'."), call. = FALSE)
    else by.name <- colnames(by)
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
    stop(paste0("'",by.arg, "' must be a string containing the name of the variable in data for which weighting is to occur within strata or a one-sided formula with the stratifying variable on the right-hand side."), call. = FALSE)
  }

  if (treat.type != "continuous" && any(vapply(levels(by.factor), function(x) nunique(treat) != nunique(treat[by.factor == x]), logical(1L)))) {
    stop(paste0("Not all the groups formed by '", by.arg, "' contain all treatment levels", if (is_not_null(treat.name)) paste("in", treat.name) else "", ". Consider coarsening ", by.arg, "."), call. = FALSE)
  }
  attr(by.components, "by.factor") <- by.factor
  return(by.components)
}
process.moments.int <- function(moments, int, method) {

  if (is.function(method) || method %in% c("npcbps", "ebal", "ebcw", "optweight", "energy")) {
    if (length(int) != 1 || !is.logical(int)) {
      stop("int must be a logical (TRUE/FALSE) of length 1.", call. = FALSE)
    }
    if (is_not_null(moments)) {
      if (length(moments) != 1 || !is.numeric(moments) ||
          !check_if_zero(moments - round(moments))) {
        if (method == "energy") {
          if (moments < 0) stop("'moments' must be a nonnegative integer of length 1.", call. = FALSE)
        }
        else if (method %in% c("npcbps", "ebal", "ebcw", "optweight")) {
          if (moments < 1) stop("'moments' must be a positive integer of length 1.", call. = FALSE)
        }
        moments <- as.integer(moments)
      }
    }
    else {
      if (!is.function(method) && method == "energy") moments <- 0L
      else moments <- 1L
    }
  }
  else if (is_not_null(moments) && any(mi0 <- c(as.integer(moments) != 1L, int))) {
    warning(paste0(word_list(c("moments", "int")[mi0], and.or = "and", is.are = TRUE, quotes = 1),
                   " not compatible with ", method.to.phrase(method), ". Ignoring ", word_list(c("moments", "int")[mi0], and.or = "and", quotes = 1), "."), call. = FALSE, immediate. = TRUE)
    moments <- NULL
    int <- FALSE
  }
  moments <- as.integer(moments)

  return(list(moments = moments, int = int))
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
      message(paste0(method.to.phrase(method), " can be used with a single model when multiple time points are present.\nUsing a seperate model for each time point. To use a single model, set 'is.MSM.method' to TRUE."))
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
  AE <- list(binary = list(ps = c("ind"
                                  , "saem"
  )
  , gbm = c("ind", "surr")
  , cbps = c("ind")
  , npcbps = c("ind")
  , ebal = c("ind")
  , ebcw = c("ind")
  , optweight = c("ind")
  , super = c("ind")
  , bart = c("ind")
  , energy = c("ind")
  # , kbal = c("ind")
  ),
  multinomial = list(ps = c("ind")
                     , gbm = c("ind", "surr")
                     , cbps = c("ind")
                     , npcbps = c("ind")
                     , ebal = c("ind")
                     , ebcw = c("ind")
                     , optweight = c("ind")
                     , super = c("ind")
                     , bart = c("ind")
                     , energy = c("ind")
                     # , kbal = c("ind")
  ),
  continuous = list(ps = c("ind"
                           , "saem"
  )
  , gbm = c("ind", "surr")
  , cbps = c("ind")
  , npcbps = c("ind")
  , ebal = c("ind")
  , ebcw = c("ind")
  , optweight = c("ind")
  , super = c("ind")
  , bart = c("ind")
  , energy = c("ind")
  # , kbal = c("ind")
  ))

  allowable.missings <- AE[[treat.type]][[method]]
  if (is_null(missing)) {
    missing <- allowable.missings[1]
    warning(paste0("Missing values are present in the covariates. See ?WeightIt::method_",
                   method, " for information on how these are handled."), call. = FALSE, immediate. = TRUE)
  }
  else {
    if (!is.character(missing) || length(missing) != 1) stop("'missing' must be a string of length 1.", call. = FALSE)
    if (missing %pin% allowable.missings) {
      missing <- allowable.missings[pmatch(missing, allowable.missings)]
    }
    else {
      missing <- allowable.missings[1]
      warning(paste0("Only ", word_list(allowable.missings, quotes = 2, is.are = TRUE), " allowed for 'missing' with ",
                     treat.type,
                     " treatments. Using link = ", word_list(allowable.missings[1], quotes = 2), "."),
              call. = FALSE, immediate. = TRUE)
    }
  }
  return(missing)
}
process.bin.vars <- function(bin.vars, mat) {
  if (missing(bin.vars)) bin.vars <- is_binary_col(mat)
  else if (is_null(bin.vars)) bin.vars <- rep(FALSE, ncol(mat))
  else {
    if (is.logical(bin.vars)) {
      bin.vars[is.na(bin.vars)] <- FALSE
      if (length(bin.vars) != ncol(mat)) stop("If 'bin.vars' is logical, it must have length equal to the number of columns of 'mat'.")
    }
    else if (is.numeric(bin.vars)) {
      bin.vars <- bin.vars[!is.na(bin.vars) & bin.vars != 0]
      if (any(bin.vars < 0) && any(bin.vars > 0)) stop("Positive and negative indices cannot be mixed with 'bin.vars'.")
      if (any(abs(bin.vars) > ncol(mat))) stop("If 'bin.vars' is numeric, none of its values can exceed the number of columns of 'mat'.")
      logical.bin.vars <- rep(any(bin.vars < 0), ncol(mat))
      logical.bin.vars[abs(bin.vars)] <- !logical.bin.vars[abs(bin.vars)]
      bin.vars <- logical.bin.vars
    }
    else if (is.character(bin.vars)) {
      bin.vars <- bin.vars[!is.na(bin.vars) & bin.vars != ""]
      if (is_null(colnames(mat))) stop("If 'bin.vars' is character, 'mat' must have column names.")
      if (any(bin.vars %nin% colnames(mat))) stop("If 'bin.vars' is character, all its values must be column names of 'mat'.")
      bin.vars <- colnames(mat) %in% bin.vars
    }
    else stop("'bin.vars' must be a logical, numeric, or character vector.")
  }
  return(bin.vars)
}

int.poly.f <- function(d, ex = NULL, int = FALSE, poly = 1, center = TRUE, orthogonal_poly = TRUE) {
  #Adds to data frame interactions and polynomial terms
  #d=matrix input
  #ex=names of variables to exclude in interactions and polynomials; a subset of df
  #int=whether to include interactions or not; currently only 2-way are supported
  #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included

  if (is_null(ex)) ex <- rep(FALSE, ncol(d))

  binary.vars <- is_binary_col(d)

  if (center) {
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
  else out <- NULL

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

  return(s.d.denom)
}
get.s.d.denom.cont.weightit <- function(s.d.denom = NULL) {
  s.d.denom.specified <- is_not_null(s.d.denom)

  if (s.d.denom.specified) {
    allowable.s.d.denoms <- c("all", "weighted")

    try.s.d.denom <- tryCatch(match_arg(s.d.denom, allowable.s.d.denoms),
                              error = function(cond) NA_character_)
    if (anyNA(try.s.d.denom)) {
      s.d.denom <- "all"
    }
    else {
      s.d.denom <- try.s.d.denom
    }
  }
  else {
    s.d.denom <- "all"
  }

  return(s.d.denom)
}
compute_s.d.denom <- function(mat, treat, s.d.denom = "pooled", s.weights = NULL, bin.vars = NULL, subset = NULL, weighted.weights = NULL, to.sd = rep(TRUE, ncol(mat)), na.rm = TRUE) {
  denoms <- setNames(rep(1, ncol(mat)), colnames(mat))
  if (is.character(s.d.denom) && length(s.d.denom) == 1L) {
    if (is_null(bin.vars)) {
      bin.vars <- rep(FALSE, ncol(mat))
      bin.vars[to.sd] <- is_binary_col(mat[subset, to.sd,drop = FALSE])
    }
    else if (!is.atomic(bin.vars) || length(bin.vars) != ncol(mat) ||
             anyNA(as.logical(bin.vars))) {
      stop("'bin.vars' must be a logical vector with length equal to the number of columns of 'mat'.")
    }

    possibly.supplied <- c("mat", "treat", "weighted.weights", "s.weights", "subset")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
      stop(paste(word_list(possibly.supplied[supplied], quotes = 1), "must have the same number of units."))
    }

    if (lengths["weighted.weights"] == 0) weighted.weights <- rep(1, NROW(mat))
    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
    if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
    else if (anyNA(as.logical(subset))) stop("'subset' must be a logical vector.")

    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    cont.treat <- get.treat.type(treat) == "continuous"

    if (!cont.treat) {
      treat <- as.character(treat)
      unique.treats <- unique(treat)
    }
    else unique.treats <- NULL

    if (s.d.denom %in% unique.treats)
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(col.w.v(mat[treat == s.d.denom, , drop = FALSE],
                     w = s.weights[treat == s.d.denom],
                     bin.vars = bin.vars, na.rm = na.rm))
      }

    else if (s.d.denom == "pooled")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(Reduce("+", lapply(unique.treats,
                                function(t) col.w.v(mat[treat == t, , drop = FALSE],
                                                    w = s.weights[treat == t],
                                                    bin.vars = bin.vars, na.rm = na.rm))) / length(unique.treats))
      }
    else if (s.d.denom == "all")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(col.w.v(mat, w = s.weights, bin.vars = bin.vars, na.rm = na.rm))
      }
    else if (s.d.denom == "weighted")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        sqrt(col.w.v(mat, w = weighted.weights * s.weights, bin.vars = bin.vars, na.rm = na.rm))
      }
    else if (s.d.denom == "hedges")
      denom.fun <- function(mat, treat, s.weights, weighted.weights, bin.vars,
                            unique.treats, na.rm) {
        (1 - 3/(4*length(treat) - 9))^-1 * sqrt(Reduce("+", lapply(unique.treats,
                                                                   function(t) (sum(treat == t) - 1) * col.w.v(mat[treat == t, , drop = FALSE],
                                                                                                               w = s.weights[treat == t],
                                                                                                               bin.vars = bin.vars, na.rm = na.rm))) / (length(treat) - 2))
      }
    else stop("s.d.denom is not an allowed value.")

    denoms[to.sd] <- denom.fun(mat = mat[, to.sd, drop = FALSE], treat = treat, s.weights = s.weights,
                               weighted.weights = weighted.weights, bin.vars = bin.vars[to.sd],
                               unique.treats = unique.treats, na.rm = na.rm)

    if (any(zero_sds <- check_if_zero(denoms[to.sd]))) {
      denoms[to.sd][zero_sds] <- sqrt(col.w.v(mat[, to.sd, drop = FALSE][, zero_sds, drop = FALSE],
                                              w = s.weights,
                                              bin.vars = bin.vars[to.sd][zero_sds], na.rm = na.rm))
    }

    if (cont.treat) {
      treat.sd <- denom.fun(mat = treat, s.weights = s.weights,
                            weighted.weights = weighted.weights, bin.vars = FALSE,
                            na.rm = na.rm)
      denoms[to.sd] <- denoms[to.sd]*treat.sd
    }
  }
  else {
    if (is.numeric(s.d.denom)) {
      if (is_not_null(names(s.d.denom)) && any(colnames(mat) %in% names(s.d.denom))) {
        denoms[colnames(mat)[colnames(mat) %in% names(s.d.denom)]] <- s.d.denom[names(s.d.denom)[names(s.d.denom) %in% colnames(mat)]]
      }
      else if (length(s.d.denom) == sum(to.sd)) {
        denoms[to.sd] <- s.d.denom
      }
      else if (length(s.d.denom) == ncol(mat)) {
        denoms[] <- s.d.denom
      }
      else {
        stop("'s.d.denom' must be an allowable value or a numeric vector of with length equal to the number of columns of 'mat'. See ?cobalt::col_w_smd for allowable values.")
      }
    }
    else {
      stop("'s.d.denom' must be an allowable value or a numeric vector of with length equal to the number of columns of 'mat'. See ?cobalt::col_w_smd for allowable values.")
    }
  }
  return(denoms)
}
check_estimated_weights <- function(w, treat, treat.type, s.weights) {

  tw <- w*s.weights

  extreme.warn <- FALSE
  if (treat.type == "continuous") {
    if (all_the_same(w)) {
      warning(paste0("All weights are ", w[1], ", possibly indicating an estimation failure."), call. = FALSE)
    }
    else if (sd(tw, na.rm = TRUE)/mean(tw, na.rm = TRUE) > 4) extreme.warn <- TRUE
  }
  else {
    if (all_the_same(w)) {
      warning(paste0("All weights are ", w[1], ", possibly indicating an estimation failure."), call. = FALSE)
    }
    else {
      t.levels <- unique(treat)
      bad.treat.groups <- setNames(rep(FALSE, length(t.levels)), t.levels)
      for (i in t.levels) {
        ti <- which(treat == i)
        if (all(is.na(w[ti])) || all(w[ti] == 0)) bad.treat.groups[as.character(i)] <- TRUE
        else if (!extreme.warn && sum(!is.na(tw[ti])) > 1 && sd(tw[ti], na.rm = TRUE)/mean(tw[ti], na.rm = TRUE) > 4) extreme.warn <- TRUE
      }

      if (any(bad.treat.groups)) {
        n <- sum(bad.treat.groups)
        warning(paste0("All weights are NA or 0 in treatment ", ngettext(n, "group ", "groups "),
                       word_list(t.levels[bad.treat.groups], quotes = TRUE), "."), call. = FALSE)
      }
    }
  }

  if (extreme.warn) warning("Some extreme weights were generated. Examine them with summary() and maybe trim them with trim().", call. = FALSE)

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
    if ((is.matrix(ps) && is.numeric(ps)) ||
        (is.data.frame(ps) && all(vapply(ps, is.numeric, logical(1L))))) {
      if (ncol(ps) == 1) {
        if (can_str2num(treat) &&
            all(check_if_zero(binarize(treat) - str2num(treat)))) treated.level <- 1
        else if (is_not_null(treated)) {
          if (treated %in% treat) treated.level <- treated
          else stop("The argument to 'treated' must be a value in treat.", call. = FALSE)
        }
        else if (is_not_null(assumed.treated)) {
          treated.level <- assumed.treated
        }
        else if (is_not_null(colnames(ps)) && colnames(ps) %in% as.character(t.levels)) {
          treated.level <- colnames(ps)
        }
        else {
          stop("If the treatment has two non-0/1 levels and 'ps' is a vector or has only one column, an argument to 'treated' must be supplied.", call. = FALSE)
        }

        t.levels <- c(treated.level, t.levels[t.levels != treated.level])
        ps <- matrix(c(ps[,1], 1-ps[,1]), ncol = 2, dimnames = list(ps.names, as.character(t.levels)))
      }
      else if (ncol(ps) == 2) {
        if (!all(as.character(t.levels) %in% colnames(ps))) {
          stop("If 'ps' has two columns, they must be named with the treatment levels.", call. = FALSE)
        }
        # else ps <- ps
      }
      else {
        stop("'ps' cannot have more than two columns if the treatment is binary.", call. = FALSE)
      }
    }
    else {
      stop("'ps' must be a matrix, data frame, or vector of propensity scores.", call. = FALSE)
    }
  }
  else if (treat.type == "multinomial") {
    if ((is.matrix(ps) && is.numeric(ps)) ||
        (is.data.frame(ps) && all(vapply(ps, is.numeric, logical(1L))))) {
      if (ncol(ps) == 1) {
        if (toupper(estimand) == "ATE") {
          ps <- matrix(rep(ps, nunique(treat)), nrow = length(treat), dimnames = list(ps.names, t.levels))
        }
        else {
          stop("With multinomial treatments, 'ps' can be a vector or have only one column only if the estimand is the ATE.", call. = FALSE)
        }
      }
      else if (ncol(ps) == nunique(treat)) {
        if (!all(t.levels %in% colnames(ps))) {
          stop("The columns of 'ps' must be named with the treatment levels.", call. = FALSE)
        }
        else ps <- ps
      }
      else {
        stop("'ps' must have as many columns as there are treatment levels.", call. = FALSE)
      }
    }
    else {
      stop("'ps' must be a matrix or data frame of propensity scores.", call. = FALSE)
    }
  }
  return(ps)
}
subclass_ps <- function(ps_mat, treat, estimand = "ATE", focal = NULL, subclass) {
  if (length(subclass) != 1 || !is.numeric(subclass)) {
    stop("'subclass' must be a single number.", call. = FALSE)
  }
  else if (round(subclass) <= 1) {
    stop("'subclass' must be greater than 1.", call. = FALSE)
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
  return(ps_sub)
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

  return(sub)
}
stabilize_w <- function(weights, treat) {
  if (is.factor(treat)) t.levels <- levels(treat)
  else t.levels <- unique(treat)
  w.names <- names(weights)
  tab <- setNames(vapply(t.levels, function(x) mean_fast(treat == x), numeric(1L)), t.levels)
  return(setNames(weights * tab[as.character(treat)], w.names))
}
`%+%` <- function(...) {
  if (is_(..1, "atomic") && is_(..2, "atomic")) crayon::`%+%`(as.character(..1), as.character(..2))
  else ggplot2::`%+%`(...)
}