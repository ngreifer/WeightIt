method.to.proper.method <- function(method) {
  method <- tolower(method)
  if      (method %in% c("ps")) return("ps")
  else if (method %in% c("gbm", "gbr", "twang")) return("gbm")
  else if (method %in% c("cbps")) return("cbps")
  else if (method %in% c("npcbps")) return("npcbps")
  else if (method %in% c("entropy", "ebal", "ebalance")) return("ebal")
  else if (method %in% c("ebcw", "ate")) return("ebcw")
  else if (method %in% c("optweight", "opt", "sbw")) return("optweight")
  else if (method %in% c("super", "superlearner")) return("super")
  # else if (method %in% c("kbal")) return("kbal")
  # else if (method %in% c("sbps", "subgroup")) return("sbps")
  else return(method)
}
check.acceptable.method <- function(method, msm = FALSE, force = FALSE) {
  bad.method <- FALSE
  acceptable.methods <- c("ps"
                          , "gbm", "twang", "gbr"
                          , "cbps"
                          , "npcbps"
                          , "ebal", "entropy", "ebalance"
                          , "sbw"
                          , "ebcw", "ate"
                          , "optweight", "opt"
                          , "super", "superlearner"
                          # "kbal",
                          # "sbps", "subgroup"
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
    else if (method %in% c("cbps")) return("covariate balancing propensity score weighting")
    else if (method %in% c("npcbps")) return("non-parametric covariate balancing propensity score weighting")
    else if (method %in% c("ebal")) return("entropy balancing")
    else if (method %in% c("ebcw")) return("empirical balancing calibration weighting")
    else if (method %in% c("optweight")) return("targeted stable balancing weights")
    else if (method %in% c("super")) return("propensity score weighting with SuperLearner")
    # else if (method %in% c("kbal")) return("kernel balancing")
    # else if (method %in% c("sbps")) return("subgroup balancing propensity score weighting")
    else return("the chosen method of weighting")
  }
}
process.estimand <- function(estimand, method, treat.type) {
  #Allowable estimands
  AE <- list(binary = list(ps = c("ATT", "ATC", "ATE", "ATO", "ATM")
                           , gbm = c("ATT", "ATC", "ATE")
                           , cbps = c("ATT", "ATC", "ATE")
                           , npcbps = c("ATE")
                           , ebal = c("ATT", "ATC", "ATE")
                           , ebcw = c("ATT", "ATC", "ATE")
                           , optweight = c("ATT", "ATC", "ATE")
                           , super = c("ATT", "ATC", "ATE", "ATO", "ATM")
                           # , kbal = c("ATT", "ATC", "ATE")
                           # , sbps = c("ATT", "ATC", "ATE", "ATO", "ATM")
                           ),
             multinomial = list(ps = c("ATT", "ATC", "ATE", "ATO", "ATM")
                                , gbm = c("ATT", "ATC", "ATE")
                                , cbps = c("ATT", "ATC", "ATE")
                                , npcbps = c("ATE")
                                , ebal = c("ATT", "ATC", "ATE")
                                , ebcw = c("ATT", "ATC", "ATE")
                                , optweight = c("ATT", "ATC", "ATE")
                                , super = c("ATT", "ATC", "ATE", "ATO", "ATM")
                                # , kbal = c("ATT", "ATE")
                                ))

  if (treat.type != "continuous" && !is.function(method) &&
      toupper(estimand) %nin% AE[[treat.type]][[method]]) {
    stop(paste0("\"", estimand, "\" is not an allowable estimand for ", method.to.phrase(method),
                " with ", treat.type, " treatments. Only ", word_list(AE[[treat.type]][[method]], quotes = TRUE, and.or = "and", is.are = TRUE),
                " allowed."), call. = FALSE)
  }
  else {
    return(toupper(estimand))
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
  acceptable.bys <- names(data)
  by.vars <- character(0)
  by.components <- NULL
  n <- length(treat)

  if (missing(by) || is_null(by)) by.factor <- factor(rep(1, n))
  else if (!is.atomic(by)) bad.by <- TRUE
  else if (is.character(by) && all(by %in% acceptable.bys)) {
    by.components <- data[by]
    by.factor <- factor(apply(by.components, 1, paste, collapse = "|"))
    by.vars <- by
  }
  else if (length(by) == n) {
    by.components <- setNames(data.frame(by), deparse(substitute(by)))
    by.factor <- factor(by.components[[1]])
    by.vars <- acceptable.bys[vapply(acceptable.bys, function(x) equivalent.factors(by, data[[x]]), logical(1L))]
  }
  else bad.by <- TRUE

  if (bad.by) stop(paste(by.arg, "must be the quoted names of variables in data for which weighting is to occur within strata or the variable itself."), call. = FALSE)

  if (any(vapply(levels(by.factor), function(x) nunique(treat) != nunique(treat[by.factor == x]), logical(1L)))) {
    stop(paste0("Not all the groups formed by ", by.arg, " contain all treatment levels", if (is_not_null(treat.name)) paste("in", treat.name) else "", ". Consider coarsening", by.arg, "."), call. = FALSE)
  }

  return(list(by.components = by.components,
              by.factor = by.factor))
}
process.moments.int <- function(moments, int, method) {
  if (!is.function(method)) {
    if (method %in% c("ebal", "ebcw", "optweight")) {
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
check.package <- function(package.name, alternative = FALSE) {
  package.not.installed <- package.name[package.name %nin% .packages(all.available = TRUE)]
  if (is_not_null(package.not.installed) && !alternative) {
    plural <- length(package.not.installed) > 1
    stop(paste0("Package", if (plural) "s " else " ",
                word_list(package.not.installed, quotes = TRUE, is.are = TRUE),
                " needed for this function to work. Please install ",
                if (plural) "them" else "it","."),
         call. = FALSE)
  }
  else return(invisible(TRUE))
}
make.closer.to.1 <- function(x) {
  if (is.factor(x) || is.character(x)) return(x)
  else if (is_binary(x)) {
    return(as.numeric(x == x[!is.na(x)][1]))
  }
  else {
    ndigits <- round(mean(floor(log10(abs(x[!check_if_zero(x)]))), na.rm = TRUE))
    if (abs(ndigits) > 2) return(x/(10^ndigits))
    else return(x)
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

#For the user to use
make_full_rank <- function(mat, with.intercept = TRUE) {

  if (is.data.frame(mat)) {
    is.mat <- FALSE
    if (!all(vapply(mat, is.numeric, logical(1L)))) stop("All columns in mat must be numeric.", call. = FALSE)
    mat <- as.matrix(mat)
  }
  else if (is.matrix(mat)) {
    is.mat <- TRUE
    if (!is.numeric(mat)) stop("mat must be a numeric matrix.", call. = FALSE)
  }
  else {
    stop("mat must be a numeric matrix or data.frame.", call. = FALSE)
  }

  if (anyNA(mat)) stop("Missing values are not allowed in mat.", call. = FALSE)

  keep <- setNames(rep(TRUE, ncol(mat)), colnames(mat))

  #Variables that have only 1 value can be removed
  all.the.same <- apply(mat, 2, all_the_same)
  keep[all.the.same] <- FALSE

  #If intercept is to be included in check, add column of 1s
  if (with.intercept) mat1 <- cbind(mat, rep(1, nrow(mat)))
  else mat1 <- mat

  for (i in colnames(mat)[keep]) {
    #Add extra value for intercept if desired
    keep1 <- c(keep, TRUE[with.intercept])

    #Create vector of keep with ith entry FALSE to compare rank with full vector
    keep1. <- keep1
    keep1[i] <- FALSE

    #Check if rank without is the same as rank with; if so, remove variable i
    if (qr(mat1[, keep1., drop = FALSE])$rank == qr(mat1[, keep1, drop = FALSE])$rank) {
      keep[i] <- FALSE
    }
  }

  if (is.mat) return(mat[, keep, drop = FALSE])
  else return(as.data.frame(mat[, keep, drop = FALSE]))

}
get_w_from_ps <- function(ps, treat, estimand = "ATE", focal = NULL, treated = NULL) {
  #ps must be a matrix/df with columns named after treat levels

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)

  treat.type <- get.treat.type(treat)

  processed.estimand <- process.focal.and.estimand(focal, estimand, treat, treat.type, treated)
  estimand <- processed.estimand$estimand
  focal <- processed.estimand$focal
  assumed.treated <- processed.estimand$treated

  if (treat.type == "continuous") {
    stop("get_w_from_ps can only be used with binary or multinomial treatments.", call. = FALSE)
  }
  else {
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
          if (suppressWarnings(!any(is.na(as.numeric(as.character(treat))))) &&
            all(check_if_zero(binarize(treat) - as.numeric(as.character(treat))))) treated.level <- 1
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
  }

  if (nrow(ps) != length(treat)) {
    stop("ps and treat must have the same number of units.", call. = FALSE)
  }

  w <- setNames(rep(0, nrow(ps)), ps.names)

  for (i in t.levels) {
    w[treat == i] <- 1/ps[treat == i, as.character(i)]
  }

  if (toupper(estimand) == "ATE") {
    # w <- w
  }
  else if (toupper(estimand) == "ATT") {
    w <- w*ps[, as.character(focal)]
  }
  else if (toupper(estimand) == "ATO") {
    w <- w/apply(ps, 1, function(x) sum(1/x)) #Li & Li (2018)
  }
  else if (toupper(estimand) == "ATM") {
    w <- w*apply(ps, 1, min)
  }
  else w <- NULL

  return(w)
}

#To pass CRAN checks:
utils::globalVariables(c(".s.weights", "dens", "x", "y"))
