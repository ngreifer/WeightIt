method.to.proper.method <- function(method) {
  method <- tolower(method)
  if (method %in% c("ps")) return("ps")
  else if (method %in% c("gbm", "gbr", "twang")) return("gbm")
  else if (method %in% c("cbps")) return("cbps")
  else if (method %in% c("npcbps")) return("npcbps")
  else if (method %in% c("entropy", "ebal", "ebalance")) return("ebal")
  else if (method %in% c("sbw")) return("sbw")
  else if (method %in% c("ebcw", "ate")) return("ebcw")
  else if (method %in% c("optweight", "opt")) return("optweight")
  else if (method %in% c("super", "superlearner")) return("super")
  else return(method)
}
check.acceptable.method <- function(method, msm = FALSE, force = FALSE) {
  bad.method <- FALSE
  acceptable.methods <- c("ps",
                          "gbm", "twang", "gbr",
                          "cbps",
                          "npcbps",
                          "ebal", "entropy", "ebalance",
                          "sbw",
                          "ebcw", "ate",
                          "optweight", "opt",
                          "super", "superlearner")

  if (missing(method)) method <- "ps"
  else if (is_null(method) || length(method) > 1) bad.method <- TRUE
  else if (is.character(method)) {
    if (tolower(method) %nin% acceptable.methods) bad.method <- TRUE
  }
  else if (!is.function(method)) bad.method <- TRUE

  if (bad.method) stop("method must be a string of length 1 containing the name of an acceptable weighting method or a function that produces weights.", call. = FALSE)

  if (msm && !force && is.character(method)) {
    m <- method.to.proper.method(method)
    if (m %in% c("nbcbps", "ebal", "sbw", "ebcw", "optweight")) {
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
    method <- tolower(method)
    if (method %in% c("ps")) return("propensity score weighting")
    else if (method %in% c("gbm", "gbr", "twang")) return("generalized boosted modeling")
    else if (method %in% c("cbps")) return("covariate balancing propensity score weighting")
    else if (method %in% c("npcbps")) return("non-parametric covariate balancing propensity score weighting")
    else if (method %in% c("entropy", "ebal", "ebalance")) return("entropy balancing")
    else if (method %in% c("sbw")) return("stable balancing weights")
    else if (method %in% c("ebcw", "ate")) return("empirical balancing calibration weighting")
    else if (method %in% c("optweight", "opt")) return("targeted stable balancing weights")
    else if (method %in% c("super", "superlearner")) return("propensity score weighting with SuperLearner")
    else return("the chosen method of weighting")
  }
}
process.estimand <- function(estimand, method, treat.type) {
  #Allowable estimands
  AE <- list(binary = list(ps = c("ATT", "ATC", "ATE", "ATO", "ATM"),
                           gbm = c("ATT", "ATC", "ATE"),
                           cbps = c("ATT", "ATC", "ATE"),
                           npcbps = c("ATE"),
                           ebal = c("ATT", "ATC", "ATE"),
                           sbw = c("ATT", "ATC", "ATE"),
                           ebcw = c("ATT", "ATC", "ATE"),
                           optweight = c("ATT", "ATC", "ATE"),
                           super = c("ATT", "ATC", "ATE", "ATO", "ATM")),
             multinomial = list(ps = c("ATT", "ATE", "ATO"),
                                gbm = c("ATT", "ATE"),
                                cbps = c("ATT", "ATE"),
                                npcbps = c("ATE"),
                                ebal = c("ATT", "ATE"),
                                sbw = c("ATT", "ATE"),
                                ebcw = c("ATT", "ATE"),
                                optweight = c("ATT", "ATE"),
                                super = c("ATT", "ATE")))

  if (treat.type != "continuous" && !is.function(method) &&
      toupper(estimand) %nin% AE[[treat.type]][[method]]) {
    stop(paste0("\"", estimand, "\" is not an allowable estimand for ", method.to.phrase(method),
                " with ", treat.type, " treatments. Only ", word.list(AE[[treat.type]][[method]], quotes = TRUE, and.or = "and", is.are = TRUE),
                " allowed."), call. = FALSE)
  }
  else {
    return(toupper(estimand))
  }
}
process.focal.and.estimand <- function(focal, estimand, treat, treat.type) {
  reported.estimand <- estimand

  #Check focal
  if (treat.type %in% c("binary", "multinomial")) {
    if (estimand == "ATT") {
      if (is_null(focal)) {
        if (treat.type == "multinomial") {
          stop("When estimand = \"ATT\" for multinomial treatments, an argument must be supplied to focal.", call. = FALSE)
        }
      }
      else if (length(focal) > 1L || !any(unique(treat) == focal)) {
        stop("The argument supplied to focal must be the name of a level of treat.", call. = FALSE)
      }
    }
    else {
      if (is_not_null(focal)) {
        warning(paste(estimand, "is not compatible with focal. Ignoring focal."), call. = FALSE)
        focal <- NULL
      }
    }
  }

  #Get focal, estimand, and reported estimand
  if (isTRUE(treat.type == "binary")) {
    unique.treat <- unique(treat, nmax = 2)
    unique.treat.bin <- unique(binarize(treat), nmax = 2)
    if (estimand == "ATT") {
      if (is_null(focal)) {
        focal <- unique.treat[unique.treat.bin == 1]
      }
      else if (focal == unique.treat[unique.treat.bin == 0]){
        reported.estimand <- "ATC"
      }
    }
    else if (estimand == "ATC") {
      focal <- unique.treat[unique.treat.bin == 0]
      estimand <- "ATT"
    }
  }
  return(list(focal = focal,
              estimand = estimand,
              reported.estimand = reported.estimand))
}
process.by <- function(by, data, treat, treat.name = NULL) {

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

  if (bad.by) stop("by must be the quoted names of variables in data for which weighting is to occur within strata or the variable itself.", call. = FALSE)

  if (any(vapply(levels(by.factor), function(x) nunique(treat) != nunique(treat[by.factor == x]), logical(1L)))) {
    stop("Not all the groups formed by by contain all treatment levels", if (is_not_null(treat.name)) paste("in", treat.name) else "", ". Consider coarsening by.", call. = FALSE)
  }

  return(list(by.components = by.components,
              by.factor = by.factor))
}
process.moments.int <- function(moments, int, method) {
  if (!is.function(method)) {
    if (method %in% c("ebal", "ebcw", "sbw", "optweight")) {
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
      warning(paste0(word.list(c("moments", "int")[mi0], and.or = "and", is.are = TRUE),
                     " not compatible with ", method.to.phrase(method), ". Ignoring ", word.list(c("moments", "int")[mi0], and.or = "and"), "."), call. = FALSE)
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
  package.is.installed <- any(.packages(all.available = TRUE) == package.name)
  if (!package.is.installed && !alternative) {
    stop(paste0("Package \"", package.name, "\" needed for this function to work. Please install it."),
         call. = FALSE)
  }
  return(invisible(package.is.installed))
}
make.closer.to.1 <- function(x) {
  if (is.factor(x) || is.character(x)) return(x)
  else if (is_binary(x)) {
    return(as.numeric(x == x[1]))
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
make.full.rank <- function(mat, with.intercept = TRUE) {

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

#To pass CRAN checks:
utils::globalVariables(c(".s.weights"))
