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
check.user.method <- function(method) {
  #Check to make sure it accepts treat and covs
  if (all(c("covs", "treat") %in% names(formals(method)))) {
    #attr(method, "is.MSM.method") <- FALSE
  }
  else if (all(c("covs.list", "treat.list") %in% names(formals(method)))) {
    #attr(method, "is.MSM.method") <- TRUE
  }
  else {
    stop("The user-provided function to method must contain \"covs\" and \"treat\" as named parameters.", call. = FALSE)
  }
  return(method)
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
process.s.weights <- function(s.weights, data = NULL) {
  #Process s.weights
  if (is_not_null(s.weights)) {
    if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
      stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
    }
    if (is.character(s.weights) && length(s.weights)==1) {
      if (is_null(data)) {
        stop("s.weights was specified as a string but there was no argument to data.", call. = FALSE)
      }
      else if (s.weights %in% names(data)) {
        s.weights <- data[[s.weights]]
      }
      else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
    }
  }
  return(s.weights)
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
get.treat.type <- function(treat) {
  #Returns treat with treat.type attribute
  nunique.treat <- nunique(treat)
  if (nunique.treat == 2) {
    treat.type <- "binary"
  }
  else if (nunique.treat < 2) {
    stop("The treatment must have at least two unique values.", call. = FALSE)
  }
  else if (is.factor(treat) || is.character(treat)) {
    treat.type <- "multinomial"
    treat <- factor(treat)
  }
  else {
    treat.type <- "continuous"
  }
  attr(treat, "treat.type") <- treat.type
  return(treat)
}
check.moments.int <- function(method, moments, int) {
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
check.MSM.method <- function(method, is.MSM.method) {
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
text.box.plot <- function(range.list, width = 12) {
  full.range <- range(unlist(range.list))
  ratio = diff(full.range)/(width+1)
  rescaled.range.list <- lapply(range.list, function(x) round(x/ratio))
  rescaled.full.range <- round(full.range/ratio)
  d <- as.data.frame(matrix(NA_character_, ncol = 3, nrow = length(range.list),
                            dimnames = list(names(range.list), c("Min", paste(rep(" ", width + 1), collapse = ""), "Max"))),
                     stringsAsFactors = FALSE)
  d[,"Min"] <- vapply(range.list, function(x) x[1], numeric(1L))
  d[,"Max"] <- vapply(range.list, function(x) x[2], numeric(1L))
  for (i in seq_len(nrow(d))) {
    spaces1 <- rescaled.range.list[[i]][1] - rescaled.full.range[1]
    #|
    dashes <- max(0, diff(rescaled.range.list[[i]]) - 2)
    #|
    spaces2 <- max(0, diff(rescaled.full.range) - (spaces1 + 1 + dashes + 1))

    d[i, 2] <- paste0(paste(rep(" ", spaces1), collapse = ""), "|", paste(rep("-", dashes), collapse = ""), "|", paste(rep(" ", spaces2), collapse = ""))
  }
  return(d)
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
make.full.rank <- function(mat, with.intercept = TRUE) {

  if (!all(apply(mat, 2, is.numeric))) stop("All columns in mat must be numeric.", call. = FALSE)

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

  return(mat[, keep, drop = FALSE])
}
int.poly.f <- function(mat, ex=NULL, int=FALSE, poly=1, center = FALSE, sep = " * ") {
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

#Shared with cobalt
word.list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  if (quotes) word.list <- vapply(word.list, function(x) paste0("\"", x, "\""), character(1L))
  if (L == 0) {
    out <- ""
    attr(out, "plural") = FALSE
  }
  else {
    word.list <- word.list[word.list %nin% c(NA, "")]
    L <- length(word.list)
    if (L == 0) {
      out <- ""
      attr(out, "plural") = FALSE
    }
    else if (L == 1) {
      out <- word.list
      if (is.are) out <- paste(out, "is")
      attr(out, "plural") = FALSE
    }
    else {
      and.or <- match.arg(and.or)
      if (L == 2) {
        out <- paste(word.list, collapse = paste0(" ", and.or," "))
      }
      else {
        out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "),
                     word.list[L], sep = paste0(", ", and.or," "))

      }
      if (is.are) out <- paste(out, "are")
      attr(out, "plural") = TRUE
    }


  }
  return(out)
}
binarize <- function(variable) {
  nas <- is.na(variable)
  if (!is_binary(variable[!nas])) stop(paste0("Cannot binarize ", deparse(substitute(variable)), ": more than two levels."))
  if (is.character(variable)) variable <- factor(variable)
  variable.numeric <- as.numeric(variable)
  if (!is.na(match(0, unique(variable.numeric)))) zero <- 0
  else zero <- min(unique(variable.numeric), na.rm = TRUE)
  newvar <- setNames(ifelse(!nas & variable.numeric==zero, 0, 1), names(variable))
  newvar[nas] <- NA
  return(newvar)
}
num_to_superscript <- function(x) {
  nums <- setNames(c("\u2070",
                     "\u00B9",
                     "\u00B2",
                     "\u00B3",
                     "\u2074",
                     "\u2075",
                     "\u2076",
                     "\u2077",
                     "\u2078",
                     "\u2079"),
                   as.character(0:9))
  x <- as.character(x)
  splitx <- strsplit(x, "")
  supx <- sapply(splitx, function(y) paste0(nums[y], collapse = ""))
  return(supx)
}
null.or.error <- function(x) {is_null(x) || class(x) == "try-error"}
get.covs.and.treat.from.formula <- function(f, data = NULL, env = .GlobalEnv, ...) {
  A <- list(...)

  tt <- terms(f, data = data)

  #Check if data exists
  if (is_not_null(data) && is.data.frame(data)) {
    data.specified <- TRUE
  }
  else data.specified <- FALSE

  #Check if response exists
  if (is.formula(tt, 2)) {
    resp.vars.mentioned <- as.character(tt)[2]
    resp.vars.failed <- vapply(resp.vars.mentioned, function(v) {
      null.or.error(try(eval(parse(text = v), c(data, env)), silent = TRUE))
    }, logical(1L))

    if (any(resp.vars.failed)) {
      if (is_null(A[["treat"]])) stop(paste0("The given response variable, \"", as.character(tt)[2], "\", is not a variable in ", word.list(c("data", "the global environment")[c(data.specified, TRUE)], "or"), "."), call. = FALSE)
      tt <- delete.response(tt)
    }
  }
  else resp.vars.failed <- TRUE

  if (any(!resp.vars.failed)) {
    treat.name <- resp.vars.mentioned[!resp.vars.failed][1]
    tt.treat <- terms(as.formula(paste0(treat.name, " ~ 1")))
    mf.treat <- quote(stats::model.frame(tt.treat, data,
                                         drop.unused.levels = TRUE,
                                         na.action = "na.pass"))

    tryCatch({mf.treat <- eval(mf.treat, c(data, env))},
             error = function(e) {stop(conditionMessage(e), call. = FALSE)})
    treat <- model.response(mf.treat)
  }
  else {
    treat <- A[["treat"]]
    treat.name <- NULL
  }

  #Check if RHS variables exist
  tt.covs <- delete.response(tt)
  rhs.vars.mentioned.lang <- attr(tt.covs, "variables")[-1]
  rhs.vars.mentioned <- vapply(rhs.vars.mentioned.lang, deparse, character(1L))
  rhs.vars.failed <- vapply(rhs.vars.mentioned.lang, function(v) {
    null.or.error(try(eval(v, c(data, env)), silent = TRUE))
  }, logical(1L))

  if (any(rhs.vars.failed)) {
    stop(paste0(c("All variables in formula must be variables in data or objects in the global environment.\nMissing variables: ",
                  paste(rhs.vars.mentioned[rhs.vars.failed], collapse=", "))), call. = FALSE)

  }

  rhs.term.labels <- attr(tt.covs, "term.labels")
  rhs.term.orders <- attr(tt.covs, "order")

  rhs.df <- vapply(rhs.vars.mentioned.lang, function(v) {
    is.data.frame(try(eval(v, c(data, env)), silent = TRUE))
  }, logical(1L))

  if (any(rhs.df)) {
    if (any(rhs.vars.mentioned[rhs.df] %in% unlist(sapply(rhs.term.labels[rhs.term.orders > 1], function(x) strsplit(x, ":", fixed = TRUE))))) {
      stop("Interactions with data.frames are not allowed in the input formula.", call. = FALSE)
    }
    addl.dfs <- setNames(lapply(rhs.vars.mentioned.lang[rhs.df], function(x) {eval(x, env)}),
                         rhs.vars.mentioned[rhs.df])

    for (i in rhs.term.labels[rhs.term.labels %in% rhs.vars.mentioned[rhs.df]]) {
      ind <- which(rhs.term.labels == i)
      rhs.term.labels <- append(rhs.term.labels[-ind],
                                values = names(addl.dfs[[i]]),
                                after = ind - 1)
    }
    new.form <- as.formula(paste("~", paste(rhs.term.labels, collapse = " + ")))

    tt.covs <- terms(new.form)
    if (is_not_null(data)) data <- do.call("cbind", unname(c(addl.dfs, list(data))))
    else data <- do.call("cbind", unname(addl.dfs))
  }

  #Get model.frame, report error
  mf.covs <- quote(stats::model.frame(tt.covs, data,
                                      drop.unused.levels = TRUE,
                                      na.action = "na.pass"))

  tryCatch({covs <- eval(mf.covs, c(data, env))},
           error = function(e) {stop(conditionMessage(e), call. = FALSE)})

  if (is_not_null(treat.name) && treat.name %in% names(covs)) stop("The variable on the left side of the formula appears on the right side too.", call. = FALSE)

  if (is_null(rhs.vars.mentioned)) {
    covs <- data.frame(Intercept = rep(1, if (is_null(treat)) 1 else length(treat)))
  }
  else attr(tt.covs, "intercept") <- 0

  #Get full model matrix with interactions too
  covs.matrix <- model.matrix(tt.covs, data = covs,
                              contrasts.arg = lapply(Filter(is.factor, covs),
                                                     contrasts, contrasts=FALSE))
  attr(covs, "terms") <- NULL

  return(list(reported.covs = covs,
              model.covs = covs.matrix,
              treat = treat,
              treat.name = treat.name))
}
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  nas <- is.na(df)
  if (!is.data.frame(df)) df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  rn <- rownames(df)
  cn <- colnames(df)
  df <- as.data.frame(lapply(df, function(col) {
    if (suppressWarnings(all(!is.na(as.numeric(as.character(col)))))) {
      as.numeric(as.character(col))
    } else {
      col
    }
  }), stringsAsFactors = FALSE)
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  o.negs <- sapply(1:ncol(df), function(x) if (nums[x]) df[[x]] < 0 else rep(FALSE, length(df[[x]])))
  df[nums] <- round(df[nums], digits = digits)
  df[nas] <- ""

  df <- as.data.frame(lapply(df, format, scientific = FALSE, justify = "none"), stringsAsFactors = FALSE)

  for (i in which(nums)) {
    if (any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      lengths <- lengths(s)
      digits.r.of.. <- vapply(seq_along(s), function(x) {
        if (lengths[x] > 1) nchar(s[[x]][lengths[x]])
        else 0 }, numeric(1L))
      df[[i]] <- sapply(seq_along(df[[i]]), function(x) {
        if (df[[i]][x] == "") ""
        else if (lengths[x] <= 1) {
          paste0(c(df[[i]][x], rep(".", pad == 0), rep(pad, max(digits.r.of..) - digits.r.of..[x] + as.numeric(pad != 0))),
                 collapse = "")
        }
        else paste0(c(df[[i]][x], rep(pad, max(digits.r.of..) - digits.r.of..[x])),
                    collapse = "")
      })
    }
  }

  df[o.negs & df == 0] <- paste0("-", df[o.negs & df == 0])

  # Insert NA placeholders
  df[nas] <- na_vals

  if (length(rn) > 0) rownames(df) <- rn
  if (length(cn) > 0) names(df) <- cn

  return(df)
}
nunique <- function(x, nmax = NA, na.rm = TRUE) {
  if (is_null(x)) return(0)
  else {
    if (na.rm) x <- x[!is.na(x)]
    if (is.factor(x)) return(nlevels(x))
    else return(length(unique(x, nmax = nmax)))
  }

}
nunique.gt <- function(x, n, na.rm = TRUE) {
  if (missing(n)) stop("n must be supplied.")
  if (n < 0) stop("n must be non-negative.")
  if (is_null(x)) FALSE
  else {
    if (na.rm) x <- x[!is.na(x)]
    if (n == 1 && is.numeric(x)) !check_if_zero(max(x) - min(x))
    else if (length(x) < 2000) nunique(x) > n
    else tryCatch(nunique(x, nmax = n) > n, error = function(e) TRUE)
  }
}
all_the_same <- function(x) !nunique.gt(x, 1)
is.formula <- function(f, sides = NULL) {
  res <- is.name(f[[1]])  && deparse(f[[1]]) %in% c( '~', '!') &&
    length(f) >= 2
  if (is_not_null(sides) && is.numeric(sides) && sides %in% c(1,2)) {
    res <- res && length(f) == sides + 1
  }
  return(res)
}
check_if_zero <- function(x) {
  # this is the default tolerance used in all.equal
  tolerance <- .Machine$double.eps^0.5
  # If the absolute deviation between the number and zero is less than
  # the tolerance of the floating point arithmetic, then return TRUE.
  # This means, to me, that I can treat the number as 0 rather than
  # -3.20469e-16 or some such.
  abs(x - 0) < tolerance
}
is_binary <- function(x) !nunique.gt(x, 2)
is_null <- function(x) length(x) == 0L
is_not_null <- function(x) !is_null(x)
probably.a.bug <- function() {
  fun <- paste(deparse(sys.call(-1)), collapse = "\n")
  stop(paste0("An error was produced and is likely a bug. Please let the maintainer know a bug was produced by the function\n",
              fun), call. = FALSE)
}
ESS <- function(w) {
  sum(w)^2/sum(w^2)
}
center <- function(x, na.rm = TRUE, at = NULL) {
  if (!is.numeric(x)) warning("x is not numeric and will not be centered.")
  else {
    if (is.matrix(x)) x <- apply(x, 2, center, na.rm = na.rm, at = at)
    else {
      if (is_null(at)) at <- mean(x, na.rm = na.rm)
      else if (!is.numeric(at)) stop("at must be numeric.")
      x <- x - at
    }
  }
  return(x)
}

#Other helpful tools
##Check if value is between other values
between <- function(x, range, inclusive = TRUE, na.action = FALSE) {
  if (!all(is.numeric(x))) stop("x must be a numeric vector.", call. = FALSE)
  if (length(range) != 2) stop("range must be of length 2.", call. = FALSE)
  if (any(is.na(range) | !is.numeric(range))) stop("range must contain numeric entries only.", call. = FALSE)
  range <- sort(range)

  if (any(is.na(x))) {
    if (length(na.action) != 1 || !is.atomic(na.action)) stop("na.action must be an atomic vector of length 1.", call. = FALSE)
  }
  if (inclusive) out <- ifelse(is.na(x), na.action, x >= range[1] & x <= range[2])
  else out <- ifelse(is.na(x), na.action, x > range[1] & x < range[2])

  return(out)
}
##Check if factors are equivalent
equivalent.factors <- function(f1, f2) {
  return(nunique(f1) == nunique(interaction(f1, f2)))
}
##Produce ordinal version of number (e.g., 1 -> 1st, etc.)
ordinal <- function(x) {
  x <- abs(x)
  if (as.integer(substring(x, seq(nchar(x)), seq(nchar(x))))[nchar(x)] == 1) {
    return(paste0(x, 'st'))
  } else if (as.integer(substring(x, seq(nchar(x)), seq(nchar(x))))[nchar(x)] == 2) {
    return(paste0(x, 'nd'))
  } else if (as.integer(substring(x, seq(nchar(x)), seq(nchar(x))))[nchar(x)] == 3) {
    return(paste0(x, 'rd'))
  } else {
    return(paste0(x, 'th'))
  }
}
## Making a "opposite of %in%" or "not %in%" function to simplify code
`%nin%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
## Quicker way to get last item of vector
last <- function(x) {return(x[length(x)])}
## Just so code reads more clearly when using last(x)
first <- function(x) {return(x[1])}

#To pass CRAN checks:
utils::globalVariables(c(".s.weights"))
