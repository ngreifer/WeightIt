word.list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  if (quotes) word.list <- sapply(word.list, function(x) paste0("\"", x, "\""))
  if (L == 0) {
    out <- ""
    attr(out, "plural") = FALSE
  }
  else {
    word.list <- word.list[!word.list %in% c(NA, "")]
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
method.to.proper.method <- function(method) {
  if (method %in% c("ps")) return("ps")
  else if (method %in% c("gbm", "gbr", "twang")) return("gbm")
  else if (method %in% c("cbps")) return("cbps")
  else if (method %in% c("npcbps")) return("npcbps")
  else if (method %in% c("entropy", "ebal", "ebalance")) return("ebal")
  else if (method %in% c("sbw")) return("sbw")
  else if (method %in% c("ebcw", "ate")) return("ebcw")
  else return(method)
}
method.to.phrase <- function(method) {
  if (method %in% c("ps")) return("propensity score weighting")
  else if (method %in% c("gbm", "gbr", "twang")) return("generalized boosted modeling")
  else if (method %in% c("cbps")) return("covariate balancing propensity score weighting")
  else if (method %in% c("npcbps")) return("non-parametric covariate balancing propensity score weighting")
  else if (method %in% c("entropy", "ebal", "ebalance")) return("entropy balancing")
  else if (method %in% c("sbw")) return("stable balancing weights")
  else if (method %in% c("ebcw", "ate")) return("empirical balancing calibration weighting")
  else return("the chosen method of weighting")
}
process.estimand <- function(estimand, method, treat.type) {
  #Allowable estimands
  AE <- list(binary = list(ps = c("ATT", "ATC", "ATE", "ATO", "ATM"),
                           gbm = c("ATT", "ATC", "ATE"),
                           cbps = c("ATT", "ATC", "ATE"),
                           npcbps = c("ATE"),
                           ebal = c("ATT", "ATC", "ATE"),
                           sbw = c("ATT", "ATC", "ATE"),
                           ebcw = c("ATT", "ATC", "ATE")),
             multinomial = list(ps = c("ATT", "ATE", "ATO"),
                                gbm = c("ATT", "ATE"),
                                cbps = c("ATT", "ATE"),
                                npcbps = c("ATE"),
                                ebal = c("ATT", "ATE"),
                                sbw = c("ATT", "ATE"),
                                ebcw = c("ATT", "ATE")))

  if (treat.type != "continuous" && !toupper(estimand) %in% AE[[treat.type]][[method]]) {
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
  if (estimand == "ATT") {
    if (length(focal) == 0L) {
      if (treat.type == "multinomial") {
        stop("When estimand = \"ATT\" for multinomial treatments, an argument must be supplied to focal.", call. = FALSE)
      }
    }
    else if (length(focal) > 1L || !any(unique(treat) == focal)) {
      stop("The argument supplied to focal must be the name of a level of treat.", call. = FALSE)
    }
  }
  else {
    if (length(focal) > 0L) {
      warning(paste(estimand, "is not compatible with focal. Ignoring focal."), call. = FALSE)
      focal <- NULL
    }
  }

  #Get focal, estimand, and reported estimand
  if (isTRUE(treat.type == "binary")) {
    unique.treat <- unique(treat, nmax = 2)
    unique.treat.bin <- unique(binarize(treat), nmax = 2)
    if (estimand == "ATT") {
      if (length(focal) == 0) {
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
check.moments.int <- function(method, moments, int) {
  if (method %in% c("ebal", "ebcw", "sbw")) {
    if (length(int) != 1 || !is.logical(int)) {
      stop("int must be a logical (TRUE/FALSE) of length 1.", call. = FALSE)
    }
    if (length(moments) != 1 || !is.numeric(moments) ||
        abs(moments - round(moments)) > sqrt(.Machine$double.eps) ||
        moments < 1) {
      stop("moments must be a positive integer of length 1.", call. = FALSE)
    }
  }
  else if (any(mi0 <- c(as.integer(moments) != 1L, int == TRUE))) {
    warning(paste0(word.list(c("moments", "int")[mi0], and.or = "and", is.are = TRUE),
                   " not compatible with ", method.to.phrase(method), ". Ignoring ", word.list(c("moments", "int")[mi0], and.or = "and"), "."), call. = FALSE)
    moments <- 1
    int <- FALSE
  }
  return(c(moments = as.integer(moments), int = int))
}
int.poly.f <- function(mat, ex=NULL, int=FALSE, poly=1, nunder=1, ncarrot=1) {
  #Adds to data frame interactions and polynomial terms; interaction terms will be named "v1_v2" and polynomials will be named "v1_2"
  #Only to be used in base.bal.tab; for general use see int.poly()
  #mat=matrix input
  #ex=matrix of variables to exclude in interactions and polynomials; a subset of df
  #int=whether to include interactions or not; currently only 2-way are supported
  #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included
  #nunder=number of underscores between variables

  if (length(ex) > 0) d <- mat[, !colnames(mat) %in% colnames(ex), drop = FALSE]
  else d <- mat
  nd <- ncol(d)
  nrd <- nrow(d)
  no.poly <- apply(d, 2, function(x) !nunique.gt(x, 2))
  npol <- nd - sum(no.poly)
  new <- matrix(ncol = (poly-1)*npol + int*(.5*(nd)*(nd-1)), nrow = nrd)
  nc <- ncol(new)
  new.names <- character(nc)
  if (poly > 1 && npol != 0) {
    for (i in 2:poly) {
      new[, (1 + npol*(i - 2)):(npol*(i - 1))] <- apply(d[, !no.poly, drop = FALSE], 2, function(x) x^i)
      new.names[(1 + npol*(i - 2)):(npol*(i - 1))] <- paste0(colnames(d)[!no.poly], paste0(replicate(ncarrot, "_"), collapse = ""), i)
    }
  }
  if (int && nd > 1) {
    new[,(nc - .5*nd*(nd-1) + 1):nc] <- matrix(t(apply(d, 1, combn, 2, prod)), nrow = nrd)
    new.names[(nc - .5*nd*(nd-1) + 1):nc] <- combn(colnames(d), 2, paste, collapse=paste0(replicate(nunder, "_"), collapse = ""))
  }

  single.value <- apply(new, 2, function(x) abs(max(x) - min(x)) < sqrt(.Machine$double.eps))
  colnames(new) <- new.names
  #new <- setNames(data.frame(new), new.names)[!single.value]
  return(new[, !single.value, drop = FALSE])
}
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
null.or.error <- function(x) {length(x) == 0 || class(x) == "try-error"}
get.covs.and.treat.from.formula <- function(f, data, env = .GlobalEnv, ...) {
  A <- list(...)

  tt <- terms(f, data = data)
  attr(tt, "intercept") <- 0

  #Check if data exists
  if (length(data) > 0 && is.data.frame(data)) {
    data.specified <- TRUE
  }
  else data.specified <- FALSE

  #Check if response exists
  if (is.formula(tt, 2)) {
    resp.vars.mentioned <- as.character(tt)[2]
    resp.vars.failed <- sapply(resp.vars.mentioned, function(v) {
      null.or.error(try(eval(parse(text = v), c(data, env)), silent = TRUE))
    })

    if (any(resp.vars.failed)) {
      if (length(A[["treat"]]) == 0) stop(paste0("The given response variable, \"", as.character(tt)[2], "\", is not a variable in ", word.list(c("data", "the global environment")[c(data.specified, TRUE)], "or"), "."), call. = FALSE)
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
  rhs.vars.mentioned <- attr(tt.covs, "term.labels")
  rhs.vars.failed <- sapply(rhs.vars.mentioned, function(v) {
    null.or.error(try(eval(parse(text = v), c(data, env)), silent = TRUE))
  })

  if (any(rhs.vars.failed)) {
    stop(paste0(c("All variables in formula must be variables in data or objects in the global environment.\nMissing variables: ",
                  paste(rhs.vars.mentioned[rhs.vars.failed], collapse=", "))), call. = FALSE)

  }

  rhs.df <- sapply(rhs.vars.mentioned, function(v) {
    is.data.frame(try(eval(parse(text = v), c(data, env)), silent = TRUE))
  })

  if (any(rhs.df)) {
    addl.dfs <- lapply(rhs.vars.mentioned[rhs.df], function(x) {eval(parse(text = x), env)})
    new.form <- paste0("~ . - ",
                       paste(rhs.vars.mentioned[rhs.df], collapse = " - "), " + ",
                       paste(unlist(sapply(addl.dfs, names)), collapse = " + "))
    tt.covs <- update(tt.covs, as.formula(new.form))
    data <- do.call("cbind", c(addl.dfs, data))
  }

  #Get model.frame, report error
  mf.covs <- quote(stats::model.frame(tt.covs, data,
                                 drop.unused.levels = TRUE,
                                 na.action = "na.pass"))
  tryCatch({covs <- eval(mf.covs, c(data, env))},
           error = function(e) {stop(conditionMessage(e), call. = FALSE)})

  if (length(rhs.vars.mentioned) == 0) covs <- data.frame(Intercept = rep(1, if (length(treat) == 0 ) 1 else length(treat)))

  return(list(covs = covs,
              treat = treat,
              treat.name = treat.name))
}
equivalent.factors <- function(f1, f2) {
  return(nunique(f1) == nunique(interaction(f1, f2)))
}
text.box.plot <- function(range.list, width = 12) {
  full.range <- range(unlist(range.list))
  ratio = diff(full.range)/(width+1)
  rescaled.range.list <- lapply(range.list, function(x) round(x/ratio))
  rescaled.full.range <- round(full.range/ratio)
  d <- as.data.frame(matrix(NA_character_, ncol = 3, nrow = length(range.list),
                            dimnames = list(names(range.list), c("Min", paste(rep(" ", width + 1), collapse = ""), "Max"))),
                     stringsAsFactors = FALSE)
  d[,"Min"] <- sapply(range.list, function(x) x[1])
  d[,"Max"] <- sapply(range.list, function(x) x[2])
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
round_df <- function(df, digits) {
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[, nums] <- round(df[, nums], digits = digits)
  return(df)
}
round_df_char <- function(df, digits, pad = "0") {
  nas <- is.na(df)
  if (!is.data.frame(df)) df <- as.data.frame(df, stringsAsFactors = FALSE)
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
  df[nums] <- round(df[nums], digits = digits)
  df[nas] <- ""

  df <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)

  for (i in which(nums)) {
    if (any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      lengths <- lengths(s)
      digits.r.of.. <- sapply(seq_along(s), function(x) {
        if (lengths[x] > 1) nchar(s[[x]][lengths[x]])
        else 0 })
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

  if (length(rn) > 0) rownames(df) <- rn
  if (length(cn) > 0) names(df) <- cn

  return(df)
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
  if (nunique.gt(x, 2)) {
    ndigits <- round(mean(floor(log10(abs(x[abs(x) > sqrt(.Machine$double.eps)]))), na.rm = TRUE))
    if (abs(ndigits) > 2) return(x/(10^ndigits))
    else return(x)
  }
  else {
    return(as.numeric(x == x[1]))
  }
}
remove.collinearity <- function(mat) {
  keep <- rep(TRUE, ncol(mat))
  for (i in seq_along(keep)) {
    keep. <- keep; keep.[i] <- FALSE
    if (qr(mat[, keep., drop = FALSE])$rank == qr(mat[, keep, drop = FALSE])$rank) {
      keep[i] <- FALSE
    }
  }
  return(mat[,keep, drop = FALSE])
}
is.formula <- function(f, sides = NULL) {
  res <- is.name(f[[1]])  && deparse(f[[1]]) %in% c( '~', '!') &&
    length(f) >= 2
  if (length(sides) > 0 && is.numeric(sides) && sides %in% c(1,2)) {
    res <- res && length(f) == sides + 1
  }
  return(res)
}
nunique <- function(x, nmax = NA, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (is.factor(x)) return(nlevels(x))
  else return(length(unique(x, nmax = nmax)))
}
nunique.gt <- function(x, n, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (length(x) < 2000) nunique(x, na.rm = na.rm) > n
  else tryCatch(nunique(x, nmax = n, na.rm = na.rm) > n, error = function(e) TRUE)
}
binarize <- function(variable) {
  if (nunique.gt(variable, 2)) stop(paste0("Cannot binarize ", deparse(substitute(variable)), ": more than two levels."))
  variable.numeric <- as.numeric(variable)
  if (!is.na(match(0, unique(variable.numeric)))) zero <- 0
  else zero <- min(unique(variable.numeric))
  newvar <- setNames(ifelse(variable.numeric==zero, 0, 1), names(variable))
  return(newvar)
}
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

#To pass CRAN checks:
utils::globalVariables(c(".s.weights"))
