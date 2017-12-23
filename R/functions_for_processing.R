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
nunique <- function(x) {
  if (is.factor(x)) return(nlevels(x))
  else return(length(unique(x)))
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

process.estimand <- function(estimand, allowable.estimands, method, treat.type) {
  if (!toupper(estimand) %in% toupper(allowable.estimands)) {
    stop(paste0("\"", estimand, "\" is not an allowable estimand for ", method.to.phrase(method),
                " with ", treat.type, " treatments. Only ", word.list(allowable.estimands, quotes = TRUE, and.or = "and", is.are = TRUE),
                " allowed."), call. = FALSE)
  }
  else {
    return(toupper(estimand))
  }

}

process.focal <- function(focal, estimand, treat) {
  if (estimand == "ATT") {
    if (length(focal) == 0) {
      stop("When estimand = \"ATT\" for multinomial treatments, an argument must be supplied to focal.", call. = FALSE)
    }
    if (length(focal) > 1 || !any(unique(treat) == focal)) {
      stop("The argument supplied to focal must be the name of a level of treat.", call. = FALSE)
    }
  }
  else {
    if (length(focal) > 0) {
      warning(paste(estimand, "is not compatible with focal. Ignoring focal."), call. = FALSE)
    }
  }
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
  package.is.intalled <- any(.packages(all.available = TRUE) == package.name)
  if (!package.is.intalled && !alternative) {
    stop(paste0("Package \"", package.name, "\" needed for this function to work. Please install it."),
         call. = FALSE)
  }
  return(invisible(package.is.intalled))
}

make.closer.to.1 <- function(x) {
  ndigits <- round(mean(floor(log10(abs(x[abs(x) > sqrt(.Machine$double.eps)])))))
  return(x/(10^ndigits))
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

#To pass CRAN checks:
utils::globalVariables(c(".s.weights"))
