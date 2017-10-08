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

process.estimand <- function(estimand, allowable.estimands, method) {
  if (!toupper(estimand) %in% toupper(allowable.estimands)) {
    stop(paste0(estimand, " is not an allowable estimand for ", method.to.phrase(method),
                ". Please select one of ", word.list(allowable.estimands, quotes = TRUE, and.or = "or"),
                "."), call. = FALSE)
  }
  else {
    return(toupper(estimand))
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
