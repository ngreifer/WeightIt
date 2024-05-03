#This document is shared across cobalt, WeightIt, and optweight

#Strings
word_list <- function(word.list = NULL, and.or = "and", is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  word.list <- add_quotes(word.list, quotes)

  if (L == 0) {
    out <- ""
    attr(out, "plural") <- FALSE
  }
  else {
    word.list <- word.list[!word.list %in% c(NA_character_, "")]
    L <- length(word.list)
    if (L == 0) {
      out <- ""
      attr(out, "plural") <- FALSE
    }
    else if (L == 1) {
      out <- word.list
      if (is.are) out <- paste(out, "is")
      attr(out, "plural") <- FALSE
    }
    else {
      and.or <- match_arg(and.or, c("and", "or"))
      if (L == 2) {
        out <- paste(word.list, collapse = paste0(" ", and.or, " "))
      }
      else {
        out <- paste(paste(word.list[seq_len(L - 1)], collapse = ", "),
                     word.list[L], sep = paste0(", ", and.or, " "))

      }
      if (is.are) out <- paste(out, "are")
      attr(out, "plural") <- TRUE
    }

  }
  out
}
add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) return(x)

  if (isTRUE(quotes)) quotes <- 2

  if (chk::vld_string(quotes)) x <- paste0(quotes, x, quotes)
  else if (chk::vld_whole_number(quotes)) {
    if (as.integer(quotes) == 0) return(x)
    else if (as.integer(quotes) == 1) x <- paste0("\'", x, "\'")
    else if (as.integer(quotes) == 2) x <- paste0("\"", x, "\"")
    else stop("`quotes` must be boolean, 1, 2, or a string.")
  }
  else {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  x
}
firstup <- function(x) {
  #Capitalize first letter
  `substr<-`(x, 1, 1, toupper(substr(x, 1, 1)))
}
expand.grid_string <- function(..., collapse = "") {
  do.call("paste", c(expand.grid(...), sep = collapse))
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
  splitx <- strsplit(x, "", fixed = TRUE)

  vapply(splitx, function(y) paste0(nums[y], collapse = ""), character(1L))
}
ordinal <- function(x) {
  if (!is.numeric(x) || !is.vector(x) || is_null(x)) stop("'x' must be a numeric vector.")
  if (length(x) > 1) return(vapply(x, ordinal, character(1L)))

  x0 <- abs(x)
  out <- paste0(x0, switch(substring(x0, nchar(x0), nchar(x0)),
                           "1" = "st",
                           "2" = "nd",
                           "3" = "rd",
                           "th"))
  if (sign(x) == -1) out <- paste0("-", out)

  out
}
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  if (NROW(df) == 0 || NCOL(df) == 0) return(df)
  if (!is.data.frame(df)) df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  rn <- rownames(df)
  cn <- colnames(df)

  infs <- o.negs <- array(FALSE, dim = dim(df))
  nas <- is.na(df)
  nums <- vapply(df, is.numeric, logical(1))
  infs[,nums] <- vapply(which(nums), function(i) !nas[,i] & !is.finite(df[[i]]), logical(NROW(df)))

  for (i in which(!nums)) {
    if (can_str2num(df[[i]])) {
      df[[i]] <- str2num(df[[i]])
      nums[i] <- TRUE
    }
  }

  o.negs[,nums] <- !nas[,nums] & df[nums] < 0 & round(df[nums], digits) == 0
  df[nums] <- round(df[nums], digits = digits)

  for (i in which(nums)) {
    df[[i]] <- format(df[[i]], scientific = FALSE, justify = "none", trim = TRUE,
                      drop0trailing = !identical(as.character(pad), "0"))

    if (!identical(as.character(pad), "0") && any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      lengths <- lengths(s)
      digits.r.of.. <- rep(0, NROW(df))
      digits.r.of..[lengths > 1] <- nchar(vapply(s[lengths > 1], `[[`, character(1L), 2))
      max.dig <- max(digits.r.of..)

      dots <- ifelse(lengths > 1, "", if (as.character(pad) != "") "." else pad)
      pads <- vapply(max.dig - digits.r.of.., function(n) paste(rep(pad, n), collapse = ""), character(1L))

      df[[i]] <- paste0(df[[i]], dots, pads)
    }
  }

  df[o.negs] <- paste0("-", df[o.negs])

  # Insert NA placeholders
  df[nas] <- na_vals
  df[infs] <- "N/A"

  if (length(rn) > 0) rownames(df) <- rn
  if (length(cn) > 0) names(df) <- cn

  df
}
text_box_plot <- function(range.list, width = 12) {
  full.range <- range(unlist(range.list))
  if (all_the_same(full.range)) {
    for (i in seq_along(range.list)) {
      range.list[[i]][1] <- range.list[[i]][1] - 1e-6
      range.list[[i]][2] <- range.list[[i]][2] + 1e-6
    }
    full.range <- range(unlist(range.list))
  }
  ratio <- diff(full.range) / (width + 1)
  rescaled.range.list <- lapply(range.list, function(x) round(x/ratio))
  rescaled.full.range <- round(full.range/ratio)
  d <- make_df(c("Min", paste(rep(" ", width + 1), collapse = ""), "Max"),
               names(range.list),
               "character")
  d[["Min"]] <- vapply(range.list, function(x) x[1], numeric(1L))
  d[["Max"]] <- vapply(range.list, function(x) x[2], numeric(1L))
  for (i in seq_row(d)) {
    spaces1 <- rescaled.range.list[[i]][1] - rescaled.full.range[1]
    #|
    dashes <- max(c(0, diff(rescaled.range.list[[i]]) - 2))
    #|
    spaces2 <- max(c(0, diff(rescaled.full.range) - (spaces1 + 1 + dashes + 1)))

    d[i, 2] <- paste0(paste(rep(" ", spaces1), collapse = ""), "|", paste(rep("-", dashes), collapse = ""), "|", paste(rep(" ", spaces2), collapse = ""))
  }

  d
}
equivalent.factors <- function(f1, f2) {
  nu1 <- nunique(f1)
  nu2 <- nunique(f2)
  if (nu1 == nu2) {
    return(nu1 == nunique(paste.(f1, f2)))
  }
  FALSE
}
equivalent.factors2 <- function(f1, f2) {
  return(qr(cbind(1, as.numeric(f1), as.numeric(f2)))$rank == 2)
}
paste. <- function(..., collapse = NULL) {
  #Like paste0 but with sep = ".'
  paste(..., sep = ".", collapse = collapse)
}
wrap <- function(s, nchar, ...) {
  vapply(s, function(s_) {
    x <- strwrap(s_, width = nchar, ...)
    paste(x, collapse = "\n")
  }, character(1L))
}
strsplits <- function(x, splits, fixed = TRUE, ...) {
  #Link strsplit but takes multiple split values.
  #Only works for one string at a time (in x).
  for (split in splits) x <- unlist(strsplit(x, split, fixed = TRUE, ...))
  return(x[x != ""]) # Remove empty values
}
can_str2num <- function(x) {
  if (is.numeric(x) || is.logical(x)) return(TRUE)
  nas <- is.na(x)
  suppressWarnings(x_num <- as.numeric(as.character(x[!nas])))

  !anyNA(x_num)
}
str2num <- function(x) {
  nas <- is.na(x)
  if (!is.numeric(x) && !is.logical(x)) x <- as.character(x)
  suppressWarnings(x_num <- as.numeric(x))
  is.na(x_num[nas]) <- TRUE
  x_num
}
trim_string <- function(x, char = " ", symmetrical = TRUE, recursive = TRUE) {
  sw <- startsWith(x, char)
  ew <- endsWith(x, char)

  if (symmetrical) {
    if (!any(sw & ew)) return(x)

    x[sw & ew] <- gsub('^.|.$', '', x[sw & ew])
  }
  else {
    asw <- any(sw)
    aew <- any(ew)

    if (!asw && !aew) return(x)

    if (asw) x[sw] <- gsub('^.', '', x[sw])
    if (aew) x[ew] <- gsub('.$', '', x[ew])
  }

  if (!recursive) return(x)

  trim_string(x, char, symmetrical, recursive)
}

#Numbers
check_if_zero <- function(x) {
  # this is the default tolerance used in all.equal
  tolerance <- .Machine$double.eps^0.5
  abs(x) < tolerance
}
between <- function(x, range, inclusive = TRUE, na.action = FALSE) {
  if (!all(is.numeric(x))) stop("'x' must be a numeric vector.", call. = FALSE)
  if (length(range) != 2) stop("'range' must be of length 2.", call. = FALSE)
  if (anyNA(range) || !is.numeric(range)) stop("'range' must contain numeric entries only.", call. = FALSE)

  if (range[2] < range[1]) range <- c(range[2], range[1])

  if (anyNA(x)) {
    if (length(na.action) != 1 || !is.atomic(na.action)) stop("'na.action' must be an atomic vector of length 1.", call. = FALSE)
  }
  if (inclusive) out <- ifelse(is.na(x), na.action, x >= range[1] & x <= range[2])
  else out <- ifelse(is.na(x), na.action, x > range[1] & x < range[2])

  out
}
max_ <- function(..., na.rm = TRUE) {
  if (!any(is.finite(unlist(list(...))))) NA_real_
  else max(..., na.rm = na.rm)
}
min_ <- function(..., na.rm = TRUE) {
  if (!any(is.finite(unlist(list(...))))) NA_real_
  else min(..., na.rm = na.rm)
}
check_if_int <- function(x) {
  #Checks if integer-like
  if (is.integer(x)) rep(TRUE, length(x))
  else if (is.numeric(x)) check_if_zero(x - round(x))
  else rep(FALSE, length(x))
}
squish <- function(p, lo = 1e-6, hi = 1 - lo) {
  pmax(pmin(p, hi), lo)
}

#Statistics
binarize <- function(variable, zero = NULL, one = NULL) {
  var.name <- deparse1(substitute(variable))
  if (is.character(variable) || is.factor(variable)) {
    variable <- factor(variable, nmax = if (is.factor(variable)) nlevels(variable) else NA)
    unique.vals <- levels(variable)
  }
  else {
    unique.vals <- unique(variable)
  }

  if (length(unique.vals) == 1) return(setNames(rep(1, length(variable)), names(variable)))
  if (length(unique.vals) != 2) stop(paste0("Cannot binarize ", var.name, ": more than two levels."))

  if (is_null(zero)) {
    if (is_null(one)) {
      if (can_str2num(unique.vals)) {
        variable.numeric <- str2num(variable)
      }
      else {
        variable.numeric <- as.numeric(as.factor(variable))
      }

      if (0 %in% variable.numeric) zero <- 0
      else zero <- min(variable.numeric, na.rm = TRUE)

      return(setNames(as.integer(variable.numeric != zero), names(variable)))
    }
    else {
      if (one %in% unique.vals)
        return(setNames(as.integer(variable == one), names(variable)))
      stop("The argument to 'one' is not the name of a level of variable.", call. = FALSE)
    }
  }
  else {
    if (zero %in% unique.vals)
      return(setNames(as.integer(variable != zero), names(variable)))
    stop("The argument to 'zero' is not the name of a level of variable.", call. = FALSE)
  }
}
## ESS
center <- function(x, at = NULL, na.rm = TRUE) {
  if (is.data.frame(x)) {
    x <- as.matrix.data.frame(x)
    type <- "df"
  }
  if (!is.numeric(x)) stop("'x' must be numeric.")
  else if (is.array(x) && length(dim(x)) > 2) stop("'x' must be a numeric or matrix-like (not array).")
  else if (!is.matrix(x)) {
    x <- matrix(x, ncol = 1)
    type <- "vec"
  }
  else type <- "matrix"
  if (is_null(at)) at <- colMeans(x, na.rm = na.rm)
  else if (length(at) %nin% c(1, ncol(x))) stop("'at' is not the right length.")
  out <- x - matrix(at, byrow = TRUE, ncol = ncol(x), nrow = nrow(x))

  if (type == "df") out <- as.data.frame.matrix(out)
  else if (type == "vec") out <- drop(out)

  out
}
w.m <- function(x, w = NULL, na.rm = TRUE) {
  if (is_null(w)) w <- rep(1, length(x))
  if (anyNA(x)) is.na(w[is.na(x)]) <- TRUE

  sum(x * w, na.rm = na.rm) / sum(w, na.rm = na.rm)
}
col.w.m <- function(mat, w = NULL, na.rm = TRUE) {
  if (is_null(w)) w <- 1
  w.sum <- colSums(w*!is.na(mat))

  colSums(mat*w, na.rm = na.rm) / w.sum
}
col.w.v <- function(mat, w = NULL, bin.vars = NULL, na.rm = TRUE) {
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
        stop("'mat' must be a numeric matrix.")
      }

      mat <- data.matrix(mat)
    }
    else if (is.numeric(mat)) {
      mat <- matrix(mat, ncol = 1)
    }
    else stop("'mat' must be a numeric matrix.")
  }

  if (is_null(bin.vars)) bin.vars <- rep(FALSE, ncol(mat))
  else if (length(bin.vars) != ncol(mat) || anyNA(as.logical(bin.vars))) {
    stop("'bin.vars' must be a logical vector with length equal to the number of columns of 'mat'.", call. = FALSE)
  }
  bin.var.present <- any(bin.vars)
  non.bin.vars.present <- any(!bin.vars)

  var <- setNames(numeric(ncol(mat)), colnames(mat))
  if (is_null(w)) {
    if (non.bin.vars.present) {
      den <- colSums(!is.na(mat[, !bin.vars, drop = FALSE])) - 1
      var[!bin.vars] <- colSums(center(mat[, !bin.vars, drop = FALSE])^2, na.rm = na.rm)/den
    }
    if (bin.var.present) {
      means <- colMeans(mat[, bin.vars, drop = FALSE], na.rm = na.rm)
      var[bin.vars] <- means * (1 - means)
    }
  }
  else if (na.rm && anyNA(mat)) {
    # n <- nrow(mat)
    w <- array(w, dim = dim(mat))
    is.na(w[is.na(mat)]) <- TRUE
    s <- colSums(w, na.rm = na.rm)
    w <- mat_div(w, s)
    if (non.bin.vars.present) {
      x <- sqrt(w[, !bin.vars, drop = FALSE]) * center(mat[, !bin.vars, drop = FALSE],
                                                       at = colSums(w[, !bin.vars, drop = FALSE] * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
      var[!bin.vars] <- colSums(x*x, na.rm = na.rm)/(1 - colSums(w[, !bin.vars, drop = FALSE]^2, na.rm = na.rm))
    }
    if (bin.var.present) {
      means <- colSums(w[, bin.vars, drop = FALSE] * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
      var[bin.vars] <- means * (1 - means)
    }
  }
  else {
    if (is_null(w)) w <- rep(1, nrow(mat))
    w <- w/sum(w)
    if (non.bin.vars.present) {
      x <- sqrt(w) * center(mat[, !bin.vars, drop = FALSE],
                            at = colSums(w * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
      var[!bin.vars] <- colSums(x*x, na.rm = na.rm)/(1 - sum(w^2))
    }
    if (bin.var.present) {
      means <- colSums(w * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
      var[bin.vars] <- means * (1 - means)
    }
  }

  var
}
col.w.cov <- function(mat, y, w = NULL, na.rm = TRUE) {
  if (!is.matrix(mat)) {
    if (is_null(w)) return(cov(mat, y, use = if (na.rm) "pair" else "everything"))
    else mat <- matrix(mat, ncol = 1)
  }
  if (is_null(w)) {
    y <- array(y, dim = dim(mat))
    if (anyNA(mat)) is.na(y[is.na(mat)]) <- TRUE
    if (anyNA(y)) is.na(mat[is.na(y)]) <- TRUE
    den <- colSums(!is.na(mat*y)) - 1
    cov <- colSums(center(mat, na.rm = na.rm)*center(y, na.rm = na.rm), na.rm = na.rm)/den
  }
  else if (na.rm && anyNA(mat)) {
    n <- nrow(mat)
    w <- array(w, dim = dim(mat))
    is.na(w[is.na(mat)]) <- TRUE
    s <- colSums(w, na.rm = na.rm)
    w <- mat_div(w, s)
    x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
    cov <- colSums(x*y, na.rm = na.rm)/(1 - colSums(w^2, na.rm = na.rm))
  }
  else {
    n <- nrow(mat)
    w <- w/sum(w)
    x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
    cov <- colSums(x*y, na.rm = na.rm)/(1 - sum(w^2))
  }

  cov
}
col.w.r <- function(mat, y, w = NULL, s.weights = NULL, bin.vars = NULL, na.rm = TRUE) {
  if (is_null(w) && is_null(s.weights)) {
    return(cor(mat, y, use = if (na.rm) "pair" else "everything"))
  }

  cov <- col.w.cov(mat, y = y, w = w, na.rm = na.rm)
  den <- sqrt(col.w.v(mat, w = s.weights, bin.vars = bin.vars, na.rm = na.rm)) *
    sqrt(col.w.v(y, w = s.weights, na.rm = na.rm))

  cov / den
}
scale_w <- function(x, w = NULL) {
  (x - w.m(x, w)) / sqrt(col.w.v(x, w))
}
mean_abs_dev <- function(x) {
  mean_fast(abs(x - mean_fast(x, TRUE)), TRUE)
}
rms <- function(x) {
  sqrt(mean_fast(x^2))
}
mat_div <- function(mat, vec) {
  mat/vec[col(mat)]
}
abs_ <- function(x, ratio = FALSE) {
  if (ratio) pmax(x, 1/x)
  else abs(x)
}
mean_fast <- function(x, nas.possible = FALSE) {
  #Equal to mean(x, na.rm = TRUE) but faster
  #Set no.nas = FALSE if it's possible there are NAs
  if (nas.possible && anyNA(x)) {
    s <- sum(x, na.rm = TRUE)
    n <- sum(!is.na(x))
    return(s/n)
  }

  sum(x)/length(x)
}
bw.nrd <- function(x) {
  #R's bw.nrd doesn't always work, but bw.nrd0 does
  bw.nrd0(x)*1.06/.9
}
w.quantile <- function(x, probs = seq(0, 1, 0.25), w = NULL, na.rm = FALSE, ...) {

  n <- length(x)
  if (n == 0 || (!isTRUE(na.rm) && anyNA(x))) {
    return(rep(NA_real_, length(probs)))
  }

  if (!is.null(w) && all(w == 0)) {
    return(rep(0, length(probs)))
  }

  if (isTRUE(na.rm)) {
    indices <- !is.na(x)
    x <- x[indices]
    if (!is.null(w))
      w <- w[indices]
  }

  order <- order(x)
  x <- x[order]
  w <- w[order]

  rw <- {
    if (is.null(w)) (1:n)/n
    else cumsum(w)/sum(w)
  }

  q <- vapply(probs, function(p) {
    if (p == 0) return(x[1])
    if (p == 1) return(x[n])
    select <- min(which(rw > p))
    if (rw[select] == p)
      mean(x[c(select, select + 1)])
    else x[select]
  }, x[1])

  unname(q)
}

#Formulas
hasbar <- function(term) {
  any(c("|", "||") %in% all.names(term))
}
nobars <- function(term) {
  #Replace formula with version without "|"s, i.e., random effects

  isBar <- function(term) {
    is.call(term) && (term[[1]] == as.name("|") || term[[1]] == as.name("||"))
  }
  isAnyArgBar <- function(term) {
    if ((term[[1]] != as.name("~")) && (term[[1]] != as.name("("))) {
      for (i in seq_along(term)) {
        if (isBar(term[[i]])) return(TRUE)
      }
    }
    FALSE
  }

  nobars_ <- function(term) {
    if (!hasbar(term))
      return(term)
    if (isBar(term))
      return(NULL)
    if (isAnyArgBar(term))
      return(NULL)
    if (length(term) == 2) {
      nb <- nobars_(term[[2]])
      if (is.null(nb))
        return(NULL)
      term[[2]] <- nb
      return(term)
    }
    nb2 <- nobars_(term[[2]])
    nb3 <- nobars_(term[[3]])
    if (is.null(nb2))
      return(nb3)
    if (is.null(nb3))
      return(nb2)
    term[[2]] <- nb2
    term[[3]] <- nb3
    term
  }

  nb <- nobars_(term)
  if (is_(term, "formula") && length(term) == 3 && is.symbol(nb)) {
    nb <- reformulate("1", response = deparse(nb))
  }
  if (is.null(nb)) {
    nb <- if (is_(term, "formula")) {~1} else 1
  }
  nb
}

#treat/covs
get_covs_and_treat_from_formula <- function(f, data = NULL, terms = FALSE, sep = "", ...) {
  A <- list(...)

  #Check if data exists
  if (is_not_null(data)) {
    if (is.data.frame(data)) {
      data.specified <- TRUE
    }
    else {
      warning("The argument supplied to data is not a data.frame object. This may causes errors or unexpected results.", call. = FALSE)
      data <- environment(f)
      data.specified <- FALSE
    }
  }
  else {
    data <- environment(f)
    data.specified <- FALSE
  }

  env <- environment(f)

  if (!rlang::is_formula(f)) stop("'f' must be a formula.")

  eval.model.matrx <- !hasbar(f)

  tryCatch(tt <- terms(f, data = data),
           error = function(e) {
             if (conditionMessage(e) == "'.' in formula and no 'data' argument") {
               stop("'.' is not allowed in formulas.", call. = FALSE)
             }
             else stop(conditionMessage(e), call. = FALSE)
           })

  #Check if response exists
  if (rlang::is_formula(tt, lhs = TRUE)) {
    resp.var.mentioned <- attr(tt, "variables")[[2]]
    resp.var.mentioned.char <- deparse1(resp.var.mentioned)

    resp.var.failed <- {
      test <- tryCatch(eval(resp.var.mentioned, data, env), error = function(e) e)
      if (inherits(test, "simpleError")) {
        if (startsWith(conditionMessage(test), "object") &&
            endsWith(conditionMessage(test), "not found")) TRUE
        else stop(test)
      }
      else is_null(test)
    }

    if (resp.var.failed) {
      if (is_null(A[["treat"]])) stop(paste0("The given response variable, \"", resp.var.mentioned.char, "\", is not a variable in ", word_list(c("data", "the global environment")[c(data.specified, TRUE)], "or"), "."), call. = FALSE)
      tt <- delete.response(tt)
    }
  }
  else resp.var.failed <- TRUE

  if (resp.var.failed) {
    treat <- A[["treat"]]
    treat.name <- NULL
  }
  else {
    treat.name <- resp.var.mentioned.char
    treat <- eval(resp.var.mentioned, data, env)
  }

  #Check if RHS variables exist
  tt.covs <- delete.response(tt)

  rhs.vars.mentioned <- attr(tt.covs, "variables")[-1]
  rhs.vars.mentioned.char <- vapply(rhs.vars.mentioned, deparse1, character(1L))
  rhs.vars.failed <- vapply(seq_along(rhs.vars.mentioned), function(i) {
    test <- tryCatch(eval(rhs.vars.mentioned[[i]], data, env), error = function(e) e)
    if (inherits(test, "simpleError")) {
      if (startsWith(conditionMessage(test), "object") &&
          endsWith(conditionMessage(test), "not found")) return(TRUE)
      else stop(test)
    }
    else is_null(test)
  }, logical(1L))

  if (any(rhs.vars.failed)) {
    stop(paste0(c("All variables in 'formula' must be variables in 'data' or objects in the global environment.\nMissing variables: ",
                  paste(rhs.vars.mentioned.char[rhs.vars.failed], collapse= ", "))), call. = FALSE)

  }

  rhs.term.labels <- attr(tt.covs, "term.labels")
  rhs.term.orders <- attr(tt.covs, "order")

  rhs.df <- setNames(vapply(rhs.vars.mentioned, function(v) {
    length(dim(try(eval(v, data, env), silent = TRUE))) == 2L
  }, logical(1L)), rhs.vars.mentioned.char)

  rhs.term.labels.list <- setNames(as.list(rhs.term.labels), rhs.term.labels)
  if (any(rhs.df)) {
    if (any(rhs.vars.mentioned.char[rhs.df] %in% unlist(lapply(rhs.term.labels[rhs.term.orders > 1], function(x) strsplit(x, ":", fixed = TRUE))))) {
      stop("Interactions with data.frames are not allowed in the input formula.", call. = FALSE)
    }
    addl.dfs <- setNames(lapply(which(rhs.df), function(i) {
      df <- eval(rhs.vars.mentioned[[i]], data, env)
      if (is_(df, "rms")) {
        class(df) <- "matrix"
        df <- setNames(as.data.frame(as.matrix(df)), attr(df, "colnames"))
      }
      else if (can_str2num(colnames(df))) colnames(df) <- paste(rhs.vars.mentioned.char[i], colnames(df), sep = sep)
      return(as.data.frame(df))
    }),
    rhs.vars.mentioned.char[rhs.df])

    for (i in rhs.term.labels[rhs.term.labels %in% rhs.vars.mentioned.char[rhs.df]]) {
      ind <- which(rhs.term.labels == i)
      rhs.term.labels <- append(rhs.term.labels[-ind],
                                values = names(addl.dfs[[i]]),
                                after = ind - 1)
      rhs.term.labels.list[[i]] <- names(addl.dfs[[i]])
    }

    if (data.specified) data <- do.call("cbind", unname(c(addl.dfs, list(data))))
    else data <- do.call("cbind", unname(addl.dfs))
  }

  if (is_null(rhs.term.labels)) {
    new.form <- as.formula("~ 0")
    tt.covs <- terms(new.form)
    covs <- data.frame(Intercept = rep(1, if (is_null(treat)) 1 else length(treat)))[,-1, drop = FALSE]
    # if (is_not_null(treat.name) && treat.name == "Intercept") {
    #   names(covs) <- "Intercept_"
    # }
  }
  else {
    new.form.char <- paste("~", paste(vapply(names(rhs.term.labels.list), function(x) {
      if (x %in% rhs.vars.mentioned.char[rhs.df]) paste0("`", rhs.term.labels.list[[x]], "`", collapse = " + ")
      else rhs.term.labels.list[[x]]
      # try.form <- try(as.formula(paste("~", x)), silent = TRUE)
      # if (null_or_error(try.form) || (grepl("^", x, fixed = TRUE) && !startsWith(x, "I("))) {
      #     paste0("`", x, "`")
      # }
      # else x
    } , character(1L)), collapse = " + "))

    new.form <- as.formula(new.form.char)
    tt.covs <- terms(update(new.form,  ~ . - 1))

    #Get model.frame, report error
    mf.covs <- quote(stats::model.frame(tt.covs, data,
                                        drop.unused.levels = TRUE,
                                        na.action = "na.pass"))

    tryCatch({covs <- eval(mf.covs)},
             error = function(e) {stop(conditionMessage(e), call. = FALSE)})

    if (is_not_null(treat.name) && treat.name %in% names(covs)) {
      .err("the variable on the left side of the formula appears on the right side too")
    }
  }

  if (eval.model.matrx) {
    s <- !identical(sep, "")
    if (!is.character(sep) || length(sep) > 1) stop("'sep' must be a string of length 1.", call. = FALSE)
    if (s) original.covs.levels <- make_list(names(covs))
    for (i in names(covs)) {
      if (is.character(covs[[i]])) covs[[i]] <- factor(covs[[i]])
      if (is.factor(covs[[i]])) {
        if (length(unique(covs[[i]])) == 1) covs[[i]] <- 1
        else if (s) {
          original.covs.levels[[i]] <- levels(covs[[i]])
          levels(covs[[i]]) <- paste0(sep, original.covs.levels[[i]])
        }
      }
    }

    #Get full model matrix with interactions too
    covs.matrix <- model.matrix(tt.covs, data = covs,
                                contrasts.arg = lapply(Filter(is.factor, covs),
                                                       contrasts, contrasts = FALSE))

    if (s) {
      for (i in names(covs)[vapply(covs, is.factor, logical(1L))]) {
        levels(covs[[i]]) <- original.covs.levels[[i]]
      }
    }
  }
  else {
    covs.matrix <- NULL
  }

  if (!terms) attr(covs, "terms") <- NULL

  list(reported.covs = covs,
       model.covs = covs.matrix,
       treat = treat,
       treat.name = treat.name)
}
assign_treat_type <- function(treat, use.multi = FALSE) {
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

  treat
}
get_treat_type <- function(treat) {
  attr(treat, "treat.type")
}
has_treat_type <- function(treat) {
  is_not_null(get_treat_type(treat))
}
get_treated_level <- function(treat) {
  if (!is_binary(treat)) stop("'treat' must be a binary variable.")
  if (is.character(treat) || is.factor(treat)) {
    treat <- factor(treat, nmax = 2)
    unique.vals <- levels(treat)
  }
  else {
    unique.vals <- unique(treat, nmax = 2)
  }

  unique.vals.numeric <- {
    if (can_str2num(unique.vals)) str2num(unique.vals)
    else seq_along(unique.vals)
  }

  if (0 %in% unique.vals.numeric) treated <- unique.vals[unique.vals.numeric != 0]
  else treated <- unique.vals[which.max(unique.vals.numeric)]

  treated
}

#Input processing
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

  bin.vars
}
process.s.weights <- function(s.weights, data = NULL) {
  #Process s.weights
  if (is_null(s.weights)) return(NULL)

  if (is.numeric(s.weights)) return(s.weights)

  if (!is.character(s.weights) || length(s.weights) != 1) {
    .err("the argument to `s.weights` must be a vector or data frame of sampling weights or the (quoted) names of the variable in `data` that contains sampling weights")
  }

  if (is_null(data)) {
    .err("`s.weights` was specified as a string but there was no argument to `data`")
  }
  if (!s.weights %in% names(data)) {
    .err("the name supplied to `s.weights` is not the name of a variable in `data`")
  }

  data[[s.weights]]
}

#Uniqueness
nunique <- function(x, nmax = NA, na.rm = TRUE) {
  if (is_null(x)) return(0)
  if (is.factor(x)) return(nlevels(x))
  if (na.rm && anyNA(x)) x <- na.rem(x)
  length(unique(x, nmax = nmax))
}
nunique.gt <- function(x, n, na.rm = TRUE) {
  if (missing(n)) stop("'n' must be supplied.")
  if (n < 0) stop("'n' must be non-negative.")
  if (is_null(x)) FALSE
  else {
    if (n == 1) !all_the_same(x, na.rm)
    else if (length(x) < 2000) nunique(x, na.rm = na.rm) > n
    else tryCatch(nunique(x, nmax = n, na.rm = na.rm) > n, error = function(e) TRUE)
  }
}
all_the_same <- function(x, na.rm = TRUE) {
  if (anyNA(x)) {
    x <- na.rem(x)
    if (!na.rm) return(is_null(x))
  }
  if (is.double(x)) check_if_zero(max(x) - min(x))
  else all(x == x[1])
}
is_binary <- function(x, na.rm = TRUE) {
  if (na.rm && anyNA(x)) x <- na.rem(x)
  !all_the_same(x) && all_the_same(x[x != x[1]])
}
is_binary_col <- function(dat, na.rm = TRUE) {
  if (length(dim(dat)) != 2) stop("is_binary_col cannot be used with objects that don't have 2 dimensions.")
  apply(dat, 2, is_binary)
}

#R Processing
.pkgenv <- environment() #Needed for get1()
get1 <- function(name, env = .pkgenv, mode = "any", ifnotfound = NULL) {
  #Replacement for get0 that only searches within package. Provided at
  #https://stackoverflow.com/a/66897921/6348551 by user MrFlick.
  get1_ <- function(name, env, mode, ifnotfound) {
    if (identical(env, emptyenv())) {
      ifnotfound
    } else if (identical(env, globalenv())) {
      # stop at global env
      ifnotfound
    } else if (exists(name, envir = env, mode = mode, inherits = FALSE)) {
      env[[name]]
    } else {
      # try parent
      if (!identical(parent.env(env), emptyenv()) && !identical(parent.env(env), globalenv())) {
        get1_(name, parent.env(env), mode = mode, ifnotfound = ifnotfound)
      }
      else ifnotfound
    }
  }

  if (identical(env, emptyenv())) {
    ifnotfound
  } else if (identical(env, globalenv())) {
    get0(name, mode = mode, ifnotfound = ifnotfound)
  } else if (exists(name, envir = env, mode = mode, inherits = FALSE)) {
    env[[name]]
  } else {
    # try parent
    if (!identical(parent.env(env), emptyenv()) && !identical(parent.env(env), globalenv())) {
      get1_(name, parent.env(env), mode = mode, ifnotfound = ifnotfound)
    }
    else ifnotfound
  }
}
make_list <- function(n) {
  if (length(n) == 1L && is.numeric(n)) {
    vector("list", as.integer(n))
  }
  else if (length(n) > 0L && is.atomic(n)) {
    setNames(vector("list", length(n)), as.character(n))
  }
  else stop("'n' must be an integer(ish) scalar or an atomic variable.")
}
make_df <- function(ncol, nrow = 0, types = "numeric") {
  if (length(ncol) == 1L && is.numeric(ncol)) {
    col_names <- NULL
    ncol <- as.integer(ncol)
  }
  else if (is_(ncol, "atomic")) {
    col_names <- as.character(ncol)
    ncol <- length(ncol)
  }
  if (length(nrow) == 1L && is.numeric(nrow)) {
    row_names <- NULL
    nrow <- as.integer(nrow)
  }
  else if (is_(nrow, "atomic")) {
    row_names <- as.character(nrow)
    nrow <- length(nrow)
  }
  df <- as.data.frame.matrix(matrix(NA_real_, nrow = nrow, ncol = ncol))
  colnames(df) <- col_names
  rownames(df) <- row_names
  if (is_not_null(types)) {
    if (length(types) %nin% c(1, ncol)) stop("'types' must be equal to the number of columns.")
    if (any(types %nin% c("numeric", "integer", "logical", "character", NA))) {
      stop("'types' must be an acceptable type. For factors, use NA.")
    }
    if (length(types) == 1) types <- rep(types, ncol)
    for (i in seq_len(ncol)) if (!is.na(types)[i] && types[i] != "numeric") df[[i]] <- get(types[i])(nrow)
  }
  return(df)
}
ifelse_ <- function(...) {
  dotlen <- ...length()
  if (dotlen %% 2 == 0) stop("ifelse_ must have an odd number of arguments: pairs of test/yes, and one no.")
  out <- ...elt(dotlen)
  if (dotlen > 1) {
    if (!is_(out, "atomic")) stop("The last entry to ifelse_ must be atomic.")
    if (length(out) == 1) out <- rep(out, length(..1))
    n <- length(out)
    for (i in seq_len((dotlen - 1)/2)) {
      test <- ...elt(2*i - 1)
      yes <- ...elt(2*i)
      if (length(yes) == 1) yes <- rep(yes, n)
      if (length(yes) != n || length(test) != n) stop("All entries must have the same length.")
      if (!is.logical(test)) stop(paste("The", ordinal(2*i - 1), "entry to ifelse_ must be logical."))
      if (!is_(yes, "atomic")) stop(paste("The", ordinal(2*i), "entry to ifelse_ must be atomic."))
      pos <- which(test)
      out[pos] <- yes[pos]
    }
  }
  else {
    if (!is_(out, "atomic")) stop("The first entry to ifelse_ must be atomic.")
  }
  return(out)
}
is_ <- function(x, types, stop = FALSE, arg.to = FALSE) {
  s1 <- deparse1(substitute(x))
  if (is_not_null(x)) {
    for (i in types) {
      if (i == "list") it.is <- is.list(x) && !is.data.frame(x)
      else if (is_not_null(fn <- get1(paste0("is_", i), mode = "function"))) {
        it.is <- fn(x)
      }
      else if (is_not_null(fn <- get1(paste.("is", i), mode = "function"))) {
        it.is <- fn(x)
      }
      else it.is <- inherits(x, i)

      if (isTRUE(it.is)) return(TRUE)
    }
  }

  if (stop) {
    s0 <- ifelse(arg.to, "The argument to ", "")
    s2 <- ifelse(any(types %in% c("factor", "character", "numeric", "logical")),
                 "vector", "")
    stop(paste0(s0, "'", s1, "' must be a ", word_list(types, and.or = "or"), " ", s2, "."), call. = FALSE)
  }

  FALSE
}
is_null <- function(x) length(x) == 0L
is_not_null <- function(x) !is_null(x)
if_null_then <- function(x1 = NULL, x2 = NULL, ...) {
  if (is_not_null(x1)) x1
  else if (is_not_null(x2)) x2
  else if (...length() > 0) {
    for (k in seq_len(...length())) {
      if (is_not_null(...elt(k))) return(...elt(k))
    }
    return(..1)
  }
  else return(x1)
}
clear_null <- function(x) {
  x[lengths(x) == 0] <- NULL
  return(x)
}
clear_attr <- function(x, all = FALSE) {
  if (all) {
    attributes(x) <- NULL
  }
  else {
    dont_clear <- c("names", "class", "dim", "dimnames", "row.names")
    attributes(x)[names(attributes(x)) %nin% dont_clear] <- NULL
  }
  return(x)
}
probably.a.bug <- function() {
  fun <- paste(deparse1(sys.call(-1)), collapse = "\n")
  stop(paste0("An error was produced and is likely a bug. Please let the maintainer know a bug was produced by the function\n",
              fun), call. = FALSE)
}
`%nin%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
`%pin%` <- function(x, table) {
  #Partial in. TRUE if x uniquely identifies values in table.
  !is.na(pmatch(x, table))
}
`%cin%` <- function(x, table) {
  #Partial in w/ charmatch. TRUE if x at all in table.
  !is.na(charmatch(x, table))
}
null_or_error <- function(x) {is_null(x) || inherits(x, "try-error")}

#More informative and cleaner version of base::match.arg(). Uses chk.
match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg))
    stop("No argument was supplied to match_arg.")
  arg.name <- deparse1(substitute(arg), width.cutoff = 500L)

  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (length(arg) == 0) return(choices[1L])

  if (several.ok) {
    chk::chk_character(arg, add_quotes(arg.name, "`"))
  }
  else {
    chk::chk_string(arg, add_quotes(arg.name, "`"))
    if (identical(arg, choices)) return(arg[1L])
  }

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    .err(sprintf("the argument to `%s` should be %s%s.",
                 arg.name, ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                 word_list(choices, and.or = "or", quotes = 2)))
  i <- i[i > 0L]

  choices[i]
}

grab <- function(x, what) {
  lapply(x, function(z) z[[what]])
}
last <- function(x) {
  x[[length(x)]]
}
`last<-` <- function(x, value) {
  `[[<-`(x, length(x), value)
}
len <- function(x, recursive = TRUE) {
  if (is_null(x)) 0L
  else if (length(dim(x)) > 1) NROW(x)
  else if (is.list(x) && recursive) vapply(x, len, numeric(1L), recursive = FALSE)
  else length(x)
}
seq_row <- function(x) {
  if (is_null(x)) return(integer(0))
  if (length(dim(x)) != 2) stop("dim(x) must be 2")
  seq_len(NROW(x))
}
seq_col <- function(x) {
  if (is_null(x)) return(integer(0))
  if (length(dim(x)) != 2) stop("dim(x) must be 2")
  seq_len(NCOL(x))
}
na.rem <- function(x) {
  #A faster na.omit for vectors
  x[!is.na(x)]
}
anyNA_col <- function(x) {
  colSums(is.na(x)) > 0
}
check_if_call_from_fun <- function(fun) {
  # Check if called from within function f
  if (missing(fun) || !exists(deparse1(substitute(fun)), mode = "function")) return(FALSE)
  sp <- sys.parents()
  sys.funs <- lapply(sp, sys.function)
  for (x in sys.funs) {
    if (identical(fun, x)) return(TRUE)
  }
  FALSE
}

Invert <- function(f) {
  f <- match.fun(f)
  function(...) 1 / f(...)
}

.chol2 <- function(Sinv) {
  ch <- suppressWarnings(chol(Sinv, pivot = TRUE))
  p <- order(attr(ch, "pivot"))
  ch[, p, drop = FALSE]
}