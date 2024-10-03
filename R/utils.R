#Strings
word_list <- function(word.list = NULL, and.or = "and", is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately

  word.list <- setdiff(word.list, c(NA_character_, ""))

  if (is_null(word.list)) {
    out <- ""
    attr(out, "plural") <- FALSE
    return(out)
  }

  word.list <- add_quotes(word.list, quotes)

  L <- length(word.list)

  if (L == 1L) {
    out <- word.list
    if (is.are) out <- paste(out, "is")
    attr(out, "plural") <- FALSE
    return(out)
  }

  if (is_null(and.or) || isFALSE(and.or)) {
    out <- paste(word.list, collapse = ", ")
  }
  else {
    and.or <- match_arg(and.or, c("and", "or"))

    if (L == 2L) {
      out <- sprintf("%s %s %s",
                     word.list[1L],
                     and.or,
                     word.list[2L])
    }
    else {
      out <- sprintf("%s, %s %s",
                     paste(word.list[-L], collapse = ", "),
                     and.or,
                     word.list[L])
    }
  }

  if (is.are) out <- sprintf("%s are", out)

  attr(out, "plural") <- TRUE

  out
}
add_quotes <- function(x, quotes = 2L) {
  if (isFALSE(quotes)) {
    return(x)
  }

  if (isTRUE(quotes))
    quotes <- '"'

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, quotes))
  }

  if (!chk::vld_count(quotes) || quotes > 2) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0L) {
    return(x)
  }

  x <- {
    if (quotes == 1) sprintf("'%s'", x)
    else sprintf('"%s"', x)
  }

  x
}
expand_grid_string <- function(..., collapse = "") {
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
  if (is_null(x) || !is.numeric(x)) {
    stop("'x' must be a numeric vector.")
  }

  if (length(x) > 1L) {
    out <- setNames(vapply(x, ordinal, character(1L)), names(x))
    return(out)
  }

  x0 <- abs(x)
  out <- paste0(x0, switch(substring(x0, nchar(x0), nchar(x0)),
                           "1" = "st",
                           "2" = "nd",
                           "3" = "rd",
                           "th"))
  if (x < 0) out <- sprintf("-%s", out)

  setNames(out, names(x))
}
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  if (NROW(df) == 0L || NCOL(df) == 0L) {
    return(df)
  }

  if (!is.data.frame(df)) {
    df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  }

  rn <- rownames(df)
  cn <- colnames(df)

  infs <- o.negs <- array(FALSE, dim = dim(df))
  nas <- is.na(df)
  nums <- vapply(df, is.numeric, logical(1))
  for (i in which(nums)) {
    infs[,i] <- !nas[,i] & !is.finite(df[[i]])
  }

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
      digits.r.of.. <- rep.int(0, NROW(df))
      digits.r.of..[lengths > 1] <- nchar(vapply(s[lengths > 1], `[[`, character(1L), 2))
      max.dig <- max(digits.r.of..)

      dots <- ifelse(lengths > 1, "", if (as.character(pad) != "") "." else pad)
      pads <- vapply(max.dig - digits.r.of.., function(n) paste(rep.int(pad, n), collapse = ""), character(1L))

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
  d <- make_df(c("Min", space(width + 1), "Max"),
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

    d[i, 2] <- paste0(space(spaces1),
                      "|",
                      paste(rep.int("-", dashes), collapse = ""),
                      "|",
                      space(spaces2))
  }

  d
}
equivalent.factors <- function(f1, f2) {
  nu1 <- nunique(f1)
  nu2 <- nunique(f2)

  if (nu1 != nu2) {
    return(FALSE)
  }

  nu1 == nunique(paste.(f1, f2))

}
equivalent.factors2 <- function(f1, f2) {
  qr(cbind(1, as.numeric(f1), as.numeric(f2)))$rank == 2
}
paste. <- function(..., collapse = NULL) {
  #Like paste0 but with sep = ".'
  paste(..., sep = ".", collapse = collapse)
}
can_str2num <- function(x) {
  if (is.numeric(x) || is.logical(x)) {
    return(TRUE)
  }

  nas <- is.na(x)
  suppressWarnings(x_num <- as.numeric(as.character(x[!nas])))

  !anyNA(x_num)
}
str2num <- function(x) {
  nas <- is.na(x)
  if (!is.numeric(x) && !is.logical(x)) x <- as.character(x)
  suppressWarnings(x_num <- as.numeric(x))
  is.na(x_num)[nas] <- TRUE
  x_num
}
cat0 <- function(... , file = "", sep = "", fill = FALSE, labels = NULL,
                 append = FALSE) {
  cat(... , file = file, sep = sep, fill = fill, labels = labels,
      append = append)
}
space <- function(n) {
  paste0(rep.int(" ", n), collapse = "")
}

#Numbers
check_if_zero <- function(x) {
  # this is the default tolerance used in all.equal
  tolerance <- .Machine$double.eps^0.5
  abs(x) < tolerance
}
between <- function(x, range, inclusive = TRUE, na.action = FALSE) {
  if (!all(is.numeric(x))) {
    stop("'x' must be a numeric vector.", call. = FALSE)
  }

  if (length(range) != 2L) {
    stop("'range' must be of length 2.", call. = FALSE)
  }

  if (anyNA(range) || !is.numeric(range)) {
    stop("'range' must contain numeric entries only.", call. = FALSE)
  }

  if (range[2] < range[1]) {
    range <- c(range[2], range[1])
  }

  if (anyNA(x)) {
    if (length(na.action) != 1L || !is.logical(na.action)) {
      stop("'na.action' must be a logical vector of length 1.", call. = FALSE)
    }
  }

  out <- {
    if (inclusive) ifelse(is.na(x), na.action, x >= range[1] & x <= range[2])
    else ifelse(is.na(x), na.action, x > range[1] & x < range[2])
  }

  out
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

  if (length(unique.vals) == 1L) {
    return(setNames(rep.int(1L, length(variable)), names(variable)))
  }

  if (length(unique.vals) != 2L) {
    .err(sprintf("cannot binarize %s: more than two levels", var.name))
  }

  if (is_not_null(zero)) {
    if (zero %nin% unique.vals) {
      .err(sprintf("the argument to `zero` is not the name of a level of %s", var.name))
    }

    return(setNames(as.integer(variable != zero), names(variable)))
  }

  if (is_not_null(one)) {
    if (one %nin% unique.vals) {
      .err(sprintf("the argument to `one` is not the name of a level of %s", var.name))
    }

    return(setNames(as.integer(variable == one), names(variable)))
  }

  variable.numeric <- {
    if (can_str2num(unique.vals)) str2num(variable)
    else as.numeric(as.factor(variable))
  }

  zero <- {
    if (0 %in% variable.numeric) 0
    else min(variable.numeric, na.rm = TRUE)
  }

  setNames(as.integer(variable.numeric != zero), names(variable))
}
center <- function(x, at = NULL, na.rm = TRUE) {
  if (is.data.frame(x)) {
    x <- as.matrix.data.frame(x)
    type <- "df"
  }

  if (!is.numeric(x)) {
    stop("'x' must be numeric.")
  }

  if (is.array(x) && length(dim(x)) > 2) {
    stop("'x' must be a numeric or matrix-like (not array).")
  }

  if (is.matrix(x)) {
    type <- "matrix"
  }
  else {
    x <- matrix(x, ncol = 1)
    type <- "vec"
  }

  if (is_null(at)) {
    at <- colMeans(x, na.rm = na.rm)
  }
  else if (length(at) %nin% c(1, ncol(x))) {
    stop("'at' is not the right length.")
  }

  out <- x - matrix(at, byrow = TRUE, ncol = ncol(x), nrow = nrow(x))

  switch(type,
         "df" = as.data.frame.matrix(out),
         "vec" = drop(out),
         out)
}
w.m <- function(x, w = NULL, na.rm = TRUE) {
  if (is_null(w)) w <- rep.int(1, length(x))
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
      if (any(vapply(mat, chk::vld_character_or_factor, logical(1L)))) {
        stop("'mat' must be a numeric matrix.")
      }

      mat <- data.matrix(mat)
    }
    else if (is.numeric(mat)) {
      mat <- matrix(mat, ncol = 1)
    }
    else {
      stop("'mat' must be a numeric matrix.")
    }
  }

  if (is_null(bin.vars)) {
    bin.vars <- rep.int(FALSE, ncol(mat))
  }
  else if (length(bin.vars) != ncol(mat) || anyNA(as.logical(bin.vars))) {
    stop("'bin.vars' must be a logical vector with length equal to the number of columns of 'mat'.")
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
    is.na(w)[is.na(mat)] <- TRUE
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
    if (is_null(w)) {
      return(cov(mat, y, use = if (na.rm) "pair" else "everything"))
    }

    mat <- matrix(mat, ncol = 1)
  }

  if (is_null(w)) {
    y <- array(y, dim = dim(mat))
    if (anyNA(mat)) is.na(y)[is.na(mat)] <- TRUE
    if (anyNA(y)) is.na(mat)[is.na(y)] <- TRUE
    den <- colSums(!is.na(mat * y)) - 1
    cov <- colSums(center(mat, na.rm = na.rm) * center(y, na.rm = na.rm), na.rm = na.rm)/den
  }
  else if (na.rm && anyNA(mat)) {
    w <- array(w, dim = dim(mat))
    is.na(w)[is.na(mat)] <- TRUE
    s <- colSums(w, na.rm = na.rm)
    w <- mat_div(w, s)
    x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
    cov <- colSums(x*y, na.rm = na.rm)/(1 - colSums(w^2, na.rm = na.rm))
  }
  else {
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
mat_div <- function(mat, vec) {
  mat/vec[col(mat)]
}
mean_fast <- function(x, nas.possible = FALSE) {
  #Equal to mean(x, na.rm = TRUE) but faster
  #Set no.nas = FALSE if it's possible there are NAs
  if (!nas.possible || !anyNA(x)) {
    return(sum(x) / length(x))
  }

  sum(x, na.rm = TRUE) / sum(!is.na(x))
}
bw.nrd <- function(x) {
  #R's bw.nrd doesn't always work, but bw.nrd0 does
  bw.nrd0(x)*1.06/.9
}
w.quantile <- function(x, probs = seq(0, 1, 0.25), w = NULL, na.rm = FALSE, ...) {

  n <- length(x)
  if (n == 0L || (!isTRUE(na.rm) && anyNA(x))) {
    return(rep.int(NA_real_, length(probs)))
  }

  if (is_not_null(w) && all(check_if_zero(w))) {
    return(rep.int(0, length(probs)))
  }

  if (isTRUE(na.rm)) {
    indices <- !is.na(x)
    x <- x[indices]
    if (is_not_null(w))
      w <- w[indices]
  }

  order <- order(x)
  x <- x[order]
  w <- w[order]

  rw <- {
    if (is_null(w)) seq_len(n)/n
    else cumsum(w)/sum(w)
  }

  q <- vapply(probs, function(p) {
    if (p == 0) {
      return(x[1])
    }

    if (p == 1) {
      return(x[n])
    }

    select <- min(which(rw > p))

    if (rw[select] == p) mean(x[c(select, select + 1)])
    else x[select]
  }, x[1])

  unname(q)
}

#Formulas
hasbar <- function(term) {
  any(c("|", "||") %in% all.names(term))
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
      .wrn("the argument supplied to `data` is not a data.frame object. This may causes errors or unexpected results")
      data <- environment(f)
      data.specified <- FALSE
    }
  }
  else {
    data <- environment(f)
    data.specified <- FALSE
  }

  env <- environment(f)

  if (!rlang::is_formula(f)) {
    .err("`formula` must be a formula")
  }

  eval.model.matrx <- !hasbar(f)

  tryCatch(tt <- terms(f, data = data),
           error = function(e) {
             if (conditionMessage(e) == "'.' in formula and no 'data' argument") {
               .err("`.` is not allowed in formulas")
             }
             else {
               .err(conditionMessage(e), tidy = FALSE)
             }
           })

  #Check if response exists
  if (rlang::is_formula(tt, lhs = TRUE)) {
    resp.var.mentioned <- attr(tt, "variables")[[2]]
    resp.var.mentioned.char <- deparse1(resp.var.mentioned)

    resp.var.failed <- {
      test <- tryCatch(eval(resp.var.mentioned, data, env), error = function(e) e)
      if (!inherits(test, "simpleError")) {
        is_null(test)
      }
      else if (startsWith(conditionMessage(test), "object") &&
               endsWith(conditionMessage(test), "not found")) {
        TRUE
      }
      else {
        .err(conditionMessage(test), tidy = FALSE)
      }
    }

    if (resp.var.failed) {
      if (is_null(A[["treat"]])) {
        .err(sprintf("the given response variable, %s, is not a variable in %s",
                     add_quotes(resp.var.mentioned.char),
                     word_list(c("data", "the global environment")[c(data.specified, TRUE)], "or")))
      }
      tt <- delete.response(tt)
    }
  }
  else {
    resp.var.failed <- TRUE
  }

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
    if (!inherits(test, "simpleError")) {
      is_null(test)
    }
    else if (startsWith(conditionMessage(test), "object") &&
          endsWith(conditionMessage(test), "not found")) {
      TRUE
    }
    else {
      .err(conditionMessage(test), tidy = FALSE)
    }
  }, logical(1L))

  if (any(rhs.vars.failed)) {
    .err(sprintf("All variables in `formula` must be variables in `data` or objects in the global environment.\nMissing variables: %s",
                 word_list(rhs.vars.mentioned.char[rhs.vars.failed], and.or = FALSE)), tidy = FALSE)

  }

  rhs.term.labels <- attr(tt.covs, "term.labels")
  rhs.term.orders <- attr(tt.covs, "order")

  rhs.df <- setNames(vapply(rhs.vars.mentioned, function(v) {
    length(dim(try(eval(v, data, env), silent = TRUE))) == 2L
  }, logical(1L)), rhs.vars.mentioned.char)

  rhs.term.labels.list <- setNames(as.list(rhs.term.labels), rhs.term.labels)
  if (any(rhs.df)) {
    if (any(rhs.vars.mentioned.char[rhs.df] %in% unlist(lapply(rhs.term.labels[rhs.term.orders > 1],
                                                               function(x) strsplit(x, ":", fixed = TRUE))))) {
      .err("interactions with data.frames are not allowed in the input formula")
    }

    addl.dfs <- setNames(lapply(which(rhs.df), function(i) {
      df <- eval(rhs.vars.mentioned[[i]], data, env)
      if (inherits(df, "rms")) {
        class(df) <- "matrix"
        df <- setNames(as.data.frame(as.matrix(df)), attr(df, "colnames"))
      }
      else if (can_str2num(colnames(df))) {
        colnames(df) <- paste(rhs.vars.mentioned.char[i], colnames(df), sep = sep)
      }

      as.data.frame(df)
    }),
    rhs.vars.mentioned.char[rhs.df])

    for (i in rhs.term.labels[rhs.term.labels %in% rhs.vars.mentioned.char[rhs.df]]) {
      ind <- which(rhs.term.labels == i)
      rhs.term.labels <- append(rhs.term.labels[-ind],
                                values = names(addl.dfs[[i]]),
                                after = ind - 1)
      rhs.term.labels.list[[i]] <- names(addl.dfs[[i]])
    }

    data <- {
      if (data.specified) do.call("cbind", unname(c(addl.dfs, list(data))))
      else do.call("cbind", unname(addl.dfs))
    }
  }

  if (is_null(rhs.term.labels)) {
    new.form <- as.formula("~ 0")
    tt.covs <- terms(new.form)
    covs <- data.frame(Intercept = rep.int(1, if (is_null(treat)) 1L else length(treat)))[,-1, drop = FALSE]
    # if (is_not_null(treat.name) && treat.name == "Intercept") {
    #   names(covs) <- "Intercept_"
    # }
  }
  else {
    new.form.char <- sprintf("~ %s", paste(vapply(names(rhs.term.labels.list), function(x) {
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
             error = function(e) {.err(conditionMessage(e), tidy = FALSE)})

    if (is_not_null(treat.name) && treat.name %in% names(covs)) {
      .err("the variable on the left side of the formula appears on the right side too")
    }
  }

  if (eval.model.matrx) {
    if (!is.character(sep) || length(sep) > 1) {
      stop("'sep' must be a string of length 1.", call. = FALSE)
    }

    s <- !identical(sep, "")

    if (s) original.covs.levels <- make_list(names(covs))

    for (i in names(covs)) {
      if (is.character(covs[[i]])) {
        covs[[i]] <- factor(covs[[i]])
      }
      else if (!is.factor(covs[[i]])) {
        next
      }

      if (length(unique(covs[[i]])) == 1L) {
        covs[[i]] <- 1
      }
      else if (s) {
        original.covs.levels[[i]] <- levels(covs[[i]])
        levels(covs[[i]]) <- paste0(sep, original.covs.levels[[i]])
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
    .err("the treatment must have at least two unique values")
  }

  if (!use.multi && nunique.treat == 2) {
    treat.type <- "binary"
  }
  else if (use.multi || chk::vld_character_or_factor(treat)) {
    treat.type <- "multi-category"
    if (!inherits(treat, "processed.treat")) treat <- factor(treat)
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
  if (!is_binary(treat)) {
    stop("'treat' must be a binary variable.")
  }

  if (chk::vld_character_or_factor(treat)) {
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

  treated <- {
    if (any(unique.vals.numeric == 0)) unique.vals[unique.vals.numeric != 0]
    else unique.vals[which.max(unique.vals.numeric)]
  }

  treated
}

#Uniqueness
nunique <- function(x, nmax = NA, na.rm = TRUE) {
  if (is_null(x)) {
    return(0)
  }

  if (is.factor(x)) {
    return(nlevels(x))
  }

  if (na.rm && anyNA(x)) x <- na.rem(x)

  length(unique(x, nmax = nmax))
}
all_the_same <- function(x, na.rm = TRUE) {
  if (anyNA(x)) {
    x <- na.rem(x)
    if (!na.rm) {
      return(is_null(x))
    }
  }

  if (is.numeric(x)) check_if_zero(max(x) - min(x))
  else all(x == x[1])
}
is_binary <- function(x, na.rm = TRUE) {
  if (na.rm && anyNA(x)) x <- na.rem(x)
  !all_the_same(x) && all_the_same(x[x != x[1]])
}
is_binary_col <- function(dat, na.rm = TRUE) {
  if (length(dim(dat)) != 2) {
    stop("is_binary_col() cannot be used with objects that don't have 2 dimensions.")
  }

  apply(dat, 2, is_binary)
}

#R Processing
make_list <- function(n) {
  if (length(n) == 1L && is.numeric(n)) {
    vector("list", as.integer(n))
  }
  else if (length(n) > 0L && is.atomic(n)) {
    setNames(vector("list", length(n)), as.character(n))
  }
  else {
    stop("'n' must be an integer(ish) scalar or an atomic variable.")
  }
}
make_df <- function(ncol, nrow = 0L, types = "numeric") {
  if (is_null(ncol)) ncol <- 0L

  if (length(ncol) == 1L && is.numeric(ncol)) {
    col_names <- NULL
    ncol <- as.integer(ncol)
  }
  else if (is.atomic(ncol)) {
    col_names <- as.character(ncol)
    ncol <- length(ncol)
  }

  if (is_null(nrow)) nrow <- 0L

  if (length(nrow) == 1L && is.numeric(nrow)) {
    row_names <- NULL
    nrow <- as.integer(nrow)
  }
  else if (is.atomic(nrow)) {
    row_names <- as.character(nrow)
    nrow <- length(nrow)
  }

  df <- as.data.frame.matrix(matrix(NA_real_, nrow = nrow, ncol = ncol))

  names(df) <- col_names
  rownames(df) <- row_names

  if (is_null(types)) {
    return(df)
  }

  if (length(types) %nin% c(1L, ncol)) {
    stop("'types' must be equal to the number of columns.")
  }

  if (!is.character(types) ||
      !all(types %in% c("numeric", "integer", "logical", "character", NA))) {
    stop("'types' must be an acceptable type. For factors, use NA.")
  }

  if (length(types) == 1) types <- rep.int(types, ncol)

  for (i in which(!is.na(types))) {
    if (types[i] != "numeric") {
      df[[i]] <- get(types[i])(nrow)
    }
  }

  df
}
is_null <- function(x) length(x) == 0L
is_not_null <- function(x) !is_null(x)
if_null_then <- function(x1 = NULL, x2 = NULL, ...) {
  if (is_not_null(x1)) {
    return(x1)
  }

  dots_len <- ...length()

  if (is_not_null(x2) || dots_len == 0L) return(x2)

  if (dots_len > 1L) {
    for (k in seq_len(dots_len - 1L)) {
      if (is_not_null(...elt(k))) return(...elt(k))
    }
  }

  ...elt(dots_len)
}
clear_null <- function(x) {
  x[lengths(x) == 0] <- NULL
  x
}
clear_attr <- function(x, all = FALSE) {
  if (all) {
    attributes(x) <- NULL
  }
  else {
    dont_clear <- c("names", "class", "dim", "dimnames", "row.names")
    attributes(x)[names(attributes(x)) %nin% dont_clear] <- NULL
  }

  x
}
probably_a_bug <- function() {
  fun <- paste(deparse1(sys.call(-1)), collapse = "\n")
  .err(sprintf("an error was produced and is likely a bug. Please let the maintainer know a bug was produced by the function\n%s", fun))
}
`%nin%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
null_or_error <- function(x) {is_null(x) || inherits(x, "try-error")}

#More informative and cleaner version of base::match.arg(). Uses chk.
match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg)) {
    stop("No argument was supplied to match_arg.")
  }

  arg.name <- deparse1(substitute(arg), width.cutoff = 500L)

  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }

  if (is_null(arg)) {
    return(choices[1L])
  }

  if (several.ok) {
    chk::chk_character(arg, x_name = add_quotes(arg.name, "`"))
  }
  else {
    chk::chk_string(arg, x_name = add_quotes(arg.name, "`"))
    if (identical(arg, choices)) return(arg[1L])
  }

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    .err(sprintf("the argument to `%s` should be %s%s",
                 arg.name,
                 ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
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
  if (is_null(x)) {
    return(integer(0))
  }

  if (length(dim(x)) != 2) {
    stop("dim(x) must have length 2")
  }

  seq_len(NROW(x))
}
seq_col <- function(x) {
  if (is_null(x)) {
    return(integer(0))
  }

  if (length(dim(x)) != 2) {
    stop("dim(x) must have length 2")
  }

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
  if (missing(fun) || !exists(deparse1(substitute(fun)), mode = "function")) {
    return(FALSE)
  }

  sp <- sys.parents()
  sys.funs <- lapply(sp, sys.function)
  for (x in sys.funs) {
    if (identical(fun, x)) {
      return(TRUE)
    }
  }

  FALSE
}

# Return 1/f
Invert <- function(f) {
  f <- match.fun(f)
  function(...) 1 / f(...)
}

#Cholseky decomp with automatic pivoting
.chol2 <- function(Sinv) {
  ch <- suppressWarnings(chol(Sinv, pivot = TRUE))
  p <- order(attr(ch, "pivot"))
  ch[, p, drop = FALSE]
}
