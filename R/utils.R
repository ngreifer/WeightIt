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
    out <- toString(word.list)
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
                     toString(word.list[-L]),
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

  if (isTRUE(quotes)) {
    quotes <- '"'
  }

  if (chk::vld_string(quotes)) {
    return(paste0(quotes, x, str_rev(quotes)))
  }

  if (!chk::vld_count(quotes) || quotes > 2L) {
    stop("`quotes` must be boolean, 1, 2, or a string.")
  }

  if (quotes == 0L) {
    return(x)
  }

  x <- {
    if (quotes == 1L) sprintf("'%s'", x)
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

  as.character(x) |>
    strsplit("", fixed = TRUE) |>
    vapply(function(y) paste(nums[y], collapse = ""), character(1L))
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
  nums <- vapply(df, is.numeric, logical(1L))
  for (i in which(nums)) {
    infs[, i] <- !nas[, i] & !is.finite(df[[i]])
  }

  for (i in which(!nums)) {
    if (can_str2num(df[[i]])) {
      df[[i]] <- str2num(df[[i]])
      nums[i] <- TRUE
    }
  }

  o.negs[, nums] <- !nas[, nums] & df[nums] < 0 & round(df[nums], digits) == 0
  df[nums] <- round(df[nums], digits = digits)

  for (i in which(nums)) {
    df[[i]] <- format(df[[i]], scientific = FALSE, justify = "none", trim = TRUE,
                      drop0trailing = !identical(as.character(pad), "0"))

    if (!identical(as.character(pad), "0") && any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      s_lengths <- lengths(s)
      digits.r.of.. <- rep.int(0, NROW(df))
      digits.r.of..[s_lengths > 1L] <- nchar(vapply(s[s_lengths > 1L], `[[`, character(1L), 2L))
      max.dig <- max(digits.r.of..)

      dots <- ifelse(s_lengths > 1L, "", if (nzchar(as.character(pad))) "." else pad)
      pads <- vapply(max.dig - digits.r.of.., function(n) strrep(pad, n), character(1L))

      df[[i]] <- paste0(df[[i]], dots, pads)
    }
  }

  df[o.negs] <- paste0("-", df[o.negs])

  # Insert NA placeholders
  df[nas] <- na_vals
  df[infs] <- "N/A"

  if (is_not_null(rn)) rownames(df) <- rn
  if (is_not_null(cn)) names(df) <- cn

  df
}
text_box_plot <- function(range.list, width = 12L) {
  full.range <- range(unlist(range.list))
  if (all_the_same(full.range)) {
    for (i in seq_along(range.list)) {
      range.list[[i]][1L] <- range.list[[i]][1L] - 1e-6
      range.list[[i]][2L] <- range.list[[i]][2L] + 1e-6
    }
    full.range <- range(unlist(range.list))
  }
  ratio <- diff(full.range) / (width + 1)
  rescaled.range.list <- lapply(range.list, function(x) round(x / ratio))
  rescaled.full.range <- round(full.range / ratio)
  d <- make_df(c("Min", space(width + 1L), "Max"),
               names(range.list),
               "character")
  d[["Min"]] <- vapply(range.list, function(x) x[1L], numeric(1L))
  d[["Max"]] <- vapply(range.list, function(x) x[2L], numeric(1L))
  for (i in seq_row(d)) {
    spaces1 <- rescaled.range.list[[i]][1L] - rescaled.full.range[1L]
    dashes <- max(c(0L, diff(rescaled.range.list[[i]]) - 2L))
    spaces2 <- max(c(0L, diff(rescaled.full.range) - (spaces1 + 1L + dashes + 1L)))

    d[i, 2L] <- sprintf("%s|%s|%s",
                        space(spaces1),
                        strrep("-", dashes),
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
  x_num <- suppressWarnings(as.numeric(as.character(x[!nas])))

  !anyNA(x_num)
}
str2num <- function(x) {
  nas <- is.na(x)
  if (!is.numeric(x) && !is.logical(x)) x <- as.character(x)
  x_num <- suppressWarnings(as.numeric(x))
  is.na(x_num)[nas] <- TRUE
  x_num
}
cat0 <- function(..., file = "", sep = "", fill = FALSE, labels = NULL,
                 append = FALSE) {
  cat(..., file = file, sep = sep, fill = fill, labels = labels,
      append = append)
}
space <- function(n) {
  strrep(" ", n)
}
str_rev <- function(x) {
  strsplit(x, NULL) |>
    lapply(rev) |>
    vapply(paste, character(1L), collapse = "")
}
safe_str2expression <- function(text) {
  expr <- try(str2expression(text), silent = TRUE)

  if (null_or_error(expr)) {
    expr <- str2expression(add_quotes(text, "`"))
  }

  expr
}

#Numbers
check_if_zero <- function(x, tol = 1e-8) {
  abs(x) < tol
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

  if (range[2L] < range[1L]) {
    range <- c(range[2L], range[1L])
  }

  if (anyNA(x) && (length(na.action) != 1L || !is.logical(na.action))) {
    stop("'na.action' must be a logical vector of length 1.", call. = FALSE)
  }

  out <- {
    if (inclusive) ifelse(is.na(x), na.action, x >= range[1L] & x <= range[2L])
    else ifelse(is.na(x), na.action, x > range[1L] & x < range[2L])
  }

  out
}
squish <- function(p, lo = 1e-6, hi = 1 - lo) {
  if (lo > -Inf)
    p[p < lo] <- lo

  if (hi < Inf)
    p[p > hi] <- hi

  p
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
    return(rep_with(1L, variable))
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

  if (is.logical(variable)) {
    return(setNames(as.integer(variable), names(variable)))
  }

  if (is.numeric(variable)) {
    zero <- {
      if (any(unique.vals == 0)) 0
      else min(unique.vals, na.rm = TRUE)
    }

    return(setNames(as.integer(variable != zero), names(variable)))
  }

  variable.numeric <- {
    if (can_str2num(unique.vals)) setNames(str2num(unique.vals), unique.vals)[variable]
    else as.numeric(factor(variable, levels = unique.vals))
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

  if (is.array(x) && length(dim(x)) > 2L) {
    stop("'x' must be a numeric or matrix-like (not array).")
  }

  if (is.matrix(x)) {
    type <- "matrix"
  }
  else {
    x <- matrix(x, ncol = 1L)
    type <- "vec"
  }

  if (is_null(at)) {
    at <- colMeans(x, na.rm = na.rm)
  }
  else if (length(at) %nin% c(1L, ncol(x))) {
    stop("'at' is not the right length.")
  }

  out <- x - matrix(at, byrow = TRUE, ncol = ncol(x), nrow = nrow(x))

  switch(type,
         df = as.data.frame.matrix(out),
         vec = drop(out),
         out)
}
w.m <- function(x, w = NULL, na.rm = TRUE) {
  if (is_null(w)) w <- rep.int(1, length(x))
  if (anyNA(x)) is.na(w[is.na(x)]) <- TRUE

  sum(x * w, na.rm = na.rm) / sum(w, na.rm = na.rm)
}
col.w.m <- function(mat, w = NULL, na.rm = TRUE) {
  if (is_null(w)) w <- 1
  w.sum <- colSums(w * !is.na(mat))

  colSums(mat * w, na.rm = na.rm) / w.sum
}
col.w.v <- function(mat, w = NULL, bin.vars = NULL, na.rm = TRUE) {
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any_apply(mat, chk::vld_character_or_factor)) {
        stop("'mat' must be a numeric matrix.")
      }

      mat <- data.matrix(mat)
    }
    else if (is.numeric(mat)) {
      mat <- matrix(mat, ncol = 1L)
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
  non.bin.vars.present <- !all(bin.vars)

  var <- setNames(numeric(ncol(mat)), colnames(mat))
  if (is_null(w)) {
    if (non.bin.vars.present) {
      den <- colSums(!is.na(mat[, !bin.vars, drop = FALSE])) - 1L
      var[!bin.vars] <- colSums(center(mat[, !bin.vars, drop = FALSE])^2, na.rm = na.rm) / den
    }
    if (bin.var.present) {
      means <- colMeans(mat[, bin.vars, drop = FALSE], na.rm = na.rm)
      var[bin.vars] <- means * (1 - means)
    }
  }
  else if (na.rm && anyNA(mat)) {
    w <- array(w, dim = dim(mat))
    is.na(w)[is.na(mat)] <- TRUE
    s <- colSums(w, na.rm = na.rm)
    w <- mat_div(w, s)
    if (non.bin.vars.present) {
      x <- sqrt(w[, !bin.vars, drop = FALSE]) * center(mat[, !bin.vars, drop = FALSE],
                                                       at = colSums(w[, !bin.vars, drop = FALSE] * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
      var[!bin.vars] <- colSums(x * x, na.rm = na.rm) / (1 - colSums(w[, !bin.vars, drop = FALSE]^2, na.rm = na.rm))
    }
    if (bin.var.present) {
      means <- colSums(w[, bin.vars, drop = FALSE] * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
      var[bin.vars] <- means * (1 - means)
    }
  }
  else {
    w <- w / sum(w)
    if (non.bin.vars.present) {
      x <- sqrt(w) * center(mat[, !bin.vars, drop = FALSE],
                            at = colSums(w * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
      var[!bin.vars] <- colSums(x * x, na.rm = na.rm) / (1 - sum(w^2))
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

    mat <- matrix(mat, ncol = 1L)
  }

  if (is_null(w)) {
    y <- array(y, dim = dim(mat))
    if (anyNA(mat)) is.na(y)[is.na(mat)] <- TRUE
    if (anyNA(y)) is.na(mat)[is.na(y)] <- TRUE
    den <- colSums(!is.na(mat * y)) - 1
    cov <- colSums(center(mat, na.rm = na.rm) * center(y, na.rm = na.rm), na.rm = na.rm) / den
  }
  else if (na.rm && anyNA(mat)) {
    w <- array(w, dim = dim(mat))
    is.na(w)[is.na(mat)] <- TRUE
    s <- colSums(w, na.rm = na.rm)
    w <- mat_div(w, s)
    x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
    cov <- colSums(x * y, na.rm = na.rm) / (1 - colSums(w^2, na.rm = na.rm))
  }
  else {
    w <- w / sum(w)
    x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
    cov <- colSums(x * y, na.rm = na.rm) / (1 - sum(w^2))
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
  if (length(dim(x)) == 2L) {
    for (i in seq_col(x)) {
      x[, i] <- scale_w(x[, i], w)
    }

    return(x)
  }

  (x - w.m(x, w)) / sqrt(col.w.v(x, w))
}
center_w <- function(x, w = NULL) {
  center(x, at = col.w.m(x, w))
}
mean_abs_dev <- function(x) {
  mean_fast(abs(x - mean_fast(x, TRUE)), TRUE)
}
mat_div <- function(mat, vec) {
  mat / vec[col(mat)]
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
  bw.nrd0(x) * 1.06 / .9
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
    if (is_not_null(w)) {
      w <- w[indices]
    }
  }

  ord <- order(x)
  x <- x[ord]
  w <- w[ord]

  rw <- {
    if (is_null(w)) seq_len(n) / n
    else cumsum(w) / sum(w)
  }

  unname(vapply(probs, function(p) {
    if (p == 0) {
      return(x[1L])
    }

    if (p == 1) {
      return(x[n])
    }

    select <- min(which(rw > p))

    if (rw[select] == p) mean(x[c(select, select + 1L)])
    else x[select]
  }, x[1L]))
}

#Formulas
hasbar <- function(term) {
  any(c("|", "||") %in% all.names(term))
}
get_varnames <- function(expr) {
  recurse <- function(e) {
    if (is.symbol(e)) {
      # bare variable like age
      return(as.character(e))
    }

    if (!is.call(e)) {
      return(NULL)
    }

    # keep as-is for $, [[, and [
    fn <- e[[1L]]
    if (fn == as.name("$") || fn == as.name("[[") || fn == as.name("[")) {
      return(deparse1(e))
    }

    # strip outer function, recurse into arguments
    lapply(as.list(e)[-1L], recurse) |>
      unlist()
  }

  recurse(expr)
}

#treat/covs
get_covs_and_treat_from_formula2 <- function(f, data = NULL, sep = "", ...) {

  if (!rlang::is_formula(f)) {
    .err("`formula` must be a formula")
  }

  chk::chk_string(sep)

  env <- environment(f)

  #Check if data exists
  if (is_null(data)) {
    data <- env
    data.specified <- FALSE
  }
  else if (is.data.frame(data)) {
    data.specified <- TRUE
  }
  else {
    .wrn("the argument supplied to `data` is not a data.frame object. This may causes errors or unexpected results")
    data <- env
    data.specified <- FALSE
  }

  tryCatch({
    tt <- terms(f, data = data)
  },
  error = function(e) {
    msg <- {
      if (identical(conditionMessage(e), "'.' in formula and no 'data' argument"))
        "`.` is not allowed in formulas"
      else
        conditionMessage(e)
    }
    .err(msg)
  })

  treat <- ...get("treat")
  treat.name <- NULL
  m <- NULL

  #Check if response exists
  if (rlang::is_formula(tt, lhs = TRUE)) {
    resp.var.mentioned <- .attr(tt, "variables")[[2L]]
    resp.var.mentioned.char <- deparse1(resp.var.mentioned)

    test <- tryCatch(eval(resp.var.mentioned, data, env),
                     error = identity)

    if (inherits(test, "simpleError")) {
      m <- conditionMessage(test)
      if (!startsWith(m, "object '") || !endsWith(m, "' not found")) {
        .err(m, tidy = FALSE)
      }

      resp.var.failed <- TRUE
    }
    else {
      resp.var.failed <- is_null(test)
    }

    if (!resp.var.failed) {
      treat.name <- resp.var.mentioned.char
      treat <- test
    }
    else if (is_not_null(treat)) {
      tt <- delete.response(tt)
    }
    else {
      .err(sprintf("the given response variable, %s, is not a variable in %s",
                   add_quotes(utils::strcapture("object '(.*)' not found", m, character(1L))[[1L]]),
                   word_list(c("data", "the global environment")[c(data.specified, TRUE)], "or")))
    }
  }

  #Check if RHS variables exist
  tt.covs <- delete.response(tt)

  rhs.term.labels <- .attr(tt.covs, "term.labels")

  if (is_null(rhs.term.labels)) {
    new.form <- as.formula("~ 0")
    tt.covs <- terms(new.form)

    if (is_null(treat)) {
      covs <- make_df(ncol = 0L, nrow = 1)
      covs.matrix <- model.matrix(tt.covs, data = covs)
    }
    else {
      covs <- make_df(ncol = 0L, nrow = length(treat))
      covs.matrix <- model.matrix(tt.covs, data = covs)

      class(treat) <- unique(c("treat", class(treat)))
      attr(treat, "treat.name") <- treat.name
    }

    return(list(reported.covs = covs,
                model.covs = covs.matrix,
                simple.covs = covs,
                treat = treat))
  }

  rhs.term.orders <- .attr(tt.covs, "order")
  rhs.vars.mentioned <- .attr(tt.covs, "variables")[-1L]
  rhs.vars.mentioned.char <- vapply(rhs.vars.mentioned, deparse1, character(1L))

  simple.covs <- rhs.vars.mentioned |>
    lapply(get_varnames) |>
    unlist() |>
    unique() |>
    lapply(function(i) {
      iexp <- safe_str2expression(i)

      test <- tryCatch(eval(iexp, data, env),
                       error = identity)

      if (inherits(test, "simpleError")) {
        m <- conditionMessage(test)

        if (!startsWith(m, "object '") || !endsWith(m, "' not found")) {
          .err(m, tidy = FALSE)
        }

        return(NULL)
      }

      if (length(dim(test)) == 2L) {
        out <- as.data.frame(test)

        if (is_null(colnames(test))) {
          names(out) <- paste(i, seq_col(out), sep = sep)
        }

        return(as.list(out))
      }

      list(test) |> setNames(i)
    }) |>
    clear_null() |>
    unlist(recursive = FALSE) |>
    list2DF()

  rhs.vars.failed <- rhs.df <- rep_with(FALSE, rhs.vars.mentioned.char)
  addl.dfs <- make_list(rhs.vars.mentioned.char)
  terms.with.interactions <- unlist(lapply(rhs.term.labels[rhs.term.orders > 1], strsplit, ":", fixed = TRUE))

  for (i in seq_along(rhs.vars.mentioned.char)) {
    test <- tryCatch(eval(rhs.vars.mentioned[[i]], data, env),
                     error = identity)

    if (inherits(test, "simpleError")) {
      m <- conditionMessage(test)
      if (!startsWith(m, "object '") || !endsWith(m, "' not found")) {
        .err(m, tidy = FALSE)
      }

      rhs.vars.failed[i] <- TRUE
      rhs.vars.mentioned.char[i] <- utils::strcapture("object '(.*)' not found", m, character(1L))[[1L]]
    }

    if (any(rhs.vars.failed)) {
      next
    }

    rhs.vars.failed[i] <- is_null(test)

    if (length(dim(test)) == 2L) {
      rhs.df[i] <- TRUE

      if (rhs.vars.mentioned.char[i] %in% terms.with.interactions) {
        .err("interactions with data.frames are not allowed in the input formula")
      }

      if (inherits(test, "rms")) {
        class(test) <- "matrix"
        test <- setNames(as.data.frame(as.matrix(test)),
                         .attr(test, "colnames"))
      }
      else if (is_not_null(colnames(test))) {
        colnames(test) <- paste(rhs.vars.mentioned.char[i], colnames(test), sep = sep)
      }
      else {
        colnames(test) <- paste(rhs.vars.mentioned.char[i], seq_col(test), sep = sep)
      }

      addl.dfs[[i]] <- as.data.frame(test)
    }
  }

  if (any(rhs.vars.failed)) {
    .err(sprintf("All variables in `formula` must be variables in `data` or objects in the global environment.\nMissing variables: %s",
                 word_list(rhs.vars.mentioned.char[rhs.vars.failed], and.or = FALSE)),
         tidy = FALSE)

  }

  rhs.term.labels.list <- setNames(as.list(rhs.term.labels), rhs.term.labels)

  if (any(rhs.df)) {
    for (i in intersect(rhs.term.labels, rhs.vars.mentioned.char[rhs.df])) {
      ind <- match(i, rhs.term.labels)
      rhs.term.labels <- append(rhs.term.labels[-ind],
                                values = names(addl.dfs[[i]]),
                                after = ind - 1L)
      rhs.term.labels.list[[i]] <- names(addl.dfs[[i]])
    }

    data <- {
      if (data.specified) do.call("cbind", unname(c(clear_null(addl.dfs), list(data))))
      else do.call("cbind", unname(clear_null(addl.dfs)))
    }
  }

  new.form <- sprintf("~ %s", paste(vapply(names(rhs.term.labels.list), function(x) {
    if (x %in% rhs.vars.mentioned.char[rhs.df]) paste(add_quotes(rhs.term.labels.list[[x]], "`"), collapse = " + ")
    else rhs.term.labels.list[[x]]
  } , character(1L)), collapse = " + ")) |>
    as.formula()

  tt.covs <- terms(update(new.form,  ~ . - 1))

  #Get model.frame, report error
  mf.covs <- quote(stats::model.frame(tt.covs, data,
                                      drop.unused.levels = TRUE,
                                      na.action = "na.pass"))

  covs <- tryCatch(eval(mf.covs),
                   error = function(e) {
                     .err(conditionMessage(e), tidy = FALSE)
                   })

  if (is_not_null(treat.name) && utils::hasName(covs, treat.name)) {
    .err("the variable on the left side of the formula appears on the right side too")
  }

  s <- nzchar(sep)

  if (s) {
    original.covs.levels <- make_list(names(covs))
  }

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

  if (is_not_null(treat)) {
    class(treat) <- unique(c("treat", class(treat)))
    attr(treat, "treat.name") <- treat.name
  }

  list(reported.covs = covs,
       model.covs = covs.matrix,
       simple.covs = simple.covs,
       treat = treat)
}
assign_treat_type <- function(treat, use.multi = FALSE) {
  #Returns treat with treat.type attribute
  nunique.treat <- nunique(treat)

  if (nunique.treat < 2L) {
    .err("the treatment must have at least two unique values")
  }

  if (!use.multi && nunique.treat == 2L) {
    treat.type <- "binary"
  }
  else if (use.multi || chk::vld_character_or_factor(treat)) {
    treat.type <- "multinomial"
    if (!inherits(treat, "processed.treat")) treat <- factor(treat)
  }
  else {
    treat.type <- "continuous"
  }

  attr(treat, "treat.type") <- treat.type

  treat
}
get_treat_type <- function(treat) {
  .attr(treat, "treat.type")
}
has_treat_type <- function(treat) {
  is_not_null(get_treat_type(treat))
}

#Uniqueness
nunique <- function(x, nmax = NA, na.rm = TRUE) {
  if (is_null(x)) {
    return(0L)
  }

  if (is.factor(x)) {
    return(nlevels(x))
  }

  if (na.rm && anyNA(x)) {
    x <- na.rem(x)
  }

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
  else all(x == x[1L])
}
is_binary <- function(x, na.rm = TRUE) {
  if (na.rm && anyNA(x)) x <- na.rem(x)
  !all_the_same(x) && all_the_same(x[x != x[1L]])
}
is_binary_col <- function(dat, na.rm = TRUE) {
  if (length(dim(dat)) != 2L) {
    stop("is_binary_col() cannot be used with objects that don't have 2 dimensions.")
  }

  apply(dat, 2L, is_binary)
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

  if (length(types) == 1L) {
    types <- rep.int(types, ncol)
  }

  for (i in which(!is.na(types))) {
    if (types[i] != "numeric") {
      df[[i]] <- get(types[i])(nrow)
    }
  }

  df
}
sq_matrix <- function(x, n, names = NULL) {
  if (missing(x)) {
    x <- NA
  }

  matrix(x, nrow = n, ncol = n, dimnames = list(names, names))
}
rep_with <- function(x, y) {
  #Helper function to fill named vectors with x and given names of y
  rep.int(x, length(y)) |>
    setNames(names(y))
}
is_null <- function(x) {length(x) == 0L}
is_not_null <- function(x) {!is_null(x)}
if_null_then <- function(x1 = NULL, x2 = NULL, ...) {
  if (is_not_null(x1)) {
    return(x1)
  }

  dots_len <- ...length()

  if (is_not_null(x2) || dots_len == 0L) {
    return(x2)
  }

  if (dots_len > 1L) {
    for (k in seq_len(dots_len - 1L)) {
      if (is_not_null(...elt(k))) {
        return(...elt(k))
      }
    }
  }

  ...elt(dots_len)
}
clear_null <- function(x) {
  x[lengths(x) == 0L] <- NULL
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
  fun <- paste(deparse1(sys.call(-1L)), collapse = "\n")
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
    sysP <- sys.parent()
    formal.args <- formals(sys.function(sysP))
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
    if (identical(arg, choices)) {
      return(arg[1L])
    }
  }

  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)

  if (all(i == 0L)) {
    .err(sprintf("the argument to `%s` should be %s%s",
                 arg.name,
                 ngettext(length(choices), "", if (several.ok) "at least one of " else "one of "),
                 word_list(choices, and.or = "or", quotes = 2L)))
  }

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
  else if (length(dim(x)) > 1L) NROW(x)
  else if (is.list(x) && recursive) vapply(x, len, numeric(1L), recursive = FALSE)
  else length(x)
}
seq_row <- function(x) {
  if (is_null(x)) {
    return(integer(0L))
  }

  if (length(dim(x)) != 2L) {
    stop("dim(x) must have length 2")
  }

  seq_len(NROW(x))
}
seq_col <- function(x) {
  if (is_null(x)) {
    return(integer(0L))
  }

  if (length(dim(x)) != 2L) {
    stop("dim(x) must have length 2")
  }

  seq_len(NCOL(x))
}
na.rem <- function(x) {
  #A faster na.omit for vectors
  x[!is.na(x)]
}
anyNA_col <- function(x) {
  colSums(is.na(x)) > 0L
}
check_if_call_from_fun <- function(fun) {
  # Check if called from within function f
  if (missing(fun) || !exists(deparse1(substitute(fun)), mode = "function")) {
    return(FALSE)
  }

  for (x in sys.parents()) {
    if (identical(fun, sys.function(x), ignore.environment = TRUE)) {
      return(TRUE)
    }
  }

  FALSE
}
.attr <- function(x, which) {
  attr(x, which, exact = TRUE)
}

#Extract variables from ..., similar to ...elt(), by name without evaluating list(...)
...get <- function(x, ifnotfound = NULL) {
  expr <- quote({
    .m1 <- match(.x, ...names())
    if (anyNA(.m1)) {
      .ifnotfound
    }
    else {
      .m2 <- ...elt(.m1[1L])
      if (is_not_null(.m2)) .m2
      else .ifnotfound
    }
  })

  eval(expr,
       pairlist(.x = x[1L], .ifnotfound = ifnotfound),
       parent.frame(1L))
}

...mget <- function(x) {
  found <- match(x, eval(quote(...names()), parent.frame(1L)))

  not_found <- is.na(found)

  if (all(not_found)) {
    return(list())
  }

  setNames(lapply(found[!not_found], function(z) {
    eval(quote(...elt(.z)),
         pairlist(.z = z),
         parent.frame(3L))
  }), x[!not_found])
}

#Evaluate a call (usually a model call) with options for ignoring and recoding
#warnings and errors.
.eval_fit <- function(call, envir = parent.frame(2L), warnings = NULL, errors = NULL,
                      from = TRUE) {
  withCallingHandlers({
    fit <- eval(call, envir = envir)
  },
  warning = function(w) {
    w <- conditionMessage(w)
    wmatch <- which(startsWith(tolower(w), tolower(names(warnings))))

    if (is_not_null(wmatch)) {
      if (!is.na(warnings[wmatch[1L]])) {
        .wrn(warnings[wmatch[1L]])
      }
    }
    else if (is_null(from) || isFALSE(from)) {
      .wrn(w, tidy = FALSE)
    }
    else if (isTRUE(from)) {
      .wrn(sprintf("(from `%s()`): %s", rlang::call_name(call), w),
           tidy = FALSE)
    }
    else {
      .wrn(sprintf("(from %s): %s", paste(from, collapse = ""), w),
           tidy = FALSE)
    }

    invokeRestart("muffleWarning")
  },
  error = function(e) {
    e <- conditionMessage(e)
    ematch <- which(startsWith(tolower(e), tolower(names(errors))))

    if (is_not_null(ematch)) {
      if (!is.na(errors[ematch[1L]])) {
        .err(errors[ematch[1L]])
      }
    }
    else if (is_null(from) || isFALSE(from)) {
      .err(e, tidy = FALSE)
    }
    else if (isTRUE(from)) {
      .err(sprintf("(from `%s()`): %s", rlang::call_name(call), e),
           tidy = FALSE)
    }
    else {
      .err(sprintf("(from %s): %s", paste(from, collapse = ""), e),
           tidy = FALSE)
    }
  })

  fit
}

# Return 1/f
Invert <- function(f) {
  f <- match.fun(f)
  function(...) {1 / f(...)}
}

#Cholseky decomp with automatic pivoting
.chol2 <- function(Sinv) {
  ch <- suppressWarnings(chol(Sinv, pivot = TRUE))
  p <- order(.attr(ch, "pivot"))
  ch[, p, drop = FALSE]
}

.block_diag <- function(mat.list) {
  n <- length(mat.list)

  if (n == 0L) {
    return(sq_matrix(n = 0L))
  }

  if (n == 1L) {
    return(mat.list[[1L]])
  }

  out <- matrix(0,
                nrow = sum(vapply(mat.list, nrow, numeric(1L))),
                ncol = sum(vapply(mat.list, ncol, numeric(1L))))

  row_start <- 1L
  col_start <- 1L

  for (i in seq_len(n)) {
    mat <- mat.list[[i]]
    row_end <- row_start + nrow(mat) - 1L
    col_end <- col_start + ncol(mat) - 1L
    out[row_start:row_end, col_start:col_end] <- mat
    row_start <- row_end + 1L
    col_start <- col_end + 1L
  }

  rownames(out) <- unlist(lapply(mat.list, rownames))
  colnames(out) <- unlist(lapply(mat.list, colnames))

  out
}

.vec2list <- function(vec, n) {
  end <- cumsum(n)
  start <- c(1L, end[-length(n)] + 1L)

  lapply(seq_along(n), function(i) {
    vec[start[i]:end[i]]
  })
}

any_apply <- function(X, FUN, ...) {
  FUN <- match.fun(FUN)
  if (!is.vector(X) || is.object(X)) {
    X <- as.list(X)
  }

  for (x in X) {
    if (isTRUE(FUN(x, ...))) {
      return(TRUE)
    }
  }

  FALSE
}
all_apply <- function(X, FUN, ...) {
  FUN <- match.fun(FUN)
  if (!is.vector(X) || is.object(X)) {
    X <- as.list(X)
  }

  for (x in X) {
    if (isFALSE(FUN(x, ...))) {
      return(FALSE)
    }
  }

  TRUE
}

#crayon utilities
.it <- function(...) crayon::italic(...)
.ul <- function(...) crayon::underline(...)
.st <- function(...) crayon::strikethrough(...)