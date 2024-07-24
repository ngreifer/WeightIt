#chk utilities

pkg_caller_call <- function(start = 1) {
  package.funs <- c(getNamespaceExports(utils::packageName()),
                    .getNamespaceInfo(asNamespace(utils::packageName()), "S3methods")[, 3])
  k <- start #skip checking pkg_caller_call()
  e_max <- start
  while (!is.null(e <- rlang::caller_call(k))) {
    if (!is.null(n <- rlang::call_name(e)) &&
        n %in% package.funs) e_max <- k
    k <- k + 1
  }
  rlang::caller_call(e_max)
}

.err <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::abort(paste(strwrap(m), collapse = "\n"),
               call = pkg_caller_call(start = 2))
}
.wrn <- function(..., n = NULL, tidy = TRUE, immediate = TRUE) {
  if (immediate && isTRUE(all.equal(getOption("warn"), 0))) {
    op <- options(warn = 1)
    on.exit(options(op))
  }
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::warn(paste(strwrap(m), collapse = "\n"))
}
.msg <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::inform(paste(strwrap(m), collapse = "\n"), tidy = FALSE)
}

.chk_null_or <- function(x, chk, ..., x_name = NULL) {
  if (is.null(x_name)) {
    x_name <- deparse1(substitute(x))
  }

  x_name <- add_quotes(x_name, "`")

  if (is.null(x)) {
    return(invisible(x))
  }

  tryCatch(chk(x, ..., x_name = x_name),
           error = function(e) {
             msg <- sub("[.]$", " or `NULL`.",
                        conditionMessage(e))
             .err(msg)
           })
}

.chk_basic_vector <- function(x, x_name = NULL) {
  if (is.atomic(x) && is.null(dim(x))) {
    return(invisible(x))
  }
  if (is.null(x_name))
    x_name <- chk::deparse_backtick_chk((substitute(x)))
  chk::abort_chk(x_name, " must be an atomic, non-matrix vector", x = x)
}