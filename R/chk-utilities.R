#chk utilities

#Note: this version seems to do better when used inside tryCatch()
pkg_caller_call <- function() {
  pn <- utils::packageName()
  package.funs <- c(getNamespaceExports(pn),
                    .getNamespaceInfo(asNamespace(pn), "S3methods")[, 3])

  for (i in seq_len(sys.nframe())) {
    e <- sys.call(i)

    if (is_null(n <- rlang::call_name(e))) {
      next
    }

    if (n %in% package.funs) {
      return(e)
    }
  }

  NULL
}

.err <- function(..., n = NULL, tidy = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)
  rlang::abort(paste(strwrap(m), collapse = "\n"),
               call = pkg_caller_call())
}
.wrn <- function(..., n = NULL, tidy = TRUE, immediate = TRUE) {
  if (immediate && isTRUE(all.equal(0, getOption("warn")))) {
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

.chk_basic_vector <- function(x, x_name = NULL) {
  if (is.atomic(x) && is.null(dim(x))) {
    return(invisible(x))
  }
  if (is.null(x_name))
    x_name <- chk::deparse_backtick_chk((substitute(x)))
  chk::abort_chk(x_name, " must be an atomic, non-matrix vector", x = x)
}