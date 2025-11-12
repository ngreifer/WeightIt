#chk utilities

#Note: this version seems to do better when used inside tryCatch()
pkg_caller_call <- function() {
  pn <- utils::packageName()
  package.funs <- c(getNamespaceExports(pn),
                    .getNamespaceInfo(asNamespace(pn), "S3methods")[, 3L])

  for (i in seq_len(sys.nframe())) {
    e <- sys.call(i)

    n <- rlang::call_name(e)

    if (is_null(n)) {
      next
    }

    if (n %in% package.funs) {
      return(e)
    }
  }

  NULL
}

.err <- function(..., n = NULL, tidy = TRUE) {
  chk::message_chk(..., n = n, tidy = tidy) |>
    strwrap() |>
    paste(collapse = "\n") |>
    rlang::abort(call = pkg_caller_call())
}
.wrn <- function(..., n = NULL, tidy = TRUE, immediate = TRUE) {
  m <- chk::message_chk(..., n = n, tidy = tidy)

  if (immediate && isTRUE(all.equal(0, getOption("warn")))) {
    rlang::with_options({
      m |> strwrap() |>
        paste(collapse = "\n") |>
        rlang::warn()
    }, warn = 1)
  }
  else {
    m |> strwrap() |>
      paste(collapse = "\n") |>
      rlang::warn()
  }
}
.msg <- function(..., n = NULL, tidy = TRUE) {
  chk::message_chk(..., n = n, tidy = tidy) |>
    strwrap() |>
    paste(collapse = "\n") |>
    rlang::inform(tidy = FALSE)
}

.chk_basic_vector <- function(x, x_name = NULL) {
  if (is.atomic(x) && is.null(dim(x))) {
    return(invisible(x))
  }
  if (is.null(x_name))
    x_name <- chk::deparse_backtick_chk((substitute(x)))
  chk::abort_chk(x_name, " must be an atomic, non-matrix vector", x = x)
}
