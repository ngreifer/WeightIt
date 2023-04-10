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

.err <- function(...) {
  chk::err(..., call = pkg_caller_call(start = 2))
}
.wrn <- function(...) {
  chk::wrn(...)
}
.msg <- function(...) {
  chk::msg(...)
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
             .err(msg, .subclass = "chk_error")
           })
}