#Treatment class
#' @exportS3Method `[` treat
`[.treat` <- function(x, ..., value) {
  y <- NextMethod("[")
  attr(y, "treat.type") <- .attr(x, "treat.type")
  attr(y, "treat.name") <- .attr(x, "treat.name")
  attr(y, "treated") <- .attr(x, "treated")
  attr(y, "control") <- .attr(x, "control")

  class(y) <- class(x)

  y
}

as.treat <- function(x, process = NULL) {
  if (is_null(process)) {
    process <- !inherits(x, "treat")
  }

  chk::chk_flag(process)

  if (process || !has_treat_type(x)) {
    x <- assign_treat_type(x)
    treat.type <- get_treat_type(x)

    if (treat.type %in% c("multinomial", "multi-category")) {
      x <- assign_treat_type(factor(x))
    }
  }

  class(x) <- unique(c("treat", class(x)))

  x
}
