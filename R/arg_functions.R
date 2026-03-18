#chk function replacements; will probably go in another package someday
.pkg_caller_call <- function() {
  pn <- utils::packageName()
  package.funs <- c(getNamespaceExports(pn),
                    .getNamespaceInfo(asNamespace(pn), "S3methods")[, 3L])

  for (i in seq_len(sys.nframe())) {
    e <- sys.call(i)

    n <- rlang::call_name(e)

    if (is_not_null(n) && n %in% package.funs) {
      return(e)
    }
  }

  NULL
}

## Use cli
.err <- function(m, call = .pkg_caller_call(), .envir = rlang::caller_env()) {
  .tidy_msg(m) |>
    cli::cli_abort(call = call, .envir = .envir)
}
.wrn <- function(m, immediate = TRUE, .envir = rlang::caller_env()) {
  if (isTRUE(immediate) && isTRUE(all.equal(0, getOption("warn")))) {
    rlang::local_options(warn = 1)
  }

  .tidy_msg(m) |>
    cli::cli_warn(.envir = .envir)
}
.msg <- function(m, .envir = rlang::caller_env()) {
  .tidy_msg(m) |>
    cli::cli_inform(.envir = .envir, tidy = FALSE)
}

.tidy_msg <- function(m) {
  if (length(m) != 1L) {
    return(m)
  }

  # Capitalize first letter
  if (grepl("^[[:alpha:]]", m)) {
    substr(m, 1, 1) <- toupper(substr(m, 1, 1))
  }

  # Add period to end
  if (!grepl("([.]|[?]|[!])$", m)) {
    m <- paste0(m, ".")
  }

  m
}

arg_supplied <- function(x, arg = rlang::caller_arg(x)) {
  if (missing(x)) {
    arg_expr <- substitute(x)

    if (!rlang::is_symbol(arg_expr)) {
      .err("{.arg x} must be an argument name")
    }

    .err("an argument to {.arg {arg}} must be supplied")
  }
}

arg_no_NA <- function(x, arg = rlang::caller_arg(x)) {
  if (anyNA(x)) {
    if (length(x) == 1L) {
      .err("{.arg {arg}} cannot be {.val {NA}}")
    }

    .err("{.arg {arg}} cannot contain {.val {NA}} values")
  }
}

arg_atomic <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.atomic(x)) {
    .err("{.arg {arg}} must be an atomic vector")
  }
}

arg_numeric <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.numeric(x)) {
    .err("{.arg {arg}} must be numeric")
  }

  arg_no_NA(x, arg = arg)
}

arg_whole_numeric <- function(x, arg = rlang::caller_arg(x)) {
  arg_numeric(x, arg = arg)

  if (!is.integer(x) && !check_if_zero(x - trunc(x))) {
    .err("{.arg {arg}} must be a whole numeric vector")
  }
}

arg_number <- function(x, arg = rlang::caller_arg(x)) {
  arg_numeric(x, arg = arg)

  if (length(x) != 1L) {
    .err("{.arg {arg}} must have length 1")
  }
}

arg_whole_number <- function(x, arg = rlang::caller_arg(x)) {
  arg_number(x, arg = arg)

  if (!is.integer(x) && !check_if_zero(x - trunc(x))) {
    .err("{.arg {arg}} must be a whole number")
  }
}

arg_count <- function(x, arg = rlang::caller_arg(x)) {
  arg_whole_number(x, arg = arg)
  arg_gte(x, 0, arg = arg)
}

arg_counts <- function(x, arg = rlang::caller_arg(x)) {
  arg_whole_numeric(x, arg = arg)
  arg_gte(x, 0, arg = arg)
}

arg_range <- function(x, range = c(0, 1), inclusive = TRUE, arg = rlang::caller_arg(x)) {
  if (length(inclusive) == 1L) {
    inclusive <- inclusive[c(1L, 1L)]
  }

  range <- sort(range)

  .gt_comp <- if (inclusive[1L]) function(a, b) {a >= b} else function(a, b) {a > b}
  .lt_comp <- if (inclusive[2L]) function(a, b) {a <= b} else function(a, b) {a < b}

  if (!all(.gt_comp(x, range[1L])) || !all(.lt_comp(x, range[2L]))) {

    .gt_str <- if (inclusive[1L]) "greater than or equal to" else "greater than"
    .lt_str <- if (inclusive[2L]) "less than or equal to" else "less than"

    if (length(x) == 1L) {
      .err("{.arg {arg}} must be {(.gt_str)} {range[1L]} and {(.lt_str)} {range[2L]}")
    }

    .err("all entries in {.arg {arg}} must be {(.gt_str)} {range[1L]} and {(.lt_str)} {range[2L]}")
  }
}

arg_gt <- function(x, bound = 0, arg = rlang::caller_arg(x)) {
  if (any(x <= bound)) {
    if (bound == 0) {
      .err("{.arg {arg}} must be positive")
    }

    .err("{.arg {arg}} must be greater than {.val {bound}}")
  }
}

arg_gte <- function(x, bound = 0, arg = rlang::caller_arg(x)) {
  if (any(x < bound)) {
    if (bound == 0) {
      .err("{.arg {arg}} must be non-negative")
    }

    .err("{.arg {arg}} must be greater than or equal to {.val {bound}}")
  }
}

arg_lt <- function(x, bound = 0, arg = rlang::caller_arg(x)) {
  if (any(x >= bound)) {
    if (bound == 0) {
      .err("{.arg {arg}} must be negative")
    }

    .err("{.arg {arg}} must be less than {.val {bound}}")
  }
}

arg_lte <- function(x, bound = 0, arg = rlang::caller_arg(x)) {
  if (any(x > bound)) {
    if (bound == 0) {
      .err("{.arg {arg}} must be non-positive")
    }

    .err("{.arg {arg}} must be less than or equal to {.val {bound}}")
  }
}

arg_character <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.character(x)) {
    .err("{.arg {arg}} must be a character vector")
  }

  arg_no_NA(x, arg = arg)
}

arg_string <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.character(x)) {
    .err("{.arg {arg}} must be a string")
  }

  if (length(x) != 1L) {
    .err("{.arg {arg}} must have length 1")
  }

  arg_no_NA(x, arg = arg)
}

arg_factor <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.factor(x)) {
    .err("{.arg {arg}} must be a factor")
  }
}

arg_flag <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.logical(x)) {
    .err("{.arg {arg}} must be a logical value ({.or {.val {c(TRUE, FALSE)}}})")
  }

  if (length(x) != 1L) {
    .err("{.arg {arg}} must have length 1")
  }

  arg_no_NA(x, arg = arg)
}

arg_logical <- function(x, .arg = rlang::caller_arg(x)) {
  if (!is.logical(x)) {
    .err("{.arg {(.arg)}} must be a logical vector")
  }
}

arg_list <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.list(x)) {
    .err("{.arg {arg}} must be a list")
  }
}

arg_is <- function(x, class, arg = rlang::caller_arg(x)) {
  if (!inherits(x, class)) {
    .err("{.arg {arg}} must inherit from class {.or {.cls {class}}}")
  }
}

arg_data <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.data.frame(x) && !is.matrix(x)) {
    .err("{.arg {arg}} must be a data frame or matrix")
  }
}

arg_data.frame <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.data.frame(x)) {
    .err("{.arg {arg}} must be a data frame")
  }
}

arg_matrix <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.matrix(x)) {
    .err("{.arg {arg}} must be a matrix")
  }
}

arg_function <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.function(x)) {
    .err("{.arg {arg}} must be a function")
  }
}

arg_subset <- function(x, values, arg = rlang::caller_arg(x)) {
  if (!all(x %in% values)) {
    if (length(x) == 1L) {
      .err("{.arg {arg}} must be one of {.or {.val {values}}}")
    }

    .err("each element of {.arg {arg}} must be one of {.or {.val {values}}}")
  }
}

arg_equal <- function(x1, x2, arg1 = rlang::caller_arg(x1), arg2 = rlang::caller_arg(x2), ...) {
  if (!isTRUE(all.equal(x1, x2, ...))) {
    .err("{.arg {arg1}} must be equal to {.arg {arg2}}")
  }
}

arg_formula <- function(x, one_sided = NULL, arg = rlang::caller_arg(x)) {

  if (length(one_sided) == 0L) {
    if (!rlang::is_formula(x)) {
      .err("{.arg {arg}} must be a formula")
    }
  }
  else {
    arg_flag(one_sided)

    if (one_sided) {
      if (!rlang::is_formula(x, lhs = FALSE)) {
        .err("{.arg {arg}} must be a one-sided formula")
      }
    }
    else if (!rlang::is_formula(x, lhs = TRUE)) {
      .err("{.arg {arg}} must be a two-sided formula")
    }
  }
}

arg_named <- function(x, arg = rlang::caller_arg(x)) {
  if (is_null(names(x)) || anyNA(names(x)) ||
      !all(nzchar(names(x)))) {
    .err("{.arg {arg}} must have non-empty and non-{.val {NA}} names")
  }
}

arg_length <- function(x, len = 1L, arg = rlang::caller_arg(x)) {
  arg_counts(len)

  if (!length(x) %in% len) {
    .err("{.arg {arg}} must have length {.or {len}}")
  }
}

arg_vector <- function(x, arg = rlang::caller_arg(x)) {
  if (!is.atomic(x) || !is_null(dim(x))) {
    .err("{.arg {arg}} must be a vector")
  }
}