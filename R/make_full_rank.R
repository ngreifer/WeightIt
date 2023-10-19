#' Make a design matrix full rank
#'
#' @description
#' When writing [user-defined methods][method_user] for use with
#' [weightit()], it may be necessary to take the potentially non-full rank
#' `covs` data frame and make it full rank for use in a downstream
#' function. This function performs that operation.
#'
#' @param mat a numeric matrix or data frame to be transformed. Typically this
#' contains covariates. `NA`s are not allowed.
#' @param with.intercept whether an intercept (i.e., a vector of 1s) should be
#' added to `mat` before making it full rank. If `TRUE`, the
#' intercept will be used in determining whether a column is linearly dependent
#' on others. Regardless, no intercept will be included in the output.
#'
#' @returns
#' An object of the same type as `mat` containing only linearly
#' independent columns.
#'
#' @details
#' `make_full_rank()` calls [qr()] to find the rank and linearly
#' independent columns of `mat`, which are retained while others are
#' dropped. If `with.intercept` is set to `TRUE`, an intercept column
#' is added to the matrix before calling `qr()`. Note that dependent
#' columns that appear later in `mat` will be dropped first.
#'
#' See example at [`method_user`].
#'
#' @note
#' Older versions would drop all columns that only had one value. With
#' `with.intercept = FALSE`, if only one column has only one value, it
#' will not be removed, and it will function as though there was an intercept
#' present; if more than only column has only one value, only the first one
#' will remain.
#'
#' @seealso
#' [`method_user`], [model.matrix()]
#'
#' @examples
#'
#' set.seed(1000)
#' c1 <- rbinom(10, 1, .4)
#' c2 <- 1-c1
#' c3 <- rnorm(10)
#' c4 <- 10*c3
#' mat <- data.frame(c1, c2, c3, c4)
#'
#' make_full_rank(mat) #leaves c2 and c4
#'
#' make_full_rank(mat, with.intercept = FALSE) #leaves c1, c2, and c4

#' @export
make_full_rank <- function(mat, with.intercept = TRUE) {

  if (is.data.frame(mat)) {
    is.mat <- FALSE
    if (!all(vapply(mat, is.numeric, logical(1L)))) {
      .err("all columns in `mat` must be numeric")
    }
    mat <- as.matrix(mat)
  }
  else if (is.matrix(mat)) {
    if (!is.numeric(mat)) .err("`mat` must be a numeric matrix")
    is.mat <- TRUE
  }
  else {
    .err("`mat` must be a numeric matrix or data.frame")
  }

  chk::chk_not_any_na(mat)

  keep <- rep(TRUE, ncol(mat))

  #If intercept is to be included in check, add column of 1s
  if (with.intercept) {
    q <- qr(cbind(1, mat))
    keep[q$pivot[-seq(q$rank)]-1] <- FALSE
  }
  else {
    q <- qr(mat)
    keep[q$pivot[-seq(q$rank)]] <- FALSE
  }

  if (is.mat) return(mat[, keep, drop = FALSE])

  as.data.frame(mat[, keep, drop = FALSE])
}
