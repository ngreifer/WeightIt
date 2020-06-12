#For the user to use
make_full_rank <- function(mat, with.intercept = TRUE) {

  if (is.data.frame(mat)) {
    is.mat <- FALSE
    if (!all(vapply(mat, is.numeric, logical(1L)))) stop("All columns in 'mat' must be numeric.", call. = FALSE)
    mat <- as.matrix(mat)
  }
  else if (is.matrix(mat)) {
    is.mat <- TRUE
    if (!is.numeric(mat)) stop("'mat' must be a numeric matrix.", call. = FALSE)
  }
  else {
    stop("'mat' must be a numeric matrix or data.frame.", call. = FALSE)
  }

  if (anyNA(mat)) stop("Missing values are not allowed in 'mat'.", call. = FALSE)

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
  else return(as.data.frame(mat[, keep, drop = FALSE]))

}

get_w_from_ps <- function(ps, treat, estimand = "ATE", focal = NULL, treated = NULL, subclass = NULL, stabilize = FALSE) {
  #ps must be a matrix/df with columns named after treat levels

  if (!has.treat.type(treat)) treat <- assign.treat.type(treat)

  treat.type <- get.treat.type(treat)

  processed.estimand <- process.focal.and.estimand(focal, estimand, treat, treat.type, treated)
  estimand <- processed.estimand$estimand
  focal <- processed.estimand$focal
  assumed.treated <- processed.estimand$treated

  if (treat.type == "continuous") {
    stop("get_w_from_ps can only be used with binary or multinomial treatments.", call. = FALSE)
  }

  ps_mat <- ps_to_ps_mat(ps, treat, assumed.treated, treat.type, treated, estimand)

  if (nrow(ps_mat) != length(treat)) {
    stop("'ps' and 'treat' must have the same number of units.", call. = FALSE)
  }

  w <- rep(0, nrow(ps_mat))

  if (is_not_null(subclass)) {
    #Get MMW subclass propensity scores
    ps_mat <- subclass_ps(ps_mat, treat, estimand, focal, subclass)
  }

  for (i in colnames(ps_mat)) {
    w[treat == i] <- 1/ps_mat[treat == i, as.character(i)]
  }

  if (toupper(estimand) == "ATE") {
    # w <- w
  }
  else if (toupper(estimand) == "ATT") {
    w <- w*ps_mat[, as.character(focal)]
  }
  else if (toupper(estimand) == "ATO") {
    w <- w*rowSums(1/ps_mat)^-1 #Li & Li (2019)
  }
  else if (toupper(estimand) == "ATM") {
    w <- w*do.call("pmin", lapply(seq_len(ncol(ps_mat)), function(x) ps_mat[,x]), quote = TRUE)
  }
  else w <- NULL

  if (stabilize) w <- stabilize_w(w, treat)

  names(w) <- rownames(ps_mat)

  attr(w, "subclass") <- attr(ps_mat, "sub_mat")

  return(w)
}

