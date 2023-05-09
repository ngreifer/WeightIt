#For the user to use
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

get_w_from_ps <- function(ps, treat, estimand = "ATE", focal = NULL, treated = NULL, subclass = NULL, stabilize = FALSE) {
  #ps must be a matrix/df with columns named after treat levels

  if (!has_treat_type(treat)) treat <- assign_treat_type(treat)
  treat.type <- get_treat_type(treat)

  if (treat.type == "continuous") {
    .err("`get_w_from_ps()` can only be used with binary or multinomial treatments")
  }

  estimand <- process.estimand(estimand, method = "glm", treat.type = treat.type)

  processed.estimand <- process.focal.and.estimand(focal, estimand, treat, treated)
  estimand <- processed.estimand$estimand
  focal <- processed.estimand$focal
  assumed.treated <- processed.estimand$treated

  ps_mat <- ps_to_ps_mat(ps, treat, assumed.treated, treat.type, treated, estimand)

  if (nrow(ps_mat) != length(treat)) {
    .err("`ps` and `treat` must have the same number of units")
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
    w <- w*do.call("pmin", lapply(seq_col(ps_mat), function(x) ps_mat[,x]), quote = TRUE)
  }
  else if (toupper(estimand) == "ATOS") {
    #Crump et al. (2009)
    ps.sorted <- sort(c(ps_mat[,2], 1 - ps_mat[,2]))
    q <- ps_mat[,2]*(1-ps_mat[,2])
    alpha.opt <- 0
    for (i in 1:sum(ps_mat[,2] < .5)) {
      if (i == 1 || !check_if_zero(ps.sorted[i] - ps.sorted[i-1])) {
        alpha <- ps.sorted[i]
        a <- alpha*(1-alpha)
        if (1/a <= 2*sum(1/q[q >= a])/sum(q >= a)) {
          alpha.opt <- alpha
          break
        }
      }
    }
    w[!between(ps_mat[,2], c(alpha.opt, 1 - alpha.opt))] <- 0
  }
  else return(numeric(0))

  if (stabilize) w <- stabilize_w(w, treat)

  names(w) <- if_null_then(rownames(ps_mat), names(treat), NULL)

  attr(w, "subclass") <- attr(ps_mat, "sub_mat")
  if (toupper(estimand) == "ATOS") attr(w, "alpha") <- alpha.opt

  w
}

