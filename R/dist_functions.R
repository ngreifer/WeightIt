#Functions to compute distance matrices
#Copied from MatchIt with edits

#Transforms covariates so that Euclidean distance computed on transforms covariates is equivalent to
#requested distance. When discarded is not NULL, statistics relevant to transformation are computed
#using retained units, but full covariate matrix is returned.
transform_covariates <- function(formula = NULL, data = NULL, method = "mahalanobis", s.weights = NULL, var = NULL,
                                 discarded = NULL) {
  X <- get.covs.matrix.for.dist(formula, data)

  X <- check_X(X)

  #If all variables have no variance, use Euclidean to avoid errors
  #If some have no variance, removes those to avoid messing up distances
  no_variance <- which(apply(X, 2, function(x) abs(max(x) - min(x)) < sqrt(.Machine$double.eps)))
  if (length(no_variance) == ncol(X)) {
    method <- "euclidean"
    X <- X[,1, drop = FALSE]
  }
  else if (length(no_variance) > 0) {
    X <- X[,-no_variance, drop = FALSE]
  }

  method <- match_arg(method, weightit_distances())

  if (is.null(discarded)) discarded <- rep.int(FALSE, nrow(X))

  if (method == "mahalanobis") {
    # X <- sweep(X, 2, colMeans(X))

    if (is.null(var)) {
      X <- scale(X)
      #NOTE: optmatch and Rubin (1980) use pooled within-group covariance matrix
      if (is.null(s.weights)) {
        var <- cov(X[!discarded,, drop = FALSE])
      }
      else {
        var <- cov.wt(X[!discarded,, drop = FALSE], s.weights[!discarded])$cov
      }
    }
    else if (!is.cov_like(var)) {
      stop("If 'var' is not NULL, it must be a covariance matrix with as many entries as supplied variables.", call. = FALSE)
    }

    inv_var <- NULL
    d <- det(var)
    if (d > 1e-8) {
      inv_var <- try(solve(var), silent = TRUE)
    }

    if (d <= 1e-8 || inherits(inv_var, "try-error")) {
      inv_var <- generalized_inverse(var)
    }

    X <- mahalanobize(X, inv_var)
  }
  else if (method == "robust_mahalanobis") {
    #Rosenbaum (2010, ch8)
    X_r <- matrix(0, nrow = sum(!discarded), ncol = ncol(X),
                  dimnames = list(rownames(X)[!discarded], colnames(X)))
    for (i in seq_col(X_r)) X_r[,i] <- rank(X[!discarded, i])

    if (is.null(s.weights)) var_r <- cov(X_r)
    else var_r <- cov.wt(X_r, s.weights[!discarded])$cov

    multiplier <- sd(seq_len(sum(!discarded)))/sqrt(diag(var_r))
    var_r <- var_r * outer(multiplier, multiplier, "*")

    inv_var <- NULL
    d <- det(var_r)
    if (d > 1e-8) {
      inv_var <- try(solve(var_r), silent = TRUE)
    }

    if (d <= 1e-8 || inherits(inv_var, "try-error")) {
      inv_var <- generalized_inverse(var_r)
    }

    if (any(discarded)) {
      X_r <- array(0, dim = dim(X), dimnames = dimnames(X))
      for (i in seq_col(X_r)) X_r[!discarded,i] <- rank(X[!discarded,i])
    }

    X <- mahalanobize(X_r, inv_var)
  }
  else if (method == "euclidean") {
    #Do nothing
  }
  else if (method == "scaled_euclidean") {
    if (is.null(var)) {
      sds <- sqrt(col.w.v(X[!discarded,, drop = FALSE], w = s.weights[!discarded]))
    }
    else if (is.cov_like(var, X)) {
      sds <- sqrt(diag(var))
    }
    else if (is.numeric(var) && is.cov_like(diag(var), X)) {
      sds <- sqrt(var)
    }
    else {
      stop("If 'var' is not NULL, it must be a covariance matrix or a vector of variances with as many entries as supplied variables.", call. = FALSE)
    }

    for (i in seq_col(X)) {
      X[,i] <- X[,i]/sds[i]
    }
  }

  X
}

#Internal function for fast(ish) Euclidean distance
eucdist_internal <- function(X, treat = NULL) {

  # if (is.null(dim(X))) X <- as.matrix(X)

  if (is.null(treat)) {
    if (is.null(dim(X))) {
      d <- abs(outer(X, X, "-"))
      dimnames(d) <- list(names(X), names(X))
    }
    else {
      if (ncol(X) == 1) {
        d <- abs(outer(drop(X), drop(X), "-"))
      }
      else {
        d <- as.matrix(dist(X))
      }
      dimnames(d) <- list(rownames(X), rownames(X))
    }
  }
  else {
    treat_l <- as.logical(treat)
    if (is.null(dim(X))) {
      d <- abs(outer(X[treat_l], X[!treat_l], "-"))
      dimnames(d) <- list(names(X)[treat_l], names(X)[!treat_l])
    }
    else {
      if (ncol(X) == 1) {
        d <- abs(outer(X[treat_l,], X[!treat_l,], "-"))
      }
      else {
        d <- as.matrix(dist(X))[treat_l, !treat_l, drop = FALSE]
      }
      dimnames(d) <- list(rownames(X)[treat_l], rownames(X)[!treat_l])
    }
  }

  d
}

#Get covariates (RHS) vars from formula; factor variable contrasts divided by sqrt(2)
#to ensure same result as when non-factor binary variable supplied (see optmatch:::contr.match_on)
get.covs.matrix.for.dist <- function(formula = NULL, data = NULL) {

  if (is.null(formula)) {
    fnames <- colnames(data)
    fnames[!startsWith(fnames, "`")] <- paste0("`", fnames[!startsWith(fnames, "`")], "`")
    data <- as.data.frame(data)
    formula <- terms(reformulate(fnames), data = data)
  }
  else {
    data <- as.data.frame(data)
    formula <- terms(formula, data = data)
  }

  mf <- model.frame(formula, data, na.action = na.pass)

  chars.in.mf <- vapply(mf, is.character, logical(1L))
  mf[chars.in.mf] <- lapply(mf[chars.in.mf], factor)

  formula <- update(formula, NULL ~ . + 1)

  X <- model.matrix(formula, data = mf,
                    contrasts.arg = lapply(Filter(is.factor, mf),
                                           function(x) contrasts(x, contrasts = FALSE)/sqrt(2)))

  if (ncol(X) > 1) {
    assign <- attr(X, "assign")[-1]
    X <- X[,-1, drop = FALSE]
  }
  attr(X, "assign") <- assign

  attr(X, "treat") <-  model.response(mf)

  X
}
check_X <- function(X) {
  if (isTRUE(attr(X, "checked"))) return(X)

  treat <- attr(X, "treat")

  if (is.data.frame(X)) X <- as.matrix(X)
  else if (is.numeric(X) && is.null(dim(X))) {
    X <- matrix(X, nrow = length(X),
                dimnames = list(names(X), NULL))
  }

  if (anyNA(X)) stop("Missing values are not allowed in the covariates.", call. = FALSE)
  else if (any(!is.finite(X))) stop("Non-finite values are not allowed in the covariates.", call. = FALSE)

  if (!is.numeric(X) || length(dim(X)) != 2) {
    stop("bad X")
  }
  attr(X, "checked") <- TRUE
  attr(X, "treat") <- treat
  X
}

is.cov_like <- function(var, X) {
  is.numeric(var) &&
    length(dim(var)) == 2 &&
    (missing(X) || all(dim(var) == ncol(X))) &&
    isSymmetric(var) &&
    all(diag(var) >= 0)
}
weightit_distances <- function() {
  c("mahalanobis", "robust_mahalanobis", "euclidean", "scaled_euclidean")
}
mahalanobize <- function(X, inv_var) {
  ## Mahalanobize covariates by computing cholesky decomp,
  ## allowing for NPD cov matrix by pivoting
  ch <- suppressWarnings(chol(inv_var, pivot = TRUE))
  p <- order(attr(ch, "pivot"))
  # r <- seq_len(attr(ch, "rank"))

  tcrossprod(X, ch[, p, drop = FALSE])
}
