weightit2ps <- function(covs, treat, s.weights, subset, estimand, focal,
                        stabilize, subclass, missing, ps, .data, verbose, ...) {
  A <- list(...)

  fit.obj <- NULL

  n <- length(treat)
  p.score <- NULL
  treat <- factor(treat)
  treat_sub <- factor(treat[subset])

  t.lev <- get_treated_level(treat)
  c.lev <- setdiff(levels(treat_sub), t.lev)

  if (is.matrix(ps) || is.data.frame(ps)) {
    if (nrow(ps) == n) {
      if (ncol(ps) == 1) {

        ps <- data.frame(ps[subset,1], 1-ps[subset,1])

        names(ps) <- c(t.lev, c.lev)

        p.score <- ps[[t.lev]]
      }
      else if (ncol(ps) == 2) {

        if (all(colnames(ps) %in% levels(treat_sub))) {
          ps <- as.data.frame(ps[subset, , drop = FALSE])
        }
        else {
          ps <- as.data.frame(ps[subset, , drop = FALSE])
          names(ps) <- levels(treat_sub)
        }

        p.score <- ps[[t.lev]]
      }
    }
  }
  else if (is.numeric(ps)) {
    if (length(ps) == n) {
      ps <- data.frame(ps[subset], 1-ps[subset])

      names(ps) <- c(t.lev, c.lev)

      p.score <- ps[[t.lev]]
    }
  }

  if (is_null(p.score)) .err("`ps` must be a numeric vector with a propensity score for each unit")

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_bin(ps = p.score, treat = as.numeric(treat_sub == t.lev), estimand,
                                   stabilize = stabilize, subclass = subclass)

  list(w = w, ps = p.score, fit.obj = fit.obj)
}

weightit2ps.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                              stabilize, subclass, missing, ps, .data, verbose, ...) {

  n <- length(treat)
  treat <- factor(treat)
  treat_sub <- factor(treat[subset])

  bad.ps <- FALSE
  if (is.matrix(ps) || is.data.frame(ps)) {
    if (all(dim(ps) == c(n, nunique(treat)))) {
      ps <- setNames(as.data.frame(ps), levels(treat))[subset, , drop = FALSE]
    }
    else if (all(dim(ps) == c(n, 1))) {
      ps <- setNames(list2DF(lapply(levels(treat), function(x) {
        p_ <- rep.int(1, length(treat))
        p_[treat == x] <- ps[treat == x, 1]
        p_
      })), levels(treat))[subset, , drop = FALSE]
    }
    else {
      bad.ps <- TRUE
    }
  }
  else if (is.numeric(ps)) {
    if (length(ps) == n) {
      ps <- setNames(list2DF(lapply(levels(treat), function(x) {
        p_ <- rep.int(1, length(treat))
        p_[treat == x] <- ps[treat == x]
        p_
      })), levels(treat))[subset, , drop = FALSE]
    }
    else {
      bad.ps <- TRUE
    }
  }
  else bad.ps <- TRUE

  if (bad.ps) .err("`ps` must be a numeric vector with a propensity score for each unit or a matrix \n\twith the probability of being in each treatment for each unit")

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_multi(ps = ps, treat = treat_sub, estimand, focal = focal,
                                     stabilize = stabilize, subclass = subclass)

  list(w = w)
}

weightit2ps.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, verbose, ...) {
  A <- list(...)

  treat <- treat[subset]
  s.weights <- s.weights[subset]

  #Process density params
  densfun <- .get_dens_fun(use.kernel = isTRUE(A[["use.kernel"]]), bw = A[["bw"]],
                          adjust = A[["adjust"]], kernel = A[["kernel"]],
                          n = A[["n"]], treat = treat, density = A[["density"]],
                          weights = s.weights)

  #Stabilization - get dens.num
  dens.num <- densfun(scale_w(treat, s.weights))

  #Get weights
  r <- treat - ps
  dens.denom <- densfun(r / sqrt(col.w.v(r, s.weights)))

  w <- dens.num / dens.denom

  if (isTRUE(A[["plot"]])) {
    d.n <- attr(dens.num, "density")
    d.d <- attr(dens.denom, "density")
    plot_density(d.n, d.d)
  }

  list(w = w)
}
