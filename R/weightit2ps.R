weightit2ps <- function(covs, treat, s.weights, subset, estimand, focal,
                        stabilize, subclass, missing, ps, .data, verbose, ...) {

  fit.obj <- NULL

  n <- length(treat)
  p.score <- NULL
  treat <- factor(treat)
  treat_sub <- factor(treat[subset])

  t.lev <- get_treated_level(treat, estimand, focal)
  c.lev <- setdiff(levels(treat_sub), t.lev)

  if (is.matrix(ps) || is.data.frame(ps)) {
    if (nrow(ps) == n) {
      if (ncol(ps) == 1L) {

        ps <- data.frame(ps[subset, 1L], 1 - ps[subset, 1L])

        names(ps) <- c(t.lev, c.lev)

        p.score <- ps[[t.lev]]
      }
      else if (ncol(ps) == 2L) {

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
  else if (is.numeric(ps) && length(ps) == n) {
    ps <- data.frame(ps[subset], 1 - ps[subset])

    names(ps) <- c(t.lev, c.lev)

    p.score <- ps[[t.lev]]
  }

  if (is_null(p.score)) {
    .err("`ps` must be a numeric vector with a propensity score for each unit")
  }

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
    else if (nrow(ps) == n && ncol(ps) == 1L) {
      ps <- setNames(list2DF(lapply(levels(treat), function(x) {
        p_ <- rep.int(1, length(treat))
        p_[treat == x] <- ps[treat == x, 1L]
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
  else {
    bad.ps <- TRUE
  }

  if (bad.ps) {
    .err("`ps` must be a numeric vector with a propensity score for each unit or a matrix \n\twith the probability of being in each treatment for each unit")
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_multi(ps = ps, treat = treat_sub, estimand, focal = focal,
                                     stabilize = stabilize, subclass = subclass)

  list(w = w)
}

weightit2ps.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, verbose, ...) {

  treat <- treat[subset]
  s.weights <- s.weights[subset]

  #Process density params
  densfun <- .get_dens_fun(use.kernel = isTRUE(...get("use.kernel")), bw = ...get("bw"),
                          adjust = ...get("adjust"), kernel = ...get("kernel"),
                          n = ...get("n"), treat = treat, density = ...get("density"),
                          weights = s.weights)

  #Stabilization - get dens.num
  log.dens.num <- densfun(scale_w(treat, s.weights), log = TRUE)

  #Get weights
  r <- treat - ps
  log.dens.denom <- densfun(r / sqrt(col.w.v(r, s.weights)), log = TRUE)

  w <- exp(log.dens.num - log.dens.denom)

  if (isTRUE(...get("plot"))) {
    d.n <- attr(log.dens.num, "density")
    d.d <- attr(log.dens.denom, "density")
    plot_density(d.n, d.d, log = TRUE)
  }

  list(w = w)
}
