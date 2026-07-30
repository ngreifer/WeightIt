weightit2ps <- function(covs, treat, s.weights, subset, estimand, focal,
                        stabilize, missing, ps, .data, verbose, ...) {

  fit.obj <- NULL

  n <- length(treat)
  p.score <- NULL
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
    arg::err("{.arg ps} must be a numeric vector with a propensity score for each unit")
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_bin(ps = p.score,
                                   treat = as.numeric(treat_sub == t.lev), estimand,
                                   stabilize = stabilize,
                                   subclass = ...get("subclass"))

  list(w = w, ps = p.score, fit.obj = fit.obj)
}

weightit2ps.multi <- function(covs, treat, s.weights, subset, estimand, focal,
                              stabilize, missing, ps, .data, verbose, ...) {

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
        p_ <- rep_with(1, treat)
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
        p_ <- rep_with(1, treat)
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
    arg::err("{.arg ps} must be a numeric vector with a propensity score for each unit or a matrix with the probability of being in each treatment for each unit")
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- .get_w_from_ps_internal_multi(ps = ps, treat = treat_sub, estimand, focal = focal,
                                     stabilize = stabilize,
                                     subclass = ...get("subclass"))

  list(w = w)
}

weightit2ps.cens <- function(covs, treat, s.weights, subset, missing, ps, verbose,
                             estimand = NULL, focal = NULL, stabilize = FALSE, ...) {

  C <- .make_cens_treat(treat)

  out <- .cens_degenerate_out(C[subset])

  if (is_not_null(out)) {
    return(out)
  }

  #`ps` is the probability of being censored, so the censored units are the focal
  #("treated") group of the equivalent ATT problem. Delegating this way inherits
  #all of `weightit2ps()`'s input parsing.
  out <- weightit2ps(covs = covs, treat = C, s.weights = s.weights,
                     subset = subset, estimand = "ATT", focal = 1,
                     stabilize = FALSE, missing = missing, ps = ps,
                     verbose = verbose, ...)

  .att_out_to_cens(out, C[subset])
}

weightit2ps.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, verbose, ...) {

  treat <- treat[subset]
  s.weights <- s.weights[subset]

  # Process density params
  make_dens_fun <- .get_make_dens_fun(density = ...get("density"),
                                      bw = ...get("bw"),
                                      adjust = ...get("adjust"),
                                      kernel = ...get("kernel"),
                                      n = ...get("n"),
                                      use.kernel = ...get("use.kernel"))

  #Get weights
  w <- .get_w_from_gps_internal_cont(mu = ps, treat = treat,
                                     s.weights = s.weights,
                                     make_dens_fun = make_dens_fun)

  list(w = w)
}
