weightit.fit <- function(covs, treat, method, treat.type, s.weights, exact.factor, estimand, focal, stabilize, ps, moments, int, ...){

  #main function of weightit that dispatches to weightit2method and returns object containing weights and ps
  out <- setNames(vector("list", 2), c("w", "ps"))

  for (i in levels(exact.factor)) {
    #Run method
    if (method == "ps") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2ps(covs = covs,
                           treat = treat,
                           s.weights = s.weights,
                           subset = exact.factor == i,
                           estimand = estimand,
                           focal = focal,
                           stabilize = stabilize,
                           ps = ps,
                           ...)
      }
      else if (treat.type == "continuous") {
        obj <- weightit2ps.cont(covs = covs,
                                treat = treat,
                                s.weights = s.weights,
                                subset = exact.factor == i,
                                stabilize = stabilize,
                                ps = ps,
                                ...)
      }
    }
    else if (method == "gbm") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2gbm(covs = covs,
                            treat = treat,
                            s.weights = s.weights,
                            estimand = estimand,
                            focal = focal,
                            subset = exact.factor == i,
                            stabilize = stabilize,
                            ...)
      }
      else {
        obj <- weightit2gbm.cont(covs = covs,
                                treat = treat,
                                s.weights = s.weights,
                                subset = exact.factor == i,
                                stabilize = stabilize,
                                ...)
      }

    }
    else if (method == "cbps") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2cbps(covs = covs,
                             treat = treat,
                             subset = exact.factor == i,
                             s.weights = s.weights,
                             stabilize = stabilize,
                             estimand = estimand,
                             focal = focal,
                             ...)
      }
      else if (treat.type == "continuous") {
        obj <- weightit2cbps.cont(covs = covs,
                                  treat = treat,
                                  subset = exact.factor == i,
                                  s.weights = s.weights,
                                  #stabilize = stabilize,
                                  ...)

      }

    }
    else if (method == "npcbps") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2npcbps(covs = covs,
                               treat = treat,
                               subset = exact.factor == i,
                               s.weights = s.weights,
                               ...)
      }
      else if (treat.type == "continuous") {
        obj <- weightit2npcbps.cont(covs = covs,
                                    treat = treat,
                                    subset = exact.factor == i,
                                    s.weights = s.weights,
                                    ...)
      }

    }
    else if (method == "ebal") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2ebal(covs = covs,
                             treat = treat,
                             s.weights = s.weights,
                             subset = exact.factor == i,
                             estimand = estimand,
                             focal = focal,
                             stabilize = stabilize,
                             moments = moments,
                             int = int,
                             ...)
      }
      else stop("Entropy balancing is not compatible with continuous treatments.", call. = FALSE)
    }
    else if (method == "sbw") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2sbw(covs = covs,
                            treat = treat,
                            s.weights = s.weights,
                            subset = exact.factor == i,
                            estimand = estimand,
                            focal = focal,
                            moments = moments,
                            int = int,
                            ...)
      }
      else {
        stop("Stable balancing weights are not compatible with continuous treatments.", call. = FALSE)
      }
    }
    else if (method == "ebcw") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2ebcw(covs = covs,
                             treat = treat,
                             s.weights = s.weights,
                             subset = exact.factor == i,
                             estimand = estimand,
                             focal = focal,
                             #stabilize = stabilize,
                             moments = moments,
                             int = int,
                             ...)
      }
      else {
        stop("Empirical balancing calibration weights are not compatible with continuous treatments.", call. = FALSE)
      }
    }

    #Extract weights
    if (!exists("obj")) stop("No object was created. This is probably a bug,\n     and you should report it at https://github.com/ngreifer/WeightIt/issues.", call = FALSE)
    if (is_null(obj$w) || all(is.na(obj$w))) stop("No weights were estimated. This is probably a bug,\n     and you should report it at https://github.com/ngreifer/WeightIt/issues.", call = FALSE)
    if (any(is.na(obj$w))) warning("Some weights were estimated as NA, which means a value was impossible to compute (e.g., Inf). Check for extreme values of the treatment or covariates and try removing them. NA values will be set to 0.", call. = FALSE)
    obj$w[is.na(obj$w)] <- 0
    out$w[exact.factor == i] <- obj$w
    if (is_not_null(obj$ps)) out$ps[exact.factor == i] <- obj$ps


  }

  return(out)
}
