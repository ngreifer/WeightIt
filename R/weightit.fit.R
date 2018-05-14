weightit.fit <- function(formula, data, method, treat.type, s.weights, exact.factor, estimand, focal, stabilize, ps, moments, int, ...){

  #main function of weightit that dispatches to weightit2method and returns object containing weights and ps
  out <- setNames(vector("list", 2), c("w", "ps"))

  for (i in levels(exact.factor)) {
    #Run method
    if (method == "ps") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2ps(formula = formula,
                           data = data,
                           s.weights = s.weights,
                           subset = exact.factor == i,
                           estimand = estimand,
                           focal = focal,
                           stabilize = stabilize,
                           ps = ps,
                           ...)
      }
      else if (treat.type == "continuous") {
        obj <- weightit2ps.cont(formula = formula,
                                data = data,
                                s.weights = s.weights,
                                subset = exact.factor == i,
                                stabilize = stabilize,
                                ps = ps,
                                ...)
      }
    }
    else if (method == "gbm") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2gbm(formula = formula,
                            data = data,
                            s.weights = s.weights,
                            estimand = estimand,
                            focal = focal,
                            subset = exact.factor == i,
                            stabilize = stabilize,
                            ...)
      }
      else stop("Generalized boosted modeling is not compatible with continuous treatments.", call. = FALSE)

    }
    else if (method == "cbps") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2cbps(formula = formula,
                             data = data,
                             subset = exact.factor == i,
                             s.weights = s.weights,
                             stabilize = stabilize,
                             estimand = estimand,
                             focal = focal,
                             ...)
      }
      else if (treat.type == "continuous") {
        obj <- weightit2cbps.cont(formula = formula,
                                  data = data,
                                  subset = exact.factor == i,
                                  s.weights = s.weights,
                                  #stabilize = stabilize,
                                  ...)

      }

    }
    else if (method == "npcbps") {
      if (s.weights.specified) stop(paste0("Sampling weights cannot be used with ", method.to.phrase(method), "."),
                                    call. = FALSE)
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2npcbps(formula = formula,
                               data = data,
                               subset = exact.factor == i,
                               ...)
      }
      else if (treat.type == "continuous") {
        obj <- weightit2npcbps.cont(formula = formula,
                                    data = data,
                                    subset = exact.factor == i,
                                    ...)
      }

    }
    else if (method == "ebal") {
      if (treat.type %in% c("binary", "multinomial")) {
        obj <- weightit2ebal(formula = formula,
                             data = data,
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
        obj <- weightit2sbw(formula = formula,
                            data = data,
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
        obj <- weightit2ebcw(formula = formula,
                             data = data,
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
    out$w[exact.factor == i] <- obj$w
    if (is_not_null(obj$ps)) out$ps[exact.factor == i] <- obj$ps


  }

  return(out)
}
