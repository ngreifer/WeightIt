#Function to quickly estimate effect from glm_weightit() and weightit() output;
#similar to avg_comparisons()
est_effect <- function(fit, contrast = "diff", vcov = "unconditional", treat = NULL,
                       estimand = NULL, focal = NULL) {

  model_data <- insight::get_data(fit, verbose = FALSE)

  if (is_null(fit$weightit) || is_null(fit$weightit$formula)) {
    chk::chk_string(treat)
    predictors <- unlist(insight::find_predictors(fit))

    chk::chk_subset(treat, predictors)

    treat_name <- treat
  }
  else {
    treat_name <- deparse1(fit$weightit$formula[[2]])
  }

  treat <- model_data[[treat_name]]
  treat_vals <- sort(unique(treat))

  if (is_null(fit$weightit) || is_null(fit$weightit$estimand)) {
    chk::chk_string(estimand)
    estimand <- toupper(estimand)
  }
  else {
    estimand <- fit$weightit$estimand
  }

  if (estimand %in% c("ATT", "ATC")) {
    if (is_null(fit$weightit) || is_null(fit$weightit$focal)) {
      chk::chk_subset(focal, treat_vals)
    }
    else {
      focal <- fit$weightit$focal
    }
  }

  if (is_null(fit$weightit) || is_null(fit$weightit$s.weights)) {
    s.weights <- rep(1, length(treat))
  }
  else {
    s.weights <- fit$weightit$s.weights
  }

  if (is_null(fit$weightit) || is_null(fit$weightit$weights)) {
    weights <- rep(1, length(treat))
  }
  else {
    weights <- fit$weightit$weights * s.weights
  }

  if (estimand %in% c("ATT", "ATC")) {
    s.weights <- s.weights[treat == focal]
    weights <- weights[treat == focal]
    model_data <- model_data[treat == focal,]
  }

  chk::chk_string(vcov)
  vcov <- match_arg(vcov, c("unconditional", "conditional", "bootstrap", "fwb"))
  bs <- vcov %in% c("bootstrap", "fwb")

  pred_fun <- function(b, fit, t, data) {
    fit$coefficients <- b
    data[[treat_name]] <- t

    predict(fit, newdata = data, type = "response")
  }

  preds <- lapply(treat_vals, pred_fun,
                  b = fit$coefficients,
                  fit = fit, data = model_data)

  contrast_fun <- switch(contrast,
                         "diff" = function(b0, b1) b1 - b0,
                         "rr" = function(b0, b1) b1 / b0,
                         "or" = function(b0, b1) (b1 / (1 - b1)) / (b0 / (1 - b0)),
                         "sr" = function(b0, b1) (1 - b1) / (1 - b0))

  mean_fun <- {
    if (estimand %in% c("ATM", "ATO", "ATOS"))
      function(x) w.m(x, weights)
    else if (!all_the_same(s.weights))
      function(x) w.m(x, s.weights)
    else
      function(x) mean_fast(x)
  }

  Epreds <- vapply(preds, mean_fun, numeric(1L))

  names(Epreds) <- sprintf("E[Y(%s)]", treat_vals)

  combos <- combn(seq_along(treat_vals), 2, simplify = FALSE)

  est_fun <- function(mu) {
    out <- vapply(combos, function(combo) {
      contrast_fun(mu[combo[1]], mu[combo[2]])
    }, numeric(1L))

    names(out) <- vapply(combos, function(combo) {
      sprintf("%s[%s,%s]", toupper(contrast), treat_vals[combo[1]], treat_vals[combo[2]])
    }, character(1L))

    c(mu, out)
  }

  if (bs) {
    bootfun <- function(data, w) {
      boot_fit <- update(fit, weights = w * s.weights)

      boot_preds <- lapply(treat_vals, pred_fun,
                           b = boot_fit$coefficients,
                           fit = boot_fit, data = data)

      boot_Epreds <- vapply(boot_preds, w.m, numeric(1L),
                            w = {
                              if (estimand %in% c("ATT", "ATC", "ATE")) boot_fit$weightit$s.weights * w
                              else boot_fit$weightit$weights * w
                            })

      names(boot_Epreds) <- sprintf("E[Y(%s)]", treat_vals)

      est_fun(boot_Epreds)
    }

    boot_out <- fwb::fwb(model_data, bootfun, R = 999,
                         wtype = switch(vcov, "bootstrap" = "multinom", "exp"))

    est <- coef(boot_out)

    V <- vcov(boot_out)
  }
  else {
    dEpdb <- lapply(treat_vals, function(t) {
      apply(.gradient(pred_fun, .x = fit$coefficients,
                      fit = fit, t = t, data = model_data), 2, mean_fun)
    })
browser()
    beta_dot <- tcrossprod(sandwich::bread(fit),
                           sandwich::estfun(fit))

    if (vcov == "unconditional") {
      S1 <- do.call("rbind", lapply(seq_along(preds), function(i) {
        preds[[i]] - Epreds[i] + dEpdb[[i]] %*% beta_dot
      }))

      V <- tcrossprod(S1 / ncol(beta_dot))
    }
    else if (vcov == "conditional") {
      V <- do.call("rbind", dEpdb) %*% stats::vcov(fit) %*% t(do.call("rbind", dEpdb))
      # S1 <- do.call("rbind", lapply(seq_along(preds), function(i) {
      #   dEpdb[[i]] %*% beta_dot
      # }))
    }

    print(V)

    est_gr <- .gradient(est_fun, Epreds)

    vcov <- est_gr %*% V %*% t(est_gr)

    est <- est_fun(Epreds)
  }

  list(est = est,
       vcov = vcov,
       estimand = estimand)
}