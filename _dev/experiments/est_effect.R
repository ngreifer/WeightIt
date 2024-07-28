#Function to quickly estimate effect from glm_weightit() and weightit() output;
#similar to avg_comparisons()
est_effect <- function(fit, w, contrast = "diff", type = "response", treat = NULL) {
  treat_name <- deparse1(w$formula[[2]])

  preds <- unlist(insight::find_predictors(fit))
  model_data <- insight::get_data(fit, verbose = FALSE)

  if (treat_name %nin% preds) {
    .err("treatment not found")
  }

  estimand <- w$estimand
  s.weights <- w$s.weights
  weights <- w$weights * s.weights

  treat <- model_data[[treat_name]]
  treat_vals <- sort(unique(treat))

  if (estimand %in% c("ATT", "ATC")) {
    s.weights <- s.weights[treat == w$focal]
    weights <- weights[treat == w$focal]
    model_data <- model_data[treat == w$focal,]
  }

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

  combos <- combn(seq_along(treat_vals), 2, simplify = FALSE)

  est_fun <- function(b, fit) {
    fit$coefficients <- b

    preds <- lapply(treat_vals, function(t) {
      model_data[[treat_name]] <- t
      p <- predict(fit, newdata = model_data, type = type)

      if (is.matrix(p)) {
        apply(p, 2, mean_fun)
      }
      else {
        mean_fun(p)
      }
    })

    contrasts <- lapply(combos, function(co) {
      contrast_fun(preds[[co[1]]], preds[[co[2]]])
    })

    out <- c(unlist(preds), unlist(contrasts))

    attr(out, "preds") <- seq_along(unlist(preds))
    attr(out, "contrasts") <- seq_along(out)[-attr(out, "preds")]

    if (length(preds[[1]]) == 1L) {
      names(out)[attr(out, "preds")] <- sprintf("E[Y(%s)]", treat_vals)
      names(out)[attr(out, "contrasts")] <- unlist(lapply(combos, function(co) {
        sprintf("%s(%s, %s)", toupper(contrast),
                names(out)[co[1]], names(out)[co[2]])
      }))
    }
    else {

    }

    out
  }

  est <- est_fun(fit$coefficients, fit)

  est_grad <- .gradient(est_fun, fit$coefficients, fit = fit)

  vcov <- est_grad %*% tcrossprod(vcov(fit), est_grad)

  out <- list(est = est,
              vcov = vcov,
              estimand = estimand)
}