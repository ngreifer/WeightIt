weightitMSM <- function(formula.list, data = NULL, method = "glm", stabilize = FALSE, by = NULL, s.weights = NULL,
                        num.formula = NULL, moments = NULL, int = FALSE, missing = NULL,
                        verbose = FALSE, include.obj = FALSE, is.MSM.method, weightit.force = FALSE, ...) {

  A <- list(...)

  call <- match.call()

  ## Checks and processing ----

  #Checks

  ##Process method
  check.acceptable.method(method, msm = TRUE, force = weightit.force)

  if (is.character(method)) {
    method <- method.to.proper.method(method)
    attr(method, "name") <- method
    if (missing(is.MSM.method)) is.MSM.method <- NULL
    is.MSM.method <- process.MSM.method(is.MSM.method, method)
  }
  else if (is.function(method)) {
    method.name <- paste(deparse(substitute(method)))
    check.user.method(method)
    if (missing(is.MSM.method)) is.MSM.method <- NULL
    is.MSM.method <- process.MSM.method(is.MSM.method, method)
    attr(method, "name") <- method.name
  }

  #Process moments and int
  moments.int <- process.moments.int(moments, int, method)
  moments <- moments.int[["moments"]]; int <- moments.int[["int"]]

  s.weights <- process.s.weights(s.weights, data)

  if (is_not_null(num.formula)) {
    if (!stabilize) {
      .msg("setting `stabilize` to `TRUE` based on `num.formula` input")
    }
    stabilize <- TRUE
  }
  if (stabilize) {
    if (is_not_null(num.formula)) {
      if (rlang::is_formula(num.formula)) {
        if (!rlang::is_formula(num.formula, lhs = FALSE)) {
          .err("the argument to `num.formula` must have right hand side variables but not a response variable (e.g., ~ V1 + V2)")
        }

        rhs.vars.mentioned.lang <- attr(terms(num.formula), "variables")[-1]
        rhs.vars.mentioned <- vapply(rhs.vars.mentioned.lang, deparse1, character(1L))
        rhs.vars.failed <- vapply(rhs.vars.mentioned.lang, function(v) {
          null_or_error(try(eval(v, c(data, .GlobalEnv)), silent = TRUE))
        }, logical(1L))

        if (any(rhs.vars.failed)) {
          .err(paste0(c("All variables in `num.formula` must be variables in `data` or objects in the global environment.\nMissing variables: ",
                        paste(rhs.vars.mentioned[rhs.vars.failed], collapse = ", "))), tidy = FALSE)
        }
      }
      else if (is.list(num.formula)) {
        if (length(num.formula) != length(formula.list)) {
          .err("when supplied as a list, `num.formula` must have as many entries as `formula.list`", call. = FALSE)
        }
        if (!all(vapply(num.formula, rlang::is_formula, logical(1L), lhs = FALSE))) {
          .err("'num.formula' must be a single formula with no response variable and with the stabilization factors on the right hand side or a list thereof")
        }
        rhs.vars.mentioned.lang.list <- lapply(num.formula, function(nf) attr(terms(nf), "variables")[-1])
        rhs.vars.mentioned <- unique(unlist(lapply(rhs.vars.mentioned.lang.list, function(r) vapply(r, deparse1, character(1L)))))
        rhs.vars.failed <- vapply(rhs.vars.mentioned, function(v) {
          null_or_error(try(eval(parse(text=v), c(data, .GlobalEnv)), silent = TRUE))
        }, logical(1L))

        if (any(rhs.vars.failed)) {
          .err(paste0(c("All variables in `num.formula` must be variables in `data` or objects in the global environment.\nMissing variables: ",
                        paste(rhs.vars.mentioned[rhs.vars.failed], collapse=", "))), tidy = FALSE)
        }

      }
      else {
        .err("`num.formula` must be a single formula with no response variable and with the stabilization factors on the right hand side or a list thereof")
      }
    }
  }

  ##Process by
  if (is_not_null(A[["exact"]])) {
    .msg("`by` has replaced `exact` in the `weightit()` syntax, but `exact` will always work")
    # by.name <- deparse(A[["exact"]])
    by <- A[["exact"]]
    by.arg <- "exact"
  }
  else by.arg <- "by"

  reported.covs.list <- covs.list <- treat.list <- w.list <- ps.list <-
    stabout <- sw.list <- make_list(length(formula.list))

  if (is_null(formula.list) || !is.list(formula.list) || !all(vapply(formula.list, rlang::is_formula, logical(1L), lhs = TRUE))) {
    .err("`formula.list` must be a list of formulas")
  }

  for (i in seq_along(formula.list)) {

    #Process treat and covs from formula and data
    t.c <- get.covs.and.treat.from.formula(formula.list[[i]], data)
    reported.covs.list[[i]] <- t.c[["reported.covs"]]
    covs.list[[i]] <- t.c[["model.covs"]]
    treat.list[[i]] <- t.c[["treat"]]
    treat.name <- t.c[["treat.name"]]
    names(treat.list)[i] <- treat.name
    names(reported.covs.list)[i] <- treat.name

    if (is_null(covs.list[[i]])) .err(sprintf("no covariates were specified in the %s formula", ordinal(i)))
    if (is_null(treat.list[[i]])) .err(sprintf("no treatment variable was specified in the %s formula", ordinal(i)))

    n <- length(treat.list[[i]])

    if (nrow(covs.list[[i]]) != n) {
      .err("treatment and covariates must have the same number of units")
    }
    if (anyNA(treat.list[[i]])) {
      .err(sprintf("no missing values are allowed in the treatment variable. Missing values found in %s", treat.name))
    }
    if (anyNA(reported.covs.list[[i]])) {
      .wrn("missing values are present in the covariates. See `?weightit` for information on how these are handled")
    }

    treat.list[[i]] <- assign.treat.type(treat.list[[i]])

    #By is processed each for each time, but only last time is used for by.factor.
    processed.by <- process.by(by, data = data,
                               treat = treat.list[[i]],
                               treat.name = treat.name,
                               by.arg = by.arg)

    #Process missing
    if (anyNA(reported.covs.list[[i]])) {
      missing <- process.missing(missing, method, get.treat.type(treat.list[[i]]))
    }
    else if (i == length(formula.list)) missing <- ""

  }

  if (is_null(s.weights)) s.weights <- rep(1, n)

  if (is.MSM.method) {
    #Returns weights (w)

    A[["covs"]] <- covs.list
    A[["treat"]] <- treat.list
    A[["s.weights"]] <- s.weights
    A[["by.factor"]] <- attr(processed.by, "by.factor")
    A[["focal"]] <- character()
    A[["stabilize"]] <- stabilize
    A[["method"]] <- method
    A[["moments"]] <- moments
    A[["int"]] <- int
    A[["subclass"]] <- numeric()
    A[["ps"]] <- numeric()
    A[["missing"]] <- missing
    A[["verbose"]] <- verbose
    A[["is.MSM.method"]] <- TRUE
    A[["include.obj"]] <- include.obj

    obj <- do.call("weightit.fit", A)

    w <- obj[["weights"]]
    stabout <- NULL
    obj.list <- obj[["fit.obj"]]
  }
  else {
    if (length(A[["link"]]) %nin% c(0, 1, length(formula.list))) {
      .err(sprintf("the argument to `link` must have length 1 or %s", length(formula.list)))
    }
    if (length(A[["link"]]) == 1) A[["link"]] <- rep(A[["link"]], length(formula.list))
    # if (length(A[["family"]]) %nin% c(0, 1, length(formula.list))) stop(paste0("The argument to link must have length 1 or ", length(formula.list), "."), call. = FALSE)
    # if (length(A[["family"]]) == 1) A[["family"]] <- rep(A[["family"]], length(formula.list))

    obj.list <- make_list(length(formula.list))

    A[["s.weights"]] <- s.weights
    A[["by.factor"]] <- attr(processed.by, "by.factor")
    A[["estimand"]] <- "ATE"
    A[["focal"]] <- character()
    A[["stabilize"]] <- FALSE
    A[["method"]] <- method
    A[["moments"]] <- moments
    A[["int"]] <- int
    A[["subclass"]] <- numeric()
    A[["ps"]] <- numeric()
    A[["missing"]] <- missing
    A[["verbose"]] <- verbose
    A[["is.MSM.method"]] <- FALSE
    A[["include.obj"]] <- include.obj

    for (i in seq_along(formula.list)) {
      A_i <- A
      if (length(A[["link"]]) == length(formula.list)) A_i[["link"]] <- A[["link"]][[i]]

      A_i[["covs"]] <- covs.list[[i]]
      A_i[["treat"]] <- treat.list[[i]]
      A_i[["treat.type"]] <- get.treat.type(treat.list[[i]])
      A_i[[".data"]] <- data
      A_i[[".covs"]] <- reported.covs.list[[i]]

      ## Running models ----

      #Returns weights (w) and propensty score (ps)
      obj <- do.call("weightit.fit", A_i)

      w.list[[i]] <- obj[["weights"]]
      ps.list[[i]] <- obj[["ps"]]
      obj.list[[i]] <- obj[["fit.obj"]]

      if (stabilize) {
        #Process stabilization formulas and get stab weights
        if (rlang::is_formula(num.formula)) {
          if (i == 1) {
            stab.f <- update.formula(as.formula(paste(names(treat.list)[i], "~ 1")), as.formula(paste(paste(num.formula, collapse = ""), "+ .")))
          }
          else {
            stab.f <- update.formula(as.formula(paste(names(treat.list)[i], "~", paste(names(treat.list)[seq_along(names(treat.list)) < i], collapse = " * "))), as.formula(paste(num.formula, "+ .")))
          }
        }
        else if (is.list(num.formula)) {
          stab.f <- update.formula(as.formula(paste(names(treat.list)[i], "~ 1")), as.formula(paste(paste(num.formula[[i]], collapse = ""), "+ .")))
        }
        else {
          if (i == 1) {
            stab.f <- as.formula(paste(names(treat.list)[i], "~ 1"))
          }
          else {
            stab.f <- as.formula(paste(names(treat.list)[i], "~", paste(names(treat.list)[seq_along(names(treat.list)) < i], collapse = " * ")))
          }
        }
        stab.t.c_i <- get.covs.and.treat.from.formula(stab.f, data)

        A_i[["covs"]] <- stab.t.c_i[["model.covs"]]
        A_i[["method"]] <- "glm"
        A_i[["moments"]] <- numeric()
        A_i[["int"]] <- FALSE

        sw_obj <- do.call("weightit.fit", A_i)

        sw.list[[i]] <- 1/sw_obj[["weights"]]
        stabout[[i]] <- stab.f[-2]

      }

    }

    w <- Reduce("*", w.list)

    if (stabilize) {
      sw <- Reduce("*", sw.list)
      w <- w*sw

      unique.stabout <- unique(stabout)
      if (length(unique.stabout) <= 1) stabout <- unique.stabout
    }
    else stabout <- NULL
  }


  if (all_the_same(w)) .err(sprintf("All weights are %s", w[1]))
  if (all(vapply(ps.list, is_null, logical(1L)))) ps.list <- NULL
  else names(ps.list) <- names(treat.list)

  if (include.obj) names(obj.list) <- names(treat.list)

  ## Assemble output object----
  out <- list(weights = w,
              treat.list = treat.list,
              covs.list = reported.covs.list,
              #data = data,
              estimand = "ATE",
              method = method,
              ps.list = ps.list,
              s.weights = s.weights,
              #discarded = NULL,
              by = processed.by,
              call = call,
              stabilization = stabout,
              obj = obj.list
  )

  out <- clear_null(out)

  class(out) <- c("weightitMSM", "weightit")

  out
  ####----
}

print.weightitMSM <- function(x, ...) {
  treat.types <- vapply(x[["treat.list"]], get.treat.type, character(1L))
  trim <- attr(x[["weights"]], "trim")

  cat("A " %+% italic("weightitMSM") %+% " object\n")
  cat(paste0(" - method: \"", attr(x[["method"]], "name"), "\" (", method.to.phrase(x[["method"]]), ")\n"))
  cat(paste0(" - number of obs.: ", length(x[["weights"]]), "\n"))
  cat(paste0(" - sampling weights: ", ifelse(all_the_same(x[["s.weights"]]), "none", "present"), "\n"))
  cat(paste0(" - number of time points: ", length(x[["treat.list"]]), " (", paste(names(x[["treat.list"]]), collapse = ", "), ")\n"))
  cat(paste0(" - treatment: \n",
             paste0(vapply(seq_along(x$covs.list), function(i) {
               paste0("    + time ", i, ": ", if (treat.types[i] == "continuous") "continuous"
                      else paste0(nunique(x[["treat.list"]][[i]]), "-category",
                                  if (treat.types[i] == "multinomial") paste0(" (", paste(levels(x[["treat.list"]][[i]]), collapse = ", "), ")")
                                  else ""
                      ), "\n")
             }, character(1L)), collapse = ""), collapse = "\n"))
  cat(paste0(" - covariates: \n",
             paste0(vapply(seq_along(x$covs.list), function(i) {
               if (i == 1) {
                 paste0("    + baseline: ", if (is_null(x$covs.list[[i]])) "(none)" else paste(names(x$covs.list[[i]]), collapse = ", "), "\n")
               }
               else {
                 paste0("    + after time ", i-1, ": ", paste(names(x$covs.list[[i]]), collapse = ", "), "\n")
               }
             }, character(1L)), collapse = ""), collapse = "\n"))
  if (is_not_null(x[["by"]])) {
    cat(paste0(" - by: ", paste(names(x[["by"]]), collapse = ", "), "\n"))
  }
  if (is_not_null(x$stabilization)) {
    cat(" - stabilized")
    if (any(vapply(x$stabilization, function(s) is_not_null(all.vars(s)), logical(1L)))) {
      cat(paste0("; stabilization factors:\n", if (length(x$stabilization) == 1) paste0("      ", paste0(attr(terms(x[["stabilization"]][[1]]), "term.labels"), collapse = ", "))
                 else {
                   paste0(vapply(seq_along(x$stabilization), function(i) {
                     if (i == 1) {
                       paste0("    + baseline: ", if (is_null(attr(terms(x[["stabilization"]][[i]]), "term.labels"))) "(none)" else paste(attr(terms(x[["stabilization"]][[i]]), "term.labels"), collapse = ", "))
                     }
                     else {
                       paste0("    + after time ", i-1, ": ", paste(attr(terms(x[["stabilization"]][[i]]), "term.labels"), collapse = ", "))
                     }
                   }, character(1L)), collapse = "\n")
                 }))
    }
  }

  if (is_not_null(trim)) {
    if (trim < 1) {
      if (attr(x[["weights"]], "trim.lower")) trim <- c(1 - trim, trim)
      cat(paste(" - weights trimmed at", word_list(paste0(round(100*trim, 2), "%")), "\n"))
    }
    else {
      if (attr(x[["weights"]], "trim.lower")) t.b <- "top and bottom" else t.b <- "top"
      cat(paste(" - weights trimmed at the", t.b, trim, "\n"))
    }
  }
  invisible(x)
}
summary.weightitMSM <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  outnames <- c("weight.range", "weight.top",
                "coef.of.var", "scaled.mad", "negative.entropy",
                "weight.mean",
                "effective.sample.size")

  out.list <- make_list(names(object$treat.list))

  if (ignore.s.weights || is_null(object$s.weights)) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights
  w <- setNames(object$weights*sw, seq_along(sw))
  treat.types <- vapply(object[["treat.list"]], get.treat.type, character(1L))
  stabilized <- is_not_null(object[["stabilization"]])

  for (ti in seq_along(object$treat.list)) {
    out <- make_list(outnames)
    if (treat.types[ti] == "continuous") {
      out$weight.range <- list(all = c(min(w[w != 0]),
                                       max(w[w != 0])))
      out$weight.top <- list(all = rev(w[order(abs(w), decreasing = TRUE)][seq_len(top)]))
      out$coef.of.var <- c(all = sd(w)/mean_fast(w))
      out$scaled.mad <- c(all = mean.abs.dev(w/mean_fast(w)))
      out$negative.entropy <- c(all = neg_ent(w))
      out$num.zeros <- c(overall = sum(check_if_zero(w)))
      out$weight.mean <- if (stabilized) mean_fast(w) else NULL

      nn <- make_df("Total", c("Unweighted", "Weighted"))
      nn["Unweighted", ] <- ESS(sw)
      nn["Weighted", ] <- ESS(w)

      out$effective.sample.size <- nn

      out.list[[ti]] <- out

    }
    else if (treat.types[ti] == "binary") {
      t <- object$treat.list[[ti]]
      top0 <- c(treated = min(top, sum(t == 1)),
                control = min(top, sum(t == 0)))
      out$weight.range <- list(treated = c(min(w[w != 0 & t == 1]),
                                           max(w[w != 0 & t == 1])),
                               control = c(min(w[w != 0 & t == 0]),
                                           max(w[w != 0 & t == 0])))
      out$weight.top <- list(treated = rev(w[t == 1][order(abs(w[t == 1]), decreasing = TRUE)][seq_len(top0["treated"])]),
                             control = rev(w[t == 0][order(abs(w[t == 0]), decreasing = TRUE)][seq_len(top0["control"])]))
      out$coef.of.var <- c(treated = sd(w[t==1])/mean_fast(w[t==1]),
                           control = sd(w[t==0])/mean_fast(w[t==0]))
      out$scaled.mad <- c(treated = mean.abs.dev(w[t==1]/mean_fast(w[t==1])),
                          control = mean.abs.dev(w[t==0]/mean_fast(w[t==0])))
      out$negative.entropy <- c(treated = neg_ent(w[t==1]),
                                control = neg_ent(w[t==0]))
      out$num.zeros <- c(treated = sum(check_if_zero(w[t==1])),
                         control = sum(check_if_zero(w[t==0])))
      out$weight.mean <- if (stabilized) mean_fast(w) else NULL

      nn <- make_df(c("Control", "Treated"), c("Unweighted", "Weighted"))
      nn["Unweighted", ] <- c(ESS(sw[t==0]),
                              ESS(sw[t==1]))
      nn["Weighted", ] <- c(ESS(w[t==0]),
                            ESS(w[t==1]))

      out$effective.sample.size <- nn
      out.list[[ti]] <- out

    }
    else if (treat.types[ti] == "multinomial") {
      t <- object$treat.list[[ti]]
      top0 <- setNames(lapply(levels(t), function(x) min(top, sum(t == x))), levels(t))
      out$weight.range <- setNames(lapply(levels(t), function(x) c(min(w[w != 0 & t == x]),
                                                                   max(w[w != 0 & t == x]))),
                                   levels(t))
      out$weight.top <- setNames(lapply(levels(t), function(x) rev(w[t == x][order(abs(w[t == x]), decreasing = TRUE)][seq_len(top0[[x]])])),
                                 levels(t))
      out$coef.of.var <- c(vapply(levels(t), function(x) sd(w[t==x])/mean_fast(w[t==x]), numeric(1L)),
                           overall = sd(w)/mean_fast(w))
      out$scaled.mad <- c(vapply(levels(t), function(x) mean.abs.dev(w[t==x]/mean_fast(w[t==x])), numeric(1L)),
                          overall = mean.abs.dev(w)/mean_fast(w))
      out$negative.entropy <- c(vapply(levels(t), function(x) neg_ent(w[t==x]), numeric(1L)),
                                overall = sum(w[w>0]*log(w[w>0]))/sum(w[w>0]))
      out$num.zeros <- c(vapply(levels(t), function(x) sum(check_if_zero(w[t==x])), numeric(1L)),
                         overall = sum(check_if_zero(w)))
      out$weight.mean <- if (stabilized) mean_fast(w) else NULL

      nn <- make_df(levels(t), c("Unweighted", "Weighted"))
      for (i in levels(t)) {
        nn["Unweighted", i] <- ESS(sw[t==i])
        nn["Weighted", i] <- ESS(w[t==i])
      }

      out$effective.sample.size <- nn
      out.list[[ti]] <- out
    }
    else if (treat.types[ti] == "ordinal") {
      .wrn("Sneaky, sneaky! Ordinal coming soon :)", tidy = FALSE)
    }
  }

  class(out.list) <- "summary.weightitMSM"
  attr(out.list, "weights") <- w

  out.list
}
summary.weightitMSM <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  out.list <- make_list(names(object$treat.list))

  if (ignore.s.weights || is_null(object$s.weights)) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights

  for (ti in seq_along(object$treat.list)) {
    obj <- as.weightit(weights = object$weights, treat = object$treat.list[[ti]],
                       s.weights = sw, stabilization = object$stabilization)
    out.list[[ti]] <- summary.weightit(obj, top = top, ignore.s.weights = ignore.s.weights, ...)
  }

  class(out.list) <- "summary.weightitMSM"

  out.list
}
print.summary.weightitMSM <- function(x, ...) {
  only.one <- all(vapply(x, function(y) isTRUE(all.equal(x[[1]], y)), logical(1L)))

  cat(paste(rep(" ", 17), collapse = "") %+% underline("Summary of weights") %+% "\n\n")
  for (ti in seq_along(x)) {
    if (!only.one) cat(strikethrough(paste(rep(" ", 22), collapse = "")) %+% italic(" Time " %+% ti %+% " ") %+% strikethrough(paste(rep(" ", 22), collapse = "")) %+% "\n")
    print(x[[ti]])
    cat("\n")
    if (only.one) break
  }

  invisible(x)
}
plot.summary.weightitMSM <- function(x, binwidth = NULL, bins = NULL, time = 1, ...) {
  if (!is.numeric(time) || length(time) != 1 || time %nin% seq_along(x)) {
    .err("`time` must be a number corresponding to the time point for which to display the distribution of weights")
  }
  p <- plot.summary.weightit(x[[time]], binwidth = binwidth, bins = bins, ...)
  p + labs(subtitle = paste0("For Time ", time))
}