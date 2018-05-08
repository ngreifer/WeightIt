weightitMSM <- function(formula.list, data = NULL, method = "ps", stabilize = FALSE, exact = NULL, s.weights = NULL,
                        num.formula = NULL, moments = 1L, int = FALSE,
                        verbose = FALSE, ...) {

  A <- list(...)
  A[["estimand"]] <- NULL
  estimand <- "ATE"

  #Process s.weights
  if (is_not_null(s.weights)) {
    if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
      stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
    }
    if (is.character(s.weights) && length(s.weights)==1) {
      if (is_null(data)) {
        stop("s.weights was specified as a string but there was no argument to data.", call. = FALSE)
      }
      else if (s.weights %in% names(data)) {
        s.weights <- data[[s.weights]]
      }
      else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
    }
    s.weights.specified <- TRUE
  }
  else {
    s.weights.specified <- FALSE
  }

  call <- match.call()

  covs.list <- reported.covs.list <- treat.list <- data.list <- w.list <- ps.list <- vector("list", length(formula.list))
  for (i in seq_along(formula.list)) {
    #Process treat and covs from formula and data
    t.c <- get.covs.and.treat.from.formula(formula.list[[i]], data)
    reported.covs.list[[i]] <- t.c[["reported.covs"]]
    covs.list[[i]] <- t.c[["model.covs"]]
    treat.list[[i]] <- t.c[["treat"]]
    names(treat.list)[i] <- t.c[["treat.name"]]

    if (is_null(covs.list[[i]])) stop(paste0("No covariates were specified in the ", ordinal(i), "formula."), call. = FALSE)
    if (is_null(treat.list[[i]])) stop(paste0("No treatment variable was specified in the ", ordinal(i), "formula."), call. = FALSE)

    data.list[[i]] <- data.frame(treat.list[[i]], covs.list[[i]])
    formula.list[[i]] <- formula(data.list[[i]])

    if (any(is.na(treat.list[[i]]))) {
      stop(paste0("No missing values are allowed in the treatment variable. Missing values found in ", names(treat.list[i])), call. = FALSE)
    }

    if (!s.weights.specified) s.weights <- rep(1, length(treat.list[[i]]))

    #Get weights into a list
    weightit_obj <- do.call("weightit", c(list(formula.list[[i]],
                             data = data.list[[i]],
                             method = method,
                             estimand = estimand,
                             stabilize = FALSE,
                             exact = exact,
                             s.weights = s.weights,
                             verbose = verbose,
                             moments = moments,
                             int = int),
                             A))
    w.list[[i]] <- weightit_obj[["weights"]]
    if (is_not_null(weightit_obj[["ps"]])) ps.list[[i]] <- weightit_obj[["ps"]]
    attr(treat.list[[i]], "treat.type") <- attr(weightit_obj[["treat"]], "treat.type")
  }

  w <- Reduce("*", w.list)

  if (is_not_null(num.formula)) {
    if (!stabilize) {
      message("Setting stabilize to TRUE based on num.formula input.")
    }
    stabilize <- TRUE
  }
  if (stabilize) {
    stabilization.formula.list <- sw.list <- vector("list", length(formula.list))
    if (is_not_null(num.formula)) {
      bad.num.formula <- FALSE
      if (is.formula(num.formula)) {
        if (is.formula(num.formula, sides = 1)) {
          rhs.vars.mentioned <- attr(terms(num.formula), "term.labels")
          rhs.vars.failed <- sapply(rhs.vars.mentioned, function(v) {
            null.or.error(try(eval(parse(text = v), c(data, .GlobalEnv)), silent = TRUE))
          })
          if (any(rhs.vars.failed)) {
            stop(paste0(c("All variables in num.formula must be variables in data or objects in the global environment.\nMissing variables: ",
                          paste(rhs.vars.mentioned[rhs.vars.failed], collapse=", "))), call. = FALSE)
          }
          else {
            for (i in seq_along(formula.list)) {
              if (i == 1) {
                stabilization.formula.list[[i]] <- update.formula(as.formula(paste(names(treat.list)[i], "~ 1")), as.formula(paste(paste(num.formula, collapse = ""), "+ .")))
              }
              else {
                stabilization.formula.list[[i]] <- update.formula(as.formula(paste(names(treat.list)[i], "~", paste(names(treat.list)[seq_along(names(treat.list)) < i], collapse = " * "))), as.formula(paste(num.formula, "+ .")))
              }
            }
          }
        }
        else {
          stop("The argument to num.formula must have right hand side variables but not a response variable (e.g., ~ V1 + V2).", call. = FALSE)
        }
      }
      else if (is.list(num.formula)) {
        bad.num.formula <- TRUE
        # if (!all(sapply(num.formula, is.formula, sides = 2)) || length(num.formula) != length(formula.list)) {
        #   bad.num.formula <- TRUE
        # }
        # else if (!all(sapply(seq_along(num.formula), function(i) all.vars(num.formula[[i]][[2]]) == names(treat.list)[i]))) {
        #   stop("The response variable in each formula in num.formula must be the treatment variable at the same time point.", call. = FALSE)
        # }
        # else {
        #   for (i in seq_along(num.formula)) {
        #     t.c <- get.covs.and.treat.from.formula(num.formula[[i]], data)
        #     num.covs.list[[i]] <- t.c[["covs"]]
        #     num.treat.list[[i]] <- t.c[["treat"]]
        #
        #     if (is_null(covs.list[[i]])) stop(paste0("No covariates were specified in the ", ordinal(i), "num.formula."), call. = FALSE)
        #     if (is_null(treat.list[[i]])) stop(paste0("No treatment variable was specified in the ", ordinal(i), "num.formula."), call. = FALSE)
        #
        #     num.data_i <- data.frame(num.treat.list[[i]], num.covs.list[[i]])
        #     stabilization.formula.list[[i]] <- formula(num.data_i)
        #
        #   }
        # }
      }
      else {
        bad.num.formula <- TRUE
      }

      if (bad.num.formula) {
        #stop("The argument to num.formula must be a list of stabilization formulas for each time point or a single formula with no response variable and with the stabilization factors on the right hand side.", call. = FALSE)
        stop("The argument to num.formula must be a single formula with no response variable and with the stabilization factors on the right hand side.", call. = FALSE)
      }
    }
    else {
      for (i in seq_along(formula.list)) {
        if (i == 1) {
          stabilization.formula.list[[i]] <- as.formula(paste(names(treat.list)[i], "~ 1"))
        }
        else {
          stabilization.formula.list[[i]] <- as.formula(paste(names(treat.list)[i], "~", paste(names(treat.list)[seq_along(names(treat.list)) < i], collapse = " * ")))
        }
      }
    }
    for (i in seq_along(formula.list)) {
      sw_obj <- do.call("weightit", c(list(stabilization.formula.list[[i]],
                         data = data,
                         method = "ps",
                         estimand = estimand,
                         stabilize = FALSE,
                         exact = exact,
                         s.weights = s.weights,
                         verbose = verbose),
                         A))
      sw.list[[i]] <- 1/cobalt::get.w(sw_obj)

    }
    sw <- Reduce("*", sw.list)
    w <- w*sw

    stabout <- lapply(stabilization.formula.list, function(x) x[-2])
    unique.stabout <- unique(stabout)
    if (length(unique.stabout) <= 1) stabout <- unique.stabout
  }
  else stabout <- NULL

  #if (nunique(treat.type.vec) == 1) treat.type.vec <- unique(treat.type.vec)

  ## Assemble output object----
  out <- list(weights = w,
              treat.list = treat.list,
              covs.list = reported.covs.list,
              data = data,
              estimand = estimand,
              method = method,
              ps.list = ps.list,
              s.weights = s.weights,
              #discarded = NULL,
              call = call,
              stabilization = stabout
  )

  class(out) <- c("weightitMSM", "weightit")

  return(out)
}

print.weightitMSM <- function(x, ...) {
  treat.types <- sapply(x[["treat.list"]], function(y) attr(y, "treat.type"))
  trim <- attr(x[["weights"]], "trim")

  cat("A weightitMSM object\n")
  cat(paste0(" - method: \"", x[["method"]], "\" (", method.to.phrase(x[["method"]]), ")\n"))
  cat(paste0(" - number of obs.: ", length(x[["weights"]]), "\n"))
  cat(paste0(" - sampling weights: ", ifelse(all_the_same(x[["s.weights"]]), "none", "present"), "\n"))
  cat(paste0(" - number of time points: ", length(x[["treat.list"]]), " (", paste(names(x[["treat.list"]]), collapse = ", "), ")\n"))
  cat(paste0(" - treatment: \n",
             paste0(sapply(1:length(x$covs.list), function(i) {

               paste0("    + time ", i, ": ", ifelse(treat.types[i] == "continuous", "continuous", paste0(nunique(x[["treat.list"]][[i]]), "-category", ifelse(treat.types[i] == "multinomial", paste0(" (", paste(levels(x[["treat.list"]][[i]]), collapse = ", "), ")"), ""))), "\n")

             }), collapse = ""), collapse = "\n"))
  cat(paste0(" - covariates: \n",
             paste0(sapply(1:length(x$covs.list), function(i) {
               if (i == 1) {
                 paste0("    + baseline: ", paste(names(x$covs.list[[i]]), collapse = ", "), "\n")
               }
               else {
                 paste0("    + after time ", i-1, ": ", paste(names(x$covs.list[[i]])[!names(x$covs.list[[i]]) %in% c(names(x$treat.list)[i-1], names(x$covs.list[[i-1]]))], collapse = ", "), "\n")
               }
             }), collapse = ""), collapse = "\n"))
  if (is_not_null(x$stabilization)) {
    cat(" - stabilized")
    if (any(sapply(x$stabilization, function(s) is_not_null(all.vars(s))))) {
      cat(paste0("; stabilization factors:\n", if (length(x$stabilization) == 1) paste0("      ", paste0(attr(terms(x[["stabilization"]][[1]]), "term.labels"), collapse = ", "))
                 else {
                   paste0("\n", sapply(1:length(x$stabilization), function(i) {
                     if (i == 1) {
                       paste0("    + baseline: ", paste(attr(terms(x[["stabilization"]][[i]]), "term.labels"), collapse = ", "))
                     }
                     else {
                       paste0("    + after time ", i-1, ": ", paste(attr(terms(x[["stabilization"]][[i]]), "term.labels"), collapse = ", "))
                     }
                   }), collapse = "")
                 }))
    }
  }

  if (is_not_null(trim)) {
    if (trim < 1) {
      if (attr(x[["weights"]], "trim.lower")) trim <- c(1 - trim, trim)
      cat(paste(" - weights trimmed at", word.list(paste0(round(100*trim, 2), "%")), "\n"))
    }
    else {
      if (attr(x[["weights"]], "trim.lower")) t.b <- "top and bottom" else t.b <- "top"
      cat(paste(" - weights trimmed at the", t.b, trim, "\n"))
    }
  }
  invisible(x)
}
summary.weightitMSM <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  outnames <- c("weight.range", "weight.top","weight.ratio",
                "coef.of.var", "weight.mean",
                "effective.sample.size")
  out.list <- setNames(vector("list", length(object$treat.list)),
                       names(object$treat.list))

  if (ignore.s.weights) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights
  w <- object$weights*sw
  treat.types <- sapply(object[["treat.list"]], function(y) attr(y, "treat.type"))
  stabilized <- is_not_null(object[["stabilization"]])

  for (ti in seq_along(object$treat.list)) {
    if (treat.types[ti] == "continuous") {
      out <- setNames(vector("list", length(outnames)), outnames)
      out$weight.range <- list(all = c(min(w[w > 0]),
                                       max(w[w > 0])))
      out$weight.ratio <- c(all = out$weight.range[["all"]][2]/out$weight.range[["all"]][1])
      top.weights <- sort(w, decreasing = TRUE)[seq_len(top)]
      out$weight.top <- list(all = sort(setNames(top.weights, which(w %in% top.weights)[seq_len(top)])))
      out$coef.of.var <- c(all = sd(w)/mean(w))
      out$weight.mean <- if (stabilized) mean(w) else NULL

      nn <- as.data.frame(matrix(0, ncol = 1, nrow = 2))
      nn[1, ] <- (sum(sw)^2)/sum(sw^2)
      nn[2, ] <- (sum(w)^2)/sum((w)^2)
      dimnames(nn) <- list(c("Unweighted", "Weighted"),
                           c("Total"))
      out$effective.sample.size <- nn

      out.list[[ti]] <- out

    }
    else if (treat.types[ti] == "binary") {
      out <- setNames(vector("list", length(outnames)), outnames)
      t <- object$treat.list[[ti]]
      out$weight.range <- list(treated = c(min(w[w > 0 & t == 1]),
                                           max(w[w > 0 & t == 1])),
                               control = c(min(w[w > 0 & t == 0]),
                                           max(w[w > 0 & t == 0])))
      out$weight.ratio <- c(treated = out$weight.range$treated[2]/out$weight.range$treated[1],
                            control = out$weight.range$control[2]/out$weight.range$control[1],
                            overall = max(unlist(out$weight.range)/min(unlist(out$weight.range))))
      top.weights <- list(treated = sort(w[t == 1], decreasing = TRUE)[seq_len(top)],
                          control = sort(w[t == 0], decreasing = TRUE)[seq_len(top)])
      out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w[t == ifelse(x == "control", 0, 1)] %in% top.weights[[x]])[seq_len(top)]))),
                                 names(top.weights))
      out$coef.of.var <- c(treated = sd(w[t==1])/mean(w[t==1]),
                           control = sd(w[t==0])/mean(w[t==0]),
                           overall = sd(w)/mean(w))
      out$weight.mean <- if (stabilized) mean(w) else NULL

      #dc <- weightit$discarded

      nn <- as.data.frame(matrix(0, nrow = 2, ncol = 2))
      nn[1, ] <- c((sum(sw[t==0])^2)/sum(sw[t==0]^2),
                   (sum(sw[t==1])^2)/sum(sw[t==1]^2))
      nn[2, ] <- c((sum(w[t==0])^2)/sum((w[t==0])^2),
                   (sum(w[t==1])^2)/sum((w[t==1])^2))
      # nn[3, ] <- c(sum(t==0 & dc==1), #Discarded
      #              sum(t==1 & dc==1))
      dimnames(nn) <- list(c("Unweighted", "Weighted"),
                           c("Control", "Treated"))
      out$effective.sample.size <- nn
      out.list[[ti]] <- out

    }
    else if (treat.types[ti] == "multinomial") {

      out <- setNames(vector("list", length(outnames)), outnames)
      t <- object$treat.list[[ti]]
      out$weight.range <- setNames(lapply(levels(t), function(x) c(min(w[w > 0 & t == x]),
                                                                   max(w[w > 0 & t == x]))),
                                   levels(t))
      out$weight.ratio <- setNames(c(sapply(out$weight.range, function(x) x[2]/x[1]),
                                     max(unlist(out$weight.range)/min(unlist(out$weight.range)))),
                                   c(levels(t), "overall"))
      top.weights <- setNames(lapply(levels(t), function(x) sort(w[t == x], decreasing = TRUE)[seq_len(top)]),
                              levels(t))
      out$weight.top <- setNames(lapply(names(top.weights), function(x) sort(setNames(top.weights[[x]], which(w[t == x] %in% top.weights[[x]])[seq_len(top)]))),
                                 names(top.weights))
      out$coef.of.var <- c(sapply(levels(t), function(x) sd(w[t==x])/mean(w[t==x])),
                           overall = sd(w)/mean(w))
      out$weight.mean <- if (stabilized) mean(w) else NULL

      nn <- as.data.frame(matrix(0, nrow = 2, ncol = nunique(t)))
      for (i in seq_len(nunique(t))) {
        nn[1, i] <- (sum(sw[t==levels(t)[i]])^2)/sum(sw[t==levels(t)[i]]^2)
        nn[2, i] <- (sum(w[t==levels(t)[i]])^2)/sum((w[t==levels(t)[i]])^2)
        # nn[3, i] <- sum(t==levels(t)[i] & dc==1) #Discarded
      }
      dimnames(nn) <- list(c("Unweighted", "Weighted"),
                           levels(t))
      out$effective.sample.size <- nn
      out.list[[ti]] <- out
    }
    else if (treat.types[ti] == "ordinal") {
      warning("Sneaky, sneaky! Ordinal coming soon :)", call. = FALSE)
    }
  }

  class(out.list) <- "summary.weightitMSM"
  return(out.list)
}
print.summary.weightitMSM <- function(x, ...) {
  if (all(sapply(x, function(y) isTRUE(all.equal(x[[1]], y))))) {
    only.one <- TRUE
  }
  else only.one <- FALSE

  cat("Summary of weights:\n\n")
  for (ti in seq_along(x)) {
    if (!only.one) cat(paste0(" - - - - - - - - - - Time ", ti, " - - - - - - - - - -\n"))
    cat("- Weight ranges:\n")
    print.data.frame(round_df(text.box.plot(x[[ti]]$weight.range, 28), 4))

    df <- setNames(data.frame(do.call("c", lapply(names(x[[ti]]$weight.top), function(y) c(" ", y))),
                              matrix(do.call("c", lapply(x[[ti]]$weight.top, function(y) c(names(y), round(y, 4)))),
                                     byrow = TRUE, nrow = 2*length(x[[ti]]$weight.top))),
                   rep("", 1 + length(x[[ti]]$weight.top[[1]])))
    cat(paste("\n- Units with", length(x[[ti]]$weight.top[[1]]), "greatest weights by group:\n"))
    print.data.frame(df, row.names = FALSE)
    cat("\n")
    print.data.frame(as.data.frame(round(matrix(c(x[[ti]]$weight.ratio, x[[ti]]$coef.of.var), ncol = 2,
                                                dimnames = list(names(x[[ti]]$weight.ratio),
                                                                c("Ratio", "Coef of Var"))), 4)))
    if (is_not_null(x[[ti]][["weight.mean"]])) cat("\n- Mean of Weights =", round(x[[ti]][["weight.mean"]], 4), "\n")

    cat("\n- Effective Sample Sizes:\n")
    print.data.frame(round(x[[ti]]$effective.sample.size, 3))
    cat("\n")
    if (only.one) break
  }

  invisible(x)
}
