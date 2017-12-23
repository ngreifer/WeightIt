weightitMSM <- function(formula.list, data, method = "ps", stabilize = FALSE, exact = NULL, s.weights = NULL,
                        num.formula = NULL, verbose = FALSE, ...) {

  A <- list(...)
  estimand <- "ATE"

  n <- nrow(data)

  #Process s.weights
  if (length(s.weights) > 0) {
    if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
      stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
    }
    if (is.character(s.weights) && length(s.weights)==1) {
      if (s.weights %in% names(data)) {
        s.weights <- data[[s.weights]]
      }
      else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
    }
    s.weights.specified <- TRUE
  }
  else {
    s.weights.specified <- FALSE
    s.weights <- rep(1, n)
  }

  call <- match.call()

  treat.type.vec <- character(length(formula.list))
  covs.list <- treat.list <- w.list <- ps.list <- vector("list", length(formula.list))
  for (i in seq_along(formula.list)) {
    #Process treat and covs from formula and data
    tt <- terms(formula.list[[i]])
    attr(tt, "intercept") <- 0
    if (is.na(match(all.vars(tt[[2]]), names(data)))) {
      stop(paste0("The given response variable, \"", all.vars(tt[[2]]), "\", is not a variable in data."))
    }
    vars.mentioned <- all.vars(tt)
    tryCatch({mf <- model.frame(tt, data, na.action = na.pass)}, error = function(e) {
      stop(paste0(c("All variables of formula ", i, " in formula.list must be variables in data.\nVariables not in data: ",
                    paste(vars.mentioned[is.na(match(vars.mentioned, names(data)))], collapse=", "))), call. = FALSE)})

    treat.list[[i]] <- model.response(mf)
    names(treat.list)[i] <- all.vars(tt[[2]])
    covs.list[[i]] <- data[!is.na(match(names(data), vars.mentioned[vars.mentioned != all.vars(tt[[2]])]))]

    if (any(is.na(covs.list[[i]]))) {
      stop("No missing values are allowed in the covariates.", call. = FALSE)
    }
    if (any(is.na(treat.list[[i]]))) {
      stop(paste0("No missing values are allowed in the treatment variable. Missing values found in ", names(treat.list[i])), call. = FALSE)
    }

    #Get weights into a list
    weightit_obj <- weightit(formula.list[[i]],
                             data = data,
                             method = method,
                             estimand = estimand,
                             stabilize = FALSE,
                             exact = exact,
                             s.weights = s.weights,
                             verbose = verbose, ...)
    w.list[[i]] <- cobalt::get.w(weightit_obj)
    if (length(weightit_obj$ps) > 0) ps.list[[i]] <- weightit_obj$ps
    treat.type.vec[i] <- weightit_obj$treat.type
  }

  w <- Reduce("*", w.list)

  if (length(num.formula) > 0) {
    if (!stabilize) {
      message("Setting stabilize to TRUE based on num.formula input.")
    }
    stabilize <- TRUE
  }
  if (stabilize) {
    stabilization.formula.list <- sw.list <- vector("list", length(formula.list))
    if (length(num.formula) > 0) {
      bad.num.formula <- FALSE
      if (is.list(num.formula)) {
        if (!all(sapply(num.formula, is.formula, sides = 2)) || length(num.formula) != length(formula.list)) {
          bad.num.formula <- TRUE
        }
        else if (!all(sapply(seq_along(num.formula), function(i) all.vars(num.formula[[i]][[2]]) == names(treat.list)[i]))) {
          stop("The response variable in each formula in num.formula must be the treatment variable at the same time point.", call. = FALSE)
        }
        else {
          for (i in seq_along(num.formula)) {
            tt <- terms(num.formula[[i]])
            attr(tt, "intercept") <- 0
            if (is.na(match(all.vars(tt[[2]]), names(data)))) {
              stop(paste0("The given response variable of formula ", i, " in num.formula, \"", all.vars(tt[[2]]), "\", is not a variable in data."))
            }
            vars.mentioned <- all.vars(tt)
            tryCatch({mf <- model.frame(tt, data)}, error = function(e) {
              stop(paste0(c("All variables of formula ", i, " in num.formula must be variables in data.\nVariables not in data: ",
                            paste(vars.mentioned[is.na(match(vars.mentioned, names(data)))], collapse=", "))), call. = FALSE)})
          }
          stabilization.formula.list <- num.formula
        }
      }
      else if (is.formula(num.formula)) {
        if (is.formula(num.formula, sides = 1)) {
          if (!all(all.vars(num.formula) %in% names(data))) {
            stop(paste0(c("All variables in num.formula must be variables in data.\nVariables not in data: ",
                          paste(all.vars(num.formula)[is.na(match(all.vars(num.formula), names(data)))], collapse=", "))), call. = FALSE)
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
          stop("The argument to num.formula must have rigt hand side variables but not a response variable (e.g., ~ V1 + V2).", call. = FALSE)
        }
      }
      else {
        bad.num.formula <- TRUE
      }

      if (bad.num.formula) {
        stop("The argument to num.formula must be a list of stabilization formulas for each time point or a single formula with no response variable and with the stabilization factors on the right hand side.", call. = FALSE)
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

      sw_obj <- weightit(stabilization.formula.list[[i]],
                         data = data,
                         method = "ps",
                         estimand = estimand,
                         stabilize = FALSE,
                         exact = exact,
                         s.weights = s.weights,
                         verbose = verbose, ...)
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
              covs.list = covs.list,
              data = data,
              estimand = estimand,
              method = method,
              ps.list = ps.list,
              s.weights = s.weights,
              #discarded = NULL,
              treat.type = treat.type.vec,
              call = call,
              stabilization = stabout
  )
  class(out) <- c("weightitMSM", "weightit")

  return(out)
}

print.weightitMSM <- function(x, ...) {
  cat("A weightitMSM object\n")
  cat(paste0(" - method: \"", x$method, "\" (", method.to.phrase(x$method), ")\n"))
  cat(paste0(" - number of obs.: ", length(x$weights), "\n"))
  cat(paste0(" - sampling weights: ", ifelse(max(x$s.weights) - min(x$s.weights) < sqrt(.Machine$double.eps),
                                             "none", "present"), "\n"))
  cat(paste0(" - number of time points: ", length(x$treat.list), " (", paste(names(x$treat.list), collapse = ", "), ")\n"))
  cat(paste0(" - treatment: \n",
             paste0(sapply(1:length(x$covs.list), function(i) {

               paste0("    + time ", i, ": ", ifelse(x$treat.type[i] == "continuous", "continuous", paste0(nunique(x$treat.list[[i]]), "-category", ifelse(x$treat.type[i] == "multinomial", paste0(" (", paste(levels(x$treat.list[[i]]), collapse = ", "), ")"), ""))), "\n")

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
  if (length(x$stabilization) > 0) cat(paste0(" - stabilization factors: ", if (length(x$stabilization) == 1) paste0(all.vars(x$stabilization[[1]]), collapse = ", ")
                                              else {
                                                paste0("\n", sapply(1:length(x$stabilization), function(i) {
                                                  if (i == 1) {
                                                    paste0("    + baseline: ", paste(all.vars(x$stabilization[[i]]), collapse = ", "))
                                                  }
                                                  else {
                                                    paste0("    + after time ", i-1, ": ", paste(all.vars(x$stabilization[[i]]), collapse = ", "))
                                                  }
                                                }), collapse = "")
                                              }))
  invisible(x)
}
summary.weightitMSM <- function(object, top = 5, ignore.s.weights = FALSE, ...) {
  outnames <- c("weight.range", "weight.top","weight.ratio",
                "coef.of.var",
                "effective.sample.size")
  out.list <- setNames(vector("list", length(object$treat.list)),
                       names(object$treat.list))

  if (ignore.s.weights) sw <- rep(1, length(object$weights))
  else sw <- object$s.weights
  w <- object$weights*sw

  for (ti in seq_along(object$treat.list)) {
    if (object$treat.type[ti] == "continuous") {
      out <- setNames(vector("list", length(outnames)), outnames)
      out$weight.range <- list(all = c(min(w[w > 0]),
                                       max(w[w > 0])))
      out$weight.ratio <- c(all = out$weight.range[["all"]][2]/out$weight.range[["all"]][1])
      top.weights <- sort(w, decreasing = TRUE)[seq_len(top)]
      out$weight.top <- list(all = sort(setNames(top.weights, which(w %in% top.weights)[seq_len(top)])))
      out$coef.of.var <- c(all = sd(w)/mean(w))

      nn <- as.data.frame(matrix(0, ncol = 1, nrow = 2))
      nn[1, ] <- (sum(sw)^2)/sum(sw^2)
      nn[2, ] <- (sum(w)^2)/sum((w)^2)
      dimnames(nn) <- list(c("Unweighted", "Weighted"),
                           c("Total"))
      out$effective.sample.size <- nn

      out.list[[ti]] <- out

    }
    else if (object$treat.type[ti] == "binary") {
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
    else if (object$treat.type[ti] == "multinomial") {

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
    else if (object$treat.type[ti] == "ordinal") {
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
    cat("\n- Effective Sample Sizes:\n")
    print.data.frame(round(x[[ti]]$effective.sample.size, 3))
    cat("\n")
    if (only.one) break
  }

  invisible(x)
}
