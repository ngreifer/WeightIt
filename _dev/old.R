#Old

#Stable balancing weights with sbw
weightit2sbw <- function(...) {
  stop("method = \"sbw\" has been deprecated in place of method = \"optweight\". Please use that instead.", call. = FALSE)
  #stop("Stable balancing weights are not currently supported. Please choose another method.\n        The github version of WeightIt may allow stable balancing weights.\n        Install it with devtools::install_github(\"ngreifer/WeightIt\").", call. = FALSE)
}

#SBW--------
weightit2sbw <- function(covs, treat, s.weights, subset, estimand, focal, moments, int, ...) {
  A <- list(...)

  if (check.package("sbw")) {
    check.package("slam")
    if (!"package:slam" %in% search()) {
      need.to.detach.slam <- TRUE
      attachNamespace("slam")
    } else {
      need.to.detach.slam <- FALSE
    }

    if (is_null(A$l_norm)) A$l_norm <- "l_2"
    if (is_null(A$solver)) A$solver <- "quadprog"
    if (is_null(A$max_iter)) A$max_iter <- 100000
    if (is_null(A$rel_tol)) A$rel_tol <- 1e-4
    if (is_null(A$abs_tol)) A$abs_tol <- 1e-4
    if (is_null(A$gap_stop)) A$gap_stop <- TRUE
    if (is_null(A$adaptive_rho)) A$adaptive_rho <- TRUE

    solve.package = switch(A$solver, quadprog = "quadprog",
                           cplex = "Rcplex",
                           gurobi = "gurobi",
                           pogs = "pogs")
    check.package(solve.package)
    if (!paste0("package:", solve.package) %in% search()) {
      need.to.detach.solver <- TRUE
      attachNamespace(solve.package)
    } else {
      need.to.detach.solver <- FALSE
    }

    covs <- covs[subset, , drop = FALSE]
    treat <- factor(treat[subset])

    if (anyNA(covs)) {
      stop("Stable balancing weights are not compatible with missing values.", call. = FALSE)
    }
    covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int))
    covs <- apply(covs, 2, make.closer.to.1)

    #new.data <- setNames(data.frame(treat, covs), as.character(seq_len(1+ncol(covs))))

    binary.vars <- apply(covs, 2, is_binary)

    if (is_null(A$bal_tols)) bal_tols <- .0001
    else {
      bal_tols <- A$bal_tols
      if (length(bal_tols) != 1 && length(bal_tols) != ncol(covs)) {
        stop(paste0("bal_tols needs to be of length 1 or equal to the number of covariates (", ncol(covs),
                    ").\nThe covariates (in order) are:\n   ", paste0(colnames(covs), collapse = " ")), call.= FALSE)
      }
    }
    if (is_null(A$bal_tols_sd)) bal_tols_sd <- TRUE
    else bal_tols_sd <- A$bal_tols_sd

    if (is_null(A$bal_tols) && is_null(A$bal_tols_sd)) {
      message("Using bal_tols = 0.0001 and bal_tols_sd = TRUE.")
    }
    else if (is_null(A$bal_tols)) {
      message("Using bal_tols = 0.0001.")
    }
    else if (is_null(A$bal_tols_sd)) {
      message("Using bal_tols_sd = TRUE.")
    }

    if (estimand == "ATT") {
      control.levels <- levels(treat)[levels(treat) != focal]

      w <- rep(1, length(treat))

      if (bal_tols_sd) {
        bal_tols <- bal_tols * sapply(1:ncol(covs), function(x) {if (binary.vars[x]) 1 else sqrt(cov.wt(covs[treat == focal, x, drop = FALSE], s.weights[subset][treat == focal])$cov[1,1])})
      }

      covs[treat == focal,] <- covs[treat == focal,] * s.weights[subset][treat == focal] * sum(treat == focal)/sum(s.weights[subset][treat == focal])

      for (i in control.levels) {

        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- ifelse(treat[treat.in.i.focal] == i, 0, 1)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]

        new.data_ <- data.frame(treat_, covs_)
        t_ind <- names(new.data_)[1]
        bal_covs = names(new.data_)[-1]

        sbw.fit <- sbw::sbw(new.data_,
                            t_ind = t_ind,
                            bal_covs = bal_covs,
                            bal_tols = bal_tols,
                            bal_tols_sd = FALSE,
                            target = "treated",
                            l_norm = A[["l_norm"]],
                            w_min = 0,
                            normalize = TRUE,
                            solver = A[["solver"]],
                            display = 1,
                            max_iter = A[["max_iter"]],
                            rel_tol = A[["rel_tol"]],
                            abs_tol = A[["abs_tol"]],
                            gap_stop = A[["gap_stop"]],
                            adaptive_rho = A[["adaptive_rho"]])

        w[treat==i] <- sbw.fit$data_frame_weights$weights[sbw.fit$data_frame_weights[[t_ind]] == 0]*sum(treat == i) / s.weights[subset][treat == i]
      }
      w[w < 0] <- 0
    }
    else if (estimand == "ATE") {
      if (bal_tols_sd) {
        bal_tols <- bal_tols * sapply(1:ncol(covs), function(x) {if (binary.vars[x]) 1 else sqrt(mean(sapply(unique(treat), function(t) cov.wt(covs[treat == t, x, drop = FALSE], s.weights[subset][treat == t])$cov[1,1])))})
      }

      bal_tols <- bal_tols/nunique(treat)

      w <- rep(1, length(treat))

      for (i in levels(treat)) {
        covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
        treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

        covs_i[treat_i == 1,] <- covs_i[treat_i == 1, , drop = FALSE] * s.weights[subset] * sum(treat_i == 1) / sum(s.weights[subset])

        new.data_i <- data.frame(treat_i, covs_i)
        t_ind <- names(new.data_i)[1]
        bal_covs = names(new.data_i)[-1]

        sbw.fit_i <- sbw::sbw(new.data_i, t_ind = t_ind,
                              bal_covs = bal_covs,
                              bal_tols = bal_tols,
                              bal_tols_sd = FALSE,
                              target = "treated",
                              l_norm = A$l_norm,
                              w_min = 0,
                              normalize = TRUE,
                              solver = A[["solver"]],
                              display = 1,
                              max_iter = A[["max_iter"]],
                              rel_tol = A[["rel_tol"]],
                              abs_tol = A[["abs_tol"]],
                              gap_stop = A[["gap_stop"]],
                              adaptive_rho = A[["adaptive_rho"]])

        w[treat==i] <- sbw.fit_i$data_frame_weights$weights[sbw.fit_i$data_frame_weights[[t_ind]] == 0]*sum(treat == i) / s.weights[subset][treat == i]

      }
      w[w < 0] <- 0
    }
  }

  if (need.to.detach.slam) detach("package:slam", character.only = TRUE)
  if (need.to.detach.solver) detach(paste0("package:", solve.package), character.only = TRUE)

  obj <- list(w = w)
  return(obj)

}

#BART using BART package (now using dbarts)----
.weightit2bart <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, subclass, missing, ...) {
  A <- list(...)

  check.package("BART")

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (!all_the_same(s.weights)) stop("Sampling weights cannot be used with method = \"bart\".",
                                     call. = FALSE)

  if (!has_treat_type(treat)) treat <- assign_treat_type(treat)
  treat.type <- get_treat_type(treat)

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  if (is_null(A$link)) A$link <- "probit"
  else {
    acceptable.links <- c("probit", "logit")

    which.link <- acceptable.links[pmatch(A$link, acceptable.links, nomatch = 0)][1]
    if (is.na(which.link)) {
      A$link <- acceptable.links[1]
      warning(paste0("Only ", word_list(acceptable.links, quotes = TRUE, is.are = TRUE), " allowed as the link for ",
                     if (treat.type == "binary") "binary" else "multinomial",
                     " treatments. Using link = ", word_list(acceptable.links[1], quotes = TRUE), "."),
              call. = FALSE, immediate. = TRUE)
    }
    else A$link <- which.link
  }


  if (is_not_null(A[["n.trees"]]) && is_null(A[["ntree"]])) A[["ntree"]] <- A[["n.trees"]]

  ps <- make_df(levels(treat), nrow = length(treat))

  if (treat.type == "binary") {

    if (is_not_null(A[["mc.cores"]]) && A[["mc.cores"]] > 1) fun <- BART::mc.gbart
    else fun <- BART::gbart

    fit.list <- do.call(fun, c(list(covs,
                                    y.train = binarize(treat),
                                    type = switch(A$link, "logit" = "lbart", "pbart")),
                               A[intersect(names(A), setdiff(names(formals(fun)),
                                                             c("x.train", "y.train", "x.test", "type", "ntype")))]))

    p.score <- fit.list$prob.train.mean

    ps[[2]] <- p.score
    ps[[1]] <- 1 - ps[[2]]

    info <- list()

  }
  else {

    if (is_not_null(A[["mc.cores"]]) && A[["mc.cores"]] > 1) {
      if (isTRUE(A[["use.mbart2"]])) fun <- BART::mc.mbart2
      else fun <- BART::mc.mbart
    }
    else {
      if (isTRUE(A[["use.mbart2"]])) fun <- BART::mbart2
      else fun <- BART::mbart
    }

    fit.list <- do.call(fun, c(list(covs,
                                    y.train = as.integer(treat),
                                    x.test = covs,
                                    type = switch(A$link, "logit" = "lbart", "pbart")),
                               A[intersect(names(A), setdiff(names(formals(fun)), c("x.train", "y.train", "x.test", "type", "ntype")))]))

    ps[,] <- matrix(fit.list$prob.test.mean, nrow = length(treat), byrow = TRUE)

    info <- list()
    p.score <- NULL
  }

  #ps should be matrix of probs for each treat
  #Computing weights
  w <- get_w_from_ps(ps = ps, treat = treat, estimand, focal, stabilize = stabilize, subclass = subclass)

  obj <- list(w = w, ps = p.score, info = info, fit.obj = fit.list)
  return(obj)
}
.weightit2bart.cont <- function(covs, treat, s.weights, subset, stabilize, missing, ps, ...) {
  A <- list(...)

  check.package("BART")

  covs <- covs[subset, , drop = FALSE]
  treat <- treat[subset]
  s.weights <- s.weights[subset]

  if (!all_the_same(s.weights)) stop("Sampling weights cannot be used with method = \"bart\".",
                                     call. = FALSE)

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  #Process density params
  if (isTRUE(A[["use.kernel"]])) {
    if (is_null(A[["bw"]])) A[["bw"]] <- "nrd0"
    if (is_null(A[["adjust"]])) A[["adjust"]] <- 1
    if (is_null(A[["kernel"]])) A[["kernel"]] <- "gaussian"
    if (is_null(A[["n"]])) A[["n"]] <- 10*length(treat)
    use.kernel <- TRUE
    densfun <- NULL
  }
  else {
    if (is_null(A[["density"]])) densfun <- dnorm
    else if (is.function(A[["density"]])) densfun <- A[["density"]]
    else if (is.character(A[["density"]]) && length(A[["density"]] == 1)) {
      splitdens <- strsplit(A[["density"]], "_", fixed = TRUE)[[1]]
      if (exists(splitdens[1], mode = "function", envir = parent.frame())) {
        if (length(splitdens) > 1 && !can_str2num(splitdens[-1])) {
          stop(paste(A[["density"]], "is not an appropriate argument to density because",
                     word_list(splitdens[-1], and.or = "or", quotes = TRUE), "cannot be coerced to numeric."), call. = FALSE)
        }
        densfun <- function(x) {
          tryCatch(do.call(get(splitdens[1]), c(list(x), as.list(str2num(splitdens[-1])))),
                   error = function(e) stop(paste0("Error in applying density:\n  ", conditionMessage(e)), call. = FALSE))
        }
      }
      else {
        stop(paste(A[["density"]], "is not an appropriate argument to density because",
                   splitdens[1], "is not an available function."), call. = FALSE)
      }
    }
    else stop("The argument to density cannot be evaluated as a density function.", call. = FALSE)
    use.kernel <- FALSE
  }

  #Stabilization - get dens.num
  p.num <- treat - mean(treat)

  if (use.kernel) {
    d.n <- density(p.num, n = A[["n"]],
                   weights = s.weights/sum(s.weights), give.Rkern = FALSE,
                   bw = A[["bw"]], adjust = A[["adjust"]], kernel = A[["kernel"]])
    dens.num <- with(d.n, approxfun(x = x, y = y))(p.num)
  }
  else {
    dens.num <- densfun(p.num/sd(treat))
    if (is_null(dens.num) || !is.atomic(dens.num) || anyNA(dens.num)) {
      stop("There was a problem with the output of density. Try another density function or leave it blank to use the normal density.", call. = FALSE)
    }
    else if (any(dens.num <= 0)) {
      stop("The input to density may not accept the full range of treatment values.", call. = FALSE)
    }
  }

  if (is_not_null(A[["n.trees"]]) && is_null(A[["ntree"]])) A[["ntree"]] <- A[["n.trees"]]

  #Estimate GPS

  if (is_not_null(A[["mc.cores"]]) && A[["mc.cores"]] > 1) fun <- BART::mc.gbart
  else fun <- BART::gbart

  fit <- do.call(fun, c(list(covs,
                             y.train = treat,
                             type = "wbart"),
                        A[intersect(names(A), setdiff(names(formals(fun)),
                                                      c("x.train", "y.train", "x.test", "type", "ntype")))]))

  gp.score <- fit$yhat.train.mean

  #Get weights
  w <- get_cont_weights(gp.score, treat = treat, s.weights = s.weights,
                        dens.num = dens.num, densfun = densfun,
                        use.kernel = use.kernel, densControl = A)

  if (use.kernel && isTRUE(A[["plot"]])) {
    d.d <- density(treat - gp.score, n = A[["n"]],
                   weights = s.weights/sum(s.weights), give.Rkern = FALSE,
                   bw = A[["bw"]], adjust = A[["adjust"]],
                   kernel = A[["kernel"]])
    plot_density(d.n, d.d)
  }

  info <- list()

  obj <- list(w = w, info = info, fit.obj = fit)
  return(obj)
}

#ebal using ebal package (now using optim)----
.weightit2ebal <- function(covs, treat, s.weights, subset, estimand, focal, stabilize, missing, moments, int, ...) {
  A <- list(...)

  covs <- covs[subset, , drop = FALSE]
  treat <- factor(treat[subset])
  s.weights <- s.weights[subset]

  if (missing == "ind") {
    missing.ind <- apply(covs[, anyNA_col(covs), drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    if (is_not_null(missing.ind)) {
      covs[is.na(covs)] <- 0
      covs <- cbind(covs, missing.ind)
    }
  }

  covs <- cbind(covs, int.poly.f(covs, poly = moments, int = int, center = TRUE))
  for (i in seq_col(covs)) covs[,i] <- make.closer.to.1(covs[,i])

  for (f in names(formals(ebal::ebalance))) {
    if (is_null(A[[f]])) A[[f]] <- formals(ebal::ebalance)[[f]]
  }
  if (stabilize) {
    for (f in names(formals(ebal::ebalance.trim))) {
      if (is_null(A[[f]])) A[[f]] <- formals(ebal::ebalance.trim)[[f]]
    }
  }

  if (check.package("ebal")) {
    if (estimand == "ATT") {
      w <- rep(1, length(treat))
      control.levels <- levels(treat)[levels(treat) != focal]
      fit.list <- make_list(control.levels)

      covs[treat == focal,] <- covs[treat == focal, , drop = FALSE] * s.weights[treat == focal] * sum(treat == focal)/sum(s.weights[treat == focal])

      for (i in control.levels) {
        treat.in.i.focal <- treat %in% c(focal, i)
        treat_ <- as.integer(treat[treat.in.i.focal] != i)
        covs_ <- covs[treat.in.i.focal, , drop = FALSE]

        colinear.covs.to.remove <- colnames(covs_)[colnames(covs_) %nin% colnames(make_full_rank(covs_[treat_ == 0, , drop = FALSE]))]

        covs_ <- covs_[, colnames(covs_) %nin% colinear.covs.to.remove, drop = FALSE]

        if (is_not_null(A[["base.weight"]])) {
          A[["base.weight"]] <- A[["base.weight"]][treat == i]
        }

        ebal.out <- ebal::ebalance(Treatment = treat_, X = covs_,
                                   base.weight = A[["base.weight"]],
                                   norm.constant = A[["norm.constant"]],
                                   coefs = A[["coefs"]],
                                   max.iterations = A[["max.iterations"]],
                                   constraint.tolerance = A[["constraint.tolerance"]],
                                   print.level = 3)
        if (stabilize) ebal.out <- ebal::ebalance.trim(ebalanceobj = ebal.out,
                                                       max.weight = A[["max.weight"]],
                                                       min.weight = A[["min.weight"]],
                                                       max.trim.iterations = A[["max.trim.iterations"]],
                                                       max.weight.increment = A[["max.weight.increment"]],
                                                       min.weight.increment = A[["min.weight.increment"]],
                                                       print.level = 3)
        w[treat == i] <- ebal.out$w / s.weights[treat.in.i.focal][treat_ == 0]
        fit.list[[i]] <- ebal.out
      }
    }
    else if (estimand == "ATE") {
      w <- rep(1, length(treat))
      fit.list <- make_list(levels(treat))

      for (i in levels(treat)) {
        covs_i <- rbind(covs, covs[treat==i, , drop = FALSE])
        treat_i <- c(rep(1, nrow(covs)), rep(0, sum(treat==i)))

        colinear.covs.to.remove <- colnames(covs_i)[colnames(covs_i) %nin% colnames(make_full_rank(covs_i[treat_i == 0, , drop = FALSE]))]

        covs_i <- covs_i[, colnames(covs_i) %nin% colinear.covs.to.remove, drop = FALSE]

        covs_i[treat_i == 1,] <- covs_i[treat_i == 1,] * s.weights * sum(treat_i == 1) / sum(s.weights)

        if (is_not_null(A[["base.weight"]])) {
          A[["base.weight"]] <- A[["base.weight"]][treat == i]
        }

        ebal.out_i <- ebal::ebalance(Treatment = treat_i, X = covs_i,
                                     base.weight = A[["base.weight"]],
                                     norm.constant = A[["norm.constant"]],
                                     coefs = A[["coefs"]],
                                     max.iterations = A[["max.iterations"]],
                                     constraint.tolerance = A[["constraint.tolerance"]],
                                     print.level = 3)

        if (stabilize) ebal.out_i <- ebal::ebalance.trim(ebalanceobj = ebal.out_i,
                                                         max.weight = A[["max.weight"]],
                                                         min.weight = A[["min.weight"]],
                                                         max.trim.iterations = A[["max.trim.iterations"]],
                                                         max.weight.increment = A[["max.weight.increment"]],
                                                         min.weight.increment = A[["min.weight.increment"]],
                                                         print.level = 3)

        w[treat == i] <- ebal.out_i$w / s.weights[treat == i]
        fit.list[[i]] <- ebal.out_i
      }
    }
  }

  obj <- list(w = w, fit.obj = fit.list)
  return(obj)
}
