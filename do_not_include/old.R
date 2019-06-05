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
    treat <- factor(treat)[subset]

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
