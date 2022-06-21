#Criterion
#Criteria to use in methods that involve a tuning parameter. Also called "stop.method" in twang.
#Separate functions for different treatment types.

#init functions
init_es <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, pairwise = TRUE, ...) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")

    bin.vars <- process.bin.vars(mat = covs)

    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call("splitfactor", list(covs, drop.first ="if2",
                                            split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }

    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector.")

    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }

    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))

    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)

    if (treat.type == "continuous") stop("treat must be a binary or multinomial variable.")

    f.e.r <- process.focal.and.estimand(focal, estimand, treat)
    focal <- f.e.r[["focal"]]
    estimand <- f.e.r[["estimand"]]

    unique.treats <- unique(treat)

    if (treat.type == "multinomial") {
        if (is_null(focal) && !pairwise) {
            treat.all <- last(make.unique(unique.treats, "All"))
            treat <- factor(c(as.character(treat), rep(treat.all, length(treat))),
                            levels = c(unique.treats, treat.all))
            covs <- rbind(covs, covs)
            s.weights <- rep(s.weights, 2)
            focal <- treat.all
        }
        else pairwise <- TRUE

        if (is_null(focal)) {
            treatment.pairs <- combn(unique.treats, 2, simplify = FALSE)
        }
        else {
            treatment.pairs <- lapply(setdiff(unique.treats, focal), function(x) c(x, focal))
        }
    }
    else {
        treatment.pairs <- list(unique.treats)
        pairwise <- TRUE
    }

    s.d.denom <- get.s.d.denom.weightit(estimand = estimand, treat = treat, focal = focal)

    denoms <- compute_s.d.denom(covs, treat = treat,
                                s.d.denom = s.d.denom, s.weights = s.weights,
                                bin.vars = bin.vars, to.sd = rep(TRUE, ncol(covs)))

    out <- list(treat = treat,
                covs = covs,
                bin.vars = bin.vars,
                s.weights = s.weights,
                s.d.denom = denoms,
                focal = focal,
                pairwise = pairwise,
                treatment.pairs = treatment.pairs)
    class(out) <- "init_es"
    out
}
init_ks <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, pairwise = TRUE, ...) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")

    bin.vars <- process.bin.vars(mat = covs)

    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call("splitfactor", list(covs, drop.first ="if2",
                                            split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }

    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector.")

    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }

    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))

    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)

    if (treat.type == "continuous") stop("treat must be a binary or multinomial variable.")

    f.e.r <- process.focal.and.estimand(focal, estimand, treat)
    focal <- f.e.r[["focal"]]
    estimand <- f.e.r[["estimand"]]

    unique.treats <- unique(treat)

    if (treat.type == "multinomial") {
        if (is_null(focal) && !pairwise) {
            treat.all <- last(make.unique(unique.treats, "All"))
            treat <- factor(c(as.character(treat), rep(treat.all, length(treat))),
                            levels = c(unique.treats, treat.all))
            covs <- rbind(covs, covs)
            s.weights <- rep(s.weights, 2)
            focal <- treat.all
        }
        else pairwise <- TRUE

        if (is_null(focal)) {
            treatment.pairs <- combn(unique.treats, 2, simplify = FALSE)
        }
        else {
            treatment.pairs <- lapply(setdiff(unique.treats, focal), function(x) c(x, focal))
        }
    }
    else {
        treatment.pairs <- list(unique.treats)
        pairwise <- TRUE
    }

    out <- list(treat = treat,
                covs = covs,
                bin.vars = bin.vars,
                s.weights = s.weights,
                focal = focal,
                pairwise = pairwise,
                treatment.pairs = treatment.pairs)
    class(out) <- "init_ks"
    out
}
init_mahalanobis <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, ...) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")

    bin.vars <- process.bin.vars(mat = covs)

    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call("splitfactor", list(covs, drop.first = TRUE,
                                            split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }

    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector.")

    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }

    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))

    if (!is_binary(treat)) stop("treat must be a binary variable.")

    s.d.denom <- get.s.d.denom.weightit(estimand = estimand, treat = treat, focal = focal)

    if (any(!bin.vars)) covs[,!bin.vars] <- scale(covs[,!bin.vars])

    if (s.d.denom %in% as.character(treat)) {
        sigma <- cov.wt(covs[treat == s.d.denom,,drop = FALSE], s.weights[treat == s.d.denom])$cov
        if (any(zeros <- diag(sigma) == 0)) {
            sigma_all <- cov.wt(covs, s.weights)$cov
            sigma[zeros,] <- sigma_all[zeros,]
            sigma[,zeros] <- sigma_all[,zeros]
            if (any(zeros <- diag(sigma) == 0)) diag(sigma)[zeros]  <- 1
        }
    }
    else if (s.d.denom == "pooled") {
        sigma <- .5 * (cov.wt(covs[treat == treat[1],,drop = FALSE], s.weights[treat == treat[1]])$cov +
                           cov.wt(covs[treat != treat[1],,drop = FALSE], s.weights[treat != treat[1]])$cov)
        if (any(zeros <- diag(sigma) == 0)) diag(sigma)[zeros] <- 1
    }
    else {
        sigma <- cov.wt(covs, s.weights)$cov
        if (any(zeros <- diag(sigma) == 0)) diag(sigma)[zeros] <- 1
    }

    #MASS::ginv
    sigmasvd <- svd(sigma)
    pos <- sigmasvd$d > max(1e-8 * sigmasvd$d[1L], 0)
    sigma_inv <- sigmasvd$v[, pos, drop = FALSE] %*% ((1/sigmasvd$d[pos]) *
                                                          t(sigmasvd$u[, pos, drop = FALSE]))

    out <- list(treat = treat,
                covs = covs,
                bin.vars = bin.vars,
                s.weights = s.weights,
                s.d.denom = s.d.denom,
                sigma_inv = sigma_inv)
    class(out) <- "init_mahalanobis"
    out
}
init_energy.dist <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, improved = TRUE, ...) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")

    if (needs.splitting) {
        covs <- do.call("splitfactor", list(covs, drop.first ="if2"))
    }

    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector.")

    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }

    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))

    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)

    if (treat.type == "continuous") stop("treat must be a binary or multinomial variable.")

    f.e.r <- process.focal.and.estimand(focal, estimand, treat)
    focal <- f.e.r[["focal"]]
    estimand <- f.e.r[["estimand"]]

    dist.covs <- transform_covariates(data = covs, method = "scaled_euclidean",
                                      s.weights = s.weights)
    d <- unname(eucdist_internal(dist.covs))

    n <- length(treat)
    unique.treats <- unique(treat)
    treat.seq <- seq_along(unique.treats)
    diagn <- diag(n)

    for (t in unique.treats) s.weights[treat == t] <- s.weights[treat == t]/mean_fast(s.weights[treat == t])

    tmat <- vapply(unique.treats, function(t) treat == t, logical(n))
    nt <- colSums(tmat)

    J <- setNames(lapply(treat.seq, function(t) s.weights*tmat[,t]/nt[t]), unique.treats)

    if (is_null(focal)) {
        J0 <- as.matrix(s.weights/n)

        M2_array <- vapply(treat.seq, function(t) -2 * tcrossprod(J[[t]]) * d, diagn)
        M1_array <- vapply(treat.seq, function(t) 2 * J[[t]] * d %*% J0, J0)

        M2 <- rowSums(M2_array, dims = 2)
        M1 <- rowSums(M1_array)

        if (improved) {
            all_pairs <- combn(treat.seq, 2, simplify = FALSE)
            M2_pairs_array <- vapply(all_pairs, function(p) -2 * tcrossprod(J[[p[1]]]-J[[p[2]]]) * d, diagn)
            M2 <- M2 + rowSums(M2_pairs_array, dims = 2)
        }
    }
    else {
        J0_focal <- as.matrix(J[[focal]])
        clevs <- treat.seq[unique.treats != focal]

        M2_array <- vapply(clevs, function(t) -2 * tcrossprod(J[[t]]) * d, diagn)
        M1_array <- vapply(clevs, function(t) 2 * J[[t]] * d %*% J0_focal, J0_focal)

        M2 <- rowSums(M2_array, dims = 2)
        M1 <- rowSums(M1_array)
    }

    out <- list(M1 = M1,
                M2 = M2,
                s.weights = s.weights,
                treat = treat,
                unique.treats = unique.treats)
    class(out) <- "init_energy.dist"
    out
}
init_p <- function(covs, treat, s.weights = NULL, ...) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")

    bin.vars <- process.bin.vars(mat = covs)

    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call("splitfactor", list(covs, drop.first ="if2",
                                            split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }

    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector.")

    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }

    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))

    if (get.treat.type(assign.treat.type(treat)) != "continuous") stop("treat must be a numeric non-binary variable.")

    s.d.denom <- get.s.d.denom.cont.weightit()

    denoms <- compute_s.d.denom(covs, treat = treat,
                                s.d.denom = s.d.denom, s.weights = s.weights,
                                bin.vars = bin.vars, to.sd = rep(TRUE, ncol(covs)))

    out <- list(treat = treat,
                covs = covs,
                bin.vars = bin.vars,
                s.weights = s.weights,
                s.d.denom = denoms)
    class(out) <- "init_p"
    out
}
init_s <- function(covs, treat, s.weights = NULL, ...) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")

    bin.vars <- process.bin.vars(mat = covs)

    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call("splitfactor", list(covs, drop.first ="if2",
                                            split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }

    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector.")

    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }

    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))

    if (get.treat.type(assign.treat.type(treat)) != "continuous") stop("treat must be a numeric non-binary variable.")

    for (i in seq_col(covs)) if (!bin.vars[i]) covs[,i] <- rank(covs[,i], na.last = "keep")
    treat <- rank(treat, na.last = "keep")

    s.d.denom <- get.s.d.denom.cont.weightit()

    denoms <- compute_s.d.denom(covs, treat = treat,
                                s.d.denom = s.d.denom, s.weights = s.weights,
                                bin.vars = bin.vars, to.sd = rep(TRUE, ncol(covs)))

    out <- list(treat = treat,
                covs = covs,
                bin.vars = bin.vars,
                s.weights = s.weights,
                s.d.denom = denoms)
    class(out) <- "init_s"
    out
}
init_r2 <- function(covs, treat, s.weights = NULL, ...) {
    needs.splitting <- FALSE
    if (!is.matrix(covs)) {
        if (is.data.frame(covs)) {
            if (any(to.split <- vapply(covs, is_, logical(1L), types = c("factor", "character")))) {
                needs.splitting <- TRUE
            }
            else covs <- as.matrix(covs)
        }
        else if (is.numeric(covs)) covs <- matrix(covs, ncol = 1)
        else stop("covs must be a data.frame or numeric matrix.")
    }
    else if (!is.numeric(covs)) stop("covs must be a data.frame or numeric matrix.")

    bin.vars <- process.bin.vars(mat = covs)

    if (needs.splitting) {
        bin.vars[to.split] <- TRUE
        covs <- do.call("splitfactor", list(covs, drop.first ="if2",
                                            split.with = bin.vars))
        bin.vars <- attr(covs, "split.with")[[1]]
    }

    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector.")

    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }

    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))

    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)

    if (treat.type == "multinomial") stop("treat must be a binary or continuous variable.")

    if (treat.type == "binary") treat <- as.numeric(treat == treat[1])

    x <- cbind(covs)

    out <- list(treat = treat,
                x = x,
                s.weights = s.weights)
    class(out) <- "init_r2"
    out
}
init_L1.med <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, .covs = NULL, L1.min.bin = 2, L1.max.bin = 12, L1.n = 101, ...) {

    if (is_not_null(.covs)) covs <- .covs
    if (!is.data.frame(covs)) {
        if (is.atomic(covs) && is_null(dim(covs))) covs <- data.frame(covs)
        else if (!is.matrix(covs)) stop("covs must be a data.frame or matrix.")
    }
    covs <- as.data.frame(covs)

    if (missing(treat) || !is.atomic(treat)) stop("treat must be an atomic vector.")

    possibly.supplied <- c("covs", "treat", "s.weights")
    lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                        possibly.supplied)
    supplied <- lengths > 0
    if (!all_the_same(lengths[supplied])) {
        stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
    }

    if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(covs))

    if (!has.treat.type(treat)) treat <- assign.treat.type(treat)
    treat.type <- get.treat.type(treat)

    coarsen <- function(covs, cutpoints = NULL, grouping = NULL) {
        is.numeric.cov <- setNames(vapply(covs, is.numeric, logical(1L)), names(covs))
        for (i in names(cutpoints)) {
            if (cutpoints[[i]] == 0) is.numeric.cov[i] <- FALSE #Will not be binned
        }

        #Process grouping
        if (!is.null(grouping) && !is.null(names(grouping))) {
            covs[names(grouping)] <- lapply(names(grouping), function(g) {
                x <- covs[[g]]
                groups <- grouping[[g]]

                for (i in seq_along(groups)) {
                    x[x %in% groups[[i]]] <- groups[[i]][1]
                }
                x
            })
            cutpoints[names(cutpoints) %in% names(grouping)] <- NULL
        }

        #Create bins for numeric variables
        for (i in names(covs)[is.numeric.cov]) {
            bins <- cutpoints[[i]]

            #cutpoints is number of bins, unlike in cem
            breaks <- seq(min(covs[[i]]), max(covs[[i]]), length = bins + 1)
            breaks[c(1, bins + 1)] <- c(-Inf, Inf)

            covs[[i]] <- findInterval(covs[[i]], breaks)
        }

        #Reduce to strata
        factor(do.call("paste", c(covs, sep = " | ")))
    }

    f.e.r <- process.focal.and.estimand(focal, estimand, treat)
    focal <- f.e.r[["focal"]]
    estimand <- f.e.r[["estimand"]]

    if (treat.type == "continuous") {
        unique.treats <- focal <- NULL
        nunique.treat <- nunique(treat)
    }
    else {
        unique.treats <- unique(treat)
        for (t in unique.treats) s.weights[treat == t] <- s.weights[treat == t]/sum(s.weights[treat == t])
    }

    is.numeric.cov <- setNames(vapply(covs, is.numeric, logical(1L)), names(covs))
    nunique.covs <- vapply(covs, nunique, integer(1L))

    if (treat.type == "continuous") {
        coarsenings <- lapply(1:L1.n, function(i) {
            cutpoints <- setNames(lapply(nunique.covs[is.numeric.cov], function(nu) {
                sample(seq(min(L1.min.bin, nu), min(L1.max.bin, nu)), 1)
            }), names(covs)[is.numeric.cov])
            grouping <- setNames(lapply(seq_along(covs)[!is.numeric.cov], function(x) {
                nu <- nunique.covs[x]
                u <- unique(covs[[x]], nmax = nu)

                #Randomly select number of bins
                nbins <- sample(seq(min(L1.min.bin, nu), min(L1.max.bin, nu)), 1)

                #Randomly assign bin numbers to levels of covariate
                bin.assignments <- sample(seq_len(nbins), nu, replace = TRUE)

                #Group levels with same bin number
                lapply(unique(bin.assignments, nmax = nbins),
                       function(b) u[bin.assignments == b])

            }), names(covs)[!is.numeric.cov])
            treat_cutpoints <- list(.treat = sample(seq(min(L1.min.bin, nunique.treat), min(L1.max.bin, nunique.treat)), 1))
            list(cutpoints = cutpoints, grouping = grouping, treat_cutpoints = treat_cutpoints)
        })
    }
    else {
        coarsenings <- lapply(1:L1.n, function(i) {
            cutpoints <- setNames(lapply(nunique.covs[is.numeric.cov], function(nu) {
                sample(seq(min(L1.min.bin, nu), min(L1.max.bin, nu)), 1)
            }), names(covs)[is.numeric.cov])
            grouping <- setNames(lapply(seq_along(covs)[!is.numeric.cov], function(x) {
                nu <- nunique.covs[x]
                u <- unique(covs[[x]], nmax = nu)

                #Randomly select number of bins
                nbins <- sample(seq(min(L1.min.bin, nu), min(L1.max.bin, nu)), 1)

                #Randomly assign bin numbers to levels of covariate
                bin.assignments <- sample(seq_len(nbins), nu, replace = TRUE)

                #Group levels with same bin number
                lapply(unique(bin.assignments, nmax = nbins),
                       function(b) u[bin.assignments == b])

            }), names(covs)[!is.numeric.cov])
            list(cutpoints = cutpoints, grouping = grouping, treat_cutpoints = NULL)
        })
    }

    L1s <- unlist(lapply(coarsenings, function(co) {
        x <- coarsen(covs, cutpoints = co[["cutpoints"]], grouping = co[["grouping"]])

        if (treat.type == "continuous") {
            treat <- coarsen(data.frame(.treat = treat), cutpoints = co[["treat_cutpoints"]])
            unique.treats <- unique(treat, nmax = L1.max.bin)
            for (t in unique.treats) s.weights[treat == t] <- s.weights[treat == t]/sum(s.weights[treat == t])
        }

        if (treat.type == "binary" || is_null(focal)) {
            sum(vapply(levels(x), function(l) {
                in_l <- which(x == l)
                abs(diff(range(vapply(unique.treats, function(t) sum(s.weights[in_l][treat[in_l] == t]), numeric(1L)))))
            }, numeric(1L))) / length(unique.treats)
        }
        else {
            sum(vapply(levels(x), function(l) {
                in_l <- which(x == l)
                sum.s.weights.focal <- sum(s.weights[in_l][treat[in_l] == focal])
                max(abs(vapply(unique.treats[unique.treats != focal],
                               function(t) sum(s.weights[in_l][treat[in_l] == t]) - sum.s.weights.focal, numeric(1L))))
            }, numeric(1L))) / length(unique.treats)
        }
    }))

    L1.med <- sort(L1s, partial = ceiling(L1.n/2))[ceiling(L1.n/2)]

    L1.med.coarsening <- coarsenings[[which(L1s == L1.med)[1]]]

    if (treat.type == "continuous") {
        treat <- coarsen(data.frame(.treat = treat), cutpoints = L1.med.coarsening[["treat_cutpoints"]])
        unique.treats <- unique(treat, nmax = L1.max.bin)
    }

    out <- list(coarsened.covs = coarsen(covs,
                                         cutpoints = L1.med.coarsening[["cutpoints"]],
                                         grouping = L1.med.coarsening[["grouping"]]),
                s.weights = s.weights,
                treat = treat,
                unique.treats = unique.treats,
                focal = focal)
    class(out) <- "init_L1.med"
    out

}

#Statistics
es.binary <- function(init, weights = NULL) {
    check_init(init, "init_es")
    cobalt::col_w_smd(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
                      bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE)
}
ks.binary <- function(init, weights = NULL) {
    check_init(init, "init_ks")
    cobalt::col_w_ks(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
                     bin.vars = init$bin.vars)
}
mahalanobis.binary <- function(init, weights = NULL) {
    check_init(init, "init_mahalanobis")
    mean.diffs <- cobalt::col_w_smd(init$covs, init$treat, weights, s.weights = init$s.weights,
                                    bin.vars = init$bin.vars, std = FALSE)
    drop(sqrt(t(mean.diffs) %*% init$sigma_inv %*% mean.diffs))
}
energy.dist.binary <- function(init, weights = NULL) {
    check_init(init, "init_energy.dist")

    if (is_null(weights)) weights <- rep(1, nrow(init[["M2"]]))

    weights <- weights * init[["s.weights"]]

    for (t in init[["unique.treats"]]) weights[init[["treat"]] == t] <- weights[init[["treat"]] == t]/mean_fast(weights[init[["treat"]] == t])

    return(drop(.5 * t(weights) %*% init[["M2"]] %*% weights + t(init[["M1"]]) %*% weights))
}
r2.binary <- function(init, weights = NULL) {
    check_init(init, "init_r2")
    if (is_null(weights)) weights <- rep(1, length(init$treat))

    weights <- weights * init$s.weights

    fit <- glm.fit(init$x, init$treat, weights, family = quasibinomial())

    wmtreat <- sum(weights*fit$linear.predictors)/sum(weights)

    SSmodel <- sum(weights * (fit$linear.predictors - wmtreat)^2)

    SSmodel / (sum(weights) * pi^2/3 + SSmodel)
}
L1.med.binary <- function(init, weights = NULL) {
    check_init(init, "init_L1.med")

    if (is_null(weights)) weights <- rep(1, length(init[["treat"]]))

    weights <- weights * init[["s.weights"]]

    for (t in init[["unique.treats"]]) weights[init[["treat"]] == t] <- weights[init[["treat"]] == t]/sum(weights[init[["treat"]] == t])

    x <- init[["coarsened.covs"]]

    sum(vapply(levels(x), function(l) {
        in_l <- which(x == l)
        abs(diff(vapply(init[["unique.treats"]], function(t) sum(weights[in_l][init[["treat"]][in_l] == t]), numeric(1L))))
    }, numeric(1L))) / 2

}
es.multinomial <- function(init, weights = NULL) {
    check_init(init, "init_es")

    if (!init$pairwise) {
        weights <- c(weights, rep(1, length(weights)))
    }

    unlist(lapply(init$treatment.pairs, function(x) {
        cobalt::col_w_smd(init$covs[init$treat %in% x,,drop = FALSE],
                          treat = init$treat[init$treat %in% x],
                          weights = weights[init$treat %in% x],
                          s.weights = init$s.weights[init$treat %in% x],
                          bin.vars = init$bin.vars,
                          s.d.denom = init$s.d.denom, abs = TRUE)
    }))
}
ks.multinomial <- function(init, weights = NULL) {
    check_init(init, "init_ks")

    if (!init$pairwise) {
        weights <- c(weights, rep(1, length(weights)))
    }

    unlist(lapply(init$treatment.pairs, function(x) {
        cobalt::col_w_ks(init$covs[init$treat %in% x,,drop = FALSE],
                         treat = init$treat[init$treat %in% x],
                         weights = weights[init$treat %in% x],
                         s.weights = init$s.weights[init$treat %in% x],
                         bin.vars = init$bin.vars)
    }))
}
energy.dist.multinomial <- energy.dist.binary
L1.med.multinomial <- function(init, weights = NULL) {
    check_init(init, "init_L1.med")

    if (is_null(weights)) weights <- rep(1, length(init[["treat"]]))

    weights <- weights * init[["s.weights"]]

    for (t in init[["unique.treats"]]) weights[init[["treat"]] == t] <- weights[init[["treat"]] == t]/sum(weights[init[["treat"]] == t])

    x <- init[["coarsened.covs"]]

    if (is_null(init[["focal"]])) {
        sum(vapply(levels(x), function(l) {
            in_l <- which(x == l)
            abs(diff(range(vapply(init[["unique.treats"]], function(t) sum(weights[in_l][init[["treat"]][in_l] == t]), numeric(1L)))))
        }, numeric(1L))) / length(init[["unique.treats"]])
    }
    else {
        sum(vapply(levels(x), function(l) {
            in_l <- which(x == l)
            sum.weights.focal <- sum(weights[in_l][init[["treat"]][in_l] == init[["focal"]]])
            max(abs(vapply(init[["unique.treats"]][init[["unique.treats"]] != init[["focal"]]],
                           function(t) sum(weights[in_l][init[["treat"]][in_l] == t]) - sum.weights.focal, numeric(1L))))
        }, numeric(1L))) / length(init[["unique.treats"]])
    }

}
pearson.corr.continuous <- function(init, weights = NULL) {
    check_init(init, "init_p")
    cobalt::col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
                      bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
                      std = TRUE)
}
spearman.corr.continuous <- function(init, weights = NULL) {
    check_init(init, "init_s")
    cobalt::col_w_cov(init$covs, treat = init$treat, weights = weights, s.weights = init$s.weights,
                      bin.vars = init$bin.vars, s.d.denom = init$s.d.denom, abs = TRUE,
                      std = TRUE)
}
r2.continuous <- function(init, weights = NULL) {
    check_init(init, "init_r2")

    if (is_null(weights)) weights <- rep(1, length(init$treat))

    weights <- weights * init$s.weights

    fit <- lm.wfit(init$x, init$treat, weights)
    wmtreat <- sum(weights*init$treat)/sum(weights)
    SSmodel <- sum(weights * (fit$fitted.values - wmtreat)^2)
    SStreat <- sum(weights * (init$treat - wmtreat)^2)
    SSmodel / SStreat
}
L1.med.continuous <- function(init, weights = NULL) {
    check_init(init, "init_L1.med")

    if (is_null(weights)) weights <- rep(1, length(init[["treat"]]))

    weights <- weights * init[["s.weights"]]

    for (t in init[["unique.treats"]]) weights[init[["treat"]] == t] <- weights[init[["treat"]] == t]/sum(weights[init[["treat"]] == t])

    x <- init[["coarsened.covs"]]

    sum(vapply(levels(x), function(l) {
        in_l <- which(x == l)
        abs(diff(range(vapply(init[["unique.treats"]], function(t) sum(weights[in_l][init[["treat"]][in_l] == t]), numeric(1L)))))
    }, numeric(1L))) / 2

    sum(vapply(levels(x), function(l) {
        in_l <- which(x == l)
        abs(diff(range(vapply(init[["unique.treats"]], function(t) sum(weights[in_l][init[["treat"]][in_l] == t]), numeric(1L)))))
    }, numeric(1L))) / length(init[["unique.treats"]])

}

bal_criterion <- function(treat.type = "binary", criterion, list = FALSE) {
    treat.type <- match_arg(treat.type, c("binary", "multinomial", "continuous"))
    criteria <- switch(
        treat.type,
        binary = list(
            es.mean = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_es(covs, treat, s.weights, estimand)
                    }
                    mean_fast(es.binary(init, weights))
                },
                init = init_es
            ),
            es.max = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_es(covs, treat, s.weights, estimand)
                    }
                    max(es.binary(init, weights))
                },
                init = init_es
            ),
            es.rms = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_es(covs, treat, s.weights, estimand)
                    }
                    rms(es.binary(init, weights))
                },
                init = init_es
            ),
            ks.mean = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_ks(covs, treat, s.weights)
                    }
                    mean_fast(ks.binary(init, weights))
                },
                init = init_ks
            ),
            ks.max = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_ks(covs, treat, s.weights)
                    }
                    max(ks.binary(init, weights))
                },
                init = init_ks
            ),
            ks.rms = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_ks(covs, treat, s.weights)
                    }
                    rms(ks.binary(init, weights))
                },
                init = init_ks
            ),
            mahalanobis = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_mahalanobis(covs, treat, s.weights, estimand)
                    }
                    mahalanobis.binary(init, weights)
                },
                init = init_mahalanobis
            ),
            energy.dist = list(
                fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, improved = TRUE, init = NULL, ...) {

                    if (is_null(init)) {
                        init <- init_energy.dist(covs, treat, s.weights, estimand, focal, improved)
                    }
                    energy.dist.binary(init, weights)
                },
                init = init_energy.dist
            ),
            # kernel.dist = list(
            #     fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, init = NULL, ...) {
            #
            #         if (is_null(init)) {
            #             init <- init_kernel.dist(covs, treat, s.weights, estimand, focal)
            #         }
            #         kernel.dist.binary(init, weights)
            #     },
            #     init = init_kernel.dist
            # ),
            L1.med = list(
                fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_L1.med(covs, treat, s.weights, estimand, focal, ...)
                    }
                    L1.med.binary(init, weights)
                },
                init = init_L1.med
            ),
            r2 = list(
                fun = function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_r2(covs, treat, s.weights)
                    }
                    r2.binary(init, weights)
                },
                init = init_r2
            )
        ),
        multinomial = list(
            es.mean = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", focal = NULL, pairwise = TRUE, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_es(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                    }
                    mean_fast(es.multinomial(init, weights))
                },
                init = init_es
            ),
            es.max = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", focal = NULL, pairwise = TRUE, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_es(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                    }
                    max(es.multinomial(init, weights))
                },
                init = init_es
            ),
            es.rms = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, estimand = "ATE", focal = NULL, pairwise = TRUE, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_es(covs, treat, s.weights, estimand = estimand, focal = focal, pairwise = pairwise)
                    }
                    rms(es.multinomial(init, weights))
                },
                init = init_es
            ),
            ks.mean = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL, pairwise = TRUE, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_ks(covs, treat, s.weights, focal = focal, pairwise = pairwise)
                    }
                    mean_fast(ks.multinomial(init, weights))
                },
                init = init_ks
            ),
            ks.max = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL, pairwise = TRUE, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_ks(covs, treat, s.weights, focal = focal, pairwise = pairwise)
                    }
                    max(ks.multinomial(init, weights))
                },
                init = init_ks
            ),
            ks.rms = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, focal = NULL, pairwise = TRUE, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_ks(covs, treat, s.weights, focal = focal, pairwise = pairwise)
                    }
                    rms(ks.multinomial(init, weights))
                },
                init = init_ks
            ),
            energy.dist = list(
                fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, improved = TRUE, init = NULL, ...) {

                    if (is_null(init)) {
                        init <- init_energy.dist(covs, treat, s.weights, estimand, focal, improved)
                    }
                    energy.dist.multinomial(init, weights)
                },
                init = init_energy.dist
            ),
            L1.med = list(
                fun = function(covs, treat, weights = NULL, s.weights = NULL, estimand = "ATE", focal = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_L1.med(covs, treat, s.weights, estimand, focal, ...)
                    }
                    L1.med.multinomial(init, weights)
                },
                init = init_L1.med
            )
        ),

        continuous = list(
            p.mean = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_p(covs, treat, s.weights)
                    }
                    mean_fast(pearson.corr.continuous(init, weights))
                },
                init = init_p
            ),
            p.max = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_p(covs, treat, s.weights)
                    }
                    max(pearson.corr.continuous(init, weights))
                },
                init = init_p
            ),
            p.rms = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_p(covs, treat, s.weights)
                    }
                    rms(pearson.corr.continuous(init, weights))
                },
                init = init_p
            ),
            s.mean = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_s(covs, treat, s.weights)
                    }
                    mean_fast(spearman.corr.continuous(init, weights))
                },
                init = init_s
            ),
            s.max = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_s(covs, treat, s.weights)
                    }
                    max(spearman.corr.continuous(init, weights))
                },
                init = init_s
            ),
            s.rms = list(
                fun = function(covs, treat, weights = NULL, bin.vars, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_s(covs, treat, s.weights)
                    }
                    rms(spearman.corr.continuous(init, weights))
                },
                init = init_s
            ),
            L1.med = list(
                fun = function(covs, treat, weights = NULL, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_L1.med(covs, treat, s.weights, ...)
                    }
                    L1.med.continuous(init, weights)
                },
                init = init_L1.med
            ),
            r2 = list(
                fun <- function(covs, treat, weights, s.weights = NULL, init = NULL, ...) {
                    if (is_null(init)) {
                        init <- init_r2(covs, treat, s.weights)
                    }
                    r2.continuous(init, weights)
                },
                init = init_r2
            )
        )
    )
    if (isTRUE(list)) return(names(criteria))
    if (missing(criterion)) stop("'criterion' must be provided.", call. = FALSE)
    criterion <- match_arg(criterion, names(criteria))

    out <- criteria[[criterion]]
    attr(out, "treat.type") <- treat.type
    attr(out, "criterion") <- criterion
    class(out) <- "bal_criterion"
    out
}

print.bal_criterion <- function(x, ...) {
    cat("A bal_criterion object\n")
    cat(paste0("  treatment type: ", attr(x, "treat.type"), "\n"))
    cat(paste0("  criterion: ", attr(x, "criterion"), " (", bal_criterion.to.phrase(attr(x, "criterion")), ")\n"))
    cat("With components:\n - $init (for initialization)\n - $fun (to compute the criterion value)\n")
    invisible(x)
}
bal_criterion.to.phrase <- function(criterion) {

    phrase <- switch(criterion,
                     "es.mean" = "average absolute standardized mean difference",
                     "es.max" = "maximum absolute standardized mean difference",
                     "es.rms" = "root-mean-square absolute standardized mean difference",
                     "ks.mean" = "average Kolmogorov-Smirnov statistic",
                     "ks.max" = "maximum Kolmogorov-Smirnov statistic",
                     "ks.rms" = "root-mean-square Kolmogorov-Smirnov statistic",
                     "mahalanobis" = "sample Mahalanobis distance",
                     "energy.dist" = "energy distance",
                     # "kernel.dist" = "kernel distance",
                     "L1.med" = "L1 median",
                     "r2" = "post-weighting treatment R-squared",
                     NA_character_
    )
    if (anyNA(phrase)) stop(paste0("\"", criterion, "\" is not an allowed criterion."))
    return(phrase)
}

check_init <- function(init, init_class) {
    if (missing(init)) stop("'init' must be specified.")
    if (!inherits(init, init_class)) stop(paste0("'init' must be of class \"", init_class, "\"."))
}