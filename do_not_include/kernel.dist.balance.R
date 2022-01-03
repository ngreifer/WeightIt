init_kernel.dist <- function(covs, treat, s.weights = NULL, estimand = "ATE", focal = NULL, improved = TRUE, ...) {
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

  covs <- mat_div(center(covs, at = col.w.m(covs, s.weights)),
                  sqrt(col.w.v(covs, s.weights)))

  d <- as.matrix(dist(covs))

  unique.treats <- unique(treat)

  out <- list(expd = exp(-d),
              s.weights = s.weights,
              treat = treat,
              unique.treats = unique.treats)
  class(out) <- "init_kernel.dist"
  out
}

kernel.dist.binary <- function(init, weights = NULL) {
  check_init(init, "init_kernel.dist")

  if (is_null(weights)) weights <- rep(1, length(init[["treat"]]))

  weights <- weights * init[["s.weights"]]

  for (t in init[["unique.treats"]]) weights[init[["treat"]] == t] <- weights[init[["treat"]] == t]/sum(weights[init[["treat"]] == t])

  Tstar <- weights * (2 * (init[["treat"]] == init[["unique.treats"]][1]) - 1)
  sqrt(sum(outer(Tstar, Tstar) * init[["expd"]]))
}
