weightit <- function(formula, data, method, estimand = "ATE", stabilize = FALSE, exact = NULL, s.weights = NULL, ...) {

  #Checks
  if (length(data) == 0) {
    stop("Data must be specified.", call. = FALSE)}
  if (!is.data.frame(data)) {
    stop("Data must be a data.frame.", call. = FALSE)}


  #Process method
  #methods: ps, gbm, CBPS, entropy balancing, stable balancing weights, ATE weights (using ATE package)
  bad.method <- FALSE
  acceptable.methods <- c("ps",
                          "gbm", "twang", "gbr",
                          "cbps",
                          "ebal", "entropy", "ebalance",
                          "sbw",
                          "ebcw", "ate")
  if (missing(method)) method <- "ps"
  else if (!is.character(method)) bad.method <- TRUE
  else if (length(method) != 1) bad.method <- TRUE
  else if (!tolower(method) %in% acceptable.methods) bad.method <- TRUE

  if (bad.method) stop("1")
  else method <- tolower(method)

  #Transform data (get treat and covs)
  ##!!!!!! Must change !!!!!!!!!

  #Initializing variables
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  if (is.na(match(rownames(attr(tt, "factors"))[1], names(data)))) {
    stop(paste0("The given response variable, \"", rownames(attr(tt, "factors"))[1], "\", is not a variable in data."))
  }
  m.try <- try({mf <- model.frame(tt, data)}, TRUE)
  if (class(m.try) == "try-error") {
    stop(paste0(c("All variables of formula must be variables in data.\nVariables not in data: ",
                  paste(attr(tt, "term.labels")[is.na(match(attr(tt, "term.labels"), names(data)))], collapse=", "))), call. = FALSE)}
  treat <- model.response(mf)
  covs <- data[, !is.na(match(names(data), attr(tt, "term.labels"))), drop = FALSE]

  n <- nrow(data)

  ##Process exact
  bad.exact <- FALSE
  acceptable.exacts <- names(data)
  if (missing(exact)) exact.factor <- factor(rep(1, n))
  else if (!is.atomic(exact)) bad.exact <- TRUE
  else if (is.character(exact) && all(exact %in% acceptable.exacts)) {
    exact.factor <- factor(apply(data[, exact, drop = FALSE], 1, paste, collapse = "|"))
  }
  else if (length(exact) == n) exact.factor <- factor(exact)
  else bad.exact <- TRUE

  if (bad.exact) stop("2")

  if (any(sapply(levels(exact), function(x) nunique(treat) != nunique(treat[exact == x])))) {
    stop("Not all the groups formed by exact contain all treatment levels. Consider reducing exact.", call. = FALSE)
  }

  ##Process s.weights
  if (length(s.weights) == 0) s.weights <- rep(1, nrow(data))

  w <- ps <- numeric(n)

  for (i in levels(exact.factor)) {
    #Run method
    if (method == "ps") {
      estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE", "ATO"), method)
      obj <- weightit2ps(formula = formula,
                    data = data,
                    s.weights = s.weights,
                    subset = exact.factor == i,
                    ...)
      ps[exact.factor == i] <- obj$ps
    }
    else if (method %in% c("gbm", "gbr", "twang")) {

    }
    else if (method %in% c("cbps")) {

      estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method)
      obj <- weightit2cbps(formula = formula,
                           data = data,
                           #s.weights = s.weights,
                           subset = exact.factor == i,
                           estimand = estimand,
                           #stabilize = stabilize,
                           ...)
    }
    else if (method %in% c("entropy", "ebal", "ebalance")) {

      estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method)
      obj <- weightit2ebal(formula = formula,
                      data = data,
                      #s.weights = s.weights,
                      subset = exact.factor == i,
                      estimand = estimand,
                      stabilize = stabilize,
                      ...)
    }
    else if (method %in% c("sbw")) {

    }
    else if (method %in% c("ebcw", "ate")) {
      estimand <- process.estimand(estimand, c("ATT", "ATC", "ATE"), method)
      obj <- weightit2ecbw(formula = formula,
                           data = data,
                           #s.weights = s.weights,
                           subset = exact.factor == i,
                           estimand = estimand,
                           #stabilize = stabilize,
                           ...)
    }

    #Extract weights (with get.w_)
    w[exact.factor == i] <- get.w_(obj, estimand = estimand)

  }

  #Assemble output object
  out <- list(weights = w,
             treat = treat,
             covs = covs,
             data = data,
             estimand = estimand,
             ps = ps,
             s.weights = s.weights,
             discarded = NULL)
  class(out) <- "weightit"

  return(out)
}
