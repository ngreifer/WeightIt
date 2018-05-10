weightitG <- function(Sformula = NULL, data = NULL, targetdata = NULL, targetmargins  = NULL, method = "ps", exact = NULL,
                      verbose = FALSE, ...) {

  call <- match.call()
  s.data.names <- names(data)
  X <- sapply(c("Sformula", "data", "targetdata", "targetmargins"),
              function(x) length(get(x)) > 0)
  # X <- setNames(rep(FALSE,4), c("Sformula", "data", "targetdata", "targetmargins"))
  # for (i in names(X)) {
  #   X[i] <- length(get0(i)) > 0
  # }

  #selection formula + full data
  if (X["Sformula"] && X["data"] && !X["targetdata"] && !X["targetmargins"]) {
    formula <- Sformula
    full.data <- data
    allowable.methods = NULL #all
  }

  #sample data + target margins
  if (!X["Sformula"] && X["data"] && !X["targetdata"] && X["targetmargins"]) {
    t.data.names <- names(targetmargins)
    if (length(t.data.names) == 0) {
      t.data.names <- names(data)[seq_along(targetmargins)]
    }
    names.in.both <- intersect(s.data.names, t.data.names)
    if (length(names.in.both) == 0) {
      stop("no names in common")
    }

    full.data <- rbind(data.frame(.SAMPLE = rep(1, nrow(data)), data[names.in.both]),
                       setNames(data.frame(rep(0, 2), matrix(rep(targetmargins[names.in.both], 2), byrow = TRUE, nrow = 2)),
                                c(".SAMPLE", names.in.both)))

    formula <- formula(full.data)

    allowable.methods = c("ebal", "sbw", "ebcw")
    if (!method.to.proper.method(method) %in% allowable.methods) {
      stop(paste("Only methods", word.list(allowable.methods, quotes = TRUE, and.or = "and"),
                 "are allowed with targetmargins."), call. = FALSE)
    }
  }

  #sample data + target data (+ rhs formula)
  if (X["data"] && X["targetdata"] && !X["targetmargins"]) {
    t.data.names <- names(targetdata)
    names.in.both <- intersect(s.data.names, t.data.names)
    if (length(names.in.both) == 0) {
      stop("no names in common")
    }

    full.data <- rbind(data.frame(.SAMPLE = rep(1, nrow(data)), data[names.in.both]),
                       data.frame(.SAMPLE = rep(0, nrow(targetdata)), targetdata[names.in.both]))
    if (X["Sformula"] > 0) {
      formula <- update.formula(Sformula, .SAMPLE ~ .)
    }
    else {
      formula <- formula(full.data)
    }

    allowable.methods = NULL #all
  }

  #Check formula
  tt <- terms(formula)
  attr(tt, "intercept") <- 0
  sample.var <- all.vars(tt[[2]])
  if (is.na(match(sample.var, names(full.data)))) {
    stop(paste0("The given response variable, \"", all.vars(tt[[2]]), "\", is not a variable in data."))
  }
  vars.mentioned <- all.vars(tt)
  tryCatch({mf <- model.frame(tt, full.data, na.action = na.pass)}, error = function(e) {
    stop(paste0(c("All variables in Sformula must be variables in data.\nVariables not in data: ",
                  paste(vars.mentioned[is.na(match(vars.mentioned, names(full.data)))], collapse=", "))), call. = FALSE)})

  estimand = "ATC"

  W <- weightit(formula = formula, data = full.data,
                method = method, estimand = estimand,
                verbose = verbose, ...)
  w <- W$weights[W$treat == 1]
  p.score <- W$ps[W$treat == 1]

  ## Assemble output object----
  out <- list(g.weights = w,
              data = W$data[W$treat == 1, names(W$data) != sample.var],
              method = method,
              ps = if (length(p.score) == 0 || all(is.na(p.score))) NULL else p.score,
              call = call)

  class(out) <- c("weightitG", "weightit")

  return(out)
}
