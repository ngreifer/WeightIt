weightitMSM <- function(formula.list, data, method = "ps", stabilize = FALSE, exact = NULL, s.weights = NULL,
                        verbose = FALSE, ...) {

  estimand <- "ATE"
  covs.list <- treat.list <- w.list <- ps.list <- vector("list", length(formula.list))
  for (i in seq_along(formula.list)) {
    #Process treat and covs from formula and data
    tt <- terms(formula.list[[i]])
    attr(tt, "intercept") <- 0
    if (is.na(match(rownames(attr(tt, "factors"))[1], names(data)))) {
      stop(paste0("The given response variable, \"", rownames(attr(tt, "factors"))[1], "\", is not a variable in data."))
    }
    vars.mentioned <- unlist(lapply(attr(tt, "variables")[-1], function(x) if (length(x) > 1) paste(x[-1]) else paste(x)))
    m.try <- try({mf <- model.frame(tt, data)}, TRUE)
    if (class(m.try) == "try-error") {
      stop(paste0(c("All variables of formula must be variables in data.\nVariables not in data: ",
                    paste(attr(tt, "term.labels")[is.na(match(attr(tt, "term.labels"), names(data)))], collapse=", "))), call. = FALSE)}
    treat.list[[i]] <- model.response(mf)
    covs.list[[i]] <- data[!is.na(match(names(data), vars.mentioned[vars.mentioned != rownames(attr(tt, "factors"))[1]]))]

    #Get weights into a list
    weightit_obj <- weightit(formula.list[[i]], data = data,
                             method = method,
                             estimand = estimand,
                             stabilize = stabilize,
                             exact = exact, s.weights = s.weights,
                             verbose = verbose, ...)
    w.list[[i]] <- cobalt::get.w(weightit_obj)
    ps.list[[i]] <- weightit_obj$ps
  }

  w <- Reduce("*", w.list)

  ## Assemble output object----
  out <- list(weights = w,
              treat.list = treat.list,
              covs.list = covs.list,
              data = data,
              estimand = NULL,
              method = method,
              ps.list = ps.list,
              s.weights = s.weights,
              #discarded = NULL,
              treat.type = NULL
              )
  class(out) <- c("weightitMSM", "weightit")

  return(out)
}
