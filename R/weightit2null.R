weightit2null <- function(covs, treat, s.weights, subset, ...) {
  list(w = rep.int(1, length(treat[subset])))
}

weightit2null.multi <- weightit2null

weightit2null.cont <- weightit2null

weightitMSM2null <- function(covs.list, treat.list, s.weights, subset, ...) {
  list(w = rep.int(1, length(treat.list[[1L]][subset])))
}
