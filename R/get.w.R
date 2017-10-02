get.w_ <- function(obj, ...) {
  UseMethod("get.w_", obj)
}
get.w_.default <- function(obj, estimand, ...) {
  t <- obj$t
  w <- obj$w
  ps <- obj$ps

  if (length(w) > 0) {

  }
  else if (length(ps) > 0 && length(t) > 0) {
    #May need to do some checks to ensure t is compatible with ps
    if (toupper(estimand) == "ATE") {
      w <- t/ps + (1-t)/(1-ps)
    }
    else if (toupper(estimand) == "ATT") {
      w <- t + (1-t)*ps/(1-ps)
    }
    else if (toupper(estimand) == "ATC") {
      w <- (1-t) + t*(1-ps)/ps
    }
    else if (toupper(estimand) == "ATO") {
      w <- t*(1-ps) + (1-t)*ps
    }
  }
  return(w)
}
