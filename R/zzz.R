.onLoad <- function(libname, pkgname) {
  .on_cran <- !interactive() && !isTRUE(as.logical(Sys.getenv("NOT_CRAN", "false")))

  if (.on_cran) {
    # CRAN OMP THREAD LIMIT
    Sys.setenv("OMP_THREAD_LIMIT" = 1L)
  }

}