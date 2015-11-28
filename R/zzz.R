.onAttach <- function(...) {
  ver <- utils::packageVersion("loo")
  packageStartupMessage("This is loo version ", ver)
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  loo.op <- list(
    loo.cores = parallel::detectCores()
  )
  set_ops <- !(names(loo.op) %in% names(op))
  if (any(set_ops)) options(loo.op[set_ops])
  invisible()
}
