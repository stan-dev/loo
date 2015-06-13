.onAttach <- function(...) {
  ver <- utils::packageVersion("loo")
  msg <- paste("loo version", ver)
  packageStartupMessage(msg)
}
