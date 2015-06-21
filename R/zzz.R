.onAttach <- function(...) {
  ver <- utils::packageVersion("loo")
  msg <- paste("This is loo version", ver)
  packageStartupMessage(msg)
}
