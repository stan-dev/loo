.onAttach <- function(...) {
  ver <- utils::packageVersion("loo")
  packageStartupMessage(
    "This is loo version ", ver, ". ",
    "NOTE: as of version 2.0.0 loo defaults to 1 core ",
    "but we recommend using as many as possible. ",
    "Use the 'cores' argument or set options(loo.cores = VALUE) ",
    "for an entire session."
  )
}
