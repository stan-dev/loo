.onAttach <- function(...) {
  ver <- utils::packageVersion("loo")
  packageStartupMessage(
    "This is loo version ", ver, ".\n",
    "**NOTE: As of version 2.0.0 loo defaults to 1 core ",
    "but we recommend using as many as possible. ",
    "Use the 'cores' argument or set options(mc.cores = NUM_CORES) for an entire session. ",
    "Windows 10 can be super slow if \code{mc.cores} is set in \code{.rprofile}",
    "(see https://github.com/stan-dev/loo/issues/94)."
    "Visit mc-stan.org/loo/news for details on other changes."
  )
}
