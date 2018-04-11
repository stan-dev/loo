.onAttach <- function(...) {
  ver <- utils::packageVersion("loo")
  packageStartupMessage(
    "This is loo version ", ver, ".\n",
    "**NOTE: As of version 2.0.0 loo defaults to 1 core ",
    "but we recommend using as many as possible. ",
    "Use the 'cores' argument or set options(mc.cores = NUM_CORES) for an entire session. ",
    "Visit mc-stan.org/loo/news for details on other changes."
  )
}
