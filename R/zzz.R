.onAttach <- function(...) {
  ver <- utils::packageVersion("loo")
  packageStartupMessage("This is loo version ", ver)
  packageStartupMessage(
    "- Online documentation and vignettes at mc-stan.org/loo"
  )
  packageStartupMessage(
    "- As of v2.0.0 loo defaults to 1 core ",
    "but we recommend using as many as possible. ",
    "Use the 'cores' argument or set options(mc.cores = NUM_CORES) ",
    "for an entire session. "
  )
  if (os_is_windows()) {
    packageStartupMessage(
      "- Windows 10 users: loo may be very slow if 'mc.cores' ",
      "is set in your .Rprofile file (see https://github.com/stan-dev/loo/issues/94)."
    )
  }
}




