.onAttach <- function(...) {
  ver <- utils::packageVersion("loo")
  packageStartupMessage("This is loo version ", ver)
  packageStartupMessage(
    "- Online documentation and vignettes at mc-stan.org/loo"
  )
  packageStartupMessage(
    "- Parallelization uses mirai. The 'cores' argument and options('mc.cores') ",
    "are deprecated; see vignette('loo2-parallel') for mirai::daemons() and ",
    "loo_mirai()."
  )
  if (os_is_windows()) {
    packageStartupMessage(
      "- Windows 10 users: loo may be very slow if 'mc.cores' ",
      "is set in your .Rprofile file (see https://github.com/stan-dev/loo/issues/94)."
    )
  }
}




