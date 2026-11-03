.onAttach <- function(...) {
  ver <- utils::packageVersion("loo")
  packageStartupMessage("This is loo version ", ver)
  packageStartupMessage(
    "- Online documentation and vignettes at mc-stan.org/loo"
  )
  packageStartupMessage(
    "- Parallel execution uses mirai: set up daemons with mirai::daemons() ",
    "before calling loo functions, or pass parallel = FALSE."
  )
}

