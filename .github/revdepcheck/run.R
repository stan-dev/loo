options(repos = c(CRAN = "https://cloud.r-project.org"))

if (Sys.getenv("PREVIOUS_RUN_ID") == "") {
  revdepcheck::revdep_reset()
}

revdepcheck::revdep_check(
  num_workers = 4,
  timeout = as.difftime(
    as.integer(Sys.getenv("TIMEOUT_MINUTES")),
    units = "mins"
  )
)
