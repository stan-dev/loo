options(repos = c(CRAN = "https://cloud.r-project.org"))

if (Sys.getenv("PREVIOUS_RUN_ID") == "") {
  revdepcheck::revdep_reset()
}

revdepcheck::revdep_check(
  num_workers = 2,
  timeout = as.difftime(
    as.integer(Sys.getenv("TIMEOUT_MINUTES")),
    units = "mins"
  ),
  quiet = FALSE,
  env = c(
    revdepcheck::revdep_env_vars(),
    R_LIBS_SITE = normalizePath(.libPaths()[1], mustWork = TRUE)
  )
)

