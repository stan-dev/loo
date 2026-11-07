if (requireNamespace("mirai", quietly = TRUE)) {
  mirai::daemons(1, dispatcher = FALSE)
  withr::defer(mirai::daemons(0), testthat::teardown_env())
}
