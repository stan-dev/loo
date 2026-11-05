mirai::daemons(2)
withr::defer(mirai::daemons(0), testthat::teardown_env())

if (identical(Sys.getenv("LOO_TEST_BACKEND"), "mizu")) {
  mizu::mizu_set_default_pool(mizu::mizu_pool(2L))
  withr::defer(mizu::mizu_set_default_pool(NULL), testthat::teardown_env())
}
