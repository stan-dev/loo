test_that("package startup message works", {
  skip_on_cran()
  source("../../R/zzz.R")
  expect_message(.onAttach(), "As of version 2.0.0 loo defaults to 1 core")
})
