library(loo)
library(testthat)
Sys.setenv("R_TESTS" = "")
test_check("loo")
