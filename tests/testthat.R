library(testthat)
library(loo)
Sys.setenv("R_TESTS" = "")
test_check("loo")
