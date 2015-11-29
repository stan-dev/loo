library(loo)
options(loo.cores = 1)

context("compare")
test_that("compare throws appropriate errors", {
  set.seed(123)
  x1 <- waic(matrix(rnorm(5000), 100, 50))
  x2 <- waic(matrix(rnorm(4000), 100, 40))
  x3 <- waic(matrix(rnorm(3000), 100, 30))
  x4 <- waic(matrix(rnorm(5000), 100, 50))

  expect_error(loo::compare(x1, list(1,2,3)), "class 'loo'")
  expect_error(loo::compare(x1),
               regexp = "requires at least two models")
  expect_error(loo::compare(x1, x2),
               regexp = "same number of data points")
  expect_error(loo::compare(x1, x2, x3),
               regexp = "same number of data points")
  expect_silent(loo::compare(x1, x4))
  expect_silent(comp <- loo::compare(x1, x4, x4))
  expect_equal(comp[2, ], comp[3, ])
})
test_that("compare returns expected result", {
  set.seed(123)
  x1 <- loo(matrix(rnorm(5000), 100, 50))
  x2 <- loo(matrix(rnorm(5000), 100, 50))
  diff_val <- compare(x1, x2)
  diff_ans <- structure(c(-0.49758869857943333059, 1.22288821000758018975,
                          0.62189249931370793600, 0.37810750068629206400),
                        .Names = c("elpd_diff", "se", "weight1", "weight2"),
                        class = "compare.loo")
  expect_equal(diff_val, diff_ans)
})
