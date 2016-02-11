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
  expect_output(print(loo::compare(x1, x4)), regexp = "elpd_diff")
  expect_equal(comp[2, ], comp[3, ])
})
test_that("compare returns expected result (2 models)", {
  set.seed(123)
  x1 <- loo(matrix(rnorm(5000), 100, 50))
  x2 <- loo(matrix(rnorm(5000), 100, 50))
  diff_val <- compare(x1, x2)
  # diff_ans <- structure(c(-0.49758869857943333059, 1.22288821000758018975,
  #                         0.62189249931370793600, 0.37810750068629206400),
  #                       .Names = c("elpd_diff", "se", "weight1", "weight2"),
  #                       class = "compare.loo")
  diff_ans <- structure(c(-0.49758869857943333059, 1.22288821000758018975),
                        .Names = c("elpd_diff", "se"),
                        class = "compare.loo")
  expect_equal(diff_val, diff_ans)
})

test_that("compare returns expected result (3 models)", {
  set.seed(123)
  x1 <- loo(matrix(rnorm(5000), 100, 50))
  x2 <- loo(matrix(rnorm(5000), 100, 50))
  x3 <- loo(matrix(rnorm(5000), 100, 50))
  diff_val <- compare(x1, x2, x3)

  diff_ans <-
    structure(c(48.4679657321382, 49.3402690427582, 49.4631431292971,
                1.4135455787813, 1.6461246810119, 2.00130950996329, -24.2339828660691,
                -24.6701345213791, -24.7315715646486, 0.706772789390649, 0.823062340505949,
                1.00065475498164, 48.6639600634871, 49.2628565017408, 49.3678983452134,
                0.990389359906154, 0.948238536607406, 0.986848742971081),
              .Dim = c(3L, 6L),
              .Dimnames = list(c("x1", "x3", "x2"),
                               c("looic", "se_looic", "elpd_loo", "se_elpd_loo", "p_loo", "se_p_loo")),
              class = c("compare.loo", "matrix"))
  expect_equal(diff_val, diff_ans)
})
