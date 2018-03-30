library(loo)
set.seed(123)
SW <- suppressWarnings

context("compare")

LLarr <- example_loglik_array()
LLarr2 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 0.5), dim = dim(LLarr))
LLarr3 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 1), dim = dim(LLarr))
w1 <- SW(waic(LLarr))
w2 <- SW(waic(LLarr2))

test_that("compare throws appropriate errors", {
  w3 <- SW(waic(LLarr[,, -1]))
  w4 <- SW(waic(LLarr[,, -(1:2)]))

  expect_error(loo::compare(w1, w2, x = list(w1, w2)),
               regexp = "If 'x' is specified then '...' should not be specified")
  expect_error(loo::compare(w1, list(1,2,3)),
               regexp = "class 'loo'")
  expect_error(loo::compare(w1),
               regexp = "requires at least two models")
  expect_error(loo::compare(x = list(w1)),
               regexp = "requires at least two models")
  expect_error(loo::compare(w1, w3),
               regexp = "same number of data points")
  expect_error(loo::compare(w1, w2, w3),
               regexp = "same number of data points")
  expect_silent(loo::compare(w1, w2))
  expect_silent(loo::compare(w1, w1, w2))
})

test_that("compare returns expected result (2 models)", {
  comp1 <- compare(w1, w1)
  expect_output(print(comp1), "elpd_diff")
  expect_equal(comp1[1:2], c(elpd_diff = 0, se = 0))

  comp2 <- compare(w1, w2)
  expect_equal_to_reference(comp2, "compare_two_models.rds")
  expect_named(comp2, c("elpd_diff", "se"))
  expect_s3_class(comp2, "compare.loo")

  # specifying objects via ... and via arg x gives equal results
  expect_equal(comp2, compare(x = list(w1, w2)))
})

test_that("compare returns expected result (3 models)", {
  w3 <- SW(waic(LLarr3))
  comp1 <- compare(w1, w2, w3)

  expect_equal(
    colnames(comp1),
    c(
      "elpd_diff", "elpd_waic", "se_elpd_waic",
      "p_waic", "se_p_waic", "waic", "se_waic"
    )
    )
  expect_equal(rownames(comp1), c("w1", "w2", "w3"))
  expect_equal(comp1[1,1], 0)
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp1, "matrix")
  expect_equal_to_reference(comp1, "compare_three_models.rds")

  # specifying objects via '...' gives equivalent results (equal
  # except rownames) to using 'x' argument
  expect_equivalent(comp1, compare(x = list(w1, w2, w3)))
})
