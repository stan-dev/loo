library(loo)
options(loo.cores = 1)
set.seed(123)
SW <- suppressWarnings

context("compare_models")

LLarr <- source(test_path("LL_array_data.R"))$value
LLarr2 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 0.5), dim = dim(LLarr))
LLarr3 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 1), dim = dim(LLarr))
w1 <- SW(waic(LLarr))
w2 <- SW(waic(LLarr2))

test_that("compare throws deprecation warning but still works", {
  expect_warning(loo::compare(w1, w2), "deprecated")
  expect_equal(SW(loo::compare(w1, w2)),
               loo::compare_models(w1, w2))
})

test_that("compare_models throws appropriate errors", {
  w3 <- SW(waic(LLarr[,, -1]))
  w4 <- SW(waic(LLarr[,, -(1:2)]))

  expect_error(loo::compare_models(w1, w2, loos = list(w1, w2)),
               regexp = "If 'loos' is specified then '...' should not be specified")
  expect_error(loo::compare_models(w1, list(1,2,3)),
               regexp = "class 'loo'")
  expect_error(loo::compare_models(w1),
               regexp = "requires at least two models")
  expect_error(loo::compare_models(loos = list(w1)),
               regexp = "requires at least two models")
  expect_error(loo::compare_models(w1, w3),
               regexp = "same number of data points")
  expect_error(loo::compare_models(w1, w2, w3),
               regexp = "same number of data points")
  expect_silent(loo::compare_models(w1, w2))
  expect_silent(loo::compare_models(w1, w1, w2))
})

test_that("compare_models returns expected result (2 models)", {
  comp1 <- compare_models(w1, w1)
  expect_equal(comp1[1:2], c(elpd_diff = 0, se = 0))

  comp2 <- compare_models(w1, w2)
  expect_equal_to_reference(comp2, "compare_two_models.rds")
  expect_named(comp2, c("elpd_diff", "se"))
  expect_s3_class(comp2, "compare.loo")

  # specifying objects via ... and via arg loos gives equal results
  expect_equal(comp2, compare_models(loos = list(w1, w2)))
})

test_that("compare_models returns expected result (3 models)", {
  w3 <- SW(waic(LLarr3))
  comp1 <- compare_models(w1, w2, w3)

  expect_equal(
    colnames(comp1),
    c("waic", "se_waic",
      "elpd_waic", "se_elpd_waic",
      "p_waic", "se_p_waic")
    )
  expect_equal(rownames(comp1), c("w1", "w2", "w3"))
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp1, "matrix")
  expect_equal_to_reference(comp1, "compare_three_models.rds")

  # specifying objects via '...' gives equivalent results (equal
  # except rownames) to using 'loos' argument
  expect_equivalent(comp1, compare_models(loos = list(w1, w2, w3)))
})
