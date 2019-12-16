library(loo)
set.seed(123)
SW <- suppressWarnings

context("compare models")

LLarr <- example_loglik_array()
LLarr2 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 0.5), dim = dim(LLarr))
LLarr3 <- array(rnorm(prod(dim(LLarr)), c(LLarr), 1), dim = dim(LLarr))
w1 <- SW(waic(LLarr))
w2 <- SW(waic(LLarr2))

test_that("loo_compare throws appropriate errors", {
  w3 <- SW(waic(LLarr[,, -1]))
  w4 <- SW(waic(LLarr[,, -(1:2)]))

  expect_error(loo_compare(2, 3), "must be a list if not a 'loo' object")
  expect_error(loo_compare(w1, w2, x = list(w1, w2)),
               "If 'x' is a list then '...' should not be specified")
  expect_error(loo_compare(w1, list(1,2,3)), "class 'loo'")
  expect_error(loo_compare(w1), "requires at least two models")
  expect_error(loo_compare(x = list(w1)), "requires at least two models")
  expect_error(loo_compare(w1, w3), "same number of data points")
  expect_error(loo_compare(w1, w2, w3), "same number of data points")
})

test_that("loo_compare throws appropriate warnings", {
  w3 <- w1; w4 <- w2
  class(w3) <- class(w4) <- c("kfold", "loo")
  attr(w3, "K") <- 2
  attr(w4, "K") <- 3
  expect_warning(loo_compare(w3, w4), "Not all kfold objects have the same K value")

  class(w4) <- c("psis_loo", "loo")
  attr(w4, "K") <- NULL
  expect_warning(loo_compare(w3, w4), "Comparing LOO-CV to K-fold-CV")

  w3 <- w1; w4 <- w2
  attr(w3, "yhash") <- "a"
  attr(w4, "yhash") <- "b"
  expect_warning(loo_compare(w3, w4), "Not all models have the same y variable")
})



comp_colnames <- c(
  "elpd_diff", "se_diff", "elpd_waic", "se_elpd_waic",
  "p_waic", "se_p_waic", "waic", "se_waic"
)

test_that("loo_compare returns expected results (2 models)", {
  comp1 <- loo_compare(w1, w1)
  expect_s3_class(comp1, "compare.loo")
  expect_equal(colnames(comp1), comp_colnames)
  expect_equal(rownames(comp1), c("model1", "model2"))
  expect_output(print(comp1), "elpd_diff")
  expect_equivalent(comp1[1:2,1], c(0, 0))
  expect_equivalent(comp1[1:2,2], c(0, 0))

  comp2 <- loo_compare(w1, w2)
  expect_s3_class(comp2, "compare.loo")

  # expect_equal_to_reference(comp2, "reference-results/loo_compare_two_models.rds")
  comp2_ref <- readRDS(test_path("reference-results/loo_compare_two_models.rds"))
  expect_equivalent(comp2, comp2_ref)
  expect_equal(colnames(comp2), comp_colnames)

  # specifying objects via ... and via arg x gives equal results
  expect_equal(comp2, loo_compare(x = list(w1, w2)))
})


test_that("loo_compare returns expected result (3 models)", {
  w3 <- SW(waic(LLarr3))
  comp1 <- loo_compare(w1, w2, w3)

  expect_equal(colnames(comp1), comp_colnames)
  expect_equal(rownames(comp1), c("model1", "model2", "model3"))
  expect_equal(comp1[1,1], 0)
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp1, "matrix")

  # expect_equal_to_reference(comp1, "reference-results/loo_compare_three_models.rds")
  comp1_ref <- readRDS(test_path("reference-results/loo_compare_three_models.rds"))
  expect_equivalent(comp1, comp1_ref)

  # specifying objects via '...' gives equivalent results (equal
  # except rownames) to using 'x' argument
  expect_equivalent(comp1, loo_compare(x = list(w1, w2, w3)))
})

# Tests for deprecated compare() ------------------------------------------

test_that("compare throws deprecation warnings", {
  expect_warning(loo::compare(w1, w2), "Deprecated")
  expect_warning(loo::compare(w1, w1, w2), "Deprecated")
})

test_that("compare returns expected result (2 models)", {
  comp1 <- expect_warning(loo::compare(w1, w1), "Deprecated")
  expect_output(print(comp1), "elpd_diff")
  expect_equal(comp1[1:2], c(elpd_diff = 0, se = 0))

  comp2 <- expect_warning(loo::compare(w1, w2), "Deprecated")
  # expect_equal_to_reference(comp2, "reference-results/compare_two_models.rds")
  expect_named(comp2, c("elpd_diff", "se"))
  expect_s3_class(comp2, "compare.loo")

  # specifying objects via ... and via arg x gives equal results
  comp_via_list <- expect_warning(loo::compare(x = list(w1, w2)), "Deprecated")
  expect_equal(comp2, comp_via_list)
})

test_that("compare returns expected result (3 models)", {
  w3 <- SW(waic(LLarr3))
  comp1 <- expect_warning(loo::compare(w1, w2, w3), "Deprecated")

  expect_equal(
    colnames(comp1),
    c(
      "elpd_diff", "se_diff", "elpd_waic", "se_elpd_waic",
      "p_waic", "se_p_waic", "waic", "se_waic"
    ))
  expect_equal(rownames(comp1), c("w1", "w2", "w3"))
  expect_equal(comp1[1,1], 0)
  expect_s3_class(comp1, "compare.loo")
  expect_s3_class(comp1, "matrix")
  # expect_equal_to_reference(comp1, "reference-results/compare_three_models.rds")

  # specifying objects via '...' gives equivalent results (equal
  # except rownames) to using 'x' argument
  comp_via_list <- expect_warning(loo::compare(x = list(w1, w2, w3)), "Deprecated")
  expect_equivalent(comp1, comp_via_list)
})

test_that("compare throws appropriate errors", {
  expect_error(suppressWarnings(loo::compare(w1, w2, x = list(w1, w2))),
               "should not be specified")
  expect_error(suppressWarnings(loo::compare(x = 2)),
               "must be a list")
  expect_error(suppressWarnings(loo::compare(x = list(2))),
               "should have class 'loo'")
  expect_error(suppressWarnings(loo::compare(x = list(w1))),
               "requires at least two models")

  w3 <- SW(waic(LLarr2[,,-1]))
  expect_error(suppressWarnings(loo::compare(x = list(w1, w3))),
               "same number of data points")
  expect_error(suppressWarnings(loo::compare(x = list(w1, w2, w3))),
               "same number of data points")
})
