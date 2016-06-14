library(loo)
context("helper functions")

set.seed(123)
x <- matrix(rnorm(100), 20, 5)

test_that("logColMeansExp ok ", {
  expect_equal(logColMeansExp(x), log(colMeans(exp(x))))
})

test_that("cbind_list ok ", {
  xlist <- list()
  for (i in 1:ncol(x))
    xlist[[i]] <- x[,i]
  expect_identical(cbind_list(xlist), x)
})

test_that("qgpd ok ", {
  probs <- seq(from = 0, to = 1, by = 0.25)
  expect_equal(qgpd(probs), c(0, 1/3, 1, 3, Inf))
  expect_true(all(is.nan(qgpd(probs, sigma = 0))))
})

test_that("nlist ok", {
  a <- 1; b <- 2; c <- 3;
  nlist_val <- list(nlist(a, b, c), nlist(a, b, c = "tornado"))
  nlist_ans <- list(list(a = 1, b = 2, c = 3), list(a = 1, b = 2, c = "tornado"))
  expect_equal(nlist_val, nlist_ans)
})

test_that("totals ok", {
  xlist <- as.list(as.data.frame(x))
  totals_val <- unlist(totals(xlist))
  totals_ans <- c(2.83247604949197, -1.02514322943742, 2.12970459700787,
                  -2.39834129191402, 7.50189473847225, 4.34989153078281,
                  3.71159858705083, 4.28135710233064, 4.3518632175463,
                  3.70738038006634)
  expect_equal(totals_val, totals_ans)
})
