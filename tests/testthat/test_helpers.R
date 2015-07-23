library(loo)


# test helpers ------------------------------------------------------------
context("helpers and utils")
set.seed(123)
x <- matrix(rnorm(100), 20, 5)

# logColMeansExp
logColMeansExp_val <- logColMeansExp(x)
logColMeansExp_ans <- log(colMeans(exp(x)))

# cbind_list
xlist <- list()
for (i in 1:ncol(x)) xlist[[i]] <- x[,i]
cbind_list_val <- cbind_list(xlist)
cbind_list_ans <- x

# qgpd
probs <- seq(from = 0, to = 1, by = 0.25)
qgpd_val <- qgpd(probs)
qgpd_ans <- c(0, 1/3, 1, 3, Inf)

# nlist
a <- 1; b <- 2; c <- 3;
nlist_val <- list(nlist(a, b, c), nlist(a, b, c = "tornado"))
nlist_ans <- list(list(a = 1, b = 2, c = 3), list(a = 1, b = 2, c = "tornado"))

# totals
xlist <- as.list(as.data.frame(x))
totals_val <- unlist(totals(xlist))
totals_ans <- c(2.83247604949197, -1.02514322943742, 2.12970459700787,
                -2.39834129191402, 7.50189473847225, 4.34989153078281,
                3.71159858705083, 4.28135710233064, 4.3518632175463,
                3.70738038006634)

test_that("helpers and utils return correct values", {
  expect_equivalent(logColMeansExp_val, logColMeansExp_ans)
  expect_equal(qgpd_val, qgpd_ans)
  expect_identical(cbind_list_val, cbind_list_ans)
  expect_equivalent(nlist_val, nlist_ans)
  expect_equivalent(totals_val, totals_ans)
})
