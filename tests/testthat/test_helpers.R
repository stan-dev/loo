library(loo)


# test helpers ------------------------------------------------------------
context("helpers")
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
qgpd_val <- qgpd(0.5)
qgpd_ans <- 1

# nlist
a <- 1; b <- 2; c <- 3;
nlist_val1 <- nlist(a, b, c)
nlist_ans1 <- list(a = 1, b = 2, c = 3)
nlist_val2 <- nlist(a, b, c = "tornado")
nlist_ans2 <- list(a = 1, b = 2, c = "tornado")

test_that("helpers return correct values", {
  expect_equivalent(logColMeansExp_val, logColMeansExp_ans)
  expect_equal(qgpd_val, qgpd_ans)
  expect_identical(cbind_list_val, cbind_list_ans)
  expect_equivalent(nlist_val1, nlist_ans1)
  expect_equivalent(nlist_val2, nlist_ans2)
})
