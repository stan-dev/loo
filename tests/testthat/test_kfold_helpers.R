library(loo)
set.seed(14014)

context("kfold helper functions")

test_that("kfold_split_random works", {
  fold_rand <- kfold_split_random(10, 100)
  expect_length(fold_rand, 100)
  expect_equal(sort(unique(fold_rand)), 1:10)
  expect_equal(sum(fold_rand == 2), sum(fold_rand == 9))
})

test_that("kfold_split_stratified works", {
  y <- rep(c(0, 1), times = c(10, 190))
  fold_strat <- kfold_split_stratified(5, y)
  expect_true(all(table(fold_strat) == 40))

  y <- rep(c(1, 2, 3), times = c(15, 33, 42))
  fold_strat <- kfold_split_stratified(7, y)
  expect_equal(range(table(fold_strat)), c(12, 13))

  y <- mtcars$cyl
  fold_strat <- kfold_split_stratified(10, y)
  expect_equal(range(table(fold_strat)), c(3, 4))
})

test_that("kfold_split_grouped works", {
  grp <- gl(n = 50, k = 15, labels = state.name)
  fold_group <- kfold_split_grouped(x = grp)
  expect_true(all(table(fold_group) == 75))
  expect_equal(sum(table(fold_group)), length(grp))

  fold_group <- kfold_split_grouped(K = 9, x = grp)
  expect_false(all(table(fold_group) == 75))
  expect_equal(sum(table(fold_group)), length(grp))

  grp <- gl(n = 50, k = 4, labels = state.name)
  grp[grp == "Montana"] <- "Utah"
  fold_group <- kfold_split_grouped(K = 10, x = grp)
  expect_equal(sum(table(fold_group)), length(grp) - 4)

  grp <- rep(c("A","B"), each = 20)
  fold_group <- kfold_split_grouped(K = 2, x = grp)
  expect_equal(fold_group, as.integer(as.factor(grp)))
})

test_that("kfold helpers throw correct errors", {
  expect_error(kfold_split_random(10), "!is.null(N) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(10.5, 100), "K == as.integer(K) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(10, 100.5), "N == as.integer(N) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(K = c(1,1), N = 100), "length(K) == 1 is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(N = c(100, 100)), "length(N) == 1 is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(K = 5, N = 4), "K <= N is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(K = 1, N = 4), "K > 1 is not TRUE", fixed = TRUE)

  y <- sample(c(0, 1), size = 200, replace = TRUE, prob = c(0.05, 0.95))
  expect_error(kfold_split_stratified(10), "!is.null(x) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_stratified(10.5, y), "K == as.integer(K) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_stratified(K = c(1,1), y), "length(K) == 1 is not TRUE", fixed = TRUE)
  expect_error(kfold_split_stratified(K = 201, y), "K <= length(x) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_stratified(K = 1, y), "K > 1 is not TRUE", fixed = TRUE)

  grp <- gl(n = 50, k = 15)
  expect_error(kfold_split_grouped(10), "!is.null(x) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_grouped(3, c(1,1,1)), "'K' must not be bigger than the number of levels/groups in 'x'", fixed = TRUE)
  expect_error(kfold_split_grouped(10.5, grp), "K == as.integer(K) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_grouped(K = c(1,1), grp), "length(K) == 1 is not TRUE", fixed = TRUE)
  expect_error(kfold_split_grouped(K = 1, grp), "K > 1 is not TRUE", fixed = TRUE)
})


test_that("print_dims.kfold works", {
  xx <- structure(list(), K = 17, class = c("kfold", "loo"))
  expect_output(print_dims(xx), "Based on 17-fold cross-validation")

  attr(xx, "K") <- NULL
  expect_silent(print_dims(xx))
})

