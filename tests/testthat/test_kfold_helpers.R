library(loo)
set.seed(14014)

context("kfold helper functions")

test_that("kfold helpers work", {
  fold_rand <- kfold_split_random(10, 100)
  expect_length(fold_rand, 100)
  expect_equal(sort(unique(fold_rand)), 1:10)
  expect_equal(sum(fold_rand == 2), sum(fold_rand == 9))

  y <- rep(c(0, 1), times = c(10, 190))
  fold_bal <- kfold_split_balanced(5, y)
  expect_true(all(table(fold_bal) == 40))

  y <- rep(c(0, 1), times = c(11, 189))
  fold_bal <- kfold_split_balanced(5, y)
  expect_equal(range(table(fold_bal)), c(39, 41))


  grp <- gl(n = 50, k = 15, labels = state.name)
  fold_strat <- kfold_split_stratified(x = grp)
  expect_true(all(table(fold_strat) == 75))
  expect_equal(sum(table(fold_strat)), length(grp))

  fold_strat <- kfold_split_stratified(K = 9, x = grp)
  expect_false(all(table(fold_strat) == 75))
  expect_equal(sum(table(fold_strat)), length(grp))

  grp <- gl(n = 50, k = 4, labels = state.name)
  grp[grp == "Montana"] <- "Utah"
  fold_strat <- kfold_split_stratified(K = 10, x = grp)
  expect_equal(sum(table(fold_strat)), length(grp) - 4)
})

test_that("kfold helper errors", {
  expect_error(kfold_split_random(10), "!is.null(N) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(10.5, 100), "K == as.integer(K) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(10, 100.5), "N == as.integer(N) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(K = c(1,1), N = 100), "length(K) == 1 is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(N = c(100, 100)), "length(N) == 1 is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(K = 5, N = 4), "K <= N is not TRUE", fixed = TRUE)
  expect_error(kfold_split_random(K = 1, N = 4), "K > 1 is not TRUE", fixed = TRUE)

  y <- sample(c(0, 1), size = 200, replace = TRUE, prob = c(0.05, 0.95))
  expect_error(kfold_split_balanced(10), "!is.null(x) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_balanced(3, c(1,1,1)), "length(unique(x)) == 2 is not TRUE", fixed = TRUE)
  expect_error(kfold_split_balanced(10.5, y), "K == as.integer(K) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_balanced(K = c(1,1), y), "length(K) == 1 is not TRUE", fixed = TRUE)
  expect_error(kfold_split_balanced(K = 201, y), "K <= length(x) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_balanced(K = 1, y), "K > 1 is not TRUE", fixed = TRUE)

  grp <- gl(n = 50, k = 15)
  expect_error(kfold_split_stratified(10), "!is.null(x) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_stratified(3, c(1,1,1)), "'K' must not be bigger than the number of levels/groups in 'x'", fixed = TRUE)
  expect_error(kfold_split_stratified(10.5, grp), "K == as.integer(K) is not TRUE", fixed = TRUE)
  expect_error(kfold_split_stratified(K = c(1,1), grp), "length(K) == 1 is not TRUE", fixed = TRUE)
  expect_error(kfold_split_stratified(K = 1, grp), "K > 1 is not TRUE", fixed = TRUE)
})

