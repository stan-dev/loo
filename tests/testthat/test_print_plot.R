library(loo)
options(loo.cores = 2)

context("print, plot, diagnostics")

LLarr <- source(test_path("LL_array_data.R"))$value
waic1 <- suppressWarnings(waic(LLarr))
loo1 <- suppressWarnings(loo(LLarr))
psis1 <- suppressWarnings(psis(LLarr))

test_that("plot methods don't error", {
  plot(loo1, label_points = TRUE)
  plot(psis1, label_points = TRUE)
})

test_that("plot methods throw appropriate errors/warnings", {
  expect_error(plot(waic1), regexp = "No Pareto k estimates found")

  loo1$diagnostics$pareto_k[1:5] <- Inf
  psis1$pareto_k[1:5] <- Inf
  expect_warning(plot(loo1), regexp = "estimates are Inf/NA/NaN and not plotted.")
  expect_warning(plot(psis1), regexp = "estimates are Inf/NA/NaN and not plotted.")
})

test_that("print.waic output is ok",{
  expect_output(
    suppressWarnings(print(waic1)),
    "Computed from 100 by 32 log-likelihood matrix"
  )
  expect_warning(
    capture.output(print(waic1)),
    "p_waic estimates greater than 0.4. We recommend trying loo instead.",
    fixed = TRUE
  )
})

test_that("print.psis_loo and print.psis output ok",{
  expect_output(print(psis1),
                "Computed from 100 by 32 log-likelihood matrix")
  expect_output(print(psis1), "Pareto k diagnostic values")
  expect_output(print(loo1),
                "Computed from 100 by 32 log-likelihood matrix")
  expect_output(print(loo1), "Pareto k diagnostic values")

  loo1$diagnostics$pareto_k <- runif(50, 0, .49)
  expect_output(print(loo1), regexp = "All Pareto k estimates are good")
  loo1$diagnostics$pareto_k[1] <- 0.71
  expect_output(print(loo1), regexp = "Pareto k diagnostic values")
  loo1$diagnostics$pareto_k[1] <- 1.1
  expect_output(print(loo1), regexp = "Pareto k diagnostic values")
})


# pareto_k_[ids,values,table] ---------------------------------------------
test_that("pareto_k_values works for psis_loo and psis objects", {
  kpsis <- pareto_k_values(psis1)
  kloo <- pareto_k_values(loo1)
  expect_identical(kpsis, kloo)
  expect_identical(kpsis, psis1$pareto_k)
})

test_that("pareto_k_ids identifies correct observations", {
  expect_identical(
    pareto_k_ids(loo1, threshold = 0.5),
    pareto_k_ids(psis1, threshold = 0.5)
  )
  expect_identical(
    pareto_k_ids(loo1, threshold = 0.5),
    which(pareto_k_values(loo1) > 0.5)
  )
  expect_identical(
    pareto_k_ids(psis1, threshold = 0.7),
    which(pareto_k_values(psis1) > 0.7)
  )
})

test_that("pareto_k_table gives correct output", {
  k <- pareto_k_values(psis1)
  tab <- pareto_k_table(psis1)

  expect_equal(sum(k <= 0.5), tab[1,1])
  expect_equal(sum(k > 0.5 & k <= 0.7), tab[2,1])
  expect_equal(sum(k > 0.7 & k <= 1), tab[3,1])
  expect_equal(sum(k > 1), tab[4,1])
})
