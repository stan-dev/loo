library(loo)

context("psis_approximate_posterior")

# load("tests/testthat/test_data_psis_approximate_posterior.rda")
load("test_data_psis_approximate_posterior.rda")

test_that("Laplace approximation, independent posterior", {
  log_p <- test_data_psis_approximate_posterior$laplace_independent$log_p
  log_q <- test_data_psis_approximate_posterior$laplace_independent$log_q
  ll <- test_data_psis_approximate_posterior$laplace_independent$log_liks
  expect_silent(
    psis_lap <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_lap, "psis")
  expect_lt(pareto_k_values(psis_lap), 0.7)

  expect_silent(
    psis_lap_ll <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, log_liks = ll , cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_lap_ll, "loo")
  expect_true(all(pareto_k_values(psis_lap_ll) < 0.7))
})


test_that("Laplace approximation, correlated posterior", {
  log_p <- test_data_psis_approximate_posterior$laplace_correlated$log_p
  log_q <- test_data_psis_approximate_posterior$laplace_correlated$log_q
  ll <- test_data_psis_approximate_posterior$laplace_correlated$log_liks
  expect_silent(
    psis_lap <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_lap, "psis")
  expect_lt(pareto_k_values(psis_lap), 0.7)

  expect_silent(
    psis_lap_ll <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, log_liks = ll , cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_lap_ll, "loo")
  expect_true(all(pareto_k_values(psis_lap_ll) < 0.7))
})

test_that("Laplace approximation, normal model", {
  log_p <- test_data_psis_approximate_posterior$laplace_normal$log_p
  log_q <- test_data_psis_approximate_posterior$laplace_normal$log_q
  ll <- test_data_psis_approximate_posterior$laplace_normal$log_liks
  expect_warning(
    psis_lap <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_lap, "psis")
  expect_gt(pareto_k_values(psis_lap), 0.5)

  expect_warning(
    psis_lap_ll <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, log_liks = ll , cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_lap_ll, "loo")
  expect_true(all(pareto_k_values(psis_lap_ll) > 0.5))
})



test_that("ADVI fullrank approximation, independent posterior", {
  log_p <- test_data_psis_approximate_posterior$fullrank_independent$log_p
  log_q <- test_data_psis_approximate_posterior$fullrank_independent$log_q
  ll <- test_data_psis_approximate_posterior$fullrank_independent$log_liks
  expect_silent(
    psis_advi <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi, "psis")
  expect_lt(pareto_k_values(psis_advi), 0.7)

  expect_silent(
    psis_advi_ll <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, log_liks = ll , cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi_ll, "loo")
  expect_true(all(pareto_k_values(psis_advi_ll) < 0.7))
})


test_that("ADVI fullrank approximation, correlated posterior", {
  log_p <- test_data_psis_approximate_posterior$fullrank_correlated$log_p
  log_q <- test_data_psis_approximate_posterior$fullrank_correlated$log_q
  ll <- test_data_psis_approximate_posterior$fullrank_correlated$log_liks
  expect_silent(
    psis_advi <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi, "psis")
  expect_lt(pareto_k_values(psis_advi), 0.7)

  expect_silent(
    psis_advi_ll <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, log_liks = ll , cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi_ll, "loo")
  expect_true(all(pareto_k_values(psis_advi_ll) < 0.7))
})

test_that("ADVI fullrank approximation, correlated posterior", {
  log_p <- test_data_psis_approximate_posterior$fullrank_normal$log_p
  log_q <- test_data_psis_approximate_posterior$fullrank_normal$log_q
  ll <- test_data_psis_approximate_posterior$fullrank_normal$log_liks
  expect_warning(
    psis_advi <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi, "psis")
  expect_gt(pareto_k_values(psis_advi), 0.7)

  expect_warning(
    psis_advi_ll <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, log_liks = ll , cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi_ll, "loo")
  expect_true(all(pareto_k_values(psis_advi_ll) > 0.7))
})


test_that("ADVI meanfield approximation, independent posterior", {
  log_p <- test_data_psis_approximate_posterior$meanfield_independent$log_p
  log_q <- test_data_psis_approximate_posterior$meanfield_independent$log_q
  ll <- test_data_psis_approximate_posterior$meanfield_independent$log_liks
  expect_silent(
    psis_advi <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi, "psis")
  expect_lt(pareto_k_values(psis_advi), 0.7)

  expect_silent(
    psis_advi_ll <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, log_liks = ll , cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi_ll, "loo")
  expect_true(all(pareto_k_values(psis_advi_ll) < 0.7))
})


test_that("ADVI meanfield approximation, correlated posterior", {
  log_p <- test_data_psis_approximate_posterior$meanfield_correlated$log_p
  log_q <- test_data_psis_approximate_posterior$meanfield_correlated$log_q
  ll <- test_data_psis_approximate_posterior$meanfield_correlated$log_liks
  expect_warning(
    psis_advi <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi, "psis")
  expect_gt(pareto_k_values(psis_advi), 0.7)

  expect_warning(
    psis_advi_ll <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, log_liks = ll , cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi_ll, "loo")
  expect_true(all(pareto_k_values(psis_advi_ll) > 0.5))
  expect_true(any(pareto_k_values(psis_advi_ll) > 0.7))
})

test_that("ADVI meanfield approximation, normal model", {
  log_p <- test_data_psis_approximate_posterior$meanfield_normal$log_p
  log_q <- test_data_psis_approximate_posterior$meanfield_normal$log_q
  ll <- test_data_psis_approximate_posterior$meanfield_normal$log_liks
  expect_warning(
    psis_advi <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi, "psis")
  expect_gt(pareto_k_values(psis_advi), 0.7)

  expect_warning(
    psis_advi_ll <-
      psis_approximate_posterior(log_p = log_p, log_q = log_q, log_liks = ll , cores = 1, save_psis = FALSE)
  )
  expect_s3_class(psis_advi_ll, "loo")
  expect_true(all(pareto_k_values(psis_advi_ll) > 0.7))
})
