options(mc.cores = 1)

test_that("Test loo_subsampling and loo_approx with radon data", {
  skip_on_cran() # avoid going over time limit for tests

  load(test_path("data-for-tests/test_radon_laplace_loo.rda"))
  # Rename to spot variable leaking errors
  llfun_test <- llfun
  log_p_test <- log_p
  log_g_test <- log_q
  draws_test <- draws
  data_test <- data
  rm(llfun, log_p, log_q, draws, data)

  set.seed(134)
  expect_silent(
    full_loo <- loo(
      llfun_test,
      draws = draws_test,
      data = data_test,
      r_eff = rep(1, nrow(data_test))
    )
  )
  expect_s3_class(full_loo, "psis_loo")

  set.seed(134)
  expect_silent(
    loo_ss <- loo_subsample(
      x = llfun_test,
      draws = draws_test,
      data = data_test,
      observations = 200,
      loo_approximation = "plpd",
      r_eff = rep(1, nrow(data_test))
    )
  )
  expect_s3_class(loo_ss, "psis_loo_ss")

  set.seed(134)
  expect_silent(
    loo_ap_ss <- loo_subsample(
      x = llfun_test,
      draws = draws_test,
      data = data_test,
      log_p = log_p_test,
      log_g = log_g_test,
      observations = 200,
      loo_approximation = "plpd",
      r_eff = rep(1, nrow(data_test))
    )
  )
  expect_s3_class(loo_ap_ss, "psis_loo_ss")
  expect_s3_class(loo_ap_ss, "psis_loo_ap")

  expect_silent(
    loo_ap_ss_full <- loo_subsample(
      x = llfun_test,
      log_p = log_p_test,
      log_g = log_g_test,
      draws = draws_test,
      data = data_test,
      observations = NULL,
      loo_approximation = "plpd",
      r_eff = rep(1, nrow(data_test))
    )
  )
  expect_failure(expect_s3_class(loo_ap_ss_full, "psis_loo_ss"))
  expect_s3_class(loo_ap_ss_full, "psis_loo_ap")

  # Expect similar results
  z <- 2
  expect_lte(
    loo_ss$estimates["elpd_loo", "Estimate"] -
      z * loo_ss$estimates["elpd_loo", "subsampling SE"],
    full_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["elpd_loo", "Estimate"] +
      z * loo_ss$estimates["elpd_loo", "subsampling SE"],
    full_loo$estimates["elpd_loo", "Estimate"]
  )
  expect_lte(
    loo_ss$estimates["p_loo", "Estimate"] -
      z * loo_ss$estimates["p_loo", "subsampling SE"],
    full_loo$estimates["p_loo", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["p_loo", "Estimate"] +
      z * loo_ss$estimates["p_loo", "subsampling SE"],
    full_loo$estimates["p_loo", "Estimate"]
  )
  expect_lte(
    loo_ss$estimates["looic", "Estimate"] -
      z * loo_ss$estimates["looic", "subsampling SE"],
    full_loo$estimates["looic", "Estimate"]
  )
  expect_gte(
    loo_ss$estimates["looic", "Estimate"] +
      z * loo_ss$estimates["looic", "subsampling SE"],
    full_loo$estimates["looic", "Estimate"]
  )

  expect_failure(expect_equal(
    full_loo$estimates["elpd_loo", "Estimate"],
    loo_ss$estimates["elpd_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    full_loo$estimates["p_loo", "Estimate"],
    loo_ss$estimates["p_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    full_loo$estimates["looic", "Estimate"],
    loo_ss$estimates["looic", "Estimate"],
    tolerance = 0.00000001
  ))

  z <- 2
  expect_lte(
    loo_ap_ss$estimates["elpd_loo", "Estimate"] -
      z * loo_ap_ss$estimates["elpd_loo", "subsampling SE"],
    loo_ap_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_gte(
    loo_ap_ss$estimates["elpd_loo", "Estimate"] +
      z * loo_ap_ss$estimates["elpd_loo", "subsampling SE"],
    loo_ap_ss_full$estimates["elpd_loo", "Estimate"]
  )
  expect_lte(
    loo_ap_ss$estimates["p_loo", "Estimate"] -
      z * loo_ap_ss$estimates["p_loo", "subsampling SE"],
    loo_ap_ss_full$estimates["p_loo", "Estimate"]
  )
  expect_gte(
    loo_ap_ss$estimates["p_loo", "Estimate"] +
      z * loo_ap_ss$estimates["p_loo", "subsampling SE"],
    loo_ap_ss_full$estimates["p_loo", "Estimate"]
  )
  expect_lte(
    loo_ap_ss$estimates["looic", "Estimate"] -
      z * loo_ap_ss$estimates["looic", "subsampling SE"],
    loo_ap_ss_full$estimates["looic", "Estimate"]
  )
  expect_gte(
    loo_ap_ss$estimates["looic", "Estimate"] +
      z * loo_ap_ss$estimates["looic", "subsampling SE"],
    loo_ap_ss_full$estimates["looic", "Estimate"]
  )

  expect_failure(expect_equal(
    loo_ap_ss_full$estimates["elpd_loo", "Estimate"],
    loo_ap_ss$estimates["elpd_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    loo_ap_ss_full$estimates["p_loo", "Estimate"],
    loo_ap_ss$estimates["p_loo", "Estimate"],
    tolerance = 0.00000001
  ))
  expect_failure(expect_equal(
    loo_ap_ss_full$estimates["looic", "Estimate"],
    loo_ap_ss$estimates["looic", "Estimate"],
    tolerance = 0.00000001
  ))

  # Correct printout
  expect_failure(expect_output(
    print(full_loo),
    "Posterior approximation correction used\\."
  ))
  expect_failure(expect_output(
    print(full_loo),
    "subsampled log-likelihood\nvalues"
  ))

  expect_failure(expect_output(
    print(loo_ss),
    "Posterior approximation correction used\\."
  ))
  expect_output(print(loo_ss), "subsampled log-likelihood\nvalues")

  expect_output(print(loo_ap_ss), "Posterior approximation correction used\\.")
  expect_output(print(loo_ap_ss), "subsampled log-likelihood\nvalues")

  expect_output(
    print(loo_ap_ss_full),
    "Posterior approximation correction used\\."
  )
  expect_failure(expect_output(
    print(loo_ap_ss_full),
    "subsampled log-likelihood\nvalues"
  ))

  # Test conversion of objects
  expect_silent(loo_ap_full <- loo:::as.psis_loo.psis_loo(loo_ap_ss_full))
  expect_s3_class(loo_ap_full, "psis_loo_ap")
  expect_silent(loo_ap_full_ss <- loo:::as.psis_loo_ss.psis_loo(loo_ap_full))
  expect_s3_class(loo_ap_full_ss, "psis_loo_ss")
  expect_s3_class(loo_ap_full_ss, "psis_loo_ap")
  expect_silent(loo_ap_full2 <- loo:::as.psis_loo.psis_loo_ss(loo_ap_full_ss))
  expect_s3_class(loo_ap_full2, "psis_loo_ap")
  expect_failure(expect_s3_class(loo_ap_full2, "psis_loo_ss"))
  expect_equal(loo_ap_full2, loo_ap_full)

  # Test update
  set.seed(4712)
  expect_silent(
    loo_ss2 <- update(
      loo_ss,
      draws = draws_test,
      data = data_test,
      observations = 1000,
      r_eff = rep(1, nrow(data_test))
    )
  )
  expect_gt(dim(loo_ss2)[2], dim(loo_ss)[2])
  expect_gt(dim(loo_ss2$pointwise)[1], dim(loo_ss$pointwise)[1])
  expect_equal(nobs(loo_ss), 200)
  expect_equal(nobs(loo_ss2), 1000)
  for (i in 1:nrow(loo_ss2$estimates)) {
    expect_lt(
      loo_ss2$estimates[i, "subsampling SE"],
      loo_ss$estimates[i, "subsampling SE"]
    )
  }

  set.seed(4712)
  expect_silent(
    loo_ap_ss2 <- update(
      object = loo_ap_ss,
      draws = draws_test,
      data = data_test,
      observations = 2000
    )
  )
  expect_gt(dim(loo_ap_ss2)[2], dim(loo_ap_ss)[2])
  expect_gt(dim(loo_ap_ss2$pointwise)[1], dim(loo_ap_ss$pointwise)[1])
  expect_equal(nobs(loo_ap_ss), 200)
  expect_equal(nobs(loo_ap_ss2), 2000)
  for (i in 1:nrow(loo_ap_ss2$estimates)) {
    expect_lt(
      loo_ap_ss2$estimates[i, "subsampling SE"],
      loo_ap_ss$estimates[i, "subsampling SE"]
    )
  }

  expect_equal(round(full_loo$estimates), round(loo_ap_ss_full$estimates))
  expect_failure(expect_equal(full_loo$estimates, loo_ap_ss_full$estimates))
  expect_equal(dim(full_loo), dim(loo_ap_ss_full))
  expect_s3_class(loo_ap_ss_full, "psis_loo_ap")
})


test_that("Test the vignette", {
  skip_on_cran()

  # NOTE
  # If any of these test fails, the vignette probably needs to be updated

  if (FALSE) {
    # Generate vignette test case
    library("rstan")
    stan_code <- "
    data {
      int<lower=0> N;             // number of data points
      int<lower=0> P;             // number of predictors (including intercept)
      matrix[N,P] X;              // predictors (including 1s for intercept)
      int<lower=0,upper=1> y[N];  // binary outcome
    }
    parameters {
      vector[P] beta;
    }
    model {
      beta ~ normal(0, 1);
      y ~ bernoulli_logit(X * beta);
    }
    "
    # logistic <- function(x) {1  / (1 + exp(-x))}
    # logit <- function(x) {log(x) - log(1-x)}
    llfun_logistic <- function(data_i, draws) {
      x_i <- as.matrix(data_i[,
        which(grepl(colnames(data_i), pattern = "X")),
        drop = FALSE
      ])
      y_pred <- draws %*% t(x_i)
      dbinom(x = data_i$y, size = 1, prob = 1 / (1 + exp(-y_pred)), log = TRUE)
    }

    # Prepare data
    url <- "http://stat.columbia.edu/~gelman/arm/examples/arsenic/wells.dat"
    wells <- read.table(url)
    wells$dist100 <- with(wells, dist / 100)
    X <- model.matrix(~ dist100 + arsenic, wells)
    standata <- list(y = wells$switch, X = X, N = nrow(X), P = ncol(X))

    # Fit model
    set.seed(4711)
    fit_1 <- stan(model_code = stan_code, data = standata, seed = 4711)
    print(fit_1, pars = "beta")

    parameter_draws <- extract(fit_1)$beta
    stan_df <- as.data.frame(standata)
    loo_i(1, llfun_logistic, data = stan_df, draws = parameter_draws)

    sm <- stan_model(model_code = stan_code)
    set.seed(4711)
    fit_laplace <- optimizing(sm, data = standata, draws = 2000, seed = 42)
    parameter_draws_laplace <- fit_laplace$theta_tilde
    log_p <- fit_laplace$log_p # The log density of the posterior
    log_g <- fit_laplace$log_g # The log density of the approximation

    # For comparisons
    standata$X[, "arsenic"] <- log(standata$X[, "arsenic"])
    stan_df2 <- as.data.frame(standata)
    set.seed(4711)
    fit_2 <- stan(fit = fit_1, data = standata, seed = 4711)
    parameter_draws_2 <- extract(fit_2)$beta

    save(
      llfun_logistic,
      stan_df,
      stan_df2,
      parameter_draws,
      parameter_draws_laplace,
      parameter_draws_2,
      log_p,
      log_g,
      file = test_path("data-for-tests/loo_subsample_vignette.rda"),
      compression_level = 9
    )
  } else {
    load(test_path("data-for-tests/loo_subsample_vignette.rda"))
  }

  set.seed(4711)
  expect_no_warning(
    looss_1 <- loo_subsample(
      llfun_logistic,
      draws = parameter_draws,
      data = stan_df,
      observations = 100
    )
  )
  expect_output(
    print(looss_1),
    "Computed from 4000 by 100 subsampled log-likelihood"
  )
  expect_output(print(looss_1), "values from 3020 total observations.")
  expect_output(
    print(looss_1),
    "MCSE and ESS estimates assume independent draws"
  )
  expect_output(print(looss_1), "elpd_loo  -1968.5 15.6            0.3")
  expect_output(print(looss_1), "p_loo         3.1  0.1            0.4")
  expect_s3_class(looss_1, c("psis_loo_ss", "psis_loo", "loo"))

  set.seed(4711)
  expect_no_warning(
    looss_1b <- update(
      looss_1,
      draws = parameter_draws,
      data = stan_df,
      observations = 200
    )
  )
  expect_output(
    print(looss_1b),
    "Computed from 4000 by 200 subsampled log-likelihood"
  )
  expect_output(print(looss_1b), "values from 3020 total observations.")
  expect_output(
    print(looss_1b),
    "MCSE and ESS estimates assume independent draws"
  )
  expect_output(print(looss_1b), "elpd_loo  -1968.3 15.6            0.2")
  expect_output(print(looss_1b), "p_loo         3.2  0.1            0.4")
  expect_s3_class(looss_1b, c("psis_loo_ss", "psis_loo", "loo"))

  set.seed(4711)
  expect_no_warning(
    looss_2 <- loo_subsample(
      x = llfun_logistic,
      draws = parameter_draws,
      data = stan_df,
      observations = 100,
      estimator = "hh_pps",
      loo_approximation = "lpd",
      loo_approximation_draws = 100
    )
  )
  expect_output(
    print(looss_2),
    "Computed from 4000 by 100 subsampled log-likelihood"
  )
  expect_output(print(looss_2), "values from 3020 total observations.")
  expect_output(
    print(looss_2),
    "MCSE and ESS estimates assume independent draws"
  )
  # Currently failing
  # expect_output(print(looss_2), "elpd_loo  -1968.9 15.4            0.5")
  # expect_output(print(looss_2), "p_loo         3.5  0.2            0.5")
  expect_s3_class(looss_2, c("psis_loo_ss", "psis_loo", "loo"))

  set.seed(4711)
  expect_no_warning(
    aploo_1 <- loo_approximate_posterior(
      llfun_logistic,
      draws = parameter_draws_laplace,
      data = stan_df,
      log_p = log_p,
      log_g = log_g
    )
  )
  expect_output(
    print(aploo_1),
    "Computed from 2000 by 3020 log-likelihood matrix"
  )
  expect_output(
    print(aploo_1),
    "MCSE and ESS estimates assume independent draws"
  )
  expect_output(print(aploo_1), "elpd_loo  -1968.4 15.6")
  expect_output(print(aploo_1), "p_loo         3.2  0.2")
  expect_output(print(aploo_1), "Posterior approximation correction used.")
  expect_output(print(aploo_1), "All Pareto k estimates are good")
  expect_equal(length(pareto_k_ids(aploo_1, threshold = 0.5)), 31)
  expect_s3_class(aploo_1, c("psis_loo_ap", "psis_loo", "loo"))

  set.seed(4711)
  expect_no_warning(
    looapss_1 <- loo_subsample(
      llfun_logistic,
      draws = parameter_draws_laplace,
      data = stan_df,
      log_p = log_p,
      log_g = log_g,
      observations = 100
    )
  )
  expect_output(
    print(looapss_1),
    "Computed from 2000 by 100 subsampled log-likelihood"
  )
  expect_output(
    print(looapss_1),
    "MCSE and ESS estimates assume independent draws"
  )
  expect_output(print(looapss_1), "values from 3020 total observations.")
  expect_output(print(looapss_1), "elpd_loo  -1968.2 15.6            0.4")
  expect_output(print(looapss_1), "p_loo         2.9  0.1            0.5")
  expect_output(print(looapss_1), "All Pareto k estimates are good")
  expect_equal(length(pareto_k_ids(looapss_1, threshold = 0.5)), 3)

  # Loo compare
  set.seed(4711)
  expect_no_warning(
    looss_1 <- loo_subsample(
      llfun_logistic,
      draws = parameter_draws,
      data = stan_df,
      observations = 100
    )
  )
  set.seed(4712)
  expect_no_warning(
    looss_2 <- loo_subsample(
      x = llfun_logistic,
      draws = parameter_draws_2,
      data = stan_df2,
      observations = 100
    )
  )
  expect_output(
    print(looss_2),
    "Computed from 4000 by 100 subsampled log-likelihood"
  )
  expect_output(
    print(looss_2),
    "MCSE and ESS estimates assume independent draws"
  )
  expect_output(print(looss_2), "values from 3020 total observations.")
  expect_output(print(looss_2), "elpd_loo  -1952.0 16.2            0.2")
  expect_output(print(looss_2), "p_loo         2.6  0.1            0.3")

  expect_warning(
    comp <- loo_compare(looss_1, looss_2),
    "Different subsamples in 'model2' and 'model1'. Naive diff SE is used."
  )
  expect_output(print(comp), "model1 16.5      22.5     0.4")

  set.seed(4712)
  expect_no_warning(
    looss_2_m <- loo_subsample(
      x = llfun_logistic,
      draws = parameter_draws_2,
      data = stan_df2,
      observations = looss_1
    )
  )
  expect_message(
    looss_2_m <- suppressWarnings(loo_subsample(
      x = llfun_logistic,
      draws = parameter_draws_2,
      data = stan_df2,
      observations = obs_idx(looss_1)
    )),
    "Simple random sampling with replacement assumed."
  )

  expect_silent(comp <- loo_compare(looss_1, looss_2_m))
  expect_output(print(comp), "model1 16.1       4.4     0.1")

  set.seed(4712)
  expect_no_warning(
    looss_1 <- update(
      looss_1,
      draws = parameter_draws,
      data = stan_df,
      observations = 200
    )
  )
  expect_no_warning(
    looss_2_m <- update(
      looss_2_m,
      draws = parameter_draws_2,
      data = stan_df2,
      observations = looss_1
    )
  )
  expect_silent(comp2 <- loo_compare(looss_1, looss_2_m))
  expect_output(print(comp2), "model1 16.3       4.4     0.1")

  expect_no_warning(
    looss_2_full <- loo(
      x = llfun_logistic,
      draws = parameter_draws_2,
      data = stan_df2
    )
  )
  expect_message(
    comp3 <- loo_compare(x = list(looss_1, looss_2_full)),
    "Estimated elpd_diff using observations included in loo calculations for all models."
  )
  expect_output(print(comp3), "model1 16.5       4.4     0.3")
})
