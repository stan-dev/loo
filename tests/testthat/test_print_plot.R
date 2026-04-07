set.seed(1414)

LLarr <- example_loglik_array()
waic1 <- suppressWarnings(waic(LLarr))
loo1 <- suppressWarnings(loo(LLarr))
psis1 <- suppressWarnings(psis(-LLarr))
r_eff_arr <- relative_eff(exp(LLarr))
loo1_r_eff <- suppressWarnings(loo(LLarr, r_eff = r_eff_arr))


# plotting ----------------------------------------------------------------
test_that("plot methods don't error", {
  expect_silent(plot(loo1, label_points = FALSE))
  expect_silent(plot(psis1, label_points = TRUE))
  expect_silent(plot(psis1, diagnostic = "n_eff", label_points = FALSE))

  loo1$diagnostics$pareto_k[1] <- 10
  expect_silent(plot(loo1, label_points = TRUE))

  expect_output(print(loo1, plot_k = TRUE))
  expect_output(print(psis1, plot_k = TRUE))
})

test_that("plot methods throw appropriate errors/warnings", {
  expect_error(plot(waic1), regexp = "No Pareto k estimates found")

  loo1$diagnostics$pareto_k[1:5] <- Inf
  psis1$diagnostics$pareto_k[1:5] <- Inf
  expect_warning(
    plot(loo1),
    regexp = "estimates are Inf/NA/NaN and not plotted."
  )
  expect_warning(
    plot(psis1),
    regexp = "estimates are Inf/NA/NaN and not plotted."
  )
})


# printing ----------------------------------------------------------------
lldim_msg <- paste0(
  "Computed from ",
  prod(dim(LLarr)[1:2]),
  " by ",
  dim(LLarr)[3],
  " log-likelihood matrix"
)
lwdim_msg <- paste0(
  "Computed from ",
  prod(dim(LLarr)[1:2]),
  " by ",
  dim(LLarr)[3],
  " log-weights matrix"
)

test_that("print.waic output is ok", {
  expect_output(print(waic1), lldim_msg)
  expect_output(
    print(waic1),
    "p_waic estimates greater than 0.4. We recommend trying loo instead."
  )
})

test_that("print.psis_loo and print.psis output ok", {
  expect_output(print(psis1), lwdim_msg)
  expect_output(print(psis1), "Pareto k estimates are good")
  expect_output(print(loo1), lldim_msg)
  expect_output(print(loo1), "MCSE and ESS estimates assume independent draws")
  expect_output(print(loo1_r_eff), "MCSE and ESS estimates assume MCMC draws")
  expect_output(print(loo1), "Pareto k estimates are good")

  loo1$diagnostics$pareto_k <- psis1$diagnostics$pareto_k <- runif(32, 0, .49)
  expect_output(print(loo1), regexp = "Pareto k estimates are good")
  expect_output(print(psis1), regexp = "Pareto k estimates are good")

  loo1$diagnostics$pareto_k[1] <- psis1$diagnostics$pareto_k[1] <- 0.71
  expect_output(print(loo1), regexp = "Pareto k diagnostic")
  loo1$diagnostics$pareto_k[1] <- psis1$diagnostics$pareto_k[1] <- 1.1
  expect_output(print(loo1), regexp = "Pareto k diagnostic")
})


# pareto_k_[ids,values,table] ---------------------------------------------
test_that("pareto_k_values works for psis_loo and psis objects, errors for waic", {
  kpsis <- pareto_k_values(psis1)
  kloo <- pareto_k_values(loo1)
  expect_identical(kpsis, kloo)
  expect_identical(kpsis, psis1$diagnostics$pareto_k)

  expect_error(pareto_k_values(waic1), "No Pareto k estimates found")
})

test_that("pareto_k_influence_values works for psis_loo objects, errors for psis waic", {
  kloo <- pareto_k_influence_values(loo1)
  kloo2 <- pareto_k_values(loo1)
  expect_identical(kloo, kloo2)

  expect_error(
    pareto_k_influence_values(psis1),
    "No Pareto k influence estimates found"
  )
  expect_error(
    pareto_k_influence_values(waic1),
    "No Pareto k influence estimates found"
  )
})

test_that("pareto_k_ids identifies correct observations", {
  for (j in 1:5) {
    loo1$diagnostics$pareto_k <- psis1$diagnostics$pareto_k <- runif(
      32,
      .25,
      1.25
    )
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
  }
})

test_that("pareto_k_table gives correct output", {
  threshold <- ps_khat_threshold(dim(psis1)[1])
  psis1$diagnostics$pareto_k[1:10] <- runif(10, 0, threshold)
  psis1$diagnostics$pareto_k[11:20] <- runif(10, threshold + 0.01, 0.99)
  psis1$diagnostics$pareto_k[21:32] <- runif(12, 1, 10)
  k <- pareto_k_values(psis1)
  tab <- pareto_k_table(psis1)

  expect_output(print(tab), "Pareto k diagnostic values")
  expect_identical(colnames(tab), c("Count", "Proportion", "Min. n_eff"))
  expect_equal(sum(tab[, "Count"]), length(k))
  expect_equal(sum(tab[, "Proportion"]), 1)

  expect_equal(sum(k <= threshold), tab[1, 1])
  expect_equal(sum(k > threshold & k <= 1), tab[2, 1])
  expect_equal(sum(k > 1), tab[3, 1])

  # if n_eff is NULL
  psis1$diagnostics$n_eff <- NULL
  tab2 <- pareto_k_table(psis1)
  expect_output(print(tab2), "<NA>")
  expect_equal(unname(tab2[, "Min. n_eff"]), rep(NA_real_, 3))

  psis1$diagnostics$pareto_k[1:32] <- 0.4
  expect_output(
    print(pareto_k_table(psis1)),
    paste0("All Pareto k estimates are good (k < ", round(threshold, 2), ")"),
    fixed = TRUE
  )
})


# psis_neff and mcse_loo --------------------------------------------------
test_that("psis_n_eff_values extractor works", {
  n_eff_psis <- psis1$diagnostics$n_eff
  expect_type(n_eff_psis, "double")
  expect_identical(psis_n_eff_values(psis1), n_eff_psis)
  expect_identical(psis_n_eff_values(psis1), psis_n_eff_values(loo1))

  psis1$diagnostics$n_eff <- NULL
  expect_error(psis_n_eff_values(psis1), "No PSIS ESS estimates found")
})

test_that("mcse_loo extractor gives correct value", {
  mcse <- mcse_loo(loo1)
  expect_type(mcse, "double")
  expect_snapshot_value(mcse, style = "serialize")
})

test_that("mcse_loo returns NA when it should", {
  loo1$diagnostics$pareto_k[1] <- 1.5
  mcse <- mcse_loo(loo1)
  expect_equal(mcse, NA)
})

test_that("mcse_loo errors if not psis_loo object", {
  expect_error(mcse_loo(psis1), "psis_loo")
})
