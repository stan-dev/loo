options(mc.cores = 1)
set.seed(123)

LLarr <- example_loglik_array()
r_eff_arr <- relative_eff(exp(LLarr))
loo1 <- suppressWarnings(loo(LLarr, r_eff = r_eff_arr))
waic1 <- suppressWarnings(waic(LLarr))
test_that("extracting estimates by name is deprecated for loo objects", {
  # $ method
  expect_snapshot(loo1$elpd_loo)
  expect_equal(
    suppressWarnings(loo1$elpd_loo),
    loo1$estimates["elpd_loo", "Estimate"]
  )
  expect_snapshot(loo1$se_elpd_loo)
  expect_equal(
    suppressWarnings(loo1$se_elpd_loo),
    loo1$estimates["elpd_loo", "SE"]
  )
  expect_snapshot(loo1$p_loo)
  expect_equal(
    suppressWarnings(loo1$p_loo),
    loo1$estimates["p_loo", "Estimate"]
  )
  expect_snapshot(loo1$se_p_loo)
  expect_equal(
    suppressWarnings(loo1$se_p_loo),
    loo1$estimates["p_loo", "SE"]
  )
  expect_snapshot(loo1$looic)
  expect_equal(
    suppressWarnings(loo1$looic),
    loo1$estimates["looic", "Estimate"]
  )
  expect_snapshot(loo1$se_looic)
  expect_equal(
    suppressWarnings(loo1$se_looic),
    loo1$estimates["looic", "SE"]
  )

  # [ method
  expect_snapshot(loo1["elpd_loo"])
  expect_equal(
    suppressWarnings(loo1["elpd_loo"][[1]]),
    loo1$estimates["elpd_loo", "Estimate"]
  )
  expect_snapshot(loo1["se_elpd_loo"])
  expect_equal(
    suppressWarnings(loo1["se_elpd_loo"][[1]]),
    loo1$estimates["elpd_loo", "SE"]
  )
  expect_snapshot(loo1["p_loo"])
  expect_equal(
    suppressWarnings(loo1["p_loo"][[1]]),
    loo1$estimates["p_loo", "Estimate"]
  )
  expect_snapshot(loo1["se_p_loo"])
  expect_equal(
    suppressWarnings(loo1["se_p_loo"][[1]]),
    loo1$estimates["p_loo", "SE"]
  )
  expect_snapshot(loo1["looic"])
  expect_equal(
    suppressWarnings(loo1["looic"][[1]]),
    loo1$estimates["looic", "Estimate"]
  )
  expect_snapshot(loo1["se_looic"])
  expect_equal(
    suppressWarnings(loo1["se_looic"][[1]]),
    loo1$estimates["looic", "SE"]
  )

  # [[ method
  expect_snapshot(loo1[["elpd_loo"]])
  expect_equal(
    suppressWarnings(loo1[["elpd_loo"]]),
    loo1$estimates["elpd_loo", "Estimate"]
  )
  expect_snapshot(loo1[["se_elpd_loo"]])
  expect_equal(
    suppressWarnings(loo1[["se_elpd_loo"]]),
    loo1$estimates["elpd_loo", "SE"]
  )
  expect_snapshot(loo1[["p_loo"]])
  expect_equal(
    suppressWarnings(loo1[["p_loo"]]),
    loo1$estimates["p_loo", "Estimate"]
  )
  expect_snapshot(loo1[["se_p_loo"]])
  expect_equal(
    suppressWarnings(loo1[["se_p_loo"]]),
    loo1$estimates["p_loo", "SE"]
  )
  expect_snapshot(loo1[["looic"]])
  expect_equal(
    suppressWarnings(loo1[["looic"]]),
    loo1$estimates["looic", "Estimate"]
  )
  expect_snapshot(loo1[["se_looic"]])
  expect_equal(
    suppressWarnings(loo1[["se_looic"]]),
    loo1$estimates["looic", "SE"]
  )
})

test_that("accessing diagnostics$n_eff is deprecated, diagnostics$ess works", {
  expect_s3_class(loo1$diagnostics, "loo_diagnostics")

  # ess should work without warning
  expect_no_warning(loo1$diagnostics$ess)
  expect_no_warning(loo1$diagnostics[["ess"]])
  expect_no_warning(loo1$diagnostics["ess"])

  # n_eff should warn
  expect_warning(loo1$diagnostics$n_eff, "deprecated")
  expect_warning(loo1$diagnostics[["n_eff"]], "deprecated")
  expect_warning(loo1$diagnostics["n_eff"], "deprecated")

  # both should return the same values
  ess_vals <- loo1$diagnostics$ess
  n_eff_vals <- suppressWarnings(loo1$diagnostics$n_eff)
  expect_equal(ess_vals, n_eff_vals)

  # psis_n_eff_values should return ess without warning
  expect_no_warning(psis_n_eff_values(loo1))
  expect_equal(psis_n_eff_values(loo1), ess_vals)
})

test_that("assigning to diagnostics$ess keeps n_eff in sync", {
  loo_copy <- loo1
  new_vals <- rep(100, length(loo_copy$diagnostics$ess))
  loo_copy$diagnostics$ess <- new_vals
  expect_equal(loo_copy$diagnostics$ess, new_vals)
  n_eff_vals <- suppressWarnings(loo_copy$diagnostics$n_eff)
  expect_equal(n_eff_vals, new_vals)
})

test_that("extracting estimates by name is deprecated for waic objects", {
  expect_snapshot(waic1$elpd_waic)
  expect_equal(
    suppressWarnings(waic1$elpd_waic),
    waic1$estimates["elpd_waic", "Estimate"]
  )
  expect_snapshot(waic1$se_elpd_waic)
  expect_equal(
    suppressWarnings(waic1$se_elpd_waic),
    waic1$estimates["elpd_waic", "SE"]
  )
  expect_snapshot(waic1$p_waic)
  expect_equal(
    suppressWarnings(waic1$p_waic),
    waic1$estimates["p_waic", "Estimate"]
  )
  expect_snapshot(waic1$se_p_waic)
  expect_equal(
    suppressWarnings(waic1$se_p_waic),
    waic1$estimates["p_waic", "SE"]
  )
  expect_snapshot(waic1$waic)
  expect_equal(
    suppressWarnings(waic1$waic),
    waic1$estimates["waic", "Estimate"]
  )
  expect_snapshot(waic1$se_waic)
  expect_equal(
    suppressWarnings(waic1$se_waic),
    waic1$estimates["waic", "SE"]
  )

  # [ method
  expect_snapshot(waic1["elpd_waic"])
  expect_equal(
    suppressWarnings(waic1["elpd_waic"][[1]]),
    waic1$estimates["elpd_waic", "Estimate"]
  )
  expect_snapshot(waic1["se_elpd_waic"])
  expect_equal(
    suppressWarnings(waic1["se_elpd_waic"][[1]]),
    waic1$estimates["elpd_waic", "SE"]
  )
  expect_snapshot(waic1["p_waic"])
  expect_equal(
    suppressWarnings(waic1["p_waic"][[1]]),
    waic1$estimates["p_waic", "Estimate"]
  )
  expect_snapshot(waic1["se_p_waic"])
  expect_equal(
    suppressWarnings(waic1["se_p_waic"][[1]]),
    waic1$estimates["p_waic", "SE"]
  )
  expect_snapshot(waic1["waic"])
  expect_equal(
    suppressWarnings(waic1["waic"][[1]]),
    waic1$estimates["waic", "Estimate"]
  )
  expect_snapshot(waic1["se_waic"])
  expect_equal(
    suppressWarnings(waic1["se_waic"][[1]]),
    waic1$estimates["waic", "SE"]
  )

  # [[ method
  expect_snapshot(waic1[["elpd_waic"]])
  expect_equal(
    suppressWarnings(waic1[["elpd_waic"]]),
    waic1$estimates["elpd_waic", "Estimate"]
  )
  expect_snapshot(waic1[["se_elpd_waic"]])
  expect_equal(
    suppressWarnings(waic1[["se_elpd_waic"]]),
    waic1$estimates["elpd_waic", "SE"]
  )
  expect_snapshot(waic1[["p_waic"]])
  expect_equal(
    suppressWarnings(waic1[["p_waic"]]),
    waic1$estimates["p_waic", "Estimate"]
  )
  expect_snapshot(waic1[["se_p_waic"]])
  expect_equal(
    suppressWarnings(waic1[["se_p_waic"]]),
    waic1$estimates["p_waic", "SE"]
  )
  expect_snapshot(waic1[["waic"]])
  expect_equal(
    suppressWarnings(waic1[["waic"]]),
    waic1$estimates["waic", "Estimate"]
  )
  expect_snapshot(waic1[["se_waic"]])
  expect_equal(
    suppressWarnings(waic1[["se_waic"]]),
    waic1$estimates["waic", "SE"]
  )
})
