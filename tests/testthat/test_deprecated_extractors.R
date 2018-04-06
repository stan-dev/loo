library(loo)
options(mc.cores = 1)
set.seed(123)

context("Depracted extractors")

LLarr <- example_loglik_array()
r_eff_arr <- relative_eff(exp(LLarr))
loo1 <- suppressWarnings(loo(LLarr, r_eff = r_eff_arr))
waic1 <- suppressWarnings(waic(LLarr))

expect_warning_fixed <- function(object, regexp = NULL) {
  expect_warning(object, regexp = regexp, fixed = TRUE)
}

test_that("extracting estimates by name is deprecated for loo objects", {
  # $ method
  expect_equal(
    expect_warning_fixed(loo1$elpd_loo, "elpd_loo using '$' is deprecated"),
    loo1$estimates["elpd_loo", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(loo1$se_elpd_loo, "se_elpd_loo using '$' is deprecated"),
    loo1$estimates["elpd_loo", "SE"]
  )
  expect_equal(
    expect_warning_fixed(loo1$p_loo, "p_loo using '$' is deprecated"),
    loo1$estimates["p_loo", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(loo1$se_p_loo, "se_p_loo using '$' is deprecated"),
    loo1$estimates["p_loo", "SE"]
  )
  expect_equal(
    expect_warning_fixed(loo1$looic, "looic using '$' is deprecated"),
    loo1$estimates["looic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(loo1$se_looic, "se_looic using '$' is deprecated"),
    loo1$estimates["looic", "SE"]
  )

  # [ method
  expect_equal(
    expect_warning_fixed(loo1["elpd_loo"], "elpd_loo using '[' is deprecated")[[1]],
    loo1$estimates["elpd_loo", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(loo1["se_elpd_loo"], "se_elpd_loo using '[' is deprecated")[[1]],
    loo1$estimates["elpd_loo", "SE"]
  )
  expect_equal(
    expect_warning_fixed(loo1["p_loo"], "p_loo using '[' is deprecated")[[1]],
    loo1$estimates["p_loo", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(loo1["se_p_loo"], "se_p_loo using '[' is deprecated")[[1]],
    loo1$estimates["p_loo", "SE"]
  )
  expect_equal(
    expect_warning_fixed(loo1["looic"], "looic using '[' is deprecated")[[1]],
    loo1$estimates["looic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(loo1["se_looic"], "se_looic using '[' is deprecated")[[1]],
    loo1$estimates["looic", "SE"]
  )


  # [[ method
  expect_equal(
    expect_warning_fixed(loo1[["elpd_loo"]], "elpd_loo using '[[' is deprecated"),
    loo1$estimates["elpd_loo", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(loo1[["se_elpd_loo"]], "se_elpd_loo using '[[' is deprecated"),
    loo1$estimates["elpd_loo", "SE"]
  )
  expect_equal(
    expect_warning_fixed(loo1[["p_loo"]], "p_loo using '[[' is deprecated"),
    loo1$estimates["p_loo", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(loo1[["se_p_loo"]], "se_p_loo using '[[' is deprecated"),
    loo1$estimates["p_loo", "SE"]
  )
  expect_equal(
    expect_warning_fixed(loo1[["looic"]], "looic using '[[' is deprecated"),
    loo1$estimates["looic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(loo1[["se_looic"]], "se_looic using '[[' is deprecated"),
    loo1$estimates["looic", "SE"]
  )
})

test_that("extracting estimates by name is deprecated for waic objects", {
  expect_equal(
    expect_warning_fixed(waic1$elpd_waic, "elpd_waic using '$' is deprecated"),
    waic1$estimates["elpd_waic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(waic1$se_elpd_waic, "se_elpd_waic using '$' is deprecated"),
    waic1$estimates["elpd_waic", "SE"]
  )
  expect_equal(
    expect_warning_fixed(waic1$p_waic, "p_waic using '$' is deprecated"),
    waic1$estimates["p_waic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(waic1$se_p_waic, "se_p_waic using '$' is deprecated"),
    waic1$estimates["p_waic", "SE"]
  )
  expect_equal(
    expect_warning_fixed(waic1$waic, "waic using '$' is deprecated"),
    waic1$estimates["waic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(waic1$se_waic, "se_waic using '$' is deprecated"),
    waic1$estimates["waic", "SE"]
  )


  expect_equal(
    expect_warning_fixed(waic1["elpd_waic"], "elpd_waic using '[' is deprecated")[[1]],
    waic1$estimates["elpd_waic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(waic1["se_elpd_waic"], "se_elpd_waic using '[' is deprecated")[[1]],
    waic1$estimates["elpd_waic", "SE"]
  )
  expect_equal(
    expect_warning_fixed(waic1["p_waic"], "p_waic using '[' is deprecated")[[1]],
    waic1$estimates["p_waic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(waic1["se_p_waic"], "se_p_waic using '[' is deprecated")[[1]],
    waic1$estimates["p_waic", "SE"]
  )
  expect_equal(
    expect_warning_fixed(waic1["waic"], "waic using '[' is deprecated")[[1]],
    waic1$estimates["waic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(waic1["se_waic"], "se_waic using '[' is deprecated")[[1]],
    waic1$estimates["waic", "SE"]
  )


  expect_equal(
    expect_warning_fixed(waic1[["elpd_waic"]], "elpd_waic using '[[' is deprecated"),
    waic1$estimates["elpd_waic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(waic1[["se_elpd_waic"]], "se_elpd_waic using '[[' is deprecated"),
    waic1$estimates["elpd_waic", "SE"]
  )
  expect_equal(
    expect_warning_fixed(waic1[["p_waic"]], "p_waic using '[[' is deprecated"),
    waic1$estimates["p_waic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(waic1[["se_p_waic"]], "se_p_waic using '[[' is deprecated"),
    waic1$estimates["p_waic", "SE"]
  )
  expect_equal(
    expect_warning_fixed(waic1[["waic"]], "waic using '[[' is deprecated"),
    waic1$estimates["waic", "Estimate"]
  )
  expect_equal(
    expect_warning_fixed(waic1[["se_waic"]], "se_waic using '[[' is deprecated"),
    waic1$estimates["waic", "SE"]
  )
})
