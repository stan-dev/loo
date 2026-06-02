# load data -----------------------------
res <- readRDS("data-for-tests/res_roaches.Rds")

supported_measures_list <- getFromNamespace("supported_measures_list", "loo")

# helpers for tests ------------------------------------------------
.builtin_entry <- function(name) {
  list(name = name, type = "builtin", key = name)
}

# .normalize_measure() ----------------------------------------------

test_that(".normalize_measure() handles NULL and character input", {
  expect_equal(.normalize_measure(NULL), list())
  entries <- .normalize_measure(c("mse", "rps"))
  expect_length(entries, 2)
  expect_equal(entries[[1]], .builtin_entry("mse"))
})

test_that(".normalize_measure() handles a custom function", {
  f <- function(y, mupred) list(estimate = 1, se = 0, pointwise = y)
  attr(f, "measure_name") <- "custom_mae"
  entries <- .normalize_measure(f)
  expect_length(entries, 1)
  expect_equal(entries[[1]]$type, "custom")
  expect_equal(entries[[1]]$name, "custom_mae")
})

test_that(".normalize_measure() handles a mixed list", {
  f <- function(y, mupred) list(estimate = 1, se = 0, pointwise = y)
  attr(f, "measure_name") <- "custom_mae"
  entries <- .normalize_measure(list("r2", custom_mae = f))
  expect_length(entries, 2)
  expect_equal(entries[[1]]$name, "r2")
  expect_equal(entries[[2]]$name, "custom_mae")
})

test_that(".normalize_measure() errors on duplicate names", {
  expect_error(
    .normalize_measure(c("mse", "mse")),
    regexp = "Duplicate measure"
  )
})

test_that(".normalize_measure() errors on unnamed list function", {
  f <- function(y, mupred) list(estimate = 1, se = 0, pointwise = y)
  expect_error(
    .normalize_measure(list(f)),
    regexp = "must be named"
  )
})

# .prepare_measures() -----------------------------------------------

test_that(".prepare_measures() errors on invalid built-in names", {
  expect_error(
    .prepare_measures("pps", res$predperf, supported_measures_list),
    regexp = "Invalid measure"
  )
})

test_that(".prepare_measures() filters measures already in predperf", {  
  expect_warning(
    .prepare_measures(c("mse", "elpd"), 
    predperf = res$predperf, supported_measures_list),
    regexp = "already present in"
  )

  expect_warning(
    .prepare_measures(c("mse", "elpd", "r2"), 
    predperf = res$predperf, supported_measures_list),
    regexp = "already present in"
  )

  entries <- .prepare_measures(
    c("mse", "rps"),
    predperf = res$predperf,
    supported_measures_list
  )
  expect_equal(vapply(entries, `[[`, "", "name"), c("mse", "rps"))

  entries <- .prepare_measures(
    c("mse"),
    predperf = res$predperf,
    supported_measures_list
  )
  expect_equal(vapply(entries, `[[`, "", "name"), "mse")
})

# .validate_measure_result() ----------------------------------------

test_that(".validate_measure_result() accepts standard and CRPS-style output", {
  res_std <- list(estimate = 1, se = 0.1, pointwise = c(1, 2))
  expect_invisible(.validate_measure_result(res_std, "m", n_obs = 2))

  res_crps <- list(estimates = c(1, 0.1), pointwise = c(1, 2))
  expect_invisible(.validate_measure_result(res_crps, "m", n_obs = 2))
})

test_that(".validate_measure_result() errors on invalid output", {
  expect_error(
    .validate_measure_result(list(estimate = 1), "m"),
    regexp = "Missing"
  )
  expect_error(
    .validate_measure_result(
      list(estimate = 1, se = 0.1, pointwise = c(1, 2, 3)),
      "m",
      n_obs = 2
    ),
    regexp = "length 2, not 3"
  )
})
