# load data -----------------------------
res <- readRDS("data-for-tests/test_data_roaches.Rds")

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

# .validate_control() ---------------------------------------

test_that(".validate_control() accepts valid control silently", {
  expect_invisible(.validate_control(list()))
  expect_invisible(.validate_control(list(rps = list())))
  expect_invisible(.validate_control(list(rps = list(scaled = TRUE))))
  expect_invisible(.validate_control(list(
    rps = list(scaled = TRUE),
    srps = list(revert_sign = TRUE)
  )))
})

test_that(".validate_control() warns on invalid measure args", {
  expect_warning(
    .validate_control(list(rps = list(size = 10))),
    regexp = "Ignoring `size` as it is not a valid argument"
  )

  expect_warning(
    .validate_control(list(rps = list(foo = 1, bar = 2))),
    regexp = "Ignoring `foo` and `bar` as it is not a valid argument"
  )

  expect_warning(
    .validate_control(list(rps = list(scaled = TRUE, bad = 1))),
    regexp = "Ignoring `bad` as it is not a valid argument"
  )

  expect_warning(
    expect_warning(
      .validate_control(list(rps = list(foo = 1), mse = list(bar = 2))),
      regexp = "Ignoring `foo` as it is not a valid argument"
    ),
    regexp = "Ignoring `bar` as it is not a valid argument"
  )
})

test_that(".validate_control() errors on malformed control", {
  expect_error(
    .validate_control("rps"),
    regexp = "must be a named list of named lists."
  )
  expect_error(
    .validate_control(list(list(scaled = TRUE))),
    regexp = "must be a named list of named lists."
  )
  expect_error(
    .validate_control(list(rps = c(scaled = TRUE))),
    regexp = "must be a named list of named lists."
  )
  expect_error(
    .validate_control(list(not_a_function = list(x = 1))),
    regexp = "not_a_function"
  )
})

# subset_measures() -----------------------------------------

.make_measure_result <- function() {
  list(
    estimates = matrix(
      1:4, 2, 2,
      dimnames = list(c("a", "b"), c("Estimate", "SE"))
    ),
    pointwise = matrix(
      1:6, 3, 2,
      dimnames = list(NULL, c("a", "b"))
    ),
    diagnostics = list(pareto_k = c(0.1, 0.2, 0.3)),
    psis_object = list(foo = 1)
  )
}

test_that("subset_measures() subsets kfold and loo base measures", {
  kfold_sub <- subset_measures(
    res$kfold,
    measures = c("elpd_kfold", "p_kfold"),
    components = c("estimates", "pointwise")
  )
  expect_equal(names(kfold_sub), c("estimates", "pointwise"))
  expect_equal(rownames(kfold_sub$estimates), c("elpd_kfold", "p_kfold"))
  expect_equal(colnames(kfold_sub$pointwise), c("elpd_kfold", "p_kfold"))

  loo_sub <- subset_measures(
    res$loo,
    measures = c("elpd_loo", "p_loo"),
    components = c("estimates", "pointwise", "diagnostics")
  )
  expect_equal(names(loo_sub), c("estimates", "pointwise", "diagnostics"))
  expect_equal(rownames(loo_sub$estimates), c("elpd_loo", "p_loo"))
  expect_equal(colnames(loo_sub$pointwise), c("elpd_loo", "p_loo"))
  expect_identical(loo_sub$diagnostics, res$loo$diagnostics)
})

test_that("subset_measures() respects components argument", {
  x <- .make_measure_result()

  estimates_only <- subset_measures(x, measures = c("a", "b"), components = "estimates")
  expect_equal(names(estimates_only), "estimates")
  expect_equal(rownames(estimates_only$estimates), c("a", "b"))

  pointwise_only <- subset_measures(x, measures = "a", components = "pointwise")
  expect_equal(names(pointwise_only), "pointwise")
  expect_equal(colnames(pointwise_only$pointwise), "a")

  diagnostics_only <- subset_measures(x, measures = "a", components = "diagnostics")
  expect_equal(names(diagnostics_only), "diagnostics")
  expect_identical(diagnostics_only$diagnostics, x$diagnostics)
})

test_that("subset_measures() drops unknown measures and components", {
  x <- .make_measure_result()

  expect_error(
    subset_measures(
      x,
      measures = c("a", "missing", "b"),
      components = c("estimates", "pointwise")
    ),
    regexp = "contains invalid value:"
  )
  
  expect_error(
    subset_measures(
      x,
      measures = c("a", "b"),
      components = c("estimates", "pointwise", "measure")
    ),
    regexp = "contains invalid value:"
  )

  empty_measures <- subset_measures(
    x,
    measures = character(0),
    components = c("estimates", "pointwise")
  )
  expect_equal(dim(empty_measures$estimates), c(0, 2))
  expect_equal(dim(empty_measures$pointwise), c(3, 0))
})

test_that("subset_measures() preserves requested measure order", {
  x <- .make_measure_result()

  sub <- subset_measures(
    x,
    measures = c("b", "a"),
    components = c("estimates", "pointwise")
  )

  expect_equal(rownames(sub$estimates), c("b", "a"))
  expect_equal(colnames(sub$pointwise), c("b", "a"))
})