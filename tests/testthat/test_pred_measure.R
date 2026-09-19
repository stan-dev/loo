# load data -----------------------------
res <- readRDS("data-for-tests/test_data_roaches.Rds")
res_sleep_test <- readRDS("data-for-tests/test_data_sleep_cv.Rds")
n_test <- length(res_sleep_test$y_test)

# unit tests ----------------------
test_that("group_ids errors as not yet implemented", {
  expect_error(
    insample_pred_measure(
      ylp = res$ylp, measure = "elpd", group_ids = rep(1:2, 131)
    ),
    "not yet implemented"
  )
})

## .compute_measure() --------------------

.builtin_entry <- function(name) {
  list(name = name, type = "builtin", key = name)
}

test_that(".compute_measure() with elpd works as expected", {
  lppd_i <- .elpd_pointwise(
    source = "insample", ylp = res$ylp, ylp_test = NULL,
    log_weights = NULL, loo = NULL, kfold = NULL, predperf = NULL
  )
  measure_res <- .compute_measure(
    y = NULL,
    ypred = NULL,
    mupred = NULL,
    ylp = res$ylp,
    measure_entry = .builtin_entry("elpd"),
    log_weights = NULL,
    lppd_i = lppd_i
  )

  expect_equal(names(measure_res), c("estimates", "pointwise"))
  expect_equal(measure_res$estimates, measure_elpd(res$ylp)$estimates)
})

test_that(".compute_measure() with rps works as expected", {
  measure_res <- .compute_measure(
    y = res$y,
    ypred = res$ypred,
    mupred = NULL,
    ylp = res$ylp,
    measure_entry = .builtin_entry("rps"),
    log_weights = NULL
  )

  expect_equal(names(measure_res), c("estimates", "pointwise"))
  expect_equal(colnames(measure_res$estimates), c("Estimate", "SE"))
})

test_that(".compute_measure() fails if insufficient input is provided", {
  expect_error(
    .compute_measure(
      y = res$y,
      ypred = res$ypred,
      mupred = NULL,
      ylp = res$ylp,
      measure_entry = .builtin_entry("r2"),
      log_weights = NULL
    ),
    regexp = "`mupred` must be a numeric matrix."
  )
})

## .elpd_pointwise() -------------------------

.elpd_pw <- function(source, ylp = NULL, ylp_test = NULL, log_weights = NULL,
                     loo = NULL, kfold = NULL, predperf = NULL) {
  .elpd_pointwise(source, ylp, ylp_test, log_weights, loo, kfold, predperf)
}

test_that(".elpd_pointwise() reuses the elpd column of predperf", {
  expect_equal(
    .elpd_pw("insample", predperf = res$predperf),
    res$predperf$pointwise[, "elpd"]
  )
})

test_that(".elpd_pointwise() takes elpd from loo and kfold objects", {
  expect_equal(
    .elpd_pw("loo", loo = res$loo), res$loo$pointwise[, "elpd_loo"]
  )
  expect_equal(
    .elpd_pw("kfold", kfold = res$kfold), res$kfold$pointwise[, "elpd_kfold"]
  )
})

test_that(".elpd_pointwise() computes elpd from ylp and ylp_test", {
  expect_equal(
    .elpd_pw("loo", ylp = res$ylp,
             log_weights = res$loo$psis_object$log_weights),
    res$loo$pointwise[, "elpd_loo"],
    ignore_attr = TRUE
  )
  expect_length(
    .elpd_pw("test", ylp_test = res_sleep_test$ylp_test), n_test
  )
})

test_that(".elpd_pointwise() errors if the input is missing", {
  expect_error(.elpd_pw("insample"), regexp = "`ylp` is required")
  expect_error(.elpd_pw("test"), regexp = "`ylp_test` is required")
  expect_error(.elpd_pw("kfold"), regexp = "not stored in")
})

## .get_psis_object() -------------------------

test_that(".get_psis_object() accepts loo and psis_object together", {
  expect_identical(
    .get_psis_object(
      ylp = res$ylp,
      loo = res$loo,
      predperf = NULL,
      psis_object = res$loo$psis_object
    ),
    res$loo$psis_object
  )
})

## .merge_matrix() ---------------------------

test_that(".merge_matrix() works as expected", {
  mat <- matrix(c(1, 2, 3, 4), ncol = 2, nrow = 2)
  val <- c(5, 6)
  expected_name <- "test"

  res <- .merge_matrix(
    source = "insample", mat = mat, name = expected_name,
    values = val, margin = 1
  )

  expect_equal(rownames(res)[[3]], "test")
  expect_equal(dim(res), c(3, 2))
  expect_equal(colnames(res), c("Estimate", "SE"))

  res <- .merge_matrix(
    source = "insample", mat = mat, name = expected_name,
    values = val, margin = 2
  )

  expect_equal(dim(res), c(2, 3))
  expect_equal(colnames(res)[3], "test")
})

test_that(".merge_matrix() with mat = NULL works as expected", {
  res <- .merge_matrix(
    source = "insample", mat = NULL, name = "test",
    values = c(1, 2), margin = 1
  )

  expect_equal(rownames(res), "test")
  expect_equal(colnames(res), c("Estimate", "SE"))
  expect_equal(dim(res), c(1, 2))
})

test_that(".merge_matrix() with loo showes correct names", {
  res <- .merge_matrix(
    source = "loo", mat = NULL, name = "test",
    values = c(1, 2), margin = 1
  )

  expect_equal(rownames(res), "test_loo")
  expect_equal(colnames(res), c("Estimate", "SE"))
  expect_equal(dim(res), c(1, 2))
})

test_that(".merge_matrix() with kfold showes correct names", {
  res <- .merge_matrix(
    source = "kfold", mat = NULL, name = "test",
    values = c(1, 2), margin = 1
  )

  expect_equal(rownames(res), "test_kfold")
  expect_equal(colnames(res), c("Estimate", "SE"))
  expect_equal(dim(res), c(1, 2))
})

test_that("duplicate measure on update warns once and keeps the results", {
  predperf <- insample_pred_measure(ylp = res$ylp)
  warnings <- character()
  withCallingHandlers(
    updated <- pred_measure(
      ylp = res$ylp, predperf = predperf, measure = "elpd"
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  expect_length(warnings, 1L)
  expect_match(warnings[[1]], "already present in")
  expect_equal(rownames(updated$estimates), "elpd")
  expect_equal(colnames(updated$pointwise), "elpd")
})

# integration tests ------------------------------
## loo_pred_measure() / pred_measure() / kfold_pred_measure() ---------

test_that("control scaled = TRUE stores the result as srps, not rps", {
  out <- insample_pred_measure(
    y = res$y,
    ypred = res$ypred,
    ylp = res$ylp,
    measure = "rps",
    control = list(rps = list(scaled = TRUE))
  )

  expect_true("srps" %in% rownames(out$estimates))
  expect_false("rps" %in% rownames(out$estimates))
  expect_true("srps" %in% colnames(out$pointwise))
})

test_that("pred_measure() updates loo results as expected", {
  predperf_loo <- loo_pred_measure(
    loo = res$loo,
    y = res$y,
    mupred = res$mupred,
    ylp = res$ylp,
    measure = c("elpd", "r2", "mse"),
    save_psis = TRUE
  )

  updated_predperf <- pred_measure(
    y = res$y,
    mupred = res$mupred,
    predperf = predperf_loo,
    measure = "mae"
  )

  expect_equal(
    rownames(updated_predperf$estimates),
    c("elpd_loo", "p_loo", "r2_loo", "mse_loo", "mae_loo")
  )
  expect_equal(dim(updated_predperf$estimates), c(5, 2))
})

test_that("pred_measure() keeps dims when the update has no matrix input", {
  predperf_loo <- loo_pred_measure(
    loo = res$loo, y = res$y, mupred = res$mupred, ylp = res$ylp,
    measure = c("elpd", "r2"), save_psis = TRUE
  )
  updated <- pred_measure(predperf = predperf_loo, measure = "mlpd")

  expect_false(is.null(attr(updated, "dims")))
  expect_equal(attr(updated, "dims"), attr(predperf_loo, "dims"))
})

test_that("pred_measure() reuses stored log_weights when save_psis = FALSE", {
  predperf_loo <- loo_pred_measure(
    loo = res$loo, y = res$y, mupred = res$mupred, ylp = res$ylp,
    measure = "r2"
  )
  expect_null(predperf_loo$psis_object)
  expect_false(is.null(predperf_loo$log_weights))

  updated <- pred_measure(
    y = res$y, mupred = res$mupred, predperf = predperf_loo, measure = "mae"
  )
  expect_true("mae_loo" %in% rownames(updated$estimates))
})

test_that("pred_measure() provides warning for duplicate measure", {
  predperf_loo <- loo_pred_measure(
    loo = res$loo,
    y = res$y,
    mupred = res$mupred,
    ylp = res$ylp,
    measure = "r2",
    save_psis = TRUE
  )

  expect_warning(
    pred_measure(
      y = res$y,
      mupred = res$mupred,
      predperf = predperf_loo,
      measure = "r2"
    ),
    regexp = "already present in .* and will be skipped"
  )

  expect_error(
    loo_pred_measure(
      y = res$y,
      mupred = res$mupred,
      ylp = res$ylp,
      loo = res$loo,
      measure = c("mse", "r2", "r2")
    ),
    regexp = "Duplicate measure"
  )
})

test_that("loo_pred_measure() computes expected measures", {
  predperf1 <- loo_pred_measure(
    loo = res$loo,
    y = res$y,
    mupred = res$mupred,
    ylp = res$ylp,
    measure = c("r2", "mse")
  )

  expect_equal(
    rownames(predperf1$estimates),
    c("r2_loo", "mse_loo")
  )
  expect_equal(dim(predperf1$estimates), c(2, 2))
  expect_true(is.loo(predperf1))
})

test_that("loo_pred_measure() has class 'loo' for all input patterns", {
  predperf_loo <- loo_pred_measure(
    loo = res$loo,
    y = res$y,
    ylp = res$ylp
  )
  predperf_ylp_psis <- suppressMessages(loo_pred_measure(
    ylp = res$ylp,
    psis_object = res$loo$psis_object
  ))
  predperf_ylp <- suppressMessages(loo_pred_measure(ylp = res$ylp))

  expect_true(is.loo(predperf_loo))
  expect_true(is.loo(predperf_ylp_psis))
  expect_true(is.loo(predperf_ylp))
  expect_true(is.psis_loo(predperf_loo))
  expect_false(is.psis_loo(predperf_ylp_psis))
  expect_false(is.psis_loo(predperf_ylp))
})

test_that("do_pred_measure() warns if control args are invalid", {
  expect_warning(
    kfold_pred_measure(
      y = res$y,
      ypred = res$ypred,
      mupred = res$mupred,
      ylp = res$ylp,
      measure = c("rps", "srps"),
      kfold = res$kfold,
      control = list(
        rps = list(size = 10)
      )
    ),
    regexp = "Ignoring `size` as it is not a valid argument"
  )
})

test_that("kfold_pred_measure() requires kfold argument", {
  expect_error(
    kfold_pred_measure(
      y = res$y,
      mupred = res$mupred,
      measure = "rmse"
    ),
    regexp = "`kfold` is required"
  )
})

test_that("kfold_pred_measure() works with rps as expected", {
  kfold_res <- kfold_pred_measure(
    y = res$y,
    ypred = res$ypred,
    mupred = res$mupred,
    ylp = res$ylp,
    measure = c("mlpd", "ic", "rps", "srps"),
    kfold = res$kfold
  )

  expect_equal(
    rownames(kfold_res$estimates),
    c("mlpd_kfold", "ic_kfold", "rps_kfold", "srps_kfold")
  )
})

## test_pred_measure() -------------------------------------------------

test_that("test_pred_measure() computes holdout measures as expected", {
  test_res <- test_pred_measure(
    y = res_sleep_test$y_test,
    ypred = res_sleep_test$ypred_test,
    mupred = res_sleep_test$mupred_test,
    ylp_test = res_sleep_test$ylp_test,
    measure = c("elpd", "rmse", "r2")
  )

  expect_s3_class(test_res, "test_pred_measure")
  expect_s3_class(test_res, "pred_measure")
  expect_equal(attr(test_res, "source"), "test")
  expect_equal(
    rownames(test_res$estimates),
    c("elpd_test", "rmse_test", "r2_test")
  )
  expect_equal(dim(test_res$estimates), c(3, 2))
  expect_equal(attr(test_res, "dims"), c(400L, n_test))
  expect_equal(dim(test_res$pointwise), c(n_test, 3L))
})

test_that("test_pred_measure() works with ylp_test only for base summary", {
  test_res <- test_pred_measure(
    y = res_sleep_test$y_test,
    mupred = res_sleep_test$mupred_test,
    ylp_test = res_sleep_test$ylp_test,
    measure = c("elpd", "mae")
  )

  expect_equal(rownames(test_res$estimates), c("elpd_test", "mae_test"))
  expect_equal(nrow(test_res$pointwise), length(res_sleep_test$y_test))
})

test_that("pred_measure() updates test_pred_measure results as expected", {
  test_res <- test_pred_measure(
    y = res_sleep_test$y_test,
    ypred = res_sleep_test$ypred_test,
    mupred = res_sleep_test$mupred_test,
    ylp_test = res_sleep_test$ylp_test,
    measure = "rmse"
  )

  updated <- pred_measure(
    y = res_sleep_test$y_test,
    mupred = res_sleep_test$mupred_test,
    predperf = test_res,
    measure = "mae"
  )

  expect_equal(rownames(updated$estimates), c("rmse_test", "mae_test"))
  expect_equal(attr(updated, "source"), "test")
  expect_equal(dim(updated$pointwise), c(n_test, 2L))
})

# pred_measure() with custom function ------------------------------
test_that("insample_pred_measure() accepts a custom measure function", {
  set.seed(42)
  S <- 4L
  n <- 8L
  y <- rnorm(n)
  mupred <- matrix(rnorm(S * n), nrow = S, ncol = n)
  ylp <- matrix(rnorm(S * n), nrow = S, ncol = n)

  custom_rmse <- function(y, mupred, log_weights = NULL) {
    measure_rmse(y, mupred, log_weights = log_weights)
  }
  attr(custom_rmse, "measure_name") <- "custom_rmse"

  res <- insample_pred_measure(
    y = y,
    mupred = mupred,
    ylp = ylp,
    measure = custom_rmse
  )

  expect_true("custom_rmse" %in% rownames(res$estimates))
  expect_true("custom_rmse" %in% colnames(res$pointwise))
})

test_that("insample_pred_measure() accepts mixed built-in and custom measures", {
  set.seed(1)
  S <- 4L
  n <- 8L
  y <- rnorm(n)
  mupred <- matrix(rnorm(S * n), nrow = S, ncol = n)
  ylp <- matrix(rnorm(S * n), nrow = S, ncol = n)

  custom_rmse <- function(y, mupred, log_weights = NULL) {
    measure_rmse(y, mupred, log_weights = log_weights)
  }
  attr(custom_rmse, "measure_name") <- "custom_rmse"

  res <- insample_pred_measure(
    y = y,
    mupred = mupred,
    ylp = ylp,
    measure = list("r2", custom_rmse = custom_rmse)
  )

  expect_true(all(c("r2", "custom_rmse") %in% rownames(res$estimates)))
})

## elpd on demand -------------------------------------------------------

test_that("measure = NULL reports elpd and p for every source", {
  expect_equal(rownames(insample_pred_measure(ylp = res$ylp)$estimates), "elpd")
  expect_equal(
    rownames(loo_pred_measure(loo = res$loo)$estimates), c("elpd_loo", "p_loo")
  )
  expect_equal(
    rownames(kfold_pred_measure(kfold = res$kfold)$estimates),
    c("elpd_kfold", "p_kfold")
  )
  expect_equal(
    rownames(test_pred_measure(ylp_test = res_sleep_test$ylp_test)$estimates),
    "elpd_test"
  )
})

test_that("elpd and p from loo and kfold objects equal the object estimates", {
  expect_equal(
    loo_pred_measure(loo = res$loo)$estimates,
    res$loo$estimates[c("elpd_loo", "p_loo"), ],
    ignore_attr = TRUE
  )
  expect_equal(
    kfold_pred_measure(kfold = res$kfold)$estimates,
    res$kfold$estimates[c("elpd_kfold", "p_kfold"), ],
    ignore_attr = TRUE
  )
})

test_that("insample_pred_measure() does not need ylp without elpd", {
  x <- insample_pred_measure(y = res$y, mupred = res$mupred, measure = "rmse")
  expect_equal(rownames(x$estimates), "rmse")
})

test_that("loo_pred_measure() keeps diagnostics without elpd", {
  x <- loo_pred_measure(
    loo = res$loo, y = res$y, mupred = res$mupred, measure = "rmse"
  )
  expect_equal(x$diagnostics$pareto_k, res$loo$diagnostics$pareto_k)
})

test_that("pred_measure() recomputes elpd for loo but aborts for kfold", {
  loo_res <- loo_pred_measure(
    loo = res$loo, y = res$y, mupred = res$mupred, measure = "r2",
    save_psis = TRUE
  )
  updated <- suppressMessages(
    pred_measure(ylp = res$ylp, predperf = loo_res, measure = "ic")
  )
  expect_equal(rownames(updated$estimates), c("r2_loo", "ic_loo"))
  expect_equal(
    updated$estimates["ic_loo", "Estimate"],
    -2 * res$loo$estimates["elpd_loo", "Estimate"],
    ignore_attr = TRUE
  )

  kfold_res <- kfold_pred_measure(
    y = res$y, mupred = res$mupred, kfold = res$kfold, measure = "r2"
  )
  expect_error(
    pred_measure(ylp = res$ylp, predperf = kfold_res, measure = "ic"),
    regexp = "not stored in"
  )
})
