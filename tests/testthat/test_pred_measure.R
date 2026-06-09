# load data -----------------------------
res <- readRDS("data-for-tests/res_roaches.Rds")

y <- res$y
ypred <- res$ypred
mupred <- res$mupred
ylp <- res$ylp
kfold1 <- res$kfold
loo1 <- res$loo
predperf <- res$predperf

# unit tests ----------------------
## .compute_measure() --------------------

.builtin_entry <- function(name) {
  list(name = name, type = "builtin", key = name)
}

test_that(".compute_measure() with r2 works as expected", {
  measure_res <- .compute_measure(
    y = y,
    ypred = NULL,
    mupred = mupred,
    ylp = ylp,
    measure_entry = .builtin_entry("r2"),
    log_weights = NULL
  )

  expect_equal(names(measure_res), c("estimates", "pointwise"))
})

test_that(".compute_measure() with rps works as expected", {
  measure_res <- .compute_measure(
    y = y,
    ypred = ypred,
    mupred = NULL,
    ylp = ylp,
    measure_entry = .builtin_entry("rps"),
    log_weights = NULL
  )

  expect_equal(names(measure_res), c("estimates", "pointwise"))
  expect_equal(names(measure_res$estimates), c("Estimate", "SE"))
})

test_that(".compute_measure() with elpd works as expected", {
  measure_res <- .compute_measure(
    y = NULL,
    ypred = NULL,
    mupred = NULL,
    ylp = ylp,
    measure_entry = .builtin_entry("elpd"),
    log_weights = NULL
  )

  expect_equal(names(measure_res), c("estimates", "pointwise"))
})

test_that(".compute_measure() fails if insufficient input is provided", {
  expect_error(
    .compute_measure(
      y = y,
      ypred = ypred,
      mupred = NULL,
      ylp = ylp,
      measure_entry = .builtin_entry("r2"),
      log_weights = NULL
    ),
    regexp = "`mupred` must be a numeric matrix."
  )
})

## .compute_base_measure() -------------------------

test_that(".compute_base_measure() returns predperf if already existent", {
  res <- .compute_base_measure(
    ylp = NULL,
    ylp_test = NULL,
    loo = NULL,
    kfold = NULL,
    predperf = predperf,
    psis_object = NULL,
    source = "insample"
  )
  expect_equal(res, predperf)
})

test_that(".compute_base_measure() errors if missing input", {
  expect_error(
    .compute_base_measure(
      ylp = NULL,
      ylp_test = NULL,
      loo = NULL,
      kfold = NULL,
      predperf = NULL,
      psis_object = NULL,
      source = "insample"
    ),
    regexp = "`ylp` must be a numeric matrix."
  )
})

test_that(".compute_base_measure() computes elpd as expected", {
  base_measure <- .compute_base_measure(
    ylp = ylp,
    ylp_test = NULL,
    loo = NULL,
    kfold = NULL,
    predperf = NULL,
    psis_object = NULL,
    source = "insample"
  )

  expect_equal(names(base_measure), c("estimates", "pointwise", "diagnostics"))
  expect_null(base_measure$diagnostics)
  expect_equal(rownames(base_measure$estimates), "elpd")
  expect_equal(colnames(base_measure$estimates), c("Estimate", "SE"))
  expect_equal(dimnames(base_measure$pointwise)[[2]], "elpd")
})

test_that(".compute_base_measure() computes elpd_loo as expected", {
  base_measure <- .compute_base_measure(
    ylp = ylp,
    ylp_test = NULL,
    loo = NULL,
    kfold = NULL,
    predperf = NULL,
    psis_object = loo1$psis_object,
    source = "loo"
  )

  expect_equal(
    rownames(base_measure$estimates),
    c("elpd_loo", "p_eff_loo")
  )
  expect_equal(
    dimnames(base_measure$pointwise)[[2]],
    c("elpd_loo", "p_eff_loo")
  )
})

test_that(".compute_base_measure() computes elpd_kfold as expected", {
  base_measure <- .compute_base_measure(
    ylp = ylp,
    ylp_test = NULL,
    loo = NULL,
    kfold = kfold1,
    predperf = NULL,
    psis_object = NULL,
    source = "kfold"
  )

  expect_equal(
    rownames(base_measure$estimates),
    c("elpd_kfold", "p_kfold")
  )
  expect_equal(
    dimnames(base_measure$pointwise)[[2]],
    c("elpd_kfold", "p_kfold")
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

test_that(".merge_matrix() with duplicate naming is skipped", {
  mat <- matrix(c(1, 2, 3, 4), ncol = 2, nrow = 2)
  rownames(mat) <- c("elpd", "mse")

  res <- .merge_matrix(
    source = "insample", mat = mat, name = "rmse",
    values = c(5, 6), margin = 1
  )

  expect_equal(c("elpd", "mse", "rmse"), rownames(res))
  
  expect_warning(
    .merge_matrix(
      source = "insample", mat = mat, name = "elpd",
      values = c(5, 6), margin = 1
    ),
    regexp = "already present in results. Skipping the update."
  )
})

# integration tests ------------------------------
## loo_pred_measure() / pred_measure() / kfold_pred_measure() ---------

test_that("pred_measure() updates loo results as expected", {
  predperf_loo <- loo_pred_measure(
    loo = loo1,
    y = y,
    mupred = mupred,
    ylp = ylp,
    measure = c("r2", "mse"),
    save_psis = TRUE
  )

  updated_predperf <- pred_measure(
    y = y,
    mupred = mupred,
    predperf = predperf_loo,
    measure = "mae"
  )

  expect_equal(
    rownames(updated_predperf$estimates),
    c("elpd_loo", "p_loo", "r2_loo", "mse_loo", "mae_loo")
  )
  expect_equal(dim(updated_predperf$estimates), c(5, 2))
})

test_that("pred_measure() provides warning for duplicate measure", {
  predperf_loo <- loo_pred_measure(
    loo = loo1,
    y = y,
    mupred = mupred,
    ylp = ylp,
    measure = "r2",
    save_psis = TRUE
  )

  expect_warning(
    pred_measure(
      y = y,
      mupred = mupred,
      predperf = predperf_loo,
      measure = "r2"
    ),
    regexp = "already present in results. Skipping the update."
  )

  expect_error(
    loo_pred_measure(
      y = y,
      mupred = mupred,
      ylp = ylp,
      loo = loo1,
      measure = c("mse", "r2", "r2")
    ),
    regexp = "Duplicate measure"
  )
})

test_that("loo_pred_measure() computes expected measures", {
  predperf1 <- loo_pred_measure(
    loo = loo1,
    y = y,
    mupred = mupred,
    ylp = ylp,
    measure = c("r2", "mse")
  )

  expect_equal(
    rownames(predperf1$estimates),
    c("elpd_loo", "p_loo", "r2_loo", "mse_loo")
  )
  expect_equal(dim(predperf1$estimates), c(4, 2))
})

test_that("pred_measure_engine() warns if control args are invalid", {
  expect_warning(
    kfold_pred_measure(
      y = y,
      ypred = ypred,
      mupred = mupred,
      ylp = ylp,
      measure = c("rps", "srps"),
      kfold = kfold1,
      control = list(
        rps = list(size = 10)
      )
    ),
    regexp = "Ignoring `size` as it is not a valid argument"
  )
})

test_that("kfold_pred_measure() works with rps as expected", {
  res <- kfold_pred_measure(
    y = y,
    ypred = ypred,
    mupred = mupred,
    ylp = ylp,
    measure = c("mlpd", "ic" ,"rps", "srps"),
    kfold = kfold1
  )

  expect_equal(
    rownames(res$estimates),
    c("elpd_kfold", "p_kfold", "mlpd_kfold", "ic_kfold", "rps_kfold", "srps_kfold")
  )
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
    rmse(y, mupred, log_weights = log_weights)
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
    rmse(y, mupred, log_weights = log_weights)
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