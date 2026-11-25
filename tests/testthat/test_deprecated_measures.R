test_that("elpd() is deprecated", {
  ll <- example_loglik_matrix()
  lifecycle::expect_deprecated(elpd(ll))
})

test_that("elpd() still returns elpd_generic", {
  ll <- example_loglik_matrix()
  out <- suppressWarnings(elpd(ll))
  expect_s3_class(out, "elpd_generic")
})

test_that("elpd(array) warns only once", {
  ll <- example_loglik_array()
  msgs <- character()
  withCallingHandlers(
    elpd(ll),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(msgs, 1L)
  expect_match(msgs[[1]], "deprecated", ignore.case = TRUE)
})

test_that("crps() is deprecated", {
  set.seed(1)
  y <- rnorm(5)
  x1 <- matrix(rnorm(50), nrow = 10)
  x2 <- matrix(rnorm(50), nrow = 10)
  lifecycle::expect_deprecated(crps(x1, x2, y))
  lifecycle::expect_deprecated(scrps(x1, x2, y))
})

test_that("loo_crps() is deprecated", {
  set.seed(1)
  y <- rnorm(5)
  x1 <- matrix(rnorm(50), nrow = 10)
  x2 <- matrix(rnorm(50), nrow = 10)
  ll <- matrix(rnorm(50) * 0.1 - 1, nrow = 10)
  lifecycle::expect_deprecated(
    tryCatch(
      loo_crps(x1, x2, y, ll),
      error = function(e) invisible(NULL)
    )
  )
  lifecycle::expect_deprecated(
    tryCatch(
      loo_scrps(x1, x2, y, ll),
      error = function(e) invisible(NULL)
    )
  )
})

test_that("loo_predictive_metric() is deprecated", {
  LL <- example_loglik_matrix()
  chain_id <- rep(1:2, each = nrow(LL) / 2)
  r_eff <- relative_eff(exp(LL), chain_id)
  x <- matrix(rnorm(length(LL)), nrow = nrow(LL), ncol = ncol(LL))
  y <- rnorm(ncol(LL))
  lifecycle::expect_deprecated(
    loo_predictive_metric(x, y, LL, metric = "mae", r_eff = r_eff)
  )
})

# Equivalence of deprecated and new functions ------------------------------

.equiv_data <- function() {
  set.seed(3)
  S <- 4000
  n <- 20
  y <- rnorm(n)
  mu <- matrix(rnorm(S, mean(y), 1 / sqrt(n)), S, n)
  list(
    y = y,
    x1 = mu + matrix(rnorm(S * n), S),
    x2 = mu + matrix(rnorm(S * n), S),
    ll = dnorm(matrix(y, S, n, byrow = TRUE), mu, 1, log = TRUE)
  )
}

test_that("elpd() matches measure_elpd()", {
  ll <- example_loglik_matrix()
  old <- suppressWarnings(elpd(ll))
  new <- measure_elpd(ylp = ll)
  expect_equal(unname(new$estimates["elpd", ]), unname(old$estimates["elpd", ]))
  expect_equal(unname(new$pointwise[, "elpd"]), unname(old$pointwise[, "elpd"]))
})

test_that("loo_predictive_metric() matches the measure functions", {
  d <- .equiv_data()
  lw <- weights(psis(-d$ll, r_eff = 1), log = TRUE)
  new_fun <- list(mae = measure_mae, rmse = measure_rmse, mse = measure_mse)
  for (m in names(new_fun)) {
    old <- suppressWarnings(loo_predictive_metric(d$x1, d$y, d$ll, metric = m))
    new <- new_fun[[m]](d$y, d$x1, log_weights = lw)
    expect_equal(unname(new$estimates[1, ]), c(old$estimate, old$se))
  }
  yb <- as.integer(d$y > 0)
  xb <- plogis(d$x1)
  old <- suppressWarnings(loo_predictive_metric(xb, yb, d$ll, metric = "acc"))
  new <- measure_acc(yb, xb, log_weights = lw)
  expect_equal(unname(new$estimates[1, ]), c(old$estimate, old$se))
  old <- suppressWarnings(
    loo_predictive_metric(xb, yb, d$ll, metric = "balanced_acc")
  )
  new <- measure_bacc(yb, xb, log_weights = lw)
  expect_equal(new$estimates[1, "Estimate"], old$estimate)
})

test_that("crps() and scrps() match measure_rps() and measure_srps()", {
  d <- .equiv_data()
  old_crps <- suppressWarnings(crps(d$x1, d$x2, d$y))
  old_scrps <- suppressWarnings(scrps(d$x1, d$x2, d$y))
  expect_equal(
    measure_rps(d$y, d$x1)$estimates[1, "Estimate"],
    -old_crps$estimates[["Estimate"]],
    tolerance = 0.01
  )
  expect_equal(
    measure_srps(d$y, d$x1)$estimates[1, "Estimate"],
    old_scrps$estimates[["Estimate"]],
    tolerance = 0.01
  )
})

test_that("loo_crps() and loo_scrps() match loo_pred_measure()", {
  d <- .equiv_data()
  new <- suppressMessages(loo_pred_measure(
    y = d$y, ypred = d$x1, ylp = d$ll, measures = c("rps", "srps")
  ))
  old_crps <- suppressWarnings(loo_crps(d$x1, d$x2, d$y, d$ll))
  old_scrps <- suppressWarnings(loo_scrps(d$x1, d$x2, d$y, d$ll))
  expect_equal(
    new$estimates["rps_loo", "Estimate"],
    -old_crps$estimates[["Estimate"]],
    tolerance = 0.01
  )
  expect_equal(
    new$estimates["srps_loo", "Estimate"],
    old_scrps$estimates[["Estimate"]],
    tolerance = 0.01
  )
})
