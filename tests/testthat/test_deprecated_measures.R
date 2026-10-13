expect_deprecated <- function(object) {
  testthat::expect_warning(object, "deprecated", ignore.case = TRUE)
}

test_that("elpd() is deprecated", {
  ll <- example_loglik_matrix()
  expect_deprecated(elpd(ll))
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
  expect_deprecated(crps(x1, x2, y))
  expect_deprecated(scrps(x1, x2, y))
})

test_that("loo_crps() is deprecated", {
  set.seed(1)
  y <- rnorm(5)
  x1 <- matrix(rnorm(50), nrow = 10)
  x2 <- matrix(rnorm(50), nrow = 10)
  ll <- matrix(rnorm(50) * 0.1 - 1, nrow = 10)
  expect_warning(
    tryCatch(
      loo_crps(x1, x2, y, ll),
      error = function(e) invisible(NULL)
    ),
    "deprecated",
    ignore.case = TRUE
  )
  expect_warning(
    tryCatch(
      loo_scrps(x1, x2, y, ll),
      error = function(e) invisible(NULL)
    ),
    "deprecated",
    ignore.case = TRUE
  )
})

test_that("loo_predictive_metric() is deprecated", {
  LL <- example_loglik_matrix()
  chain_id <- rep(1:2, each = nrow(LL) / 2)
  r_eff <- relative_eff(exp(LL), chain_id)
  x <- matrix(rnorm(length(LL)), nrow = nrow(LL), ncol = ncol(LL))
  y <- rnorm(ncol(LL))
  expect_deprecated(
    loo_predictive_metric(x, y, LL, metric = "mae", r_eff = r_eff)
  )
})
