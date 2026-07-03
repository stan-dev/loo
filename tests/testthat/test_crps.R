set.seed(123456789)
n <- 10
S <- 100
y <- rnorm(n)
x1 <- matrix(rnorm(n * S), nrow = S)
x2 <- matrix(rnorm(n * S), nrow = S)
ll <- matrix(rnorm(n * S) * 0.1 - 1, nrow = S)

with_seed <- function(seed, code) {
  code <- substitute(code)
  orig.seed <- .Random.seed
  on.exit(.Random.seed <<- orig.seed)
  set.seed(seed)
  eval.parent(code)
}

test_that("crps computation is correct", {
  expect_equal(.crps_fun(2.0, 1.0), 0.0)
  expect_equal(.crps_fun(1.0, 2.0), -1.5)
  expect_equal(.crps_fun(pi, pi^2), 0.5 * pi - pi^2)

  expect_equal(.crps_fun(1.0, 0.0, scale = TRUE), 0.0)
  expect_equal(.crps_fun(1.0, 2.0, scale = TRUE), -2.0)
  expect_equal(.crps_fun(pi, pi^2, scale = TRUE), -pi^2/pi - 0.5 * log(pi))
})

test_that("crps matches snapshots", {
  expect_snapshot_value(
    with_seed(1, suppressWarnings(crps(x1, x2, y))),
    style = "serialize"
  )
  expect_snapshot_value(
    with_seed(1, suppressWarnings(scrps(x1, x2, y))),
    style = "serialize"
  )
  expect_snapshot_value(
    with_seed(1, suppressWarnings(loo_crps(x1, x2, y, ll))),
    style = "serialize"
  )
  expect_snapshot_value(
    with_seed(1, suppressWarnings(loo_scrps(x1, x2, y, ll))),
    style = "serialize"
  )
})

test_that("input validation throws correct errors", {
  expect_error(validate_crps_input(as.character(x1), x2, y),
               "is.numeric(x) is not TRUE",
               fixed = TRUE)
  expect_error(validate_crps_input(x1, as.character(x2), y),
               "is.numeric(x2) is not TRUE",
               fixed = TRUE)
  expect_error(validate_crps_input(x1, x2, c('a', 'b')),
               "is.numeric(y) is not TRUE",
               fixed = TRUE)
  expect_error(validate_crps_input(x1, t(x2), y),
               "identical(dim(x), dim(x2)) is not TRUE",
               fixed = TRUE)
  expect_error(validate_crps_input(x1, x2, c(1, 2)),
               "ncol(x) == length(y) is not TRUE",
               fixed = TRUE)
  expect_error(validate_crps_input(x1, x2, y, t(ll)),
               "ifelse(is.null(log_lik), TRUE, identical(dim(log_lik), dim(x))) is not TRUE",
               fixed = TRUE)
})

test_that("methods for single data point don't error", {
  expect_silent(suppressWarnings(crps(x1[, 1], x2[, 1], y[1])))
  expect_silent(suppressWarnings(scrps(x1[, 1], x2[, 1], y[1])))
})

# -------------------------------------------------------------------------
# Comparison with measure_rps() / measure_srps() (PWM vs permutation EXX)
# See notes/developer-notes.Rmd ("CRPS / RPS numerical comparison") for details.
# -------------------------------------------------------------------------

.exx_pwm <- function(ypred) {
  n_draws <- nrow(ypred)
  ypred_sorted <- apply(ypred, 2, sort)
  colMeans(ypred_sorted * ((seq_len(n_draws) * (4 / (n_draws - 1))) - 2))
}

.exy_crps <- function(ypred, y) {
  colMeans(abs(sweep(ypred, 2, y)))
}

.crps_draws <- function(seed = 123456789L, n = 10L, S = 100L) {
  set.seed(seed)
  y <- rnorm(n)
  x1 <- matrix(rnorm(n * S), nrow = S)
  x2 <- matrix(rnorm(n * S), nrow = S)
  list(y = y, x1 = x1, x2 = x2)
}

test_that("measure_rps(revert_sign = TRUE) matches deprecated crps() sign convention", {
  d <- .crps_draws()
  old <- suppressWarnings(crps(d$x1, d$x2, d$y))
  new_rev <- measure_rps(d$y, d$x1, revert_sign = TRUE)

  expect_equal(
    as.vector(new_rev$pointwise),
    -as.vector(measure_rps(d$y, d$x1)$pointwise)
  )
  expect_gt(cor(old$pointwise, as.vector(new_rev$pointwise)), 0.98)
  expect_false(isTRUE(all.equal(
    old$pointwise,
    as.vector(new_rev$pointwise),
    tolerance = 1e-6
  )))
})

test_that("crps() and measure_rps() share EXy; differences come from EXX estimator", {
  d <- .crps_draws()
  EXy <- .exy_crps(d$x1, d$y)
  EXX_perm <- suppressWarnings({
    set.seed(1)
    colMeans(abs(d$x1 - d$x2[sample(nrow(d$x1)), , drop = FALSE]))
  })
  EXX_pwm <- .exx_pwm(d$x1)

  expect_equal(EXy, colMeans(abs(sweep(d$x1, 2, d$y))))

  pw_perm <- 0.5 * EXX_perm - EXy
  pw_pwm <- 0.5 * EXX_pwm - EXy
  expect_false(isTRUE(all.equal(EXX_perm, EXX_pwm, tolerance = 1e-6)))
  expect_gt(cor(pw_perm, pw_pwm), 0.98)
  expect_equal(
    max(abs(pw_perm - pw_pwm)),
    max(abs((0.5 * EXX_perm - EXy) - (0.5 * EXX_pwm - EXy)))
  )
})

test_that("scrps() and measure_srps() share sign convention but differ numerically", {
  d <- .crps_draws()
  old <- suppressWarnings(scrps(d$x1, d$x2, d$y))
  new <- measure_srps(d$y, d$x1)

  expect_gt(cor(old$pointwise, as.vector(new$pointwise)), 0.98)
  expect_false(isTRUE(all.equal(
    old$pointwise,
    as.vector(new$pointwise),
    tolerance = 1e-6
  )))
})

test_that("PWM and permutation EXX estimates converge with more draws", {
  set.seed(42)
  n <- 20L
  results <- vapply(
    c(100L, 1000L, 5000L),
    function(S) {
      d <- .crps_draws(seed = 99L, n = n, S = S)
      EXX_perm <- colMeans(abs(d$x1 - d$x2[sample(S), , drop = FALSE]))
      EXX_pwm <- .exx_pwm(d$x1)
      mean(abs(EXX_perm - EXX_pwm) / EXX_pwm)
    },
    numeric(1)
  )
  expect_lt(results[3], results[1])
})

test_that("loo_crps() and loo_pred_measure(..., measure = 'rps') differ like in-sample APIs", {
  set.seed(1)
  d <- .crps_draws(seed = 1L, n = 10L, S = 100L)
  ll <- matrix(rnorm(10 * 100) * 0.1 - 1, nrow = 100)

  old <- suppressWarnings(loo_crps(d$x1, d$x2, d$y, ll))
  psis_obj <- psis(-ll, r_eff = 1)
  new <- loo_pred_measure(
    y = d$y,
    ypred = d$x1,
    ylp = ll,
    measure = "rps",
    psis_object = psis_obj
  )

  expect_gt(
    cor(old$pointwise, -new$pointwise[, "rps_loo"]),
    0.98
  )
  expect_false(isTRUE(all.equal(
    old$pointwise,
    -new$pointwise[, "rps_loo"],
    tolerance = 1e-6
  )))
})

test_that("loo_scrps() and loo_pred_measure(..., measure = 'srps') differ like in-sample APIs", {
  set.seed(1)
  d <- .crps_draws(seed = 1L, n = 10L, S = 100L)
  ll <- matrix(rnorm(10 * 100) * 0.1 - 1, nrow = 100)

  old <- suppressWarnings(loo_scrps(d$x1, d$x2, d$y, ll))
  psis_obj <- psis(-ll, r_eff = 1)
  new <- loo_pred_measure(
    y = d$y,
    ypred = d$x1,
    ylp = ll,
    measure = "srps",
    psis_object = psis_obj
  )

  expect_gt(
    cor(old$pointwise, new$pointwise[, "srps_loo"]),
    0.98
  )
  expect_false(isTRUE(all.equal(
    old$pointwise,
    new$pointwise[, "srps_loo"],
    tolerance = 1e-6
  )))
})
