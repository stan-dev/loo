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

test_that("crps matches references", {
  expect_equal_to_reference(with_seed(1, crps(x1, x2, y)), 'reference-results/crps.rds')
  expect_equal_to_reference(with_seed(1, scrps(x1, x2, y)), 'reference-results/scrps.rds')
  # Suppress warnings for the missing r_eff
  suppressWarnings(expect_equal_to_reference(
    with_seed(1, loo_crps(x1, x2, y, ll)),
    'reference-results/loo_crps.rds'))
  suppressWarnings(expect_equal_to_reference(
    with_seed(1, loo_scrps(x1, x2, y, ll)),
    'reference-results/loo_scrps.rds'))
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
  expect_silent(crps(x1[,1], x2[,1], y[1]))
  expect_silent(scrps(x1[,1], x2[,1], y[1]))
})
