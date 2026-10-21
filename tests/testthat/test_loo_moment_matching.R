options(mc.cores = 1)

set.seed(123)
S <- 4000

# helper functions for sampling from the posterior distribution
rinvchisq <- function(n, df, scale = 1 / df, ...) {
  if ((length(scale) != 1) & (length(scale) != n)) {
    stop("scale should be a scalar or a vector of the same length as x")
  }
  if (df <= 0) {
    stop("df must be greater than zero")
  }
  if (any(scale <= 0)) {
    stop("scale must be greater than zero")
  }
  return((df * scale) / rchisq(n, df = df))
}

dinvchisq <- function(x, df, scale = 1 / df, log = FALSE, ...) {
  if (df <= 0) {
    stop("df must be greater than zero")
  }
  if (scale <= 0) {
    stop("scale must be greater than zero")
  }
  nu <- df / 2
  if (log) {
    return(ifelse(
      x > 0,
      nu *
        log(nu) -
        log(gamma(nu)) +
        nu * log(scale) -
        (nu + 1) * log(x) -
        (nu * scale / x),
      NA
    ))
  } else {
    return(ifelse(
      x > 0,
      (((nu)^(nu)) / gamma(nu)) *
        (scale^nu) *
        (x^(-(nu + 1))) *
        exp(-nu * scale / x),
      NA
    ))
  }
}


# generate toy data
# normally distributed data with known variance
data_sd <- 1.1
data_mean <- 1.3
n <- as.integer(30)
y <- rnorm(n = n, mean = data_mean, sd = data_sd)
y_tilde <- 11
y[1] <- y_tilde

ymean <- mean(y)
s2 <- sum((y - ymean)^2) / (n - 1)

# draws from the posterior distribution when including all observations
draws_full_posterior_sigma2 <- rinvchisq(S, n - 1, s2)
draws_full_posterior_mu <- rnorm(
  S,
  ymean,
  sqrt(draws_full_posterior_sigma2 / n)
)


# create a dummy model object
x <- list()
x$data <- list()
x$data$y <- y
x$data$n <- n
x$data$ymean <- ymean
x$data$s2 <- s2

x$draws <- data.frame(
  mu = draws_full_posterior_mu,
  sigma = sqrt(draws_full_posterior_sigma2)
)


# implement functions for moment matching loo

# extract original posterior draws
post_draws_test <- function(x, ...) {
  as.matrix(x$draws)
}

# extract original log lik draws
log_lik_i_test <- function(x, i, ...) {
  -0.5 *
    log(2 * pi) -
    log(x$draws$sigma) -
    1.0 / (2 * x$draws$sigma^2) * (x$data$y[i] - x$draws$mu)^2
}


loglik <- matrix(0, S, n)
for (j in seq(n)) {
  loglik[, j] <- log_lik_i_test(x, j)
}


# mu, log(sigma)
unconstrain_pars_test <- function(x, pars, ...) {
  upars <- as.matrix(pars)
  upars[, 2] <- log(upars[, 2])
  upars
}

log_prob_upars_test <- function(x, upars, ...) {
  dinvchisq(exp(upars[, 2])^2, x$data$n - 1, x$data$s2, log = TRUE) +
    dnorm(
      upars[, 1],
      x$data$ymean,
      exp(upars[, 2]) / sqrt(x$data$n),
      log = TRUE
    )
}

# compute log_lik_i values based on the unconstrained parameters
log_lik_i_upars_test <- function(x, upars, i, ...) {
  -0.5 *
    log(2 * pi) -
    upars[, 2] -
    1.0 / (2 * exp(upars[, 2])^2) * (x$data$y[i] - upars[, 1])^2
}


upars <- unconstrain_pars_test(x, x$draws)
lwi_1 <- -loglik[, 1]
lwi_1 <- lwi_1 - matrixStats::logSumExp(lwi_1)


test_that("log_prob_upars_test works", {
  upars <- unconstrain_pars_test(x, x$draws)
  xloo <- list()
  xloo$data <- list()
  xloo$data$y <- y[-1]
  xloo$data$n <- n - 1
  xloo$data$ymean <- mean(y[-1])
  xloo$data$s2 <- sum((y[-1] - mean(y[-1]))^2) / (n - 2)

  post1 <- log_prob_upars_test(x, upars)
  post1 <- post1 - matrixStats::logSumExp(post1)
  post2 <- log_prob_upars_test(xloo, upars) + loglik[, 1]
  post2 <- post2 - matrixStats::logSumExp(post2)
  expect_equal(post1, post2)
})


test_that("loo_moment_match.default warnings work", {
  # loo object
  loo_manual <- suppressWarnings(loo(loglik))
  loo_manual_tis <- suppressWarnings(loo(loglik, is_method = "tis"))

  expect_warning(
    loo_moment_match(
      x,
      loo_manual,
      post_draws_test,
      log_lik_i_test,
      unconstrain_pars_test,
      log_prob_upars_test,
      log_lik_i_upars_test,
      max_iters = 30L,
      k_thres = 0.5,
      split = FALSE,
      cov = TRUE,
      cores = 1
    ),
    "The accuracy of self-normalized importance sampling"
  )

  expect_warning(
    loo_moment_match(
      x,
      loo_manual,
      post_draws_test,
      log_lik_i_test,
      unconstrain_pars_test,
      log_prob_upars_test,
      log_lik_i_upars_test,
      max_iters = 30L,
      split = FALSE,
      cov = TRUE,
      cores = 1
    ),
    "The accuracy of self-normalized importance sampling"
  )

  expect_no_warning(loo_moment_match(
    x,
    loo_manual,
    post_draws_test,
    log_lik_i_test,
    unconstrain_pars_test,
    log_prob_upars_test,
    log_lik_i_upars_test,
    max_iters = 30L,
    k_thres = 100,
    split = TRUE,
    cov = TRUE,
    cores = 1
  ))

  expect_snapshot(loo_moment_match(
    x,
    loo_manual,
    post_draws_test,
    log_lik_i_test,
    unconstrain_pars_test,
    log_prob_upars_test,
    log_lik_i_upars_test,
    max_iters = 1,
    k_thres = 0.5,
    split = TRUE,
    cov = TRUE,
    cores = 1
  ))

  expect_error(
    loo_moment_match(
      x,
      loo_manual_tis,
      post_draws_test,
      log_lik_i_test,
      unconstrain_pars_test,
      log_prob_upars_test,
      log_lik_i_upars_test,
      max_iters = 30L,
      k_thres = 0.5,
      split = TRUE,
      cov = TRUE,
      cores = 1
    ),
    "loo_moment_match currently supports only"
  )
})


test_that("loo_moment_match.default works", {
  # allow -Inf
  lwi_x <- lwi_1
  lwi_x[which.min(lwi_1)] <- -Inf
  expect_no_error(suppressWarnings(importance_sampling.default(
    lwi_x,
    method = "psis",
    r_eff = 1,
    cores = 1
  )))

  # loo object
  loo_manual <- suppressWarnings(loo(loglik))

  loo_moment_match_object <- suppressWarnings(loo_moment_match(
    x,
    loo_manual,
    post_draws_test,
    log_lik_i_test,
    unconstrain_pars_test,
    log_prob_upars_test,
    log_lik_i_upars_test,
    max_iters = 30L,
    k_thres = 0.8,
    split = FALSE,
    cov = TRUE,
    cores = 1
  ))

  # diagnostic Pareto k decreases but influence pareto k stays the same
  expect_lt(
    loo_moment_match_object$diagnostics$pareto_k[1],
    loo_moment_match_object$pointwise[1, "influence_pareto_k"]
  )
  expect_equal(
    loo_moment_match_object$pointwise[, "influence_pareto_k"],
    loo_manual$pointwise[, "influence_pareto_k"]
  )
  expect_equal(
    loo_moment_match_object$pointwise[, "influence_pareto_k"],
    loo_manual$diagnostics$pareto_k
  )
  expect_equal(
    unname(loo_moment_match_object$pointwise[1, "mcse_elpd_loo"]),
    0.1904162,
    tolerance = 1e-6
  )

  expect_snapshot_value(loo_moment_match_object, style = "serialize")

  loo_moment_match_object2 <- suppressWarnings(loo_moment_match(
    x,
    loo_manual,
    post_draws_test,
    log_lik_i_test,
    unconstrain_pars_test,
    log_prob_upars_test,
    log_lik_i_upars_test,
    max_iters = 30L,
    k_thres = 0.5,
    split = FALSE,
    cov = TRUE,
    cores = 1
  ))

  expect_snapshot_value(loo_moment_match_object2, style = "serialize")

  loo_moment_match_object3 <- suppressWarnings(loo_moment_match(
    x,
    loo_manual,
    post_draws_test,
    log_lik_i_test,
    unconstrain_pars_test,
    log_prob_upars_test,
    log_lik_i_upars_test,
    max_iters = 30L,
    k_thres = 0.5,
    split = TRUE,
    cov = TRUE,
    cores = 1
  ))

  expect_snapshot_value(loo_moment_match_object3, style = "serialize")

  loo_moment_match_object4 <- suppressWarnings(loo_moment_match(
    x,
    loo_manual,
    post_draws_test,
    log_lik_i_test,
    unconstrain_pars_test,
    log_prob_upars_test,
    log_lik_i_upars_test,
    max_iters = 30L,
    k_thres = 100,
    split = FALSE,
    cov = TRUE,
    cores = 1
  ))

  expect_equal(loo_manual, loo_moment_match_object4)

  loo_manual_with_psis <- suppressWarnings(loo(loglik, save_psis = TRUE))
  loo_moment_match_object5 <- suppressWarnings(loo_moment_match(
    x,
    loo_manual_with_psis,
    post_draws_test,
    log_lik_i_test,
    unconstrain_pars_test,
    log_prob_upars_test,
    log_lik_i_upars_test,
    max_iters = 30L,
    k_thres = 0.8,
    split = FALSE,
    cov = TRUE,
    cores = 1
  ))

  expect_equal(
    loo_moment_match_object5$diagnostics,
    loo_moment_match_object5$psis_object$diagnostics
  )
})

test_that("variance and covariance transformations work", {
  S <- 2000

  set.seed(8493874)
  draws_full_posterior_sigma2 <- rinvchisq(S, n - 1, s2)
  draws_full_posterior_mu <- rnorm(
    S,
    ymean,
    sqrt(draws_full_posterior_sigma2 / n)
  )

  x$draws <- data.frame(
    mu = draws_full_posterior_mu,
    sigma = sqrt(draws_full_posterior_sigma2)
  )
  loglik <- matrix(0, S, n)
  for (j in seq(n)) {
    loglik[, j] <- log_lik_i_test(x, j)
  }

  upars <- unconstrain_pars_test(x, x$draws)
  lwi_1 <- -loglik[, 1]
  lwi_1 <- lwi_1 - matrixStats::logSumExp(lwi_1)

  loo_manual <- suppressWarnings(loo(loglik))

  loo_moment_match_object <- suppressWarnings(loo_moment_match(
    x,
    loo_manual,
    post_draws_test,
    log_lik_i_test,
    unconstrain_pars_test,
    log_prob_upars_test,
    log_lik_i_upars_test,
    max_iters = 30L,
    k_thres = 0.0,
    split = FALSE,
    cov = TRUE,
    cores = 1
  ))

  expect_snapshot_value(loo_moment_match_object, style = "serialize")
})

test_that("loo_moment_match.default works with multiple cores", {
  # forked processes can segfault on CRAN's macOS machines if R is linked to
  # Apple's Accelerate/vecLib BLAS, which is not fork-safe
  skip_on_cran()
  # loo object
  loo_manual <- suppressWarnings(loo(loglik))

  loo_moment_match_manual3 <- suppressWarnings(loo_moment_match(
    x,
    loo_manual,
    post_draws_test,
    log_lik_i_test,
    unconstrain_pars_test,
    log_prob_upars_test,
    log_lik_i_upars_test,
    max_iters = 30L,
    k_thres = 0.5,
    split = FALSE,
    cov = TRUE,
    cores = 1
  ))

  loo_moment_match_manual4 <- suppressWarnings(loo_moment_match(
    x,
    loo_manual,
    post_draws_test,
    log_lik_i_test,
    unconstrain_pars_test,
    log_prob_upars_test,
    log_lik_i_upars_test,
    max_iters = 30L,
    k_thres = 0.5,
    split = FALSE,
    cov = TRUE,
    cores = 2
  ))

  expect_equal(
    loo_moment_match_manual3$diagnostics$pareto_k,
    loo_moment_match_manual4$diagnostics$pareto_k
  )
  expect_equal(
    loo_moment_match_manual3$diagnostics$n_eff,
    loo_moment_match_manual4$diagnostics$n_eff,
    tolerance = 5e-4
  )

  expect_equal(
    loo_moment_match_manual3$estimates,
    loo_moment_match_manual4$estimates
  )

  expect_equal(
    loo_moment_match_manual3$pointwise,
    loo_moment_match_manual4$pointwise,
    tolerance = 5e-4
  )
})


test_that("loo_moment_match_split works", {
  # skip on M1 Mac until we figure out why this test fails only on M1 Mac
  skip_if(Sys.info()[["sysname"]] == "Darwin" && R.version$arch == "aarch64")

  is_obj_1 <- suppressWarnings(importance_sampling.default(
    lwi_1,
    method = "psis",
    r_eff = 1,
    cores = 1
  ))
  lwi_1_ps <- as.vector(weights(is_obj_1))

  split <- loo_moment_match_split(
    x,
    upars,
    cov = FALSE,
    total_shift = c(0, 0),
    total_scaling = c(1, 1),
    total_mapping = diag(c(1, 1)),
    i = 1,
    log_prob_upars = log_prob_upars_test,
    log_lik_i_upars = log_lik_i_upars_test,
    cores = 1,
    r_eff_i = 1,
    is_method = "psis"
  )

  expect_named(split, c("lwi", "lwfi", "log_liki", "r_eff_i"))

  expect_equal(lwi_1_ps, split$lwi)

  # with cov = TRUE, loo_moment_match() includes the scaling in the mapping
  # and leaves total_scaling at one
  split2 <- loo_moment_match_split(
    x,
    upars,
    cov = TRUE,
    total_shift = c(-0.1, -0.2),
    total_scaling = c(1, 1),
    total_mapping = 0.7 * matrix(c(1, 0.1, 0.1, 1), 2, 2),
    i = 1,
    log_prob_upars = log_prob_upars_test,
    log_lik_i_upars = log_lik_i_upars_test,
    cores = 1,
    r_eff_i = 1,
    is_method = "psis"
  )

  expect_snapshot_value(split2, style = "serialize")

  log_prob_with_inf <- function(x, upars, ...) {
    out <- log_prob_upars_test(x, upars, ...)
    out[1] <- -Inf
    out
  }
  expect_no_error(split3 <- loo_moment_match_split(
    x,
    upars,
    cov = FALSE,
    total_shift = c(0, 0),
    total_scaling = c(1, 1),
    total_mapping = diag(c(1, 1)),
    i = 1,
    log_prob_upars = log_prob_with_inf,
    log_lik_i_upars = log_lik_i_upars_test,
    cores = 1,
    r_eff_i = 1,
    is_method = "psis"
  ))
  expect_false(anyNA(split3$lwi))
})

test_that("loo_moment_match_split evaluates the inverse of its transformation", {
  # The second half of the draws is untransformed. For these draws theta,
  # the mixture proposal needs the density at T^-1(theta), where
  # T(theta) = M diag(s) (theta - m) + m + shift.
  shift <- c(-0.1, -0.2)
  scaling <- c(0.7, 1.3)
  m <- colMeans(upars)
  second <- seq_len(S)[-seq_len(as.integer(0.5 * S))]
  for (cov in c(FALSE, TRUE)) {
    mapping <- if (cov) matrix(c(1, 0.1, 0.1, 1), 2, 2) else diag(2)
    seen <- list()
    log_prob_record <- function(x, upars, ...) {
      seen[[length(seen) + 1]] <<- upars
      log_prob_upars_test(x, upars, ...)
    }
    loo_moment_match_split(
      x,
      upars,
      cov = cov,
      total_shift = shift,
      total_scaling = scaling,
      total_mapping = mapping,
      i = 1,
      log_prob_upars = log_prob_record,
      log_lik_i_upars = log_lik_i_upars_test,
      cores = 1,
      r_eff_i = 1,
      is_method = "psis"
    )
    expect_length(seen, 2)
    first <- seq_len(as.integer(0.5 * S))
    # first call: second half untransformed; second call: first half
    # untransformed
    expect_equal(seen[[1]][second, , drop = FALSE], upars[second, , drop = FALSE])
    expect_equal(seen[[2]][first, , drop = FALSE], upars[first, , drop = FALSE])
    # the second call evaluates the inverse-transformed draws
    inv <- seen[[2]][second, , drop = FALSE]
    fwd <- sweep(inv, 2, m, "-")
    fwd <- sweep(fwd, 2, scaling, "*")
    if (cov) {
      fwd <- tcrossprod(fwd, mapping)
    }
    fwd <- sweep(fwd, 2, shift + m, "+")
    expect_equal(unname(fwd), unname(upars[second, , drop = FALSE]))
  }
})

test_that("the total transformation composes the accepted steps in order", {
  # The moment matching loop starts again with the shift step after every
  # accepted step, so a scaling step can follow a covariance step. The
  # transformation that loo_moment_match_split() rebuilds from the totals
  # must map the original draws to the draws that the steps produced one
  # after another.
  first <- seq_len(as.integer(0.5 * S))
  sequences <- list(
    list(cov = TRUE, steps = c("cov", "scale", "shift", "cov", "scale")),
    list(cov = FALSE, steps = c("scale", "shift", "scale"))
  )
  for (sq in sequences) {
    total <- list(shift = c(0, 0), scaling = c(1, 1), mapping = diag(2))
    u <- upars
    for (step in sq$steps) {
      trans <- switch(step,
        shift = shift(x, u, lwi_1),
        scale = shift_and_scale(x, u, lwi_1),
        cov = shift_and_cov(x, u, lwi_1)
      )
      u <- trans$upars
      total <- add_step_to_total(total, trans, sq$cov)
    }
    seen <- list()
    log_prob_record <- function(x, upars, ...) {
      seen[[length(seen) + 1]] <<- upars
      log_prob_upars_test(x, upars, ...)
    }
    loo_moment_match_split(
      x,
      upars,
      cov = sq$cov,
      total_shift = total$shift,
      total_scaling = total$scaling,
      total_mapping = total$mapping,
      i = 1,
      log_prob_upars = log_prob_record,
      log_lik_i_upars = log_lik_i_upars_test,
      cores = 1,
      r_eff_i = 1,
      is_method = "psis"
    )
    # the first call evaluates T(theta) for the first half of the draws
    expect_equal(
      unname(seen[[1]][first, , drop = FALSE]),
      unname(u[first, , drop = FALSE])
    )
  }
})

test_that("shift_and_scale is stable for large parameter values", {
  upars <- matrix(1e8 + c(-1, 0, 1), ncol = 1)
  result <- shift_and_scale(NULL, upars, rep(-log(3), 3))

  expect_equal(result$scaling, 1)
  expect_equal(result$upars, upars)
})

test_that("passing arguments works", {
  log_lik_i_upars_test_additional_argument <- function(
    x,
    upars,
    i,
    passed_arg = FALSE,
    ...
  ) {
    if (!passed_arg) {
      warning("passed_arg was not passed here")
    }
    -0.5 *
      log(2 * pi) -
      upars[, 2] -
      1.0 / (2 * exp(upars[, 2])^2) * (x$data$y[i] - upars[, 1])^2
  }
  unconstrain_pars_test_additional_argument <- function(
    x,
    pars,
    passed_arg = FALSE,
    ...
  ) {
    if (!passed_arg) {
      warning("passed_arg was not passed here")
    }
    upars <- as.matrix(pars)
    upars[, 2] <- log(upars[, 2])
    upars
  }

  log_prob_upars_test_additional_argument <- function(
    x,
    upars,
    passed_arg = FALSE,
    ...
  ) {
    if (!passed_arg) {
      warning("passed_arg was not passed here")
    }
    dinvchisq(exp(upars[, 2])^2, x$data$n - 1, x$data$s2, log = TRUE) +
      dnorm(
        upars[, 1],
        x$data$ymean,
        exp(upars[, 2]) / sqrt(x$data$n),
        log = TRUE
      )
  }
  post_draws_test_additional_argument <- function(x, passed_arg = FALSE, ...) {
    if (!passed_arg) {
      warning("passed_arg was not passed here")
    }
    as.matrix(x$draws)
  }
  log_lik_i_test_additional_argument <- function(
    x,
    i,
    passed_arg = FALSE,
    ...
  ) {
    if (!passed_arg) {
      warning("passed_arg was not passed here")
    }
    -0.5 *
      log(2 * pi) -
      log(x$draws$sigma) -
      1.0 / (2 * x$draws$sigma^2) * (x$data$y[i] - x$draws$mu)^2
  }

  # loo object
  loo_manual <- suppressWarnings(loo(loglik))
  expect_silent(loo_moment_match(
    x,
    loo_manual,
    post_draws_test_additional_argument,
    log_lik_i_test_additional_argument,
    unconstrain_pars_test_additional_argument,
    log_prob_upars_test_additional_argument,
    log_lik_i_upars_test_additional_argument,
    max_iters = 30L,
    k_thres = 0.5,
    split = TRUE,
    cov = TRUE,
    cores = 1,
    passed_arg = TRUE
  ))
})

# log_lik_i_upars() that fails in the calls listed in `fail_calls`, as a model
# can fail for the extreme draws of a transformation
failing_log_lik_i_upars <- function(fail_calls) {
  calls <- 0
  function(x, upars, i, ...) {
    calls <<- calls + 1
    if (calls %in% fail_calls) {
      stop("overflow in log_lik_i_upars")
    }
    log_lik_i_upars_test(x, upars, i)
  }
}

mm_with <- function(loo_x, log_lik_i = log_lik_i_test,
                    log_lik_i_upars = log_lik_i_upars_test, split = TRUE,
                    k_threshold = 0.7, cores = 1) {
  loo_moment_match(
    x, loo = loo_x, post_draws = post_draws_test, log_lik_i = log_lik_i,
    unconstrain_pars = unconstrain_pars_test,
    log_prob_upars = log_prob_upars_test, log_lik_i_upars = log_lik_i_upars,
    k_threshold = k_threshold, split = split, cores = cores
  )
}

test_that("an error in the split step keeps the result without the split", {
  loo_manual <- suppressWarnings(loo(loglik))
  # observation 1 needs two transformations; the third call is the split step
  warnings <- character(0)
  mm <- withCallingHandlers(
    mm_with(loo_manual, log_lik_i_upars = failing_log_lik_i_upars(3)),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_match(warnings, "stopped early for observation(s) 1", fixed = TRUE,
               all = FALSE)
  expect_match(warnings, "The split transformation failed", all = FALSE)
  mm_no_split <- suppressWarnings(mm_with(loo_manual, split = FALSE))
  expect_equal(mm$pointwise, mm_no_split$pointwise)
  expect_equal(mm$diagnostics, mm_no_split$diagnostics)
})

test_that("the accuracy warning names a failed split transformation", {
  expect_warning(throw_large_kf_warning(0.8, 0.7), "Setting the argument 'split'")
  expect_warning(
    throw_large_kf_warning(0.8, 0.7, split_failed = TRUE),
    "The split transformation failed"
  )
  expect_no_warning(throw_large_kf_warning(0.5, 0.7, split_failed = TRUE))
})

test_that("an error in the first transformation keeps the PSIS results", {
  loo_manual <- suppressWarnings(loo(loglik))
  # k_threshold = 0 moment-matches all observations; the first call of
  # log_lik_i_upars is the first transformation of observation 1
  warnings <- character(0)
  mm <- withCallingHandlers(
    mm_with(loo_manual, log_lik_i_upars = failing_log_lik_i_upars(1),
            k_threshold = 0),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  # no transformation was accepted, so the observation counts as failed
  expect_match(warnings, "failed for observation(s) 1.", fixed = TRUE,
               all = FALSE)
  expect_match(warnings, "overflow in log_lik_i_upars", all = FALSE)
  expect_equal(mm$pointwise[1, ], loo_manual$pointwise[1, ])
  expect_equal(mm$diagnostics$pareto_k[1], loo_manual$diagnostics$pareto_k[1])
})

test_that("moment matching stops if it fails for all observations", {
  loo_manual <- suppressWarnings(loo(loglik))
  # only observation 1 is above k_threshold = 0.7
  expect_error(
    suppressWarnings(mm_with(loo_manual,
                             log_lik_i_upars = failing_log_lik_i_upars(1))),
    "Moment matching failed for all observations.*overflow in log_lik_i_upars"
  )
  # a log_lik_i with the wrong arguments fails for every observation
  wrong_log_lik_i <- function(x) x
  expect_error(
    suppressWarnings(mm_with(loo_manual, log_lik_i = wrong_log_lik_i,
                             k_threshold = 0)),
    "Moment matching failed for all observations"
  )
})

test_that("errors in loo itself are not reported as failed observations", {
  loo_manual <- suppressWarnings(loo(loglik))
  # log_lik_i returns without an error, and loo's own code then fails
  non_numeric_log_lik_i <- function(x, i, ...) {
    if (i == 1) "not numeric" else log_lik_i_test(x, i)
  }
  cores_list <- if (os_is_windows()) 1 else c(1, 2)
  for (cores in cores_list) {
    expect_error(
      # with cores = 2, mclapply() also warns that a worker had an error
      suppressWarnings(mm_with(loo_manual, log_lik_i = non_numeric_log_lik_i,
                               k_threshold = 0, cores = cores)),
      "non-numeric argument to mathematical function"
    )
  }
})

test_that("a parallel worker that dies gives failed observations", {
  skip_on_cran()
  skip_on_os("windows")  # the worker is a forked process
  loo_manual <- suppressWarnings(loo(loglik))
  parent <- Sys.getpid()
  dying_log_lik_i <- function(x, i, ...) {
    if (i == 1 && Sys.getpid() != parent) {
      tools::pskill(Sys.getpid(), tools::SIGKILL)
    }
    log_lik_i_test(x, i)
  }
  warnings <- character(0)
  mm <- withCallingHandlers(
    mm_with(loo_manual, log_lik_i = dying_log_lik_i, k_threshold = 0,
            cores = 2),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_match(warnings, "a parallel worker did not return a result",
               all = FALSE)
  expect_equal(mm$pointwise[1, ], loo_manual$pointwise[1, ])
})

test_that("an error for one observation keeps its PSIS results", {
  skip_on_cran()
  loo_manual <- suppressWarnings(loo(loglik))
  failing_log_lik_i <- function(x, i, ...) {
    if (i == 1) stop("log_lik_i failed")
    log_lik_i_test(x, i)
  }
  mm_ok <- suppressWarnings(mm_with(loo_manual, k_threshold = 0))
  # On Windows the workers are PSOCK processes that load the installed loo, so
  # with devtools::test() they would not run this version of the code
  cores_list <- if (os_is_windows()) 1 else c(1, 2)
  for (cores in cores_list) {
    # observation 1 keeps its high Pareto k, which gives a second warning
    suppressWarnings(expect_warning(
      mm <- mm_with(loo_manual, log_lik_i = failing_log_lik_i,
                    k_threshold = 0, cores = cores),
      "failed for observation(s) 1.",
      fixed = TRUE
    ))
    expect_equal(mm$pointwise[1, ], loo_manual$pointwise[1, ])
    expect_equal(mm$pointwise[-1, ], mm_ok$pointwise[-1, ])
  }
})
