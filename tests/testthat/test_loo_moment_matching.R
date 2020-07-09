library(loo)
options(mc.cores = 1)


context("moment matching")
set.seed(123)
S <- 4000

# helper functions for sampling from the posterior distribution
rinvchisq <- function(n, df, scale = 1/df, ...)
{
  if ((length(scale) != 1) & (length(scale) != n))
    stop("scale should be a scalar or a vector of the same length as x")
  if (df <= 0)
    stop("df must be greater than zero")
  if (any(scale <= 0))
    stop("scale must be greater than zero")
  return((df*scale)/rchisq(n, df = df))
}

dinvchisq <- function(x, df, scale=1/df, log = FALSE, ...)
{
  if (df <= 0)
    stop("df must be greater than zero")
  if (scale <= 0)
    stop("scale must be greater than zero")
  nu <- df/2
  if (log)
    return(ifelse(x > 0, nu*log(nu) - log(gamma(nu)) + nu*log(scale) -
                    (nu + 1)*log(x) - (nu*scale/x), NA))
  else
    return(ifelse(x > 0,
                  (((nu)^(nu))/gamma(nu)) * (scale^nu) *
                    (x^(-(nu + 1))) * exp(-nu*scale/x), NA))
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
s2 <- sum((y - ymean)^2)/(n - 1)

# draws from the posterior distribution when including all observations
draws_full_posterior_sigma2 <- rinvchisq(S, n - 1, s2)
draws_full_posterior_mu <- rnorm(S, ymean, sqrt(draws_full_posterior_sigma2/n))


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
  -0.5*log(2*pi) - log(x$draws$sigma) - 1.0/(2*x$draws$sigma^2)*(x$data$y[i] - x$draws$mu)^2
}


loglik <- matrix(0,S,n)
for (j in seq(n)) {
  loglik[,j] <- log_lik_i_test(x, j)
}



# mu, log(sigma)
unconstrain_pars_test <- function(x, pars, ...) {
  upars <- as.matrix(pars)
  upars[,2] <- log(upars[,2])
  upars
}

log_prob_upars_test <- function(x, upars, ...) {
  dinvchisq(exp(upars[,2])^2,x$data$n - 1,x$data$s2, log = TRUE) +
    dnorm(upars[,1],x$data$ymean,exp(upars[,2])/sqrt(x$data$n), log = TRUE)
}

# compute log_lik_i values based on the unconstrained parameters
log_lik_i_upars_test <- function(x, upars, i, ...) {
  -0.5*log(2*pi) - upars[,2] - 1.0/(2*exp(upars[,2])^2)*(x$data$y[i] - upars[,1])^2
}



upars <- unconstrain_pars_test(x, x$draws)
lwi_1 <- -loglik[,1]
lwi_1 <- lwi_1 - matrixStats::logSumExp(lwi_1)



test_that("log_prob_upars_test works", {
  upars <- unconstrain_pars_test(x, x$draws)
  xloo <- list()
  xloo$data <- list()
  xloo$data$y <- y[-1]
  xloo$data$n <- n - 1
  xloo$data$ymean <- mean(y[-1])
  xloo$data$s2 <- sum((y[-1] - mean(y[-1]))^2)/(n - 2)

  post1 <- log_prob_upars_test(x,upars)
  post1 <- post1 - matrixStats::logSumExp(post1)
  post2 <- log_prob_upars_test(xloo,upars) + loglik[,1]
  post2 <- post2 - matrixStats::logSumExp(post2)
  expect_equal(post1,post2)
})


test_that("loo_moment_match.default warnings work", {
  # loo object
  loo_manual <- suppressWarnings(loo(loglik))
  loo_manual_tis <- suppressWarnings(loo(loglik, is_method = "tis"))


  expect_warning(loo_moment_match(x, loo_manual, post_draws_test, log_lik_i_test,
                              unconstrain_pars_test, log_prob_upars_test,
                              log_lik_i_upars_test, max_iters = 30L,
                              k_thres = 100, split = FALSE,
                              cov = TRUE, cores = 1), "Some Pareto k")

  expect_warning(loo_moment_match(x, loo_manual, post_draws_test, log_lik_i_test,
                              unconstrain_pars_test, log_prob_upars_test,
                              log_lik_i_upars_test, max_iters = 30L,
                              k_thres = 0.5, split = FALSE,
                              cov = TRUE, cores = 1), "The accuracy of self-normalized importance sampling")

  expect_warning(loo_moment_match(x, loo_manual, post_draws_test, log_lik_i_test,
                              unconstrain_pars_test, log_prob_upars_test,
                              log_lik_i_upars_test, max_iters = 1,
                              k_thres = 0.5, split = TRUE,
                              cov = TRUE, cores = 1), "The maximum number of moment matching iterations")

  expect_error(loo_moment_match(x, loo_manual_tis, post_draws_test, log_lik_i_test,
                       unconstrain_pars_test, log_prob_upars_test,
                       log_lik_i_upars_test, max_iters = 30L,
                       k_thres = 0.5, split = TRUE,
                       cov = TRUE, cores = 1), "loo_moment_match currently supports only")
})




test_that("loo_moment_match.default works", {

  # loo object
  loo_manual <- suppressWarnings(loo(loglik))

  loo_moment_match_object <- suppressWarnings(loo_moment_match(x, loo_manual, post_draws_test, log_lik_i_test,
                                                unconstrain_pars_test, log_prob_upars_test,
                                                log_lik_i_upars_test, max_iters = 30L,
                                                k_thres = 0.8, split = FALSE,
                                                cov = TRUE, cores = 1))

  # diagnostic pareto k decreases but influence pareto k stays the same
  expect_lt(loo_moment_match_object$diagnostics$pareto_k[1], loo_moment_match_object$pointwise[1,"influence_pareto_k"])
  expect_equal(loo_moment_match_object$pointwise[,"influence_pareto_k"],loo_manual$pointwise[,"influence_pareto_k"])
  expect_equal(loo_moment_match_object$pointwise[,"influence_pareto_k"],loo_manual$diagnostics$pareto_k)

  expect_equal_to_reference(loo_moment_match_object, "reference-results/moment_match_loo_1.rds")

  loo_moment_match_object2 <- suppressWarnings(loo_moment_match(x, loo_manual, post_draws_test, log_lik_i_test,
                                                unconstrain_pars_test, log_prob_upars_test,
                                                log_lik_i_upars_test, max_iters = 30L,
                                                k_thres = 0.5, split = FALSE,
                                                cov = TRUE, cores = 1))

  expect_equal_to_reference(loo_moment_match_object2, "reference-results/moment_match_loo_2.rds")

  loo_moment_match_object3 <- suppressWarnings(loo_moment_match(x, loo_manual, post_draws_test, log_lik_i_test,
                                                unconstrain_pars_test, log_prob_upars_test,
                                                log_lik_i_upars_test, max_iters = 30L,
                                                k_thres = 0.5, split = TRUE,
                                                cov = TRUE, cores = 1))

  expect_equal_to_reference(loo_moment_match_object3, "reference-results/moment_match_loo_3.rds")

  loo_moment_match_object4 <- suppressWarnings(loo_moment_match(x, loo_manual, post_draws_test, log_lik_i_test,
                                                unconstrain_pars_test, log_prob_upars_test,
                                                log_lik_i_upars_test, max_iters = 30L,
                                                k_thres = 100, split = FALSE,
                                                cov = TRUE, cores = 1))

  expect_equal(loo_manual,loo_moment_match_object4)

  loo_manual_with_psis <- suppressWarnings(loo(loglik, save_psis = TRUE))
  loo_moment_match_object5 <- suppressWarnings(loo_moment_match(x, loo_manual_with_psis, post_draws_test, log_lik_i_test,
                                          unconstrain_pars_test, log_prob_upars_test,
                                          log_lik_i_upars_test, max_iters = 30L,
                                          k_thres = 0.8, split = FALSE,
                                          cov = TRUE, cores = 1))

  expect_equal(loo_moment_match_object5$diagnostics,loo_moment_match_object5$psis_object$diagnostics)


})

test_that("variance and covariance transformations work", {
  S <- 2000

  set.seed(8493874)
  draws_full_posterior_sigma2 <- rinvchisq(S, n - 1, s2)
  draws_full_posterior_mu <- rnorm(S, ymean, sqrt(draws_full_posterior_sigma2/n))

  x$draws <- data.frame(
    mu = draws_full_posterior_mu,
    sigma = sqrt(draws_full_posterior_sigma2)
  )
  loglik <- matrix(0,S,n)
  for (j in seq(n)) {
    loglik[,j] <- log_lik_i_test(x, j)
  }

  upars <- unconstrain_pars_test(x, x$draws)
  lwi_1 <- -loglik[,1]
  lwi_1 <- lwi_1 - matrixStats::logSumExp(lwi_1)

  loo_manual <- suppressWarnings(loo(loglik))

  loo_moment_match_object <- suppressWarnings(loo_moment_match(x, loo_manual, post_draws_test, log_lik_i_test,
                                         unconstrain_pars_test, log_prob_upars_test,
                                         log_lik_i_upars_test, max_iters = 30L,
                                         k_thres = 0.0, split = FALSE,
                                         cov = TRUE, cores = 1))

  expect_equal_to_reference(loo_moment_match_object, "reference-results/moment_match_var_and_cov.rds")

})

test_that("loo_moment_match.default works with multiple cores", {

  # loo object
  loo_manual <- suppressWarnings(loo(loglik))

  loo_moment_match_manual3 <- suppressWarnings(loo_moment_match(x, loo_manual, post_draws_test, log_lik_i_test,
                                                 unconstrain_pars_test, log_prob_upars_test,
                                                 log_lik_i_upars_test, max_iters = 30L,
                                                 k_thres = 0.5, split = FALSE,
                                                 cov = TRUE, cores = 1))

  loo_moment_match_manual4 <- suppressWarnings(loo_moment_match(x, loo_manual, post_draws_test, log_lik_i_test,
                                                 unconstrain_pars_test, log_prob_upars_test,
                                                 log_lik_i_upars_test, max_iters = 30L,
                                                 k_thres = 0.5, split = FALSE,
                                                 cov = TRUE, cores = 2))

  expect_equal(loo_moment_match_manual3$diagnostics$pareto_k, loo_moment_match_manual4$diagnostics$pareto_k)
  expect_equal(loo_moment_match_manual3$diagnostics$n_eff, loo_moment_match_manual4$diagnostics$n_eff)

  expect_equal(loo_moment_match_manual3$estimates, loo_moment_match_manual4$estimates)

  expect_equal(loo_moment_match_manual3$pointwise, loo_moment_match_manual4$pointwise, tolerance=5e-4)

})



test_that("loo_moment_match_split works", {

  is_obj_1 <- suppressWarnings(importance_sampling.default(lwi_1, method = "psis", r_eff = 1, cores = 1))
  lwi_1_ps <- as.vector(weights(is_obj_1))

  split <- loo_moment_match_split(
    x, upars, cov = FALSE, total_shift = c(0,0), total_scaling = c(1,1), total_mapping = diag(c(1,1)), i = 1,
    log_prob_upars = log_prob_upars_test, log_lik_i_upars = log_lik_i_upars_test,
    cores = 1, r_eff_i = 1, is_method = "psis")

  expect_named(split,c("lwi", "lwfi", "log_liki", "r_eff_i"))

  expect_equal(lwi_1_ps,split$lwi)

  split2 <- loo_moment_match_split(
    x, upars, cov = FALSE, total_shift = c(-0.1,-0.2), total_scaling = c(0.7,0.7),
    total_mapping = matrix(c(1,0.1,0.1,1),2,2), i = 1,
    log_prob_upars = log_prob_upars_test, log_lik_i_upars = log_lik_i_upars_test,
    cores = 1, r_eff_i = 1, is_method = "psis")

  expect_equal_to_reference(split2, "reference-results/moment_match_split.rds")

})

test_that("passing arguments works", {
  log_lik_i_upars_test_additional_argument <- function(x, upars, i, passed_arg = FALSE, ...) {
    if (!passed_arg) {
      warning("passed_arg was not passed here")
    }
    -0.5*log(2*pi) - upars[,2] - 1.0/(2*exp(upars[,2])^2)*(x$data$y[i] - upars[,1])^2

  }
  unconstrain_pars_test_additional_argument <- function(x, pars, passed_arg = FALSE, ...) {
    if (!passed_arg) {
      warning("passed_arg was not passed here")
    }
    upars <- as.matrix(pars)
    upars[,2] <- log(upars[,2])
    upars
  }

  log_prob_upars_test_additional_argument <- function(x, upars, passed_arg = FALSE, ...) {
    if (!passed_arg) {
      warning("passed_arg was not passed here")
    }
    dinvchisq(exp(upars[,2])^2,x$data$n - 1,x$data$s2, log = TRUE) +
      dnorm(upars[,1],x$data$ymean,exp(upars[,2])/sqrt(x$data$n), log = TRUE)
  }
  post_draws_test_additional_argument <- function(x, passed_arg = FALSE, ...) {
    if (!passed_arg) {
      warning("passed_arg was not passed here")
    }
    as.matrix(x$draws)
  }
  log_lik_i_test_additional_argument <- function(x, i, passed_arg = FALSE, ...) {
    if (!passed_arg) {
      warning("passed_arg was not passed here")
    }
    -0.5*log(2*pi) - log(x$draws$sigma) - 1.0/(2*x$draws$sigma^2)*(x$data$y[i] - x$draws$mu)^2
  }

  # loo object
  loo_manual <- suppressWarnings(loo(loglik))
  expect_silent(loo_moment_match(x, loo_manual, post_draws_test_additional_argument, log_lik_i_test_additional_argument,
                                 unconstrain_pars_test_additional_argument, log_prob_upars_test_additional_argument,
                                 log_lik_i_upars_test_additional_argument, max_iters = 30L,
                                 k_thres = 0.5, split = TRUE,
                                 cov = TRUE, cores = 1, passed_arg = TRUE))
})
