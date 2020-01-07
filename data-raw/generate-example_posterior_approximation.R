# This file create data to test the loo using approximative posteriors
# Three datasets are used:
# 1. Linear regression with normal (independent) posterior - Laplace, ADVI-mf, ADVI-fr should all pass
# 2. Linear regression with normal (correlated) posterior - Laplace, ADVI-fr pass, ADVI-mf fail
# 3. Normal model unknown mu and sigma - Laplace, ADVI-fr, ADVI-mf all should fail

# Function used to generate data

#' Generate simulated data for a linear regression model with a correlated posterior
#' @param D dimension size
#' @param N Number of observations
#' @param off_diag_cov Off diagonal element covariance
#' @param seed The numeric seed to use for simulation
#'
#' @export
generate_lr_data <- function(D = 10, N = 1000, off_diag_cov = 0.7, seed = 4711, scaled = FALSE, sigma_noise = 0.01, round.digits = 7){
  set.seed(seed)
  x <- mvtnorm::rmvnorm(n = N, mean = rep(0, D), sigma = solve(diag(D) + off_diag_cov - diag(D)*off_diag_cov)/N)
  if (scaled) x <- scale(x); attr(x,"scaled:scale") <- attr(x,"scaled:center") <- NULL
  x <- round(x, digits = round.digits)
  b <- as.matrix(rep(1.0, D))
  y <- as.vector(x%*%b + rnorm(N, sd = sigma_noise))
  y <- round(y, digits = round.digits)
  list(x = x, D = D, N = N, y = y)
}

#' Function to compute log likelihoods for a linear regression model
#'
#' Computes the log likelihood in either the constrained or unconstrained space.
#'
#' @param b matrix of size (\eqn{S \times P}) with (S) draws from posteror of dimension P ?
#' @param x matrix of size (\eqn{N \times P}) with N observations from posteror of dimension P ?
#' @param y a vector with observed values of size (\eqn{N}) with N observations
#' @param sigma a vector with draws from the posterior for sigma.
#' @param par_unconstrained is the sigma on unconstrained space?
#'
#' @export
log_lik_blinreg <- function(b, x, y, sigma, par_unconstrained = FALSE){
  checkmate::assert_class(b, classes = "matrix")
  checkmate::assert_class(x, classes = "matrix")
  checkmate::assert_numeric(y, len = nrow(x))
  checkmate::assert_numeric(sigma, len = nrow(b))
  checkmate::assert_flag(par_unconstrained)

  bx <- b %*% t(x)
  log_lik <- matrix(0, nrow = nrow(b), ncol = nrow(x))
  if (par_unconstrained) sigma <- exp(sigma)
  for(s in 1:nrow(b)){
    for(i in 1:nrow(x)){
      log_lik[s,i] <-
        dnorm(x = y[i], mean = bx[s,i], sd = sigma[s], log = TRUE)
    }
  }
  colnames(log_lik) <- paste0("log_lik.", 1:nrow(x))
  log_lik
}


# Generate datasets ----
lm_cor_post <- generate_lr_data(D = 10, N = 1000, scaled = TRUE, sigma_noise = 1, off_diag_cov = 0.7, seed = 4711)
lm_indep_post <- generate_lr_data(D = 10, N = 1000, scaled = TRUE, sigma_noise = 1, off_diag_cov = 0, seed = 4712)


# Generate datasets ----
# Train ADVI and Laplace models
stan_lm <- "
data {
int <lower=0> N;
int <lower=0> D;
matrix [N,D] x ;
vector [N] y;
}
parameters {
vector [D] b;
real <lower=0> sigma;
}
model {
target += normal_lpdf(y | x * b, sigma);
target += normal_lpdf(b | 0, 1);
}

generated quantities{
vector[N] log_lik;
// Compute the log likelihoods needed to compute the loo
for (n in 1:N) {
log_lik[n] = normal_lpdf(y[n] | x[n,] * b, sigma);
}
}
"
stm <- stan_model(model_code = stan_lm)

# Laplace test data
opt_independent <- rstan::optimizing(stm, data = lm_indep_post, draws = 5000, constrained = FALSE, hessian = TRUE, seed = 4711)
opt_correlated <- rstan::optimizing(stm, data = lm_indep_post, draws = 5000, constrained = FALSE, hessian = TRUE, seed = 4711)
npar_independent <- ncol(opt_independent$theta_tilde)
npar_correlated <- ncol(opt_correlated$theta_tilde)
b_independent <- as.matrix(opt_independent$theta_tilde[,1:(npar_independent-1)])
b_correlated <- as.matrix(opt_correlated$theta_tilde[,1:(npar_correlated-1)])
sigma_independent <- opt_independent$theta_tilde[,npar_independent]
sigma_correlated <- opt_correlated$theta_tilde[,npar_correlated]
ll_independent <- log_lik_blinreg(b = b_independent, x = lm_indep_post$x, y = lm_indep_post$y, sigma = sigma_independent, par_unconstrained = TRUE)
ll_correlated <- log_lik_blinreg(b = b_correlated, x = lm_cor_post$x, y = lm_cor_post$y, sigma = sigma_correlated, par_unconstrained = TRUE)

# ADVI test data
# The ADVI was run using Yuling Yaos code to get both log_p and log_q.
# When this is included in Stan this code should be possible to run from R directly

advi_fullrank_correlated <- read.table("data-raw/raw_data/output_lr_loo_fullrank_grad100_lr.cor.data1k.R.csv", header = TRUE, sep = ",", fill = TRUE)
advi_fullrank_independent <- read.table("data-raw/raw_data/output_lr_loo_fullrank_grad100_lr.nocor.data1k.R.csv", header = TRUE, sep = ",", fill = TRUE)
advi_fullrank_normal <- read.table("data-raw/raw_data/output_normal_model_fullrank_grad100_normal_data.R.csv", header = TRUE, sep = ",", fill = TRUE)
advi_meanfield_correlated <- read.table("data-raw/raw_data/output_lr_loo_meanfield_grad100_lr.cor.data1k.R.csv", header = TRUE, sep = ",", fill = TRUE)
advi_meanfield_independent <- read.table("data-raw/raw_data/output_lr_loo_meanfield_grad100_lr.nocor.data1k.R.csv", header = TRUE, sep = ",", fill = TRUE)
advi_meanfield_normal <- read.table("data-raw/raw_data/output_normal_model_fullrank_grad100_normal_data.R.csv", header = TRUE, sep = ",", fill = TRUE)
log_lik_idx <- which(grepl(colnames(advi_meanfield_independent), pattern = "log_lik"))
log_lik_idx_normal <- which(grepl(colnames(advi_meanfield_normal), pattern = "log_lik"))

stan_normal <- "
data {
int <lower=0> N;
vector[N] y;
}
parameters {
real mu;
real sigma;
}
model {
target += normal_lpdf(y | mu, sigma);
target += normal_lpdf(mu | 0, 10);
target += gamma_lpdf(sigma | 1, 1);
}

generated quantities{
vector[N] log_lik;
// Compute the log likelihoods needed to compute the loo
for (n in 1:N) {
log_lik[n] = normal_lpdf(y[n] | mu , sigma);
}
}
"
stnm <- stan_model(model_code = stan_normal)

stan_normal_data <- list(N = 3, y = c(-15,0,15))
opt_normal <- rstan::optimizing(stnm, data = stan_normal_data, draws = 1000, constrained = FALSE, hessian = TRUE, seed = 4712)
# s_normal <- rstan::sampling(stnm, data = stan_normal_data)

log_lik <- matrix(0, nrow = 1000, ncol = stan_normal_data$N)
for(i in 1:ncol(log_lik)){
  log_lik[,i] <-
    dnorm(x = stan_normal_data$y[i], mean = opt_normal$theta_tilde[, "mu"], sd = opt_normal$theta_tilde[, "sigma"], log = TRUE)
}

# hist(opt_normal$theta_tilde[, "sigma"])
# hist(unlist(extract(s_normal, "sigma")))
# hist(opt_normal$theta_tilde[, "mu"])
# hist(unlist(extract(s_normal, "mu")))
# hist(log_lik[,1])
# hist(extract(s_normal, "log_lik")$log_lik[,1])


# Create test data list
# We only keep 1000 samples and compute likelihoods for 10 observations to minimize memory burden.
test_data_psis_approximate_posterior <-
             list(laplace_independent = list(log_p = opt_independent$log_p[1:1000],
                                             log_q = opt_independent$log_g[1:1000],
                                             log_liks = ll_independent[1:1000,1:10]),
                  laplace_correlated = list(log_p = opt_correlated$log_p[1:1000],
                                            log_q = opt_correlated$log_g[1:1000],
                                            log_liks = ll_correlated[1:1000,1:10]),
                  laplace_normal = list(log_p = opt_normal$log_p[1:1000],
                                            log_q = opt_normal$log_g[1:1000],
                                            log_liks = log_lik[1:1000,]),
                  meanfield_independent = list(log_p = advi_meanfield_independent$log_p[2:1001],
                                             log_q = advi_meanfield_independent$log_q[2:1001],
                                             log_liks = as.matrix(advi_meanfield_independent[2:1001, log_lik_idx[1:10]])),
                  meanfield_correlated = list(log_p = advi_meanfield_correlated$log_p[2:1001],
                                            log_q = advi_meanfield_correlated$log_q[2:1001],
                                            log_liks = as.matrix(advi_meanfield_correlated[2:1001, log_lik_idx[1:10]])),
                  meanfield_normal = list(log_p = advi_meanfield_normal$log_p[2:1001],
                                          log_q = advi_meanfield_normal$log_q[2:1001],
                                          log_liks = as.matrix(advi_meanfield_normal[2:1001, log_lik_idx_normal])),
                  fullrank_independent = list(log_p = advi_fullrank_independent$log_p[2:1001],
                                               log_q = advi_fullrank_independent$log_q[2:1001],
                                               log_liks = as.matrix(advi_fullrank_independent[2:1001, log_lik_idx[1:10]])),
                  fullrank_correlated = list(log_p = advi_fullrank_correlated$log_p[2:1001],
                                              log_q = advi_fullrank_correlated$log_q[2:1001],
                                              log_liks = as.matrix(advi_fullrank_correlated[2:1001, log_lik_idx[1:10]])),
                  fullrank_normal = list(log_p = advi_fullrank_normal$log_p[2:1001],
                                          log_q = advi_fullrank_normal$log_q[2:1001],
                                          log_liks = as.matrix(advi_fullrank_normal[2:1001, log_lik_idx_normal])))

save(test_data_psis_approximate_posterior, file = "tests/testthat/test_data_psis_approximate_posterior.rda")
