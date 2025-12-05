# Avoiding model refits in leave-one-out cross-validation with moment matching

## Introduction

This vignette demonstrates how to improve the Monte Carlo sampling
accuracy of leave-one-out cross-validation with the **loo** package and
Stan. The **loo** package automatically monitors the sampling accuracy
using Pareto \\k\\ diagnostics for each observation. Here, we present a
method for quickly improving the accuracy when the Pareto diagnostics
indicate problems. This is done by performing some additional
computations using the existing posterior sample. If successful, this
will decrease the Pareto \\k\\ values, making the model assessment more
reliable. **loo** also stores the original Pareto \\k\\ values with the
name `influence_pareto_k` which are not changed. They can be used as a
diagnostic of how much each observation influences the posterior
distribution.

The methodology presented is based on the paper

- Paananen, T., Piironen, J., Buerkner, P.-C., Vehtari, A. (2020).
  Implicitly Adaptive Importance Sampling. [arXiv preprint
  arXiv:1906.08850](https://arxiv.org/abs/1906.08850).

More information about the Pareto \\k\\ diagnostics is given in the
following papers

- Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian
  model evaluation using leave-one-out cross-validation and WAIC.
  *Statistics and Computing*. 27(5), 1413–1432.
  :10.1007/s11222-016-9696-4. Links:
  [published](https://link.springer.com/article/10.1007/s11222-016-9696-4)
  \| [arXiv preprint](https://arxiv.org/abs/1507.04544).

- Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
  Pareto smoothed importance sampling. *Journal of Machine Learning
  Research*, 25(72):1-58. [PDF](https://jmlr.org/papers/v25/19-556.html)

## Example: Eradication of Roaches

We will use the same example as in the vignette [*Using the loo package
(version \>=
2.0.0)*](https://mc-stan.org/loo/articles/loo2-example.html). See the
demo for a description of the problem and data. We will use the same
Poisson regression model as in the case study.

### Coding the Stan model

Here is the Stan code for fitting the Poisson regression model, which we
will use for modeling the number of roaches.

``` r
# Note: some syntax used in this Stan program requires RStan >= 2.26 (or CmdStanR)
# To use an older version of RStan change the line declaring `y` to: int y[N];
stancode <- "
data {
  int<lower=1> K;
  int<lower=1> N;
  matrix[N,K] x;
  array[N] int y;
  vector[N] offset;

  real beta_prior_scale;
  real alpha_prior_scale;
}
parameters {
  vector[K] beta;
  real intercept;
}
model {
  y ~ poisson(exp(x * beta + intercept + offset));
  beta ~ normal(0,beta_prior_scale);
  intercept ~ normal(0,alpha_prior_scale);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = poisson_lpmf(y[n] | exp(x[n] * beta + intercept + offset[n]));
}
"
```

Following the usual approach recommended in [*Writing Stan programs for
use with the loo
package*](http://mc-stan.org/loo/articles/loo2-with-rstan.md), we
compute the log-likelihood for each observation in the
`generated quantities` block of the Stan program.

### Setup

In addition to **loo**, we load the **rstan** package for fitting the
model, and the **rstanarm** package for the data.

``` r
library("rstan")
library("loo")
seed <- 9547
set.seed(seed)
```

### Fitting the model with RStan

Next we fit the model in Stan using the **rstan** package:

``` r
# Prepare data
data(roaches, package = "rstanarm")
roaches$roach1 <- sqrt(roaches$roach1)
y <- roaches$y
x <- roaches[, c("roach1", "treatment", "senior")]
offset <- log(roaches[, "exposure2"])
n <- dim(x)[1]
k <- dim(x)[2]

standata <- list(
  N = n,
  K = k,
  x = as.matrix(x),
  y = y,
  offset = offset,
  beta_prior_scale = 2.5,
  alpha_prior_scale = 5.0
)

# Compile
stanmodel <- stan_model(model_code = stancode)
```

    Running /opt/R/4.5.2/lib/R/bin/R CMD SHLIB foo.c
    using C compiler: ‘gcc (Ubuntu 13.3.0-6ubuntu2~24.04) 13.3.0’
    gcc -std=gnu2x -I"/opt/R/4.5.2/lib/R/include" -DNDEBUG   -I"/home/runner/work/_temp/Library/Rcpp/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/"  -I"/home/runner/work/_temp/Library/RcppEigen/include/unsupported"  -I"/home/runner/work/_temp/Library/BH/include" -I"/home/runner/work/_temp/Library/StanHeaders/include/src/"  -I"/home/runner/work/_temp/Library/StanHeaders/include/"  -I"/home/runner/work/_temp/Library/RcppParallel/include/"  -I"/home/runner/work/_temp/Library/rstan/include" -DEIGEN_NO_DEBUG  -DBOOST_DISABLE_ASSERTS  -DBOOST_PENDING_INTEGER_LOG2_HPP  -DSTAN_THREADS  -DUSE_STANC3 -DSTRICT_R_HEADERS  -DBOOST_PHOENIX_NO_VARIADIC_EXPRESSION  -D_HAS_AUTO_PTR_ETC=0  -include '/home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp'  -D_REENTRANT -DRCPP_PARALLEL_USE_TBB=1   -I/usr/local/include    -fpic  -g -O2  -c foo.c -o foo.o
    In file included from /home/runner/work/_temp/Library/RcppEigen/include/Eigen/Core:19,
                     from /home/runner/work/_temp/Library/RcppEigen/include/Eigen/Dense:1,
                     from /home/runner/work/_temp/Library/StanHeaders/include/stan/math/prim/fun/Eigen.hpp:22,
                     from <command-line>:
    /home/runner/work/_temp/Library/RcppEigen/include/Eigen/src/Core/util/Macros.h:679:10: fatal error: cmath: No such file or directory
      679 | #include <cmath>
          |          ^~~~~~~
    compilation terminated.
    make: *** [/opt/R/4.5.2/lib/R/etc/Makeconf:202: foo.o] Error 1

``` r
# Fit model
fit <- sampling(stanmodel, data = standata, seed = seed, refresh = 0)
print(fit, pars = "beta")
```

    Inference for Stan model: anon_model.
    4 chains, each with iter=2000; warmup=1000; thin=1; 
    post-warmup draws per chain=1000, total post-warmup draws=4000.

             mean se_mean   sd  2.5%   25%   50%   75% 97.5% n_eff Rhat
    beta[1]  0.16       0 0.00  0.16  0.16  0.16  0.16  0.16  2318    1
    beta[2] -0.57       0 0.02 -0.62 -0.59 -0.57 -0.55 -0.52  2467    1
    beta[3] -0.31       0 0.04 -0.38 -0.34 -0.32 -0.29 -0.24  2000    1

    Samples were drawn using NUTS(diag_e) at Fri Dec  5 18:14:01 2025.
    For each parameter, n_eff is a crude measure of effective sample size,
    and Rhat is the potential scale reduction factor on split chains (at 
    convergence, Rhat=1).

Let us now evaluate the predictive performance of the model using
[`loo()`](https://mc-stan.org/loo/dev/reference/loo.md).

``` r
loo1 <- loo(fit)
```

    Replacing NAs in `r_eff` with 1s

    Warning: Some Pareto k diagnostic values are too high. See help('pareto-k-diagnostic') for details.

``` r
loo1
```

    Computed from 4000 by 262 log-likelihood matrix.

             Estimate     SE
    elpd_loo  -5462.4  695.2
    p_loo       259.6   55.7
    looic     10924.8 1390.4
    ------
    MCSE of elpd_loo is NA.
    MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.0]).

    Pareto k diagnostic values:
                             Count Pct.    Min. ESS
    (-Inf, 0.7]   (good)     247   94.3%   104     
       (0.7, 1]   (bad)        8    3.1%   <NA>    
       (1, Inf)   (very bad)   7    2.7%   <NA>    
    See help('pareto-k-diagnostic') for details.

The [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) function
output warnings that there are some observations which are highly
influential, and thus the accuracy of importance sampling is compromised
as indicated by the large Pareto \\k\\ diagnostic values (\> 0.7). As
discussed in the vignette [*Using the loo package (version \>=
2.0.0)*](https://mc-stan.org/loo/articles/loo2-example.html), this may
be an indication of model misspecification. Despite that, it is still
beneficial to be able to evaluate the predictive performance of the
model accurately.

### Moment matching correction for importance sampling

To improve the accuracy of the
[`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) result above, we
could perform leave-one-out cross-validation by explicitly leaving out
single observations and refitting the model using MCMC repeatedly.
However, the Pareto \\k\\ diagnostics indicate that there are 19
observations which are problematic. This would require 19 model refits
which may require a lot of computation time.

Instead of refitting with MCMC, we can perform a faster moment matching
correction to the importance sampling for the problematic observations.
This can be done with the
[`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md)
function in the **loo** package, which takes our existing `loo` object
as input and modifies it. The moment matching requires some evaluations
of the model posterior density. For models fitted with **rstan**, this
can be conveniently done by using the existing `stanfit` object.

First, we show how the moment matching can be used for a model fitted
using **rstan**. It only requires setting the argument `moment_match` to
`TRUE` in the [`loo()`](https://mc-stan.org/loo/dev/reference/loo.md)
function. Optionally, you can also set the argument `k_threshold` which
determines the Pareto \\k\\ threshold, above which moment matching is
used. By default, it operates on all observations whose Pareto \\k\\
value is larger than the sample size (\\S\\) specific threshold
\\\min(1 - 1 / \log\_{10}(S), 0.7)\\ (which is \\0.7\\ for \\S\>2200\\).

``` r
# available in rstan >= 2.21
loo2 <- loo(fit, moment_match = TRUE)
```

    Replacing NAs in `r_eff` with 1s

``` r
loo2
```

    Computed from 4000 by 262 log-likelihood matrix.

             Estimate     SE
    elpd_loo  -5479.1  700.2
    p_loo       276.3   62.6
    looic     10958.2 1400.4
    ------
    MCSE of elpd_loo is NA.
    MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.0]).

    All Pareto k estimates are good (k < 0.7).
    See help('pareto-k-diagnostic') for details.

After the moment matching, all observations have the diagnostic Pareto
\\k\\ less than 0.7, meaning that the estimates are now reliable. The
total `elpd_loo` estimate also changed from `-5457.8` to `-5478.5`,
showing that before moment matching,
[`loo()`](https://mc-stan.org/loo/dev/reference/loo.md) overestimated
the predictive performance of the model.

The updated Pareto \\k\\ values stored in `loo2$diagnostics$pareto_k`
are considered algorithmic diagnostic values that indicate the sampling
accuracy. The original Pareto \\k\\ values are stored in
`loo2$pointwise[,"influence_pareto_k"]` and these are not modified by
the moment matching. These can be considered as diagnostics for how big
influence each observation has on the posterior distribution. In
addition to the Pareto \\k\\ diagnostics, moment matching also updates
the effective sample size estimates.

## Using `loo_moment_match()` directly

The moment matching can also be performed by explicitly calling the
function
[`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md).
This enables its use also for models that are not using **rstan** or
another package with built-in support for
[`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md).
To use
[`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md),
the user must give the model object `x`, the `loo` object, and 5 helper
functions as arguments to
[`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md).
The helper functions are

- `post_draws`
  - A function the takes `x` as the first argument and returns a matrix
    of posterior draws of the model parameters, `pars`.
- `log_lik_i`
  - A function that takes `x` and `i` and returns a matrix (one column
    per chain) or a vector (all chains stacked) of log-likeliood draws
    of the ith observation based on the model `x`. If the draws are
    obtained using MCMC, the matrix with MCMC chains separated is
    preferred.
- `unconstrain_pars`
  - A function that takes arguments `x` and `pars`, and returns
    posterior draws on the unconstrained space based on the posterior
    draws on the constrained space passed via `pars`.
- `log_prob_upars`
  - A function that takes arguments `x` and `upars`, and returns a
    matrix of log-posterior density values of the unconstrained
    posterior draws passed via `upars`.
- `log_lik_i_upars`
  - A function that takes arguments `x`, `upars`, and `i` and returns a
    vector of log-likelihood draws of the `i`th observation based on the
    unconstrained posterior draws passed via `upars`.

Next, we show how the helper functions look like for RStan objects, and
show an example of using
[`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md)
directly. For stanfit objects from **rstan** objects, the functions look
like this:

``` r
# create a named list of draws for use with rstan methods
.rstan_relist <- function(x, skeleton) {
  out <- utils::relist(x, skeleton)
  for (i in seq_along(skeleton)) {
    dim(out[[i]]) <- dim(skeleton[[i]])
  }
  out
}

# rstan helper function to get dims of parameters right
.create_skeleton <- function(pars, dims) {
  out <- lapply(seq_along(pars), function(i) {
    len_dims <- length(dims[[i]])
    if (len_dims < 1) {
      return(0)
    }
    return(array(0, dim = dims[[i]]))
  })
  names(out) <- pars
  out
}

# extract original posterior draws
post_draws_stanfit <- function(x, ...) {
  as.matrix(x)
}

# compute a matrix of log-likelihood values for the ith observation
# matrix contains information about the number of MCMC chains
log_lik_i_stanfit <- function(x, i, parameter_name = "log_lik", ...) {
  loo::extract_log_lik(x, parameter_name, merge_chains = FALSE)[,, i]
}

# transform parameters to the unconstraint space
unconstrain_pars_stanfit <- function(x, pars, ...) {
  skeleton <- .create_skeleton(x@sim$pars_oi, x@par_dims[x@sim$pars_oi])
  upars <- apply(pars, 1, FUN = function(theta) {
    rstan::unconstrain_pars(x, .rstan_relist(theta, skeleton))
  })
  # for one parameter models
  if (is.null(dim(upars))) {
    dim(upars) <- c(1, length(upars))
  }
  t(upars)
}

# compute log_prob for each posterior draws on the unconstrained space
log_prob_upars_stanfit <- function(x, upars, ...) {
  apply(
    upars,
    1,
    rstan::log_prob,
    object = x,
    adjust_transform = TRUE,
    gradient = FALSE
  )
}

# compute log_lik values based on the unconstrained parameters
log_lik_i_upars_stanfit <- function(
  x,
  upars,
  i,
  parameter_name = "log_lik",
  ...
) {
  S <- nrow(upars)
  out <- numeric(S)
  for (s in seq_len(S)) {
    out[s] <- rstan::constrain_pars(x, upars = upars[s, ])[[parameter_name]][i]
  }
  out
}
```

Using these function, we can call
[`loo_moment_match()`](https://mc-stan.org/loo/dev/reference/loo_moment_match.md)
to update the existing `loo` object.

``` r
loo3 <- loo::loo_moment_match.default(
  x = fit,
  loo = loo1,
  post_draws = post_draws_stanfit,
  log_lik_i = log_lik_i_stanfit,
  unconstrain_pars = unconstrain_pars_stanfit,
  log_prob_upars = log_prob_upars_stanfit,
  log_lik_i_upars = log_lik_i_upars_stanfit
)
loo3
```

    Computed from 4000 by 262 log-likelihood matrix.

             Estimate     SE
    elpd_loo  -5479.1  700.2
    p_loo       276.3   62.6
    looic     10958.2 1400.4
    ------
    MCSE of elpd_loo is NA.
    MCSE and ESS estimates assume MCMC draws (r_eff in [0.4, 1.0]).

    All Pareto k estimates are good (k < 0.7).
    See help('pareto-k-diagnostic') for details.

As expected, the result is identical to the previous result of
`loo2 <- loo(fit, moment_match = TRUE)`.

## References

Gelman, A., and Hill, J. (2007). *Data Analysis Using Regression and
Multilevel Hierarchical Models.* Cambridge University Press.

Stan Development Team (2020) *RStan: the R interface to Stan, Version
2.21.1* <https://mc-stan.org>

Paananen, T., Piironen, J., Buerkner, P.-C., Vehtari, A. (2021).
Implicitly adaptive importance sampling. *Statistics and Computing*, 31,
16. :10.1007/s11222-020-09982-2. arXiv preprint arXiv:1906.08850.

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*. 27(5), 1413–1432. :10.1007/s11222-016-9696-4. Links:
[published](https://link.springer.com/article/10.1007/s11222-016-9696-4)
\| [arXiv preprint](https://arxiv.org/abs/1507.04544).

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research*, 25(72):1-58. [PDF](https://jmlr.org/papers/v25/19-556.html)
