# Holdout validation and K-fold cross-validation of Stan programs with the loo package

## Introduction

This vignette demonstrates how to do holdout validation and K-fold
cross-validation with **loo** for a Stan program.

## Example: Eradication of Roaches using holdout validation approach

This vignette uses the same example as in the vignettes [*Using the loo
package (version \>=
2.0.0)*](http://mc-stan.org/loo/articles/loo2-example.md) and [*Avoiding
model refits in leave-one-out cross-validation with moment
matching*](https://mc-stan.org/loo/articles/loo2-moment-matching.html).

### Coding the Stan model

Here is the Stan code for fitting a Poisson regression model:

``` r
# Note: some syntax used in this Stan program requires RStan >= 2.26 (or CmdStanR)
# To use an older version of RStan change the line declaring `y` to: int y[N];
stancode <- "
data {
  int<lower=1> K;
  int<lower=1> N;
  matrix[N,K] x;
  array[N] int y;
  vector[N] offset_; // offset is reserved keyword in Stan so use offset_

  real beta_prior_scale;
  real alpha_prior_scale;
}
parameters {
  vector[K] beta;
  real intercept;
}
model {
  y ~ poisson(exp(x * beta + intercept + offset_));
  beta ~ normal(0,beta_prior_scale);
  intercept ~ normal(0,alpha_prior_scale);
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = poisson_lpmf(y[n] | exp(x[n] * beta + intercept + offset_[n]));
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
model. We will also need the **rstanarm** package for the data.

``` r
library("rstan")
library("loo")
seed <- 9547
set.seed(seed)
```

## Holdout validation

For this approach, the model is first fit to the “train” data and then
is evaluated on the held-out “test” data.

### Splitting the data between train and test

The data is divided between train (80% of the data) and test (20%):

``` r
# Prepare data
data(roaches, package = "rstanarm")
roaches$roach1 <- sqrt(roaches$roach1)
roaches$offset <- log(roaches[,"exposure2"])
# 20% of the data goes to the test set:
roaches$test <- 0
roaches$test[sample(.2 * seq_len(nrow(roaches)))] <- 1
# data to "train" the model
data_train <- list(y = roaches$y[roaches$test == 0],
                   x = as.matrix(roaches[roaches$test == 0,
                                         c("roach1", "treatment", "senior")]),
                   N = nrow(roaches[roaches$test == 0,]),
                   K = 3,
                   offset_ = roaches$offset[roaches$test == 0],
                   beta_prior_scale = 2.5,
                   alpha_prior_scale = 5.0
                   )
# data to "test" the model
data_test <- list(y = roaches$y[roaches$test == 1],
                   x = as.matrix(roaches[roaches$test == 1,
                                         c("roach1", "treatment", "senior")]),
                   N = nrow(roaches[roaches$test == 1,]),
                   K = 3,
                   offset_ = roaches$offset[roaches$test == 1],
                   beta_prior_scale = 2.5,
                   alpha_prior_scale = 5.0
                   )
```

### Fitting the model with RStan

Next we fit the model to the “test” data in Stan using the **rstan**
package:

``` r
# Compile
stanmodel <- stan_model(model_code = stancode)
```

    Trying to compile a simple C file

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
fit <- sampling(stanmodel, data = data_train, seed = seed, refresh = 0)
```

We recompute the generated quantities using the posterior draws
conditional on the training data, but we now pass in the held-out data
to get the log predictive densities for the test data. Because we are
using independent data, the log predictive density coincides with the
log likelihood of the test data.

``` r
gen_test <- gqs(stanmodel, draws = as.matrix(fit), data= data_test)
log_pd <- extract_log_lik(gen_test)
```

### Computing holdout elpd:

Now we evaluate the predictive performance of the model on the test data
using [`elpd()`](https://mc-stan.org/loo/dev/reference/elpd.md).

``` r
(elpd_holdout <- elpd(log_pd))
```

    Computed from 4000 by 52 log-likelihood matrix using the generic elpd function

         Estimate    SE
    elpd  -1741.0 290.7
    ic     3482.1 581.5

When one wants to compare different models, the function
[`loo_compare()`](https://mc-stan.org/loo/dev/reference/loo_compare.md)
can be used to assess the difference in performance.

## K-fold cross validation

For this approach the data is divided into folds, and each time one fold
is tested while the rest of the data is used to fit the model (see
Vehtari et al., 2017).

### Splitting the data in folds

We use the data that is already pre-processed and we divide it in 10
random folds using `kfold_split_random`

``` r
# Prepare data
roaches$fold <- kfold_split_random(K = 10, N = nrow(roaches))
```

### Fitting and extracting the log pointwise predictive densities for each fold

We now loop over the 10 folds. In each fold we do the following. First,
we fit the model to all the observations except the ones belonging to
the left-out fold. Second, we compute the log pointwise predictive
densities for the left-out fold. Last, we store the predictive density
for the observations of the left-out fold in a matrix. The output of
this loop is a matrix of the log pointwise predictive densities of all
the observations.

``` r
# Prepare a matrix with the number of post-warmup iterations by number of observations:
log_pd_kfold <- matrix(nrow = 4000, ncol = nrow(roaches))
# Loop over the folds
for(k in 1:10){
  data_train <- list(y = roaches$y[roaches$fold != k],
                   x = as.matrix(roaches[roaches$fold != k,
                                         c("roach1", "treatment", "senior")]),
                   N = nrow(roaches[roaches$fold != k,]),
                   K = 3,
                   offset_ = roaches$offset[roaches$fold != k],
                   beta_prior_scale = 2.5,
                   alpha_prior_scale = 5.0
                   )
  data_test <- list(y = roaches$y[roaches$fold == k],
                   x = as.matrix(roaches[roaches$fold == k,
                                         c("roach1", "treatment", "senior")]),
                   N = nrow(roaches[roaches$fold == k,]),
                   K = 3,
                   offset_ = roaches$offset[roaches$fold == k],
                   beta_prior_scale = 2.5,
                   alpha_prior_scale = 5.0
                   )
  fit <- sampling(stanmodel, data = data_train, seed = seed, refresh = 0)
  gen_test <- gqs(stanmodel, draws = as.matrix(fit), data= data_test)
  log_pd_kfold[, roaches$fold == k] <- extract_log_lik(gen_test)
}
```

### Computing K-fold elpd:

Now we evaluate the predictive performance of the model on the 10 folds
using [`elpd()`](https://mc-stan.org/loo/dev/reference/elpd.md).

``` r
(elpd_kfold <- elpd(log_pd_kfold))
```

    Computed from 4000 by 262 log-likelihood matrix using the generic elpd function

         Estimate     SE
    elpd  -5560.1  730.0
    ic    11120.2 1460.1

If one wants to compare several models (with `loo_compare`), one should
use the same folds for all the different models.

## References

Gelman, A., and Hill, J. (2007). *Data Analysis Using Regression and
Multilevel Hierarchical Models.* Cambridge University Press.

Stan Development Team (2020) *RStan: the R interface to Stan, Version
2.21.1* <https://mc-stan.org>

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*. 27(5), 1413–1432. :10.1007/s11222-016-9696-4. Links:
[published](https://link.springer.com/article/10.1007/s11222-016-9696-4)
\| [arXiv preprint](https://arxiv.org/abs/1507.04544).
