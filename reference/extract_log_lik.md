# Extract pointwise log-likelihood from a Stan model

Convenience function for extracting the pointwise log-likelihood matrix
or array from a `stanfit` object from the rstan package. Note: recent
versions of rstan now include a
[`loo()`](https://mc-stan.org/loo/reference/loo.md) method for `stanfit`
objects that handles this internally.

## Usage

``` r
extract_log_lik(stanfit, parameter_name = "log_lik", merge_chains = TRUE)
```

## Arguments

- stanfit:

  A `stanfit` object (rstan package).

- parameter_name:

  A character string naming the parameter (or generated quantity) in the
  Stan model corresponding to the log-likelihood.

- merge_chains:

  If `TRUE` (the default), all Markov chains are merged together (i.e.,
  stacked) and a matrix is returned. If `FALSE` they are kept separate
  and an array is returned.

## Value

If `merge_chains=TRUE`, an \\S\\ by \\N\\ matrix of (post-warmup)
extracted draws, where \\S\\ is the size of the posterior sample and
\\N\\ is the number of data points. If `merge_chains=FALSE`, an \\I\\ by
\\C\\ by \\N\\ array, where \\I \times C = S\\.

## Details

Stan does not automatically compute and store the log-likelihood. It is
up to the user to incorporate it into the Stan program if it is to be
extracted after fitting the model. In a Stan model, the pointwise log
likelihood can be coded as a vector in the transformed parameters block
(and then summed up in the model block) or it can be coded entirely in
the generated quantities block. We recommend using the generated
quantities block so that the computations are carried out only once per
iteration rather than once per HMC leapfrog step.

For example, the following is the `generated quantities` block for
computing and saving the log-likelihood for a linear regression model
with `N` data points, outcome `y`, predictor matrix `X`, coefficients
`beta`, and standard deviation `sigma`:

`vector[N] log_lik;`

`for (n in 1:N) log_lik[n] = normal_lpdf(y[n] | X[n, ] * beta, sigma);`

## References

Stan Development Team (2017). The Stan C++ Library, Version 2.16.0.
<https://mc-stan.org/>

Stan Development Team (2017). RStan: the R interface to Stan, Version
2.16.1. <https://mc-stan.org/>
