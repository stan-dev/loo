# Pareto smoothed importance sampling (PSIS) using approximate posteriors

Pareto smoothed importance sampling (PSIS) using approximate posteriors

## Usage

``` r
ap_psis(log_ratios, log_p, log_g, ...)

# S3 method for class 'array'
ap_psis(log_ratios, log_p, log_g, ..., cores = getOption("mc.cores", 1))

# S3 method for class 'matrix'
ap_psis(log_ratios, log_p, log_g, ..., cores = getOption("mc.cores", 1))

# Default S3 method
ap_psis(log_ratios, log_p, log_g, ...)
```

## Arguments

- log_ratios:

  The log-likelihood ratios (ie -log_liks)

- log_p:

  The log-posterior (target) evaluated at S samples from the proposal
  distribution (g). A vector of length S.

- log_g:

  The log-density (proposal) evaluated at S samples from the proposal
  distribution (g). A vector of length S.

- ...:

  Currently not in use.

- cores:

  The number of cores to use for parallelization. This defaults to the
  option `mc.cores` which can be set for an entire R session by
  `options(mc.cores = NUMBER)`. The old option `loo.cores` is now
  deprecated but will be given precedence over `mc.cores` until
  `loo.cores` is removed in a future release. **As of version 2.0.0 the
  default is now 1 core if `mc.cores` is not set**, but we recommend
  using as many (or close to as many) cores as possible.

  - Note for Windows 10 users: it is **strongly**
    [recommended](https://github.com/stan-dev/loo/issues/94) to avoid
    using the `.Rprofile` file to set `mc.cores` (using the `cores`
    argument or setting `mc.cores` interactively or in a script is
    fine).

## Methods (by class)

- `ap_psis(array)`: An \\I\\ by \\C\\ by \\N\\ array, where \\I\\ is the
  number of MCMC iterations per chain, \\C\\ is the number of chains,
  and \\N\\ is the number of data points.

- `ap_psis(matrix)`: An \\S\\ by \\N\\ matrix, where \\S\\ is the size
  of the posterior sample (with all chains merged) and \\N\\ is the
  number of data points.

- `ap_psis(default)`: A vector of length \\S\\ (posterior sample size).
