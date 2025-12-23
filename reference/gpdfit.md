# Estimate parameters of the Generalized Pareto distribution

Given a sample \\x\\, Estimate the parameters \\k\\ and \\\sigma\\ of
the generalized Pareto distribution (GPD), assuming the location
parameter is 0. By default the fit uses a prior for \\k\\, which will
stabilize estimates for very small sample sizes (and low effective
sample sizes in the case of MCMC samples). The weakly informative prior
is a Gaussian prior centered at 0.5.

## Usage

``` r
gpdfit(x, wip = TRUE, min_grid_pts = 30, sort_x = TRUE)
```

## Arguments

- x:

  A numeric vector. The sample from which to estimate the parameters.

- wip:

  Logical indicating whether to adjust \\k\\ based on a weakly
  informative Gaussian prior centered on 0.5. Defaults to `TRUE`.

- min_grid_pts:

  The minimum number of grid points used in the fitting algorithm. The
  actual number used is `min_grid_pts + floor(sqrt(length(x)))`.

- sort_x:

  If `TRUE` (the default), the first step in the fitting algorithm is to
  sort the elements of `x`. If `x` is already sorted in ascending order
  then `sort_x` can be set to `FALSE` to skip the initial sorting step.

## Value

A named list with components `k` and `sigma`.

## Details

Here the parameter \\k\\ is the negative of \\k\\ in Zhang & Stephens
(2009).

## References

Zhang, J., and Stephens, M. A. (2009). A new and efficient estimation
method for the generalized Pareto distribution. *Technometrics* **51**,
316-325.

## See also

[`psis()`](https://mc-stan.org/loo/reference/psis.md),
[pareto-k-diagnostic](https://mc-stan.org/loo/reference/pareto-k-diagnostic.md)
