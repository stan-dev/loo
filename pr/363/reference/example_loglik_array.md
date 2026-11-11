# Objects to use in examples and tests

Example pointwise log-likelihood objects to use in demonstrations and
tests. See the **Value** and **Examples** sections below.

## Usage

``` r
example_loglik_array()

example_loglik_matrix()
```

## Value

`example_loglik_array()` returns a 500 (draws) x 2 (chains) x 32
(observations) pointwise log-likelihood array.

`example_loglik_matrix()` returns the same pointwise log-likelihood
values as `example_loglik_array()` but reshaped into a 1000
(draws\*chains) x 32 (observations) matrix.

## Examples

``` r
LLarr <- example_loglik_array()
(dim_arr <- dim(LLarr))
#> [1] 500   2  32
LLmat <- example_loglik_matrix()
(dim_mat <- dim(LLmat))
#> [1] 1000   32

all.equal(dim_mat[1], dim_arr[1] * dim_arr[2])
#> [1] TRUE
all.equal(dim_mat[2], dim_arr[3])
#> [1] TRUE

all.equal(LLarr[, 1, ], LLmat[1:500, ])
#> [1] TRUE
all.equal(LLarr[, 2, ], LLmat[501:1000, ])
#> [1] TRUE
```
