# Define a custom predictive measure

Attaches the name, the orientation, and the standard error of the
difference to a measure function. Pass the result to the `measures`
argument of the `*_pred_measure()` functions.

## Usage

``` r
custom_measure(fun, name, se_diff_fun = NULL, loss = FALSE)
```

## Arguments

- fun:

  A function that returns `estimate`, `se`, and `pointwise`. Only the
  arguments `y`, `ypred`, `mupred`, `ylp`, and `log_weights` that are in
  its signature are supplied.

- name:

  A character string. The measure is reported under this name.

- se_diff_fun:

  How
  [`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
  computes the standard error of a difference: a function `(ref, cmp)`,
  `"sum"`, `"mean"`, or `NULL`. With `NULL` (the default), the standard
  error of the difference is `NA`.

- loss:

  `TRUE` if lower values are better. The default `FALSE` treats the
  measure as a utility.

## Value

`fun` with the attributes `measure_name`, `measure_loss`, and
`measure_se_diff`.

## See also

[`pred_measure()`](https://mc-stan.org/loo/pr/363/reference/pred_measure.md),
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)

## Examples

``` r
my_abs_err <- custom_measure(
  function(y, mupred) {
    pw <- abs(y - colMeans(mupred))
    list(estimate = mean(pw), se = sd(pw) / sqrt(length(pw)), pointwise = pw)
  },
  name = "my_abs_err", se_diff_fun = "mean", loss = TRUE
)
```
