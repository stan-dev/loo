# Model comparison: Explanation of \`model_compare()\`

## Introduction

This vignette introduces the
[`loo::model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
function. Use
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
instead of the deprecated
[`loo::loo_compare()`](https://mc-stan.org/loo/pr/363/reference/loo_compare.md)
(see the note below).

This vignette shows how to use the function. It does not explain how to
interpret the results. It also does not explain what to do in the
Bayesian workflow when you get diagnostic warnings, for example high
Pareto-k values. For these topics, see the “TODO-case
study-model_compare vignette” \[Link-TODO\].

This vignette does not derive the formulas behind `<measure>_diff` and
`<measure>_se_diff` (introduced below). For the derivations, see the
supplement (TODO-Vehtari et al., 2026).

[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
accepts the output of these functions:

- [`loo()`](https://mc-stan.org/loo/pr/363/reference/loo.md),
  [`kfold()`](https://mc-stan.org/loo/pr/363/reference/kfold-generic.md),
  and [`waic()`](https://mc-stan.org/loo/pr/363/reference/waic.md)
- `<type>_pred_measure()`, where `<type>` is `loo`, `kfold`, `test`, or
  `insample`

[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
is new in loo v3.0.0. It compares models on several predictive
performance measures, for example `elpd`, `mae`, `rmse`, and `acc`. For
the full list, see the [Predictive model performance: Overview of
measures](https://mc-stan.org/loo/pr/363/articles/articles-online-only/overview-measures.md)
vignette. You can also define custom measures.

This vignette uses the output of `<type>_pred_measure()`. It also
briefly shows how it compares to the deprecated
[`loo_compare()`](https://mc-stan.org/loo/pr/363/reference/loo_compare.md)
output. The vignette explains:

1.  The arguments and output of
    [`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md).
2.  Which standard error
    [`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
    uses for each built-in measure.
3.  How to compute the difference and its standard error for a
    **custom** measure.

**Relation between `model_compare` and `loo_compare`**

The `loo_compare` function has traditionally been used for model
comparison in loo, accepting any object that inherits from the `loo`
class (such as those produced by
[`loo()`](https://mc-stan.org/loo/pr/363/reference/loo.md),
[`kfold()`](https://mc-stan.org/loo/pr/363/reference/kfold-generic.md),
and [`waic()`](https://mc-stan.org/loo/pr/363/reference/waic.md)).

With the introduction of the *`pred_measure` family*, we have deprecated
`loo_compare` and introduced `model_compare`. While `model_compare`
retains full backward compatibility, it also expands functionality to
handle outputs from the `pred_measure` family.

`loo_compare` still works and is still a generic, so methods registered
by other packages keep dispatching, but it now warns and is frozen at
its previous behavior: `loo` objects compared on ELPD only. The
`pred_measure` sources and custom measures described in this article
require `model_compare`. `loo_compare` will be removed in a future
release. For more details, see
[`vignette("migration-guide", package = "loo")`](https://mc-stan.org/loo/pr/363/articles/migration-guide.md).

## Running example (count data)

To illustrate the functionality of
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
we first introduce a *running example*. In this example, we compare
three Poisson regression models fitted to the **`roaches` data** from
`rstanarm` (Gelman and Hill, 2007). This data set records a
pest-management trial in 262 apartment buildings. The outcome `y` is the
number of roaches caught in traps after treatment and ranges from 0 to
357. The predictors are:

- `roach1`: pre-treatment roach count (used as
  `sqrt_roach1 = sqrt(roach1)`)
- `treatment`: whether the building received pest management
- `senior`: whether the building is restricted to elderly residents

The **three models** used in the following comparison tasks are defined
as follows:

      m1:      y ~ treatment + senior
      m2:      y ~ sqrt_roach1 + senior
      m3:      y ~ sqrt_roach1 + treatment

As **performance criteria** we select six predictive performance
measures: `elpd`, `mae`, `mse`, `rmse`, `r2`, and `rps` computed by
`loo_pred_measure`.

Note that `rps` is the only one of these that needs the posterior
*predictive* draws `ypred`; the error measures are computed from the
posterior expectation `mupred` (see the input requirements table in
[`?loo_pred_measure`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md)).

**Running example**

The running example shows how to call
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md).
It does not show how to interpret the results. For that, see
\[TODO-LINK\]. We therefore do not discuss diagnostic flags or
interpretation of the output here.

``` r

library(loo)
#> This is loo version 3.0.0
#> - Online documentation and vignettes at mc-stan.org/loo
#> - As of v2.0.0 loo defaults to 1 core but we recommend using as many as possible. Use the 'cores' argument or set options(mc.cores = NUM_CORES) for an entire session.
library(brms)
#> Loading required package: Rcpp
#> Loading 'brms' package (version 2.23.2). Useful instructions
#> can be found by typing help('brms'). A more detailed introduction
#> to the package is available through vignette('brms_overview').
#> 
#> Attaching package: 'brms'
#> The following object is masked from 'package:stats':
#> 
#>     ar
```

``` r

measures <- c("elpd", "mae", "mse", "rmse", "r2", "rps")
```

## The classical `loo_compare` output

As the actual model building workflow is not the focus of this tutorial,
we fit the three models here without further comment.
[`add_criterion()`](https://paulbuerkner.com/brms/reference/add_criterion.html)
attaches a `loo` object (as returned by
[`loo::loo()`](https://mc-stan.org/loo/pr/363/reference/loo.md)) to each
fit.

``` r

SEED <- 42

data(roaches, package = "rstanarm")
roaches$sqrt_roach1 <- sqrt(roaches$roach1)

fit_m1 <- brms::brm(
  y ~ treatment + senior,
  data = roaches,
  family = poisson,
  prior = brms::prior(normal(0, 1), class = b),
  chains = 2, iter = 400, seed = SEED, refresh = 0
)
fit_m2 <- update(fit_m1, formula = y ~ sqrt_roach1 + senior,
                 newdata = roaches, refresh = 0)
fit_m3 <- update(fit_m1, formula = y ~ sqrt_roach1 + treatment,
                 newdata = roaches, refresh = 0)

add_loo <- function(fit) {
  brms::add_criterion(fit, criterion = "loo", moment_match = TRUE,
                      save_psis = TRUE)
}
fit_m1 <- add_loo(fit_m1)
fit_m2 <- add_loo(fit_m2)
fit_m3 <- add_loo(fit_m3)

loos <- list(m1 = fit_m1$criteria$loo, m2 = fit_m2$criteria$loo,
             m3 = fit_m3$criteria$loo)
```

Passing these `loo` objects to the classical
[`loo_compare()`](https://mc-stan.org/loo/pr/363/reference/loo_compare.md)
reproduces the familiar ELPD-only comparison, and warns once per session
that
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
is now the preferred function:

``` r

loo_compare(loos)
#> Warning: 'loo_compare' is deprecated. Use 'model_compare' instead.
#>  model elpd_diff se_diff p_worse diag_diff       diag_elpd
#>     m3       0.0     0.0      NA           2 k_psis > 0.62
#>     m2    -211.0   226.1    0.82           2 k_psis > 0.62
#>     m1   -2752.6   637.6    1.00           8 k_psis > 0.62
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
```

When we use `model_compare` instead, we get the same output for `loo`
objects:

``` r

model_compare(loos)
#>  model elpd_diff se_diff p_worse diag_diff       diag_elpd
#>     m3       0.0     0.0      NA           2 k_psis > 0.62
#>     m2    -211.0   226.1    0.82           2 k_psis > 0.62
#>     m1   -2752.6   637.6    1.00           8 k_psis > 0.62
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
```

It also supports the `simplify` argument in the print method:

``` r

print(model_compare(loos), simplify = FALSE)
#>  model elpd_diff se_diff p_worse diag_diff       diag_elpd elpd_loo se_elpd_loo
#>     m3       0.0     0.0      NA           2 k_psis > 0.62  -5465.3       669.8
#>     m2    -211.0   226.1    0.82           2 k_psis > 0.62  -5676.3       772.3
#>     m1   -2752.6   637.6    1.00           8 k_psis > 0.62  -8218.0       876.3
#>  p_loo se_p_loo   looic se_looic
#>  216.4     55.4 10930.7   1339.5
#>  181.1     41.8 11352.6   1544.5
#>  271.1     56.4 16436.0   1752.6
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
```

## Basic usage

We now use
[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md)
to compute the predictive performance for each model (`m1, m2, m3`) and
measure (see list above). Finally, we compare the performance of the
three models using
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md).
By default, only the ranking measure is shown. The ranking measure is
`elpd` when all models share it. Otherwise, it is the first shared
measure in alphabetical order.

``` r

fit_measure <- function(fit, measure) {
  loo_pred_measure(
    loo = fit$criteria$loo,
    y = fit$data$y,
    ypred = brms::posterior_predict(fit),
    mupred = brms::posterior_epred(fit),
    measures = measure
  )
}

set.seed(SEED)
m1 <- fit_measure(fit_m1, measures)
#> Warning: The `loo` object was corrected with "moment_match".
#> ! Only "elpd", "mlpd" and "ic" include this correction.
#> ℹ "mae", "mse", "rmse", "r2", and "rps" use the uncorrected PSIS weights.
m2 <- fit_measure(fit_m2, measures)
#> Warning: The `loo` object was corrected with "moment_match".
#> ! Only "elpd", "mlpd" and "ic" include this correction.
#> ℹ "mae", "mse", "rmse", "r2", and "rps" use the uncorrected PSIS weights.
m3 <- fit_measure(fit_m3, measures)
#> Warning: The `loo` object was corrected with "moment_match".
#> ! Only "elpd", "mlpd" and "ic" include this correction.
#> ℹ "mae", "mse", "rmse", "r2", and "rps" use the uncorrected PSIS weights.

(comp <- model_compare(list(m1 = m1, m2 = m2, m3 = m3)))
#> Each measure compared against its own best model (elpd: m3, ...).
#> PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
#>  model bad_k
#>     m1     8
#>     m3     2
#>     m2     2
#> 
#>  model elpd_diff se_diff p_worse diag_diff
#>     m3       0.0     0.0      NA          
#>     m2    -211.0   226.1    0.82          
#>     m1   -2752.6   637.6    1.00
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
#> 
#> Use print(x, measures = "all") to see all measures.
```

When we want to see a **specific measure**, we can change this by
passing the information to the `print` method:

``` r

print(comp, measures = "r2")
#> Each measure compared against its own best model (elpd: m3, ...).
#> PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
#>  model bad_k
#>     m1     8
#>     m3     2
#>     m2     2
#> 
#> -- r2 (vs m2) --
#>  model r2_diff r2_se_diff
#>     m2   0.000      0.000
#>     m3  -0.005      0.071
#>     m1  -0.219      0.092
```

If you want to see **all columns related to a measure** specified in
`measures` use `simplify = FALSE`:

``` r

print(comp, measures = "r2", simplify = FALSE)
#> Each measure compared against its own best model (elpd: m3, ...).
#> PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
#>  model bad_k
#>     m1     8
#>     m3     2
#>     m2     2
#> 
#> -- r2 (vs m2) --
#>  model r2_diff r2_se_diff    r2 se_r2
#>     m2   0.000      0.000 0.249 0.094
#>     m3  -0.005      0.071 0.243 0.131
#>     m1  -0.219      0.092 0.030 0.024
```

If you want to see **all measures** use `measures = "all"`:

``` r

print(comp, measures = "all")
#> Printing 6 measure comparisons; consider `measures = c(...)`.
#> Each measure compared against its own best model (elpd: m3, ...).
#> PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
#>  model bad_k
#>     m1     8
#>     m3     2
#>     m2     2
#> 
#> -- elpd (vs m3) --
#>  model elpd_diff se_diff p_worse diag_diff
#>     m3       0.0     0.0      NA          
#>     m2    -211.0   226.1    0.82          
#>     m1   -2752.6   637.6    1.00          
#> 
#> -- mae (vs m3, sign flipped) --
#>  model mae_diff mae_se_diff
#>     m3     0.00        0.00
#>     m2    -0.12        0.70
#>     m1    -7.71        1.67
#> 
#> -- mse (vs m2, sign flipped) --
#>  model mse_diff mse_se_diff
#>     m2      0.0         0.0
#>     m3    -13.5       180.6
#>     m1   -563.4       289.6
#> 
#> -- r2 (vs m2) --
#>  model r2_diff r2_se_diff
#>     m2   0.000      0.000
#>     m3  -0.005      0.071
#>     m1  -0.219      0.092
#> 
#> -- rmse (vs m2, sign flipped) --
#>  model rmse_diff rmse_se_diff
#>     m2       0.0          0.0
#>     m3      -0.2          2.0
#>     m1      -6.0          2.9
#> 
#> -- rps (vs m3, sign flipped) --
#>  model rps_diff rps_se_diff
#>     m3     0.00        0.00
#>     m2    -0.12        0.69
#>     m1    -7.34        1.65
#> 
#> All differences: 0 = best model, negative = worse.
#> Signs flipped for loss measures: mae, mse, rmse, rps.
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
```

### Output of `model_compare()`

Using `print(., measures = "all")` returns a subsection for each
measure. The *rows* represent the different models, with the best model
always in the first row, followed by the second best, and so on. The
*columns*:

- `model`: Name of each model
- `<measure>_diff`: The point estimate representing the difference in
  the measure between the current model and the reference model
  (calculated as model - reference). For the reference model itself,
  this value is always zero.
- `<measure>_se_diff`\\^{(1)}\\: The standard error of the difference
- `p_worse`: (only for `elpd`) Estimated probability that a model has
  worse predictive performance than the reference model
- `diag_diff`: (only for `elpd`) Diagnostic column flagging
  miscalibration of the normal approximation behind `p_worse` and
  `se_diff`

Above the measure subsections,
[`print()`](https://rdrr.io/r/base/print.html) names the reference model
each set of differences was computed against, and, for
[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md)
inputs, reports the PSIS-LOO diagnostic once per model:

- `bad_k`: Number of observations whose Pareto \\\hat{k}\\ exceeds the
  PSIS reliability threshold. Pareto \\\hat{k}\\ describes a model’s
  PSIS-LOO approximation.

This block is omitted when no model is flagged, and for `kfold`, `test`,
and `insample` sources, which have no Pareto \\\hat{k}\\ diagnostic.

\\^{(1)}\\ The only measure that deviates from this naming convention is
`elpd`, where the column is called `se_diff` instead of `elpd_se_diff`
for backwards compatibility.

### The `measure_info` attribute

Each model has a `measure_info` attribute, which provides meta
information for each performance measure.

``` r

meta <- attr(m1, "measure_info")

# measure info for each single measure (e.g., R2)
names(meta$r2)
#> [1] "loss"        "diff_method" "se_diff_fun" "extra"
```

The information is stored in different *slots*:

| slot | description | functionality |
|:---|:---|:---|
| `loss` | Whether the measure is a *loss* (`TRUE`; lower values of the measure are better) or a *utility* (`FALSE`, higher values are better) | During model comparison, loss measures have their sign flipped so that all differences are reported on a utility scale (higher is better) |
| `diff_method` | How paired pointwise differences are aggregated within model comparison: `"sum"`, `"mean"`, `"measure_specific"`, or `"custom"` | Selects the formula for `<measure>_se_diff` (see the table in Section [Computation of the standard error](#computation-of-the-standard-error-measure_se_diff)) |
| `se_diff_fun` | Name of the measure’s own standard-error function for computing `<measure>_se_diff`, only set when `diff_method = "measure_specific"` (`"rmse"`, `"r2"`, `"bacc"`) | Points to the built-in function that computes `<measure>_se_diff` for measures whose estimate is not a plain sum or mean of pointwise values |
| `extra` | Optional auxiliary data the measure stored for its `se_diff_fun`, on the measure’s natural scale (e.g., `mse_y_i` for `r2`, `class_id` for `bacc`); `NULL` for measures that need none | Carries the quantities the standard error needs but the pointwise values do not contain; passed to `se_diff_fun()` as part of each model’s input |

## How `model_compare()` computes the output

This section relates `<measure>_diff` and `<measure>_se_diff` to the
model objects. All measures use the same point estimate. The standard
error depends on the measure. For the formulas, see the supplement
(TODO-Vehtari et al., 2026, Section TODO-ARXIV-SECTION).

### Computation of the point estimate `<measure>_diff`

The estimates `elpd_diff`, `mae_diff`, `rmse_diff`, and `r2_diff` in the
output of the running example represent the point estimates of the
pairwise model comparisons. They are computed as the measure-specific
difference between two model estimates (comparison model - reference
model).

For example, the difference in `elpd` between model `m1` (comparison)
and model `m3` (reference) is computed as follows:

``` r

m1_elpd <- m1$estimates["elpd_loo", "Estimate"]
m3_elpd <- m3$estimates["elpd_loo", "Estimate"]

m1_elpd - m3_elpd
```

This computation is done for each measure and each pairwise model
comparison.

By default, however, the reference model is chosen *per measure* as the
best model for that measure. The reference used for each measure is
recorded in the `compare_reference` attribute:

``` r

models <- list(m1 = m1, m2 = m2, m3 = m3)
ref <- attr(comp, "compare_reference")
ref
#> elpd  mae  mse   r2 rmse  rps 
#> "m3" "m3" "m2" "m2" "m2" "m3"
```

The following code block shows, for model `m1`, the hand-computed
differences against each measure’s own reference model and the
differences reported in the `comp` object:

``` r

est <- function(x, measure) x$estimates[measure, "Estimate"]
measures_all <- c("elpd", "mse", "mae", "rmse", "r2")

data.frame(
  reference = ref[measures_all],
  by_hand   = sapply(measures_all, function(m) {
    est(m1, paste0(m, "_loo")) - est(models[[ref[[m]]]], paste0(m, "_loo"))
  }),
  reported  = sapply(measures_all, function(m) {
    comp[[paste0(m, "_diff")]][comp$model == "m1"]
  })
)
#>      reference       by_hand      reported
#> elpd        m3 -2752.6484331 -2752.6484331
#> mse         m2   563.4409737  -563.4409737
#> mae         m3     7.7100076    -7.7100076
#> rmse        m2     5.9957668    -5.9957668
#> r2          m2    -0.2187693    -0.2187693
```

Note that in the above output, the sign for some measures is flipped
between the `by_hand` and `reported` column, specifically for `mse`,
`mae`, and `rmse`. These measures are *losses* (see callout box below).

However, when we compare all measures together, we want them all to be
on the same scale, which in our case is the **utility scale**.
Therefore, losses are converted to the utility scale. As such, higher
values are always better. Consequently, the best model has the highest
value, and when the best model is also the reference model, the
estimated difference `<measure>_diff` (computed as comp - reference) is
always non-positive.

When a measure is converted to the utility scale,
[`print()`](https://rdrr.io/r/base/print.html) marks it with “sign
flipped” in the table header. It also names it below the tables:

``` r

comp <- model_compare(list(m1 = m1, m2 = m2, m3 = m3))
print(comp, measures = c("r2", "mae"))
#> Each measure compared against its own best model (elpd: m3, ...).
#> PSIS-LOO unreliable for all 3 models (k_psis > 0.62); measures may be biased.
#>  model bad_k
#>     m1     8
#>     m3     2
#>     m2     2
#> 
#> -- r2 (vs m2) --
#>  model r2_diff r2_se_diff
#>     m2   0.000      0.000
#>     m3  -0.005      0.071
#>     m1  -0.219      0.092
#> 
#> -- mae (vs m3, sign flipped) --
#>  model mae_diff mae_se_diff
#>     m3     0.00        0.00
#>     m2    -0.12        0.70
#>     m1    -7.71        1.67
#> 
#> All differences: 0 = best model, negative = worse.
#> Signs flipped for loss measures: mae.
```

**Utility vs. loss scale**

We call a measure a **loss** when it has, by default, the orientation
*lower is better*.

By contrast, we call a measure a **utility** if it has, by default, the
orientation *higher is better*.

Which orientation a measure has can be inspected from the `loss` slot of
the `measure_info` attribute of the model object:

``` r

sapply(meta, function(x) x$loss)
#>  elpd   mae   mse  rmse    r2   rps 
#> FALSE  TRUE  TRUE  TRUE FALSE  TRUE
```

### Computation of the standard error `<measure>_se_diff`

The standard error summarizes the uncertainty in the estimated
difference. Depending on the measure, different approaches for computing
the standard error are used. We summarize them in four different
categories:

| `diff_method` | SE of the difference | Measures | Supplement |
|:---|:---|:---|:---|
| `"sum"` | `sqrt(N) * sd(d_i)` | `elpd`, `ic` | Eq. TODO-ARXIV |
| `"mean"` | `sd(d_i) / sqrt(N)` | `mlpd`, `mae`, `mse`, `acc`, `brier`, `rps`, `srps` | Eq. TODO-ARXIV |
| `"measure_specific"` | the measure’s own `se_diff_fun` | `rmse`, `r2`, `bacc` | Eq. TODO-ARXIV |
| `"custom"` | declared via `custom_measure(se_diff_fun = )` | custom measures | – |

where `N` is the number of observations and `d_i` the paired difference
of the pointwise estimates of the two models (`comparison - reference`).

The user can check which approach is used for each measure via the
`diff_method` slot in the `measure_info` attribute of the model object:

``` r

sapply(meta, function(x) x$diff_method)
#>               elpd                mae                mse               rmse 
#>              "sum"             "mean"             "mean" "measure_specific" 
#>                 r2                rps 
#> "measure_specific"             "mean"
```

## Custom measures

A custom measure is any function returning a list with `estimate`, `se`,
and `pointwise`. Define it with
[`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md),
which attaches three declarations to the function:

| argument of [`custom_measure()`](https://mc-stan.org/loo/pr/363/reference/custom_measure.md) | attribute | meaning |
|:---|:---|:---|
| `name` | `measure_name` | the name under which the measure is reported |
| `loss` | `measure_loss` | `TRUE` if lower values are better (default `FALSE`) |
| `se_diff_fun` | `measure_se_diff` | how [`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md) computes the SE of a difference |

A custom measure always gets `diff_method = "custom"`. Setting the
attributes directly with [`attr()`](https://rdrr.io/r/base/attr.html)
gives the same result.

The argument `se_diff_fun` accepts one of three things:

| value | meaning |
|:---|:---|
| `"mean"` | estimate is the mean of the pointwise values. SE is computed as `sd(d_i) / sqrt(N)` |
| `"sum"` | estimate is the sum of the pointwise values. SE is computed as `sqrt(N) * sd(d_i)` |
| a function | custom function for computing the SE, called as `fn(ref, cmp)` |

With `se_diff_fun = NULL` (the default), the measure declares nothing.
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
then reports the difference with an `NA` standard error and a message.
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
needs no further argument in either case.

We elaborate on these points in the following examples.

### Example (one custom measure)

In the following, we provide a custom measure that is not part of the
built-in measures. The Huber loss is quadratic for small residuals and
linear for large residuals. We define it as follows:

``` r

huber_fn <- function(y, mupred) {
  delta <- 10
  r <- y - colMeans(mupred)
  l <- ifelse(abs(r) <= delta, 0.5 * r^2, delta * (abs(r) - 0.5 * delta))
  list(estimate = mean(l), se = sd(l) / sqrt(length(l)), pointwise = l)
}
huber_se_fn <- function(ref, cmp) {
  d <- cmp$pointwise - ref$pointwise
  sd(d) / sqrt(length(d))
}

huber_measure <- custom_measure(
  fun = huber_fn,
  name = "huber",
  se_diff_fun = huber_se_fn,
  loss = TRUE
)
```

With `name`, we report the measure as “huber”. With `loss = TRUE`, we
declare that it is a loss rather than a utility (i.e., lower values are
better). Note that a custom measure is by default considered to be a
utility. The argument `se_diff_fun` declares the standard error of the
difference. It is called as `huber_se_fn(ref, cmp)`, where each argument
is a list describing **one** model with elements `estimate`, `se`,
`pointwise`, and `extra`. This is the general form. The shorthands
`"mean"` and `"sum"` cover the two common cases.

A custom measure is by default considered to be a utility (i.e., higher
values are better).

This can be changed with `custom_measure(loss = TRUE)` (see example
above).

First, we compute the predictive performance measures for two different
models using `"rmse"` and the custom Huber loss as performance criteria.

``` r

h1 <- fit_measure(fit_m1, measure = list("rmse", huber_measure))
#> Warning: The `loo` object was corrected with "moment_match".
#> ! Only "elpd", "mlpd" and "ic" include this correction.
#> ℹ "rmse" and "huber" use the uncorrected PSIS weights.
h3 <- fit_measure(fit_m3, measure = list("rmse", huber_measure))
#> Warning: The `loo` object was corrected with "moment_match".
#> ! Only "elpd", "mlpd" and "ic" include this correction.
#> ℹ "rmse" and "huber" use the uncorrected PSIS weights.
```

(Remember that we have created the `fit_measure` helper at the
[beginning of this tutorial](#example). It simply wraps
[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md).)

Inspecting the attribute `measure_info` shows that the Huber loss uses
`diff_method = "custom"`, and that our declaration is recorded as
`se_diff_fun`.

``` r

str(attr(h1, "measure_info"))
#> List of 3
#>  $ elpd :List of 3
#>   ..$ loss       : logi FALSE
#>   ..$ diff_method: chr "sum"
#>   ..$ se_diff_fun: NULL
#>  $ rmse :List of 3
#>   ..$ loss       : logi TRUE
#>   ..$ diff_method: chr "measure_specific"
#>   ..$ se_diff_fun: chr "rmse"
#>  $ huber:List of 3
#>   ..$ loss       : logi TRUE
#>   ..$ diff_method: chr "custom"
#>   ..$ se_diff_fun:function (ref, cmp)
```

[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
reads the declaration from the objects, so the comparison needs no
further argument:

``` r

comp_h <- model_compare(list(m3 = h3, m1 = h1))
print(comp_h, measures = "huber")
#> Each measure compared against its own best model (huber: m3, rmse: m3).
#> PSIS-LOO unreliable for both models (k_psis > 0.62); measures may be biased.
#>  model bad_k
#>     m1     8
#>     m3     2
#> 
#> -- huber (vs m3, sign flipped) --
#>  model huber_diff huber_se_diff
#>     m3        0.0           0.0
#>     m1      -75.0          16.2
#> 
#> All differences: 0 = best model, negative = worse.
#> Signs flipped for loss measures: huber.
```

As discussed above, the difference estimates are always reported on the
utility scale. Therefore `huber_diff` for `m1` is negative.

### Extended example (multiple custom measures)

In the following, we extend the example above by adding one further
custom measure and show how this is implemented. As an additional custom
measure we use RMSE normalized by the standard deviation of the outcome:

``` r

nrmse_fn <- function(y, mupred) {
  sqe <- (y - colMeans(mupred))^2
  list(
    estimate = sqrt(mean(sqe)) / sd(y),
    se = sqrt(var(sqe) / length(sqe)) / (2 * sqrt(mean(sqe))) / sd(y),
    pointwise = sqe,
    extra = list(sd_y = sd(y))
  )
}

nrmse_se_fn <- function(ref, cmp) {
  n <- length(ref$pointwise)
  mse_ref <- mean(ref$pointwise)
  mse_cmp <- mean(cmp$pointwise)
  se_ref <- sqrt(var(ref$pointwise) / n)
  se_cmp <- sqrt(var(cmp$pointwise) / n)
  cov_rc <- sum((cmp$pointwise - mse_cmp) * (ref$pointwise - mse_ref)) /
    (n * (n - 1))

  rel_ref <- se_ref / sqrt(mse_ref)
  rel_cmp <- se_cmp / sqrt(mse_cmp)
  rho <- cov_rc / (se_ref * se_cmp)

  # the RMSE formula of the supplement (Eq. TODO-ARXIV), written so
  # that a model compared against itself cancels to exactly zero
  se_rmse <- 0.5 * sqrt(
    (rel_cmp - rel_ref)^2 + 2 * rel_cmp * rel_ref * max(1 - rho, 0)
  )

  se_rmse / ref$extra$sd_y
}

nrmse_measure <- custom_measure(
  fun = nrmse_fn,
  name = "nrmse",
  se_diff_fun = nrmse_se_fn,
  loss = TRUE
)
```

Subsequently, we compute the predictive performance measures for two
different models using `"rmse"`, the custom Huber loss, and the custom
normalized RMSE (`"nrmse"`) as performance criteria.

``` r

n1 <- fit_measure(fit_m1, list("rmse", huber_measure, nrmse_measure))
#> Warning: The `loo` object was corrected with "moment_match".
#> ! Only "elpd", "mlpd" and "ic" include this correction.
#> ℹ "rmse", "huber", and "nrmse" use the uncorrected PSIS weights.
n3 <- fit_measure(fit_m3, list("rmse", huber_measure, nrmse_measure))
#> Warning: The `loo` object was corrected with "moment_match".
#> ! Only "elpd", "mlpd" and "ic" include this correction.
#> ℹ "rmse", "huber", and "nrmse" use the uncorrected PSIS weights.
```

Finally, we perform model comparison with multiple custom measures. Both
measures declare their own `se_diff_fun`, so no argument is needed:

``` r

comp_n <- model_compare(list(m3 = n3, m1 = n1))
print(comp_n, measures = "all")
#> Each measure compared against its own best model (huber: m3, nrmse: m3, rmse:
#> m3).
#> PSIS-LOO unreliable for both models (k_psis > 0.62); measures may be biased.
#>  model bad_k
#>     m1     8
#>     m3     2
#> 
#> -- huber (vs m3, sign flipped) --
#>  model huber_diff huber_se_diff
#>     m3        0.0           0.0
#>     m1      -75.0          16.2
#> 
#> -- nrmse (vs m3, sign flipped) --
#>  model nrmse_diff nrmse_se_diff
#>     m3      0.000         0.000
#>     m1     -0.124         0.075
#> 
#> -- rmse (vs m3, sign flipped) --
#>  model rmse_diff rmse_se_diff
#>     m3       0.0          0.0
#>     m1      -5.8          4.0
#> 
#> All differences: 0 = best model, negative = worse.
#> Signs flipped for loss measures: huber, nrmse, rmse.
```

## Comparing other predictive sources

[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
treats every `*_pred_measure()` result the same way. Only the
constructor changes. All models in one call must share one source. A
mixed call is an error.

### K-fold cross-validation

[`brms::kfold()`](https://mc-stan.org/loo/reference/kfold-generic.html)
refits each model `K` times, so this is the most expensive step in this
article.

``` r

set.seed(SEED)
folds <- kfold_split_random(K = 5, N = nrow(roaches))

kfold_measure <- function(fit, measure) {
  kf <- brms::kfold(fit, folds = folds, save_fits = TRUE)

  kfold_pred_measure(
    y = fit$data$y,
    ypred = brms::kfold_predict(kf, method = "predict")$yrep,
    mupred = brms::kfold_predict(kf, method = "fitted")$yrep,
    kfold = kf,
    measures = measure
  )
}

k1 <- kfold_measure(fit_m1, measures)
#> Warning: Found 21 observations with a pareto_k > 0.61568910657988 in model
#> 'fit'.
k2 <- kfold_measure(fit_m2, measures)
#> Warning: Found 13 observations with a pareto_k > 0.61568910657988 in model
#> 'fit'.
k3 <- kfold_measure(fit_m3, measures)
#> Warning: Found 21 observations with a pareto_k > 0.61568910657988 in model
#> 'fit'.

model_compare(list(m1 = k1, m2 = k2, m3 = k3))
#> Each measure compared against its own best model (elpd: m3, ...).
#> Predictive measures evaluated on 5-fold cross-validation.
#>  model elpd_diff se_diff p_worse diag_diff
#>     m3       0.0     0.0      NA          
#>     m2    -273.7   222.1    0.89          
#>     m1   -2847.8   610.6    1.00
```

[`print()`](https://rdrr.io/r/base/print.html) names the source above
the table. The PSIS-LOO block is gone, because Pareto \\\hat{k}\\
describes the LOO approximation only.

### Held-out test data

Here we hold out 60 buildings, refit each model on the rest, and score
the held-out ones. `ylp_test` supplies `elpd_test`; `ylp` stays the
training log-likelihood.

``` r

set.seed(SEED)
test_id <- sample(nrow(roaches), 60)
roaches_train <- roaches[-test_id, ]
roaches_test <- roaches[test_id, ]

test_measure <- function(fit, measure) {
  fit_train <- update(fit, newdata = roaches_train, refresh = 0)
  test_pred_measure(
    y = roaches_test$y,
    ypred = brms::posterior_predict(fit_train, newdata = roaches_test),
    mupred = brms::posterior_epred(fit_train, newdata = roaches_test),
    ylp = brms::log_lik(fit_train),
    ylp_test = brms::log_lik(fit_train, newdata = roaches_test),
    measures = measure
  )
}

t1 <- test_measure(fit_m1, measures)
#> Start sampling
t2 <- test_measure(fit_m2, measures)
#> Start sampling
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
t3 <- test_measure(fit_m3, measures)
#> Start sampling
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

model_compare(list(m1 = t1, m2 = t2, m3 = t3))
#> Each measure compared against its own best model (elpd: m2, ...).
#> Predictive measures evaluated on held-out test data (N = 60).
#>  model elpd_diff se_diff p_worse diag_diff
#>     m2       0.0     0.0      NA          
#>     m3     -46.1    55.1    0.80   N < 100
#>     m1    -666.2   212.0    1.00   N < 100
#> 
#> Diagnostic flags present.
#> See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
#> or https://mc-stan.org/loo/reference/loo-glossary.html.
#> 
#> Use print(x, measures = "all") to see all measures.
```

### In-sample

In-sample scores need no refit. They are optimistically biased and favor
the more complex model, so
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
warns:

``` r

insample_measure <- function(fit, measure) {
  insample_pred_measure(
    y = fit$data$y,
    ypred = brms::posterior_predict(fit),
    mupred = brms::posterior_epred(fit),
    ylp = brms::log_lik(fit),
    measures = measure
  )
}

i1 <- insample_measure(fit_m1, measures)
i2 <- insample_measure(fit_m2, measures)
i3 <- insample_measure(fit_m3, measures)

model_compare(list(m1 = i1, m2 = i2, m3 = i3))
#> Warning: Comparing in-sample predictive measures. In-sample scores are
#> optimistically biased and favor more complex models. For out-of-sample
#> comparison use loo_pred_measure(), kfold_pred_measure(), or
#> test_pred_measure().
#> Each measure compared against its own best model (elpd: m3, ...).
#> Predictive measures evaluated on in-sample (training) data.
#>  model elpd_diff se_diff p_worse diag_diff
#>     m3       0.0     0.0      NA          
#>     m2    -246.3   210.0    0.88          
#>     m1   -2698.0   610.3    1.00
#> Use print(x, measures = "all") to see all measures.
```

Everything else in this article, custom measures, the `print` arguments
`measures` and `simplify`, works unchanged for all four sources.

## Summary

| Measure | diff_method | Standard error |
|:---|:---|:---|
| elpd, ic | sum | sqrt(N) \* sd(d_i) |
| mlpd, mae, mse, acc, brier, rps, srps | mean | sd(d_i) / sqrt(N) |
| rmse | measure_specific | delta method on the MSE scale |
| r2 | measure_specific | delta method, including baseline uncertainty |
| bacc | measure_specific | stratified by class, paired within class |
| custom | custom | from se_diff_fun: your function, “sum”, “mean”, or NA |

The things worth remembering:

1.  The difference is always `estimate_cmp - estimate_ref`, reported so
    that negative means worse.
2.  The standard error is always paired, because the models are
    evaluated on the same observations.
3.  `NA` in a standard error column means “not available”, not “zero”
    and not “an error occurred”.
4.  A custom measure is assumed to be a utility. If yours is a loss, say
    so with `custom_measure(loss = TRUE)`, or its differences and its
    ranking will come out reversed.

## References

Sivula, T., Magnusson, M., Matamoros, A. A., and Vehtari, A. (2025).
Uncertainty in Bayesian leave-one-out cross-validation based model
comparison. *Bayesian Analysis*. <doi:10.1214/25-BA1569>

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*, 27(5):1413–1432. <doi:10.1007/s11222-016-9696-4>

Vehtari, A., Bockting, F., Axen, S., and Weber, F. (2026).
*TODO-ARXIV-TITLE*. arXiv preprint arXiv:TODO-ARXIV-ID.

## See also

- [Predictive model performance: Predictive schemes with
  `pred_measure()`](https://mc-stan.org/loo/pr/363/articles/articles-online-only/pred-measure-workflow.md)
- [Predictive model performance: Overview of
  measures](https://mc-stan.org/loo/pr/363/articles/articles-online-only/overview-measures.md)
- [`?"loo-glossary"`](https://mc-stan.org/loo/pr/363/reference/loo-glossary.md)
- [Predictive measures: A formula
  reference](https://mc-stan.org/loo/pr/363/articles/articles-online-only/)
