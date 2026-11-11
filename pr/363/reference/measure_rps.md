# Ranked Probability Score (RPS, SRPS, CRPS, SCRPS)

Computes proper scoring rules based on the ranked probability score
family, covering both discrete and continuous outcomes, with optional
scaling. The specific score computed depends on the type of `y` and
`ypred` and the value of `scaled`:

## Usage

``` r
measure_rps(y, ypred, log_weights = NULL, pointwise = NULL, scaled = FALSE)
```

## Arguments

- y:

  A numeric vector of \\n\\ observed outcomes. May be integer-valued
  (for RPS/SRPS) or continuous (for CRPS/SCRPS).

- ypred:

  A numeric matrix of posterior predictive draws with dimensions \\S
  \times n\\ (draws × observations).

- log_weights:

  Optional numeric matrix of unnormalized log-importance weights with
  dimensions \\S \times n\\. Weights are column-normalized before
  computing each per-observation contribution.

- pointwise:

  Optional numeric vector of precomputed pointwise values, on the scale
  of the measure that is returned: RPS/CRPS when `scaled = FALSE` and
  SRPS/SCRPS when `scaled = TRUE`. The values are used as they are, with
  no further transformation; in particular, unscaled values are not
  converted to the scaled variant. If provided, `y`, `ypred`, and
  `log_weights` are ignored.

- scaled:

  Logical; if `TRUE`, computes the scaled variant (SRPS for discrete
  outcomes, SCRPS for continuous outcomes). Default is `FALSE`.

## Value

An object of class `"measure"`: a list with:

- `estimates`:

  A 1 \\\times\\ 2 matrix. The row is named after the measure, the
  columns are `Estimate` and `SE` (standard error).

- `pointwise`:

  An \\n \times 1\\ matrix of per-observation values. The column is
  named after the measure.

- `extra`:

  Optional. Auxiliary data for the standard error of a difference. Most
  measures do not set it.

Attributes `measure` (i.e., measure name) and `dims` (draws \\\times\\
observations) are also set. Use
[`print()`](https://rdrr.io/r/base/print.html) for a readable summary.

## Details

|                  |                  |                 |
|------------------|------------------|-----------------|
| `y`/`ypred` type | `scaled = FALSE` | `scaled = TRUE` |
| Discrete         | RPS              | SRPS            |
| Continuous       | CRPS             | SCRPS           |

**Scoring rules:**

- **RPS** (Epstein, 1969): Compares predictive and observed cumulative
  distributions for ordered discrete outcomes.

- **CRPS** (Matheson & Winkler, 1976; Gneiting & Raftery, 2007):
  Generalizes RPS to continuous outcomes. Defined as
  \$\$\mathrm{CRPS}(X; y) = E\[\|X - y\|\] - \frac{1}{2} E\[\|X -
  X'\|\],\$\$ where \\X, X'\\ are independent draws from the predictive
  distribution.

- **SRPS/SCRPS** (Bolin & Wallin, 2023): Scaled variants that are
  invariant to the scale of the predictive distribution. Defined as
  \$\$\mathrm{SCRPS}(X; y) = -\frac{E\[\|X - y\|\]}{E\[\|X - X'\|\]} -
  \frac{1}{2} \log E\[\|X - X'\|\].\$\$

**Estimation:**

Scores are estimated using the probability-weighted moment (PWM)
estimator (Taillardat et al., 2016; Zamo & Naveau, 2018), which requires
only a single set of predictive draws unlike permutation-based
estimators, which require two independent draw sets. The PWM estimator
is unbiased and generally more accurate than single-permutation
estimators. The same estimator is used for both discrete and continuous
outcomes; see Hosking (1990, 1996) for theoretical justification in the
discrete case.

The term \\E\[\|X - X'\|\]\\ is computed as a weighted average over all
pairs of draws, \$\$E\[\|X - X'\|\] = \frac{\sum_i \sum\_{j \neq i} w_i
w_j \|x_i - x_j\|}{1 - \sum_i w_i^2},\$\$ which is evaluated in \\O(S
\log S)\\ using the sorted draws and the cumulative sums of their
weights. If log-weights (`log_weights`) are provided (e.g., PSIS weights
for LOO cross-validation), the \\w_i\\ are the normalized importance
weights and \\E\[\|X - y\|\]\\ is likewise a weighted mean. Otherwise
the weights are equal, \\w_i = 1/S\\, and the expression reduces to the
usual unbiased PWM estimator with the \\1 / (S (S - 1))\\ normalization.
Because the estimate is a convex combination of \\\|x_i - x_j\|\\, it is
always non-negative and invariant to shifts of `ypred`.

**Sign convention:**

Unscaled scores (RPS, CRPS) are returned on the loss scale of Gneiting &
Raftery (2007): lower is better. Scaled scores (SRPS, SCRPS) follow
Bolin & Wallin (2023) and are returned on the utility scale: higher is
better.
[`model_compare()`](https://mc-stan.org/loo/pr/363/reference/model_compare.md)
reports differences for both on the utility scale, so no manual
conversion is needed there. The deprecated
[`crps()`](https://mc-stan.org/loo/pr/363/reference/crps.md) returned
the negated unscaled score; `-measure_rps(...)$pointwise` reproduces it.

## References

Bolin, D. and Wallin, J. (2023). Local scale invariance and robustness
of proper scoring rules. *Statistical Science*, 38(1):140–159.

Epstein, E. S. (1969). A scoring system for probability forecasts of
ranked categories. *Journal of Applied Meteorology*, 8(6):985–987.

Gneiting, T. and Raftery, A. E. (2007). Strictly proper scoring rules,
prediction, and estimation. *Journal of the American Statistical
Association*, 102(477):359–378.

Hosking, J. R. M. (1990). L-moments: Analysis and estimation of
distributions using linear combinations of order statistics. *Journal of
the Royal Statistical Society Series B*, 52(1):105–124.

Hosking, J. R. M. (1996). Some theoretical results concerning L-moments.
Research report RC 14492. IBM Thomas J. Watson Research Division.

Matheson, J. E. and Winkler, R. L. (1976). Scoring rules for continuous
probability distributions. *Management Science*, 22(10):1087–1096.

Taillardat, M., Mestre, O., Zamo, M., and Naveau, P. (2016). Calibrated
ensemble forecasts using quantile regression forests and ensemble model
output statistics. *Monthly Weather Review*, 144(6):2375–2393.

Zamo, M. and Naveau, P. (2018). Estimation of the continuous ranked
probability score with limited information and applications to ensemble
weather forecasts. *Mathematical Geosciences*, 50:209–234.

## Examples

``` r
# Discrete outcomes: RPS
y <- c(2L, 1L, 3L)
ypred <- matrix(c(2, 1, 2, 3, 1, 3), nrow = 2)
measure_rps(y, ypred)
#> 
#> Computed from 2 draws by 3 observations.
#> 
#>     Estimate   SE
#> rps     0.33 0.33

# Discrete outcomes: SRPS (scaled)
measure_rps(y, ypred, scaled = TRUE)
#> 
#> Computed from 2 draws by 3 observations.
#> 
#>      Estimate    SE
#> srps    -0.95  0.29

# Continuous outcomes: CRPS
y_cont <- c(0.5, -1.2, 2.3)
ypred_cont <- matrix(rnorm(200), nrow = 100, ncol = 3)
#> Warning: data length [200] is not a sub-multiple or multiple of the number of columns [3]
measure_rps(y_cont, ypred_cont)
#> 
#> Computed from 100 draws by 3 observations.
#> 
#>     Estimate   SE
#> rps     0.97 0.47

# With importance weights: LOO-CRPS
log_weights <- matrix(rnorm(200), nrow = 100, ncol = 3)
#> Warning: data length [200] is not a sub-multiple or multiple of the number of columns [3]
measure_rps(y_cont, ypred_cont, log_weights = log_weights)
#> 
#> Computed from 100 draws by 3 observations.
#> 
#>     Estimate   SE
#> rps     1.03 0.46
```
