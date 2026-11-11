# Continuously ranked probability score (deprecated)

As of loo 3.0.0, `crps()`, `scrps()`, `loo_crps()`, and `loo_scrps()`
are **deprecated**. Please use
[`measure_rps()`](https://mc-stan.org/loo/pr/363/reference/measure_rps.md)
and
[`measure_srps()`](https://mc-stan.org/loo/pr/363/reference/measure_srps.md)
instead, or
[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md)
for leave-one-out predictive performance. See
[`vignette("migration-guide", package = "loo")`](https://mc-stan.org/loo/pr/363/articles/migration-guide.md)
for a full mapping table.

## Usage

``` r
crps(x, ...)

scrps(x, ...)

loo_crps(x, ...)

loo_scrps(x, ...)

# S3 method for class 'matrix'
crps(x, x2, y, ..., permutations = 1)

# S3 method for class 'numeric'
crps(x, x2, y, ..., permutations = 1)

# S3 method for class 'matrix'
loo_crps(
  x,
  x2,
  y,
  log_lik,
  ...,
  permutations = 1,
  r_eff = 1,
  cores = getOption("mc.cores", 1)
)

# S3 method for class 'matrix'
scrps(x, x2, y, ..., permutations = 1)

# S3 method for class 'numeric'
scrps(x, x2, y, ..., permutations = 1)

# S3 method for class 'matrix'
loo_scrps(
  x,
  x2,
  y,
  log_lik,
  ...,
  permutations = 1,
  r_eff = 1,
  cores = getOption("mc.cores", 1)
)
```

## Arguments

- x:

  A `S` by `N` matrix (draws by observations), or a vector of length `S`
  when only single observation is provided in `y`.

- ...:

  Passed on to
  [`E_loo()`](https://mc-stan.org/loo/pr/363/reference/E_loo.md) in the
  `loo_*()` version of these functions.

- x2:

  Independent draws from the same distribution as draws in `x`. Should
  be of the identical dimension.

- y:

  A vector of observations or a single value.

- permutations:

  An integer, with default value of 1, specifying how many times the
  expected value of \|X - X'\| (`|x - x2|`) is computed. The row order
  of `x2` is shuffled as elements `x` and `x2` are typically drawn given
  the same values of parameters. This happens, e.g., when one calls
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  twice for a fitted rstanarm or brms model. Generating more
  permutations is expected to decrease the variance of the computed
  expected value.

- log_lik:

  A log-likelihood matrix the same size as `x`.

- r_eff:

  An optional vector of relative effective sample size estimates
  containing one element per observation. See
  [`psis()`](https://mc-stan.org/loo/pr/363/reference/psis.md) for
  details.

- cores:

  The number of cores to use for parallelization of `[psis()]`. See
  [`psis()`](https://mc-stan.org/loo/pr/363/reference/psis.md) for
  details.

## Value

A list containing two elements: `estimates` and `pointwise`. The former
reports estimator and standard error and latter the pointwise values.
Following Bolin & Wallin (2023), a larger value is better.

## Details

The new functions use a single matrix of posterior predictive draws
(`ypred`) and a PWM estimator. The deprecated functions require two
independent draw matrices (`x` and `x2`) and use a permutation-based
estimator. See
[`measure_rps()`](https://mc-stan.org/loo/pr/363/reference/measure_rps.md)
for migration details.

The `crps()` and `scrps()` functions and their `loo_*()` counterparts
can be used to compute the continuously ranked probability score (CRPS)
and scaled CRPS (SCRPS) (as defined by Bolin and Wallin, 2023). CRPS is
a proper scoring rule, and strictly proper when the first moment of the
predictive distribution is finite. Both can be expressed in terms of
samples form the predictive distribution. See, for example, a paper by
Gneiting and Raftery (2007) for a comprehensive discussion on CRPS.

To compute (S)CRPS, the user needs to provide two sets of draws, `x` and
`x2`, from the predictive distribution. This is due to the fact that
formulas used to compute CRPS involve an expectation of the absolute
difference of `x` and `x2`, both having the same distribution. See the
`permutations` argument, as well as Gneiting and Raftery (2007) for
details.

## References

Bolin, D., & Wallin, J. (2023). Local scale invariance and robustness of
proper scoring rules. Statistical Science, 38(1):140-159.

Gneiting, T., & Raftery, A. E. (2007). Strictly Proper Scoring Rules,
Prediction, and Estimation. Journal of the American Statistical
Association, 102(477), 359–378.

## See also

[`measure_rps()`](https://mc-stan.org/loo/pr/363/reference/measure_rps.md),
[`measure_srps()`](https://mc-stan.org/loo/pr/363/reference/measure_srps.md),
[`loo_pred_measure()`](https://mc-stan.org/loo/pr/363/reference/loo_pred_measure.md)

## Examples

``` r
# \dontrun{
# Deprecated:
library(rstanarm)
data("kidiq")
fit <- stan_glm(kid_score ~ mom_hs + mom_iq, data = kidiq)
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 2.2e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 1: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 1: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 1: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.024 seconds (Warm-up)
#> Chain 1:                0.042 seconds (Sampling)
#> Chain 1:                0.066 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 8e-06 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.08 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 2: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 2: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 2: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.025 seconds (Warm-up)
#> Chain 2:                0.044 seconds (Sampling)
#> Chain 2:                0.069 seconds (Total)
#> Chain 2: 
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 3).
#> Chain 3: 
#> Chain 3: Gradient evaluation took 8e-06 seconds
#> Chain 3: 1000 transitions using 10 leapfrog steps per transition would take 0.08 seconds.
#> Chain 3: Adjust your expectations accordingly!
#> Chain 3: 
#> Chain 3: 
#> Chain 3: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 3: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 3: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 3: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 3: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 3: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 3: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 3: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 3: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 3: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 3: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 3: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 3: 
#> Chain 3:  Elapsed Time: 0.028 seconds (Warm-up)
#> Chain 3:                0.045 seconds (Sampling)
#> Chain 3:                0.073 seconds (Total)
#> Chain 3: 
#> 
#> SAMPLING FOR MODEL 'continuous' NOW (CHAIN 4).
#> Chain 4: 
#> Chain 4: Gradient evaluation took 8e-06 seconds
#> Chain 4: 1000 transitions using 10 leapfrog steps per transition would take 0.08 seconds.
#> Chain 4: Adjust your expectations accordingly!
#> Chain 4: 
#> Chain 4: 
#> Chain 4: Iteration:    1 / 2000 [  0%]  (Warmup)
#> Chain 4: Iteration:  200 / 2000 [ 10%]  (Warmup)
#> Chain 4: Iteration:  400 / 2000 [ 20%]  (Warmup)
#> Chain 4: Iteration:  600 / 2000 [ 30%]  (Warmup)
#> Chain 4: Iteration:  800 / 2000 [ 40%]  (Warmup)
#> Chain 4: Iteration: 1000 / 2000 [ 50%]  (Warmup)
#> Chain 4: Iteration: 1001 / 2000 [ 50%]  (Sampling)
#> Chain 4: Iteration: 1200 / 2000 [ 60%]  (Sampling)
#> Chain 4: Iteration: 1400 / 2000 [ 70%]  (Sampling)
#> Chain 4: Iteration: 1600 / 2000 [ 80%]  (Sampling)
#> Chain 4: Iteration: 1800 / 2000 [ 90%]  (Sampling)
#> Chain 4: Iteration: 2000 / 2000 [100%]  (Sampling)
#> Chain 4: 
#> Chain 4:  Elapsed Time: 0.028 seconds (Warm-up)
#> Chain 4:                0.043 seconds (Sampling)
#> Chain 4:                0.071 seconds (Total)
#> Chain 4: 
ypred1 <- posterior_predict(fit)
ypred2 <- posterior_predict(fit)
crps(ypred1, ypred2, y = fit$y)
#> Warning: 'crps.matrix' is deprecated.
#> Use 'measure_rps' instead.
#> See help("Deprecated")
#> $estimates
#>    Estimate          SE 
#> -10.2201033   0.3462796 
#> 
#> $pointwise
#>          1          2          3          4          5          6          7 
#> -24.900751  -9.429146  -7.566819  -4.757940 -21.776399  -7.188831 -30.804403 
#>          8          9         10         11         12         13         14 
#>  -4.757819 -15.637096  -6.468492  -6.359427 -17.768061  -4.239904 -10.673001 
#>         15         16         17         18         19         20         21 
#> -11.561414  -6.692402  -5.390397  -4.318643  -4.445514 -12.452083  -5.665604 
#>         22         23         24         25         26         27         28 
#>  -7.378751  -7.754736  -4.345221 -11.053228  -4.832656  -8.728251  -5.861742 
#>         29         30         31         32         33         34         35 
#>  -4.313555  -4.057937 -15.823206 -30.627095  -6.826580 -12.101544  -5.089308 
#>         36         37         38         39         40         41         42 
#>  -9.514070 -13.078258  -4.691303  -5.944815 -18.215365  -9.966671  -5.378535 
#>         43         44         45         46         47         48         49 
#>  -8.406113  -4.406818 -20.135966  -5.749869 -24.510404  -4.460993 -17.431957 
#>         50         51         52         53         54         55         56 
#>  -4.765400  -6.663401  -5.582078 -19.759953 -14.754404  -5.325623  -8.904427 
#>         57         58         59         60         61         62         63 
#>  -6.197788  -4.937223  -5.284561  -5.616381  -9.268595  -8.761120  -4.327131 
#>         64         65         66         67         68         69         70 
#>  -5.930424 -23.483479 -25.138686 -13.046133  -4.409787  -5.123954  -5.424425 
#>         71         72         73         74         75         76         77 
#>  -5.830134  -5.192773  -4.421525 -14.751313 -11.591811 -16.786968  -4.808163 
#>         78         79         80         81         82         83         84 
#>  -4.873056 -20.363815  -4.659020 -11.978782 -16.056477  -8.690020 -10.145217 
#>         85         86         87         88         89         90         91 
#>  -6.371287  -4.215216 -34.773780  -5.828478  -4.379739  -7.326905  -6.964549 
#>         92         93         94         95         96         97         98 
#> -10.216889  -8.668180  -5.046801 -11.936571  -4.861803  -4.321718  -6.345546 
#>         99        100        101        102        103        104        105 
#> -13.690298 -16.683949  -4.866824 -15.809229 -10.871605  -4.586989  -8.808631 
#>        106        107        108        109        110        111        112 
#>  -4.481927 -12.945530  -7.816259  -5.494882  -8.391380 -29.463484 -13.190499 
#>        113        114        115        116        117        118        119 
#> -26.475558  -7.355715  -5.609233  -5.046325 -25.125862 -29.227117  -5.179128 
#>        120        121        122        123        124        125        126 
#> -15.897256 -11.517018  -5.355114  -4.743005  -4.542042 -19.157431 -13.700954 
#>        127        128        129        130        131        132        133 
#> -11.060062  -5.684123  -4.306502  -8.767889 -23.833348  -8.716688  -6.633172 
#>        134        135        136        137        138        139        140 
#>  -5.609098  -7.071051 -27.677424 -12.652698  -8.112550  -4.550719  -6.658283 
#>        141        142        143        144        145        146        147 
#>  -4.024724  -4.785121  -4.542477  -5.888854  -4.471665  -4.260616  -6.749968 
#>        148        149        150        151        152        153        154 
#> -23.320616 -15.247624 -12.348571  -5.178561 -28.495674  -4.154838  -9.125980 
#>        155        156        157        158        159        160        161 
#>  -8.651110  -4.695778  -7.205551  -7.762747 -10.042472  -4.724869  -5.419278 
#>        162        163        164        165        166        167        168 
#>  -7.321711  -5.710599  -5.684938  -5.157638  -8.308155  -8.532738  -4.083349 
#>        169        170        171        172        173        174        175 
#> -15.950781  -6.366487  -4.809156  -5.605162  -5.755546  -6.754661 -12.383820 
#>        176        177        178        179        180        181        182 
#>  -4.333933  -4.414126  -6.459915  -4.825279  -4.501981  -5.024166  -5.309689 
#>        183        184        185        186        187        188        189 
#>  -4.292834 -13.955145 -19.676055 -14.354093  -5.286629 -20.570589  -4.099824 
#>        190        191        192        193        194        195        196 
#> -10.327451  -8.586051  -4.618217  -7.028630  -6.801230  -7.308523  -5.014850 
#>        197        198        199        200        201        202        203 
#>  -6.764636 -11.187975  -6.635365  -4.580954  -4.725686  -5.049591  -7.791539 
#>        204        205        206        207        208        209        210 
#>  -4.811319  -8.697904  -7.435702  -7.453368  -4.270350 -15.245878  -4.791924 
#>        211        212        213        214        215        216        217 
#>  -7.842196 -12.926809 -39.772207 -19.405911 -11.434198  -6.739385  -6.149844 
#>        218        219        220        221        222        223        224 
#>  -7.846718  -5.139589  -6.247213  -7.814136 -22.926283  -5.032341  -5.200068 
#>        225        226        227        228        229        230        231 
#>  -5.528244  -5.002408  -4.802062  -6.562702 -15.494864  -9.674893  -4.705059 
#>        232        233        234        235        236        237        238 
#>  -4.274015  -8.972172  -7.139933  -6.944002  -4.425469  -7.772937  -4.217749 
#>        239        240        241        242        243        244        245 
#>  -4.504299  -6.113204  -4.403393 -21.564528  -4.113582 -13.595483  -6.663990 
#>        246        247        248        249        250        251        252 
#> -15.946724  -4.254775 -26.205819  -6.977809  -4.696165  -5.512515  -8.670659 
#>        253        254        255        256        257        258        259 
#> -10.980411  -4.761930  -4.379472  -5.854158 -10.457894 -10.903596  -8.402074 
#>        260        261        262        263        264        265        266 
#>  -4.369433  -4.794512  -6.834050  -6.691491 -15.167920  -6.513032  -6.754408 
#>        267        268        269        270        271        272        273 
#>  -6.094088 -21.070727 -19.046839  -8.256524  -4.761618 -27.479445 -43.050347 
#>        274        275        276        277        278        279        280 
#> -10.008444  -5.791209  -8.147274 -16.712370  -4.356783 -17.768102  -6.208954 
#>        281        282        283        284        285        286        287 
#> -13.991192 -19.411188 -21.483631  -4.081200  -4.213080 -43.131797 -18.057646 
#>        288        289        290        291        292        293        294 
#> -22.560500  -4.346636 -17.934602  -9.454650 -11.018643 -21.950537  -5.720397 
#>        295        296        297        298        299        300        301 
#>  -7.864779 -22.960938  -4.465751 -16.873918  -7.333854 -15.125827  -5.576326 
#>        302        303        304        305        306        307        308 
#> -23.463821  -4.656249  -5.861642  -4.114459  -4.525957 -37.019695  -5.745356 
#>        309        310        311        312        313        314        315 
#>  -4.251819 -28.584794 -12.486891 -29.241154 -10.702439  -4.357488  -4.326050 
#>        316        317        318        319        320        321        322 
#>  -4.240167  -5.449423  -6.822403 -10.114359  -6.980591  -4.855649  -4.462275 
#>        323        324        325        326        327        328        329 
#> -12.718060 -18.244126  -8.433017  -8.983453  -7.264701 -27.258033  -4.433268 
#>        330        331        332        333        334        335        336 
#>  -5.595724  -4.215245 -18.285949 -26.972903 -14.868302  -6.930501 -12.082793 
#>        337        338        339        340        341        342        343 
#> -14.953421 -19.599673  -8.106585  -5.695342 -21.703369  -4.586396 -12.748782 
#>        344        345        346        347        348        349        350 
#>  -3.955999 -16.971795 -16.534197 -29.610012 -16.826868  -5.513721  -6.937369 
#>        351        352        353        354        355        356        357 
#>  -8.680805 -10.166254 -11.578243  -4.163982 -19.910044 -27.248136 -15.752426 
#>        358        359        360        361        362        363        364 
#> -20.500953  -4.108723  -5.362101  -4.605677 -10.187313  -4.427260  -4.733038 
#>        365        366        367        368        369        370        371 
#> -11.077258  -4.469114 -12.863830 -30.380626 -11.802507  -4.564728  -8.907167 
#>        372        373        374        375        376        377        378 
#>  -4.336490  -8.981592  -4.676666 -15.353319 -15.790499 -18.287881  -4.549392 
#>        379        380        381        382        383        384        385 
#>  -4.367500  -6.379122 -13.979059  -9.809787  -7.731674  -9.709499 -11.793193 
#>        386        387        388        389        390        391        392 
#> -13.517956 -11.014441  -4.417755  -6.439535 -10.247766  -9.259853  -4.896456 
#>        393        394        395        396        397        398        399 
#> -19.425236 -20.799539  -4.224040  -5.224484 -20.871493 -24.302692  -4.486746 
#>        400        401        402        403        404        405        406 
#>  -6.671885  -4.364008  -4.453198  -5.858647  -4.548513 -13.924773  -5.371274 
#>        407        408        409        410        411        412        413 
#> -16.381598  -5.336469 -28.237185  -4.300187  -4.449397 -12.718340 -10.111967 
#>        414        415        416        417        418        419        420 
#>  -4.856456  -6.795591  -7.134434  -4.319025 -17.680057  -5.012004 -17.640512 
#>        421        422        423        424        425        426        427 
#> -18.777509  -4.649689 -11.628672  -9.772419 -23.594912  -4.159432  -4.310980 
#>        428        429        430        431        432        433        434 
#>  -7.447002 -16.413943 -12.962150  -6.031866 -19.480648  -4.140881  -8.171842 
#> 
# ->
measure_rps(y = fit$y, ypred = ypred1)
#> 
#> Computed from 4000 draws by 434 observations.
#> 
#>     Estimate    SE
#> rps    10.22  0.35
# }
```
