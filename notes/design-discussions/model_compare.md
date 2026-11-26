# Design discussion (`model_compare`)

> **Participants:** Florence Bockting, Jonah Gabry
>
> **Last update:** 2026-09-28

## Design of the `rank_by` argument

**Scope**

- accepting both **model name** and **measure name** might be confusing
- `rank_by` might be most reasonable for **measure name**
- perhaps we want to have two arguments:
  - `rank_by` = measure name
  - `reference` = model name

**Behavior**

- Currently: 
  - Best model is always in the first row; followed by second-best, etc.
  - If reference model is not the "best" model then this will look as follows (reference = m3):
  ```r
  #> -- rmse (vs m3) --
  #>  model rmse_diff rmse_se_diff
  #>     m2       0.6          2.1
  #>     m3       0.0          0.0
  #>     m1      -5.5          4.0
  ```

After discussing this aspect, we realized that it makes a lot of things very complex (what should be the rule for ordering the rows? What does it mean to rank by a measure?, etc.)
We were also not sure in how far users might really have use cases where they want to
order according to a specific model or measure. Therefore we decided to remove
the `rank_by` argument initially. If users say they would like to have such an
argument we can add it later.

## Sign flipping of losses

> I think it would be helpful to include the sign flipping message a bit more prominently e.g. (-- mse (vs m2) — utility scale, higher is better --) . Right now only model_compare() prints the explanation so I think printing a saved comparison object later loses that info.

Sounds like a good idea. However, providing this information only for losses might look confusing as also utilities are on a "utility scale" but there we would not show the information. Might be obvious for users that know the difference between utility and loss but less clear otherwise.

The following output-snippet shows the behavior of the current proposed
messaging.

```r
#   -- elpd (vs m2) --
#    model elpd_diff se_diff p_worse diag_diff
#       m2      0.00    0.00      NA          
#       m3    -25.47  129.10    0.58          
#       m1   -850.29  372.31    0.99          
    
#   -- r2 (vs m2) --
#    model r2_diff r2_se_diff
#       m2    0.00       0.00
#       m3   -0.09       0.18
#       m1   -0.10       0.22
    
#   -- mae (vs m3, sign flipped) --
#    model mae_diff mae_se_diff
#       m3     0.00        0.00
#       m2    -0.07        1.24
#       m1    -6.34        3.08
    
#   All differences: 0 = best model, negative = worse.
#   Signs flipped for loss measures: mae.
```

## Custom measures (incl. difference SE estimate)

> measure_name and measure_loss are attributes of the measure function, but the SE method is an argument to model_compare(). I think (although I could be wrong), that means that if a package author wants to include a custom measure in their package, users would have to remember to set custom_se_fn themselves. In other words, it can’t be fully self contained the way it’s currently designed. Is that right, or am I wrong about this? Should we instead use attr(fn, "measure_se_diff") and allow the custom_se_fn to override it?

Thank you for pointing this out. This is indeed a flaw in the design.
I refactored the design such that a custom measure can have now the attribute
`measure_se_diff`. Furthermore, I added an exported wrapper `custom_measure(fun, name, se_diff_fun = NULL, loss = FALSE)` which sets the three attributes `measure_name`, `measure_loss`, and `measure_se_diff`. 

With the declaration in place, the `custom_se_fn` argument of `model_compare()`
was redundant, so it is removed. For a custom measure that declares
nothing, `model_compare()` reports the difference with an `NA` standard error
and a message.

```r
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

h1 <- fit_measure(fit_m1, measure = list("rmse", huber_measure))
h3 <- fit_measure(fit_m3, measure = list("rmse", huber_measure))

comp_h <- model_compare(list(m3 = h3, m1 = h1))
```

## Using `model_compare` with `kfold`

> I think the brms kfold example has a mistake. It uses brms::kfold(fit, K = 5) separately for each model but I think that means they’re using different folds, which means se_diff is wrong? I think we need to do something like this: folds <- loo::kfold_split_random(K = 5, N = nrow(roaches)) and then pass that to kfold.

Yes, indeed. I changed the corresponding cell in the notebook and added a warning (`throw_kfold_folds_mismatch_warning`) when folds are not equal.

## The helper `add_loo()`

> The add_loo() helper is using moment_match = TRUE and r_eff (since brms does). But loo_pred_measure doesn’t. So the displayed loo_compare and model_compare results don’t actually match for loo objects. 

Yes, that's actually a tricky one.
Currently, we accept for `loo_pred_measure` three input schemes:

+ `loo`: able to reproduce loo_moment_match results
+ `ylp` + `psis_object`: the weights are the moment-matched, so it would work for measures where we only use the weights; but for `elpd` is does not work, as it is recomputed from `ylp`
+ `ylp`: nothing from the moment-matching reaches the computation of measures

So, the question is how we want to handle this and pass the diagnostic information to pred_measure

## Printing
### Diagnostic flags

> Some of the print output says Diagnostic flags present but doesn’t actually show any diagnostic flags in the output 

I assume you refer here to the missing p_worse and diag_diff column for measures that are not elpd.

The reason is that I was not sure where the normal approximation is reasonable for these measures as well. We want indeed to include the diagnostic columns here as well but we have to check first whether the normal approximation and thus the diagnostics are valid for all other measures.

I updated the printing method such that the message "diagnostic flags present" is only printed when `elpd` is present.

### Number of digits per measure

> I think the default number of digits to display is tricky. We might need different defaults per measure or use significant digits or something? I’m not sure, but I think it’s going to be annoying/confusing for users. Especially for constrained measures like R2, acc, brier, etc.

Yes, I agree. I updated the digits rule and set different default formatting for the measures:

+ 1 digit: elpd, ic
+ 3 digits: mlpd, r2, acc, bacc, brier
+ dependent on SE: mae, rmse, mse, rps, srps, custom_measure

```
before            after
 model r2_diff    model r2_diff r2_se_diff
    m2     0.0       m2   0.000      0.000
    m3    -0.1       m3  -0.093      0.178
    m1    -0.1       m1  -0.105      0.217
```

## Improving vignette

> There’s a ton of explanation about the computation but not really very much explanation about how to interpret the results 
> We already have this problem in the pre-existing loo vignette I think, but it occurs to me that the example has many bad K values, which means we don’t even recommend trusting it! I wonder how big of a problem it is for a tutorial vignette? Not sure

Is a TODO and refers to [Issue #401](https://github.com/stan-dev/loo/issues/401)