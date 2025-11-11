# Bayesian Stacking and Pseudo-BMA weights using the loo package

## Introduction

This vignette demonstrates the new functionality in **loo** v2.0.0 for
Bayesian stacking and Pseudo-BMA weighting. In this vignette we can’t
provide all of the necessary background on this topic, so we encourage
readers to refer to the paper

- Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018). Using
  stacking to average Bayesian predictive distributions. In Bayesian
  Analysis, :10.1214/17-BA1091.
  [Online](https://projecteuclid.org/euclid.ba/1516093227)

which provides important details on the methods demonstrated in this
vignette. Here we just quote from the abstract of the paper:

> **Abstract**: Bayesian model averaging is flawed in the
> \\\mathcal{M}\\-open setting in which the true data-generating process
> is not one of the candidate models being fit. We take the idea of
> stacking from the point estimation literature and generalize to the
> combination of predictive distributions. We extend the utility
> function to any proper scoring rule and use Pareto smoothed importance
> sampling to efficiently compute the required leave-one-out posterior
> distributions. We compare stacking of predictive distributions to
> several alternatives: stacking of means, Bayesian model averaging
> (BMA), Pseudo-BMA, and a variant of Pseudo-BMA that is stabilized
> using the Bayesian bootstrap. Based on simulations and real-data
> applications, we recommend stacking of predictive distributions, with
> bootstrapped-Pseudo-BMA as an approximate alternative when computation
> cost is an issue.

Ideally, we would avoid the Bayesian model combination problem by
extending the model to include the separate models as special cases, and
preferably as a continuous expansion of the model space. For example,
instead of model averaging over different covariate combinations, all
potentially relevant covariates should be included in a predictive model
(for causal analysis more care is needed) and a prior assumption that
only some of the covariates are relevant can be presented with
regularized horseshoe prior (Piironen and Vehtari, 2017a). For variable
selection we recommend projective predictive variable selection
(Piironen and Vehtari, 2017a; [**projpred**
package](https://cran.r-project.org/package=projpred)).

To demonstrate how to use **loo** package to compute Bayesian stacking
and Pseudo-BMA weights, we repeat two simple model averaging examples
from Chapters 6 and 10 of *Statistical Rethinking* by Richard McElreath.
In *Statistical Rethinking* WAIC is used to form weights which are
similar to classical “Akaike weights”. Pseudo-BMA weighting using
PSIS-LOO for computation is close to these WAIC weights, but named after
the Pseudo Bayes Factor by Geisser and Eddy (1979). As discussed below,
in general we prefer using stacking rather than WAIC weights or the
similar pseudo-BMA weights.

## Setup

In addition to the **loo** package we will also load the **rstanarm**
package for fitting the models.

``` r
library(rstanarm)
library(loo)
```

## Example: Primate milk

In *Statistical Rethinking*, McElreath describes the data for the
primate milk example as follows:

> A popular hypothesis has it that primates with larger brains produce
> more energetic milk, so that brains can grow quickly. … The question
> here is to what extent energy content of milk, measured here by
> kilocalories, is related to the percent of the brain mass that is
> neocortex. … We’ll end up needing female body mass as well, to see the
> masking that hides the relationships among the variables.

``` r
data(milk)
d <- milk[complete.cases(milk),]
d$neocortex <- d$neocortex.perc /100
str(d)
```

    'data.frame':   17 obs. of  9 variables:
     $ clade         : Factor w/ 4 levels "Ape","New World Monkey",..: 4 2 2 2 2 2 2 2 3 3 ...
     $ species       : Factor w/ 29 levels "A palliata","Alouatta seniculus",..: 11 2 1 6 27 5 3 4 21 19 ...
     $ kcal.per.g    : num  0.49 0.47 0.56 0.89 0.92 0.8 0.46 0.71 0.68 0.97 ...
     $ perc.fat      : num  16.6 21.2 29.7 53.4 50.6 ...
     $ perc.protein  : num  15.4 23.6 23.5 15.8 22.3 ...
     $ perc.lactose  : num  68 55.2 46.9 30.8 27.1 ...
     $ mass          : num  1.95 5.25 5.37 2.51 0.68 0.12 0.47 0.32 1.55 3.24 ...
     $ neocortex.perc: num  55.2 64.5 64.5 67.6 68.8 ...
     $ neocortex     : num  0.552 0.645 0.645 0.676 0.688 ...

We repeat the analysis in Chapter 6 of *Statistical Rethinking* using
the following four models (here we use the default weakly informative
priors in **rstanarm**, while flat priors were used in *Statistical
Rethinking*).

``` r
fit1 <- stan_glm(kcal.per.g ~ 1, data = d, seed = 2030)
fit2 <- update(fit1, formula = kcal.per.g ~ neocortex)
fit3 <- update(fit1, formula = kcal.per.g ~ log(mass))
fit4 <- update(fit1, formula = kcal.per.g ~ neocortex + log(mass))
```

McElreath uses WAIC for model comparison and averaging, so we’ll start
by also computing WAIC for these models so we can compare the results to
the other options presented later in the vignette. The **loo** package
provides `waic` methods for log-likelihood arrays, matrices and
functions. Since we fit our model with rstanarm we can use the `waic`
method provided by the **rstanarm** package (a wrapper around `waic`
from the **loo** package), which allows us to just pass in our fitted
model objects instead of first extracting the log-likelihood values.

``` r
waic1 <- waic(fit1)
waic2 <- waic(fit2)
waic3 <- waic(fit3)
```

    Warning: 
    1 (5.9%) p_waic estimates greater than 0.4. We recommend trying loo instead.

``` r
waic4 <- waic(fit4)
```

    Warning: 
    2 (11.8%) p_waic estimates greater than 0.4. We recommend trying loo instead.

``` r
waics <- c(
  waic1$estimates["elpd_waic", 1],
  waic2$estimates["elpd_waic", 1],
  waic3$estimates["elpd_waic", 1],
  waic4$estimates["elpd_waic", 1]
)
```

We get some warnings when computing WAIC for models 3 and 4, indicating
that we shouldn’t trust the WAIC weights we will compute later.
Following the recommendation in the warning, we next use the `loo`
methods to compute PSIS-LOO instead. The **loo** package provides `loo`
methods for log-likelihood arrays, matrices, and functions, but since we
fit our model with **rstanarm** we can just pass the fitted model
objects directly and **rstanarm** will extract the needed values to pass
to the **loo** package. (Like **rstanarm**, some other R packages for
fitting Stan models, e.g. **brms**, also provide similar methods for
interfacing with the **loo** package.)

``` r
# note: the loo function accepts a 'cores' argument that we recommend specifying
# when working with bigger datasets

loo1 <- loo(fit1)
loo2 <- loo(fit2)
loo3 <- loo(fit3)
loo4 <- loo(fit4)
lpd_point <- cbind(
  loo1$pointwise[,"elpd_loo"], 
  loo2$pointwise[,"elpd_loo"],
  loo3$pointwise[,"elpd_loo"], 
  loo4$pointwise[,"elpd_loo"]
)
```

With `loo` we don’t get any warnings for models 3 and 4, but for
illustration of good results, we display the diagnostic details for
these models anyway.

``` r
print(loo3)
```

    Computed from 4000 by 17 log-likelihood matrix.

             Estimate  SE
    elpd_loo      4.5 2.3
    p_loo         2.1 0.5
    looic        -9.1 4.6
    ------
    MCSE of elpd_loo is 0.0.
    MCSE and ESS estimates assume independent draws (r_eff=1).

    All Pareto k estimates are good (k < 0.7).
    See help('pareto-k-diagnostic') for details.

``` r
print(loo4)
```

    Computed from 4000 by 17 log-likelihood matrix.

             Estimate  SE
    elpd_loo      8.4 2.8
    p_loo         3.3 0.9
    looic       -16.8 5.5
    ------
    MCSE of elpd_loo is 0.1.
    MCSE and ESS estimates assume independent draws (r_eff=1).

    All Pareto k estimates are good (k < 0.7).
    See help('pareto-k-diagnostic') for details.

One benefit of PSIS-LOO over WAIC is better diagnostics. Here for both
models 3 and 4 all \\k\<0.7\\ and the Monte Carlo SE of `elpd_loo` is
0.1 or less, and we can expect the model comparison to be reliable.

Next we compute and compare 1) WAIC weights, 2) Pseudo-BMA weights
without Bayesian bootstrap, 3) Pseudo-BMA+ weights with Bayesian
bootstrap, and 4) Bayesian stacking weights.

``` r
waic_wts <- exp(waics) / sum(exp(waics))
pbma_wts <- pseudobma_weights(lpd_point, BB=FALSE)
pbma_BB_wts <- pseudobma_weights(lpd_point) # default is BB=TRUE
stacking_wts <- stacking_weights(lpd_point)
round(cbind(waic_wts, pbma_wts, pbma_BB_wts, stacking_wts), 2)
```

           waic_wts pbma_wts pbma_BB_wts stacking_wts
    model1     0.01     0.02        0.07         0.01
    model2     0.01     0.01        0.04         0.00
    model3     0.02     0.02        0.04         0.00
    model4     0.96     0.96        0.86         0.99

With all approaches Model 4 with `neocortex` and `log(mass)` gets most
of the weight. Based on theory, Pseudo-BMA weights without Bayesian
bootstrap should be close to WAIC weights, and we can also see that
here. Pseudo-BMA+ weights with Bayesian bootstrap provide more cautious
weights further away from 0 and 1 (see Yao et al. (2018) for a
discussion of why this can be beneficial and results from related
experiments). In this particular example, the Bayesian stacking weights
are not much different from the other weights.

One of the benefits of stacking is that it manages well if there are
many similar models. Consider for example that there could be many
irrelevant covariates that when included would produce a similar model
to one of the existing models. To emulate this situation here we simply
copy the first model a bunch of times, but you can imagine that instead
we would have ten alternative models with about the same predictive
performance. WAIC weights for such a scenario would be close to the
following:

``` r
waic_wts_demo <- 
  exp(waics[c(1,1,1,1,1,1,1,1,1,1,2,3,4)]) /
  sum(exp(waics[c(1,1,1,1,1,1,1,1,1,1,2,3,4)]))
round(waic_wts_demo, 3)
```

     [1] 0.013 0.013 0.013 0.013 0.013 0.013 0.013 0.013 0.013 0.013 0.006 0.016
    [13] 0.847

Notice how much the weight for model 4 is lowered now that more models
similar to model 1 (or in this case identical) have been added. Both
WAIC weights and Pseudo-BMA approaches first estimate the predictive
performance separately for each model and then compute weights based on
estimated relative predictive performances. Similar models share similar
weights so the weights of other models must be reduced for the total sum
of the weights to remain the same.

On the other hand, stacking optimizes the weights *jointly*, allowing
for the very similar models (in this toy example repeated models) to
share their weight while more unique models keep their original weights.
In our example we can see this difference clearly:

``` r
stacking_weights(lpd_point[,c(1,1,1,1,1,1,1,1,1,1,2,3,4)])
```

    Method: stacking
    ------
            weight
    model1  0.001 
    model2  0.001 
    model3  0.001 
    model4  0.001 
    model5  0.001 
    model6  0.001 
    model7  0.001 
    model8  0.001 
    model9  0.001 
    model10 0.001 
    model11 0.000 
    model12 0.000 
    model13 0.987 

Using stacking, the weight for the best model stays essentially
unchanged.

## Example: Oceanic tool complexity

Another example we consider is the Kline oceanic tool complexity data,
which McElreath describes as follows:

> Different historical island populations possessed tool kits of
> different size. These kits include fish hooks, axes, boats, hand
> plows, and many other types of tools. A number of theories predict
> that larger populations will both develop and sustain more complex
> tool kits. … It’s also suggested that contact rates among populations
> effectively increases population \[sic, probably should be tool kit\]
> size, as it’s relevant to technological evolution.

We build models predicting the total number of tools given the log
population size and the contact rate (high vs. low).

``` r
data(Kline)
d <- Kline
d$log_pop <- log(d$population)
d$contact_high <- ifelse(d$contact=="high", 1, 0)
str(d)
```

    'data.frame':   10 obs. of  7 variables:
     $ culture     : Factor w/ 10 levels "Chuuk","Hawaii",..: 4 7 6 10 3 9 1 5 8 2
     $ population  : int  1100 1500 3600 4791 7400 8000 9200 13000 17500 275000
     $ contact     : Factor w/ 2 levels "high","low": 2 2 2 1 1 1 1 2 1 2
     $ total_tools : int  13 22 24 43 33 19 40 28 55 71
     $ mean_TU     : num  3.2 4.7 4 5 5 4 3.8 6.6 5.4 6.6
     $ log_pop     : num  7 7.31 8.19 8.47 8.91 ...
     $ contact_high: num  0 0 0 1 1 1 1 0 1 0

We start with a Poisson regression model with the log population size,
the contact rate, and an interaction term between them (priors are
informative priors as in *Statistical Rethinking*).

``` r
fit10 <-
  stan_glm(
    total_tools ~ log_pop + contact_high + log_pop * contact_high,
    family = poisson(link = "log"),
    data = d,
    prior = normal(0, 1, autoscale = FALSE),
    prior_intercept = normal(0, 100, autoscale = FALSE),
    seed = 2030
  )
```

Before running other models, we check whether Poisson is good choice as
the conditional observation model.

``` r
loo10 <- loo(fit10)
```

    Warning: Found 2 observation(s) with a pareto_k > 0.7. We recommend calling 'loo' again with argument 'k_threshold = 0.7' in order to calculate the ELPD without the assumption that these observations are negligible. This will refit the model 2 times to compute the ELPDs for the problematic observations directly.

``` r
print(loo10)
```

    Computed from 4000 by 10 log-likelihood matrix.

             Estimate   SE
    elpd_loo    -40.9  6.1
    p_loo         5.7  2.0
    looic        81.8 12.2
    ------
    MCSE of elpd_loo is NA.
    MCSE and ESS estimates assume independent draws (r_eff=1).

    Pareto k diagnostic values:
                             Count Pct.    Min. ESS
    (-Inf, 0.7]   (good)     8     80.0%   288     
       (0.7, 1]   (bad)      2     20.0%   <NA>    
       (1, Inf)   (very bad) 0      0.0%   <NA>    
    See help('pareto-k-diagnostic') for details.

We get at least one observation with \\k\>0.7\\ and the estimated
effective number of parameters `p_loo` is larger than the total number
of parameters in the model. This indicates that Poisson might be too
narrow. A negative binomial model might be better, but with so few
observations it is not so clear.

We can compute LOO more accurately by running Stan again for the
leave-one-out folds with high \\k\\ estimates. When using **rstanarm**
this can be done by specifying the `k_threshold` argument:

``` r
loo10 <- loo(fit10, k_threshold=0.7)
```

    2 problematic observation(s) found.
    Model will be refit 2 times.

    Fitting model 1 out of 2 (leaving out observation 6)

    Fitting model 2 out of 2 (leaving out observation 10)

``` r
print(loo10)
```

    Computed from 4000 by 10 log-likelihood matrix.

             Estimate   SE
    elpd_loo    -41.1  6.0
    p_loo         5.9  2.0
    looic        82.2 12.1
    ------
    MCSE of elpd_loo is 0.2.
    MCSE and ESS estimates assume independent draws (r_eff=1).

    All Pareto k estimates are good (k < 0.7).
    See help('pareto-k-diagnostic') for details.

In this case we see that there is not much difference, and thus it is
relatively safe to continue.

As a comparison we also compute WAIC:

``` r
waic10 <- waic(fit10)
```

    Warning: 
    4 (40.0%) p_waic estimates greater than 0.4. We recommend trying loo instead.

``` r
print(waic10)
```

    Computed from 4000 by 10 log-likelihood matrix.

              Estimate   SE
    elpd_waic    -40.2  6.0
    p_waic         5.0  1.8
    waic          80.4 12.0

    4 (40.0%) p_waic estimates greater than 0.4. We recommend trying loo instead. 

The WAIC computation is giving warnings and the estimated ELPD is
slightly more optimistic. We recommend using the PSIS-LOO results
instead.

To assess whether the contact rate and interaction term are useful, we
can make a comparison to models without these terms.

``` r
fit11 <- update(fit10, formula = total_tools ~ log_pop + contact_high)
fit12 <- update(fit10, formula = total_tools ~ log_pop)
```

``` r
(loo11 <- loo(fit11))
```

    Computed from 4000 by 10 log-likelihood matrix.

             Estimate   SE
    elpd_loo    -39.7  5.8
    p_loo         4.4  1.6
    looic        79.4 11.6
    ------
    MCSE of elpd_loo is 0.1.
    MCSE and ESS estimates assume independent draws (r_eff=1).

    All Pareto k estimates are good (k < 0.7).
    See help('pareto-k-diagnostic') for details.

``` r
(loo12 <- loo(fit12))
```

    Warning: Found 1 observation(s) with a pareto_k > 0.7. We recommend calling 'loo' again with argument 'k_threshold = 0.7' in order to calculate the ELPD without the assumption that these observations are negligible. This will refit the model 1 times to compute the ELPDs for the problematic observations directly.

    Computed from 4000 by 10 log-likelihood matrix.

             Estimate  SE
    elpd_loo    -42.5 4.7
    p_loo         4.0 1.1
    looic        85.0 9.4
    ------
    MCSE of elpd_loo is NA.
    MCSE and ESS estimates assume independent draws (r_eff=1).

    Pareto k diagnostic values:
                             Count Pct.    Min. ESS
    (-Inf, 0.7]   (good)     9     90.0%   1251    
       (0.7, 1]   (bad)      1     10.0%   <NA>    
       (1, Inf)   (very bad) 0      0.0%   <NA>    
    See help('pareto-k-diagnostic') for details.

``` r
loo11 <- loo(fit11, k_threshold=0.7)
```

    All pareto_k estimates below user-specified threshold of 0.7. 
    Returning loo object.

``` r
loo12 <- loo(fit12, k_threshold=0.7)
```

    1 problematic observation(s) found.
    Model will be refit 1 times.

    Fitting model 1 out of 1 (leaving out observation 10)

``` r
lpd_point <- cbind(
  loo10$pointwise[, "elpd_loo"], 
  loo11$pointwise[, "elpd_loo"], 
  loo12$pointwise[, "elpd_loo"]
)
```

For comparison we’ll also compute WAIC values for these additional
models:

``` r
waic11 <- waic(fit11)
```

    Warning: 
    3 (30.0%) p_waic estimates greater than 0.4. We recommend trying loo instead.

``` r
waic12 <- waic(fit12)
```

    Warning: 
    5 (50.0%) p_waic estimates greater than 0.4. We recommend trying loo instead.

``` r
waics <- c(
  waic10$estimates["elpd_waic", 1], 
  waic11$estimates["elpd_waic", 1], 
  waic12$estimates["elpd_waic", 1]
)
```

The WAIC computation again gives warnings, and we recommend using
PSIS-LOO instead.

Finally, we compute 1) WAIC weights, 2) Pseudo-BMA weights without
Bayesian bootstrap, 3) Pseudo-BMA+ weights with Bayesian bootstrap, and
4) Bayesian stacking weights.

``` r
waic_wts <- exp(waics) / sum(exp(waics))
pbma_wts <- pseudobma_weights(lpd_point, BB=FALSE)
pbma_BB_wts <- pseudobma_weights(lpd_point) # default is BB=TRUE
stacking_wts <- stacking_weights(lpd_point)
round(cbind(waic_wts, pbma_wts, pbma_BB_wts, stacking_wts), 2)
```

           waic_wts pbma_wts pbma_BB_wts stacking_wts
    model1     0.32     0.19        0.18          0.0
    model2     0.64     0.79        0.69          0.8
    model3     0.04     0.02        0.14          0.2

All weights favor the second model with the log population and the
contact rate. WAIC weights and Pseudo-BMA weights (without Bayesian
bootstrap) are similar, while Pseudo-BMA+ is more cautious and closer to
stacking weights.

It may seem surprising that Bayesian stacking is giving zero weight to
the first model, but this is likely due to the fact that the estimated
effect for the interaction term is close to zero and thus models 1 and 2
give very similar predictions. In other words, incorporating the model
with the interaction (model 1) into the model average doesn’t improve
the predictions at all and so model 1 is given a weight of 0. On the
other hand, models 2 and 3 are giving slightly different predictions and
thus their combination may be slightly better than either alone. This
behavior is related to the repeated similar model illustration in the
milk example above.

## Simpler coding using `loo_model_weights` function

Although in the examples above we called the `stacking_weights` and
`pseudobma_weights` functions directly, we can also use the
`loo_model_weights` wrapper, which takes as its input either a list of
pointwise log-likelihood matrices or a list of precomputed loo objects.
There are also `loo_model_weights` methods for stanreg objects (fitted
model objects from **rstanarm**) as well as fitted model objects from
other packages (e.g. **brms**) that do the preparation work for the user
(see, e.g., the examples at
[`help("loo_model_weights", package = "rstanarm")`](https://mc-stan.org/loo/reference/loo_model_weights.html)).

``` r
# using list of loo objects
loo_list <- list(loo10, loo11, loo12)
loo_model_weights(loo_list)
```

    Method: stacking
    ------
          weight
    fit10 0.000 
    fit11 0.803 
    fit12 0.197 

``` r
loo_model_weights(loo_list, method = "pseudobma")
```

    Method: pseudo-BMA+ with Bayesian bootstrap
    ------
          weight
    fit10 0.169 
    fit11 0.644 
    fit12 0.187 

``` r
loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)
```

    Method: pseudo-BMA
    ------
          weight
    fit10 0.189 
    fit11 0.792 
    fit12 0.019 

## References

McElreath, R. (2016). *Statistical rethinking: A Bayesian course with
examples in R and Stan*. Chapman & Hall/CRC.
<http://xcelab.net/rm/statistical-rethinking/>

Piironen, J. and Vehtari, A. (2017a). Sparsity information and
regularization in the horseshoe and other shrinkage priors. In
Electronic Journal of Statistics, 11(2):5018-5051.
[Online](https://projecteuclid.org/euclid.ejs/1513306866).

Piironen, J. and Vehtari, A. (2017b). Comparison of Bayesian predictive
methods for model selection. Statistics and Computing, 27(3):711-735.
:10.1007/s11222-016-9649-y.
[Online](https://link.springer.com/article/10.1007/s11222-016-9649-y).

Vehtari, A., Gelman, A., and Gabry, J. (2017). Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*. 27(5), 1413–1432. :10.1007/s11222-016-9696-4.
[online](https://link.springer.com/article/10.1007/s11222-016-9696-4),
[arXiv preprint arXiv:1507.04544](https://arxiv.org/abs/1507.04544).

Vehtari, A., Simpson, D., Gelman, A., Yao, Y., and Gabry, J. (2024).
Pareto smoothed importance sampling. *Journal of Machine Learning
Research*, 25(72):1-58. [PDF](https://jmlr.org/papers/v25/19-556.html)

Yao, Y., Vehtari, A., Simpson, D., and Gelman, A. (2018). Using stacking
to average Bayesian predictive distributions. In Bayesian Analysis,
:10.1214/17-BA1091.
[Online](https://projecteuclid.org/euclid.ba/1516093227).
