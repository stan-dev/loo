## Test data for testing print method of `kfold` object

### Case 1: All pareto-k values are good

```{r}
set.seed(123)
dat <- dplyr::tibble(
  x = rnorm(200),
  y = 2 + 1.5 * x + rnorm(200, sd = 1)
)

fit <- brm(y ~ x, data = dat, seed = 42)
kfold1 <- kfold(fit)
saveRDS(kfold, "kfold-calibrated.Rds")
```

### Case 2: Some pareto-k values are problematic

```{r}
data(roaches, package = "rstanarm")
roaches$sqrt_roach1 <- sqrt(roaches$roach1)

fit_p <- brm(y ~ sqrt_roach1 + treatment + senior + offset(log(exposure2)),
                data = roaches,
                family = poisson,
                prior = prior(normal(0,1), class = b),
                refresh = 0)

kfold2 <- kfold(fit_p)
saveRDS(kfold2, "kfold-miscalibrated.Rds")
```