library(brms)
library(loo)
options(brms.backend = "cmdstanr")
options(mc.cores = 4)

fitos <- read.delim("data-raw/winequality-red.csv", sep = ";") |>
  unique() |>
  scale() |>
  as.data.frame() |>
  brm(
    ordered(quality) ~ .,
    family = cumulative("logit"),
    prior = prior(R2D2(mean_R2 = 1 / 3, prec_R2 = 3)),
    data = _,
    seed = 1,
    silent = 2,
    refresh = 0
  )

saveRDS(log_lik(fitos), "touchstone/wine.rds")
