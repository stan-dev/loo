# This creates the file sysdata.rda in the 'R' directory
# Users are given access via example_loglik() functions, see loo/R/example_loglik_objects.R


# .example_loglik_array
library(rstanarm)
fit <- stan_glm(mpg ~ wt, data = mtcars, chains = 2, iter = 1000, warmup = 500)
ll <- log_lik(fit)
.example_loglik_array <- loo:::llmatrix_to_array(ll, chain_id = rep(1:2, each = 500))


# .example_wine_loglik_matrix
library(dplyr)
library(brms)
options(brms.backend = "cmdstanr")
options(mc.cores = 4)

wine <- read.delim("data-raw/winequality-red.csv", sep = ";") |> distinct()

wine_scaled <- as.data.frame(scale(wine))

fitos <- brm(ordered(quality) ~ .,
             family = cumulative("logit"),
             prior = prior(R2D2(mean_R2 = 1/3, prec_R2 = 3)),
             data = wine_scaled,
             seed = 1,
             silent = 2,
             refresh = 0)

.example_wine_loglik_matrix <- log_lik(fitos)


# Store in R/sysdata.rda
devtools::use_data(
  .example_loglik_array,
  .example_wine_loglik_matrix,
  internal = TRUE,
  overwrite = TRUE
)
