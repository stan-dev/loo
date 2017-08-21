# this creates the file sysdata.rda in the 'R' directory
# users are given access via example_loglik()
library(rstanarm)
fit <- stan_glm(mpg ~ wt, data = mtcars, chains = 2, iter = 1000, warmup = 500)
ll <- log_lik(fit)
.example_loglik_array <- loo:::llmatrix_to_array(ll, chain_id = rep(1:2, each = 500))
devtools::use_data(.example_loglik_array, internal = TRUE, overwrite = TRUE)
