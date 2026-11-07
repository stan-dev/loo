args <- commandArgs(TRUE)
book <- args[1]
out <- args[2]
suppressMessages(library(rstanarm))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
as_array <- function(fit) {
  ll <- log_lik(fit)
  n_chain <- length(fit$stanfit@sim$samples)
  array(ll, c(nrow(ll) / n_chain, n_chain, ncol(ll)))
}
save_ll <- function(name, fit) saveRDS(as_array(fit), file.path(out, paste0(name, ".rds")))
wells <- read.csv(file.path(book, "misc/chapter_08/section_08_04/wells.csv"))
save_ll("arsenic_3", stan_glm(switch ~ dist100 + arsenic, family = binomial, data = wells, refresh = 0, seed = 1))
save_ll("arsenic_3l", stan_glm(switch ~ dist100 + log10(arsenic), family = binomial, data = wells, refresh = 0, seed = 1))
data(roaches, package = "rstanarm")
roaches$sqrt_roach1 <- sqrt(roaches$roach1)
f <- y ~ sqrt_roach1 + treatment + senior + offset(log(exposure2))
save_ll("roaches_poisson", stan_glm(f, family = poisson, data = roaches, refresh = 0, seed = 1))
save_ll("roaches_nb", stan_glm(f, family = neg_binomial_2, data = roaches, refresh = 0, seed = 1))
save_ll("roaches_poisson_vi", stan_glmer(y ~ sqrt_roach1 + treatment + senior + offset(log(exposure2)) + (1 | id), family = poisson, data = transform(roaches, id = seq_len(nrow(roaches))), refresh = 0, seed = 1))
