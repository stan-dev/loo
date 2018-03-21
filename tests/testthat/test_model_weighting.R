library(loo)
options(loo.cores = 2)

context("loo_model_weights")

# generate fake data
set.seed(123)
y<-rnorm(50,0,1)
sd_sim1<- abs(rnorm(500,1.5, 0.1))
sd_sim2<- abs(rnorm(500,1.2, 0.1))
sd_sim3<- abs(rnorm(500,1, 0.05))
log_lik1 <- log_lik2 <- log_lik3 <- matrix(NA, 500, 50)
for(s in 1:500) {
  log_lik1[s,] <- dnorm(y,-1,sd_sim1[s], log=T)
  log_lik2[s,] <-  dnorm(y,0.7,sd_sim2[s], log=T)
  log_lik3[s,] <-  dnorm(y,1,sd_sim3[s], log=T)
}

ll_list <- list(log_lik1, log_lik2,log_lik3)
r_eff_list <- list(rep(0.9,50), rep(0.9,50), rep(0.9,50))

tol <- 0.01 # absoulte tolerance of weights

test_that("loo_model_weights throws correct errors", {
  expect_error(loo_model_weights(log_lik1), "is.list")
  expect_error(loo_model_weights(list(log_lik1)), "At least two models")
  expect_error(loo_model_weights(list(log_lik1, log_lik2[-1, ])), "same dimensions")
  expect_error(loo_model_weights(list(log_lik1, log_lik2, log_lik3[, -1])), "same dimensions")

  expect_error(loo_model_weights(ll_list, r_eff_list = r_eff_list[-1]),
               "one component for each model")

  r_eff_list[[3]] <- rep(0.9, 51)
  expect_error(loo_model_weights(ll_list, r_eff_list = r_eff_list),
               "same length as the number of columns")
})


test_that("loo_model_weights (stacking and pseudo-BMA) gives expected result", {
  w1 <- loo_model_weights(ll_list, method = "stacking", r_eff_list = r_eff_list)
  expect_type(w1,"double")
  expect_s3_class(w1, "stacking_weights")
  expect_length(w1, 3)
  expect_named(w1, paste0("model"  ,c(1:3)))
  expect_equal_to_reference(as.numeric(w1), "model_weights_stacking.rds",
                            tolerance  = tol, scale=1)

  w2 <- loo_model_weights(ll_list, r_eff_list=r_eff_list,
                      method = "pseudobma", BB = TRUE)
  expect_type(w2, "double")
  expect_s3_class(w2, "pseudobma_bb_weights")
  expect_length(w2, 3)
  expect_named(w2, paste0("model", c(1:3)))
  expect_equal_to_reference(as.numeric(w2), "model_weights_pseudobma.rds",
                            tolerance  = tol, scale=1)

  w3 <- loo_model_weights(ll_list, r_eff_list=r_eff_list,
                      method = "pseudobma", BB = FALSE)
  expect_type(w3,"double")
  expect_length(w3, 3)
  expect_named(w3, paste0("model"  ,c(1:3)))
  expect_equal(as.numeric(w3), c(5.365279e-05, 9.999436e-01, 2.707028e-06),
               tolerance  = tol, scale = 1)
})