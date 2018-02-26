# Unit test for model_weights() and model_select()
library(loo)
options(loo.cores = 2)
context("stacking")
# generate fake data
set.seed(123)
y<-rnorm(50,0,1)
sd_sim1<- abs(rnorm(500,1.5, 0.1))
sd_sim2<- abs(rnorm(500,1.2, 0.1))
sd_sim3<- abs(rnorm(500,1, 0.05))
log_lik1 <- log_lik2<-log_lik3<- matrix(NA, 500, 50)
for( s in 1:500){
  log_lik1[s,] <- dnorm(y,-1,sd_sim1[s], log=T)
  log_lik2[s,] <-  dnorm(y,0.7,sd_sim2[s], log=T)
  log_lik3[s,] <-  dnorm(y,1,sd_sim3[s], log=T)
}

tol<-0.01 ## absoulte tolerance of weights
test_that("model_weights (stacking and pseudo-BMA) gives expected result", {
  w1 <- model_weights(list(log_lik1, log_lik2,log_lik3),method="stacking", r_eff_list=list(rep(0.9,50),rep(0.9,50),rep(0.9,50)))
  expect_type(w1,"double" )
  expect_length(w1,3 )
  expect_named(w1, paste("Model"  ,c(1:3)))
  expect_equal_to_reference(unname(w1), "model_weights_stacking.rds", tolerance  = tol, scale=1)

  w2 <- model_weights(list(log_lik1, log_lik2,log_lik3),r_eff_list=list(rep(0.9,50),rep(0.9,50),rep(0.9,50)),
                      method="pseudobma",BB=T,seed=123)
  expect_type(w2,"double" )
  expect_length(w2,3 )
  expect_named(w2, paste("Model"  ,c(1:3)))

  expect_equal_to_reference(unname(w2), "model_weights_pseudobma.rds", tolerance  = tol, scale=1)

  w3 <- model_weights(list(log_lik1, log_lik2,log_lik3),r_eff_list=list(rep(0.9,50),rep(0.9,50),rep(0.9,50)),
                      method="pseudobma",BB=F)
  expect_type(w3,"double" )
  expect_length(w3,3 )
  expect_named(w3, paste("Model"  ,c(1:3)))
  expect_equal(unname(w3), c(5.365279e-05, 9.999436e-01, 2.707028e-06), tolerance  = tol, scale=1)
})

test_that("model_select  gives expected result", {
select_prob1=model_select(list(log_lik1, log_lik2,log_lik3),r_eff_list=list(rep(0.9,50),rep(0.9,50),rep(0.9,50)),  visualise = FALSE, seed=123, BB=F)
expect_type(select_prob1,"integer")
expect_length(select_prob1,1 )
expect_equal(select_prob1,2)
select_prob2<-model_select(list(log_lik1, log_lik2,log_lik3),r_eff_list=list(rep(0.9,50),rep(0.9,50),rep(0.9,50)),  visualise = FALSE, seed=123, BB=TRUE)
expect_type(select_prob2,"double")
expect_length(select_prob2,3 )
expect_named(select_prob2, paste("Model"  ,c(1:3)))
expect_equal(unname(select_prob2),c(0.054, 0.946, 0.000))

})
#saveRDS(w1,"model_weights_stacking.rds")
#saveRDS(w2,"model_weights_pseudobma.rds")
