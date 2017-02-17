##### Function for Stacking and Loo Weighting
library(loo)
library(Rsolnp)    ## package used for optimization
library(MCMCpack)  ## package used for generating dirichlet distribution  

model_average=function(N, model_list,   BB_sample=1000)  ### N: number of sample points. model_list: list of stan fit for different models
{
  ## sum to 1 condition ###
  constrain=function(w){
    sum(w)-1
  }
  ## stacking on log score ###
  stacking_score=function(lpd_point){
    log_score_loo=function(w)
    {
      sum=0
      for(i in 1:N)
        sum=sum+log(exp(lpd_point[i,])%*%w)
      return(as.numeric(sum))
    }
    negative_log_score_loo=function(w)
    {
      return(-log_score_loo(w))
    }
    return( solnp(pars =rep(1/K,K),fun=negative_log_score_loo,  control=list( trace=0),
                  eqfun=constrain,eqB = 0, LB=rep(0,K),UB=rep(1,K) ) $pars)
  }
  
  K=length(model_list)
  lpd_point=matrix(NA,N,K)            #point wise log likelihood
  elpd_loo=rep(NA,K)
  for( k in 1:K){
    log_likelihood=extract( model_list[[k]] )[['log_lik']]
    L=loo(log_likelihood)
    lpd_point[,k] = L$pointwise[,1]    ## log(p_k (y_i | y_-i))
    elpd_loo[k]=L$elpd_loo
  }
  w_stacking = stacking_score(lpd_point)
  
  ## 2) loo weights
  uwts <- exp( elpd_loo - max( elpd_loo))
  w_loo1 <- uwts / sum(uwts)            ## original loo weights
  ## 3) loo withs using BB sample
  temp=matrix(NA,BB_sample,K)
  BB_weighting=rdirichlet(BB_sample, rep(1,N))
  for( bb in 1:BB_sample){
    z_bb=   BB_weighting[bb,] %*% lpd_point *N
    uwts  <- exp( z_bb - max(z_bb  ))
    temp[bb,] <- uwts / sum(uwts)
  }
  w_loo2  <- colMeans(temp)
  ## 4) best single model selection using loo
  k_best=which(elpd_loo==max(elpd_loo))
  w_selection=rep(0,K)
  w_selection[k_best]=1
  cat("The stacking weights are:\n")
  print(w_stacking)
  cat("The LOO weights are:\n")
  print(w_loo1)
  cat("The LOO weights using BB sample are:\n")
  print(w_loo2)
  cat("The best model selection is:\n")
  print(paste("Model ",k_best, ""))
  return(matrix( c(w_stacking,w_loo1,w_loo2,w_selection),N,4)   )
}
 

#### examples:###########
y=c(4,9.5)
x=c(2,3)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
#model 1:
stan_code_1<-"  
data{
int N;
vector[N] y;
vector[N] x;
}
parameters{
real beta;
real<lower=0> sigma;
}
transformed parameters{
vector[N] theta;
theta= x * beta ;
}
model{
y~normal(theta, sigma);          
beta~normal(0, 2);
}
generated quantities {
vector[N] log_lik; 
for (n in 1:N)
log_lik[n] = normal_lpdf(y[n] |theta[n], sigma);
}
"
fit1=stan(model_code = stan_code_1, data=list(N=2,y=y,x=x), iter=2000,  chains =3)

#model 2:
stan_code_2<-"  
data{
int N;
vector[N] y;
vector[N] x;
}
parameters{
real beta;
real<lower=0> sigma;
}
transformed parameters{
vector[N] theta;
for(n in 1:N)
theta[n]= x[n] * x[n] * beta ;
}
model{
y~normal(theta, sigma);          
beta~normal(0,2);
}
generated quantities {
vector[N] log_lik; 
for (n in 1:N)
log_lik[n] = normal_lpdf(y[n] |theta[n], sigma);
}
"
fit2=stan(model_code = stan_code_2, data=list(N=2,y=y,x=x), iter=2000,  chains =3)

Model_list=list(fit1, fit2)
model_average (N=2, model_list=Model_list,   BB_sample=1000)
