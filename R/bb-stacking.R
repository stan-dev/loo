#' Model averaging via stacking or pseudo-BMA weighting
#'
#' @description Model averaging via stacking or pseudo-BMA weighting. See
#'  Yao et al. (2017) and  Vehtari, Gelman, and Gabry (2017) for
#'   background.
#' @export
#' @param log_lik_list A list of pointwise log likelihood simulation matrixes
#' (\eqn{S} by \eqn{N}), where \eqn{S} is the size of the posterior sample
#' (with all chains merged) and \eqn{N} is the number of data points.
#' The \eqn{i}-th element corresponds to the \eqn{i}-th model.
#' @param method One of  \code{"stacking"} or \code{"pseudobma"}, indicating
#'   which method is to use for obtaining the optimal weights. \code{"stacking"}
#'   refers to stacking of predictive distributions and  \code{"pseudobma"}
#'   refers to pseudo-BMA weighting (by setting \code{"BB"=FALSE}) or pseudo-BMA+
#'   weighting (by setting \code{"BB"=TRUE}).
#' @param BB Logical used when \code{"method"}=\code{"pseudobma"}. If
#'   \code{TRUE}(default), the Bayesian bootstrap will be used to adjust the
#'   pseudo-BMA weighting, which is called pseudo-BMA+ weighting. It helps
#'   regularize the weight away from 0 and 1, so as to reduce the variance.
#' @param BB_n A positive integer indicating the number of samples in
#'   Bayesian bootstrap. It is necessary when  \code{BB}=\code{TRUE}. The
#'   default number is 1000.
#' @param alpha A positive scalar; the shape parameter in the Dirichlet
#'   distribution of Bayesian bootstrap. The default is 1.
#' @param seed An integer; optional. If specified, it will fix the random seed
#'   when dong Bayesian bootstrap sampling.
#' @param optim_method	 The optimization method to be used in stacking. It can
#'   be chosen from "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent".
#'   The default method is "BFGS".
#' @param optim_control A list of control parameters in optimization. See \code{\link{constrOptim}}.
#' @param r_eff_list (Optional)  A list of relative effective sample size estimates
#'  for the likelihood \code{(exp(log_lik))} of each observation in each model. See
#'   \code{\link{psis}} and  \code{\link{relative_eff}} helper function for
#'   computing r_eff. The default is \code{NULL}.
#' @param cores	 The number of cores to use for parallelization. Same as for
#' \code{\link{psis}}.  The default for an entire R session can be set with
#' \code{options(loo.cores = NUMBER)}. As of version 2.0.0 the default is now
#' 1 core, but we recommend using as many (or close to as many) cores as possible.
#'
#' @return A vector of optimal model weights.
#' @details This function implements  stacking of predictive distributions,
#'   pseudo-BMA and pseudo-BMA+ weighting for combining multiple predictive
#'   distributions.
#'
#' For either method, we can use  Leave-one-out cross-validation (LOO) to
#' estimate the expected log predictive density(elpd).  \code{Stacking} combines
#' all model by maximizing the leave-one-out predictive density of the
#' combination distribution. \code{Pseudo-BMA} finds the relative weights
#' proportional to elpd of each model. \code{Pseudo-BMA+} takes into account the
#' uncertainty resulted from having a finite number of proxy samples from the
#' future data distribution through Bayesian bootstrap (set \code{BB=TRUE}),
#' which will keep weights further away from 0 and 1.
#'
#' In general, we recommend  \code{stacking} for averaging predictive
#' distributions, while \code{pseudo-BMA+} can serve as a computationally
#' easier alternative.
#'
#' @seealso
#' \itemize{
#' \item \code{\link{pseudobma_weight}} for details on pseudo-BMA and pseudo-BMA+ weights.
#' \item \code{\link{stacking_weight}} for details on stacking weighs.
#' \item \code{\link{constrOptim}} for the choice of optimization methods and control-parameters.
#' \item \code{\link{relative_eff}}  for  computing r_eff.
#' \item \code{\link{model_select}} for single-model selection.
#' }
#'
#'@template stacking-references
#'
#' @examples
#' \dontrun{
#' ### Usage with stanfit objects.
#' library(rstan)
#' log_lik1 <- extract(stan(model=model_1,data=data))[['log_lik']]
#' log_lik2 <- extract(stan(model=model_2,data=data))[['log_lik']]
#' w1=model_weights(list(log_lik1, log_lik2),method="stacking")
#' w2=model_weights(list(log_lik1, log_lik2),method="pseudobma",BB=TRUE)}
#'
#' \dontrun{
#' ### A top example
#' # generate fake data from N(0,1).
#' set.seed(100)
#' N <- 100
#' y <- rnorm(N, 0, 1)
#' # Suppose we have three models: N(-1, sigma), N(0.5, sigma) and N(0.6,sigma).
#' stan_code <- ' data { int n; vector[n] y; real mu_fixed; }
#' parameters { real<lower=0> sigma;}model {y ~ normal(mu_fixed, sigma);}
#' generated quantities {vector[n] log_lik;
#' for (i in 1:n) log_lik[i] = normal_lpdf(y[i]| mu_fixed, sigma); }'
#' fit_sample_1 <- stan(model_code=stan_code, data=list(n=N, y=y, mu_fixed=-1))
#' fit_sample_2 <- stan(model_code=stan_code, data=list(n=N, y=y, mu_fixed=0.5))
#' fit_sample_3 <- stan(model_code=stan_code, data=list(n=N, y=y, mu_fixed=0.6))
#' log_lik_list <- list(extract( fit_sample_1)[['log_lik']],extract( fit_sample_2)[['log_lik']],
#' extract( fit_sample_3)[['log_lik']])
#' r_eff_list <- list( relative_eff(exp( extract_log_lik(fit_sample_1, merge_chains = FALSE))),
#' relative_eff( exp(extract_log_lik(fit_sample_2, merge_chains = FALSE))),
#' relative_eff( exp(extract_log_lik(fit_sample_3, merge_chains = FALSE))))
#' # Three model-weighting methods:
#' weights1 <- model_weights(log_lik_list,method="stacking", r_eff_list=r_eff_list,
#' optim_control = list(reltol=1e-10))
#' weights2 <- model_weights(log_lik_list,method="pseudobma",BB=FALSE,r_eff_list=r_eff_list)
#' weights3 <- model_weights(log_lik_list,method="pseudobma",BB=TRUE,r_eff_list=r_eff_list)}


model_weights <-function(log_lik_list, method="stacking",BB=TRUE,BB_n=1000, alpha=1, seed=NULL, optim_method="BFGS", optim_control=list(), r_eff_list = NULL, cores=1)
{
  if (!method %in%c("stacking","pseudobma") )
    stop("Must specify a method in stacking or pseudobma .")
  K<-length(log_lik_list)                #number of models
  if (!is.null(r_eff_list) & length(r_eff_list)!=K)
    stop("Dimensions do not match. If specified, r_eff_list should have the same length as log_lik_list. You can set r_eff_list=NULL instead.")
  if (K==1)
    stop("Only one model is found.")
  if(length(unique(unlist(lapply(log_lik_list,ncol))))!=1 |length(unique(unlist(lapply(log_lik_list,nrow))))!=1)
    stop("Dimensions do not match. Each element of log_lik_list should have same dimensions.")
  N<-ncol(log_lik_list[[1]])             #number of data points
  lpd_point<-matrix(NA,N,K)            #point wise log likelihood
  elpd_loo<-rep(NA,K)
  for( k in 1:K){
    if(!is.null(r_eff_list) ){
       r_eff_k <- r_eff_list[[k]]
       if(  length(r_eff_k) != N )
           stop("Dimensions of r_eff_list do not match. Each r_eff should have the same length as the number of columns in the likelihood matrix.")
    }
    log_likelihood<- log_lik_list[[k]]
    L <- loo(log_likelihood, r_eff = r_eff_list[[k]],cores=cores)
    lpd_point[,k] <- L$pointwise[,1]    #calculate log(p_k (y_i | y_-i))
    elpd_loo[k] <- L$estimates["elpd_loo", 1]
  }
  ## 1) stacking on log score
  if (method =="stacking"){
    w_stacking <- stacking_weight(lpd_point, optim_method=optim_method, optim_control=optim_control)
    cat("The stacking weights are:\n")
    print_weight_vector(w_stacking)
    return(setNames( w_stacking,paste("Model"  ,c(1:K))))
  }
  else
    if (method =="pseudobma"){
      uwts <- exp( elpd_loo - max( elpd_loo))
      w_loo1 <- uwts / sum(uwts)
      if(BB==FALSE){
        cat("The pseudo-BMA weights are:\n")
        print_weight_vector(w_loo1)
        return(setNames( w_loo1,paste("Model"  ,c(1:K))))
      }
      if(BB==TRUE){
        w_loo2  <- pseudobma_weight(lpd_point, BB_n,alpha, seed)
        cat("The pseudo-BMA+ weights using Bayesian bootstrap are:\n ")
        print_weight_vector(w_loo2)
        return(setNames( w_loo2,paste("Model"  ,c(1:K))))
      }
    }
}

#' Model selection via Leave-one-out log predictive density.
#' @description Model selection via Leave-one-out log predictive density estimation and Bayesian bootstrap adjustment

#' @importFrom graphics barplot
#' @export
#' @param log_lik_list A list of pointwise log likelihood simulation matrixes.
#'   The \eqn{i}-th element corresponds to the \eqn{i}-th model. Each row of the
#'   matrix is the log likelihood vector evaluated using a simulated parameter.
#' @param BB Logical. If \code{TRUE}(default), the Bayesian bootstrap will be used to adjust the LOO estimator.
#' @param BB_n A positive integer indicating the number of samples in Bayesian
#'   bootstrap. It is necessary  when \code{BB}=\code{True}. The  default number
#'   is 1000.
#' @param alpha A positive scalar; the shape parameter in the Dirichlet
#'   distribution when doing the bootstrap.
#' @param seed An integer; optional. If specified, it will fix the random seed
#'   when dong the Bayesian bootstrap.
#' @param r_eff_list (Optional)  A list of relative effective sample size estimates
#'  for the likelihood (exp(log_lik)) of each observation in each model. See
#'   \code{\link{psis}} and  \code{\link{relative_eff}} helper function for
#'   computing r_eff. The default is \code{NULL}.
#' @param cores	 The number of cores to use for parallelization. Same as for
#' \code{\link{psis}}.  The default for an entire R session can be set with
#' \code{options(loo.cores = NUMBER)}. As of version 2.0.0 the default is now
#' 1 core, but we recommend using as many (or close to as many) cores as possible.
#'
#' @param visualise Logical, whether to visualise the selection probability.
#' @return When \code{BB=FALSE}, it returns an integer indicating the index of
#'   the best model.  When\code{BB=TRUE}, it returns a vector indicating the
#'   probability of each model being selected to be the best model.
#' @details \code{\link{loo}} gives an estimation of the expected log predictive
#'   density of each model, we can pick the best model by picking the model with
#'   the largest elpd estimation. Just like \code{\link{pseudobma_weight}}, to
#'   make the elpd estimation more reliable, we can use the Bayesian bootstrap
#'   adjustment. With each sample in the Bayesian bootstrap, we compare the
#'   adjusted elpd estimation and finally compute the probability of that model
#'   being the optimal one. If none of the probability is close to 1, then it is
#'   better to do model averaging rather than model selection.
#' @seealso
#' \code{\link{model_weights}} for model combination.
#'
#' \code{\link{compare}} for two-model comparison.
#'
#' \code{\link{pseudobma_weight}} for details on pseudo-BMA weighting.
#'

#' @examples
#' \dontrun{
#' ### Usage with stanfit objects
#' library(rstan)
#' log_lik1 <- extract(stan(model=model_1,data=data))[['log_lik']]
#' log_lik2 <- extract(stan(model=model_1,data=data))[['log_lik']]
#' select_prob <- model_select(list(log_lik1,log_lik2),BB=TRUE)}
#'
model_select <-function(log_lik_list, BB=TRUE,BB_n=1000, alpha=1,seed=NULL,visualise=FALSE,
                        r_eff_list=NULL, cores=1)
{
  if (!is.logical(BB))
    stop("BB must be logical.")
  K<-length(log_lik_list)                #number of models
  if (!is.null(r_eff_list) & length(r_eff_list)!=K)
    stop("Dimensions do not match. If specified, r_eff_list should have the same length as log_lik_list. You can set r_eff_list=NULL instead.")

  if (K==1)
    stop("Only one model is found.")
  N<-ncol(log_lik_list[[1]])             #number of data points
  lpd_point<-matrix(NA,N,K)            #point wise log likelihood
  elpd_loo<-rep(NA,K)
  for( k in 1:K){
    if(!is.null(r_eff_list) ){
      r_eff_k=r_eff_list[[k]]
      if(  length(r_eff_k) != N )
        stop("Dimensions of r_eff_list do not match. Each r_eff should have the same length as the number of columns in the likelihood matrix.")
    }
    log_likelihood<- log_lik_list[[k]]
    L<-loo(log_likelihood, r_eff=r_eff_list[[k]], cores=cores)
    lpd_point[,k] <- L$pointwise[,1]    #calculate log(p_k (y_i | y_-i))
    elpd_loo[k]<-L$estimates[1,1]  #calculate elpd
  }
  if(BB==FALSE){
    k_best_loo<-which(elpd_loo==max(elpd_loo))
    cat("The best model selected by LOO  elpd (without BB) is:\n")
    print(k_best_loo)
    return(k_best_loo)
  }
  else{
    temp<-matrix(NA,BB_n, K)
    if(!is.null(seed))
      set.seed(seed)
    BB_weighting <- dirichlet_rng(BB_n, rep(alpha,N))
    best_count<-rep(0,K)
    for( bb in 1:BB_n){
      z_bb <-   BB_weighting[bb,] %*% lpd_point
      index_best=which( z_bb==max(z_bb) )
      best_count[index_best]=best_count[index_best]+1/length(index_best)
    }
    prob<- best_count/sum(best_count)
    prob_order<-cbind(paste( "Model"  ,c(1:K) ),  prob) [order(prob, decreasing=TRUE),]
    cat("The probability of each model being selected to be best model are :\n")
    print_weight_vector(prob)
    if(max(prob)<=0.5)
      warning("The highest probability of any single model being the best one is smaller than 0.5.
              It is better to do model averaging rather than model selection")
    if(visualise)
    barplot(as.numeric(prob_order[,2]), main="The probability of each single model \n being the best model", names.arg=prob_order[,1], ylim=c(0,1) )
    return(setNames( prob,paste("Model"  ,c(1:K)))  )
  }
}

#' Stacking of predictive distributions
#'
#' @importFrom stats constrOptim
#' @export
#' @param lpd_point A  matrix of pointwise leave-one-out likelihood evaluated in
#'   different models. Each column corresponds to one model.
#' @param optim_method 	The optimization method to be used in stacking; it can
#'   be chosen from "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent".
#'   The default one is "BFGS".
#' @param optim_control A list of control parameters in optimization. See \code{\link{constrOptim}}.

#' @return A vector of best model weights that maximize the leave-one-out log
#'   score of the combination.
#' @details \code{lpd_point} is a matrix of pointwise log leave-one-out
#'   likelihood, which can be calculated from \code{\link{loo}} or through
#'   running exact LOO in advance. It should be a \eqn{N} by \eqn{K}  matrix
#'   when Sample size is \eqn{N} and model number is \eqn{K}. \code{stacking} is
#'   an approach that finds the optimal linear combining weight which maximizes
#'   the leave-one-out log score.
#'
#'@seealso
#' \itemize{
#' \item \code{\link{loo}} for details on leave-one-out elpd estimation.
#' \item \code{\link{model_weights}} for model weighting from stanfit objects.
#' \item \code{\link{pseudobma_weight}} for model weighting by pseudo-BMA and Bayesian bootstrap adjustment.
#' \item \code{\link{optim}} for the choice of optimization methods and other control-parameters.
#' }
#'
#' @examples
#' \dontrun{
#' loo1 <- loo(log_lik1)$pointwise[,1]
#' loo2 <- loo(log_lik2)$pointwise[,1]
#' stacking_weight(cbind(loo1, loo2))}
#'
stacking_weight<-function( lpd_point, optim_method="BFGS", optim_control=list()){
  K<-ncol(lpd_point)                #number of models
  if (K==1)
    stop("Only one model is found.")
  N<-nrow(lpd_point)
  exp_lpd_point<-exp(lpd_point)
  negative_log_score_loo<-function(w)  #objective function: log score
  {
    if(length(w)!=K-1)
      break
    w_full<-c(w, 1-sum(w))
    sum<-0
    for(i in 1:N)
      sum<-sum+log(exp(lpd_point[i,])%*%w_full)
    return(-as.numeric(sum))
  }
  gradient<-function(w)  #gradient of the objective function
  {
    if(length(w)!=K-1)
      break
    w_full<-c(w, 1-sum(w))
    grad<-rep(0,K-1)
    for(k in 1:(K-1))
      for(i in 1:N)
        grad[k]<-grad[k]+ (exp_lpd_point[i,k] -exp_lpd_point[i,K]) /( exp_lpd_point[i,]  %*%w_full)
    return(-grad)
  }
  ui<-rbind(rep(-1,K-1),diag(K-1))  # K-1 simplex constraint matrix
  ci<-c(-1,rep(0,K-1))
  w<-constrOptim (theta =rep(1/K,K-1),f=negative_log_score_loo,
                  grad=gradient, ui=ui, ci=ci, method=optim_method, control=optim_control)$par   # constrOptim: function from {base}.
  return (c(w, 1-sum(w)))
}


#' Pseudo-BMA weighting using Bayesian bootstrap samples adjustment
#'
#' @export
#' @param lpd_point A matrix of pointwise leave-one-out likelihood evaluated in
#'   different models. Each column corresponds to one model.
#' @param BB_n  A positive integer indicating the number of samples in Bayesian
#'   bootstrap. Default is 1000.
#' @param alpha A positive scalar; the shape parameter in the Dirichlet
#'   distribution when doing Bayesian bootstrap. Default is 1.
#' @param seed An integer; optional. If specified, it fixes the random seed when
#'   dong Bayesian bootstrap.
#' @return A vector of model weights.
#' @details \code{lpd_point} is a matrix of pointwise leave-one-out likelihood,
#'   which can be calculated from  \code{\link{loo}} or through running LOO in
#'   advance. It should be a \eqn{N} by \eqn{K} matrix when Sample size is
#'   \eqn{N} and model number is \eqn{K}.  The Bayesian bootstrap  takes into
#'   account the uncertainty of finite data points and regularize the weights
#'   making them go further away from 0 and 1. The shape parameter in the
#'   Dirichlet distribution is \code{alpha}. When \code{alpha}=1, the
#'   distribution is uniform on the simplex space.
#'
#' @seealso
#' \itemize{
#' \item \code{\link{loo}} for details on leave-one-out elpd estimation.
#' \item \code{\link{model_weights}} for model weighting from stanfit objects.
#' \item \code{\link{stacking_weight}} for model weighting by stacking on log score.
#' }
#'
#' @examples
#' \dontrun{
#' loo1 <- loo(log_lik1)$pointwise[,1]
#' loo2 <- loo(log_lik2)$pointwise[,1]
#' pseudobma_weight(cbind(loo1, loo2))}
#'
pseudobma_weight<-function(lpd_point, BB_n=1000, alpha=1, seed=NULL) {
  K<-ncol(lpd_point)                #number of models
  if (K==1)
    stop("Only one model is found.")
  temp<-matrix(NA,BB_n, K)
  N <- nrow(lpd_point)
  if(!is.null(seed))
    set.seed(seed)
  BB_weighting <- dirichlet_rng(BB_n, rep(alpha,N))
  for( bb in 1:BB_n){
    z_bb <-   BB_weighting[bb,] %*% lpd_point * N
    uwts <- exp( z_bb - max(z_bb) )
    temp[bb,] <- uwts / sum(uwts)
  }
  return(colMeans(temp))
}

#' Generate dirichlet simulations, rewritten version
#' @importFrom stats rgamma
#' @noRd
dirichlet_rng <- function(n, alpha) {
  K <- length(alpha)
  gamma_sim<-matrix(rgamma(K * n, alpha), ncol = K, byrow = TRUE)
  return(gamma_sim/rowSums( gamma_sim))
}

#' Print model weights
#' @noRd
print_weight_vector <- function(weights_vector) {
  K=length(weights_vector)
  print(setNames( round(weights_vector,digits=3 ),paste("Model"  ,c(1:K))) ,    quote=FALSE)
}
