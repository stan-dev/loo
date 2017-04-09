#' Model averaging via stacking or pseudo-BMA weighting.
#' @export
#' @param log_lik_list   A list of pointwise log likelihood simulation matrixes. The \eqn{i} th element corresponds to the \eqn{i} th model. Each row of the matrix is the log likelihood vector evaluated using a simulated parameter
#' @param method    One of  \code{"stacking"} or \code{"pseudobma"}, indicating which method is to use for obtaining the optimal weights. \code{"stacking"}  refers to stacking of predictive distributions and  \code{"pseudobma"} refers to pseudo-BMA weighting (by setting \code{"BB"=F}) or pseudo-BMA+ weighting (by setting \code{"BB"=T}).
#' @param BB    Logicals used when \code{"method"}=\code{"pseudobma"}. If \code{True}(default), Bayesian Bootstrap will be used to adjust the pseudo-BMA weighting, which is called pseudo-BMA+ weighting. It helps regularize the weight away from 0 and 1, so as to reduce the variance.
#' @param BB_n    A positive integer indicating the number of samples in Bayesian Bootstrap. It is necessary when  \code{BB}=\code{TRUE}. The  default number is 1000.
#' @param alpha A positive scaler; the shape parameter in the Dirichlet distribution when doing Bootstrap. The default is \eqn{1}.
#' @param  seed An integer; optional. If specified, it will fix the random seed when dong Bayesian Bootstrap sampling.
#' @param  optim_method	The optimization method to be used in stacking. It can be chosen from "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent". The default method is "BFGS".
#'
#' @return A vector of optimal model weights.
#' @details This function implements  stacking of predictive distributions, pseudo-BMA and pseudo-BMA+ weighting for combining multiple predictive distributions.
#'
#' For either method, we can use  Leave-one-out cross-validation (LOO) to estimate the expected log predictive density(elpd).  \code{Stacking} combines all model by maximizing the leave-one-out predictive density of the combination distribution. \code{Pseudo-BMA} finds the relative weights proportional to elpd of each model. \code{Pseudo-BMA+} takes into account the uncertainty resulted from having a finite number of proxy samples from the future data distribution through Bayesian bootstrap (set \code{"BB=T"}), which will keep weights further away from 0 and 1.
#'
#' In general, we recommend  \code{stacking} for averaging predictive distributions, while  \code{Pseudo-BMA+} can serve as a computationally easier alternative.
#' @seealso
#' \code{\link{model_select}} for single-model selection.
#'
#' \code{\link{pseudobma_weight}} for details on pseudo-BMA and pseudo-BMA+ weights.
#'
#' \code{\link{stacking_weight}} for details on stacking weighs.
#'
#' \code{\link{optim}} for the choice of optimization methods.
#'
#' @examples
#' \dontrun{
### Usage with stanfit objects
#' log_lik1 <- extract(stan(model=model_1,data=data))[['log_lik']]
#' log_lik2 <- extract(stan(model=model_1,data=data))[['log_lik']]
#' w1=model_weights(list(log_lik1, log_lik2),method="stacking")
#' w2=model_weights(list(log_lik1, log_lik2),method="pseudobma",BB=T)
#' }
#'
model_weights <-function(log_lik_list, method="stacking",BB=T,BB_n=1000, alpha=1, seed=NULL, optim_method="BFGS")
{
  if (!method %in%c("stacking","pseudobma") )
    stop("Must specify a method in stacking or pseudobma .")
  K<-length(log_lik_list)                #number of models
  if (K==1)
    stop("Only one model is fould.")
  if(length(unique(unlist(lapply(log_lik_list,ncol))))!=1 |length(unique(unlist(lapply(log_lik_list,nrow))))!=1)
    stop("Dimensions do not match. Each element of log_lik_list should have same dimensions.")
  N<-ncol(log_lik_list[[1]])             #number of data points
  lpd_point<-matrix(NA,N,K)            #point wise log likelihood
  elpd_loo<-rep(NA,K)
  for( k in 1:K){
    log_likelihood<- log_lik_list[[k]]
    L<-loo(log_likelihood)
    lpd_point[,k] <- L$pointwise[,1]    #calculate log(p_k (y_i | y_-i))
    elpd_loo[k]<-L$elpd_loo
  }
  ## 1) stacking on log score
  if (method =="stacking"){
    w_stacking <- stacking_weight(lpd_point, optim_method=optim_method)
    cat("The stacking weights are:\n")
    print(rbind(paste("Model"  ,c(1:K) ), round(w_stacking*100 )/100))
    return(w_stacking)
  }
  else
    if (method =="pseudobma"){
      uwts <- exp( elpd_loo - max( elpd_loo))
      w_loo1 <- uwts / sum(uwts)
      if(BB==F){
        cat("The Pseudo-BMA weights are:\n")
        print(rbind(paste("Model"  ,c(1:K) ),  round(w_loo1*100 )/100))
        return(w_loo1)
      }
      if(BB==T){
        w_loo2  <- pseudobma_weight(lpd_point, BB_n,alpha, seed)   #3) loo withs using BB sample
        cat("The Pseudo-BMA+ weights using Bayesian Bootstrap  are:\n ")
        print(rbind(paste("Model",c(1:K) ),  round(w_loo2*100 )/100))
        return (w_loo2 )
      }
    }
}

#' Model selection via Leave-one-out log predictive density estimation and Bayesian Bootstrap adjustment.
#' @export
#' @param log_lik_list A list of pointwise log likelihood simulation matrixes The \eqn{i} th element corresponds to the \eqn{i} th model. Each row of the matrix is the log likelihood vector evaluated using a simulated parameter.
#' @param BB Logicals. If \code{True}(default), Bayesian Bootstrap will be used to adjust the LOO estimator.
#' @param BB_n A positive integer indicating the number of samples in Bayesian Bootstrap. It is necessary  when \code{BB}=\code{True}. The  default number is 1000.
#' @param alpha A positive scaler; the shape parameter in the Dirichlet distribution when doing Bootstrap.
#' @param seed An integer; optional. If specified, it will fix the random seed when dong Bayesian Bootstrap.
#'@param visualise Logical, whether to visualise the selection probability.
#' @return   When \code{BB}=\code{False}, it returns an integer indicating the  index of the best model.  When\code{BB}=\code{TRUE}, it return a vector indicating the probability of each model being selected to be the best model.
#' @details \code{\link{loo}} gives an estimation of the expected log predictive density of each model, we can pick the best model by picking the model with the largest elpd estimation. Just like \code{\link{pseudobma_weight}}, to make the elpd estimation more reliable, we can use Bayesian Bootstrap adjustment. With each sample in the Bayesian Bootstrap, we compare the adjusted elpd estimation and finally compute the probability of that model being the optimal one. If none of the probability is close to 1, then it is better to do model averaging rather than model selection.
#' @seealso
#' \code{\link{Combine}} for model combination.
#'
#' \code{\link{compare}} for two-model comparison.
#'
#' \code{\link{loo_weight}} for details on LOO weighting.
#'

#' @examples
#' \dontrun{
#' ### Usage with stanfit objects
#' library(rstan)
#' log_lik1 <- extract(stan(model=model_1,data=data))[['log_lik']]
#' log_lik2 <- extract(stan(model=model_1,data=data))[['log_lik']]
#' k=model_select(list(log_lik1,log_lik2),BB=T)
#' }

model_select <-function(log_lik_list, BB=T,BB_n=1000, alpha=1,seed=NULL,visualise=T)
{
  if (!is.logical(BB))
    stop("BB must be logical.")
  K<-length(log_lik_list)                #number of models
  if (K==1)
    stop("Only one model is fould.")
  N<-ncol(log_lik_list[[1]])             #number of data points
  lpd_point<-matrix(NA,N,K)            #point wise log likelihood
  elpd_loo<-rep(NA,K)
  for( k in 1:K){
    log_likelihood<- log_lik_list[[k]]
    L<-loo(log_likelihood)
    lpd_point[,k] <- L$pointwise[,1]    #calculate log(p_k (y_i | y_-i))
    elpd_loo[k]<-L$elpd_loo
  }
  if(BB==F){
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
    prob_order<-cbind(paste( "Model"  ,c(1:K) ),  prob) [order(prob, decreasing=T),]
    cat("The probability of each model being selected to be best model are :\n")
    print( prob_order)
    if(max(prob)<=0.5)
      warning("The highest probability of any single model being the best one is smaller than 0.5. It is better to do model averaging rather than model selection")
    if(visualise)
      barplot(as.numeric(prob_order[,2]), main="The probability of each single model being the best", names.arg=prob_order[,1], ylim=c(0,1) )
    return(prob_order)
  }
}

#' Stacking of predictive distributions
#' @export
#' @param lpd_point A  matrix of pointwise leave-one-out likelihood evaluated in different models. Each column corresponds to one model.
#' @param  optim_method	The optimization method to be used in stacking; it can be chosen from "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent". The default one is "BFGS".
#'
#' @return A vector of best model weights that maximize the leave-one-out log score of the combination.
#' @details \code{lpd_point} is a matrix of pointwise log leave-one-out likelihood, which can be calculated from \code{\link{loo}} or through running exact LOO in advance. It should be a \eqn{N} by \eqn{K}  matrix when Sample size is \eqn{N} and model number is \eqn{K}. \code{stacking} is an approach that finds the optimal linear combining weight which maximizes the leave-one-out log score.
#'
#'@seealso
#' \code{\link{loo}} for details on leave-one-out elpd estimation.
#'
#' \code{\link{pseudobma_weight}} for model weighting by pseudo-BMA and Bayesian bootstrap adjustment.
#'
#' \code{\link{optim}} for the choice of optimization methods.
#'

#'@examples
#' \dontrun{
#' loo1 <- loo(log_lik1)$pointwise[,1]
#' loo2 <- loo(log_lik2)$pointwise[,1]
#' stacking_weight(cbind(loo1, loo2))
#' }
stacking_weight<-function( lpd_point, optim_method="BFGS"){
  K<-ncol(lpd_point)                #number of models
  if (K==1)
    stop("Only one model is fould.")
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
                  grad=gradient, ui=ui, ci=ci, method=optim_method )$par   # constrOptim: function from {base}.
  return (c(w, 1-sum(w)))
}


#' Pseudo-BMA weighting using Bayesian bootstrap samples adjustment
#' @export
#' @param lpd_point A matrix of pointwise leave-one-out likelihood evaluated in different models. Each column corresponds to one model.
#' @param  BB_n  A positive integer indicating the number of samples in Bayesian Bootstrap. Default is 1000.
#' @param alpha A postive scaler; the shape parameter in the Dirichlet distribution when doing Bayesian Bootstrap. Default is 1.
#' @param  seed An integer; optional. If specified, it fixes the random seed when dong Bayesian Bootstrap.
#' @return A vector of model weights.
#' @details \code{lpd_point} is a matrix of pointwise leave-one-out likelihood, which can be calculated from  \code{\link{loo}} or through running LOO in advance. It should be a \eqn{N} by \eqn{K} matrix when Sample size is \eqn{N} and model number is \eqn{K}.   Bayesian bootstrap  takes into account the uncertainty of finite data points and regularize the weights making them go further away from 0 and 1. The shape parameter in the Dirichlet distribution is \code{alpha}. When \code{alpha}=1, the distribution is uniform on the simplex space. A smaller \code{alpha} will keeps the final weights more away from 0 and 1.
#'
#' @seealso
#' \code{\link{loo}} for details on leave-one-out elpd estimation.
#'
#' \code{\link{stacking_weight}} for model weighting by stacking on log score.
#'
#' @examples
#' \dontrun{
#' loo1 <- loo(log_lik1)$pointwise[,1]
#' loo2 <- loo(log_lik2)$pointwise[,1]
#' pseudobma_weight(cbind(loo1, loo2))
#' }
#'
pseudobma_weight<-function(lpd_point, BB_n=1000,alpha=1, seed=NULL) {
  K<-ncol(lpd_point)                #number of models
  if (K==1)
    stop("Only one model is fould.")
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

#  generate dirichlet simulations, rewritten version
dirichlet_rng <- function(n, alpha) {
  K <- length(alpha)
  gamma_sim<-matrix(rgamma(K * n, alpha), ncol = K, byrow = TRUE)
  return(gamma_sim/rowSums( gamma_sim))
}


