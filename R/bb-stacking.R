#' Model averaging via stacking and loo weighting
#' @export
#' @param log_lik_list   A list of pointwise log likelihood simlation matrixs. The \eqn{i} th element corresponds to the \eqn{i} th model. Each row of the matrix is the log likelihood evaluted using a simulated paramter.
#' @param method    One of  \code{"stacking"} or \code{"loo"} indicating which method is to use in obtaining the optimal weights.
#' @param BB    Logicals used in method \code{"loo"}. If \code{Ture}(defalt), Bayesian Bootstrap will be used to adjust the LOO estimator. It helps regularize the weight away from 0 and 1, so as to reduce the variance.
#' @param BB_n    A positivs integer indicating the number of samples in Bayesian Bootstrap. It is necessary when \code{method}=\code{"loo"} and \code{BB}=\code{Ture}. The  default number is 1000.
#' @param alpha A postive scaler; the shape paramter in the Dirichlet distribution when doing Bootstrap.
#' @param  seed An integer; optional. If specified, it will fix the random seed when dong Bayesian Bootstrap.
#' @return A vector of optimal model weights.
#' @details This function implements  \code{stacking} and \code{loo} weighting  for combining multiple predictive distributions. In either method, we can use  Leave-one-out cross-validation (LOO) to estimate the expected log predictive density(elpd).  \code{stacking}  finds the optimal linear weight in combining elpd. \code{loo} finds the relative weights proportional to elpd. However, this estimation doesn’t take into account the uncertainty resulted from having a finite number of proxy samples from the future data distribution. By Bayesian bootstrap, we can take into account the uncertainty and regularize the weights making them go further away from 0 and 1.
#' @seealso
#' \code{\link{Select}} for single-model selection.
#'
#' \code{\link{bb_weight}} for details on Bayesian Bootstrap adjustment.
#' @examples
#' \dontrun{
### Usage with stanfit objects
#' log_lik1 <- extract(stan(model=model_1,data=data))[['log_lik']]
#' log_lik2 <- extract(stan(model=model_1,data=data))[['log_lik']]
#' w1=combine(list(log_lik1, log_lik2),method="stacking")
#' w2=combine(list(log_lik1, log_lik2),method="loo",BB=T)
#' }
#'
Combine <-function(log_lik_list, method, BB=T,BB_n=1000, alpha=1, seed=NULL)
{
  if (!method %in%c("stacking","loo") )
    stop("Must specify a method in stacking or loo.")
  if ( method =="loo" & !is.logical(BB))
    stop("BB must be logical.")
  K<-length(log_lik_list)                #number of models
  if (K==1)
    stop("Only one model is fould.")
  if(length(unique(unlist(lapply(log_lik_list,ncol))))!=1 |length(unique(unlist(lapply(log_lik_list,nrow))))!=1)
    stop("dimensions do not match.")
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
    w_stacking <- stacking_weight(lpd_point)
    cat("The stacking weights are:\n")
    print(rbind(paste("Model"  ,c(1:K) ), round(w_stacking*100 )/100))
    return(w_stacking)
  }
  else
    if (method =="loo"){
      ## 2) loo weights w/o BB sample
      uwts <- exp( elpd_loo - max( elpd_loo))
      w_loo1 <- uwts / sum(uwts)
      if(BB==F){
        cat("The LOO weights are:\n")
        print(rbind(paste("Model"  ,c(1:K) ),  round(w_loo1*100 )/100))
        return(w_loo1)
      }
      if(BB==T){
        w_loo2  <- bb_weight(lpd_point, BB_n,alpha, seed)   #3) loo withs using BB sample
        cat("The LOO weights using BB sample are:\n")
        print(rbind(paste("Model"  ,c(1:K) ),  round(w_loo2*100 )/100))
        return (w_loo2 )
      }
    }
}

#' Model selection via Leave-one-out log predictive density estimation and Bayesian Bootstrap adjustment.
#' @export
#' @param log_lik_list A list of pointwise log likelihood simlation matrixs. The \eqn{i} th element corresponds to the \eqn{i} th model. Each row of the matrix is the log likelihood evaluted using a simulated paramter.
#' @param BB Logicals. If \code{Ture}(defalt), Bayesian Bootstrap will be used to adjust the LOO estimator. It helps regularize the weight away from 0 and 1, so as to reduce the variance.
#' @param BB_n A positivs integer indicating the number of samples in Bayesian Bootstrap. It is necessary  when \code{BB}=\code{Ture}. The  default number is 1000.
#' @param alpha A postive scaler; the shape paramter in the Dirichlet distribution when doing Bootstrap.
#' @param seed An integer; optional. If specified, it fix the random seed when dong Bayesian Bootstrap.
#'@param visualise Logical, whether to visualise the selection probability.
#' @return   When \code{BB}=\code{False}, it returns an integer indicating the  index of the best model.  When\code{BB}=\code{TRUE}, it return a vector indicating the probability of each model being selected to be the best model.
#' @details \code{\link{loo}} gives an estimation of the expected log predictive density of each model, we can pick the best model by picking the model with the largest elpd estimation. Just like \code{\link{Combine}}, to make the elpd estimation more reliable, we can use Bayesian Bootstrap adjustment. With each sample in the Bayesian Bootstrap, we compare the adjusted elpd estimation and finally compute the probability of that model being the optimal one. If none of the probability is close to 1, then it is better to do model averaging rather than model selection.
#' @seealso
#' \code{\link{Combine}} for model combination.
#'
#' \code{\link{compare}} for two-model comparison.
#'
#' \code{\link{bb_weight}} for details on Bayesian Bootstrap adjustment.
#'
#' @examples
#' \dontrun{
#' ### Usage with stanfit objects
#' log_lik1 <- extract(stan(model=model_1,data=data))[['log_lik']]
#' log_lik2 <- extract(stan(model=model_1,data=data))[['log_lik']]
#' k=Select(list(log_lik1,log_lik2),BB=T)
#' }

Select <-function(log_lik_list, BB=T,BB_n=1000, alpha=1,seed=NULL,visualise=T)
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

#' stacking on log score
#' @export
#' @param lpd_point A  matrix of pointwise leave-one-out likelihood evaluated in different models. Each column corresponds to one model.
#' @return A vector of best model weights that minimize the levae-one-out log score.
#' @details \code{lpd_point} is a matrix of pointwise log leave-one-out likelihood, which can be calculeted from  \code{\link{loo}}. It should be a \eqn{N} by \eqn{K}  matrix when Sample size is \eqn{N} and model number is \eqn{K}. \code{stacking} is an approach that finds the optimal linear combining weight so as to minimize the leave-one-out log score.
#'@seealso
#' \code{\link{loo}} for details on leave-one-out elpd estimation.
#'
#' \code{\link{bb-weight}} for model weighting by Leave-one-out and Bayesian bootstrap adjustment.
#'@examples
#' \dontrun{
#' loo1 <- loo(log_lik1)$pointwise[,1]
#' loo2 <- loo(log_lik2)$pointwise[,1]
#' stacking_weight(cbind(loo1, loo2))
#' }
stacking_weight<-function( lpd_point){
  K<-ncol(lpd_point)                #number of models
  if (K==1)
    stop("Only one model is fould.")
  N<-nrow(lpd_point)
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
    exp_lpd_point<-exp(lpd_point)
    for(k in 1:(K-1))
     for(i in 1:N)
       grad[k]<-grad[k]+ (exp_lpd_point[i,k] -exp_lpd_point[i,K]) /( exp_lpd_point[i,]  %*%w_full)
     return(-grad)
  }
  ui<-rbind(rep(-1,K-1),diag(K-1))  # K-1 simplex constraint matrix
  ci<-c(-1,rep(0,K-1))
  w<-constrOptim (theta =rep(1/K,K-1),f=negative_log_score_loo,
                       grad=gradient, ui=ui, ci=ci  )$par   # constrOptim: function from {base}.
  return (c(w, 1-sum(w)))
}


#' LOO weighting using Bayesian bootstrap samples adjustment
#' @export
#' @param lpd_point A matrix of pointwise leave-one-out likelihood evaluated in different models. Each column corresponds to one model.
#' @param  BB_n  A positivs integer indicating the number of samples in Bayesian Bootstrap. Default is 1000.
#' @param alpha A postive scaler; the shape paramter in the Dirichlet distribution when doing Bayesian Bootstrap.
#' @param  seed An integer; optional. If specified, it fix the random seed when dong Bayesian Bootstrap.
#' @return A vector of LOO weights.
#' @details \code{lpd_point} is a matrix of pointwise leave-one-out likelihood, which can be calculeted from  \code{\link{loo}}. It should be a \eqn{N} by \eqn{K}  matrix when Sample size is \eqn{N} and model number is \eqn{K}.   Bayesian bootstrap  takes into account the uncertainty of finite data points and regularize the weights making them go further away from 0 and 1. The shape parameter in the Dirichlet distribution is \code{alpha}. When \code{alpha}=1, the distribution is uniform on the simplex space. A smaller {alpha} will keeps the final weights more away from 0 and 1.
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
#' bb_weight(cbind(loo1, loo2))
#' }
#'
bb_weight <- function(lpd_point, BB_n=1000,alpha=1, seed=NULL) {
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


