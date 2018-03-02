#' Model averaging/weighting via stacking or pseudo-BMA weighting
#'
#' @description Model averaging via stacking of predictive distributions,
#'   pseudo-BMA weighting or pseudo-BMA+ weighting with the Bayesian bootstrap.
#'   See Yao et al. (2017) and  Vehtari, Gelman, and Gabry (2017) for
#'   background.
#'
#' @export
#' @param x A list of pointwise log-likelihood matrices, one for each model.
#'   Each matrix should have dimensions \eqn{S} by \eqn{N}, where \eqn{S} is the
#'   size of the posterior sample (with all chains merged) and \eqn{N} is the
#'   number of data points.
#' @param method Either \code{"stacking"} or \code{"pseudobma"}, indicating
#'   which method to use for obtaining the weights. \code{"stacking"} refers to
#'   stacking of predictive distributions and  \code{"pseudobma"} refers to
#'   pseudo-BMA weighting (by setting \code{BB=FALSE}) or pseudo-BMA+ weighting
#'   (by leaving the default \code{BB=TRUE}).
#' @param BB Logical used when \code{"method"}=\code{"pseudobma"}. If
#'   \code{TRUE} (the default), the Bayesian bootstrap will be used to adjust the
#'   pseudo-BMA weighting, which is called pseudo-BMA+ weighting. It helps
#'   regularize the weight away from 0 and 1, so as to reduce the variance.
#' @param BB_n When \code{BB}=\code{TRUE}, a positive integer indicating the
#'   number of samples for the Bayesian bootstrap. The default is 1000.
#' @param alpha A positive scalar; the shape parameter in the Dirichlet
#'   distribution used for the Bayesian bootstrap. The default is 1, which
#'   corresponds to a uniform distribution on the simplex space.
#' @param seed When \code{BB}=\code{TRUE}, an optional integer seed for the
#'   Bayesian bootstrap sampling.
#' @param optim_method The optimization method to use if
#'   \code{method="stacking"}. It can be chosen from "Nelder-Mead", "BFGS",
#'   "CG", "L-BFGS-B", "SANN" and "Brent". The default method is "BFGS".
#' @param optim_control If \code{method="stacking"}, a list of control
#'   parameters for optimization. See \code{\link{constrOptim}} for details.
#' @param r_eff_list Optionally, a list of relative effective sample size
#'   estimates for the likelihood \code{(exp(log_lik))} of each observation in
#'   each model. See \code{\link{psis}} and  \code{\link{relative_eff}} helper
#'   function for computing \code{r_eff}.
#' @param cores	 The number of cores to use for parallelization. Same as for
#'   \code{\link{psis}}. The default for an entire R session can be set with
#'   \code{options(loo.cores = NUMBER)}. \strong{As of version 2.0.0 the default
#'   is now 1 core, but we recommend using as many (or close to as many) cores
#'   as possible.}
#' @param ... Unused, except for the generic to pass arguments to individual
#'   methods.
#'
#' @return A numeric vector containing one weight for each model.
#'
#' @details
#' \code{model_weights} implements stacking, pseudo-BMA, and pseudo-BMA+
#' weighting for combining multiple predictive distributions. In all cases, we
#' can use leave-one-out cross-validation (LOO) to estimate the expected log
#' predictive density (ELPD).
#'
#' The stacking method (\code{method="stacking"}) combines all models by
#' maximizing the leave-one-out predictive density of the combination
#' distribution. That is, it finds the optimal linear combining weights for
#' maximizing the leave-one-out log score.
#'
#' The pseudo-BMA method (\code{method="pseudobma"}) finds the relative weights
#' proportional to the ELPD of each model. However, when
#' \code{method="pseudobma"}, the default is to also use the Bayesian bootstrap
#' (\code{BB=TRUE}), which corresponds to the pseudo-BMA+ method. The Bayesian
#' bootstrap  takes into account the uncertainty of finite data points and
#' regularizes the weights away from the extremes of 0 and 1.
#'
#' In general, we recommend stacking for averaging predictive distributions,
#' while pseudo-BMA+ can serve as a computationally easier alternative.
#'
#' @seealso
#' \itemize{
#' \item \code{\link{loo}} for details on leave-one-out ELPD estimation.
#' \item \code{\link{constrOptim}} for the choice of optimization methods and control-parameters.
#' \item \code{\link{relative_eff}} for computing \code{r_eff}.
#' }
#'
#'@template stacking-references
#'
#' @examples
#' \dontrun{
#' ### Demonstrating usage after fitting models with RStan
#' library(rstan)
#'
#' # generate fake data from N(0,1).
#' set.seed(100)
#' N <- 100
#' y <- rnorm(N, 0, 1)
#'
#' # Suppose we have three models: N(-1, sigma), N(0.5, sigma) and N(0.6,sigma).
#' stan_code <- "
#'   data {
#'     int N;
#'     vector[N] y;
#'     real mu_fixed;
#'   }
#'   parameters {
#'     real<lower=0> sigma;
#'   }
#'   model {
#'     sigma ~ exponential(1);
#'     y ~ normal(mu_fixed, sigma);
#'   }
#'   generated quantities {
#'     vector[N] log_lik;
#'     for (n in 1:N) log_lik[n] = normal_lpdf(y[n]| mu_fixed, sigma);
#'   }"
#'
#' mod <- stan_model(model_code = stan_code)
#' fit1 <- sampling(mod, data=list(N=N, y=y, mu_fixed=-1))
#' fit2 <- sampling(mod, data=list(N=N, y=y, mu_fixed=0.5))
#' fit3 <- sampling(mod, data=list(N=N, y=y, mu_fixed=0.6))
#' log_lik_list <- lapply(c(fit1, fit2, fit3), extract_log_lik)
#'
#' # optional but recommended
#' r_eff_list <- lapply(c(fit1, fit2, fit3), function(x) {
#'   relative_eff(exp(extract_log_lik(x, merge_chains = FALSE)))
#' })
#'
#' # stacking method:
#' model_weights(
#'   log_lik_list,
#'   method="stacking",
#'   r_eff_list = r_eff_list,
#'   optim_control = list(reltol=1e-10)
#' )
#'
#' # pseudo-BMA+ method:
#' model_weights(
#'   log_lik_list,
#'   method = "pseudobma",
#'   r_eff_list=r_eff_list
#'  )
#'
#' # pseudo-BMA method (set BB = FALSE):
#' model_weights(
#'   log_lik_list,
#'   method = "pseudobma",
#'   BB = FALSE,
#'   r_eff_list=r_eff_list
#'  )
#'
#' # calling stacking_weights or pseudobma_weights directly
#' lpd1 <- loo(log_lik_list[[1]], r_eff = r_eff_list[[1]])$pointwise[,1]
#' lpd2 <- loo(log_lik_list[[2]], r_eff = r_eff_list[[2]])$pointwise[,1]
#' lpd3 <- loo(log_lik_list[[3]], r_eff = r_eff_list[[3]])$pointwise[,1]
#' stacking_weights(cbind(lpd1, lpd2, lpd3))
#' pseudobma_weights(cbind(lpd1, lpd2, lpd3))
#' pseudobma_weights(cbind(lpd1, lpd2, lpd3), BB = FALSE)
#' }
#'
model_weights <- function(x, ...) {
  UseMethod("model_weights")
}

#' @rdname model_weights
#' @export
model_weights.default <-
  function(x,
           ...,
           method = c("stacking", "pseudobma"),
           optim_method = "BFGS",
           optim_control = list(),
           BB = TRUE,
           BB_n = 1000,
           alpha = 1,
           seed = NULL,
           r_eff_list = NULL,
           cores = getOption("loo.cores", 1)) {

    method <- match.arg(method)
    K <- length(x) # number of models
    N <- ncol(x[[1]]) # number of data points
    validate_log_lik_list(x)
    validate_r_eff_list(r_eff_list, K, N)

    lpd_point <- matrix(NA, N, K)
    elpd_loo <- rep(NA, K)
    for (k in 1:K) {
      r_eff_k <- r_eff_list[[k]] # possibly NULL
      log_likelihood <- x[[k]]
      L <- loo(log_likelihood, r_eff = r_eff_k, cores = cores)
      lpd_point[, k] <- L$pointwise[, 1]    #calculate log(p_k (y_i | y_-i))
      elpd_loo[k] <- L$estimates["elpd_loo", 1]
    }

    ## 1) stacking on log score
    if (method =="stacking"){
      wts <- stacking_weights(
        lpd_point = lpd_point,
        optim_method = optim_method,
        optim_control = optim_control
      )

    } else {
      # method =="pseudobma"
      wts <- pseudobma_weights(
        lpd_point = lpd_point,
        BB = BB,
        BB_n = BB_n,
        alpha = alpha,
        seed = seed
      )
    }

    return(wts)
  }


#' @rdname model_weights
#' @export
#' @param lpd_point A matrix of pointwise log leave-one-out likelihoods
#'   evaluated for different models. It should be a \eqn{N} by \eqn{K}  matrix
#'   where \eqn{N} is sample size and \eqn{K} is the number of models. Each
#'   column corresponds to one model. These values can be calculated
#'   approximately using \code{\link{loo}} or by running exact leave-one-out
#'   cross-validation.
#'
#' @importFrom stats constrOptim
#'
stacking_weights <-
  function(lpd_point,
           optim_method = "BFGS",
           optim_control = list()) {

    stopifnot(is.matrix(lpd_point))
    N <- nrow(lpd_point)
    K <- ncol(lpd_point)
    if (K < 2) {
      stop("At least two models are required for stacking weights.")
    }

    exp_lpd_point <- exp(lpd_point)
    negative_log_score_loo <- function(w) {
      # objective function: log score
      stopifnot(length(w) == K - 1)
      w_full <- c(w, 1 - sum(w))
      sum <- 0
      for (i in 1:N) {
        sum <- sum + log(exp(lpd_point[i, ]) %*% w_full)
      }
      return(-as.numeric(sum))
    }

    gradient <- function(w) {
      # gradient of the objective function
      stopifnot(length(w) == K - 1)
      w_full <- c(w, 1 - sum(w))
      grad <- rep(0, K - 1)
      for (k in 1:(K - 1)) {
        for (i in 1:N) {
          grad[k] <- grad[k] +
            (exp_lpd_point[i, k] - exp_lpd_point[i, K]) / (exp_lpd_point[i,]  %*% w_full)
        }
      }
      return(-grad)
    }

    ui <- rbind(rep(-1, K - 1), diag(K - 1))  # K-1 simplex constraint matrix
    ci <- c(-1, rep(0, K - 1))
    w <- constrOptim(
      theta = rep(1 / K, K - 1),
      f = negative_log_score_loo,
      grad = gradient,
      ui = ui,
      ci = ci,
      method = optim_method,
      control = optim_control
    )$par

    wts <- structure(
      c(w, 1 - sum(w)),
      names = paste0("model", 1:K),
      class = c("stacking_weights")
    )

    return(wts)
  }


#' @rdname model_weights
#' @export
#'
pseudobma_weights <-
  function(lpd_point,
           BB = TRUE,
           BB_n = 1000,
           alpha = 1,
           seed = NULL) {
    stopifnot(is.matrix(lpd_point))
    N <- nrow(lpd_point)
    K <- ncol(lpd_point)
    if (K < 2) {
      stop("At least two models are required for pseudo-BMA weights.")
    }
    if (!is.null(seed)) {
      set.seed(seed)
    }

    if (!BB) {
      elpd <- colSums2(lpd_point)
      uwts <- exp(elpd - max(elpd))
      wts <- structure(
        uwts / sum(uwts),
        names = paste0("model", 1:K),
        class = "pseudobma_weights"
      )
      return(wts)
    }

    temp <- matrix(NA, BB_n, K)
    BB_weighting <- dirichlet_rng(BB_n, rep(alpha, N))
    for (bb in 1:BB_n) {
      z_bb <- BB_weighting[bb, ] %*% lpd_point * N
      uwts <- exp(z_bb - max(z_bb))
      temp[bb, ] <- uwts / sum(uwts)
    }
    wts <- structure(
      colMeans(temp),
      names = paste0("model", 1:K),
      class = "pseudobma_bb_weights"
    )
    return(wts)
  }


#' Generate dirichlet simulations, rewritten version
#' @importFrom stats rgamma
#' @noRd
dirichlet_rng <- function(n, alpha) {
  K <- length(alpha)
  gamma_sim <- matrix(rgamma(K * n, alpha), ncol = K, byrow = TRUE)
  return(gamma_sim / rowSums(gamma_sim))
}

#' @export
print.stacking_weights <- function(x, digits = 3, ...) {
  cat("Method: stacking\n------\n")
  print_weight_vector(x, digits = digits)
}

#' @export
print.pseudobma_weights <- function(x, digits = 3, ...) {
  cat("Method: pseudo-BMA\n------\n")
  print_weight_vector(x, digits = digits)
}

#' @export
print.pseudobma_bb_weights <- function(x, digits = 3, ...) {
  cat("Method: pseudo-BMA+ with Bayesian bootstrap\n------\n")
  print_weight_vector(x, digits = digits)
}

print_weight_vector <- function(x, digits) {
  z <- cbind(x)
  colnames(z) <- "weight"
  print(.fr(z, digits = digits), quote = FALSE)
  invisible(x)
}

#' Validate r_eff_list argument if provided
#'
#' @noRd
#' @param r_eff_list User's r_eff_list argument
#' @param K Required length of r_eff_list (number of models).
#' @param N Required length of each element of r_eff_list (number of data points).
#' @return Either throws an error or returns \code{TRUE} invisibly.
#'
validate_r_eff_list <- function(r_eff_list, K, N) {
  if (is.null(r_eff_list)) return(invisible(TRUE))

  if (length(r_eff_list) != K) {
    stop("If r_eff_list is specified then it must contain ",
         "one component for each model being compared.",
         call. = FALSE)
  }
  if (any(sapply(r_eff_list, length) != N)) {
    stop("Each component of r_eff list must have the same length ",
         "as the number of columns in the log-likelihood matrix.",
         call. = FALSE)
  }
  invisible(TRUE)
}


#' Validate log-likelihood list argument
#'
#' Checks that log-likelihood list has at least 2 elements and that each element
#' has the same dimensions.
#'
#' @noRd
#' @param log_lik_list User's list of log-likelihood matrices (the 'x' argument
#'   to model_weights).
#' @return Either throws an error or returns \code{TRUE} invisibly.
#'
validate_log_lik_list <- function(log_lik_list) {
  stopifnot(is.list(log_lik_list))
  if (length(log_lik_list) < 2) {
    stop("At least two models are required.", call. = FALSE)
  }
  if (length(unique(sapply(log_lik_list, ncol))) != 1 |
     length(unique(sapply(log_lik_list, nrow))) != 1) {
    stop("Each log-likelihood matrix must have the same dimensions.", call. = FALSE)
  }
  return(invisible(TRUE))
}

