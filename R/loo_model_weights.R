#' Model averaging/weighting via stacking or pseudo-BMA weighting
#'
#' Model averaging via stacking of predictive distributions, pseudo-BMA
#' weighting or pseudo-BMA+ weighting with the Bayesian bootstrap. See Yao et
#' al. (2018), Vehtari, Gelman, and Gabry (2017), and Vehtari, Simpson,
#' Gelman, Yao, and Gabry (2024) for background.
#'
#' @export
#' @param x A list of `"psis_loo"` objects (objects returned by [loo()]) or
#'   pointwise log-likelihood matrices or , one for each model. If the list
#'   elements are named the names will be used to label the models in the
#'   results. Each matrix/object should have dimensions \eqn{S} by \eqn{N},
#'   where \eqn{S} is the size of the posterior sample (with all chains merged)
#'   and \eqn{N} is the number of data points. If `x` is a list of
#'   log-likelihood matrices then [loo()] is called internally on each matrix.
#'   Currently the `loo_model_weights()` function is not implemented to be used
#'   with results from K-fold CV, but you can still obtain weights using K-fold
#'   CV results by calling the `stacking_weights()` or `pseudobma_weights()`
#'   function directly.
#' @param method Either `"stacking"` (the default) or `"pseudobma"`, indicating which method
#'   to use for obtaining the weights. `"stacking"` refers to stacking of
#'   predictive distributions and  `"pseudobma"` refers to pseudo-BMA+ weighting
#'   (or plain pseudo-BMA weighting if argument `BB` is `FALSE`).
#' @param BB Logical used when `"method"`=`"pseudobma"`. If
#'   `TRUE` (the default), the Bayesian bootstrap will be used to adjust
#'   the pseudo-BMA weighting, which is called pseudo-BMA+ weighting. It helps
#'   regularize the weight away from 0 and 1, so as to reduce the variance.
#' @param BB_n For pseudo-BMA+ weighting only, the number of samples to use for
#'   the Bayesian bootstrap. The default is `BB_n=1000`.
#' @param alpha Positive scalar shape parameter in the Dirichlet distribution
#'   used for the Bayesian bootstrap. The default is `alpha=1`, which
#'   corresponds to a uniform distribution on the simplex space.
#' @param optim_method If `method="stacking"`, a string passed to the `method`
#'   argument of [stats::constrOptim()] to specify the optimization algorithm.
#'   The default is `optim_method="BFGS"`, but other options are available (see
#'   [stats::optim()]).
#' @param optim_control If `method="stacking"`, a list of control parameters for
#'   optimization passed to the `control` argument of [stats::constrOptim()].
#' @param r_eff_list Optionally, a list of relative effective sample size
#'   estimates for the likelihood `(exp(log_lik))` of each observation in
#'   each model. See [psis()] and  [relative_eff()] helper
#'   function for computing `r_eff`. If `x` is a list of `"psis_loo"`
#'   objects then `r_eff_list` is ignored.
#' @template cores
#' @param ... Unused, except for the generic to pass arguments to individual
#'   methods.
#'
#' @return A numeric vector containing one weight for each model.
#'
#' @details
#' `loo_model_weights()` is a wrapper around the `stacking_weights()` and
#' `pseudobma_weights()` functions that implements stacking, pseudo-BMA, and
#' pseudo-BMA+ weighting for combining multiple predictive distributions. We can
#' use approximate or exact leave-one-out cross-validation (LOO-CV) or K-fold CV
#' to estimate the expected log predictive density (ELPD).
#'
#' The stacking method (`method="stacking"`), which is the default for
#' `loo_model_weights()`, combines all models by maximizing the leave-one-out
#' predictive density of the combination distribution. That is, it finds the
#' optimal linear combining weights for maximizing the leave-one-out log score.
#'
#' The pseudo-BMA method (`method="pseudobma"`) finds the relative weights
#' proportional to the ELPD of each model. However, when
#' `method="pseudobma"`, the default is to also use the Bayesian bootstrap
#' (`BB=TRUE`), which corresponds to the pseudo-BMA+ method. The Bayesian
#' bootstrap  takes into account the uncertainty of finite data points and
#' regularizes the weights away from the extremes of 0 and 1.
#'
#' In general, we recommend stacking for averaging predictive distributions,
#' while pseudo-BMA+ can serve as a computationally easier alternative.
#'
#' @seealso
#' * The __loo__ package [vignettes](https://mc-stan.org/loo/articles/), particularly
#'   [Bayesian Stacking and Pseudo-BMA weights using the __loo__ package](https://mc-stan.org/loo/articles/loo2-weights.html).
#' * [loo()] for details on leave-one-out ELPD estimation.
#' * [constrOptim()] for the choice of optimization methods and control-parameters.
#' * [relative_eff()] for computing `r_eff`.
#'
#' @template loo-and-psis-references
#' @template stacking-references
#'
#' @examples
#' \dontrun{
#' ### Demonstrating usage after fitting models with RStan
#' library(rstan)
#'
#' # generate fake data from N(0,1).
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
#' model_list <- list(fit1, fit2, fit3)
#' log_lik_list <- lapply(model_list, extract_log_lik)
#'
#' # optional but recommended
#' r_eff_list <- lapply(model_list, function(x) {
#'   ll_array <- extract_log_lik(x, merge_chains = FALSE)
#'   relative_eff(exp(ll_array))
#' })
#'
#' # stacking method:
#' wts1 <- loo_model_weights(
#'   log_lik_list,
#'   method = "stacking",
#'   r_eff_list = r_eff_list,
#'   optim_control = list(reltol=1e-10)
#' )
#' print(wts1)
#'
#' # can also pass a list of psis_loo objects to avoid recomputing loo
#' loo_list <- lapply(1:length(log_lik_list), function(j) {
#'   loo(log_lik_list[[j]], r_eff = r_eff_list[[j]])
#' })
#'
#' wts2 <- loo_model_weights(
#'   loo_list,
#'   method = "stacking",
#'   optim_control = list(reltol=1e-10)
#' )
#' all.equal(wts1, wts2)
#'
#' # can provide names to be used in the results
#' loo_model_weights(setNames(loo_list, c("A", "B", "C")))
#'
#'
#' # pseudo-BMA+ method:
#' set.seed(1414)
#' loo_model_weights(loo_list, method = "pseudobma")
#'
#' # pseudo-BMA method (set BB = FALSE):
#' loo_model_weights(loo_list, method = "pseudobma", BB = FALSE)
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
loo_model_weights <- function(x, ...) {
  UseMethod("loo_model_weights")
}

#' @rdname loo_model_weights
#' @export
#' @export loo_model_weights.default
loo_model_weights.default <-
  function(x,
           ...,
           method = c("stacking", "pseudobma"),
           optim_method = "BFGS",
           optim_control = list(),
           BB = TRUE,
           BB_n = 1000,
           alpha = 1,
           r_eff_list = NULL,
           cores = getOption("mc.cores", 1)) {

    cores <- loo_cores(cores)
    method <- match.arg(method)
    K <- length(x) # number of models

    if (is.matrix(x[[1]])) {
      N <- ncol(x[[1]]) # number of data points
      validate_log_lik_list(x)
      validate_r_eff_list(r_eff_list, K, N)
      lpd_point <- matrix(NA, N, K)
      elpd_loo <- rep(NA, K)
      for (k in 1:K) {
        r_eff_k <- r_eff_list[[k]] # possibly NULL
        log_likelihood <- x[[k]]
        loo_object <- loo(log_likelihood, r_eff = r_eff_k, cores = cores)
        lpd_point[, k] <- loo_object$pointwise[, "elpd_loo"]    #calculate log(p_k (y_i | y_-i))
        elpd_loo[k] <- loo_object$estimates["elpd_loo", "Estimate"]
      }
    } else if (is.psis_loo(x[[1]])) {
      validate_psis_loo_list(x)
      lpd_point <- do.call(cbind, lapply(x, function(obj) obj$pointwise[, "elpd_loo"]))
      elpd_loo <- sapply(x, function(obj) obj$estimates["elpd_loo", "Estimate"])
    } else {
      stop("'x' must be a list of matrices or a list of 'psis_loo' objects.")
    }

    ## 1) stacking on log score
    if (method =="stacking") {
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
        alpha = alpha
      )
    }

    if (is.matrix(x[[1]])) {
      if (!is.null(names(x)) && all(nzchar(names(x)))) {
        wts <- setNames(wts, names(x))
      }
    } else { # list of loo objects
      wts <- setNames(wts, find_model_names(x))
    }
    wts
  }


#' @rdname loo_model_weights
#' @export
#' @param lpd_point If calling `stacking_weights()` or `pseudobma_weights()`
#'   directly, a matrix of pointwise leave-one-out (or K-fold) log likelihoods
#'   evaluated for different models. It should be a \eqn{N} by \eqn{K}  matrix
#'   where \eqn{N} is sample size and \eqn{K} is the number of models. Each
#'   column corresponds to one model. These values can be calculated
#'   approximately using [loo()] or by running exact leave-one-out or K-fold
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

    negative_log_score_loo <- function(w) {
      # objective function: log score
      stopifnot(length(w) == K - 1)
      w_full <- c(w, 1 - sum(w))
      # avoid over- and underflows using log weights and rowLogSumExps
      sum <- sum(matrixStats::rowLogSumExps(sweep(lpd_point[1:N,], 2, log(w_full), '+')))
      return(-as.numeric(sum))
    }

    gradient <- function(w) {
      # gradient of the objective function
      stopifnot(length(w) == K - 1)
      w_full <- c(w, 1 - sum(w))
      grad <- rep(0, K - 1)
      # avoid over- and underflows using log weights, rowLogSumExps,
      # and by subtracting the row maximum of lpd_point
      mlpd <- matrixStats::rowMaxs(lpd_point)
      for (k in 1:(K - 1)) {
        grad[k] <- sum((exp(lpd_point[, k] - mlpd) - exp(lpd_point[, K] - mlpd)) / exp(matrixStats::rowLogSumExps(sweep(lpd_point, 2, log(w_full), '+')) - mlpd))
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


#' @rdname loo_model_weights
#' @export
#'
pseudobma_weights <-
  function(lpd_point,
           BB = TRUE,
           BB_n = 1000,
           alpha = 1) {
    stopifnot(is.matrix(lpd_point))
    N <- nrow(lpd_point)
    K <- ncol(lpd_point)
    if (K < 2) {
      stop("At least two models are required for pseudo-BMA weights.")
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
  gamma_sim / rowSums(gamma_sim)
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
#' @param r_eff_list User's `r_eff_list` argument
#' @param K Required length of `r_eff_list` (number of models).
#' @param N Required length of each element of `r_eff_list` (number of data points).
#' @return Either throws an error or returns `TRUE` invisibly.
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
#' @param log_lik_list User's list of log-likelihood matrices (the `x` argument
#'   to loo_model_weights).
#' @return Either throws an error or returns `TRUE` invisibly.
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
  invisible(TRUE)
}

validate_psis_loo_list <- function(psis_loo_list) {
  stopifnot(is.list(psis_loo_list))
  if (length(psis_loo_list) < 2) {
    stop("At least two models are required.", call. = FALSE)
  }
  if (!all(sapply(psis_loo_list, is.psis_loo))) {
    stop("List elements must all be 'psis_loo' objects or log-likelihood matrices.")
  }

  dims <- sapply(psis_loo_list, dim)
  if (length(unique(dims[1, ])) != 1 |
      length(unique(dims[2, ])) != 1) {
    stop("Each object in the list must have the same dimensions.", call. = FALSE)
  }
  invisible(TRUE)
}
