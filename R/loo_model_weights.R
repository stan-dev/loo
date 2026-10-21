#' Model averaging/weighting via stacking or pseudo-BMA weighting
#'
#' Model averaging via stacking of predictive distributions, pseudo-BMA
#' weighting or pseudo-BMA+ weighting with the Bayesian bootstrap. See Yao et
#' al. (2018), Vehtari, Gelman, and Gabry (2017), and Vehtari, Simpson,
#' Gelman, Yao, and Gabry (2024) for background.
#'
#' @export
#' @param x A list of `"psis_loo"` objects (objects returned by [loo()]) or
#'   pointwise log-likelihood matrices, one for each model. If the list
#'   elements are named the names will be used to label the models in the
#'   results. Each matrix/object should have dimensions \eqn{S} by \eqn{N},
#'   where \eqn{S} is the size of the posterior sample (with all chains merged)
#'   and \eqn{N} is the number of data points. All models must have the same
#'   \eqn{N}, but \eqn{S} can differ between models. If `x` is a list of
#'   log-likelihood matrices then [loo()] is called internally on each matrix.
#'   Currently the `loo_model_weights()` function is not implemented to be used
#'   with results from K-fold CV, but you can still obtain weights using K-fold
#'   CV results by calling the `stacking_weights()` or `pseudobma_weights()`
#'   function directly. Objects from [loo_subsample()] are not supported,
#'   because they have pointwise values only for the subsampled observations.
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
#' @param optim_method,optim_control Deprecated and ignored. Earlier versions
#'   passed them to [stats::constrOptim()]. The stacking weights are now
#'   computed with an active-set Newton method, see Details.
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
#' `stacking_weights()` maximizes the log score with an active-set method. It
#' starts with the model with the largest total log predictive density, finds
#' the optimal weights of the models in the active set with Newton steps, and
#' adds the models with the largest gradients of the log score (at most ten at a
#' time), until no model can improve it. A model whose weight reaches zero
#' leaves the active set. At the returned weights, the log score is at most
#' \eqn{N \times 10^{-8}} below its maximum. If the method cannot reach this
#' bound, it gives a warning. The optimal weights are usually sparse, so the
#' method is fast also for many models. The method is a variant of the
#' constrained Newton method for mixture proportions (Wang, 2007), and the
#' stopping rule uses the optimality condition of Lindsay (1983).
#'
#' Wang, Y. (2007). On fast computation of the non-parametric maximum
#' likelihood estimate of a mixing distribution. *Journal of the Royal
#' Statistical Society: Series B*, 69(2), 185--198.
#'
#' Lindsay, B. G. (1983). The geometry of mixture likelihoods: a general
#' theory. *The Annals of Statistics*, 11(1), 86--94.
#'
#' @seealso
#' * The __loo__ package [vignettes](https://mc-stan.org/loo/articles/), particularly
#'   [Bayesian Stacking and Pseudo-BMA weights using the __loo__ package](https://mc-stan.org/loo/articles/loo2-weights.html).
#' * [loo()] for details on leave-one-out ELPD estimation.
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
#'   r_eff_list = r_eff_list
#' )
#' print(wts1)
#'
#' # can also pass a list of psis_loo objects to avoid recomputing loo
#' loo_list <- lapply(1:length(log_lik_list), function(j) {
#'   loo(log_lik_list[[j]], r_eff = r_eff_list[[j]])
#' })
#'
#' wts2 <- loo_model_weights(loo_list, method = "stacking")
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
           optim_method = NULL,
           optim_control = NULL,
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
stacking_weights <-
  function(lpd_point,
           optim_method = NULL,
           optim_control = NULL) {

    stopifnot(is.matrix(lpd_point))
    N <- nrow(lpd_point)
    K <- ncol(lpd_point)
    if (K < 2) {
      stop("At least two models are required for stacking weights.")
    }
    validate_lpd_point(lpd_point)
    if (any(rowSums(is.finite(lpd_point)) == 0)) {
      stop("Each observation must have a finite predictive density for at least one model.")
    }

    if (!is.null(optim_method) || !is.null(optim_control)) {
      warning(
        "The arguments 'optim_method' and 'optim_control' are deprecated and ",
        "ignored. The stacking weights are now computed with an active-set ",
        "Newton method.",
        call. = FALSE
      )
    }

    structure(
      stacking_weights_active_set(lpd_point),
      names = paste0("model", 1:K),
      class = c("stacking_weights")
    )
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
    validate_lpd_point(lpd_point)
    elpd <- colSums2(lpd_point)
    if (!any(is.finite(elpd))) {
      stop("At least one model must have a finite total predictive density.")
    }

    if (!BB) {
      uwts <- exp(elpd - max(elpd))
      wts <- structure(
        uwts / sum(uwts),
        names = paste0("model", 1:K),
        class = "pseudobma_weights"
      )
      return(wts)
    }

    BB_weighting <- dirichlet_rng(BB_n, rep(alpha, N))
    z <- BB_weighting %*% lpd_point * N
    uwts <- exp(z - rowMaxs(z))
    structure(
      colMeans(uwts / rowSums(uwts)),
      names = paste0("model", 1:K),
      class = "pseudobma_bb_weights"
    )
  }


#' Stacking weights with an active-set Newton method
#'
#' Maximizes `sum_i log(sum_k w_k exp(lpd_point[i, k]))` over the simplex. The
#' problem is concave in `w`, and the optimal weights are usually sparse. The
#' method optimizes the weights of the active models with Newton steps and adds
#' the models with the largest `g_k = mean_i exp(lpd_point[i, k]) / p_i`, where
#' `p_i` is the predictive density of the combination, until `max_k g_k` is at
#' most `exp(tol)`. Then the log score is at most `N * tol` below its maximum,
#' by Jensen's inequality.
#'
#' @noRd
#' @param lpd_point N by K matrix of pointwise log predictive densities. Each
#'   row must have at least one finite value.
#' @param tol Tolerance for `max_k log(g_k)`.
#' @param max_add Largest number of models added in one outer iteration. Each
#'   added model has `log(g_k) > tol` and starts at weight 0.001.
#' @return Vector of K weights that sums to 1.
#'
stacking_weights_active_set <- function(lpd_point, tol = 1e-8, max_add = 10) {
  N <- nrow(lpd_point)
  K <- ncol(lpd_point)

  # start with the best model, and add the best model of each observation
  # that no model in the set can predict, so that every row has a finite value
  active <- which.max(colSums(pmax(lpd_point, -1e300)))
  uncovered <- which(rowSums(is.finite(lpd_point[, active, drop = FALSE])) == 0)
  while (length(uncovered) > 0) {
    active <- c(active, which.max(lpd_point[uncovered[1], ]))
    uncovered <- which(
      rowSums(is.finite(lpd_point[, active, drop = FALSE])) == 0
    )
  }
  w <- rep(1 / length(active), length(active))

  converged <- FALSE
  for (outer in seq_len(K + 100)) {
    w_res <- stacking_newton(w, lpd_point[, active, drop = FALSE], tol = tol)
    keep <- w_res > 0
    active <- active[keep]
    w <- w_res[keep] / sum(w_res[keep])

    log_g <- stacking_log_g(lpd_point, active, w)
    k_new <- which.max(log_g)
    if (log_g[k_new] <= tol) {
      converged <- TRUE
      break
    }
    if (k_new %in% active) {
      # Newton stopped before the optimum of the active set
      break
    }
    # add the models with the largest g_k, at most max_add at a time; with many
    # nonzero weights this needs fewer outer iterations than one model at a time
    new <- setdiff(which(log_g > tol), active)
    new <- new[order(log_g[new], decreasing = TRUE)]
    new <- new[seq_len(min(length(new), max_add))]
    active <- c(active, new)
    w <- c(w * (1 - 1e-3 * length(new)), rep(1e-3, length(new)))
  }

  if (!converged) {
    # after the last iteration, the weights can include a model that was added
    # without Newton steps
    gap <- N * max(stacking_log_g(lpd_point, active, w))
    warning(
      "The stacking weights may not be optimal. ",
      "The log score can be up to ", signif(gap, 2), " below its maximum.",
      call. = FALSE
    )
  }

  out <- numeric(K)
  out[active] <- w
  out
}

#' Log of the stacking optimality measure g_k for all models
#'
#' @noRd
#' @param lpd_point N by K matrix of pointwise log predictive densities.
#' @param active Indices of the models with positive weight.
#' @param w Weights of the models in `active`, summing to 1.
#' @return Vector of K values `log(g_k)`, with
#'   `g_k = mean_i exp(lpd_point[i, k]) / p_i`. The maximum is 0 at the optimum.
#'
stacking_log_g <- function(lpd_point, active, w) {
  log_p <- matrixStats::rowLogSumExps(
    sweep(lpd_point[, active, drop = FALSE], 2, log(w), "+")
  )
  matrixStats::colLogSumExps(lpd_point - log_p) - log(nrow(lpd_point))
}

#' Negative log score of a combination of models
#'
#' @noRd
#' @param w Weights, non-negative and summing to 1.
#' @param lpd N by length(w) matrix of log predictive densities.
#' @return `-sum_i log(sum_k w_k exp(lpd[i, k]))`.
#'
stacking_neg_log_score <- function(w, lpd) {
  -sum(matrixStats::rowLogSumExps(sweep(lpd, 2, log(w), "+")))
}

#' Newton steps for the stacking weights of a fixed set of models
#'
#' @noRd
#' @param w Starting weights, positive and summing to 1.
#' @param lpd N by |active| matrix of log predictive densities. Each row has at
#'   least one finite value.
#' @param tol Tolerance of `stacking_weights_active_set()`. Newton stops when
#'   every `g_k` of the set is within `tol / 10` of 1.
#' @return Weights, non-negative and summing to 1. Models whose weight reaches
#'   zero leave the set, and their weights in the result are zero.
#'
stacking_newton <- function(w, lpd, tol, max_iter = 100) {
  out <- numeric(length(w))
  in_set <- seq_along(w)
  for (iter in seq_len(max_iter)) {
    n_active <- length(w)
    if (n_active == 1) {
      break
    }
    lpd_set <- lpd[, in_set, drop = FALSE]
    # scale each row by its maximum over the current set, so that the
    # densities are at most 1 and the largest one in each row is 1
    dens <- exp(lpd_set - matrixStats::rowMaxs(lpd_set))
    p <- drop(dens %*% w)
    grad <- -colSums(dens / p)
    hess <- crossprod(dens / p)
    # a small ridge keeps the system solvable for identical models
    hess <- hess + diag(1e-12 * mean(diag(hess)), n_active)
    kkt <- rbind(cbind(hess, 1), c(rep(1, n_active), 0))
    step <- tryCatch(
      solve(kkt, c(-grad, 0))[seq_len(n_active)],
      error = function(e) NULL
    )
    if (is.null(step)) {
      break
    }
    # at the optimum of this set, g_k = -grad_k / N is 1 for every model
    if (max(abs(-grad / nrow(lpd) - 1)) < tol / 10) {
      break
    }
    f0 <- stacking_neg_log_score(w, lpd_set)
    decrement <- -sum(grad * step)
    # close to the optimum the change of the log score is below its rounding
    # error, and the line search cannot see it; take the Newton step there
    near_optimum <- decrement < 1e-10 * max(1, abs(f0))
    # largest step that keeps the weights non-negative
    ratio <- ifelse(step < 0, -w / step, Inf)
    t_max <- min(1, ratio)
    if (t_max < 1e-10) {
      # a weight is almost zero and the step would make it negative: the
      # optimum of this set has that weight at zero
      w[ratio < 1e-10] <- 0
    } else {
      t <- t_max
      repeat {
        w_new <- pmax(w + t * step, 0)
        f_new <- stacking_neg_log_score(w_new, lpd_set)
        # a step to the boundary sets a weight to zero; accept it if the log
        # score does not get worse beyond rounding
        to_boundary <- t == t_max && t_max < 1 &&
          f_new <= f0 + 1e-12 * max(1, abs(f0))
        if (near_optimum || f_new <= f0 - 1e-4 * t * decrement ||
            to_boundary || t < 1e-10) {
          break
        }
        t <- t / 2
      }
      if (t < 1e-10) {
        break
      }
      w <- w_new
    }
    w <- w / sum(w)
    # models whose weight reached zero leave the set
    keep <- w > 0
    in_set <- in_set[keep]
    w <- w[keep] / sum(w[keep])
  }
  out[in_set] <- w
  out
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



# internal ----------------------------------------------------------------
# `-Inf` is a valid zero predictive density; `NA`, `NaN` and `+Inf` make
# stacking fail in the optimizer and pseudo-BMA return invalid weights.
validate_lpd_point <- function(lpd_point) {
  if (anyNA(lpd_point) || any(lpd_point == Inf)) {
    stop("All values in 'lpd_point' must be finite or -Inf.")
  }
  invisible(lpd_point)
}

#' Generate dirichlet simulations, rewritten version
#' @importFrom stats rgamma
#' @noRd
dirichlet_rng <- function(n, alpha) {
  K <- length(alpha)
  gamma_sim <- matrix(rgamma(K * n, alpha), ncol = K, byrow = TRUE)
  gamma_sim / rowSums(gamma_sim)
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
#' has the same number of observations (columns). The number of posterior
#' draws (rows) can differ, because the weights use only the pointwise elpd.
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
  if (length(unique(sapply(log_lik_list, ncol))) != 1) {
    stop("Each log-likelihood matrix must have the same number of ",
         "observations (columns).", call. = FALSE)
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
  # the pointwise values of a subsampled object cover only its subsample, so
  # the objects would not have values for the same observations
  if (any(sapply(psis_loo_list, inherits, "psis_loo_ss"))) {
    stop("loo_model_weights() does not support objects from loo_subsample(), ",
         "because they have pointwise values only for the subsampled ",
         "observations. Use loo() for all observations.", call. = FALSE)
  }

  # the weights use only the pointwise elpd, so the number of posterior draws
  # can differ between the objects
  n_obs <- sapply(psis_loo_list, function(x) nrow(x$pointwise))
  if (length(unique(n_obs)) != 1) {
    stop("Each object in the list must have the same number of observations.",
         call. = FALSE)
  }
  invisible(TRUE)
}
