#' Iterative Pareto smoothed importance sampling
#'
#' Iterative importance weighting using Pareto Smoothed Importance Sampling  (PSIS).
#' The number of iterations is a random variable and is determined by the
#' shape parameter \eqn{k} of a generalized Pareto distribution. The algorithm
#' stops once \eqn{k < 1/2} (see Details).
#'
#' @export
#' @param start Named list with components \code{log_p} (target), \code{log_g}
#'   (proposal), \code{draws} (parameter draws from approximate posterior).
#' @param stanfit \code{\link[=stanfit-class]{stanfit}} object from which to
#'   obtain \code{log_prob()} and \code{unconstrain_pars()} functions.
#' @param ... Optional PSIS tuning parameters passed to \code{\link{psislw}}.
#' @param control Parameters controlling the iterative process.
#'
#' @return A named list with components
#' \itemize{
#' \item \code{mu} Mean vectors
#' \item \code{Sigma} Covariance matrices
#' \item \code{khat} Pareto shape estimates
#' \item \code{coef_lg} Slopes from regression of log_p on log_g
#' }
#' Each component is a list of length equal to the number of iterations.
#'
#' @details
#' We start by computing log weights from \code{start$log_p} and
#' \code{start$log_g}, smoothing them with PSIS, and computing the weighted mean
#' vector and covariance matrix using the smoothed weights and
#' \code{start$draws}. We then take this first weighted mean vector and
#' covariance matrix pair and use them to start the iterative process:
#'
#' \enumerate{
#' \item Draw from multivariate normal approximation.
#' \item Evaluate log_p and log_g.
#' \item Apply PSIS to smooth log weights (if \code{smooth_weights} = TRUE).
#' \item Calculate weighted mean vector and covariance matrix.
#' }
#'
#' \subsection{Stopping rule}{
#' The iterative algorithm will stop once the estimate of the Pareto shape
#' parameter \eqn{k} (estimated during PSIS smoothing) is below \eqn{1/2} or
#' when \code{max_iter} is reached.
#'
#'\itemize{
#' \item If \eqn{k < 1/2} the variance of the raw importance ratios is finite,
#'  the central limit theorem holds, and the estimate converges quickly.
#' \item If \eqn{k} is between \eqn{1/2} and \eqn{1} the variance of the raw importance
#'  ratios is infinite but the mean exists, the generalized central limit theorem
#'  for stable distributions holds, and the convergence of the estimate is
#'  slower. The variance of the PSIS estimate is finite but may be large.
#' \item If \eqn{k > 1} the variance and the mean of the raw ratios distribution
#'  do not exist. The variance of the PSIS estimate is finite but may be large.
#'}
#'}
#'
iterate_psis <- function(start, stanfit, control = iter_control(), ...) {
  mu <- Sigma <- coef_lg <- khat <- list()
  psis1 <- psislw(start$log_p - start$log_g, ...)
  starting_mean_and_var <- weighted_mean_and_var(start$draws, lw = psis1$lw_smooth)

  mu[[1L]] <- starting_mean_and_var$mean
  Sigma[[1L]] <- starting_mean_and_var$var
  khat[[1L]] <- psis1$pareto_k
  coef_lg[[1L]] <- coef(lm(start$log_p ~ start$log_g))[2L]
  skeleton <- get_inits(stanfit)[[1L]]

  if (control$verbose == TRUE) {
    if (control$smooth_weights == FALSE) message("smooth_weights = FALSE, using raw weights")
    .khat_msg(khat[[1]], iter = 1)
  }
  if (control$plot == TRUE) .initialize_khat_plot(khat[[1]], max_iter = control$max_iter)

  # Begin iterations
  for (n in 2:control$max_iter) {
    mvn_draws <- mvtnorm::rmvnorm(control$ndraws, mu[[n-1]], Sigma[[n-1]])
#     upars <- t(apply(mvn_draws, 1L, FUN = function(theta) {
#       unconstrain_pars(stanfit, relist(theta, skeleton))
#     }))
    log_p <- apply(mvn_draws, 1L, FUN = function(u) {
      log_prob(stanfit, u, adjust_transform = FALSE) # does it matter if this is false or true?
    })
    log_g <- apply(mvn_draws, 1L, FUN = function(u) {
      mvtnorm::dmvnorm(u, mu[[n-1]], Sigma[[n-1]], log = TRUE)
    })
    psis_n <- psislw(lw = log_p - log_g, ...)
    lw <- if (control$smooth_weights == TRUE)
      psis_n$lw_smooth else lw_normalize(log_p - log_g)

    next_mu_sigma <- weighted_mean_and_var(mvn_draws, lw = lw)
    mu[[n]] <- next_mu_sigma$mean
    Sigma[[n]] <- next_mu_sigma$var
    coef_lg[[n]] <- coef(lm(log_p ~ log_g))["log_g"]
    khat[[n]] <- k_n <- psis_n$pareto_k

    if (control$plot == TRUE) {
      clr <- .khat_clr(k_n)
      segments(x0 = n-1, y0 = khat[[n-1]], x1 = n, y1 = k_n)
      points(n, k_n, col = clr, pch = 19)
    }
    if (control$verbose == TRUE)
      .khat_msg(k_n, iter = n)
    if (k_n < 0.5) {
      if (control$verbose == TRUE)
        message("Stopping... khat below 0.5 after ", n, " iterations.")
      break
    }
  }
  if (n == control$max_iter)
    warning("Stopping... max_iter reached. khat is still above 0.5 after ",
            control$max_iter, " iterations.")

  nlist(mu, Sigma, khat, coef_lg)
}

#' @rdname iterate_psis
#' @export
#' @param smooth_weights If \code{TRUE} (the default) \code{\link[psislw]{PSIS}}
#'   is used to smooth the raw weights. If \code{FALSE} the raw weights are
#'   used.
#' @param max_iter The maximum allowed number of iterations. The algorithm will
#'   stop after \code{max_iter} iterations even if \eqn{k > 1/2}.
#' @param ndraws The size of the multivariate normal sample to generate at each
#'   iteration.
#' @param verbose Should status updates be printed?
#' @param plot Should pareto k estimates be plotted as each iteration finishes?
#'
iter_control <- function(smooth_weights = TRUE, max_iter = 20, ndraws = 4000,
                         verbose = TRUE, plot = interactive()) {
  nlist(smooth_weights, max_iter, ndraws, verbose, plot)
}

#' @rdname iterate_psis
#' @param mu_approx,Sigma_approx,mu_true,Sigma_true Mean and covariance matrices
#'   for approximating distribution and comparison (true) distribution.
#'
#' @importFrom monomvn kl.norm
KLdivergence <- function(mu_approx, Sigma_approx, mu_true, Sigma_true) {
  monomvn::kl.norm(mu1 = mu_approx,
                   S1 = Sigma_approx,
                   mu2 = mu_true,
                   S2 = Sigma_true)
}

.initialize_khat_plot <- function(k1, max_iter) {
  plot(1, k1, xlab = "Iteration", ylab = "khat",
       xlim = c(1, max_iter), ylim = c(0, 1.5),
       col = .khat_clr(k1), pch = 19)
  abline(h = c(0.5, 1), lty = 2, col = "darkgray")
}

.khat_msg <- function(k, iter, digits = 3) {
  message("iteration: 1, khat = ", round(k, digits))
}
.khat_clr <- function(k, clrs = c("blue", "purple", "red")) {
  ifelse(k < 0.5, clrs[1], ifelse(k < 1, clrs[2], clrs[3]))
}

exp_norm <- function(lw) {
  exp(lw - max(lw))
}

efficiency <- function(w, log = FALSE) {
  S <- length(w)
  if(log) {
    w <- w - matrixStats::logSumExp(w) + log(S)
    efficiency <- exp(log(S) - matrixStats::logSumExp(2 * w))
  } else {
    w_norm <- w / mean(w)
    efficiency <- S / sum(w_norm^2)
  }
  return(efficiency)
}

## calculates reweighted mean and variance given log-weights
weighted_mean_and_var <- function(X, lw, ...) {
  UseMethod("weighted_mean_and_var")
}
weighted_mean_and_var.data.frame <- function(X, lw, ...) {
  weighted_mean_and_var.matrix(X = as.matrix(X), lw = lw, ...)
}
weighted_mean_and_var.matrix <- function(X, lw, ...) {
  S <- length(lw)
  ess <- S * efficiency(lw, log=TRUE)
  if (ess <= 1) {
    ## for values <= 1 the ess correction of the variance below
    ## collapses, hence correct for this
    warning(paste("ESS =", ess, "artifically increased to", 1.1))
    ess <- 1.1
  }
  w <- exp_norm(lw)
  if (is.matrix(X)) {
    K <- ncol(X)
    mu <- apply(X, 2L, weighted.mean, w = w)
    Sigma  <- array(NA, c(K,K))
    for (i in 1:K) {
      for (j in 1:K) {
        tmp <- (X[,i] - mu[i]) * (X[,j] - mu[j])
        Sigma[i,j] <- weighted.mean(tmp, w)
      }
    }
    Sigma <- Sigma * ess / (ess - 1)
    sds  <- sqrt(diag(Sigma))
  }
  list(mean = mu, var = Sigma, sd = sds)
}

weighted_mean_and_var.vector <- function(X, lw, ...) {
  S <- length(lw)
  ess <- S * efficiency(lw, log=TRUE)
  if (ess <= 1) {
    ## for values <= 1 the ess correction of the variance below
    ## collapses, hence correct for this
    warning(paste("ESS =", ess, "artifically increased to", 1.1))
    ess <- 1.1
  }
  w <- exp_norm(lw)
  mu <- weighted.mean(X, w)
  var  <- weighted.mean((X - mu)^2, w)
  var  <- var * ess / (ess - 1)
  sd   <- sqrt(var)
  list(mean = mu, var = var, sd = sd)
}

