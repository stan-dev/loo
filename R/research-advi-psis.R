#' Iterative Pareto smoothed importance sampling
#'
#' Iterative importance weighting using Pareto Smoothed Importance Sampling  (PSIS) \eqn{k} estimate is below 0.5.
#'
#' @export
#' @param start Named list with components \code{log_p}, \code{log_g}, \code{draws}.
#' @param stanfit \code{\link[=stanfit-class]{stanfit}} object from which to obtain \code{log_prob()} and \code{unconstrain_pars()} functions.
#' @param smooth_weights If \code{TRUE} (the default) \code{\link[psislw]{PSIS}} is used to smooth the raw weights. If FALSE the raw weights are used.
#' @param ... Optional PSIS tuning parameters passed to \code{\link{psislw}}.
#' @param control Parameters controlling the iterative process.
#'
iterate_psis <- function(start, stanfit, smooth_weights = TRUE, ...,
                         control = iter_control()) {
  mu <- Sigma <- coef_lg <- khat <- list()
  psis1 <- psislw(start$log_p - start$log_g)
  starting_mean_and_var <- weighted_mean_and_var(start$draws, lw = psis1$lw_smooth)

  mu[[1]] <- starting_mean_and_var$mean
  Sigma[[1]] <- starting_mean_and_var$var
  khat[[1]] <- psis1$pareto_k
  coef_lg[[1]] <- coef(lm(start$log_p ~ start$log_g))[2]
  skeleton <- get_inits(stanfit)[[1]]

  if (control$verbose) message("iteration: 1, khat = ", round(khat[[1]], 3))

  if (control$plot) {
    clr <- ifelse(khat[[1]] < 0.5, "blue", ifelse(khat[[1]] < 1, "purple", "red"))
    plot(1, khat[[1]], xlab = "Iteration", ylab = "khat", col = clr, pch = 19,
         xlim = c(1, control$max_iter), ylim = c(-1, 2))
    abline(h = c(0.5, 1), lty = 2, col = "maroon")
  }
  for (n in 2:control$max_iter) {
    mvn_draws <- mvtnorm::rmvnorm(control$ndraws, mu[[n-1]], Sigma[[n-1]])
    upars <- t(apply(mvn_draws, 1, FUN = function(theta) {
      unconstrain_pars(stanfit, relist(theta, skeleton))
    }))
    log_p <- apply(upars, 1, FUN = function(u) {
      log_prob(stanfit, u, adjust_transform = FALSE) # does it matter if this is false or true?
    })
    log_g <- apply(upars, 1, FUN = function(u) {
      mvtnorm::dmvnorm(u, mu[[n-1]], Sigma[[n-1]], log = TRUE)
    })
    psis_n <- psislw(lw = log_p - log_g, ...)
    if (smooth_weights) lw <- psis_n$lw_smooth
    else lw <- lw_normalize(log_p - log_g)
    next_mu_sigma <- weighted_mean_and_var(mvn_draws, lw = lw)
    mu[[n]] <- next_mu_sigma$mean
    Sigma[[n]] <- next_mu_sigma$var
    coef_lg[[n]] <- coef(lm(log_p ~ log_g))["log_g"]
    khat[[n]] <- psis_n$pareto_k

    if (control$plot) {
      clr <- ifelse(khat[[n]] < 0.5, "blue", ifelse(khat[[n]] < 1, "purple", "red"))
      segments(x0 = n-1, y0 = khat[[n-1]], x1 = n, y1 = khat[[n]])
      points(n, khat[[n]], col = clr, pch = 19)
    }
    if (control$verbose) message("iteration: ", n, ", khat = ", round(khat[[n]], 3))
    if (khat[[n]] < 0.5) {
      if (control$verbose) message("Stopping... khat below 0.5 after ", n, " iterations.")
      break
    }
  }
  if (n == control$max_iter)
    warning("khat is still above 0.5 after ", control$max_iter, " iterations.")

  nlist(mu, Sigma, khat, coef_lg)
}

#' @rdname iterate_psis
#' @export
#'
#' @param max_iter The maximum allowed number of iterations.
#' @param The number of multivariate normal draws to use at each iteration.
#' @param verbose Should status updates be printed?
#' @param plot Should pareto k estimates be plotted as each iteration finishes?
#'
iter_control <- function(max_iter = 20, ndraws = 4000,
                         verbose = TRUE, plot = interactive()) {
  nlist(max_iter, ndraws, verbose, plot)
}


exp_norm <- function(lw) {
  exp(lw - max(lw))
}

efficiency_weights <- function(w, log = FALSE) {
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
  ess <- S * efficiency_weights(lw, log=TRUE)
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
  ess <- S * efficiency_weights(lw, log=TRUE)
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

