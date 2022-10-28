#' R2
#'
#' The `loo_r2()` method can compute the pointwise R2 metric
#'
#' @export
#' @param y A numeric vector of observations. Length should be equal to the
#'          number of rows in `yrep`.
#' @param yrep A numeric matrix of predictions.
#' @param log_lik A matrix of pointwise log-likelihoods. Should be of same
#'     dimension as `yrep`.
#' @param r_eff A Vector of relative effective sample size estimates containing
#'     one element per observation. See [psis()] for more details.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{`pointwise`}{
#'   Pointwise components of the R2
#'   }
#'   \item{`estimate`}{
#'   R2 estimate
#'   }
#'  \item{`se`}{
#'   Standard error of the R2
#'   }
#'  }
#' @export
#' @examples
#' \donttest{
#' if (requireNamespace("rstanarm", quietly = TRUE)) {
#'   # Use rstanarm package to quickly fit a model and get both a log-likelihood
#'   # matrix and draws from the posterior predictive distribution
#'   library("rstanarm")
#'   data("mtcars")
#'
#'   fit <- stan_glm(mpg ~ cyl + hp + wt, data = mtcars, refresh = 0)
#'   ll <- log_lik(fit)
#'   r_eff <- relative_eff(exp(-ll), chain_id = rep(1:4, each = 1000))
#'
#'   yrep <- posterior_predict(fit)
#'   # Leave-one-out R2
#'   r2 <- loo_r2(
#'     y = mtcars$mpg,
#'     yrep = yrep,
#'     log_lik = ll,
#'     r_eff = r_eff
#'   )
#' }
#' }
loo_r2 <- function(y, yrep, log_lik, r_eff, ...) {
  UseMethod("loo_r2")
}

#' @rdname loo_r2
#' @export
loo_r2.default <-
  function(y,
           yrep,
           log_lik,
           r_eff = NULL) {
    stopifnot(
      is.numeric(y),
      is.numeric(yrep),
      identical(ncol(yrep), length(y)),
      identical(dim(yrep), dim(log_lik))
    )
    psis_object <- psis(-log_lik, r_eff = r_eff)
    pointwise <- .r2(
      y = y,
      yrep = yrep,
      weights = exp(psis_object$log_weights)
    )
    list(
      estimate = sum(pointwise),
      pointwise = pointwise,
      se = sqrt(length(pointwise) * var(pointwise)),
      diagnostics = psis_object$diagnostics
    )
  }

# internal ----------------------------------------------------------------

.r2 <- function(y, yrep, weights = NULL) {
  ss_y <- sum((y - mean(y))^2)
  pointwise_loo_r2 <- vector(mode = "numeric", length = length(y))

  if (is.null(weights)) {
    for (n in seq_along(pointwise_loo_r2)) {
      ss_e <- sum((y[n] - yrep[, n])^2)
      pointwise_loo_r2[[n]] <- 1 / length(y) - ss_e / ss_y
    }
  } else {
    for (n in seq_along(pointwise_loo_r2)) {
      ss_e <- sum((weights[, n] * (y[n] - yrep[, n])^2) / sum(weights[, n]))
      pointwise_loo_r2[[n]] <- 1 / length(y) - ss_e / ss_y
    }
  }
  pointwise_loo_r2
}
