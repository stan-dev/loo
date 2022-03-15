#' Diagnostics for Pareto smoothed importance sampling (PSIS)
#'
#' Print a diagnostic table summarizing the estimated Pareto shape parameters
#' and PSIS effective sample sizes, find the indexes of observations for which
#' the estimated Pareto shape parameter \eqn{k} is larger than some
#' `threshold` value, or plot observation indexes vs. diagnostic estimates.
#' The **Details** section below provides a brief overview of the
#' diagnostics, but we recommend consulting Vehtari, Gelman, and Gabry (2017)
#' and Vehtari, Simpson, Gelman, Yao, and Gabry (2019) for full details.
#'
#' @name pareto-k-diagnostic
#' @param x An object created by [loo()] or [psis()].
#' @param threshold For `pareto_k_ids()`, `threshold` is the minimum \eqn{k}
#'   value to flag (default is `0.5`). For `mcse_loo()`, if any \eqn{k}
#'   estimates are greater than `threshold` the MCSE estimate is returned as
#'   `NA` (default is `0.7`). See **Details** for the motivation behind these
#'   defaults.
#'
#' @details
#' The reliability and approximate convergence rate of the PSIS-based estimates
#' can be assessed using the estimates for the shape parameter \eqn{k} of the
#' generalized Pareto distribution:
#' * If \eqn{k < 0.5} then the distribution of raw importance ratios has
#'   finite variance and the central limit theorem holds. However, as \eqn{k}
#'   approaches \eqn{0.5} the RMSE of plain importance sampling (IS) increases
#'   significantly while PSIS has lower RMSE.
#'
#' * If \eqn{0.5 \leq k < 1}{0.5 <= k < 1} then the variance of the raw
#'   importance ratios is infinite, but the mean exists. TIS and PSIS estimates
#'   have finite variance by accepting some bias. The convergence of the
#'   estimate is slower with increasing \eqn{k}.
#'   If \eqn{k} is between 0.5 and approximately 0.7 then we observe practically
#'   useful convergence rates and Monte Carlo error estimates with PSIS (the
#'   bias of TIS increases faster than the bias of PSIS). If \eqn{k > 0.7} we
#'   observe impractical convergence rates and unreliable Monte Carlo error
#'   estimates.
#'
#' * If \eqn{k \geq 1}{k >= 1} then neither the variance nor the mean of
#'   the raw importance ratios exists. The convergence rate is close to
#'   zero and bias can be large with practical sample sizes.
#'
#' \subsection{What if the estimated tail shape parameter \eqn{k} exceeds
#' \eqn{0.5}}{ Importance sampling is likely to work less well if the
#' marginal posterior \eqn{p(\theta^s | y)} and LOO posterior
#' \eqn{p(\theta^s | y_{-i})} are very different, which is more likely to
#' happen with a non-robust model and highly influential observations.
#' If the estimated tail shape parameter \eqn{k} exceeds \eqn{0.5}, the user
#' should be warned. (Note: If \eqn{k} is greater than \eqn{0.5}
#' then WAIC is also likely to fail, but WAIC lacks its own diagnostic.)
#' In practice, we have observed good performance for values of \eqn{k} up to 0.7.
#' When using PSIS in the context of approximate LOO-CV, we recommend one of
#' the following actions when \eqn{k > 0.7}:
#'
#' * With some additional computations, it is possible to transform the MCMC
#'   draws from the posterior distribution to obtain more reliable importance
#'   sampling estimates. This results in a smaller shape parameter \eqn{k}.
#'   See [loo_moment_match()] for an example of this.
#'
#' * Sampling directly from \eqn{p(\theta^s | y_{-i})} for the problematic
#'   observations \eqn{i}, or using \eqn{k}-fold cross-validation will generally
#'   be more stable.
#'
#' * Using a model that is more robust to anomalous observations will
#'   generally make approximate LOO-CV more stable.
#'
#' }
#'
#' \subsection{Observation influence statistics}{ The estimated shape parameter
#' \eqn{k} for each observation can be used as a measure of the observation's
#' influence on posterior distribution of the model. These can be obtained with
#' `pareto_k_influence_values()`.
#' }
#'
#' \subsection{Effective sample size and error estimates}{ In the case that we
#' obtain the samples from the proposal distribution via MCMC the **loo**
#' package also computes estimates for the Monte Carlo error and the effective
#' sample size for importance sampling, which are more accurate for PSIS than
#' for IS and TIS (see Vehtari et al (2019) for details). However, the PSIS
#' effective sample size estimate will be
#'  **over-optimistic when the estimate of \eqn{k} is greater than 0.7**.
#' }
#'
#' @seealso
#'  * [psis()] for the implementation of the PSIS algorithm.
#'  * The [FAQ page](https://mc-stan.org/loo/articles/online-only/faq.html) on
#'    the __loo__ website for answers to frequently asked questions.
#'
#' @template loo-and-psis-references
#'
NULL

#' @rdname pareto-k-diagnostic
#' @export
#' @return `pareto_k_table()` returns an object of class
#'   `"pareto_k_table"`, which is a matrix with columns `"Count"`,
#'   `"Proportion"`, and `"Min. n_eff"`, and has its own print method.
#'
pareto_k_table <- function(x) {
  k <- pareto_k_values(x)
  n_eff <- try(psis_n_eff_values(x), silent = TRUE)
  if (inherits(n_eff, "try-error")) {
    n_eff <- rep(NA, length(k))
  }

  kcut <- k_cut(k)
  min_n_eff <- min_n_eff_by_k(n_eff, kcut)
  count <- table(kcut)
  out <- cbind(
    Count = count,
    Proportion = prop.table(count),
    "Min. n_eff" = min_n_eff
  )
  structure(out, class = c("pareto_k_table", class(out)))
}

#' @export
print.pareto_k_table <- function(x, digits = 1, ...) {
  count <- x[, "Count"]

  if (sum(count[2:4]) == 0) {
    cat("\nAll Pareto k estimates are good (k < 0.5).\n")
  } else {
    tab <- cbind(
      " " = rep("", 4),
      " " = c("(good)", "(ok)", "(bad)", "(very bad)"),
      "Count" = .fr(count, 0),
      "Pct.   " = paste0(.fr(100 * x[, "Proportion"], digits), "%"),
      "Min. n_eff" = round(x[, "Min. n_eff"])
    )
    tab2 <- rbind(tab)
    cat("Pareto k diagnostic values:\n")
    rownames(tab2) <- format(rownames(tab2), justify = "right")
    print(tab2, quote = FALSE)

    if (sum(count[3:4]) == 0)
      cat("\nAll Pareto k estimates are ok (k < 0.7).\n")

    invisible(x)
  }
}

#' @rdname pareto-k-diagnostic
#' @export
#' @return `pareto_k_ids()` returns an integer vector indicating which
#' observations have Pareto \eqn{k} estimates above `threshold`.
#'
pareto_k_ids <- function(x, threshold = 0.5) {
  k <- pareto_k_values(x)
  which(k > threshold)
}

#' @rdname pareto-k-diagnostic
#' @export
#' @return `pareto_k_values()` returns a vector of the estimated Pareto
#'   \eqn{k} parameters. These represent the reliability of sampling.
pareto_k_values <- function(x) {
  k <- x$diagnostics[["pareto_k"]]
  if (is.null(k)) {
    # for compatibility with objects from loo < 2.0.0
    k <- x[["pareto_k"]]
  }
  if (is.null(k)) {
    stop("No Pareto k estimates found.", call. = FALSE)
  }
  return(k)
}

#' @rdname pareto-k-diagnostic
#' @export
#' @return `pareto_k_influence_values()` returns a vector of the estimated Pareto
#'   \eqn{k} parameters. These represent influence of the observations on the
#'   model posterior distribution.
pareto_k_influence_values <- function(x) {
  if ("influence_pareto_k" %in% colnames(x$pointwise)) {
    k <- x$pointwise[,"influence_pareto_k"]
  }
  else {
    stop("No Pareto k influence estimates found.", call. = FALSE)
  }
  return(k)
}

#' @rdname pareto-k-diagnostic
#' @export
#' @return `psis_n_eff_values()` returns a vector of the estimated PSIS
#'   effective sample sizes.
psis_n_eff_values <- function(x) {
  n_eff <- x$diagnostics[["n_eff"]]
  if (is.null(n_eff)) {
    stop("No PSIS n_eff estimates found.", call. = FALSE)
  }
  return(n_eff)
}

#' @rdname pareto-k-diagnostic
#' @export
#' @return `mcse_loo()` returns the Monte Carlo standard error (MCSE)
#'   estimate for PSIS-LOO. MCSE will be NA if any Pareto \eqn{k} values are
#'   above `threshold`.
#'
mcse_loo <- function(x, threshold = 0.7) {
  stopifnot(is.psis_loo(x))
  if (any(pareto_k_values(x) > 0.7, na.rm = TRUE)) {
    return(NA)
  }
  mc_var <- x$pointwise[, "mcse_elpd_loo"]^2
  sqrt(sum(mc_var))
}

#' @rdname pareto-k-diagnostic
#' @aliases plot.loo
#' @export
#' @param label_points,... For the `plot()` method, if `label_points` is
#'   `TRUE` the observation numbers corresponding to any values of \eqn{k}
#'   greater than 0.5 will be displayed in the plot. Any arguments specified in
#'   `...` will be passed to [graphics::text()] and can be used
#'   to control the appearance of the labels.
#' @param diagnostic For the `plot` method, which diagnostic should be
#'   plotted? The options are `"k"` for Pareto \eqn{k} estimates (the
#'   default) or `"n_eff"` for PSIS effective sample size estimates.
#' @param main For the `plot()` method, a title for the plot.
#'
#' @return The `plot()` method is called for its side effect and does not
#'   return anything. If `x` is the result of a call to [loo()]
#'   or [psis()] then `plot(x, diagnostic)` produces a plot of
#'   the estimates of the Pareto shape parameters (`diagnostic = "k"`) or
#'   estimates of the PSIS effective sample sizes (`diagnostic = "n_eff"`).
#'
plot.psis_loo <- function(x,
                          diagnostic = c("k", "n_eff"),
                          ...,
                          label_points = FALSE,
                          main = "PSIS diagnostic plot") {
  diagnostic <- match.arg(diagnostic)
  k <- pareto_k_values(x)
  k[is.na(k)] <- 0  # FIXME when reloo is changed to make NA k values -Inf
  k_inf <- !is.finite(k)
  if (any(k_inf)) {
    warning(signif(100 * mean(k_inf), 2),
            "% of Pareto k estimates are Inf/NA/NaN and not plotted.")
  }

  if (diagnostic == "n_eff") {
    n_eff <- psis_n_eff_values(x)
  } else {
    n_eff <- NULL
  }

  plot_diagnostic(
    k = k,
    n_eff = n_eff,
    ...,
    label_points = label_points,
    main = main
  )
}

#' @export
#' @noRd
#' @rdname pareto-k-diagnostic
plot.loo <- plot.psis_loo

#' @export
#' @rdname pareto-k-diagnostic
plot.psis <- function(x, diagnostic = c("k", "n_eff"), ...,
                      label_points = FALSE,
                      main = "PSIS diagnostic plot") {
  plot.psis_loo(x, diagnostic = diagnostic, ...,
                label_points = label_points, main = main)
}



# internal ----------------------------------------------------------------

plot_diagnostic <-
  function(k,
           n_eff = NULL,
           ...,
           label_points = FALSE,
           main = "PSIS diagnostic plot") {
    use_n_eff <- !is.null(n_eff)
    graphics::plot(
      x = if (use_n_eff) n_eff else k,
      xlab = "Data point",
      ylab = if (use_n_eff) "PSIS n_eff" else "Pareto shape k",
      type = "n",
      bty = "l",
      yaxt = "n",
      main = main
    )
    graphics::axis(side = 2, las = 1)

    in_range <- function(x, lb_ub) {
      x >= lb_ub[1L] & x <= lb_ub[2L]
    }

    if (!use_n_eff) {
      krange <- range(k, na.rm = TRUE)
      breaks <- c(0, 0.5, 0.7, 1)
      hex_clrs <- c("#C79999", "#A25050", "#7C0000")
      ltys <- c(3, 4, 2, 1)
      for (j in seq_along(breaks)) {
        val <- breaks[j]
        if (in_range(val, krange))
          graphics::abline(
            h = val,
            col = ifelse(val == 0, "darkgray", hex_clrs[j - 1]),
            lty = ltys[j],
            lwd = 1
          )
      }
    }

    breaks <- c(-Inf, 0.5, 1)
    hex_clrs <- c("#6497b1", "#005b96", "#03396c")
    clrs <- ifelse(
      in_range(k, breaks[1:2]),
      hex_clrs[1],
      ifelse(in_range(k, breaks[2:3]), hex_clrs[2], hex_clrs[3])
    )
    if (all(k < 0.5) || !label_points) {
      graphics::points(x = if (use_n_eff) n_eff else k,
                       col = clrs, pch = 3, cex = .6)
      return(invisible())
    } else {
      graphics::points(x = if (use_n_eff) n_eff[k < 0.5] else k[k < 0.5],
                       col = clrs[k < 0.5], pch = 3, cex = .6)
      sel <- !in_range(k, breaks[1:2])
      dots <- list(...)
      txt_args <- c(
        list(
          x = seq_along(k)[sel],
          y = if (use_n_eff) n_eff[sel] else k[sel],
          labels = seq_along(k)[sel]
        ),
        if (length(dots)) dots
      )
      if (!("adj" %in% names(txt_args))) txt_args$adj <- 2 / 3
      if (!("cex" %in% names(txt_args))) txt_args$cex <- 0.75
      if (!("col" %in% names(txt_args))) txt_args$col <- clrs[sel]

      do.call(graphics::text, txt_args)
    }
  }


#' Convert numeric Pareto k values to a factor variable.
#'
#' @noRd
#' @param k Vector of Pareto k estimates.
#' @return A factor variable (the same length as k) with 4 levels.
#'
k_cut <- function(k) {
  cut(
    k,
    breaks = c(-Inf, 0.5, 0.7, 1, Inf),
    labels = c("(-Inf, 0.5]", "(0.5, 0.7]", "(0.7, 1]", "(1, Inf)")
  )
}

#' Calculate the minimum PSIS n_eff within groups defined by Pareto k values
#'
#' @noRd
#' @param n_eff Vector of PSIS n_eff estimates.
#' @param kcut Factor returned by the k_cut() function.
#' @return Vector of length `nlevels(kcut)` containing the minimum n_eff within
#'   each k group. If there are no k values in a group the corresponding element
#'   of the returned vector is NA.
min_n_eff_by_k <- function(n_eff, kcut) {
  n_eff_split <- split(n_eff, f = kcut)
  n_eff_split <- sapply(n_eff_split, function(x) {
    # some k groups might be empty.
    # split gives numeric(0) but replace with NA
    if (!length(x)) NA else x
  })
  sapply(n_eff_split, min)
}
