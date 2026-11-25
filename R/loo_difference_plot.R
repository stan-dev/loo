#' Compare models across domains
#'
#' The LOO difference plot shows how the ELPD of two different models
#' changes when a predictor is varied. This can be useful for identifying
#' opportunities for model stacking or expansion. Pointwise differences
#' are computed as `loo_1 - loo_2`, so positive values indicate better
#' predictive performance for `loo_1`.
#'
#' @param y A vector of observations.
#' @param loo_1,loo_2 Objects returned by [loo()].
#' @param group An optional grouping variable with the same length as `y`.
#'   Points are colored according to group membership.
#' @param size,alpha Point size and opacity passed to [ggplot2::geom_point()].
#' @param jitter Amount of horizontal jitter passed to
#'   [ggplot2::position_jitter()].
#' @param sort_by_group If `TRUE`, observations are ordered by `group`
#'   and the x-axis is replaced by a sequential index. The supplied `y` values
#'   are therefore not used as x coordinates. Plotting by index can be useful
#'   when categories have very different sample sizes. To control the group
#'   order, supply `group` as a factor with levels in the desired order.
#' @param label_threshold Optional nonnegative threshold for labeling
#'   observations. Observations for which the absolute pointwise ELPD
#'   difference exceeds this value are labeled. If `NULL`, no observations
#'   are labeled.
#' @param labels Optional vector of labels with the same length as `y`, used
#'   for observations selected by `label_threshold`. If `NULL`, observation
#'   indices are used.
#'
#' @template bayesvis-reference
#'
#' @return A [ggplot2::ggplot()] object.
#'
#' @examples
#' # Artificial example
#' log_lik <- example_loglik_matrix()
#' shift <- seq(-0.5, 0.5, length.out = ncol(log_lik))
#' log_lik_2 <- sweep(log_lik, 2, shift, FUN = "+")
#'
#' loo_1 <- loo(log_lik)
#' loo_2 <- loo(log_lik_2)
#'
#' plot_loo_difference(
#'   seq_len(ncol(log_lik)),
#'   loo_1,
#'   loo_2
#' )
#'
#' # Label observations with large pointwise ELPD differences
#' plot_loo_difference(
#'   seq_len(ncol(log_lik)),
#'   loo_1,
#'   loo_2,
#'   label_threshold = 0.3
#' )
#'
#' # Create interspersed groups, then sort them in the plot
#' group <- rep(c("A", "A", "A", "B"), length.out = ncol(log_lik))
#'
#' plot_loo_difference(
#'   seq_len(ncol(log_lik)),
#'   loo_1,
#'   loo_2,
#'   group = group,
#'   sort_by_group = TRUE
#' )
#'
#' @export
plot_loo_difference <-
  function(
    y,
    loo_1,
    loo_2,
    group = NULL,
    size = 1,
    alpha = 1,
    jitter = 0,
    sort_by_group = FALSE,
    label_threshold = NULL,
    labels = NULL
  ) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop(
        "Please install 'ggplot2' to use `plot_loo_difference()`.",
        call. = FALSE
      )
    }

    checkmate::assert_flag(sort_by_group)
    loo_compare_checks(nlist(loo_1, loo_2))

    # elpd_diffs(a, b) computes b - a
    elpd_diff <- elpd_diffs(loo_2, loo_1)

    checkmate::assert_atomic_vector(
      y,
      len = length(elpd_diff)
    )

    if (!is.null(group)) {
      checkmate::assert_atomic_vector(
        group,
        len = length(y),
        any.missing = FALSE
      )
    }

    if (!is.null(labels)) {
      checkmate::assert_atomic_vector(
        labels,
        len = length(y)
      )

      if (is.null(label_threshold)) {
        stop(
          "`label_threshold` must be supplied when `labels` is supplied.",
          call. = FALSE
        )
      }
    }

    if (!is.null(label_threshold)) {
      checkmate::assert_number(
        label_threshold,
        lower = 0,
        finite = TRUE
      )
    }

    if (!is.null(label_threshold) && is.null(labels)) {
      labels <- seq_along(y)
    }

    if (sort_by_group) {
      if (is.null(group)) {
        stop(
          "`group` must be supplied when `sort_by_group = TRUE`.",
          call. = FALSE
        )
      }

      ordering <- order(group)
      elpd_diff <- elpd_diff[ordering]
      group <- group[ordering]

      if (!is.null(labels)) {
        labels <- labels[ordering]
      }

      y <- seq_along(elpd_diff)
    }

    plot_data <- data.frame(
      y = y,
      elpd_diff = elpd_diff
    )

    if (!is.null(group)) {
      plot_data$group <- factor(group)
    }

    if (!is.null(label_threshold)) {
      plot_data$labels <- ifelse(
        abs(plot_data$elpd_diff) > label_threshold,
        as.character(labels),
        ""
      )
    }

    jitter_position <- ggplot2::position_jitter(
      width = jitter,
      height = 0,
      seed = 1
    )

    plot <- ggplot2::ggplot(
      data = plot_data,
      mapping = ggplot2::aes(x = y, y = elpd_diff)
    ) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::labs(
        x = if (sort_by_group) "Index" else NULL,
        y = "Pointwise ELPD Difference (loo_1 - loo_2)"
      )

    if (is.null(group)) {
      plot <- plot +
        ggplot2::geom_point(
          position = jitter_position,
          alpha = alpha,
          size = size
        )
    } else {
      plot <- plot +
        ggplot2::geom_point(
          ggplot2::aes(color = group),
          position = jitter_position,
          alpha = alpha,
          size = size
        ) +
        ggplot2::labs(color = "Group")
    }

    if (!is.null(label_threshold)) {
      plot <- plot +
        ggplot2::geom_text(
          ggplot2::aes(label = labels),
          position = jitter_position,
          vjust = -0.5
        )
    }

    plot
  }
