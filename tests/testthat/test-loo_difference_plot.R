skip_if_not_installed("ggplot2")

log_lik <- example_loglik_matrix()
shift <- seq(-0.5, 0.5, length.out = ncol(log_lik))
log_lik_2 <- sweep(log_lik, 2, shift, FUN = "+")

loo_1 <- loo(log_lik)
loo_2 <- loo(log_lik_2)
y <- seq_len(ncol(log_lik))

test_that("plot_loo_difference returns pointwise ELPD differences", {
  p <- plot_loo_difference(y, loo_1, loo_2)

  expect_s3_class(p, "ggplot")
  expect_equal(p$data$y, y)
  expect_equal(p$data$elpd_diff, -shift)
})

test_that("plot_loo_difference sorts observations and labels by group", {
  group <- factor(
    rep(c("A", "B"), length.out = length(y)),
    levels = c("B", "A")
  )
  labels <- paste0("obs", y)
  ordering <- order(group)

  p <- plot_loo_difference(
    y,
    loo_1,
    loo_2,
    group = group,
    sort_by_group = TRUE,
    label_threshold = 0.3,
    labels = labels
  )

  expect_equal(p$data$y, seq_along(y))
  expect_equal(p$data$elpd_diff, -shift[ordering])
  expect_equal(
    as.character(p$data$group),
    as.character(group[ordering])
  )
  expect_equal(
    p$data$labels,
    ifelse(
      abs(shift[ordering]) > 0.3,
      labels[ordering],
      ""
    )
  )
})

test_that("plot_loo_difference uses observation indices as default labels", {
  p <- plot_loo_difference(
    y,
    loo_1,
    loo_2,
    label_threshold = 0.3
  )

  expect_equal(
    p$data$labels,
    ifelse(
      abs(shift) > 0.3,
      as.character(seq_along(y)),
      ""
    )
  )
})

test_that("plot_loo_difference checks observation-level arguments", {
  expect_error(
    plot_loo_difference(y[-1], loo_1, loo_2)
  )

  expect_error(
    plot_loo_difference(
      y,
      loo_1,
      loo_2,
      group = rep("A", length(y) - 1)
    )
  )

  expect_error(
    plot_loo_difference(
      y,
      loo_1,
      loo_2,
      label_threshold = -1
    )
  )

  expect_error(
    plot_loo_difference(
      y,
      loo_1,
      loo_2,
      label_threshold = 0.3,
      labels = y[-1]
    )
  )

  expect_error(
    plot_loo_difference(
      y,
      loo_1,
      loo_2,
      labels = y
    ),
    "`label_threshold` must be supplied when `labels` is supplied.",
    fixed = TRUE
  )

  expect_error(
    plot_loo_difference(
      y,
      loo_1,
      loo_2,
      sort_by_group = TRUE
    ),
    "`group` must be supplied",
    fixed = TRUE
  )
})
