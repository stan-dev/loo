test_that("gpdfit returns correct result", {
  set.seed(123)
  x <- rexp(100)
  gpdfit_val_old <- unlist(gpdfit(x, wip = FALSE, min_grid_pts = 80))
  expect_snapshot_value(gpdfit_val_old, style = "serialize")

  gpdfit_val_wip <- unlist(gpdfit(x, wip = TRUE, min_grid_pts = 80))
  expect_snapshot_value(gpdfit_val_wip, style = "serialize")

  gpdfit_val_wip_default_grid <- unlist(gpdfit(x, wip = TRUE))
  expect_snapshot_value(gpdfit_val_wip_default_grid, style = "serialize")
})

