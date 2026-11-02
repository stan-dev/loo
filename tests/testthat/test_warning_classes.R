test_that(".warn() signals one message with the loo classes", {
  w <- capture_warning(.warn("a", c("b", "c"), 1, class = "loo_warning_x"))
  expect_equal(conditionMessage(w), "abc1")
  expect_equal(class(w)[1:2], c("loo_warning_x", "loo_warning"))
})

test_that(".warn_diagnostic() adds the diagnostic class", {
  w <- capture_warning(.warn_diagnostic("a", class = "loo_warning_x"))
  expect_equal(
    class(w)[1:3],
    c("loo_warning_x", "loo_warning_diagnostic", "loo_warning")
  )
})

test_that("warning helpers signal their classes", {
  expect_warning(throw_pareto_warnings(1, 0.7), class = "loo_warning_pareto_k_high")
  expect_warning(throw_tail_length_warnings(1), class = "loo_warning_pareto_tail_not_fitted")
  expect_warning(throw_pwaic_warnings(1), class = "loo_warning_p_waic_high")
  expect_warning(throw_moment_match_max_iters_warning(), class = "loo_warning_moment_match_max_iters")
  expect_warning(throw_moment_match_error_warning(1, 2, "e"), class = "loo_warning_moment_match_failed")
  expect_warning(throw_large_kf_warning(1, 0.7), class = "loo_warning_large_kf")
})
