options(mc.cores = 1)

LLarr <- example_loglik_array()
LLmat <- example_loglik_matrix()
chain_id <- rep(1:2, each = dim(LLarr)[1])
r_eff_mat <- relative_eff(exp(LLmat), chain_id = chain_id)

# log-likelihood of other draws that correspond one to one to the draws of
# LLmat (for example, projected draws)
set.seed(123)
LLmat_other <- LLmat + matrix(rnorm(length(LLmat), sd = 0.1), nrow(LLmat))

test_that("loo with psis_object equals loo without it for the same ratios", {
  loo_a <- loo(LLmat, r_eff = r_eff_mat, save_psis = TRUE)
  loo_b <- loo(LLmat, save_psis = TRUE,
               psis_object = psis(-LLmat, r_eff = r_eff_mat))
  expect_identical(loo_b, loo_a)

  r_eff_arr <- relative_eff(exp(LLarr))
  loo_c <- loo(LLarr, r_eff = r_eff_arr, save_psis = TRUE)
  loo_d <- loo(LLarr, save_psis = TRUE,
               psis_object = psis(-LLarr, r_eff = r_eff_arr))
  expect_identical(loo_d, loo_c)
  loo_e <- loo(LLmat, save_psis = TRUE, psis_object = loo_d$psis_object)
  expect_identical(loo_e$pointwise, loo_c$pointwise)
})

test_that("loo with psis_object uses the weights of psis_object", {
  psis_obj <- psis(-LLmat, r_eff = r_eff_mat)
  loo_x <- loo(LLmat_other, psis_object = psis_obj)
  lw <- weights(psis_obj)
  elpd <- matrixStats::colLogSumExps(LLmat_other + lw)
  expect_equal(loo_x$pointwise[, "elpd_loo"], elpd)
  expect_equal(
    loo_x$pointwise[, "elpd_loo"],
    log(E_loo(exp(LLmat_other), psis_obj, type = "mean",
              log_ratios = -LLmat)$value)
  )
  lpd <- matrixStats::colLogSumExps(LLmat_other) - log(nrow(LLmat_other))
  expect_equal(loo_x$pointwise[, "p_loo"], lpd - elpd)
  expect_equal(
    loo_x$pointwise[, "mcse_elpd_loo"],
    mcse_elpd(LLmat_other, lw, E_elpd = elpd, r_eff = r_eff_mat)
  )
  expect_identical(loo_x$diagnostics, psis_obj$diagnostics)
  expect_identical(loo_x$pointwise[, "influence_pareto_k"],
                   psis_obj$diagnostics$pareto_k)
  expect_identical(dim(loo_x), dim(psis_obj))
  expect_null(loo_x$psis_object)
  expect_s3_class(loo_x, c("psis_loo", "importance_sampling_loo", "loo"),
                  exact = TRUE)
})

test_that("r_eff with psis_object sets n_eff and the MCSE, not the smoothing", {
  psis_1 <- psis(-LLmat, r_eff = 1)
  loo_x <- loo(LLmat_other, psis_object = psis_1, r_eff = r_eff_mat,
               save_psis = TRUE)
  # the smoothing of psis_1 is kept
  expect_identical(loo_x$diagnostics$pareto_k, psis_1$diagnostics$pareto_k)
  expect_identical(loo_x$psis_object$log_weights, psis_1$log_weights)
  expect_identical(attr(loo_x$psis_object, "tail_len"),
                   attr(psis_1, "tail_len"))
  # n_eff and the MCSE use the given r_eff
  w <- weights(psis_1, log = FALSE)
  expect_equal(loo_x$diagnostics$n_eff, r_eff_mat / colSums(w^2))
  expect_identical(loo_x$diagnostics$r_eff, r_eff_mat)
  expect_identical(relative_eff(loo_x$psis_object), r_eff_mat)
  expect_identical(loo_x$psis_object$diagnostics, loo_x$diagnostics)
  lw <- weights(psis_1)
  elpd <- matrixStats::colLogSumExps(LLmat_other + lw)
  expect_equal(
    loo_x$pointwise[, "mcse_elpd_loo"],
    mcse_elpd(LLmat_other, lw, E_elpd = elpd, r_eff = r_eff_mat)
  )
  # without r_eff, the r_eff of psis_object is used
  loo_y <- loo(LLmat_other, psis_object = psis_1)
  expect_identical(loo_y$diagnostics, psis_1$diagnostics)
  # a scalar r_eff is recycled
  loo_z <- loo(LLmat_other, psis_object = psis_1, r_eff = 0.5)
  expect_identical(loo_z$diagnostics$r_eff, rep(0.5, ncol(LLmat)))
})

test_that("loo with psis_object keeps the importance sampling method", {
  # few draws, the case in which a caller uses sis() instead of psis()
  S <- 20
  sis_obj <- sis(-LLmat[1:S, ])
  loo_s <- loo(LLmat_other[1:S, ], psis_object = sis_obj, save_psis = TRUE)
  expect_s3_class(loo_s, c("sis_loo", "importance_sampling_loo", "loo"),
                  exact = TRUE)
  expect_named(loo_s, names(loo(LLmat[1:S, ], is_method = "sis")))
  expect_identical(loo_s$sis_object, sis_obj)
  tis_obj <- tis(-LLmat)
  loo_t <- loo(LLmat_other, psis_object = tis_obj, is_method = "tis")
  expect_s3_class(loo_t, "tis_loo")
})

test_that("loo with psis_object works with loo_compare and print", {
  psis_obj <- psis(-LLmat, r_eff = r_eff_mat)
  loo_a <- loo(LLmat, r_eff = r_eff_mat)
  loo_x <- loo(LLmat_other, psis_object = psis_obj)
  expect_no_error(comp <- loo_compare(loo_a, loo_x))
  expect_equal(
    abs(comp[2, "elpd_diff"]),
    abs(sum(loo_x$pointwise[, "elpd_loo"] - loo_a$pointwise[, "elpd_loo"]))
  )
  expect_output(print(loo_x), "elpd_loo")
})

test_that("loo with psis_object gives errors for wrong input", {
  psis_obj <- psis(-LLmat, r_eff = r_eff_mat)
  expect_error(loo(LLmat, psis_object = weights(psis_obj)),
               "must be an object of class")
  expect_error(loo(LLmat[, -1], psis_object = psis_obj),
               "same number of draws and observations")
  expect_error(loo(LLmat[-1, ], psis_object = psis_obj),
               "same number of draws and observations")
  expect_error(loo(LLmat, psis_object = psis_obj, is_method = "sis"),
               "'is_method' is 'sis'")
  expect_error(loo(LLmat, psis_object = psis_obj, r_eff = c(1, 1)),
               "'r_eff' must have one value or one value per observation")
  psis_arr <- psis(-LLarr, r_eff = relative_eff(exp(LLarr)))
  expect_error(loo(LLarr[, , -1], psis_object = psis_arr),
               "same number of draws and observations")
  expect_error(loo(LLarr[-1, , ], psis_object = psis_arr),
               "same number of draws and observations")
  mixed <- psis_obj
  attr(mixed, "method") <- rep(c("psis", "sis"), length.out = ncol(LLmat))
  expect_error(loo(LLmat, psis_object = mixed), "one importance sampling")
  llfun <- function(data_i, draws, ...) {
    dnorm(data_i$y, draws, log = TRUE)
  }
  expect_error(
    loo(llfun, data = data.frame(y = 1:3), draws = rnorm(10),
        psis_object = psis_obj),
    "only by the array and matrix methods"
  )
})

test_that("loo with psis_object handles r_eff = NULL and NA", {
  psis_1 <- psis(-LLmat, r_eff = r_eff_mat)
  loo_null <- loo(LLmat_other, psis_object = psis_1, r_eff = NULL)
  expect_identical(loo_null$diagnostics, psis_1$diagnostics)
  r_na <- r_eff_mat
  r_na[1:2] <- NA
  expect_message(
    loo_na <- loo(LLmat_other, psis_object = psis_1, r_eff = r_na),
    "Replacing NAs in `r_eff` with 1s"
  )
  expect_identical(loo_na$diagnostics$r_eff, c(1, 1, r_eff_mat[-(1:2)]))
  loo_na1 <- loo(LLmat_other, psis_object = psis_1, r_eff = NA)
  expect_identical(loo_na1$diagnostics$r_eff, rep(1, ncol(LLmat)))
  # an object without r_eff uses r_eff = 1
  no_r_eff <- psis_1
  attr(no_r_eff, "r_eff") <- NULL
  loo_1 <- loo(LLmat_other, psis_object = no_r_eff)
  expect_identical(loo_1$diagnostics$r_eff, rep(1, ncol(LLmat)))
  expect_false(anyNA(loo_1$pointwise[, "mcse_elpd_loo"]))
})

test_that("loo with psis_object works through a wrapper method", {
  wrapper <- function(x, psis_object, r_eff = NULL,
                      is_method = c("psis", "tis", "sis")) {
    loo(x, r_eff = r_eff, is_method = is_method, psis_object = psis_object)
  }
  sis_obj <- sis(-LLmat[1:20, ])
  loo_w <- wrapper(LLmat_other[1:20, ], psis_object = sis_obj)
  expect_s3_class(loo_w, "sis_loo")
  expect_identical(loo_w$diagnostics, sis_obj$diagnostics)
  psis_obj <- psis(-LLmat, r_eff = r_eff_mat)
  expect_no_error(loo(LLmat, psis_object = psis_obj, is_method = "psis"))
})

test_that("loo.array with psis_object uses the chain-stacked draw order", {
  r_eff_arr <- relative_eff(exp(LLarr))
  psis_mat <- psis(-llarray_to_matrix(LLarr), r_eff = r_eff_arr)
  expect_identical(
    loo(LLarr, psis_object = psis_mat)$pointwise,
    loo(LLarr, r_eff = r_eff_arr)$pointwise
  )
})
