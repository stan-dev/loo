# Data-generation pipeline for the pred_measure test fixtures.
#
# Slow: fits ~6 brms models with k-fold / LOO postprocessing (~10-20 min).
# Run from package root:
#   Rscript tests/testthat/data-for-tests/test_data_generation.R
#
# Outputs:
#   tests/testthat/data-for-tests/test_data_*.Rds
library(rstanarm)
# brms features used here require GitHub master (not yet on CRAN)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
if (!requireNamespace("brms", quietly = TRUE)) {
  remotes::install_github("paul-buerkner/brms", ref = "master", upgrade = "never")
}
suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(loo)
})

SEED <- 42

postprocess_res <- function(model, fit, chains = 2, draws = 200) {
  ypred <- brms::posterior_predict(fit)
  mupred <- brms::posterior_epred(fit)
  ylp <- log_lik(fit)
  log_ratios <- -1 * ylp
  r_eff <- relative_eff(exp(-log_ratios), chain_id = rep(1:chains, each = draws))
  psis_object <- psis(log_ratios, r_eff = r_eff, cores = 2)
  kfold2 <- brms::kfold(fit, save_fits = FALSE)
  if (model %in% c("roaches", "binary", "binomial", "sleep")) {
    kfold <- brms::kfold(fit, save_fits = TRUE)
    mupred_kfold <- brms::kfold_predict(kfold, method = "fitted")$yrep
    ypred_kfold <- brms::kfold_predict(kfold, method = "predict")$yrep
    loo <- brms::loo(fit, save_psis = TRUE)
    predperf <- insample_pred_measure(y = fit$data$y, mupred = mupred, 
      measure = "r2", ylp = ylp)
  }
  
  if (model == "roaches") {
    list(
      y = fit$data$y,
      ypred = ypred,
      mupred = mupred,
      ylp = ylp,
      log_weights = psis_object$log_weights,
      kfold = kfold2,
      loo = loo,
      predperf = predperf
    )
  } else if (model == "binomial") {
    list(
      y = fit$data$y,
      ypred = ypred,
      log_weights = psis_object$log_weights,
      ypred_kfold = ypred_kfold
    )
  } else if (model == "binary") {
    list(
      y = fit$data$y,
      ypred = ypred,
      log_weights = psis_object$log_weights,
      ypred_kfold = ypred_kfold
    )
  } else if (model == "categorical") {
    list(
      y = fit$data$y,
      mupred = mupred,
      log_weights = psis_object$log_weights
    )
  } else if (model == "sleep") {
    list(
      y = fit$data$y,
      ypred = ypred,
      log_weights = psis_object$log_weights,
      ypred_kfold = ypred_kfold,
      mupred = mupred,
      mupred_kfold = mupred_kfold
    )
  }
}

# ---- fixture shrinking ------------------------------------------------------
# These fixtures ship in the source tarball, which CRAN limits to 5 MB. Keep
# only a subset of the observations. The draws stay at 400, so the Pareto k
# threshold ps_khat_threshold(400) does not move.
N_KEEP <- c(
  roaches = 53, roaches_compare = 110, categorical = 67, sleep = 29,
  sleep_test = 20
)

# Estimate and SE of a summed pointwise column, in loo's convention.
.estimates_from_pointwise <- function(pointwise, cols) {
  pw <- pointwise[, cols, drop = FALSE]
  cbind(Estimate = colSums(pw), SE = sqrt(nrow(pw) * apply(pw, 2, var)))
}

# Index of the observations to keep. `strata` allocates proportionally, so
# every level of a categorical outcome survives.
.keep_index <- function(n, n_keep, strata = NULL) {
  set.seed(SEED)
  if (is.null(strata)) {
    return(sort(sample.int(n, n_keep)))
  }
  strata <- as.factor(strata)
  per <- pmax(1L, round(n_keep * as.vector(table(strata)) / n))
  idx <- unlist(Map(
    function(lev, k) sample(which(strata == lev), k),
    levels(strata), per
  ))
  sort(as.integer(idx))
}

# PSIS smoothing is per observation, so subsetting columns is exact. Only the
# aggregates need recomputing.
.shrink_psis_loo <- function(x, keep) {
  x$pointwise <- x$pointwise[keep, , drop = FALSE]
  x$diagnostics <- lapply(x$diagnostics, function(d) d[keep])

  po <- x$psis_object
  po$log_weights <- po$log_weights[, keep, drop = FALSE]
  po$diagnostics <- lapply(po$diagnostics, function(d) d[keep])
  for (a in c("norm_const_log", "tail_len", "r_eff")) {
    attr(po, a) <- attr(po, a)[keep]
  }
  attr(po, "dims") <- c(attr(po, "dims")[1], length(keep))
  x$psis_object <- po

  est <- .estimates_from_pointwise(x$pointwise, c("elpd_loo", "p_loo", "looic"))
  x$estimates <- est
  x$elpd_loo <- est["elpd_loo", "Estimate"]
  x$p_loo <- est["p_loo", "Estimate"]
  x$looic <- est["looic", "Estimate"]
  x$se_elpd_loo <- est["elpd_loo", "SE"]
  x$se_p_loo <- est["p_loo", "SE"]
  x$se_looic <- est["looic", "SE"]
  attr(x, "dims") <- c(attr(x, "dims")[1], length(keep))
  x
}

.shrink_kfold <- function(x, keep) {
  x$pointwise <- x$pointwise[keep, , drop = FALSE]
  x$estimates <- .estimates_from_pointwise(x$pointwise, colnames(x$pointwise))
  attr(x, "folds") <- attr(x, "folds")[keep]
  x
}

shrink_res <- function(model, res) {
  # binary and binomial hold 50 observations and 86 KB in total. Leave them.
  if (model %in% c("binary", "binomial")) {
    return(res)
  }
  if (model == "sleep_test") {
    keep_test <- .keep_index(length(res$y_test), N_KEEP[["sleep_test"]])
    res$y_test <- res$y_test[keep_test]
    for (nm in c("ypred_test", "mupred_test", "ylp_test")) {
      res[[nm]] <- res[[nm]][, keep_test, drop = FALSE]
    }
    return(res)
  }

  keep <- if (model == "categorical") {
    .keep_index(length(res$y), N_KEEP[["categorical"]], strata = res$y)
  } else {
    .keep_index(length(res$y), N_KEEP[[model]])
  }

  res$y <- res$y[keep]
  for (nm in c("ypred", "mupred", "ylp", "log_weights",
               "ypred_kfold", "mupred_kfold")) {
    if (!is.null(res[[nm]]) && length(dim(res[[nm]])) == 2L) {
      res[[nm]] <- res[[nm]][, keep, drop = FALSE]
    }
  }
  # A categorical mupred is draws x observations x categories.
  if (!is.null(res$mupred) && length(dim(res$mupred)) == 3L) {
    res$mupred <- res$mupred[, keep, , drop = FALSE]
  }
  if (!is.null(res$loo)) {
    res$loo <- .shrink_psis_loo(res$loo, keep)
  }
  if (!is.null(res$kfold)) {
    res$kfold <- .shrink_kfold(res$kfold, keep)
  }
  if (!is.null(res$predperf)) {
    res$predperf <- insample_pred_measure(
      y = res$y, mupred = res$mupred, measure = "r2", ylp = res$ylp
    )
  }
  res
}

# The model-comparison fixture holds four `psis_loo` objects and four sets of
# draws. `.keep_index()` reseeds, so this keeps the same 53 observations as
# `test_data_roaches.Rds`.
shrink_roaches_compare <- function(res) {
  keep <- .keep_index(length(res$y), N_KEEP[["roaches_compare"]])
  res$y <- res$y[keep]
  for (nm in grep("^(ypred|mupred|ylp)(_m[0-9]+)?$", names(res), value = TRUE)) {
    res[[nm]] <- res[[nm]][, keep, drop = FALSE]
  }
  for (nm in grep("^loo_p(_m[0-9]+)?$", names(res), value = TRUE)) {
    res[[nm]] <- .shrink_psis_loo(res[[nm]], keep)
  }
  res
}


get_binary_res <- function() {
  set.seed(SEED)
  df_binary <- data.frame(y = rbinom(50, 1, 0.3))

  fit_binary <- brms::brm(formula = "y ~ 1",
    data = df_binary,
    family = bernoulli,
    chains = 2,
    iter = 400,
    seed = SEED,
    refresh = 0
  )
  list(
    fit = fit_binary,
    res = postprocess_res("binary", fit_binary)
  )
}

get_roaches_res <- function() {
  data(roaches, package = "rstanarm")
  roaches$sqrt_roach1 <- sqrt(roaches$roach1)
  
  fit_roaches <- brm(
    y ~ sqrt_roach1 + treatment + senior + offset(log(exposure2)),
    data = roaches,
    family = poisson,
    prior = prior(normal(0, 1), class = b),
    chains = 2,
    iter = 400,
    refresh = 0,
    seed = SEED
  )
  list(
    fit = fit_roaches,
    res = postprocess_res("roaches", fit_roaches)
  )
}

get_roaches_compare_res <- function() {
  data(roaches, package = "rstanarm")
  roaches$sqrt_roach1 <- sqrt(roaches$roach1)
  
  fit_p <- brm(
    y ~ sqrt_roach1 + treatment + senior + offset(log(exposure2)),
    data = roaches,
    family = poisson,
    prior = prior(normal(0, 1), class = b),
    chains = 2,
    iter = 400,
    refresh = 0,
    seed = SEED
  )

  fit_p <- add_criterion(
    fit_p,
    criterion = "loo",
    moment_match = TRUE,
    save_psis = TRUE,
    overwrite = TRUE
  )

  fit_p_m1 <- update(fit_p, formula = y ~ treatment + senior) |>
    add_criterion(criterion = "loo", moment_match = TRUE, save_psis = TRUE)
  fit_p_m2 <- update(fit_p, formula = y ~ sqrt_roach1 + senior)  |>
    add_criterion(criterion = "loo", moment_match = TRUE, save_psis = TRUE)
  fit_p_m3 <- update(fit_p, formula = y ~ sqrt_roach1 + treatment) |>
    add_criterion(criterion = "loo", moment_match = TRUE, save_psis = TRUE)

  # `ypred` (posterior predictive draws) is needed by the sampling-based scores
  # such as `rps`/`srps`; `mupred` (posterior_epred) is not enough for those.
  set.seed(SEED)
  return(list(
    y = fit_p$data$y,
    loo_p = fit_p$criteria$loo,
    ypred = brms::posterior_predict(fit_p),
    mupred = brms::posterior_epred(fit_p),
    ylp = brms::log_lik(fit_p),
    loo_p_m1 = fit_p_m1$criteria$loo,
    ypred_m1 = brms::posterior_predict(fit_p_m1),
    mupred_m1 = brms::posterior_epred(fit_p_m1),
    ylp_m1 = brms::log_lik(fit_p_m1),
    loo_p_m2 = fit_p_m2$criteria$loo,
    ypred_m2 = brms::posterior_predict(fit_p_m2),
    mupred_m2 = brms::posterior_epred(fit_p_m2),
    ylp_m2 = brms::log_lik(fit_p_m2),
    loo_p_m3 = fit_p_m3$criteria$loo,
    ypred_m3 = brms::posterior_predict(fit_p_m3),
    mupred_m3 = brms::posterior_epred(fit_p_m3),
    ylp_m3 = brms::log_lik(fit_p_m3)
  ))
}

get_sleep_test_train_res <- function() {
  # specifically for testing test_pred_measure
  data("sleepstudy", package = "lme4")
  conditions <- brms::make_conditions(sleepstudy, "Subject", incl_vars = FALSE)
  sleepstudy <- sleepstudy |>
    dplyr::filter(Days >= 2) |>
    dplyr::mutate(
      Days = Days - 2,
      y = Reaction
    )
  
  test_subjects <- sample(unique(sleepstudy$Subject), size = 5)
  train_data <- sleepstudy |>
    dplyr::filter(!Subject %in% test_subjects)
  test_data <- sleepstudy |>
    dplyr::filter(Subject %in% test_subjects)

  prior_lin_base <- brms::prior(normal(200, 100), class = b, coef = "Intercept") +
    brms::prior(normal(0, 20), class = b, coef = "Days") +
    brms::prior(exponential(0.02), class = sigma)

  fit_sleep_train <- brm(
    y ~ 0 + Intercept + Days, 
    data = train_data,
    family = gaussian(),
    prior = prior_lin_base,
    chains = 2, 
    iter = 400,
    seed = SEED
  )
  
  list(
    fit = fit_sleep_train,
    res = list(
      y_test = test_data$y,
      ypred_test = brms::posterior_predict(fit_sleep_train, newdata = test_data),
      mupred_test = brms::posterior_epred(fit_sleep_train, newdata = test_data),
      ylp_test = brms::log_lik(fit_sleep_train, newdata = test_data)
    )
  )
}

get_sleep_res <- function() {
  data("sleepstudy", package = "lme4")
  conditions <- brms::make_conditions(sleepstudy, "Subject", incl_vars = FALSE)
  sleepstudy <- sleepstudy |>
    dplyr::filter(Days >= 2) |>
    dplyr::mutate(Days = Days - 2, y = Reaction)

  prior_lin_base <- brms::prior(normal(200, 100), class = b, coef = "Intercept") +
    brms::prior(normal(0, 20), class = b, coef = "Days") +
    brms::prior(exponential(0.02), class = sigma)

  fit_sleepstudy <- brms::brm(
    y ~ 0 + Intercept + Days, 
    data = sleepstudy,
    family = gaussian(),
    prior = prior_lin_base,
    chains = 2,
    iter = 400,
    refresh = 0,
    seed = SEED
  )
  list(
    fit = fit_sleepstudy,
    res = postprocess_res("sleep", fit_sleepstudy)
  )
}

get_penguins_res <- function() {
  data("penguins", package = "palmerpenguins")
  penguins <- subset(penguins, complete.cases(penguins))
  penguins$y <- penguins$species
  
  fit <- brm(
    y ~ bill_length_mm + bill_depth_mm,
    data = penguins,
    family = categorical(),
    chains = 2,
    iter = 400,
    cores = 2,
    seed = SEED
  )
  list(
    fit = fit,
    res = postprocess_res("categorical", fit)
  )
}

get_binomial_res <- function() {
  set.seed(SEED)
  df <- data.frame(
    y = rbinom(50, 10, 0.3), 
    n = 10
  )

  fit <- brms::brm(formula = "y | trials(n) ~ 1",
    data = df,
    family = binomial,
    chains = 2,
    iter = 400,
    seed = SEED,
    refresh = 0
  )
  list(
    fit = fit,
    res = postprocess_res("binomial", fit)
  )
}

generate_test_data <- function() {
  message("Generating test fixtures. This may take 10-20 minutes.")
  t0 <- proc.time()

  full_roaches <- get_roaches_res()
  full_binary <- get_binary_res()
  full_penguins <- get_penguins_res()
  full_binomial <- get_binomial_res()
  full_sleep <- get_sleep_res()
  full_sleep_test <- get_sleep_test_train_res()
  full_roaches_compare <- get_roaches_compare_res()

  test_path <- "tests/testthat/data-for-tests/"
  saveRDS(shrink_res("roaches", full_roaches$res), paste0(test_path, "test_data_roaches.Rds"))
  saveRDS(shrink_roaches_compare(full_roaches_compare), paste0(test_path, "test_data_roaches_compare.Rds"))
  saveRDS(shrink_res("binary", full_binary$res), paste0(test_path, "test_data_binary.Rds"))
  saveRDS(shrink_res("categorical", full_penguins$res), paste0(test_path, "test_data_penguins.Rds"))
  saveRDS(shrink_res("binomial", full_binomial$res), paste0(test_path, "test_data_binomial.Rds"))
  saveRDS(shrink_res("sleep", full_sleep$res), paste0(test_path, "test_data_sleep.Rds"))
  saveRDS(shrink_res("sleep_test", full_sleep_test$res), paste0(test_path, "test_data_sleep_cv.Rds"))
  message("Saved test fixtures to ", test_path)

  elapsed_min <- round((proc.time() - t0)[3] / 60, 1)
  message("Data generation finished in ", elapsed_min, " minutes.")
  invisible(NULL)
}

generate_test_data()
