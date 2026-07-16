# Data-generation pipeline for pred_measure tests and online-only vignettes.
#
# Slow: fits ~6 brms models with k-fold / LOO postprocessing (~10-20 min).
# Run from package root:
#   Rscript tests/testthat/data-for-tests/test_data_generation.R
# Vignette fits only (pkgdown CI; test RDS files stay as committed in git):
#   Rscript tests/testthat/data-for-tests/test_data_generation.R --vignettes-only
#
# Outputs:
#   tests/testthat/data-for-tests/test_data_*.Rds  (unless --vignettes-only)
#   vignettes/articles-online-only/data-for-vignettes/fit_*.Rds
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
    mupred_loo <- loo::E_loo(ypred, psis_object, type = "mean")$value
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
      mupred_loo = mupred_loo,
      mupred_kfold = mupred_kfold
    )
  }
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
      ylp_test = brms::log_lik(fit_sleep_train, newdata = test_data),
      ylp_train = brms::log_lik(fit_sleep_train)
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

generate_test_data <- function(vignettes_only = FALSE) {
  if (vignettes_only) {
    message(
      "Generating vignette brms fits only (test fixtures unchanged). ",
      "This may take 10-20 minutes."
    )
  } else {
    message(
      "Generating test fixtures and vignette brms fits. ",
      "This may take 10-20 minutes."
    )
  }
  t0 <- proc.time()

  full_roaches <- get_roaches_res()
  full_binary <- get_binary_res()
  full_penguins <- get_penguins_res()
  full_binomial <- get_binomial_res()
  full_sleep <- get_sleep_res()
  full_sleep_test <- get_sleep_test_train_res()

  if (!vignettes_only) {
    test_path <- "tests/testthat/data-for-tests/"
    saveRDS(full_roaches$res, paste0(test_path, "test_data_roaches.Rds"))
    saveRDS(full_binary$res, paste0(test_path, "test_data_binary.Rds"))
    saveRDS(full_penguins$res, paste0(test_path, "test_data_penguins.Rds"))
    saveRDS(full_binomial$res, paste0(test_path, "test_data_binomial.Rds"))
    saveRDS(full_sleep$res, paste0(test_path, "test_data_sleep.Rds"))
    saveRDS(full_sleep_test$res, paste0(test_path, "test_data_sleep_cv.Rds"))
    message("Saved test fixtures to ", test_path)
  }

  vignette_path <- "vignettes/articles-online-only/data-for-vignettes/"
  dir.create(vignette_path, recursive = TRUE, showWarnings = FALSE)
  saveRDS(full_roaches$fit, paste0(vignette_path, "fit_roaches.Rds"))
  saveRDS(full_binary$fit, paste0(vignette_path, "fit_binary.Rds"))
  saveRDS(full_penguins$fit, paste0(vignette_path, "fit_penguins.Rds"))
  saveRDS(full_binomial$fit, paste0(vignette_path, "fit_binomial.Rds"))
  saveRDS(full_sleep$fit, paste0(vignette_path, "fit_sleep.Rds"))
  message("Saved vignette fits to ", vignette_path)

  elapsed_min <- round((proc.time() - t0)[3] / 60, 1)
  message("Data generation finished in ", elapsed_min, " minutes.")
  invisible(NULL)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0L) {
  if (!("--vignettes-only" %in% args)) {
    stop(
      "Unknown argument(s): ", paste(args, collapse = ", "),
      ". Usage: [--vignettes-only]",
      call. = FALSE
    )
  }
  generate_test_data(vignettes_only = TRUE)
} else {
  generate_test_data(vignettes_only = FALSE)
}