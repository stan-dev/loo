#' Compare `psis_loo_ss` objects
#' @noRd
#' @param x A list with `psis_loo` objects.
#' @param ... Currently ignored.
#' @return A `compare.loo_ss` object.
#' @author Mans Magnusson
loo_compare.psis_loo_ss_list <- function(x, ...) {

  checkmate::assert_list(x, any.missing = FALSE, min.len = 1)
  for(i in seq_along(x)){
    if (!inherits(x[[i]], "psis_loo_ss")) x[[i]] <- as.psis_loo_ss.psis_loo(x[[i]])
  }

  loo_compare_checks.psis_loo_ss_list(x)

  comp <- loo_compare_matrix.psis_loo_ss_list(x)
  ord <- loo_compare_order(x)
  names(x) <- rownames(comp)[ord]

  rnms <- rownames(comp)
  elpd_diff_mat <- matrix(0, nrow = nrow(comp), ncol = 3,
                          dimnames = list(rnms, c("elpd_diff", "se_diff", "subsampling_se_diff")))
  for(i in 2:length(ord)){
    elpd_diff_mat[i,] <- loo_compare_ss(ref_loo = x[ord[1]], compare_loo = x[ord[i]])
  }
  comp <- cbind(elpd_diff_mat, comp)
  rownames(comp) <- rnms

  class(comp) <- c("compare.loo_ss", "compare.loo", class(comp))
  return(comp)
}

#' Compare a reference loo object with a comaprison loo object
#' @noRd
#' @param ref_loo A named list with a `psis_loo_ss` object.
#' @param compare_loo A named list with a  `psis_loo_ss` object.
#' @return A 1 by 3 elpd_diff estimation.
loo_compare_ss <- function(ref_loo, compare_loo){
  checkmate::assert_list(ref_loo, names = "named")
  checkmate::assert_list(compare_loo, names = "named")
  checkmate::assert_class(ref_loo[[1]], "psis_loo_ss")
  checkmate::assert_class(compare_loo[[1]], "psis_loo_ss")


  ref_idx <- obs_idx(ref_loo[[1]])
  compare_idx <- obs_idx(compare_loo[[1]])
  intersect_idx <- base::intersect(ref_idx, compare_idx)
  ref_subset_of_compare <- base::setequal(intersect_idx, ref_idx)
  compare_subset_of_ref <- base::setequal(intersect_idx, compare_idx)

  # Using HH estimation
  if (ref_loo[[1]]$loo_subsampling$estimator == "hh_pps" | compare_loo[[1]]$loo_subsampling$estimator == "hh_pps"){
    warning("Hansen-Hurwitz estimator used. Naive diff SE is used.", call. = FALSE)
    return(loo_compare_ss_naive(ref_loo, compare_loo))
  }

  # Same observations in both
  if (compare_subset_of_ref & ref_subset_of_compare){
    return(loo_compare_ss_diff(ref_loo, compare_loo))
  }

  # Use subset
  if (compare_subset_of_ref | ref_subset_of_compare){
    if (compare_subset_of_ref) ref_loo[[1]] <- update(object = ref_loo[[1]], observations = compare_loo[[1]])
    if (ref_subset_of_compare) compare_loo[[1]] <- update(compare_loo[[1]], observations = ref_loo[[1]])
    message("Estimated elpd_diff using observations included in loo calculations for all models.")
    return(loo_compare_ss_diff(ref_loo, compare_loo))
  }

  # If different samples
  if (!compare_subset_of_ref & !ref_subset_of_compare){
    warning("Different subsamples in '", names(ref_loo), "' and '", names(compare_loo),
            "'. Naive diff SE is used.", call. = FALSE)
    return(loo_compare_ss_naive(ref_loo, compare_loo))
  }
}

#' Compute a naive diff SE
#' @noRd
#' @inheritParams loo_compare_ss
#' @return a 1 by 3 elpd_diff estimation
loo_compare_ss_naive <- function(ref_loo, compare_loo){
  checkmate::assert_list(ref_loo, names = "named")
  checkmate::assert_list(compare_loo, names = "named")
  checkmate::assert_class(ref_loo[[1]], "psis_loo_ss")
  checkmate::assert_class(compare_loo[[1]], "psis_loo_ss")

  elpd_loo_diff <- ref_loo[[1]]$estimates["elpd_loo","Estimate"] - compare_loo[[1]]$estimates["elpd_loo","Estimate"]
  elpd_loo_diff_se <- sqrt(
    (ref_loo[[1]]$estimates["elpd_loo","SE"])^2 +
    (compare_loo[[1]]$estimates["elpd_loo","SE"])^2)
  elpd_loo_diff_subsampling_se <- sqrt(
      (ref_loo[[1]]$estimates["elpd_loo","subsampling SE"])^2 +
      (compare_loo[[1]]$estimates["elpd_loo","subsampling SE"])^2)

  c(elpd_loo_diff, elpd_loo_diff_se, elpd_loo_diff_subsampling_se)
}

#' Compare a effective diff SE
#' @noRd
#' @inheritParams loo_compare_ss
#' @return a 1 by 3 elpd_diff estimation
loo_compare_ss_diff <- function(ref_loo, compare_loo){
  checkmate::assert_list(ref_loo, names = "named")
  checkmate::assert_list(compare_loo, names = "named")
  checkmate::assert_class(ref_loo[[1]], "psis_loo_ss")
  checkmate::assert_class(compare_loo[[1]], "psis_loo_ss")
  checkmate::assert_true(identical(obs_idx(ref_loo[[1]]), obs_idx(compare_loo[[1]])))

  # Assert not none as loo approximation
  checkmate::assert_true(ref_loo[[1]]$loo_subsampling$loo_approximation != "none")
  checkmate::assert_true(compare_loo[[1]]$loo_subsampling$loo_approximation != "none")

  diff_approx <- ref_loo[[1]]$loo_subsampling$elpd_loo_approx - compare_loo[[1]]$loo_subsampling$elpd_loo_approx
  diff_sample <- ref_loo[[1]]$pointwise[,"elpd_loo"] - compare_loo[[1]]$pointwise[,"elpd_loo"]
  est <- srs_diff_est(diff_approx, y = diff_sample, y_idx = ref_loo[[1]]$pointwise[,"idx"])

  elpd_loo_diff <- est$y_hat
  elpd_loo_diff_se <- sqrt(est$hat_v_y)
  elpd_loo_diff_subsampling_se <- sqrt(est$v_y_hat)

  c(elpd_loo_diff, elpd_loo_diff_se, elpd_loo_diff_subsampling_se)
}


#' Check list of `psis_loo` objects
#' @details Similar to `loo_compare_checks()` but checks dim size rather than
#' pointwise dim since different pointwise sizes of `psis_loo_ss` will work.
#' Can probably be removed by refactoring `loo_compare_checks()`.
#' @noRd
#' @inheritParams loo_compare_ss
#' @return A 1 by 3 elpd_diff estimation.
loo_compare_checks.psis_loo_ss_list <- function(loos) {
  ## errors
  if (length(loos) <= 1L) {
    stop("'loo_compare' requires at least two models.", call.=FALSE)
  }
  if (!all(sapply(loos, is.loo))) {
    stop("All inputs should have class 'loo'.", call.=FALSE)
  }

  Ns <- sapply(loos, function(x) x$loo_subsampling$data_dim[1])
  if (!all(Ns == Ns[1L])) {
    stop("Not all models have the same number of data points.", call.=FALSE)
  }

  ## warnings

  yhash <- lapply(loos, attr, which = "yhash")
  yhash_ok <- sapply(yhash, function(x) { # ok only if all yhash are same (all NULL is ok)
    isTRUE(all.equal(x, yhash[[1]]))
  })
  if (!all(yhash_ok)) {
    warning("Not all models have the same y variable. ('yhash' attributes do not match)",
            call. = FALSE)
  }

  if (all(sapply(loos, is.kfold))) {
    Ks <- unlist(lapply(loos, attr, which = "K"))
    if (!all(Ks == Ks[1])) {
      warning("Not all kfold objects have the same K value. ",
              "For a more accurate comparison use the same number of folds. ",
              call. = FALSE)
    }
  } else if (any(sapply(loos, is.kfold)) && any(sapply(loos, is.psis_loo))) {
    warning("Comparing LOO-CV to K-fold-CV. ",
            "For a more accurate comparison use the same number of folds ",
            "or loo for all models compared.",
            call. = FALSE)
  }
}

#' @rdname loo_compare
#' @export
print.compare.loo_ss <- function(x, ..., digits = 1, simplify = TRUE) {
  xcopy <- x
  if (inherits(xcopy, "old_compare.loo")) {
    if (NCOL(xcopy) >= 2 && simplify) {
      patts <- "^elpd_|^se_diff|^p_|^waic$|^looic$"
      xcopy <- xcopy[, grepl(patts, colnames(xcopy))]
    }
  } else if (NCOL(xcopy) >= 2 && simplify) {
    xcopy <- xcopy[, c("elpd_diff", "se_diff", "subsampling_se_diff")]
  }
  print(.fr(xcopy, digits), quote = FALSE)
  invisible(x)
}


#' Compute comparison matrix for `psis_loo_ss` objects
#' @noRd
#' @keywords internal
#' @param loos List of `psis_loo_ss` objects.
#' @return A `compare.loo_ss` matrix.
loo_compare_matrix.psis_loo_ss_list <- function(loos){
  tmp <- sapply(loos, function(x) {
    est <- x$estimates
    setNames(c(est), nm = c(rownames(est),
                            paste0("se_", rownames(est)),
                            paste0("subsampling_se_", rownames(est))))
  })
  colnames(tmp) <- find_model_names(loos)
  rnms <- rownames(tmp)
  comp <- tmp
  ord <- loo_compare_order(loos)
  comp <- t(comp)[ord, ]
  patts <- c("elpd", "p_", "^waic$|^looic$", "se_waic$|se_looic$")
  col_ord <- unlist(sapply(patts, function(p) grep(p, colnames(comp))),
                    use.names = FALSE)
  comp <- comp[, col_ord]
  comp
}
