# waic and loo helpers ----------------------------------------------------

#' @importFrom matrixStats colLogSumExps
logColMeansExp <- function(x) {
  # should be more stable than log(colMeans(exp(x)))
  S <- nrow(x)
  colLogSumExps(x) - log(S)
}

logColMeansExp_ll <- function(fun, args) {
  # should be more stable than log(colMeans(exp(x)))
  logS <- log(args$S)
  clse <- vapply(seq_len(args$N), FUN = function(i) {
    logSumExp(fun(i = i, data = args$data[i,,drop=FALSE], draws = args$draws))
  }, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
  clse - logS
}
colVars_ll <- function(fun, args) {
  vapply(seq_len(args$N), FUN = function(i) {
    var(as.vector(fun(i = i, data = args$data[i,,drop=FALSE], draws = args$draws)))
  }, FUN.VALUE = numeric(1), USE.NAMES = FALSE)
}

totals <- function(pointwise) {
  N <- length(pointwise[[1L]])
  total  <- unlist_lapply(pointwise, sum)
  se <- sqrt(N * unlist_lapply(pointwise, var))
  as.list(c(total, se))
}

#' @importFrom matrixStats colVars
pointwise_waic <- function(log_lik, llfun = NULL, llargs = NULL) {
  if (!missing(log_lik)) {
    lpd <- logColMeansExp(log_lik)
    p_waic <- colVars(log_lik)
  } else {
    if (is.null(llfun) || is.null(llargs))
      stop("Either log_lik or llfun and llargs must be specified")
    lpd <- logColMeansExp_ll(llfun, llargs)
    p_waic <- colVars_ll(llfun, llargs)
  }
  elpd_waic <- lpd - p_waic
  waic <- -2 * elpd_waic
  pointwise <- nlist(elpd_waic, p_waic, waic)
  out <- totals(pointwise)
  nms <- names(pointwise)
  names(out) <- c(nms, paste0("se_", nms))
  out$pointwise <- cbind_list(pointwise)
  out
}
pointwise_loo <- function(psis, log_lik, llfun = NULL, llargs = NULL) {
  if (!missing(log_lik)) lpd <- logColMeansExp(log_lik)
  else {
    if (is.null(llfun) || is.null(llargs))
      stop("Either log_lik or llfun and llargs must be specified")
    lpd <- logColMeansExp_ll(llfun, llargs)
  }
  elpd_loo <- psis$loos
  p_loo <- lpd - elpd_loo
  looic <- -2 * elpd_loo
  pointwise <- nlist(elpd_loo, p_loo, looic)
  out <- totals(pointwise)
  nms <- names(pointwise)
  names(out) <- c(nms, paste0("se_", nms))
  out$pointwise <- cbind_list(pointwise)
  out$pareto_k <- psis$pareto_k
  out
}

# psis helpers ------------------------------------------------------------

# inverse-CDF of generalized Pareto distribution (formula from Wikipedia)
qgpd <- function(p, xi = 1, mu = 0, sigma = 1, lower.tail = TRUE) {
  if (sigma <= 0) return(rep(NaN, length(p)))
  if (!lower.tail)
    p <- 1 - p
  mu + sigma * ((1 - p)^(-xi) - 1) / xi
}

lx <- function(a,x) {
  a <- -a
  k <- sapply(a, FUN = function(y) mean(log1p(y * x)))
  log(a/k) - k - 1
}

lw_cutpoint <- function(y, wcp, min_cut) {
  if (min_cut < log(.Machine$double.xmin))
    min_cut <- -700
  cp <- quantile(y, 1 - wcp, names = FALSE)
  max(cp, min_cut)
}
lw_truncate <- function(y, wtrunc) {
  if (wtrunc == 0)
    return(y)
  logS <- log(length(y))
  lwtrunc <- wtrunc * logS - logS + logSumExp(y)
  y[y > lwtrunc] <- lwtrunc
  y
}
#' @importFrom matrixStats logSumExp
lw_normalize <- function(y) {
  y - logSumExp(y)
}

# print helpers -----------------------------------------------------------
#' @importFrom graphics abline axis plot points text
.fr <- function(x, digits) format(round(x, digits), nsmall = digits)
.warn <- function(..., call. = FALSE) warning(..., call. = call.)
k_warnings <- function(k, digits = 1) {
  brks <- c(-Inf, 0.5, 1, Inf)
  kcut <- cut(k, breaks = brks)
  count <- table(kcut)
  prop <- prop.table(count)
  if (sum(count[2:3]) == 0) {
    cat("\nAll Pareto k estimates OK (k < 0.5)")
  } else {
    if (count[2] != 0) {
      txt2 <- "%) Pareto k estimates between 0.5 and 1"
      .warn(paste0(count[2], " (", .fr(100 * prop[2], digits), txt2))
    }
    if (count[3] != 0) {
      txt3 <- "%) Pareto k estimates greater than 1"
      .warn(paste0(count[3], " (", .fr(100 * prop[3], digits), txt3))
    }
    .warn("See PSIS-LOO description (?'loo-package') for more information")
  }
  invisible(NULL)
}

plot_k <- function(k, ..., label_points = FALSE) {
  inrange <- function(a, rr) a >= rr[1L] & a <= rr[2L]
  yl <- expression(paste("Shape parameter ", italic(k)))
  xl <- expression(paste("Data ", italic(i)))
  plot(k, xlab = xl, ylab = yl, type = "n", bty = "l", yaxt = "n")
  axis(side = 2, las = 1)
  krange <- range(k)
  for (val in c(0, 0.5, 1)) {
    if (inrange(val, krange))
      abline(h = val, col = "#b17e64", lty = 2, lwd = 1)
  }
  hex_clrs <- c("#6497b1", "#005b96", "#03396c")
  brks <- c(-Inf, 0.5, 1)
  clrs <- ifelse(inrange(k, brks[1:2]), hex_clrs[1],
                 ifelse(inrange(k, brks[2:3]), hex_clrs[2L], hex_clrs[3L]))
  if (all(k < 0.5) || !label_points) {
    points(k, col = clrs, pch = 3, cex = .6)
    return(invisible())
  } else {
    points(k[k < 0.5], col = clrs[k < 0.5], pch = 3, cex = .6)
    sel <- !inrange(k, brks[1:2])
    dots <- list(...)
    txt_args <- c(list(x = seq_along(k)[sel], y = k[sel],
                       labels = seq_along(k)[sel]),
                  if (length(dots) > 0) dots)
    if (!("adj" %in% names(txt_args))) txt_args$adj <- 2/3
    if (!("cex" %in% names(txt_args))) txt_args$cex <- 0.75
    if (!("col" %in% names(txt_args))) txt_args$col <- clrs[sel]
    do.call("text", txt_args)
  }
}


