# waic and loo helpers ----------------------------------------------------

#' @importFrom matrixStats colLogSumExps
logColMeansExp <- function(x) {
  # should be more stable than log(colMeans(exp(x)))
  S <- nrow(x)
  colLogSumExps(x) - log(S)
}

#' @importFrom matrixStats colVars
pointwise_waic <- function(log_lik) {
  lpd <- logColMeansExp(log_lik)
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2 * elpd_waic
  nlist(elpd_waic, p_waic, waic)
}
pointwise_loo <- function(log_lik, vgis) {
  # vgis is output from vgisloo()
  lpd <- logColMeansExp(log_lik)
  elpd_loo <- vgis$loos
  p_loo <- lpd - elpd_loo
  looic <- -2 * elpd_loo
  nlist(elpd_loo, p_loo, looic)
}
totals <- function(pointwise) {
  N <- length(pointwise[[1L]])
  total  <- unlist_lapply(pointwise, sum)
  se <- sqrt(N * unlist_lapply(pointwise, var))
  as.list(c(total, se))
}


# VGIS helpers ------------------------------------------------------------

# inverse-CDF of generalized Pareto distribution (formula from Wikipedia)
qgpd <- function(p, xi = 1, mu = 0, beta = 1, lower.tail = TRUE) {
  if (!lower.tail)
    p <- 1 - p
  mu + beta * ((1 - p)^(-xi) - 1) / xi
}

# lx <- function(a, x) {
#   k <- mean.default(log1p(-a * x))
#   log(-a / k) - k - 1
# }
lx <- function(a, x) {
  # vectorized version
  b <- -a
  bx <- outer(x, b)
  d <- dim(bx)
  k <- .colMeans(log1p(bx), d[1], d[2])
  log(b / k) - k - 1
}

#' @importFrom matrixStats logSumExp
lw_normalize <- function(y) {
  y - logSumExp(y)
}
lw_truncate <- function(y, wtrunc) {
  if (wtrunc == 0)
    return(y)
  logS <- log(length(y))
  lwtrunc <- wtrunc * logS - logS + logSumExp(y)
  y[y > lwtrunc] <- lwtrunc
  y
}
lw_cutpoint <- function(y, wcp, min_cut) {
  cp <- quantile(y, 1 - wcp, names = FALSE)
  max(cp, min_cut)
}

# The parallelization functions mclapply and parLapply return a list of lists:
# vgis is a list of length N=ncol(lw). Each of the N elements of vgis is itself
# a list of length 2. In each of these N lists of length 2 the first component
# is a vector of length S=nrow(lw) containing the modified log weights and the
# second component is the estimate of the pareto shape parameter k. This
# function cbinds the log weight vectors into a matrix and combines the k
# estimates into a vector.
.vgis_out <- function(vgis) {
  ux <- unlist(vgis, recursive = FALSE)
  lws <- grepl("lw", names(ux))
  lw_smooth <- cbind_list(ux[lws])
  pareto_k <- unlist(ux[!lws])
  nlist(lw_smooth, pareto_k)
}


# print helpers -----------------------------------------------------------
.fr <- function(x, digits) format(round(x, digits), nsmall = digits)
.warn <- function(..., call. = FALSE) warning(..., call. = call.)
.k_warnings <- function(k, digits = 1) {
  brks <- c(-Inf, 0.5, 1, Inf)
  kcut <- cut(k, breaks = brks, right = FALSE)
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
    .warn("See VGIS-LOO description (?'loo-package') for more information")
  }
  invisible(NULL)
}

.plot_k <- function(k) {
  inrange <- function(a, rr) a >= rr[1] & a <= rr[2]
  yl <- expression(paste("Shape parameter ", italic(k)))
  xl <- expression(paste("Data ", italic(i)))
  plot(k, xlab = xl, ylab = yl, type = "n", bty = "l", yaxt = "n")
  axis(side = 2, las = 1)
  krange <- range(k)
  for (val in c(0, 0.5, 1)) {
    if (inrange(val, krange))
      abline(h = val, col = "#b17e64", lty = 2, lwd = 0.75)
  }
  hex_clrs <- c("#6497b1", "#005b96", "#03396c")
  brks <- c(-Inf, 0.5, 1)
  clrs <- ifelse(inrange(k, brks[1:2]), hex_clrs[1],
                 ifelse(inrange(k, brks[2:3]), hex_clrs[2], hex_clrs[3]))
  points(k, col = clrs, pch = 3, cex = .6)
}


