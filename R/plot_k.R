#' Plot the estimates of the shape parameter k of the generalized Pareto distribution
# #'
# #' @export
# #'
# .plot_k <- function(x, cex = 1) {
#   if (!("pareto_k" %in% names(x)))
#     stop("No k estimates found")
#
#   k <- x$pareto_k
#
#   yl <- expression(paste("Shape parameter ", italic(k)))
#   xl <- expression(paste("Data ", italic(i)))
#   par(mfrow = c(1,2))
#   plot(k, xlab = xl, ylab = yl, type = "n", bty = "l", yaxt = "n")
#   axis(side = 2, las = 1)
#   inrange <- function(a, rr) a >= rr[1] & a <= rr[2]
#   krange <- range(k)
#   for (val in c(0, 0.5, 1)) {
#     if (inrange(val, krange))
#       abline(h = val, lty = 2)
#   }
#   cnms <- c("dodgerblue", "darkorchid", "firebrick1")
#   clrs <- ifelse(inrange(k, c(-Inf, 0.5)), cnms[1],
#                  ifelse(inrange(k, c(0.5, 1)), cnms[2], cnms[3]))
#   points(k, col = clrs, pch = 20, cex = cex)
#   hist(k, col = "dodgerblue", border = "dodgerblue3")
# }
# #
