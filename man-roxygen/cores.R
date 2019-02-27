#' @param cores The number of cores to use for parallelization. This defaults to
#'   the option \code{mc.cores} which can be set for an entire R session by
#'   \code{options(mc.cores = NUMBER)}. The old option \code{loo.cores} is now
#'   deprecated but will be given precedence over \code{mc.cores} until
#'   \code{loo.cores} is removed in a future release. \strong{As of version
#'   2.0.0 the default is now 1 core if \code{mc.cores} is not set, but we
#'   recommend using as many (or close to as many) cores as possible.}
#'   Note for Windows 10 users: it is recommended to avoid using the
#'   \code{.Rprofile} file to set \code{mc.cores} (using the \code{cores}
#'   argument or setting \code{mc.cores} interactively or in a script is fine).
#'
