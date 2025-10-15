#' @param cores The number of cores to use for parallelization. This defaults to
#'   the option `mc.cores` which can be set for an entire R session by
#'   `options(mc.cores = NUMBER)`. The old option `loo.cores` is now
#'   deprecated but will be given precedence over `mc.cores` until
#'   `loo.cores` is removed in a future release. **As of version
#'   2.0.0 the default is now 1 core if `mc.cores` is not set**, but we
#'   recommend using as many (or close to as many) cores as possible.
#'   * Note for Windows 10 users: it is __strongly__
#'     [recommended](https://github.com/stan-dev/loo/issues/94) to avoid using
#'     the `.Rprofile` file to set `mc.cores` (using the `cores` argument or
#'     setting `mc.cores` interactively or in a script is fine).
#'
