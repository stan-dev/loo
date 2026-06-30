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
#'   Parallelism is implemented with the \pkg{mirai} package. There are three
#'   ways to control the backend, in increasing order of precedence:
#'   * `cores > 1` (the default behaviour): a local daemon pool is created for
#'     the duration of the call and torn down automatically when it returns, so
#'     no worker processes are left running. This is convenient for one-off
#'     calls but pays a small pool start-up/teardown cost on every call.
#'   * `options(loo.daemons = k)` or the environment variable `LOO_DAEMONS=k`
#'     (with `k >= 2`): opt in to a *persistent* local pool of `k` daemons that
#'     is created lazily on the first parallel call and kept warm for the rest
#'     of the R session, then cleaned up automatically at session exit. This
#'     avoids repeated pool start-up/teardown and is ideal for simulations,
#'     benchmarks, and batch/HPC scripts that call `loo()`/`psis()` many times.
#'     Local pools automatically use zero-copy shared memory (via \pkg{mori})
#'     for the shared posterior draws.
#'   * A pool you configure yourself with [mirai::daemons()] (including
#'     remote/SSH/HPC daemons via `mirai::daemons(url = ...)`) always takes
#'     precedence: loo reuses it and never tears it down.
#'
