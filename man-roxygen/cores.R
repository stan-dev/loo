#' @param cores **Deprecated.** The number of cores to use for parallelization.
#'   This argument will be removed in a future release. A deprecation warning is
#'   emitted when `cores` is passed explicitly (even `cores = 1`), when
#'   `options(mc.cores)` is set, or when the resolved value is `cores > 1`.
#'   When `cores > 1` and no [mirai::daemons()] pool is connected, a per-call
#'   local daemon pool of that size is started for the duration of the call.
#'   Prefer [mirai::daemons()] for session-level pools or [loo_mirai()] for
#'   map-style workflows over many models.
#'
#'   This defaults to `getOption("mc.cores", 1)`, but **`options(mc.cores)` is
#'   also deprecated** and will be removed in a future release (it still works
#'   on this branch as a bridge to per-call daemon pools). **As of version 2.0.0
#'   the default is 1 core when `mc.cores` is not set**.
#'   * Note for Windows 10 users: it is __strongly__
#'     [recommended](https://github.com/stan-dev/loo/issues/94) to avoid using
#'     the `.Rprofile` file to set `mc.cores`.
#'
#'   See `vignette("loo2-parallel")` for the mirai-based replacement workflow.
#'
