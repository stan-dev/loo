#' Evaluate parallel work with an appropriate mirai daemon pool
#'
#' @noRd
#' @keywords internal
#' @description
#' Central entry point used by loo's parallel code paths to ensure a
#' [mirai::daemons()] pool exists for the duration of a computation. It is
#' deliberately a good citizen of the user's session:
#'
#' * `cores <= 1`: runs `code` serially without touching daemons.
#' * A daemon pool is already configured (e.g. the user called
#'   [mirai::daemons()] themselves, possibly with remote/HPC daemons): `code`
#'   runs on the existing pool, which is left untouched.
#' * Otherwise: a pool of `cores` local daemons is created for the duration of
#'   `code` and automatically reset afterwards (via the scoped
#'   `with(mirai::daemons(), ...)` method), so no daemon processes are left
#'   running once the call returns.
#'
#' This keeps a single pool alive across the whole top-level computation
#' (rather than spinning daemons up and down for each unit of work) while
#' respecting any pool the user has already declared. Because it reuses an
#' existing pool, it is safe to nest: an inner call made while an outer call
#' already established a pool simply reuses it instead of creating another.
#'
#' @param cores Integer number of cores requested by the user.
#' @param code Expression to evaluate. Lazily evaluated in the calling
#'   environment, after any daemon pool has been set up.
#' @return The value of `code`.
with_loo_daemons <- function(cores, code) {
  if (cores <= 1 || loo_has_pool()) {
    # Serial work, or reuse the daemon pool the user (or an outer loo call)
    # already configured.
    return(code)
  }
  # No pool configured: create one scoped to this computation and reset it on
  # exit. `code` (including result collection via `[]`) is forced before the
  # daemons are torn down.
  with(mirai::daemons(cores), code)
}

#' Is a mirai daemon pool currently connected?
#'
#' @noRd
#' @keywords internal
#' @return `TRUE` if at least one daemon connection exists for the active
#'   compute profile, otherwise `FALSE`.
loo_has_pool <- function() {
  conns <- tryCatch(mirai::status()$connections, error = function(e) 0L)
  isTRUE(as.integer(conns) > 0L)
}

#' Number of workers available for chunking decisions
#'
#' @noRd
#' @keywords internal
#' @description
#' Returns the number of connected daemons when a pool exists (so chunking
#' matches the actual worker count, including user-supplied or remote pools),
#' otherwise falls back to the requested `cores`.
#' @param cores Integer number of cores requested by the user.
loo_n_workers <- function(cores) {
  conns <- tryCatch(mirai::status()$connections, error = function(e) 0L)
  conns <- as.integer(conns)
  if (length(conns) != 1L || is.na(conns) || conns < 1L) {
    return(as.integer(cores))
  }
  conns
}

#' Is the active daemon pool on the local machine?
#'
#' @noRd
#' @keywords internal
#' @description
#' Determines whether shared memory ([mori::share()]) can be used safely with
#' the active pool. Shared memory only works when workers run on the same
#' physical machine, so we only treat same-host transports as local:
#'
#' * `abstract://` and `ipc://` are same-machine inter-process transports used
#'   by local [mirai::daemons()] pools, so these are treated as local.
#' * `tcp://` (and anything else) may be a remote pool, or the host URL that
#'   remote SSH/HPC daemons dial back to, so it is treated as **not** local.
#'   loo then falls back to ordinary serialization instead of shared memory.
#'
#' This is intentionally conservative: an incorrect "local" classification
#' would produce wrong results on a remote pool, whereas an incorrect "remote"
#' classification merely forgoes the zero-copy optimisation.
#' @return `TRUE` if the pool is confirmed local, otherwise `FALSE`.
loo_pool_is_local <- function() {
  urls <- tryCatch(mirai::status()$daemons, error = function(e) NULL)
  if (!is.character(urls) || length(urls) == 0L) {
    return(FALSE)
  }
  all(grepl("^(abstract|ipc)://", urls))
}

#' Map a worker over elements, serially or across a mirai daemon pool
#'
#' @noRd
#' @keywords internal
#' @description
#' Single cross-platform entry point for loo's per-observation parallelism.
#' Replaces the previous platform-branching
#' [parallel::mclapply()] / [parallel::parLapply()] code paths with a single
#' [mirai::mirai_map()] path, while preserving the serial [lapply()] behaviour
#' when no parallelism is requested or available.
#'
#' Object transport is chosen automatically:
#'
#' * `broadcast` objects are reused identically by every element (e.g. the
#'   posterior `draws` matrix). On a local pool they are written once into
#'   shared memory with [mori::share()] so each daemon maps the same physical
#'   pages (zero-copy). On a remote pool, where shared memory is unavailable,
#'   they are serialized instead; chunking bounds the number of copies sent to
#'   roughly one per worker.
#' * Small per-call arguments are passed through `...`.
#'
#' @param X A vector or list to iterate over. Each element is passed as the
#'   first argument to `FUN`.
#' @param FUN Worker function. Called as `FUN(x, <broadcast>, <...>)`; the
#'   names in `broadcast` and `...` must match `FUN`'s formals.
#' @param ... Small constant arguments forwarded to `FUN` for every element.
#' @param cores Integer number of cores requested by the user. Parallelism is
#'   only used when `cores > 1` and a daemon pool is connected.
#' @param broadcast Named list of large objects reused by every element. See
#'   Description for how these are transported.
#' @param chunk Chunking strategy. `"auto"` (default) splits `X` into roughly
#'   one chunk per worker to amortise per-task overhead -- best for cheap
#'   per-element work over many elements. `"never"` dispatches one task per
#'   element for finer load balancing -- best for expensive, uneven per-element
#'   work. `"never"` is automatically promoted to `"auto"` on a remote pool
#'   that carries `broadcast` objects, to avoid re-sending them per task.
#' @return A list of `FUN` results in the same order as `X`.
loo_map <- function(X, FUN, ..., cores = 1L, broadcast = list(),
                    chunk = c("auto", "never")) {
  chunk <- match.arg(chunk)
  dots <- list(...)

  if (!(cores > 1L && loo_has_pool())) {
    # Serial path: identical behaviour to a plain lapply() with the broadcast
    # and constant arguments supplied by name.
    return(do.call(lapply, c(list(X, FUN), broadcast, dots)))
  }

  local_pool <- loo_pool_is_local()
  if (length(broadcast) > 0L) {
    if (local_pool) {
      # Zero-copy: write once to shared memory, ship tiny references.
      broadcast <- lapply(broadcast, mori::share)
    } else if (chunk == "never") {
      # Remote pool: avoid re-serializing large broadcast objects once per
      # task by collapsing to one chunk per worker instead.
      chunk <- "auto"
    }
  }
  const_args <- c(broadcast, dots)

  if (chunk == "never") {
    return(
      mirai::mirai_map(
        X,
        function(.x, .FUN, .const) do.call(.FUN, c(list(.x), .const)),
        .args = list(.FUN = FUN, .const = const_args)
      )[mirai::.stop]
    )
  }

  n_chunks <- min(loo_n_workers(cores), length(X))
  positions <- parallel::splitIndices(length(X), n_chunks)
  chunks <- lapply(positions, function(p) X[p])
  chunk_results <- mirai::mirai_map(
    chunks,
    function(.chunk, .FUN, .const) {
      lapply(.chunk, function(.x) do.call(.FUN, c(list(.x), .const)))
    },
    .args = list(.FUN = FUN, .const = const_args)
  )[mirai::.stop]
  # splitIndices() returns contiguous ascending groups, so concatenating the
  # per-chunk lists restores the original order of X.
  do.call(c, chunk_results)
}
