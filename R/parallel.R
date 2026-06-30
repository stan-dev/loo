#' Package-internal state for the parallel backend
#'
#' @noRd
#' @keywords internal
#' @description
#' Holds small bits of session-scoped state used by the parallel helpers:
#'
#' * `cleanup_registered`: guards `loo_register_daemon_cleanup()` so the
#'   session-exit finalizer is only registered once.
#' * `warned_bad_daemons`: guards the malformed-config warning in
#'   `loo_persist_config()` so it is only emitted once per session.
#' * `informed_cores_ignored`: guards the "a pool is connected so `cores` is
#'   ignored" message in `with_loo_daemons()` so it is only emitted once per
#'   session.
#'
#' It also serves as the object the daemon-cleanup finalizer is attached to.
.loo_internal <- new.env(parent = emptyenv())

#' Resolve the persistent local daemon pool size from user configuration
#'
#' @noRd
#' @keywords internal
#' @description
#' Reads the opt-in "persistent local pool" size from, in precedence order:
#'
#' 1. the R option `loo.daemons`,
#' 2. the environment variable `LOO_DAEMONS`,
#' 3. otherwise the feature is off.
#'
#' This knob enables a local [mirai::daemons()] pool that is created lazily on
#' first parallel use and kept warm for the rest of the session (see
#' `with_loo_daemons()`), which avoids paying pool spawn/teardown overhead on
#' every top-level `loo()`/`psis()` call (useful for simulations, benchmarks
#' and batch/HPC scripts).
#'
#' @return A single integer `>= 2` giving the persistent pool size, or
#'   `NA_integer_` when the feature is off (unset, `0`/`1`, or a non-integer
#'   value). Genuinely malformed (non-coercible) values warn once per session
#'   and then disable the feature.
loo_persist_config <- function() {
  raw <- getOption("loo.daemons")
  if (is.null(raw)) {
    raw <- Sys.getenv("LOO_DAEMONS", unset = NA_character_)
  }
  if (length(raw) != 1L) {
    return(NA_integer_)
  }
  if (is.na(raw) || (is.character(raw) && !nzchar(trimws(raw)))) {
    # Unset / empty -> feature off.
    return(NA_integer_)
  }
  n <- suppressWarnings(as.numeric(raw))
  if (is.na(n) || !is.finite(n)) {
    # Non-numeric garbage -> off, but tell the user once that it was ignored.
    loo_warn_bad_daemons(raw)
    return(NA_integer_)
  }
  if (n < 2 || n != trunc(n)) {
    # 0/1 (serial) or a non-integer value -> feature off, silently.
    return(NA_integer_)
  }
  as.integer(n)
}

#' Warn (once per session) about a malformed persistent-pool configuration
#'
#' @noRd
#' @keywords internal
loo_warn_bad_daemons <- function(value) {
  if (isTRUE(.loo_internal$warned_bad_daemons)) {
    return(invisible(NULL))
  }
  .loo_internal$warned_bad_daemons <- TRUE
  warning(
    "Ignoring invalid persistent-pool size ", encodeString(value, quote = "'"),
    " from 'loo.daemons'/'LOO_DAEMONS'; expected a single integer >= 2.",
    call. = FALSE
  )
  invisible(NULL)
}

#' Inform (once per session) that `cores` is ignored while a pool is connected
#'
#' @noRd
#' @keywords internal
#' @description
#' Emitted by `with_loo_daemons()` when a mirai daemon pool is connected (a
#' user-managed pool, or a persistent pool left warm by an earlier call) and
#' the current call passed `cores <= 1`. When a pool is connected loo always
#' uses it, so the `cores` argument is ignored and the work runs in parallel
#' regardless of its value. This message makes that behaviour visible the
#' first time a call looks like it asked for serial execution; the usual cause
#' is relying on the default `cores = getOption("mc.cores", 1)` after setting
#' up daemons.
loo_inform_cores_ignored <- function() {
  if (isTRUE(.loo_internal$informed_cores_ignored)) {
    return(invisible(NULL))
  }
  .loo_internal$informed_cores_ignored <- TRUE
  message(
    "A mirai daemon pool is connected, so 'cores' is ignored and this call ",
    "runs in parallel on the existing pool. Call mirai::daemons(0) to stop ",
    "the pool if you want serial execution."
  )
  invisible(NULL)
}

#' Register a one-time session-exit cleanup for the persistent daemon pool
#'
#' @noRd
#' @keywords internal
#' @description
#' Attaches a finalizer (only once per session) that resets any local daemon
#' pool with `mirai::daemons(0)` when the R session exits. mirai already
#' terminates local daemons when the host session ends; this is a
#' belt-and-suspenders guard so the lazily created persistent pool never leaves
#' orphan processes behind in batch/HPC scripts.
loo_register_daemon_cleanup <- function() {
  if (isTRUE(.loo_internal$cleanup_registered)) {
    return(invisible(NULL))
  }
  .loo_internal$cleanup_registered <- TRUE
  reg.finalizer(
    .loo_internal,
    function(e) try(mirai::daemons(0), silent = TRUE),
    onexit = TRUE
  )
  invisible(NULL)
}

#' Evaluate parallel work with an appropriate mirai daemon pool
#'
#' @noRd
#' @keywords internal
#' @description
#' Central entry point used by loo's parallel code paths to ensure a
#' [mirai::daemons()] pool exists for the duration of a computation. It is
#' deliberately a good citizen of the user's session:
#'
#' * A daemon pool is already connected (e.g. the user called
#'   [mirai::daemons()] themselves, possibly with remote/HPC daemons, or a
#'   persistent pool was left warm by an earlier call): `code` runs on the
#'   existing pool, which is left untouched. **A connected pool always wins**,
#'   so the `cores` argument is ignored entirely in this case (loo uses the
#'   pool regardless of `cores`). If `cores <= 1` here -- i.e. the call looked
#'   like it requested serial execution -- a one-time-per-session message
#'   notes that `cores` is being ignored (via `loo_inform_cores_ignored()`).
#' * `cores <= 1` with no pool connected: runs `code` serially without
#'   touching daemons. (A configured `loo.daemons` pool is created lazily only
#'   when `cores > 1`, so serial-only work never starts workers.)
#' * Otherwise, if the user opted in to a persistent session pool via the
#'   `loo.daemons` option or `LOO_DAEMONS` environment variable (see
#'   `loo_persist_config()`): a local pool of that size is created lazily on
#'   this first parallel call and left warm for the rest of the session, with a
#'   session-exit finalizer registered for cleanup. Subsequent calls reuse it
#'   via the existing-pool branch above.
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
#' @param cores Integer number of cores requested by the user. When no pool is
#'   connected this is the per-call "enable parallel" switch (and the size of
#'   the ephemeral pool). When a pool is already connected it is ignored --
#'   loo always uses the connected pool. The persistent pool size, when
#'   enabled, comes from `loo_persist_config()` rather than from `cores`.
#' @param code Expression to evaluate. Lazily evaluated in the calling
#'   environment, after any daemon pool has been set up.
#' @return The value of `code`.
with_loo_daemons <- function(cores, code) {
  if (loo_has_pool()) {
    # A pool is connected (user-managed, or a warm persistent pool): use it
    # regardless of `cores`. This always wins over the options below. If the
    # call looked like it asked for serial work, note once that `cores` is
    # being ignored.
    if (cores <= 1) {
      loo_inform_cores_ignored()
    }
    return(code)
  }
  if (cores <= 1) {
    # No pool connected and no parallelism requested: run serially. A
    # configured loo.daemons pool is created lazily only when cores > 1, so
    # serial-only work never starts workers.
    return(code)
  }
  persist <- loo_persist_config()
  if (!is.na(persist)) {
    # Opt-in persistent pool: create once, leave warm for the session, and
    # register a finalizer to tidy up at session exit. No per-call teardown.
    mirai::daemons(persist)
    loo_register_daemon_cleanup()
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
#' when no daemon pool is connected.
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
#' @param cores Integer number of cores requested by the user. Used only as a
#'   fallback worker count for chunking when (unexpectedly) no pool is
#'   connected; it does not gate execution. Parallelism is used whenever a
#'   daemon pool is connected (the connected pool always wins over `cores`).
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

  if (!loo_has_pool()) {
    # Serial path: identical behaviour to a plain lapply() with the broadcast
    # and constant arguments supplied by name. When a pool is connected loo
    # uses it regardless of `cores`; the decision to create one lives upstream
    # in with_loo_daemons().
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
