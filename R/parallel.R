#' Package-internal state for the parallel backend
#'
#' @noRd
#' @keywords internal
#' @description
#' Holds small bits of session-scoped state used by the parallel helpers:
#'
#' * `informed_cores_ignored`: guards the "a pool is connected so `cores` is
#'   ignored" message in `with_loo_daemons()` so it is only emitted once per
#'   session.
#' * `warned_cores_deprecated`: guards the `cores` argument deprecation warning
#'   in `loo_deprecate_cores()` so it is only emitted once per session.
#' * `warned_cores_in_mirai`: guards the `cores`-in-`args_list` warning in
#'   `loo_mirai()` so it is only emitted once per session.
#'
#' It also serves as the object attached to session-scoped parallel state.
.loo_internal <- new.env(parent = emptyenv())

#' Inform (once per session) that `cores` is ignored while a pool is connected
#'
#' @noRd
#' @keywords internal
#' @description
#' Emitted by `with_loo_daemons()` when a mirai daemon pool is connected (e.g.
#' the user called [mirai::daemons()] themselves) and the current call passed
#' `cores <= 1`. When a pool is connected loo always uses it, so the `cores`
#' argument is ignored and the work runs in parallel regardless of its value.
#' This message makes that behaviour visible the
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

#' Warn (once per session) that the `cores` argument is deprecated
#'
#' @noRd
#' @keywords internal
#' @description
#' Emitted from `loo_cores()` when the user passes `cores` explicitly (even
#' `cores = 1`), when `options(mc.cores)` is set, or when the resolved value
#' is `cores > 1`. The argument/option still works for now by starting a scoped
#' local [mirai::daemons()] pool when `cores > 1`, but users should migrate to
#' [mirai::daemons()] or [loo_mirai()].
#' @param cores Integer number of cores requested by the user.
#' @param explicit Whether `cores` was passed explicitly by the user.
#' @param mc_cores_set Whether `options(mc.cores)` is set and `cores` was not
#'   passed explicitly.
loo_deprecate_cores <- function(cores, explicit = FALSE, mc_cores_set = FALSE) {
  if (!explicit && !mc_cores_set && cores <= 1L) {
    return(invisible(NULL))
  }
  if (isTRUE(.loo_internal$warned_cores_deprecated)) {
    return(invisible(NULL))
  }
  .loo_internal$warned_cores_deprecated <- TRUE
  warning(
    "The 'cores' argument and options('mc.cores') are deprecated and will be ",
    "removed in a future release. For one-off parallel calls, 'cores' / ",
    "'mc.cores' currently start a per-call local daemon pool when > 1; use ",
    "mirai::daemons() before calling loo functions, or loo_mirai(..., ",
    "n_daemons = k) for map-style workflows over many models.",
    call. = FALSE
  )
  invisible(NULL)
}

#' Warn (once per session) about `cores` in `loo_mirai()` `args_list`
#'
#' @noRd
#' @keywords internal
loo_deprecate_cores_in_mirai <- function() {
  if (isTRUE(.loo_internal$warned_cores_in_mirai)) {
    return(invisible(NULL))
  }
  .loo_internal$warned_cores_in_mirai <- TRUE
  warning(
    "Passing 'cores' in args_list to loo_mirai() is deprecated and ignored. ",
    "Parallelism is controlled by loo_mirai() via n_daemons or an existing ",
    "mirai::daemons() pool.",
    call. = FALSE
  )
  invisible(NULL)
}

#' Prepare argument lists for [loo_mirai()] map jobs
#'
#' @noRd
#' @keywords internal
#' @description
#' Strips any user-supplied `cores` element (with a deprecation warning) and
#' forces `cores = 1` so each mapped job runs serially on its worker. Outer
#' parallelism comes only from `loo_mirai()` itself.
loo_mirai_prepare_args <- function(args) {
  if ("cores" %in% names(args)) {
    loo_deprecate_cores_in_mirai()
    args$cores <- NULL
  }
  args$cores <- 1L
  args
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
#'   [mirai::daemons()] themselves, possibly with remote/HPC daemons): `code`
#'   runs on the existing pool, which is left untouched. **A connected pool
#'   always wins**, so the `cores` argument is ignored entirely in this case.
#'   If `cores <= 1` here -- i.e. the call looked like it requested serial
#'   execution -- a one-time-per-session message notes that `cores` is being
#'   ignored (via `loo_inform_cores_ignored()`).
#' * `cores <= 1` with no pool connected: runs `code` serially without touching
#'   daemons.
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
#'   loo always uses the connected pool.
#' @param code Expression to evaluate. Lazily evaluated in the calling
#'   environment, after any daemon pool has been set up.
#' @return The value of `code`.
with_loo_daemons <- function(cores, code) {
  if (loo_has_pool()) {
    if (cores <= 1) {
      loo_inform_cores_ignored()
    }
    return(code)
  }
  if (cores <= 1) {
    return(code)
  }
  with(mirai::daemons(cores), code)
}

#' Run a \pkg{loo} function over argument lists in parallel with mirai
#'
#' Intended for map-style workflows such as computing [loo()], [psis()], or
#' other \pkg{loo} functions for many models in one session.
#'
#' @param fun A \pkg{loo} function (e.g. [loo()], [psis()]) to apply to each
#'   element of `args_list`.
#' @param args_list A list of lists. Each element is passed to `fun` via
#'   `do.call(fun, element)`.
#' @param n_daemons Number of local daemons to start for this call. When
#'   `NULL` (the default), an existing [mirai::daemons()] pool is reused if
#'   connected; otherwise the work runs serially. When an integer `>= 2`, a
#'   local pool of that size is created for the duration of the call and
#'   torn down automatically when it returns (unless a pool was already
#'   connected, in which case it is reused and left untouched).
#' @return A list of results in the same order as `args_list`.
#' @details Parallelism is controlled by `loo_mirai()` (via `n_daemons` or an
#'   existing daemon pool), not by a `cores` argument on the mapped calls.
#'   Any `cores` element in `args_list` is deprecated, ignored, and replaced
#'   with `cores = 1` so each job runs serially on its worker.
#' @export
#' @examples
#' r_eff <- relative_eff(exp(example_loglik_array()))
#' ll_list <- list(example_loglik_matrix(), example_loglik_matrix())
#' args_list <- lapply(ll_list, function(LL) {
#'   list(log_ratios = -LL, r_eff = r_eff)
#' })
#'
#' # Serial when no pool is connected and n_daemons is NULL
#' loo_mirai(psis, args_list)
#'
#' \dontrun{
#' # Convenience wrapper: start workers, map, tear down
#' loo_mirai(psis, args_list, n_daemons = 2)
#'
#' # Or manage the pool yourself and call loo functions directly
#' mirai::daemons(4)
#' loo_mirai(loo, lapply(ll_list, function(LL) {
#'   list(x = LL, r_eff = r_eff)
#' }))
#' mirai::daemons(0)
#' }
loo_mirai <- function(fun, args_list, n_daemons = NULL) {
  if (!is.function(fun)) {
    stop("'fun' must be a function.", call. = FALSE)
  }
  if (!identical(environment(fun), asNamespace("loo"))) {
    stop(
      "'fun' must be a loo function (e.g. loo(), psis()).",
      call. = FALSE
    )
  }
  if (!is.list(args_list)) {
    stop("'args_list' must be a list.", call. = FALSE)
  }

  run_map <- function() {
    prepared <- lapply(args_list, loo_mirai_prepare_args)
    if (!loo_has_pool()) {
      return(lapply(prepared, function(a) do.call(fun, a)))
    }
    mirai::mirai_map(
      prepared,
      function(.args, .fun) do.call(.fun, .args),
      .args = list(.fun = fun)
    )[mirai::.stop]
  }

  if (!is.null(n_daemons)) {
    n_daemons <- as.integer(n_daemons)
    if (length(n_daemons) != 1L || is.na(n_daemons) || n_daemons < 1L) {
      stop(
        "'n_daemons' must be a single positive integer or NULL.",
        call. = FALSE
      )
    }
    if (n_daemons == 1L || loo_has_pool()) {
      return(run_map())
    }
    return(with(mirai::daemons(n_daemons), run_map()))
  }

  run_map()
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
    return(do.call(lapply, c(list(X, FUN), broadcast, dots)))
  }

  local_pool <- loo_pool_is_local()
  if (length(broadcast) > 0L) {
    if (local_pool) {
      broadcast <- lapply(broadcast, mori::share)
    } else if (chunk == "never") {
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
  do.call(c, chunk_results)
}
