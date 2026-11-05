#' @param parallel Logical. If `TRUE` (the default), run on mirai daemons, which
#'   you set up beforehand with [mirai::daemons()]. This needs the \pkg{mirai}
#'   and \pkg{mori} packages, which are only suggested. Local daemons read
#'   inputs from shared memory ([mori::share()]); remote daemons receive copies.
#'   Daemons start with an empty global environment, so a log-likelihood
#'   function must be self-contained or get what it needs through `data`,
#'   `draws` or `...`.
#'
