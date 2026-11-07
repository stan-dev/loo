#' @param parallel Logical. If `TRUE`, run on mirai daemons, which you set up
#'   beforehand with [mirai::daemons()]. The default is
#'   `getOption("loo.parallel", TRUE)`. This needs the \pkg{mirai} and \pkg{mori}
#'   packages, which are only suggested. Local daemons read inputs larger than
#'   `getOption("loo.share_bytes", 4 * 2^20)` bytes (4 MiB) from shared memory
#'   ([mori::share()]) and receive copies of smaller ones; remote daemons
#'   always receive copies. Daemons start with an empty global environment, so a
#'   log-likelihood function must be self-contained or get what it needs
#'   through `data`, `draws` or `...`.
#'
