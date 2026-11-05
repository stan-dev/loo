#' @param parallel Logical. If `TRUE` (the default), run in parallel on the
#'   local machine or on remote workers. If a default pool has been set with
#'   [mizu::mizu_set_default_pool()], that pool is used: mizu workers share
#'   memory with the session and give the lowest overhead for fully local
#'   computation. Otherwise the work runs on mirai daemons, which you set up
#'   beforehand with [mirai::daemons()]. Local workers of either kind read
#'   inputs from shared memory ([mori::share()]); remote mirai daemons receive
#'   copies.
#'
