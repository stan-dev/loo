#' @param parallel Logical. If `TRUE` (the default), run on mirai daemons, which
#'   you set up beforehand with [mirai::daemons()]. Local daemons read inputs
#'   from shared memory ([mori::share()]); remote daemons receive copies.
#'
