with_two_daemons <- function(code) {
  mirai::daemons(0)
  mirai::daemons(2, dispatcher = FALSE)
  withr::defer({
    mirai::daemons(0)
    mirai::daemons(1, dispatcher = FALSE)
  })
  force(code)
}
