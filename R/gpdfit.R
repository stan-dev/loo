gpdfit <- function(x) {
  n <- length(x)
  x <- sort.int(x, method = "quick")
  prior <- 3
  m <- 80 + floor(sqrt(n)) # note: original paper used m <- 20+floor(sqrt(n))
  b <- 1/x[n] + (1 - sqrt(m/seq_min_half(m)))/prior/x[floor(n/4 + 0.5)]
  L <- vapply1m(m, function(i) n * lx(b[i], x))
  w <- vapply1m(m, function(i) 1/sum(exp(L - L[i])))
  b <- sum(b*w)
  k <- mean.default(log(1 - b*x))
  sigma <- -k/b
  list(k=k, sigma=sigma)
}