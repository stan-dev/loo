N <- 50; K <- 10; S <- 100; a0 <- 3; b0 <- 2
p <- rbeta(1, a0, b0)
y <- rbinom(N, size = K, prob = p)
a <- a0 + sum(y); b <- b0 + N * K - sum(y)
draws <- as.matrix(rbeta(S, a, b))
data <- data.frame(y,K)
llfun <- function(data_i, draws) {
  dbinom(data_i$y, size = data_i$K, prob = draws, log = TRUE)
}
llmat_from_fn <- sapply(1:N, function(i) llfun(data[i,, drop=FALSE], draws))
