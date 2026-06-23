# Synthetic multiclass data used across the multiclass tests. Returns a list
# with logits, the softmax probabilities, and integer labels drawn from those
# probabilities.
make_multiclass <- function(n = 200, k = 3, seed = 1) {
  set.seed(seed)
  logits <- matrix(stats::rnorm(n * k), ncol = k)
  prob <- t(apply(logits, 1L, function(row) {
    z <- exp(row - max(row))
    z / sum(z)
  }))
  labels <- apply(prob, 1L, function(row) sample.int(k, 1L, prob = row))
  list(logits = logits, prob = prob, labels = labels, k = k)
}
