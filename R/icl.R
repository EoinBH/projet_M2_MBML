# Complete-data log-likelihood and ICL criterion.
#   ICL = loglik_complete - (K-1)/2 log(n) - (L-1)/2 log(d) - K*L log(n*d)

# Complete-data log-likelihood and ICL for given partitions and parameters.
complete_loglik <- function(Ilist, z, w, theta, K, L) {
  m <- length(Ilist)
  N <- block_counts(Ilist, indicator(z, K), indicator(w, L))
  logP <- log_cub_table(theta$pi, theta$xi, m)
  sum(tabulate(z, K) * log(theta$alpha)) + sum(tabulate(w, L) * log(theta$beta)) +
    sum(N * logP)
}

icl_value <- function(loglik, n, d, K, L) {
  loglik - (K - 1) / 2 * log(n) - (L - 1) / 2 * log(d) - K * L * log(n * d)
}
