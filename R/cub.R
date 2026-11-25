# The CUB distribution: log-probability table and estimation of (pi, xi)
# from category counts by an EM algorithm.
#   P(R = r) = pi * C(m-1, r-1) * xi^(m-r) * (1-xi)^(r-1) + (1-pi)/m

# log of the CUB probability for every block and category: K x L x m array.
log_cub_table <- function(pi, xi, m) {
  K <- nrow(pi)
  L <- ncol(pi)
  out <- array(0, dim = c(K, L, m))
  for (r in seq_len(m)) {
    out[, , r] <- log(pi * choose(m - 1, r - 1) * xi^(m - r) * (1 - xi)^(r - 1) +
                        (1 - pi) / m)
  }
  out
}

# EM for a single CUB distribution from the counts of the m categories.
fit_cub_counts <- function(cnt, m, pi, xi, max_iter = 100, tol = 1e-8) {
  tot <- sum(cnt)
  if (tot == 0) return(c(pi, xi))
  r <- seq_len(m)
  cb <- choose(m - 1, r - 1)
  eps <- 1e-6
  for (it in seq_len(max_iter)) {
    pb <- cb * xi^(m - r) * (1 - xi)^(r - 1)
    nu <- pi * pb / (pi * pb + (1 - pi) / m)       # E step
    s <- sum(cnt * nu)
    pi_new <- min(max(s / tot, eps), 1 - eps)       # M step
    xi_new <- if (s > 0) sum(cnt * nu * (m - r)) / ((m - 1) * s) else xi
    xi_new <- min(max(xi_new, eps), 1 - eps)
    done <- abs(pi_new - pi) + abs(xi_new - xi) < tol
    pi <- pi_new
    xi <- xi_new
    if (done) break
  }
  c(pi, xi)
}
