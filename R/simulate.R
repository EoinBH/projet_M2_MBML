# Simulation of ordinal data from a CUB latent block model (useful for
# testing the package: the true partitions and parameters are returned).

#' Simulate ordinal data from a CUB latent block model
#'
#' @param n,d Number of rows and columns.
#' @param m Number of ordinal categories.
#' @param alpha,beta Row and column cluster proportions.
#' @param pi,xi K x L matrices of CUB parameters. Defaults give well
#'   separated blocks for K = L = 2.
#' @param seed Optional seed.
#'
#' @return A list with the data matrix `X`, the true partitions `z` (rows) and
#'   `w` (columns) and the true parameters.
#' @examples
#' sim <- simulate_cub_lbm(n = 40, d = 20, m = 5, seed = 1)
#' dim(sim$X)
#' @export
simulate_cub_lbm <- function(n = 100, d = 50, m = 5,
                             alpha = c(0.5, 0.5), beta = c(0.5, 0.5),
                             pi = matrix(0.9, length(alpha), length(beta)),
                             xi = matrix(c(0.2, 0.8, 0.8, 0.2), 2, 2),
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  K <- length(alpha)
  L <- length(beta)
  pi <- matrix(pi, K, L)
  xi <- matrix(xi, K, L)
  z <- sample.int(K, n, replace = TRUE, prob = alpha)
  w <- sample.int(L, d, replace = TRUE, prob = beta)
  X <- matrix(0L, n, d)
  for (i in seq_len(n)) {
    for (h in seq_len(d)) {
      k <- z[i]
      l <- w[h]
      if (stats::runif(1) < pi[k, l]) {
        X[i, h] <- 1L + stats::rbinom(1, m - 1, 1 - xi[k, l])
      } else {
        X[i, h] <- sample.int(m, 1)
      }
    }
  }
  list(X = X, z = z, w = w,
       params = list(alpha = alpha, beta = beta, pi = pi, xi = xi), m = m)
}
