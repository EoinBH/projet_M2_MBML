#' Co-clustering of ordinal data with a CUB latent block model
#'
#' Fits a latent block model (LBM) whose \eqn{K \times L} blocks follow a CUB
#' (Combination of Uniform and shifted Binomial) distribution. Inference is
#' done by a SEM-Gibbs algorithm. The algorithm is run `n_init` times from
#' random initializations and the solution with the best ICL is returned.
#'
#' @section Model:
#' Rows are partitioned into `K` clusters (proportions `alpha`) and columns
#' into `L` clusters (proportions `beta`). Given the partitions, the cells of
#' block \eqn{(k,\ell)} are i.i.d. CUB:
#' \deqn{P(X = r) = \pi_{k\ell} \binom{m-1}{r-1} \xi_{k\ell}^{m-r}
#' (1-\xi_{k\ell})^{r-1} + (1-\pi_{k\ell})/m, \quad r = 1, \dots, m.}
#' `pi` is the weight of the "feeling" component (1 - `pi` is the weight of
#' uncertainty) and `xi` drives the feeling (the larger `xi`, the lower the
#' ordinal values).
#'
#' @section Algorithm:
#' *SE step*: `n_gibbs` Gibbs sweeps sample the row partition given the
#' column partition, then the column partition given the row partition.
#' *M step*: `alpha` and `beta` are cluster frequencies; `pi` and `xi` of each
#' block are obtained by an inner EM on the block's category counts.
#' After `burn_in` iterations the sampled parameters are averaged. The
#' membership probabilities are then the frequencies of the Gibbs chain
#' (run at the averaged parameters) and the final partitions are its mode.
#'
#' @section ICL:
#' \deqn{ICL = \log p(x, \hat z, \hat w; \hat\theta) - \frac{K-1}{2}\log n -
#' \frac{L-1}{2}\log d - \frac{2KL}{2}\log(nd).}
#' The larger, the better.
#'
#' @param X Matrix or data frame (n rows, d columns) of ordinal data coded
#'   with integers `1, ..., m`. No missing values.
#' @param K Number of row clusters.
#' @param L Number of column clusters.
#' @param n_init Number of random initializations (default 20).
#' @param n_iter Number of SEM iterations per initialization (default 100).
#' @param burn_in Number of burn-in iterations (default 40).
#' @param n_gibbs Number of Gibbs sweeps per SE step (default 1).
#' @param m Number of ordinal categories. Default: `max(X)`.
#' @param init Initialization of the partitions: `"kmeans"` (random k-means,
#'   default) or `"random"` (uniform draws).
#' @param seed Optional integer seed for reproducibility.
#' @param verbose Print the ICL of each initialization?
#'
#' @return An object of class `"coclust_ordinal"`, a list with
#' \describe{
#'   \item{row_prob}{n x K matrix of row membership probabilities.}
#'   \item{col_prob}{d x L matrix of column membership probabilities.}
#'   \item{row_clusters}{Row assignments (integers in `1:K`).}
#'   \item{col_clusters}{Column assignments (integers in `1:L`).}
#'   \item{icl}{ICL value of the returned solution.}
#'   \item{loglik}{Complete-data log-likelihood at the returned solution.}
#'   \item{params}{List with `alpha` (length K), `beta` (length L), and the
#'     K x L matrices `pi` and `xi`.}
#'   \item{all_icl}{ICL of each of the `n_init` initializations.}
#'   \item{K, L, m, n, d}{Dimensions of the model and of the data.}
#' }
#'
#' @examples
#' sim <- simulate_cub_lbm(n = 60, d = 30, m = 5, seed = 1)
#' fit <- coclust_ordinal(sim$X, K = 2, L = 2, n_init = 3,
#'                        n_iter = 40, burn_in = 15, seed = 1)
#' fit
#' table(fit$row_clusters, sim$z)
#' fit$params$pi
#' @export
coclust_ordinal <- function(X, K, L, n_init = 20, n_iter = 100, burn_in = 40,
                            n_gibbs = 1, m = NULL,
                            init = c("kmeans", "random"),
                            seed = NULL, verbose = FALSE) {
  init <- match.arg(init)
  cl <- match.call()
  chk <- check_ordinal_data(X, m)
  X <- chk$X
  m <- chk$m
  n <- nrow(X)
  d <- ncol(X)
  K <- check_count(K, "K", n)
  L <- check_count(L, "L", d)
  n_init <- check_count(n_init, "n_init", Inf)
  n_iter <- check_count(n_iter, "n_iter", Inf)
  n_gibbs <- check_count(n_gibbs, "n_gibbs", Inf)
  if (!is.numeric(burn_in) || burn_in < 0 || burn_in >= n_iter) {
    stop("`burn_in` must satisfy 0 <= burn_in < n_iter.")
  }
  if (!is.null(seed)) set.seed(seed)

  Ilist <- lapply(seq_len(m), function(r) (X == r) * 1)

  best <- NULL
  all_icl <- numeric(n_init)
  for (it in seq_len(n_init)) {
    fit <- sem_gibbs_once(Ilist, n, d, m, K, L, n_iter, burn_in, n_gibbs, init, X)
    all_icl[it] <- fit$icl
    if (verbose) message(sprintf("init %d/%d: ICL = %.2f", it, n_init, fit$icl))
    if (is.null(best) || fit$icl > best$icl) best <- fit
  }
  best$all_icl <- all_icl
  best$K <- K
  best$L <- L
  best$m <- m
  best$n <- n
  best$d <- d
  best$call <- cl
  class(best) <- "coclust_ordinal"
  best
}

#' Select the numbers of row and column clusters with ICL
#'
#' Fits [coclust_ordinal()] on every pair `(K, L)` of the grid and keeps the
#' pair with the largest ICL.
#'
#' @param X Ordinal data matrix, see [coclust_ordinal()].
#' @param K_range,L_range Integer vectors of candidate numbers of clusters.
#' @param ... Further arguments passed to [coclust_ordinal()]
#'   (e.g. `n_init`, `n_iter`, `burn_in`, `seed`).
#'
#' @return A list with `table` (data frame with columns `K`, `L`, `icl`),
#'   `best` (the best `"coclust_ordinal"` fit) and `fits` (all fits).
#' @examples
#' sim <- simulate_cub_lbm(n = 50, d = 24, m = 5, seed = 2)
#' sel <- coclust_ordinal_select(sim$X, K_range = 1:3, L_range = 1:3,
#'                               n_init = 2, n_iter = 30, burn_in = 10, seed = 2)
#' sel$table
#' @export
coclust_ordinal_select <- function(X, K_range = 1:4, L_range = 1:4, ...) {
  grid <- expand.grid(K = K_range, L = L_range)
  fits <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    fits[[i]] <- coclust_ordinal(X, K = grid$K[i], L = grid$L[i], ...)
  }
  grid$icl <- vapply(fits, function(f) f$icl, numeric(1))
  list(table = grid, best = fits[[which.max(grid$icl)]], fits = fits)
}

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

#' @export
print.coclust_ordinal <- function(x, ...) {
  cat("Co-clustering of ordinal data (CUB latent block model)\n")
  cat(sprintf("  %d rows x %d columns, %d categories\n", x$n, x$d, x$m))
  cat(sprintf("  K = %d row clusters, L = %d column clusters\n", x$K, x$L))
  cat(sprintf("  ICL = %.2f\n", x$icl))
  cat("  Row cluster sizes:   ", tabulate(x$row_clusters, x$K), "\n")
  cat("  Column cluster sizes:", tabulate(x$col_clusters, x$L), "\n")
  invisible(x)
}

#' @export
summary.coclust_ordinal <- function(object, ...) {
  print(object)
  cat("\nRow proportions (alpha):\n")
  print(round(object$params$alpha, 3))
  cat("Column proportions (beta):\n")
  print(round(object$params$beta, 3))
  cat("Block parameter pi (weight of feeling):\n")
  print(round(object$params$pi, 3))
  cat("Block parameter xi (feeling is 1 - xi):\n")
  print(round(object$params$xi, 3))
  invisible(object)
}

# ---------------------------------------------------------------------------
# Internal functions
# ---------------------------------------------------------------------------

check_ordinal_data <- function(X, m) {
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("`X` must be numeric (integers 1, ..., m).")
  if (anyNA(X)) stop("`X` must not contain missing values.")
  if (any(X != round(X))) stop("`X` must contain integer values.")
  if (min(X) < 1) stop("`X` must be coded with integers starting at 1.")
  if (is.null(m)) m <- max(X)
  if (max(X) > m) stop("`X` contains values larger than `m`.")
  if (m < 2) stop("At least 2 ordinal categories are needed.")
  storage.mode(X) <- "integer"
  list(X = X, m = as.integer(m))
}

check_count <- function(v, name, upper) {
  if (length(v) != 1 || !is.numeric(v) || v < 1 || v != round(v) || v > upper) {
    stop(sprintf("`%s` must be a positive integer%s.", name,
                 if (is.finite(upper)) sprintf(" not larger than %d", upper) else ""))
  }
  as.integer(v)
}

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

indicator <- function(z, K) {
  Z <- matrix(0, length(z), K)
  Z[cbind(seq_along(z), z)] <- 1
  Z
}

# Number of cells of each category in each block: K x L x m array.
block_counts <- function(Ilist, Z, W) {
  K <- ncol(Z)
  L <- ncol(W)
  m <- length(Ilist)
  N <- array(0, dim = c(K, L, m))
  for (r in seq_len(m)) N[, , r] <- crossprod(Z, Ilist[[r]] %*% W)
  N
}

# Draw one category per row of a matrix of log-scores (Gumbel-max trick).
draw_categorical <- function(logscore) {
  g <- -log(-log(stats::runif(length(logscore))))
  max.col(logscore + matrix(g, nrow(logscore)), ties.method = "first")
}

# Make sure no cluster is empty by moving a random member of the largest one.
fix_empty <- function(z, K) {
  repeat {
    cnt <- tabulate(z, K)
    empty <- which(cnt == 0)
    if (length(empty) == 0) break
    donor <- which.max(cnt)
    members <- which(z == donor)
    z[members[sample.int(length(members), 1)]] <- empty[1]
  }
  z
}

# One Gibbs sweep: rows given columns, then columns given rows.
gibbs_sweep <- function(Ilist, z, w, logalpha, logbeta, logP, K, L, n_gibbs = 1) {
  n <- length(z)
  d <- length(w)
  m <- length(Ilist)
  for (g in seq_len(n_gibbs)) {
    W <- indicator(w, L)
    score <- matrix(logalpha, n, K, byrow = TRUE)
    for (r in seq_len(m)) {
      score <- score + (Ilist[[r]] %*% W) %*% t(matrix(logP[, , r], K, L))
    }
    z <- fix_empty(draw_categorical(score), K)

    Z <- indicator(z, K)
    score <- matrix(logbeta, d, L, byrow = TRUE)
    for (r in seq_len(m)) {
      score <- score + crossprod(Ilist[[r]], Z) %*% matrix(logP[, , r], K, L)
    }
    w <- fix_empty(draw_categorical(score), L)
  }
  list(z = z, w = w)
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

# M step: alpha, beta and (pi, xi) of every block.
m_step <- function(Ilist, z, w, K, L, pi, xi) {
  m <- length(Ilist)
  N <- block_counts(Ilist, indicator(z, K), indicator(w, L))
  for (k in seq_len(K)) {
    for (l in seq_len(L)) {
      est <- fit_cub_counts(N[k, l, ], m, pi[k, l], xi[k, l])
      pi[k, l] <- est[1]
      xi[k, l] <- est[2]
    }
  }
  list(alpha = tabulate(z, K) / length(z), beta = tabulate(w, L) / length(w),
       pi = pi, xi = xi)
}

initial_partition <- function(M, K, init) {
  n <- nrow(M)
  if (K == 1) return(rep(1L, n))
  z <- NULL
  if (init == "kmeans") {
    z <- tryCatch(suppressWarnings(stats::kmeans(M, centers = K, nstart = 1)$cluster),
                  error = function(e) NULL)
  }
  if (is.null(z)) z <- sample.int(K, n, replace = TRUE)
  fix_empty(as.integer(z), K)
}

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

# One run of SEM-Gibbs from one random initialization.
sem_gibbs_once <- function(Ilist, n, d, m, K, L, n_iter, burn_in, n_gibbs, init, X) {
  z <- initial_partition(X, K, init)
  w <- initial_partition(t(X), L, init)
  pi <- matrix(stats::runif(K * L, 0.3, 0.9), K, L)
  xi <- matrix(stats::runif(K * L, 0.1, 0.9), K, L)
  theta <- m_step(Ilist, z, w, K, L, pi, xi)

  keep <- n_iter - burn_in
  acc <- list(alpha = 0, beta = 0, pi = 0, xi = 0)
  for (it in seq_len(n_iter)) {
    logP <- log_cub_table(theta$pi, theta$xi, m)
    s <- gibbs_sweep(Ilist, z, w, log(theta$alpha), log(theta$beta), logP,
                     K, L, n_gibbs)                                  # SE step
    z <- s$z
    w <- s$w
    theta <- m_step(Ilist, z, w, K, L, theta$pi, theta$xi)           # M step
    if (it > burn_in) {
      for (nm in names(acc)) acc[[nm]] <- acc[[nm]] + theta[[nm]] / keep
    }
  }
  theta_hat <- acc

  # Gibbs chain at fixed theta_hat: membership probabilities and partitions.
  logP <- log_cub_table(theta_hat$pi, theta_hat$xi, m)
  rc <- matrix(0, n, K)
  cc <- matrix(0, d, L)
  for (it in seq_len(n_iter)) {
    s <- gibbs_sweep(Ilist, z, w, log(theta_hat$alpha), log(theta_hat$beta),
                     logP, K, L, n_gibbs)
    z <- s$z
    w <- s$w
    if (it > burn_in) {
      rc <- rc + indicator(z, K)
      cc <- cc + indicator(w, L)
    }
  }
  row_prob <- rc / keep
  col_prob <- cc / keep
  z_hat <- max.col(row_prob, ties.method = "first")
  w_hat <- max.col(col_prob, ties.method = "first")

  ll <- complete_loglik(Ilist, z_hat, w_hat, theta_hat, K, L)
  list(row_prob = row_prob, col_prob = col_prob,
       row_clusters = z_hat, col_clusters = w_hat,
       icl = icl_value(ll, n, d, K, L), loglik = ll, params = theta_hat)
}
