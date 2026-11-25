# The SEM-Gibbs algorithm for the CUB latent block model.
# sem_gibbs_once() is the main loop; the other functions are its steps.

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
