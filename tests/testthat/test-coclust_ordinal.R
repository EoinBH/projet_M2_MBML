ari <- function(a, b) {
  tab <- table(a, b)
  c2 <- function(x) x * (x - 1) / 2
  sum_ij <- sum(c2(tab))
  sum_i <- sum(c2(rowSums(tab)))
  sum_j <- sum(c2(colSums(tab)))
  expected <- sum_i * sum_j / c2(sum(tab))
  (sum_ij - expected) / ((sum_i + sum_j) / 2 - expected)
}

test_that("log CUB table sums to one over categories", {
  pi <- matrix(c(0.3, 0.9), 1, 2)
  xi <- matrix(c(0.2, 0.7), 1, 2)
  lp <- log_cub_table(pi, xi, 5)
  expect_equal(apply(exp(lp), c(1, 2), sum), matrix(1, 1, 2))
})

test_that("inner EM recovers CUB parameters", {
  m <- 7
  p <- 0.8 * choose(m - 1, 0:(m - 1)) * 0.3^(m - 1:m) * 0.7^(0:(m - 1)) + 0.2 / m
  est <- fit_cub_counts(p * 1e5, m, 0.5, 0.5)
  expect_equal(est, c(0.8, 0.3), tolerance = 1e-2)
})

test_that("output has the expected structure", {
  sim <- simulate_cub_lbm(n = 30, d = 16, m = 5, seed = 1)
  fit <- coclust_ordinal(sim$X, K = 2, L = 2, n_init = 2, n_iter = 20,
                         burn_in = 8, seed = 1)
  expect_s3_class(fit, "coclust_ordinal")
  expect_equal(dim(fit$row_prob), c(30, 2))
  expect_equal(dim(fit$col_prob), c(16, 2))
  expect_equal(rowSums(fit$row_prob), rep(1, 30))
  expect_equal(rowSums(fit$col_prob), rep(1, 16))
  expect_length(fit$row_clusters, 30)
  expect_length(fit$col_clusters, 16)
  expect_true(is.finite(fit$icl))
  expect_length(fit$all_icl, 2)
  expect_equal(sum(fit$params$alpha), 1)
  expect_equal(dim(fit$params$pi), c(2, 2))
})

test_that("well separated blocks are recovered", {
  sim <- simulate_cub_lbm(n = 80, d = 40, m = 5, seed = 3)
  fit <- coclust_ordinal(sim$X, K = 2, L = 2, n_init = 5, n_iter = 50,
                         burn_in = 20, seed = 3)
  expect_gt(ari(fit$row_clusters, sim$z), 0.8)
  expect_gt(ari(fit$col_clusters, sim$w), 0.8)
})

test_that("invalid inputs are rejected", {
  X <- matrix(sample(1:5, 40, TRUE), 8, 5)
  expect_error(coclust_ordinal(X, K = 0, L = 2))
  expect_error(coclust_ordinal(X, K = 2, L = 10))
  expect_error(coclust_ordinal(X * 0.5, K = 2, L = 2))
  X2 <- X
  X2[1, 1] <- NA
  expect_error(coclust_ordinal(X2, K = 2, L = 2))
})

test_that("K = 1 and L = 1 work", {
  sim <- simulate_cub_lbm(n = 20, d = 10, m = 4, seed = 4)
  fit <- coclust_ordinal(sim$X, K = 1, L = 1, n_init = 1, n_iter = 10,
                         burn_in = 3, seed = 4)
  expect_true(all(fit$row_clusters == 1))
})
