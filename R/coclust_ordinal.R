#' Co-clustering of ordinal data with a CUB latent block model
#'
#' Fits a latent block model (LBM) whose \eqn{K \times L}{K x L} blocks follow a CUB
#' (Combination of Uniform and shifted Binomial) distribution. Inference is
#' done by a SEM-Gibbs algorithm. The algorithm is run `n_init` times from
#' random initializations and the solution with the best ICL is returned.
#'
#' @section Model:
#' Rows are partitioned into `K` clusters (proportions `alpha`) and columns
#' into `L` clusters (proportions `beta`). Given the partitions, the cells of
#' block \eqn{(k,\ell)}{(k,l)} are i.i.d. CUB:
#' \deqn{P(X = r) = \pi_{k\ell} \binom{m-1}{r-1} \xi_{k\ell}^{m-r}
#' (1-\xi_{k\ell})^{r-1} + (1-\pi_{k\ell})/m, \quad r = 1, \dots, m.}{
#' P(X = r) = pi_kl * C(m-1, r-1) * xi_kl^(m-r) * (1-xi_kl)^(r-1)
#' + (1-pi_kl)/m,   r = 1, ..., m.}
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
#' \frac{L-1}{2}\log d - KL\log(nd).}{
#' ICL = log p(x, z_hat, w_hat; theta_hat) - (K-1)/2 log(n)
#' - (L-1)/2 log(d) - K*L log(n*d).}
#' Each of the \eqn{KL}{K*L} blocks has 2 free parameters (`pi` and `xi`)
#' and its penalty uses the \eqn{nd}{n*d} cells.
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
