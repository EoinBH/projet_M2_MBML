# Small internal helpers: input checks and matrix utilities used by the
# SEM-Gibbs algorithm.

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
