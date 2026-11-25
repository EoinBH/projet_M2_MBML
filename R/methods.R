# S3 methods for objects of class "coclust_ordinal".

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
