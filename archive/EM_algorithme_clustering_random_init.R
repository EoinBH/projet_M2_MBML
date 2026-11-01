library(mvtnorm) # Pour dmvnorm au lieu de ma fonction faite maison ci-dessous

GMM_init <- function(x, k) {
  n <- nrow(x)
  d <- ncol(x)
  # Paramètres du mélange
  # Proportions de mélange aléatoires
  p <- rgamma(k, shape = 1, rate = 1)
  p <- p / sum(p)
  # Moyennes initiales sont des x choisis aléatoirement
  mu <- x[sample(seq_len(n), k, replace = FALSE), , drop = FALSE]
  # Matrices de covariance initiales
  Sigma <- array(0, dim = c(d, d, k))
  # Covariance globale
  S <- cov(x)
  # Régularisation
  S <- S + 1e-6 * diag(d)
  for (j in seq_len(k)) {
    Sigma[, , j] <- S
  }
  return(list(
    p = p,
    mu = mu,
    Sigma = Sigma
  ))
}

GMM_log_vraisemblance <- function(x, p, mu, Sigma) {
  n <- nrow(x)
  k <- length(p)
  log_vraisemblance <- 0
  for(i in seq_len(n)) {
    densite <- 0
    for(j in seq_len(k)) {
      # Avec dmvnorm il faut faire attention aux tailles des matrices!!
      #densite <- densite + p[j] * dmvnorm(x[i, ], mean = mu[j, ], sigma = Sigma[, , j])
      densite <- densite + p[j] * dmvnorm(x[i, ], mu[j, ], Sigma[, , j], log = FALSE)
    }
    log_vraisemblance <- log_vraisemblance + log(densite)
  }
  return(log_vraisemblance)
}

GMM_E_step <- function(x, p, mu, Sigma) {
  n <- nrow(x)
  k <- length(p)
  # Matrice des t_ik
  t_ik <- matrix(0, nrow = n, ncol = k)
  for(i in seq_len(n)) {
    # Calcul du numérateur pour chaque composante
    for(j in seq_len(k)) {
      #t_ik[i, j] <- p[j] * dmvnorm(x[i, ], mean = mu[j, ], sigma = Sigma[, , j])
      t_ik[i, j] <- p[j] * dmvnorm(x[i, ], mu[j, ], Sigma[, , j], log = FALSE)
    }
    # Normalisation
    t_ik[i, ] <- t_ik[i, ] / sum(t_ik[i, ])
  }
  return(t_ik)
}

regularisation_covariance <- function(Sigma, lambda = 1e-6) {
  d <- nrow(Sigma)
  return(Sigma + lambda * diag(d))
}

GMM_M_step <- function(x, t_ik) {
  n <- nrow(x)
  d <- ncol(x)
  k <- ncol(t_ik)
  # Nombre effectif d'observations dans chaque groupe
  N_k <- numeric(k)
  # Probabilités de mélange
  p <- N_k / n
  # Moyennes
  mu <- matrix(0, nrow = k, ncol = d)
  for (j in seq_len(k)) {
    mu[j, ] <- colSums(t_ik[, j] * x) / N_k[j]
  }
  # Matrices de covariance
  Sigma <- array(0, dim = c(d, d, k))
  
  for (j in seq_len(k)) {
    centered <- sweep(x, 2, mu[j, ], "-")
    Sigma[, , j] <- t(centered) %*% (centered * t_ik[, j]) / N_k[j]
    Sigma[, , j] <- regularisation_covariance(Sigma[, , j], lambda = 1e-6)
  }
  return(list(
    p = p,
    mu = mu,
    Sigma = Sigma
  ))
}

GMM_EM <- function(x, k, max_iter = 100, tolerance = 1e-6, lambda = 1e-6) {
  # Initialisation
  parametres <- GMM_init(x, k)
  p <- parametres$p
  mu <- parametres$mu
  Sigma <- parametres$Sigma
  iterations_a_tracer <- c()
  # Log-vraisemblance initiale
  log_vraisemblance_precedente <- GMM_log_vraisemblance(x, p, mu, Sigma)
  print(paste("Itération 0 - log-vraisemblance :", log_vraisemblance_precedente))
  iterations_a_tracer <- c(iterations_a_tracer, log_vraisemblance_precedente)
  # Boucle EM
  for(iter in 1:max_iter) {
    # Expectation
    t_ik <- GMM_E_step(x, p, mu, Sigma)
    # Maximisation
    nouveaux_parametres <- GMM_M_step(x, t_ik)
    p <- nouveaux_parametres$p
    mu <- nouveaux_parametres$mu
    Sigma <- nouveaux_parametres$Sigma
    # log-vraisemblance
    log_vraisemblance <- GMM_log_vraisemblance(x, p, mu, Sigma)
    print(paste("Itération", iter, "- log-vraisemblance :", log_vraisemblance))
    iterations_a_tracer <- c(iterations_a_tracer, log_vraisemblance)
    # Convergence
    difference <- abs(log_vraisemblance - log_vraisemblance_precedente)
    if(difference < tolerance) {
      print(paste("Convergence atteinte à l'itération", iter))
      break
    }
    log_vraisemblance_precedente <- log_vraisemblance
  }
  plot(iterations_a_tracer, col = "orange", pch = 18, cex = 2, main = "Log-vraisemblance par itération")
  return(list(
    p = p,
    mu = mu,
    Sigma = Sigma,
    t_ik = t_ik,
    log_vraisemblance = log_vraisemblance,
    iterations = iter
  ))
}

main <- function() {
  x <- matrix(c(
    1.0, 1.5,
    1.2, 1.8,
    1.5, 1.3,
    1.8, 1.7,
    1.3, 1.6,
    1.7, 1.4,
    
    5.0, 5.5,
    5.3, 5.8,
    5.6, 5.2,
    5.9, 5.7,
    5.2, 5.4,
    5.7, 5.6,
    
    9.0, 1.0,
    9.3, 1.4,
    9.6, 0.8,
    9.9, 1.2,
    9.2, 1.1,
    9.7, 1.0,
    
    1.0, 8.0,
    1.3, 8.4,
    1.6, 7.8,
    1.8, 8.3,
    1.2, 8.2,
    1.7, 8.0,
    
    6.5, 1.0,
    6.8, 1.4,
    7.1, 0.8,
    7.4, 1.2,
    6.7, 1.1,
    7.2, 1.3
    
  ), ncol = 2, byrow = TRUE)
  
  k <- 5
  
  resultat <- GMM_EM(x, k)
  
  print(resultat$p)
  print(resultat$mu)
  print(resultat$Sigma)
  print(resultat$t_ik)
}

main()
