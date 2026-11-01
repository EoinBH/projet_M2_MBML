library(mvtnorm) # Pour dmvnorm au lieu de ma fonction faite maison

log_sum_exp <- function(x) {
  maximum <- max(x)
  return(maximum + log(sum(exp(x - maximum))))
}

GMM_init <- function(x, k, lambda = 1e-6) {
  n <- nrow(x)
  d <- ncol(x)
  
  if (k <= 0 || k > n) {
    stop("Le nombre de groupes k doit être compris entre 1 et n.")
  }
  
  # Proportions de mélange aléatoires
  p <- rgamma(k, shape = 1, rate = 1)
  p <- p / sum(p)
  
  # Moyennes initiales : k observations choisies aléatoirement
  mu <- x[sample(seq_len(n), k, replace = FALSE), , drop = FALSE]
  
  # Matrices de covariance initiales
  Sigma <- array(0, dim = c(d, d, k))
  
  # Covariance globale
  S <- cov(x)
  
  # Régularisation de la covariance globale
  S <- S + lambda * diag(d)
  
  for (j in seq_len(k)) {
    Sigma[, , j] <- S
  }
  
  return(list(p = p, mu = mu, Sigma = Sigma))
}

GMM_log_vraisemblance <- function(x, p, mu, Sigma) {
  n <- nrow(x)
  k <- length(p)
  log_vraisemblance <- 0
  
  for (i in seq_len(n)) {
    log_prob <- numeric(k)
    
    for (j in seq_len(k)) {
      # Calcul sur l'échelle logarithmique pour éviter les problèmes numériques
      log_prob[j] <- log(p[j]) + dmvnorm(x[i, ], mean = mu[j, ], sigma = Sigma[, , j], log = TRUE)
    }
    
    log_vraisemblance <- log_vraisemblance + log_sum_exp(log_prob)
  }
  
  return(log_vraisemblance)
}

GMM_E_step <- function(x, p, mu, Sigma) {
  n <- nrow(x)
  k <- length(p)
  
  # Matrice des t_ik
  t_ik <- matrix(0, nrow = n, ncol = k)
  
  for (i in seq_len(n)) {
    # Calcul des probabilités sur l'échelle logarithmique
    log_prob <- numeric(k)
    
    for (j in seq_len(k)) {
      log_prob[j] <- log(p[j]) + dmvnorm(x[i, ], mean = mu[j, ], sigma = Sigma[, , j], log = TRUE)
    }
    
    # Normalisation avec log-sum-exp
    log_denominator <- log_sum_exp(log_prob)
    t_ik[i, ] <- exp(log_prob - log_denominator)
  }
  
  return(t_ik)
}

regularisation_covariance <- function(Sigma, lambda = 1e-6) {
  d <- nrow(Sigma)
  return(Sigma + lambda * diag(d))
}

GMM_M_step <- function(x, t_ik, lambda = 1e-6) {
  n <- nrow(x)
  d <- ncol(x)
  k <- ncol(t_ik)
  
  # Nombre effectif d'observations dans chaque groupe
  N_k <- colSums(t_ik)
  
  if (any(N_k < 1e-8)) {
    warning("Un ou plusieurs groupes ont presque aucune observation.")
  }
  
  # Probabilités de mélange
  p <- N_k / n
  
  # Moyennes
  mu <- matrix(0, nrow = k, ncol = d)
  
  for (j in seq_len(k)) {
    if (N_k[j] < 1e-8) {
      # Réinitialisation de la moyenne d'un groupe presque vide
      mu[j, ] <- x[sample(seq_len(n), 1), , drop = TRUE]
    } else {
      mu[j, ] <- colSums(t_ik[, j] * x) / N_k[j]
    }
  }
  
  # Matrices de covariance
  Sigma <- array(0, dim = c(d, d, k))
  
  for (j in seq_len(k)) {
    if (N_k[j] < 1e-8) {
      # Utilisation de la covariance globale pour un groupe vide
      Sigma[, , j] <- cov(x) + lambda * diag(d)
    } else {
      # Centrage des observations autour de mu_j
      centered <- sweep(x, 2, mu[j, ], "-")
      
      # Calcul de la covariance pondérée
      Sigma[, , j] <- t(centered) %*% (centered * t_ik[, j]) / N_k[j]
      
      # Régularisation de la covariance
      Sigma[, , j] <- regularisation_covariance(Sigma[, , j], lambda = lambda)
    }
  }
  
  return(list(p = p, mu = mu, Sigma = Sigma))
}

GMM_EM <- function(x, k, max_iter = 100, tolerance = 1e-6, lambda = 1e-6) {
  # Vérification des données
  if (!is.matrix(x)) {
    stop("x doit être une matrice.")
  }
  
  if (any(!is.finite(x))) {
    stop("Les données contiennent des valeurs non finies.")
  }
  
  # Initialisation
  parametres <- GMM_init(x, k, lambda = lambda)
  p <- parametres$p
  mu <- parametres$mu
  Sigma <- parametres$Sigma
  iterations_a_tracer <- c()
  
  # Log-vraisemblance initiale
  log_vraisemblance_precedente <- GMM_log_vraisemblance(x, p, mu, Sigma)
  
  print(paste("Itération 0 - log-vraisemblance :", log_vraisemblance_precedente))
  iterations_a_tracer <- c(iterations_a_tracer, log_vraisemblance_precedente)
  
  # Boucle EM
  for (iter in seq_len(max_iter)) {
    # Expectation
    t_ik <- GMM_E_step(x, p, mu, Sigma)
    
    # Maximisation
    nouveaux_parametres <- GMM_M_step(x, t_ik, lambda = lambda)
    p <- nouveaux_parametres$p
    mu <- nouveaux_parametres$mu
    Sigma <- nouveaux_parametres$Sigma
    
    # Log-vraisemblance
    log_vraisemblance <- GMM_log_vraisemblance(x, p, mu, Sigma)
    
    print(paste("Itération", iter, "- log-vraisemblance :", log_vraisemblance))
    iterations_a_tracer <- c(iterations_a_tracer, log_vraisemblance)
    
    # Vérification numérique
    if (!is.finite(log_vraisemblance)) {
      stop(paste("La log-vraisemblance n'est pas finie à l'itération", iter))
    }
    
    # Vérification de la monotonie de la log-vraisemblance
    difference <- log_vraisemblance - log_vraisemblance_precedente
    
    if (difference < -1e-8) {
      warning(paste("La log-vraisemblance a diminué à l'itération", iter))
    }
    
    # Convergence
    difference_absolue <- abs(difference)
    
    if (difference_absolue < tolerance) {
      print(paste("Convergence atteinte à l'itération", iter))
      break
    }
    
    log_vraisemblance_precedente <- log_vraisemblance
  }
  
  # Graphique de convergence
  plot(iterations_a_tracer, col = "orange", pch = 18, cex = 1.5, type = "b", main = "Log-vraisemblance par itération", xlab = "Itération", ylab = "Log-vraisemblance")
  
  return(list(p = p, mu = mu, Sigma = Sigma, t_ik = t_ik, log_vraisemblance = log_vraisemblance, iterations = iter))
}

main <- function() {
  # Chargement des données Iris
  data(iris)
  
  # Sélection des quatre variables numériques
  # Species n'est pas utilisée par l'algorithme EM
  x <- as.matrix(iris[, 1:4])
  
  # Standardisation des données
  x <- scale(x)
  
  # Nombre de groupes
  k <- 3
  
  # Exécution de l'algorithme EM
  resultat <- GMM_EM(x = x, k = k, max_iter = 100, tolerance = 1e-6, lambda = 1e-6)
  
  # Affichage des paramètres estimés
  print("Proportions de mélange :")
  print(resultat$p)
  
  print("Moyennes :")
  print(resultat$mu)
  
  print("Matrices de covariance :")
  print(resultat$Sigma)
  
  # Attribution de chaque observation au groupe ayant la plus grande probabilité a posteriori
  z <- max.col(resultat$t_ik)
  
  print("Groupes attribués :")
  print(z)
  
  # Taille des groupes
  print("Taille des groupes :")
  print(table(z))
  
  # Comparaison avec les espèces Iris
  # Les noms des groupes sont arbitraires
  print("Table de contingence entre groupes et espèces :")
  print(table(Groupe = z, Espece = iris$Species))
  
  return(resultat)
}

main()