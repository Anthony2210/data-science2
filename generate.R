# GENERATION DES DONNEES POUR LE MODELE LINEAIRE

library(mvtnorm)


# Fonction pour créer la matrice de covariance

make_SigmaX <- function(p, structure = c("indep", "ar1", "block"),
                        rho = 0, K = 5) {
  structure <- match.arg(structure)
  
  # Cas 1 : variables indépendantes
  if (structure == "indep" || rho == 0) {
    return(diag(p))
  }
  
  # Cas 2 : dépendance de type AR(1)
  if (structure == "ar1") {
    Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    return(Sigma)
  }
  
  # Cas 3 : dépendance par blocs
  if (structure == "block") {
    Sigma <- diag(p)
    block_size <- p / K
    
    for (b in seq(1, p, by = block_size)) {
      idx <- b:(b + block_size - 1)
      Sigma[idx, idx] <- rho
      diag(Sigma[idx, idx]) <- 1
    }
    
    return(Sigma)
  }
}

# Fonction pour générer le vecteur de coefficients Beta_star

make_beta <- function(p, p0) {
  beta <- rep(0, p)
  idx <- sample.int(p, p0)
  values <- runif(p0, 1, 2) * sample(c(-1, 1), p0, replace = TRUE)
  beta[idx] <- values
  beta
}

# Fonction principale : generate_lm()

generate_lm <- function(n.train, p, p0, sigma2,
                        n.test = 10 * n.train,
                        structure = c("indep", "ar1", "block"),
                        rho = 0, K = 5) {
  
  structure <- match.arg(structure)
  
  # Matrice de covariance 
  SigmaX <- make_SigmaX(p, structure, rho, K)
  
  # Vrai vecteur de coefficients
  beta_star <- make_beta(p, p0)
  
  # Génération des prédicteurs X_train et X_test
  X_train <- rmvnorm(n.train, mean = rep(0, p), sigma = SigmaX)
  X_test  <- rmvnorm(n.test,  mean = rep(0, p), sigma = SigmaX)
  
  # Génération du bruit
  eps_train <- rnorm(n.train, mean = 0, sd = sqrt(sigma2))
  eps_test  <- rnorm(n.test,  mean = 0, sd = sqrt(sigma2))
  
  # Réponses
  y_train <- as.vector(X_train %*% beta_star + eps_train)
  y_test  <- as.vector(X_test  %*% beta_star + eps_test)
  
  # Fusion dans un seul jeu complet
  y_all <- c(y_train, y_test)
  X_all <- rbind(X_train, X_test)
  
  # Indices pour train/test
  train_idx <- seq_len(n.train)
  test_idx <- n.train + seq_len(n.test)
  
  # Retourne une liste complète
  list(
    y = y_all,
    X = X_all,
    beta = beta_star,
    sigmaMatrix = SigmaX,
    sigma2 = sigma2,
    structure = structure,
    rho = rho,
    K = K,
    train = train_idx,
    test = test_idx
  )
}


generate_lm_long <- function(n.train, p, p0, sigma2,
                             n.test = 10 * n.train, rho = 0.5) {
  generate_lm(n.train, p, p0, sigma2, n.test,
              structure = "ar1", rho = rho)
}

generate_lm_block <- function(n.train, p, p0, sigma2,
                              n.test = 10 * n.train, rho = 0.5, K = 5) {
  generate_lm(n.train, p, p0, sigma2, n.test,
              structure = "block", rho = rho, K = K)
}
