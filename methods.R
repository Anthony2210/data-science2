# METHODES DE SELECTION DE VARIABLES
library(MASS)
library(leaps)
library(glmnet)

# OLS (moindres carrés ordinaires)

getOLS <- function(X, y) {
  fit <- lm.fit(x = X, y = y)
  as.vector(fit$coefficients)
}


# Oracle (connait le vrai support S*)

getOracle <- function(X, y, beta_star) {
  p <- ncol(X)
  beta_hat <- rep(0, p)
  idx <- which(beta_star != 0)
  fit <- lm.fit(x = X[, idx, drop = FALSE], y = y)
  beta_hat[idx] <- fit$coefficients
  beta_hat
}

# Best Subset (meilleure combinaison possible selon AIC/BIC)
getBestSubset <- function(X, y, criterion = "bic") {
  p <- ncol(X)
  df <- as.data.frame(X)
  colnames(df) <- paste0("x", 1:p)
  df$y <- y
  
  best <- regsubsets(
    y ~ .,
    data = df,
    method = "exhaustive",
    nbest = 1
  )
  
  summ <- summary(best)
  
  if (criterion == "bic") {
    idx <- which.min(summ$bic)
  } else if (criterion == "aic") {
    idx <- which.min(summ$aic)
  } else if (criterion == "cp") {
    idx <- which.min(summ$cp)
  } else {
    stop("Critère inconnu : choisir 'bic', 'aic' ou 'cp'")
  }
  
  vars <- names(which(summ$which[idx, ]))  # garde l’intercept
  vars <- vars[vars != "(Intercept)"]       # on enlève l’intercept pour indexer
  
  beta <- rep(0, p)
  if (length(vars) > 0) {
    idx_vars <- as.numeric(gsub("x", "", vars))
    fit <- lm(y ~ 0 + ., data = df[, c(vars, "y")])  # régression sans intercept
    beta[idx_vars] <- coef(fit)
  }
  beta
}

# Sélection pas à pas (Stepwise AIC)

getStepwiseAIC <- function(X, y) {
  p <- ncol(X)
  df <- as.data.frame(X)
  colnames(df) <- paste0("x", 1:p)
  df$y <- y
  
  full <- lm(as.formula(paste0("y ~ ", paste(colnames(df)[1:p], collapse = " + "), " - 1")), data = df)
  null <- lm(y ~ -1, data = df)
  
  fit <- stepAIC(null, scope = list(lower = null, upper = full),
                 direction = "both", k = 2, trace = FALSE)
  
  beta <- rep(0, p)
  coefs <- coef(fit)
  vars <- names(coefs)
  vars <- vars[vars != "(Intercept)"]
  if (length(vars) > 0) {
    idx <- as.numeric(gsub("x", "", vars))
    beta[idx] <- coefs[vars]
  }
  beta
  
}

# Sélection pas à pas (Stepwise BIC)

getStepwiseBIC <- function(X, y) {
  p <- ncol(X)
  df <- as.data.frame(X)
  colnames(df) <- paste0("x", 1:p)
  df$y <- y
  
  full <- lm(as.formula(paste0("y ~ ", paste(colnames(df)[1:p], collapse = " + "), " - 1")), data = df)
  null <- lm(y ~ -1, data = df)
  
  fit <- stepAIC(null, scope = list(lower = null, upper = full),
                 direction = "both", k = log(nrow(X)), trace = FALSE)
  
  beta <- rep(0, p)
  coefs <- coef(fit)
  vars <- names(coefs)
  vars <- vars[vars != "(Intercept)"]
  if (length(vars) > 0) {
    idx <- as.numeric(gsub("x", "", vars))
    beta[idx] <- coefs[vars]
  }
  beta
}

# Régression LASSO

getLasso <- function(X, y) {
  fit <- cv.glmnet(X, y, alpha = 1, standardize = TRUE)
  as.vector(coef(fit, s = "lambda.min"))[-1]
}
