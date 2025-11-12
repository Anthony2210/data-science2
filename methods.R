# =========================================================
# METHODES DE SELECTION DE VARIABLES
# =========================================================
library(MASS)
library(leaps)


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


# Stepwise AIC et BIC

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
