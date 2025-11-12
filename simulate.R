library(tidyverse)

getOneSimu <- function(i, structure = "indep", rho = 0.6) {
  data <- generate_lm(n.train = 100, p = 20, p0 = 5, sigma2 = 1,
                      structure = structure, rho = rho)
  
  X_train <- data$X[data$train, ]
  y_train <- data$y[data$train]
  X_test <- data$X[data$test, ]
  y_test <- data$y[data$test]
  
  # Estimations
  beta_ols     <- getOLS(X_train, y_train)
  beta_best    <- getBestSubset(X_train, y_train, criterion = "bic")  
  beta_aic     <- getStepwiseAIC(X_train, y_train)
  beta_bic     <- getStepwiseBIC(X_train, y_train)
  beta_lasso   <- getLasso(X_train, y_train)
  beta_oracle  <- getOracle(X_train, y_train, data$beta)
  
  # Performances
  res <- rbind(
    OLS          = perf(X_test, y_test, beta_ols, data$beta),
    BestSubset   = perf(X_test, y_test, beta_best, data$beta),
    StepwiseAIC  = perf(X_test, y_test, beta_aic, data$beta),
    StepwiseBIC  = perf(X_test, y_test, beta_bic, data$beta),
    LASSO        = perf(X_test, y_test, beta_lasso, data$beta),
    Oracle       = perf(X_test, y_test, beta_oracle, data$beta)
  ) %>%
    as.data.frame() %>%
    rownames_to_column("method") %>%
    mutate(simu = i,
           structure = structure,
           rho = rho)
  
  res
}
