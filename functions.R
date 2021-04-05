
library(rlang)
library(tidyverse)
library(tidymodels)
library(parsnip)
library(dplyr)

fit_logistic_lasso <- function(x, y, lambda, beta0 = NULL, eps = 0.0001, iter_max = 100) {
  ## fit_logistic_lasso is a function to generate logistic lasso regression
  ## model for the further analysis. It will return a large list of the
  ## regression.
  ##
  ## Input: 
  ##
  ## - x: A matrix of predictors that we are wish to analyze by 
  ##      logistic lasso regression.
  ## 
  ## - y: A vector of the selected data. 
  ##
  ## - lambda: a value that indicates the penalty in lasso regression. 
  ##
  ## - beta0: the initial guess of the parameter. 
  ##
  ## - eps: parameter for stopping critereon
  ##
  ## - iter_max: maximum number of iterations. 
  ##
  ## Output:
  ##
  ## - A large list of the lasso regression, containing the members
  ##   intercept, beta and lambda
  ##
  ## Example: 
  ##    ret <- fit_logistic_lasso(x,y,lambda,beta0=NULL, eps=0.0001, iter_max=100)
  ##    ret$fit$fit$fit$beta    
  ##    ret$fit$fit$fit$intercept
  
  n <- dim(x)[1]
  m <- rep(1,n)
  x<-cbind(m,x)
  colnames(x)[1] = c('(intercept)')
  p <- dim(x)[2] 
  if (is.null(beta0)) {
    beta0 <- rep(0,p)
  } 
  fct_levels <- levels(y)  
  y <- as.numeric(y) - 1
  beta <- beta0
  x_beta <- (x %*% beta) %>% as.numeric 
  q <- 1/(1 + exp(-x_beta))  
  
  for(iter in 1:iter_max) {
    w <- q * (1 - q)  
    z <- x_beta + (y - q)/w  
    for(j in 1:p){
      rj = z - x[,-j]%*%beta[-j] 
      beta[j] = sign(sum(w*rj*x[,j])) * max(abs(sum(w*rj*x[,j]))-lambda, 0) / sum(w*x[,j]^2)
    }
    
    
    x_beta <- (x %*% beta) %>% as.numeric
    q <- 1/(1 + exp(-x_beta))
    grad <- t(x) %*% (y - q) / n 
    
    if (sqrt(sum(grad^2)) < eps) {
      return(
        list(beta = beta[-1], intercept=beta[1], 
             fct_levels = fct_levels, lambda=lambda, 
             iter = iter, converged = TRUE, error = sqrt(sum(grad^2)))
      )
    }
    
  }
  
  names(beta) <- colnames(x)
  
  warning(paste("Method did not converge in", iter_max, "iterations", sep = " "))
  return(
    list(beta = beta[-1],  intercept=beta[1],
         fct_levels = fct_levels, lambda=lambda, 
         iter = iter, converged = FALSE, error = sqrt(sum(grad^2)))
  )
}


predict_logistic_lasso <- function(fit, new_x) {
  ## predict_logistic_lasso is a function to generate logistic lasso regression
  ## predictions for the further analysis. 
  ##
  ## Input: 
  ##
  ## - fit: Output from fit_logistic_lasso()
  ## 
  ## - new_x: Data to predict at
  ##
  ##
  ## Output:
  ##
  ## - A large list containing the intercept and beta. 
  ##
  ## Example: 
  ##    ret <- fit_logistic_lasso(x,y,lambda,beta0=NULL, eps=0.0001, iter_max=100)
  ##    predict(ret, new_data = test) %>% bind_cols(test %>% select(y)) %>%
  ##            conf_mat(truth = y, estimate = .pred_class)   
  ##    
  ##  
  
  numeric_pred <- ((new_x %*% fit$beta + fit$intercept * rep(1,dim(new_x)[1])) >= 0) %>% as.numeric
  return( fit$fct_levels[numeric_pred + 1] %>% factor)
}
