rm(list=ls())
library(dplyr)

# Consider a simple case 
# The treatment assignment for the RPCT data is completely at random 
# There is no unmeasured confounder

logit_fun <- function(eta_0, eta_X, x){
  return(1/(1 + exp( - (eta_0 + eta_X*x))))
}
gen_Y <- function(mean_Y, sigma_Y){
  return( rnorm(1, mean = mean_Y, sd =  sigma_Y))
}



sim_data <- function(Nt, # sample size of treatment group in RCT
                     Nc, # sample size of control group in RCT
                     NE, # sample size of EC data
                     eta_0, # intercept in R model
                     eta_X, # coef of X in R model
                     beta_0, # intercept in OM model
                     beta_X, # coef of X in OM model 
                     alpha_0, # intercept of A in OM model
                     alpha_X, # intercept of AX in OM model
                     sigma_Y, # sd of Y's distribution
                     mu_X, # mean of X's distribution
                     sigma_X # sd of X's distribution
){
  
  # RCT group
  NR <- Nt+Nc
  Index_trt_RCT <- sample(1:NR, Nt, replace = F)
  Index_crl_RCT <- setdiff(1:NR, Index_trt_RCT)
  
  # indicator for treatment group
  A <- rep(1, NR) 
  A[Index_crl_RCT]  <- 0  
  
  X_RCT <- rnorm(NR, mean = mu_X, sd = sigma_X)
  mu_Y_RCT <- beta_0 + beta_X * X_RCT + alpha_0 * A + alpha_X * A * X_RCT
  Y_RCT <- sapply(mu_Y_RCT, gen_Y, sigma_Y = sigma_Y)
  # pi_A <- Nt / NR # propensity score for RCT data
  
  # EC grouop
  X_EC <- rnorm(NE, mean = mu_X, sd = sigma_X)
  mu_Y_EC <- beta_0 + beta_X * X_EC + alpha_0 * A + alpha_X * A * X_EC
  Y_EC <- sapply(mu_Y_EC, gen_Y, sigma_Y = sigma_Y)
  
  # R_all <- logit_fun(eta_0, eta_X, X_all)
  
  data <- data.frame(
    X = c(X_RCT, X_EC),
    Y = c(Y_RCT, Y_EC),
    A = c(A, rep(0, NE)),
    R = c(rep(1, NR), rep(0, NE))
  )
    
  return(data)
    
  
}

Nc = 50
Nt = 200
NE = 1000
beta_0 = 0
beta_X = 1
alpha_0 = 0
alpha_X = 1
sigma_Y = 1
mu_X = 1
sigma_X = 1

# Need to decide the conditional variance estimation function

Est_fun <- function(data){
  
}


beta <- c(1,1,1)
Sigma <- matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1), nrow = 3)

# smallest sample size for each group in two-sided t-test
t_sample<- function(sigma1, 
                    sigma2,
                    delta,
                    alpha,
                    beta){
  n <- (sigma1 + sigma2)* (qnorm(alpha/2) + qnorm(1-beta))^2 / delta^2
  return(n)
}

# verify that paper's sample size:
t_sample(sigma1 = 2.28^2,
         sigma2 = 2.28^2,
         delta = 0.905,
         alpha = 0.05, 
         beta = 0.8)

# In our case
t_sample(sigma1 = 5,
         sigma2 = 5,
         delta = 0.8,
         alpha = 0.05, 
         beta = 0.8)


opt_sample <- function(alpha, beta, N, ){
  
}
  
  
  
