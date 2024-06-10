rm(list=ls())
library(nloptr)

# Consider a simple case 
# The treatment assignment for the RPCT data is completely at random 
# There is no unmeasured confounders

# Generate outcome Y based on linear regression
f_Y_gen <- function(x,
                    beta0,
                    beta1, 
                    sdY){
  mu <- beta0 + beta1 * x 
  return(rnorm(1, mean = mu, sd = sdY))
}

# Generate EC data
f_EC_out <- function(seed_EC = 123,
                     N_EC = 1000,
                     muX_EC = 1, 
                     sdX_EC = 1,
                     beta0 = 1,
                     beta1 = 2,
                     sdY = 1){
  set.seed(seed_EC)
  # generate the pre-treatment covariates
  X_EC <- rnorm(N_EC, mean = muX_EC, sd = sdX_EC)
  
  Y_EC <- mapply(f_Y_gen,
                 x = X_EC, 
                 beta0, beta1, sdY)
  return(data.frame(X_EC, Y_EC))
}

data_EC <- f_EC_out(N_EC = 1000)

# Estimate variance of IF-based estimator using EC data
Est_var_simple <- function(
    N_RCT, # sample size of RCT, need to be decided
    data_EC, # available EC data
    pi = 0.5, # given propensity socre
    var_fun, # given function to estimate variance
    r = 1, # the ratio of variance in EC and variance in control of RCT
    gamma1 = 1, # the ratio of two kappa
    d_X = 1 # the ratio of f(X|R=1) with f(X|R=0)
){
  # Step 1: estimate V(Y|X,R=0) using data_EC
  var_X_EC <- var_fun(data_EC)
  
  # Step 2: estimate V(Y|X,R=1,A=0)
  var_X_RCT_0 <- r * var_X_EC
  
  # Step 3: estimate kappa 0 and kappa 1
  kappa_0 <- mean(var_X_RCT_0)
  kappa_1 <- gamma1 * kappa_0
  
  # Step 4: estimate V(Y|R=1,A=1), V(Y|R=1, A=0)
  var_EC <- var(data_EC$Y_EC)
  var_RCT_0 <- var_EC
  var_RCT_1 <- var_EC
  
  # Step 5: calculate d(X)
    
  # Step 6: plug in to get nu^2
  N <- N_EC + N_RCT
  r_R <- N_RCT / N_EC
  term1 <- 1 / N_RCT / pi * mean(var_X_EC)
  term2 <- 1 / N_RCT * (var_EC -  mean(var_X_RCT_0))
  term3 <- 1 / N_RCT * mean(r * (1 - pi)/(1 - pi + r/d_X/r_R)^2 * var_X_EC)
  term4 <- 1 / N_RCT * (var_EC - mean(var_X_RCT_0)) 
  term5 <- 1 / N_RCT * mean(r^2/r_R/(1 - pi + r/d_X/r_D)^2 * var_X_EC)
  nu2 <- N * (term1 + term2 + term3 + term4 + term5)
 
  return(nu2) 
}



# Generate RCT data based on given sample size
f_gen_sim <- function(
    seed = 123,
    NE = 1000, # sample size of EC data
    Nt = , # sample size of treatment group in RCT data
    muX = 1, 
    sdX = 1, 
    pi = 0.5, # propensity score in RCT 
    ){
  # sed random seed for each simulation
  set.seed(seed)
  
  # generate the pre-treatment covariates
  X <- rnorm(n, mean = muX, sd = sdX)
  
  # generate A by using pi
  A <- sample(c(1,0), size = n, prob = c(pi, 1-pi))
  
}





#################################################
####### Default RCT Sample Size - Schuler #######
# Need to estimate nu for IF based simple case
{
  
  # write the optimization function using Schuler's example
  eval_f0 <- function(n,
                      alpha, 
                      beta,
                      tau,
                      nu){
    # n is the sample size for the treatment group
    # pi is the proensity score for the treatment group
    return(n)
  }
  
  # the default two sample test case
  # In the form g(n) <= 0
  eval_g0 <- function(n, 
                      alpha, 
                      beta,
                      tau,
                      nu){
    res <- 1 - beta -
      pnorm(qnorm(alpha/2) + sqrt(n) * tau / 2 / nu) -
      pnorm(qnorm(alpha/2) - sqrt(n) * tau / 2 / nu)
    return(res)
  }
  
  # Use the nonlinear optim function to find the smallest sample size
  res1 <- nloptr(x0=c(100),
                 eval_f = eval_f0,
                 lb = c(0),
                 ub = c(Inf),
                 eval_g_ineq = eval_g0,
                 opts = list("algorithm"="NLOPT_LN_COBYLA",
                             "xtol_rel"=1.0e-8),
                 alpha = 0.05,
                 beta = 0.2,
                 tau = 0.8,
                 nu = sqrt(5))
  print( res1 )
  
  
  # the IF based test case
  # In the form g(n) <= 0
  eval_g1 <- function(n, 
                      alpha, 
                      beta,
                      tau,
                      nu){
    res <- 1 - beta -
      pnorm(qnorm(alpha/2) + sqrt(n) * tau / nu) -
      pnorm(qnorm(alpha/2) - sqrt(n) * tau / nu)
    return(res)
  }
  
  # Use the nonlinear optim function to find the smallest sample size
  # need to figure out the estimated nu
  res2 <- nloptr(x0=c(100),
                 eval_f = eval_f0,
                 lb = c(0),
                 ub = c(Inf),
                 eval_g_ineq = eval_g1,
                 opts = list("algorithm"="NLOPT_LN_COBYLA",
                             "xtol_rel"=1.0e-8),
                 alpha = 0.05,
                 beta = 0.2,
                 tau = 0.8,
                 nu = )
  print( res2 )
  
  
  # this power function is matched with the following case 
  power_f <- function(n,
                      alpha,
                      beta,
                      tau,
                      nu){
    return(pnorm(qnorm(alpha/2) + sqrt(n) * tau / 2 / nu) -
             pnorm(qnorm(alpha/2) - sqrt(n) * tau / 2 / nu))
    
  }
  
  power_f(n = 246,
          alpha = 0.05,
          beta = 0.2,
          tau = 0.8,
          nu = sqrt(5))
  
  
  
}


################################################
####### Default RCT Sample Size - Guo ##########
{
  # Copy Guo's default method
  # for pi = 0.5
  # Two sample T-test
  
  beta <- c(1,1,1)
  Sigma <- matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1), nrow = 3)
  
  # smallest sample size for each group in two-sided t-test
  t_sample<- function(sigma1, # variance 1
                      sigma2, # variance 2
                      delta, # expected difference
                      alpha,
                      beta # 1-beta = 0.8
  ){
    n <- (sigma1 + sigma2)* (qnorm(alpha/2) + qnorm(beta))^2 / delta^2
    return(n)
  }
  
  # verify that paper's sample size:
  t_sample(sigma1 = 2.28^2,
           sigma2 = 2.28^2,
           delta = 0.905,
           alpha = 0.05, 
           beta = 0.2)
  
  # In our case
  t_sample(sigma1 = 5,
           sigma2 = 5,
           delta = 0.8,
           alpha = 0.05, 
           beta = 0.2)
  
}


