rm(list=ls())
library(nloptr)

# Consider a simple case 
# The treatment assignment for the RPCT data is completely at random 
# There is no unmeasured confounders



#################################################################################
####### Generate EC data and define variance function for our simple case #######
{
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
  
  # Conditional variance estimation function for EC data
  var_fun <- function(data){
    res <- 1
    return(res)
  }
  
  data_EC <- f_EC_out(N_EC = 1000)
  
  
  
  # Estimate variance of IF-based estimator using EC data
  Est_var_simple <- function(
    data_EC, # available EC data
    N_RCT, # sample size of RCT, need to be decided
    r = 1, # the ratio of variance in EC and variance in control of RCT
    d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
    pi = 0.5, # given propensity score
    gamma1 = 1 # the ratio of two kappa
    # var_fun = var_fun # given function to estimate variance
  ){
    
    # Step 0: given propensity score function - skip, let's use 1/2
    
    # Step 1: estimate V(Y|X,R=0) using data_EC
    # var_X_EC <- var_fun(data_EC)
    var_X_EC <- 1
    
    # Step 2: estimate V(Y|X,R=1,A=0)
    var_X_RCT_0 <- r * var_X_EC
    
    # Step 3: estimate kappa 0 and kappa 1
    kappa_0 <- mean(var_X_RCT_0)
    kappa_1 <- gamma1 * kappa_0
    
    # Step 4: estimate V(Y|R=1,A=1), V(Y|R=1, A=0) 
    #  Step 4 relies on the d(X) function 
    var_EC <- var(data_EC$Y_EC)
    var_RCT_0 <- var_EC
    var_RCT_1 <- var_EC
    
    # Step 5: calculate d(X) - skip, let's use the true value 
    
    # Step 6: plug in to get nu^2
    N_EC <- nrow(data_EC)
    N <- N_EC + N_RCT
    r_R <- N_RCT / N_EC
    term1 <- 1 / N_RCT / pi * mean(var_X_RCT_0)
    term2 <- 1 / N_RCT * (var_EC -  mean(var_X_RCT_0))
    term3 <- 1 / N_RCT * mean(r * (1 - pi)/(1 - pi + r / d_X / r_R)^2 * var_X_EC)
    term4 <- 1 / N_RCT * (var_EC - mean(var_X_RCT_0)) 
    term5 <- 1 / N_RCT * mean(r^2/r_R/(1 - pi + r / d_X / r_R)^2 * var_X_EC)
    nu2 <- N * (term1 + term2 + term3 + term4 + term5)
    
    return(nu2) 
  } 
  
  Est_var_simple(data_EC, N_RCT = 246)
  
}


#######################################################################
####### Use Nonlinear optimization function for our simple case #######
{
  
  # This is the main optimization formula
  eval_f0 <- function(n,
                      alpha, 
                      beta,
                      tau,
                      data_EC, # available EC data
                      r, # the ratio of variance in EC and variance in control of RCT
                      d_X, # the ratio of f(X|R=1) with f(X|R=0)
                      pi, # given propensity score
                      gamma1 # the ratio of two kappa
                      ){
    # n is the sample size for the treatment group
    # pi is the proensity score for the treatment group
    return(n)
    }
  
  eval_g1 <- function(n, 
                      alpha, 
                      beta,
                      tau,
                      data_EC, # available EC data
                      r, # the ratio of variance in EC and variance in control of RCT
                      d_X, # the ratio of f(X|R=1) with f(X|R=0)
                      pi, # given propensity score
                      gamma1 # the ratio of two kappa
                      ){
    # Step 0: given propensity score function - skip, let's use 1/2
    
    # Step 1: estimate V(Y|X,R=0) using data_EC
    # var_X_EC <- var_fun(data_EC)
    var_X_EC <- 1
    
    # Step 2: estimate V(Y|X,R=1,A=0)
    var_X_RCT_0 <- r * var_X_EC
    
    # Step 3: estimate kappa 0 and kappa 1
    kappa_0 <- mean(var_X_RCT_0)
    kappa_1 <- gamma1 * kappa_0
    
    # Step 4: estimate V(Y|R=1,A=1), V(Y|R=1, A=0) 
    #  Step 4 relies on the d(X) function 
    var_EC <- var(data_EC$Y_EC)
    var_RCT_0 <- var_EC
    var_RCT_1 <- var_EC
    
    # Step 5: calculate d(X) - skip, let's use the true value 
    
    # Step 6: plug in to get nu^2
    N_EC <- nrow(data_EC)
    N_RCT <- n - N_EC
    r_R <- N_RCT / N_EC
    term1 <- 1 / N_RCT / pi * mean(var_X_RCT_0)
    term2 <- 1 / N_RCT * (var_EC -  mean(var_X_RCT_0))
    term3 <- 1 / N_RCT * mean(r * (1 - pi)/(1 - pi + r / d_X / r_R)^2 * var_X_EC)
    term4 <- 1 / N_RCT * (var_EC - mean(var_X_RCT_0)) 
    term5 <- 1 / N_RCT * mean(r^2/r_R/(1 - pi + r / d_X / r_R)^2 * var_X_EC)
    nu2 <- n * (term1 + term2 + term3 + term4 + term5)
    sqrt_nu2 <- sqrt(nu2)
    
    res <- 1 - beta -
      pnorm(qnorm(alpha/2) + sqrt(n) * tau / sqrt_nu2) -
      pnorm(qnorm(alpha/2) - sqrt(n) * tau / sqrt_nu2)
    return(res)
  }

  
  # Use the nonlinear optim function to find the smallest sample size
  # need to figure out the estimated nu
  res2 <- nloptr(x0 = c(1200),
                 eval_f = eval_f0,
                 lb = c(1010),
                 ub = c(Inf),
                 eval_g_ineq = eval_g1,
                 opts = list("algorithm"="NLOPT_LN_COBYLA",
                             "xtol_rel"=1.0e-8),
                 data_EC = data_EC,
                 r = 1,
                 d_X = 1,
                 pi = 0.5, # pi_A
                 gamma1 = 1,
                 alpha = 0.05,
                 beta = 0.2,
                 tau = 0.8)
  print( res2 )
  
  ### Use M = 2000 EC data to calcualte the smallest sample size N_RCT 
  # Since for difference EC data, the calculated sample size is different
  # We generate 2000 EC data and calculate the average RCT sample size
  {
    # pi_A <- 0.5
    # set.seed(123)
    
    # pi_A <- 0.2
    # set.seed(234)
    
    # pi_A <- 0.4
    # set.seed(345)
    
    # pi_A <- 0.6
    # set.seed(456)
    
    pi_A <- 0.8
    set.seed(567)
    
    M <- 2000
    Index_all <- sample(1:50000000, M, replace=F)
    
    N_RCT_res <- rep(NA, M)
    for(i in 1:M){
      print(i)
      data_EC <- f_EC_out(seed_EC = Index_all[i],
                          N_EC = 1000)
      res2 <- nloptr(x0 = c(1200),
                     eval_f = eval_f0,
                     lb = c(1010),
                     ub = c(Inf),
                     eval_g_ineq = eval_g1,
                     opts = list("algorithm"="NLOPT_LN_COBYLA",
                                 "xtol_rel"=1.0e-8),
                     data_EC = data_EC,
                     r = 1,
                     d_X = 1,
                     pi = pi_A,
                     gamma1 = 1,
                     alpha = 0.05,
                     beta = 0.2,
                     tau = 0.8)
      N_RCT_res[i] <- res2$objective
    }
    
    # average N_RCT for M = 1000 EC data
    mean(N_RCT_res - 1000)
    
  }
  
  
}






############################################################
####### Generate RCT data based on given sample size #######
f_gen_sim <- function(
    seed = 123,
    NR = n_RCT, # sample size of RCT data
    Nt = , # sample size of treatment group in RCT data
    muX = 1, 
    sdX = 1, 
    pi = 0.5 # propensity score in RCT 
){
  # sed random seed for each simulation
  set.seed(seed)
  
  # generate the pre-treatment covariates
  X <- rnorm(n_RCT, mean = muX, sd = sdX)
  
  # generate A by using pi
  A <- sample(1, size = n_RCT, prob = c(pi, 1-pi))
 
   
}





###############################################
####### Nonlinear optimization function #######
####### Default RCT Sample Size - Schuler #######
# Need to estimate nu for IF based simple case
{
  ### The following optimization for the traditional RCT ###
  # write the optimization function using Schuler's example
  eval_f0 <- function(n,
                      alpha, 
                      beta,
                      tau,
                      sigma){
    # n is the sample size for the treatment group
    # pi is the proensity score for the treatment group
    return(n)
  }
  
  # the default two sample test case
  # In the form g(n) <= 0
  # Check Schuler's paper for page 2
  eval_g0 <- function(n, 
                      alpha, 
                      beta,
                      tau,
                      sigma){
    res <- 1 - beta -
      pnorm(qnorm(alpha/2) + sqrt(n) * tau / 2 / sigma) -
      pnorm(qnorm(alpha/2) - sqrt(n) * tau / 2 / sigma)
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
                 sigma = sqrt(5))
  print( res1 )
  
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



