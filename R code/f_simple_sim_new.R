rm(list=ls())
library(nloptr)
library(dplyr)
library(tidyr)
library(paletteer)
# Consider a simple case 
# The treatment assignment for the RPCT data is completely at random 
# There is no unmeasured confounders



#################################################################################
####### Generate EC data and define variance function for our simple case #######
{
  # Generate outcome Y based on linear regression in EC data
  f_gen_EC_Y <- function(x,
                         beta0,
                         beta1, 
                         sdY){
    mu <- beta0 + beta1 * x 
    return(rnorm(1, mean = mu, sd = sdY))
  }
  
  # Generate EC data
  f_gen_EC <- function(seed_EC = 123,
                       N_EC = 1000,
                       muX_EC = 1, 
                       sdX_EC = 1,
                       beta0 = 1,
                       beta1 = 2,
                       sdY = 1){
    
    set.seed(seed_EC)
    
    # generate the pre-treatment covariates
    X_EC <- rnorm(N_EC, mean = muX_EC, sd = sdX_EC)
    
    Y_EC <- mapply(f_gen_EC_Y,
                   x = X_EC, 
                   beta0, 
                   beta1,
                   sdY)
    
    return(data.frame(A = rep(0, N_EC), X = X_EC, Y = Y_EC))
    
  }

  
  data_EC <- f_gen_EC(N_EC = 1000)
  
  
  # Estimate variance of IF-based estimator using EC data
  Est_var_simple <- function(
    data_EC, # available EC data
    N_RCT, # sample size of RCT, need to be decided
    r = 1, # the ratio of variance in EC and variance in control of RCT
    d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
    pi = 0.5, # given propensity score
    gamma1 = 1, # the ratio of two kappa
    gamma = 1 # the correlation 
  ){
    
    # Step 0: given propensity score function - skip, let's use 1/2
    
    # Step 1: estimate V(Y|X,R=0) using data_EC and Chenyin's method
    mu0_model <- lm(Y ~ X, data_EC)
    var_X_EC <- mean((data_EC$Y - predict(mu0_model, data_EC))^2)
    
    # Step 2: estimate V(Y|X,R=1,A=0)
    var_X_RCT_0 <- r * var_X_EC
    
    # Step 3: estimate d(X)
    # In this case, the true value of d(X) is 1, thus we use the true value.
    
    # Step 4: estimate kappa 0 and kappa 1
    kappa_0 <- mean(var_X_RCT_0)
    kappa_1 <- gamma1 * kappa_0
    
    # Step 4: estimate V(Y|R=1,A=1), V(Y|R=1, A=0) 
    #  Step 4 relies on the d(X) function 
    var_EC <- (sd(data_EC$Y))^2
    var_RCT_0 <- var_EC
    var_RCT_1 <- var_EC
    
    # Step 5: calculate d(X) - skip, let's use the true value 
    
    # Step 6: plug in to get nu^2
    N_EC <- nrow(data_EC)
    N <- N_EC + N_RCT
    r_R <- N_RCT / N_EC
    
    term1 <- N / N_RCT / pi * kappa_1
    term2 <- N / N_RCT * ((1 - pi) * var_X_RCT_0) / ((1 - pi) + r / d_X / r_R)^2 
    
    term3 <- N / N_RCT * (var_RCT_1 - kappa_1)
    term4 <- N / N_RCT * (var_RCT_0 - kappa_0) 
    term5 <- N / N_RCT * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0) )
    
    
    term6 <- N / N_RCT * r^2  * var_X_EC / r_R / ((1 - pi) + r / d_X / r_R)^2
    
    nu2 <- (term1 + term2 + term3 + term4 - 2 * term5 + term6)
    
    return(nu2) 
  } 
  
  Est_var_simple(data_EC, N_RCT = 100)
  
  
  # this power function is matched with the following case 
  power_f <- function(n,
                      alpha,
                      beta,
                      tau,
                      nu){
    # wrong formula
    # res <- pnorm(-qnorm(1 - alpha/2) + sqrt(n) * tau / sqrt(nu)) +
    #   pnorm(qnorm(1 - alpha/2) + sqrt(n) * tau / sqrt(nu))
   
    res <- pnorm(qnorm(alpha/2) + sqrt(n) * tau / sqrt(nu)) +
      pnorm(qnorm(alpha/2) - sqrt(n) * tau / sqrt(nu))
    
    return(c(res))
    
  }
  
  power_f(n = 1100,
          alpha = 0.05,
          beta = 0.8,
          tau = 0.8,
          nu = 163.8197)
  
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
                      gamma1, # the ratio of two kappa
                      gamma # correlation
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
                      gamma1, # the ratio of two kappa
                      gamma # correlation
                      ){
    
    # Step 0: given propensity score function - skip, let's use 1/2
    
    # Step 1: estimate V(Y|X,R=0) using data_EC and Chenyin's method
    mu0_model <- lm(Y ~ X, data_EC)
    var_X_EC <- mean((data_EC$Y - predict(mu0_model, data_EC))^2)
    
    # Step 2: estimate V(Y|X,R=1,A=0)
    var_X_RCT_0 <- r * var_X_EC
    
    # Step 3: estimate d(X)
    # In this case, the true value of d(X) is 1, thus we use the true value.
    
    # Step 4: estimate kappa 0 and kappa 1
    kappa_0 <- mean(var_X_RCT_0)
    kappa_1 <- gamma1 * kappa_0
    
    # Step 4: estimate V(Y|R=1,A=1), V(Y|R=1, A=0) 
    #  Step 4 relies on the d(X) function 
    var_EC <- (sd(data_EC$Y))^2
    var_RCT_0 <- var_EC
    var_RCT_1 <- var_EC
    
    # Step 5: calculate d(X) - skip, let's use the true value 
    
    # Step 6: plug in to get nu^2
    N_EC <- nrow(data_EC)
    N_RCT <- n - N_EC
    r_R <- N_RCT / N_EC
    
    term1 <- n / N_RCT / pi * kappa_1
    term2 <- n / N_RCT * ((1 - pi) * var_X_RCT_0) / ((1 - pi) + r / d_X / r_R)^2 
    
    term3 <- n / N_RCT * (var_RCT_1 - kappa_1)
    term4 <- n / N_RCT * (var_RCT_0 - kappa_0) 
    term5 <- n / N_RCT * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0) )
    
    term6 <- n / N_RCT * r^2  * var_X_EC / r_R / ((1 - pi) + r / d_X / r_R)^2
    
    nu2 <- (term1 + term2 + term3 + term4 - 2 * term5 + term6)
    
    sqrt_nu2 <- sqrt(nu2)
    
    res <- beta - pnorm(qnorm(alpha/2) + sqrt(n) * tau / sqrt_nu2) -
      pnorm(qnorm(alpha/2) - sqrt(n) * tau / sqrt_nu2)
    
    
    return(res)
  }

  
  # Use the nonlinear optim function to find the smallest sample size
  # need to figure out the estimated nu
  res2 <- nloptr(x0 = c(1200),
                 eval_f = eval_f0,
                 lb = c(1002),
                 ub = c(Inf),
                 eval_g_ineq = eval_g1,
                 opts = list("algorithm"="NLOPT_LN_COBYLA",
                             "xtol_rel"=1.0e-8),
                 data_EC = data_EC,
                 r = 1,
                 d_X = 1,
                 pi = 0.8, # pi_A
                 gamma1 = 1,
                 gamma = 1,
                 alpha = 0.05,
                 beta = 0.8,
                 tau = 0.8)
  print( res2 )
  
  ### Use M = 2000 EC data to calcualte the smallest sample size N_RCT 
  # Since for difference EC data, the calculated sample size is different
  # We generate 2000 EC data and calculate the average RCT sample size
  {
    
    # pi_A <- 0.2
    # set.seed(234)
    
    # pi_A <- 0.4
    # set.seed(345)
    
    # pi_A <- 0.5
    # set.seed(123)
    
    # pi_A <- 0.6
    # set.seed(456)
    
    pi_A <- 0.8
    set.seed(567)
    
    M <- 2000
    Index_all <- sample(1:50000000, M, replace=F)
    
    N_RCT_res <- rep(NA, M)
    
    for(i in 1:M){
      print(i)
      data_EC <- f_gen_EC(seed_EC = Index_all[i],
                          N_EC = 1000)
      res2 <- nloptr(x0 = c(1200),
                     eval_f = eval_f0,
                     lb = c(1002),
                     ub = c(Inf),
                     eval_g_ineq = eval_g1,
                     opts = list("algorithm"="NLOPT_LN_COBYLA",
                                 "xtol_rel"=1.0e-8),
                     data_EC = data_EC,
                     r = 1,
                     d_X = 1,
                     pi = pi_A, # pi_A
                     gamma1 = 1,
                     gamma = 1,
                     alpha = 0.05,
                     beta = 0.8,
                     tau = 0.8)
      N_RCT_res[i] <- res2$objective
    }
    
    
    # average N_RCT for M = 1000 EC data
    mean(N_RCT_res - 1000)
    max(N_RCT_res - 1000)
    min(N_RCT_res - 1000)
    
  }
  
  
}



######################################################################
####### Simulation to verify the smallest sample size of N_RCT #######
{
  
  # Step 1: Generate 1000 simulated data for each case
  # Step 2: Get the point estimation and variance estimation for each simulated
  # Step 3: Do simulation 
  # Step 4: Analysis
  
  
  # Step 1: 
  {
    # Generate outcome Y based on linear regression in EC data
    f_gen_EC_Y <- function(x,
                           beta0,
                           beta1, 
                           sdY){
      mu <- beta0 + beta1 * x 
      return(rnorm(1, mean = mu, sd = sdY))
    }
    
    # Generate EC data
    f_gen_EC <- function(seed_EC = 123,
                         N_EC = 1000,
                         muX_EC = 1, 
                         sdX_EC = 1,
                         beta0 = 1,
                         beta1 = 2,
                         sdY = 1){
      
      set.seed(seed_EC)
      
      # generate the pre-treatment covariates
      X_EC <- rnorm(N_EC, mean = muX_EC, sd = sdX_EC)
      
      Y_EC <- mapply(f_gen_EC_Y,
                     x = X_EC, 
                     beta0, 
                     beta1,
                     sdY)
      
      return(data.frame(A = rep(0, N_EC), X = X_EC, Y = Y_EC))
      
    }
    
    
    # Generate outcome Y based on linear regression in RCT data
    f_gen_RCT_Y <- function(x,
                            a,
                            beta0,
                            beta1,
                            tau,
                            sdY){
      mu <- beta0 + beta1 * x + tau * a
      return(rnorm(1, mean = mu, sd = sdY))
    }
    
    # Generate RCT data
    f_gen_RCT <- function(seed = 123,
                          N_RCT, # sample size of RCT data
                          pi_A, # propensity score in RCT
                          muX_RCT = 1, 
                          sdX_RCT = 1,
                          beta0 = 1,
                          beta1 = 2,
                          sdY = 1,
                          tau
                          ){
                            
        set.seed(seed)
        
        N_crl <- ceiling(N_RCT * (1 - pi_A)) 
        N_trt <- N_RCT - N_crl 
        
        # generate the pre-treatment covariates
        X_RCT <- rnorm(N_RCT, mean = muX_RCT, sd = sdX_RCT)
        
        # generate A by using pi
        A_trt <- rep(1, N_trt)
        A_crl <- rep(0, N_crl)
        A <- c(A_trt, A_crl)
        
        Y_RCT <- mapply(f_gen_RCT_Y,
                        x = X_RCT, 
                        a = A,
                        beta0, beta1, tau, sdY)
        
        return(data.frame(A = A, X = X_RCT, Y = Y_RCT))
                            
     }
    
    # combine data_RCT and data_EC
    f_gen_full <- function(data_EC,
                           data_RCT){
      
      data_full <- rbind(data_RCT, data_EC)
      data_full$R <- c(rep(1, nrow(data_RCT)), rep(0, nrow(data_EC)) )
      return(data_full)
    }
    
    
    data_EC <- f_gen_EC(N_EC = 1000)
    data_RCT <- f_gen_RCT(N_RCT = 200,  pi_A = 0.5, tau = 0.8)
    data_full <-  f_gen_full(data_EC, data_RCT)
    
  }
  
  
  # Step 2
  {
    
    # Use the formula provided in Chenyin's paper
    f_r <- function(data, 
                    N_R,
                    N_E,
                    mu1_model,
                    mu0_model){
      data_10_r <- filter(data, R == 1 & A == 0)
      epsi_10_r <- data_10_r$Y -  predict(mu0_model, data_10_r)
      
      data_0_r <- filter(data, R == 0)
      epsi_0_r <- data_0_r$Y -  predict(mu0_model, data_0_r)
      
      res <- mean(epsi_10_r^2) / mean(epsi_0_r^2) 
      
      return(res)
      
    }
      
    # Variance estimation based on EC data
    # Estimate variance of IF-based estimator using EC data
    Est_var_simple <- function(
    data_EC, # available EC data
    N_RCT, # sample size of RCT, need to be decided
    r = 1, # the ratio of variance in EC and variance in control of RCT
    d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
    pi = 0.5, # given propensity score
    gamma1 = 1, # the ratio of two kappa
    gamma = 1 # the correlation 
  ){
    
    # Step 0: given propensity score function - skip, let's use 1/2
    
    # Step 1: estimate V(Y|X,R=0) using data_EC and Chenyin's method
    mu0_model <- lm(Y ~ X, data_EC)
    var_X_EC <- mean((data_EC$Y - predict(mu0_model, data_EC))^2)
    
    # Step 2: estimate V(Y|X,R=1,A=0)
    var_X_RCT_0 <- r * var_X_EC
    
    # Step 3: estimate d(X)
    # In this case, the true value of d(X) is 1, thus we use the true value.
    
    # Step 4: estimate kappa 0 and kappa 1
    kappa_0 <- mean(var_X_RCT_0)
    kappa_1 <- gamma1 * kappa_0
    
    # Step 4: estimate V(Y|R=1,A=1), V(Y|R=1, A=0) 
    #  Step 4 relies on the d(X) function 
    var_EC <- (sd(data_EC$Y))^2
    var_RCT_0 <- var_EC
    var_RCT_1 <- var_EC
    
    # Step 5: calculate d(X) - skip, let's use the true value 
    
    # Step 6: plug in to get nu^2
    N_EC <- nrow(data_EC)
    N <- N_EC + N_RCT
    r_R <- N_RCT / N_EC
    
    term1 <- N / N_RCT / pi * kappa_1
    term2 <- N / N_RCT * ((1 - pi) * var_X_RCT_0) / ((1 - pi) + r / d_X / r_R)^2 
    
    term3 <- N / N_RCT * (var_RCT_1 - kappa_1)
    term4 <- N / N_RCT * (var_RCT_0 - kappa_0) 
    term5 <- N / N_RCT * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0) )
    
    
    term6 <- N / N_RCT * r^2  * var_X_EC / r_R / ((1 - pi) + r / d_X / r_R)^2
    
    nu2 <- (term1 + term2 + term3 + term4 - 2 * term5 + term6)
    
    return(nu2) 
  } 
      
    # Point estimation using both RCT and EC (Chenyin)
    f_est_full <- function(data = data_full){
      R <- data$R
      A <- data$A
      X <- data$X
      Y <- data$Y
      
      data_11 <- filter(data,  R == 1 & A == 1)
      mu1_model <- lm(Y ~ X, data_11)
      
      data_10 <- filter(data,  R == 1 & A == 0)
      mu0_model <- lm(Y ~ X, data_10)

      data_1 <- filter(data, R == 1)
      pi_A_model <- glm(A ~ X, data_1, family = "binomial")
      
      pi_R_model <- glm(R ~ X, data, family = "binomial")
      
      data_0 <- filter(data, R == 0)
      
      N_R <- nrow(data_1)
      N_E <- nrow(data_0)
      
      mu1_X <- predict(mu1_model, data)
      mu0_X <- predict(mu0_model, data)
      pi_A_X <- predict(pi_A_model, data, type = "response")
      pi_R_X <- predict(pi_R_model, data, type = "response")
      epsi_1 <- data$Y - mu1_X
      epsi_0 <- data$Y - mu0_X
      q_X <- pi_R_X / (1 - pi_R_X)
      
      r_X <- f_r(data, 
                 N_R,
                 N_E,
                 mu1_model,
                 mu0_model)
      r_R <- N_R / N_E
      N <- N_R + N_E
      
      ## Estimator 1 - Point estimation - based on IF
      {
        IF_part1 <- (N_R + N_E) / N_R * (data$R * (mu1_X - mu0_X + data$A * epsi_1 / pi_A_X )) 
        IF_part2 <- (N_R + N_E) / N_R * ( 
          (data$R * (1- data$A) + (1 - data$R) * r_X) * (q_X * epsi_0) / (q_X * (1 - pi_A_X) + r_X)
        )  
        
        IF <- IF_part1 - IF_part2
        
        IF_est_tau <- mean(IF)
        
        IF_est_tau
      }
      
      ## Estimator 1 - Variance estimation - based on IF
      {
        IF_est_var <- mean(IF^2)
       
        IF_est_var
      }
      
      
      return(data.frame(IF_est_tau,  
                        IF_est_var))
    }
    
    # Point estimation using only RCT (AIPW)
    f_est_RCT <- function(data = data_RCT){
      A <- data$A
      X <- data$X
      Y <- data$Y
    
      N_R <- nrow(data)
      
      pi_A_model <- glm(A ~ X, data, family = "binomial")
      
      data_1 <- filter(data, A == 1)
      mu1_model <- lm(Y ~ X, data_1)
      
      data_0 <- filter(data, A == 0)
      mu0_model <- lm(Y ~ X, data_0)
      
      mu1_X <- predict(mu1_model, data)
      mu0_X <- predict(mu0_model, data)
      pi_A_X <- predict(pi_A_model, data, type = "response")
      epsi_1 <- data$Y - mu1_X
      epsi_0 <- data$Y - mu0_X
     
      ## Estimator 2 - Point estimation - based on IF
      {
        IF_part1 <- ((data$A * epsi_1) / pi_A_X  + mu1_X)
        IF_part2 <- (((1 - data$A) * epsi_0) / (1 - pi_A_X) + mu0_X)  
        
        IF <- IF_part1 - IF_part2
        
        IF_est_tau <- 1 / N_R * sum(IF)
        
        IF_est_tau
      }
      
      ## Estimator 2 - Variance estimation - based on IF
      {
        IF_est_var <- 1 / N_R * sum(IF^2)
        
        IF_est_var
      }
      
      return(data.frame(IF_est_tau,  
                        IF_est_var))
    }
    
    
    # Naive two-sample z-test estimator
    f_est_naive <- function(data = data_RCT,
                            delta = 0){
      A <- data$A
      X <- data$X
      Y <- data$Y
      
      N_R <- nrow(data)
      
      data_1 <- filter(data, A == 1)
      N_t <- nrow(data_1)
      data_0 <- filter(data, A == 0)
      N_c <- nrow(data_0)
      
      mu_1 <- mean(data_1$Y)
      mu_0 <- mean(data_0$Y)
      sigma_1 <- sd(data_1$Y)
      sigma_0 <- sd(data_0$Y)
      
      z_test <- ((mu_1 - mu_0) - delta) / sqrt(sigma_1^2 / N_t + sigma_0^2 / N_c)
      z_test
      
      p_value <- 1 - pnorm(abs(z_test)) + pnorm(-abs(z_test))
      
      return(data.frame(z_test, p_value))
    }
    
    
    f_est_res1 <- f_est_full(data_full)
    f_est_res1
    
    data_full_RCT <- filter(data_full, R == 1)
    f_est_res2 <-  f_est_RCT(data_full_RCT)
    f_est_res2
    
    est_var_simple <- Est_var_simple(data_EC,
                                     N_RCT = 200)
    est_var_simple

   
  
  }
  
  
  # Step 3:
  {
    
    M <- 2000
    tau <- 0.8
    
    # pi = 0.5 & N_RCT = 124
    {
      pi_A <- 0.5
      n <- 30
      
      N <- 1000 + n
      set.seed(567)
      Index_all <- sample(1:50000000, M, replace=F)
      
      Result_05_124 <- c()
      
      start <- Sys.time()
      for(i in 1:M){
        
        print(i)
        data_EC <- f_gen_EC(seed_EC = Index_all[i],
                            N_EC = 1000)
        data_RCT <- f_gen_RCT(seed = Index_all[i],
                              N_RCT = n,
                              pi_A = pi_A,
                              tau = tau)
        data_full <-  f_gen_full(data_EC, data_RCT)
        
        res1 <- f_est(data_full)
        Result_05_124 <- rbind(Result_05_124, res1)
        
        
      }
      end <- Sys.time()
      end - start
    }
   
    
    # Hypothesis H0: tau = 0.8
    {
      Z_IF_var <- (Result_05_124$IF_est_tau - 0) / sqrt(Result_05_124$IF_est_var / N)
      res_IF_var <- (abs(Z_IF_var) > qnorm(0.975)) * 1
      sum(res_IF_var) / M
      
      Z_tau2 <- (Result_05_124$IF_est_tau - 0) / sqrt(Result_05_124$est_nu2 / N)
      res_tau2 <- (abs(Z_tau2) > qnorm(0.975)) * 1
      sum(res_tau2) / M
      
      Z_prosp_tau2 <- (Result_05_124$IF_est_tau - 0) / sqrt(Result_05_124$est_prosp_nu2 / N)
      res_prosp_tau2 <- (abs(Z_prosp_tau2) > qnorm(0.975)) * 1
      sum(res_prosp_tau2) / M
    }
    
    
    # Try different n
    {
     
      M <- 2000
      tau <- 0.8
      pi_A <- 0.5
      
      f_power <- function(n,
                          pi_A = pi_A,
                          tau = tau,
                          M = M,
                          N_EC = 1000,
                          alpha = 0.05
                          ){
        N <- N_EC + n
        
        set.seed(567)
        Index_all <- sample(1:50000000, M, replace=F)
        
        Result_05_124 <- c()
        # start <- Sys.time()
        for(i in 1:M){
          print(i)
          data_EC <- f_gen_EC(seed_EC = Index_all[i],
                              N_EC = N_EC)
          data_RCT <- f_gen_RCT(seed = Index_all[i],
                                N_RCT = n,
                                pi_A = pi_A,
                                tau = tau)
          data_full <-  f_gen_full(data_EC, data_RCT)
          
          res1 <- f_est(data_full)
          Result_05_124 <- rbind(Result_05_124, res1)
        }
        # end <- Sys.time()
        # end - start
        
        Z_IF_var <- (Result_05_124$IF_est_tau - 0) / sqrt(Result_05_124$IF_est_var / N)
        res_IF_var <- (abs(Z_IF_var) > qnorm(1 - alpha/2)) * 1
        
        Z_tau2 <- (Result_05_124$IF_est_tau - 0) / sqrt(Result_05_124$est_nu2 / N)
        res_tau2 <- (abs(Z_tau2) > qnorm(1 - alpha/2)) * 1
        
        Z_prosp_tau2 <- (Result_05_124$IF_est_tau - 0) / sqrt(Result_05_124$est_prosp_nu2 / N)
        res_prosp_tau2 <- (abs(Z_prosp_tau2) > qnorm(1 - alpha/2)) * 1
        
        res <- data.frame(power_IF_var = sum(res_IF_var) / M,
                          power_tau2 = sum(res_tau2) / M,
                          power_prosp_tau2 =  sum(res_prosp_tau2) / M)
        return(res)
        
      }
      
      # f_power(n = 50,
      #         pi_A = pi_A,
      #         tau = tau,
      #         M = M,
      #         N_EC = 1000,
      #         alpha = 0.05)
      
      # n_list <- seq(4, 200, by = 2) # length = 100
      n_list <- seq(4, 200, by = 2) # length = 100
      
      Res <- c()
      for(j in 1:length(n_list)){
        n <- n_list[j]
        temp <-  f_power(n = n,
                         pi_A = pi_A,
                         tau = tau,
                         M = M,
                         N_EC = 1000,
                         alpha = 0.05)
        Res <- rbind(Res, temp)
      }
      
      
      Res$n <- n_list
      setwd("C:/Users/gaoyuji/OneDrive - Merck Sharp & Dohme LLC/Documents/GitHub/hybrid_design/R code")
      saveRDS(Res, "res_05_08_part1.RData")
      
      
    }
    

    
    
    # Wald CI  - no use
    {
      true <- 0
      
      Wald_CI <- function(alpha = 0.025,
                          point_est,
                          var_est,
                          N){
        Zstar <- qnorm(1-alpha)
        temp <- sqrt(var_est / N)
        res <- c(point_est - Zstar * temp,  point_est + Zstar * temp)
        return(res)
      }
      
      Res_CI_05_124 <- c()
      Res_sign_05_124 <- c()
      
      for(i in 1:M){
        print(i)
        
        temp <- Result_05_124[i, ]
        
        res1 <- Wald_CI(point_est = temp$IF_est_tau,
                        var_est = temp$IF_est_var,
                        N = N)
        res2 <- Wald_CI(point_est = temp$IF_est_tau,
                        var_est = temp$est_nu2,
                        N = N)
        res3 <- Wald_CI(point_est = temp$IF_est_tau,
                        var_est = temp$est_prosp_nu2,
                        N = N)
        
        res <- data.frame(IF_var_CI_l = res1[1],
                          IF_var_CI_r = res1[2],
                          nu2_CI_l = res2[1],
                          nu2_CI_r = res2[2],
                          prosp_nu2_CI_l = res3[1],
                          prosp_nu2_CI_r = res3[2]
        )
        
        sign_res1 <- res1[1] <= true & true <= res1[2]
        sign_res2 <- res2[1] <= true & true <= res2[2]
        sign_res3 <- res3[1] <= true & true <= res3[2]
        sign_res <- data.frame(IF_var_sign = sign_res1,
                               nu2_sign = sign_res2,
                               prosp_nu2_sign = sign_res3)
        
        Res_CI_05_124 <- rbind(Res_CI_05_124, res)
        Res_sign_05_124 <- rbind(Res_sign_05_124, sign_res)
        
      }
      
      1 - apply(Res_sign_05_124, 2, mean)
      
    }
    
   
    
  }
  
  
  # Step 4:
  {
    
    # Power  
    {
      rm(list=ls())
      library(tidyverse)
      library(tidyr)
      library(ggplot2)
      
      # setwd("C:/Users/gaoyuji/OneDrive - Merck Sharp & Dohme LLC/Documents/GitHub/hybrid_design/R code")
      setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
      Res <- readRDS("res_08_08.RData")
      
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 700, by = 2) 
      
      
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = "power_IF_full":"power_naive", 
          names_to = "estimator",
          values_to = "power"
        )
      
      Res_plot$estimator <- factor(Res_plot$estimator, 
                                   levels = c("power_IF_full",
                                              "power_IF_RCT",
                                              "power_naive"
                                   ))
      levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3")
      
      
      # Power
      p1 <- ggplot(data = Res_plot, aes(x = n, y = power, colour = estimator)) +
        geom_line(size = 0.8)+
        geom_hline(yintercept = 0.8,
                   linetype="dashed",
                   color = "#023743FF",
                   size = 0.8)+
        geom_vline(xintercept= 95,
                   linetype="dashed",
                   color = "#011C40FF",
                   size = 0.8)+
        geom_vline(xintercept= 115,
                   linetype="dashed",
                   color = "#035AA6FF",
                   size = 0.8)+
        geom_vline(xintercept= 133,
                   linetype="dashed",
                   color = "#05C7F2FF",
                   size = 0.8)+
        theme(legend.position = "none")+
        labs(x="n", y="Power")+
        annotate("text",x=450, y=0.6,label="beta = 0.8",size = 5)+
        labs(title= expression(paste("Power When ", tau, " = 0.8, ", pi, " = 0.8")))
      p1
      
      
      ggsave(plot = p1, 
             width = 5.5, 
             height = 5, 
             dpi = 500, 
             filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/power_08.png")
      
      
    }
    
    
    # Type I error
    {
      
      rm(list=ls())
      library(tidyverse)
      library(tidyr)
      library(ggplot2)
      
      # setwd("C:/Users/gaoyuji/OneDrive - Merck Sharp & Dohme LLC/Documents/GitHub/hybrid_design/R code")
      setwd("/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation")
      Res <- readRDS("res_08_00.RData")
      
      
      Res <- bind_rows(Res)
      Res$n <- seq(6, 700, by = 2) 
      
      
      Res_plot <-  Res %>% 
        pivot_longer(
          cols = "power_IF_full":"power_naive", 
          names_to = "estimator",
          values_to = "Type_I_Error"
        )
      
      Res_plot$estimator <- factor(Res_plot$estimator, 
                                   levels = c("power_IF_full",
                                              "power_IF_RCT",
                                              "power_naive"
                                   ))
      levels(Res_plot$estimator) <- c("Estimator 1", "Estimator 2", "Estimator 3")
      
      p2 <- ggplot(data = Res_plot, aes(x = n, y = Type_I_Error, colour = estimator)) +
      geom_line(size = 0.8)+
      geom_hline(yintercept = 0.05,
                 linetype="dashed",
                 color = "#023743FF",
                 size = 0.8)+
      labs(x="n", y="Type I Error")+
      annotate("text",x=450, y=0.1,label="alpha = 0.05",size=5)+
      labs(title= expression(paste("Type I Error When ", tau, " = 0, ", pi, " = 0.8")))+
      expand_limits(y=c(0, 0.546))+
      scale_y_continuous(breaks = seq(0, 0.546, by=0.1))
      p2
    
    
    ggsave(plot = p2, 
           width = 7, 
           height = 5, 
           dpi = 500, 
           filename = "/Users/jane.gao/Desktop/Hybrid experimental design/5. simulation/plots/TypeI_08.png")
    
    
    }
    
  }
}



####################################################################
####### Verify the prospective estimation process - variance #######
{
  # Step 1 - function
  {
    # Generate outcome Y based on linear regression in EC data
    f_gen_EC_Y <- function(x,
                           beta0,
                           beta1, 
                           sdY){
      mu <- beta0 + beta1 * x 
      return(rnorm(1, mean = mu, sd = sdY))
    }
    
    # Generate EC data
    f_gen_EC <- function(seed_EC = 123,
                         N_EC = 1000,
                         muX_EC = 1, 
                         sdX_EC = 1,
                         beta0 = 1,
                         beta1 = 2,
                         sdY = 1){
      
      set.seed(seed_EC)
      
      # generate the pre-treatment covariates
      X_EC <- rnorm(N_EC, mean = muX_EC, sd = sdX_EC)
      
      Y_EC <- mapply(f_gen_EC_Y,
                     x = X_EC, 
                     beta0, 
                     beta1,
                     sdY)
      
      return(data.frame(A = rep(0, N_EC), X = X_EC, Y = Y_EC))
      
    }
    
    
    # Generate outcome Y based on linear regression in RCT data
    f_gen_RCT_Y <- function(x,
                            a,
                            beta0,
                            beta1,
                            tau,
                            sdY){
      mu <- beta0 + beta1 * x + tau * a
      return(rnorm(1, mean = mu, sd = sdY))
    }
    
    # Generate RCT data
    f_gen_RCT <- function(seed = 123,
                          N_RCT, # sample size of RCT data
                          pi_A, # propensity score in RCT
                          muX_RCT = 1, 
                          sdX_RCT = 1,
                          beta0 = 1,
                          beta1 = 2,
                          sdY = 1,
                          tau
    ){
      
      set.seed(seed)
      
      N_crl <- ceiling(N_RCT * (1 - pi_A)) 
      N_trt <- N_RCT - N_crl 
      
      # generate the pre-treatment covariates
      X_RCT <- rnorm(N_RCT, mean = muX_RCT, sd = sdX_RCT)
      
      # generate A by using pi
      A_trt <- rep(1, N_trt)
      A_crl <- rep(0, N_crl)
      A <- c(A_trt, A_crl)
      
      Y_RCT <- mapply(f_gen_RCT_Y,
                      x = X_RCT, 
                      a = A,
                      beta0, beta1, tau, sdY)
      
      return(data.frame(A = A, X = X_RCT, Y = Y_RCT))
      
    }
    
    # combine data_RCT and data_EC
    f_gen_full <- function(data_EC,
                           data_RCT){
      
      data_full <- rbind(data_RCT, data_EC)
      data_full$R <- c(rep(1, nrow(data_RCT)), rep(0, nrow(data_EC)) )
      return(data_full)
    }
    
  }
  
  # Step 2 - function
  {
    
    # Use the formula provided in Chenyin's paper
    f_r <- function(data, 
                    N_R,
                    N_E,
                    mu1_model,
                    mu0_model){
      
      data_10_r <- filter(data, R == 1 & A == 0)
      epsi_10_r <- data_10_r$Y -  predict(mu0_model, data_10_r)
      
      data_0_r <- filter(data, R == 0)
      epsi_0_r <- data_0_r$Y -  predict(mu0_model, data_0_r)
      
      res <- mean((1 - data_10_r$A) * epsi_10_r^2) / mean(epsi_0_r^2)
      
      return(res)
      
    }
    
    # Variance estimation based on EC data
    
    # Estimate variance of IF-based estimator using EC data
    Est_var_simple <- function(
    data_EC, # available EC data
    N_RCT, # sample size of RCT, need to be decided
    r = 1, # the ratio of variance in EC and variance in control of RCT
    d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
    pi = 0.5, # given propensity score
    gamma1 = 1 # the ratio of two kappa
    ){
      
      # Step 0: given propensity score function - skip, let's use 1/2
      
      # Step 1: estimate V(Y|X,R=0) using data_EC and Chenyin's method
      mu0_model <- lm(Y ~ X, data_EC)
      var_X_EC <- mean((data_EC$Y - predict(mu0_model, data_EC))^2)
      
      # Step 2: estimate V(Y|X,R=1,A=0)
      var_X_RCT_0 <- r * var_X_EC
      
      # Step 3: estimate d(X)
      # In this case, the true value of d(X) is 1, thus we use the true value.
      
      # Step 4: estimate kappa 0 and kappa 1
      kappa_0 <- mean(var_X_RCT_0)
      kappa_1 <- gamma1 * kappa_0
      
      # Step 4: estimate V(Y|R=1,A=1), V(Y|R=1, A=0) 
      #  Step 4 relies on the d(X) function 
      var_EC <- (sd(data_EC$Y))^2
      var_RCT_0 <- var_EC
      var_RCT_1 <- var_EC
      
      # Step 5: calculate d(X) - skip, let's use the true value 
      
      # Step 6: plug in to get nu^2
      N_EC <- nrow(data_EC)
      N <- N_EC + N_RCT
      r_R <- N_RCT / N_EC
      term1 <- 1 / N_RCT / pi * kappa_1
      term2 <- 1 / N_RCT * (var_EC -  kappa_1)
      term3 <- 1 / N_RCT * (r * (1 - pi) * var_X_EC) / ((1 - pi) + r / d_X / r_R)^2 
      term4 <- 1 / N_RCT * (var_EC - kappa_0) 
      term5 <- 1 / N_RCT * r^2  * var_X_EC / r_R / ((1 - pi) + r / d_X / r_R)^2
      nu2 <- N * (term1 + term2 + term3 + term4 + term5)
      
      return(nu2) 
    }  
    
    # Point estimation using both RCT and EC (Chenyin)
    f_est_full <- function(data = data_full){
      R <- data$R
      A <- data$A
      X <- data$X
      Y <- data$Y
      
      data_11 <- filter(data,  R == 1 & A == 1)
      mu1_model <- lm(Y ~ X, data_11)
      
      data_10 <- filter(data,  R == 1 & A == 0)
      mu0_model <- lm(Y ~ X, data_10)
      
      data_1 <- filter(data, R == 1)
      pi_A_model <- glm(A ~ X, data_1, family = "binomial")
      
      pi_R_model <- glm(R ~ X, data, family = "binomial")
      
      data_0 <- filter(data, R == 0)
      
      N_R <- nrow(data_1)
      N_E <- nrow(data_0)
      
      mu1_X <- predict(mu1_model, data)
      mu0_X <- predict(mu0_model, data)
      pi_A_X <- predict(pi_A_model, data, type = "response")
      pi_R_X <- predict(pi_R_model, data, type = "response")
      epsi_1 <- data$Y - mu1_X
      epsi_0 <- data$Y - mu0_X
      q_X <- pi_R_X / (1 - pi_R_X)
      
      r_X <- f_r(data, 
                 N_R,
                 N_E,
                 mu1_model,
                 mu0_model)
      r_R <- N_R / N_E
      N <- N_R + N_E
      
      ## Estimator 1 - Point estimation - based on IF
      {
        IF_part1 <- (N_R + N_E) / N_R * (data$R * (mu1_X - mu0_X + data$A * epsi_1 / pi_A_X )) 
        IF_part2 <- (N_R + N_E) / N_R * ( 
          (data$R * (1- data$A) + (1 - data$R) * r_X) * (q_X * epsi_0) / (q_X * (1 - pi_A_X) + r_X)
        )  
        
        IF <- IF_part1 - IF_part2
        
        IF_est_tau <- mean(IF)
        
        IF_est_tau
      }
      
      ## Estimator 1 - Variance estimation - based on IF
      {
        IF_est_var <- mean(IF^2)
        
        IF_est_var
      }
      
      
      return(data.frame(IF_est_tau,  
                        IF_est_var))
    }
    
    f_est_orig <- function(data = data_full){
      R <- data$R
      A <- data$A
      X <- data$X
      Y <- data$Y
      
      data_11 <- filter(data,  R == 1 & A == 1)
      mu1_model <- lm(Y ~ X, data_11)
      
      data_10 <- filter(data,  R == 1 & A == 0)
      mu0_model <- lm(Y ~ X, data_10)
      
      data_1 <- filter(data, R == 1)
      pi_A_model <- glm(A ~ X, data_1, family = "binomial")
      
      pi_R_model <- glm(R ~ X, data, family = "binomial")
      
      data_0 <- filter(data, R == 0)
      
      N_R <- nrow(data_1)
      N_E <- nrow(data_0)
      
      mu1_X <- predict(mu1_model, data)
      mu0_X <- predict(mu0_model, data)
      pi_A_X <- predict(pi_A_model, data, type = "response")
      pi_R_X <- predict(pi_R_model, data, type = "response")
      epsi_1 <- data$Y - mu1_X # error term of Y - Y|X, R=1, A=1
      epsi_0 <- data$Y - mu0_X # error term of Y - Y|X, R=1, A=0
      q_X <- pi_R_X / (1 - pi_R_X)
      
      r_X <- f_r(data, 
                 N_R,
                 N_E,
                 mu1_model,
                 mu0_model)
      r_R <- N_R / N_E
      N <- N_R + N_E
      
      
      
      # calculate the variance based on the original formula
      {
        term1 <- N / N_R * mean( data$R * 1 / pi_A_X * (epsi_1^2) )
        term2 <- N / N_R *( (sd(data_11$Y))^2 - mean( data$R * (epsi_1^2) ) )
        term3 <- N / N_R * mean( data$R * ((1 - pi_A_X) / ((1 - pi_A_X) + r_X / q_X )^2 * (epsi_0^2) ))
        term4 <- N / N_R *( (sd(data_10$Y))^2 - mean( data$R * (epsi_0^2) ) )
        term5 <- N / N_R * mean( (1 - data$R) * ((r_X^2 / r_R) / ((1 - pi_A_X) + r_X / q_X )^2 * (epsi_0^2) ))
         
        nu <- term1 + term2 + term3 + term4 + term5 
        
      }
     
      return(nu)
      
    }
    

    
  }
  
  
  # Generate data
  data_EC <- f_gen_EC(N_EC = 1000)
  data_RCT <- f_gen_RCT(N_RCT = 1000,  pi_A = 0.5, tau = 0.8)
  data_full <-  f_gen_full(data_EC, data_RCT)
  

  # Get the IF Based variance
  f_est_res1 <- f_est_full(data_full)
  f_est_res1$IF_est_var

  # Get the var using the prospective form
  f_est_pros <- Est_var_simple(data_EC = data_EC, 
                               N_RCT = 500)
  f_est_pros

  # Get the var using the original variance form
  f_est_orig <-  f_est_orig(data = data_full)
  f_est_orig
  
  
  # Let's check the influence formula
  {
    
    psi <- function(data = data_full){
      R <- data$R
      A <- data$A
      X <- data$X
      Y <- data$Y
      
      data_11 <- filter(data,  R == 1 & A == 1)
      mu1_model <- lm(Y ~ X, data_11)
      
      data_10 <- filter(data,  R == 1 & A == 0)
      mu0_model <- lm(Y ~ X, data_10)
      
      data_1 <- filter(data, R == 1)
      pi_A_model <- glm(A ~ X, data_1, family = "binomial")
      
      pi_R_model <- glm(R ~ X, data, family = "binomial")
      
      data_0 <- filter(data, R == 0)
      
      N_R <- nrow(data_1)
      N_E <- nrow(data_0)
      
      mu1_X <- predict(mu1_model, data)
      mu0_X <- predict(mu0_model, data)
      pi_A_X <- predict(pi_A_model, data, type = "response")
      pi_R_X <- predict(pi_R_model, data, type = "response")
      epsi_1 <- Y - mu1_X # error term of Y - Y|X, R=1, A=1
      epsi_0 <- Y - mu0_X # error term of Y - Y|X, R=1, A=0
      q_X <- pi_R_X / (1 - pi_R_X)
      
      r_X <- f_r(data, 
                 N_R,
                 N_E,
                 mu1_model,
                 mu0_model)
      r_R <- N_R / N_E
      
      N <- N_R + N_E
      p_R_1 <- N_R / N
      
      
      ### Old variance formula based on 3 parts
      {
        psi_1_all <- R / p_R_1 * (A / pi_A_X * (Y - mu1_X) + (mu1_X - 0.8))
        psi_2_all <- R / p_R_1 * ((1-A) * q_X / (q_X*(1-pi_A_X) + r_X) * (Y - mu0_X) + (mu0_X - 0) )
        psi_3_all <- (1-R) / p_R_1 * r_X * q_X /  (q_X*(1-pi_A_X) + r_X) * (Y - mu0_X)
        
        # psi_1_all_new <- R / p_R_1 * (A / pi_A_X * (Y - mu1_X) + (mu1_X))
        # psi_2_all_new <- R / p_R_1 * ((1-A) * q_X / (q_X*(1-pi_A_X) + r_X) * (Y - mu0_X) + (mu0_X) )
        # IF_new <-  psi_1_all_new -  psi_2_all_new -  psi_3_all
        
        # This result is matched with the previous IF's point estimation.
        point_est_3 <- mean(psi_1_all - psi_2_all - psi_3_all) + 0.8
        
        # # Next check the mean of each term - the first two terms not 0!!!
        # mean(psi_1_all) # not 0
        # mean(psi_2_all) # not 0
        # mean(psi_3_all) # 0
        # 
        # # For term 1
        # mean(R / p_R_1 * (A / pi_A_X * (Y - mu1_X))) # part 1 in term1 is 0
        # mean(R / p_R_1 * (mu1_X - 0.8))
        # term1_1 <-  R / p_R_1 * (A / pi_A_X * (Y - mu1_X))
        # term1_2 <- R / p_R_1 * (mu1_X - 0.8)
        # 
        # # For term 2
        # mean( R / p_R_1 * ((1-A) * q_X / (q_X*(1-pi_A_X) + r_X) * (Y - mu0_X) )) # part 1 in term 2 is 0
        # mean( R / p_R_1 * (mu0_X - 0))
        # term2_1 <-  R / p_R_1 * ((1-A) * q_X / (q_X*(1-pi_A_X) + r_X) * (Y - mu0_X) )
        # term2_2 <- R / p_R_1 * (mu0_X - 0)
        # 
        # 
        # 
        # # Next check the covariance of these terms
        # mean( (psi_1_all - mean(psi_1_all)) * (psi_2_all - mean(psi_2_all)) ) # not zero
        # mean( (psi_1_all - mean(psi_1_all)) * (psi_3_all - mean(psi_3_all)) ) # approximate 0
        # mean( (psi_2_all - mean(psi_2_all)) * (psi_3_all - mean(psi_3_all)) ) # approximate 0
        # 
        # 
        # # Check the covariance between term1 and term2
        # mean(psi_1_all * psi_2_all)
        
        
        # Check the final result - this is the correct result using psi_1 + psi_2 + psi_3
        # est_3 <- (sd(psi_1_all))^2 + (sd(psi_2_all))^2 + (sd(psi_3_all))^2 - 
        #   2 * mean( (psi_1_all - mean(psi_1_all)) * (psi_2_all - mean(psi_2_all)) ) -
        #   2 * mean( (psi_1_all - mean(psi_1_all)) * (psi_3_all - mean(psi_3_all)) ) +
        #   2 * mean( (psi_2_all - mean(psi_2_all)) * (psi_3_all - mean(psi_3_all)) ) 
        
        est_3 <- mean((psi_1_all - psi_2_all - psi_3_all)^2)
        
      }
      
      
      ### Estimation process using four parts of new variance form (4)
      {
        term1 <-  R / p_R_1 * (A / pi_A_X * epsi_1)
        term2 <-  R / p_R_1 * ((1 - A) * q_X * epsi_0 / (q_X * (1 - pi_A_X) + r_X)  )
        term3 <-  R / p_R_1 * ((mu1_X - 0.8) - (mu0_X - 0))
        term4 <-  (1 - R) / p_R_1 * (r_X * q_X * epsi_0 / (q_X * (1 - pi_A_X) + r_X)  )
        
        ### check covariance - matched with theory
        # cov(term1, term2)
        # cov(term1, term4)
        # cov(term2, term4)
        # cov(term3, term4)
        # cov(term1, term3)
        # cov(term2, term3)
        
        point_est_4 <- mean(term1 - term2 + term3 - term4) + 0.8
        
        # Check the final result - this is the correct result using psi_1 + psi_2 + psi_3 + psi_4
        est_4 <- var(term1) + var(term2) + var(term3) + var(term4)
        
      }
      
      
      ### Estimation process using the new variance form
      {

        data_11$X2 <- data_11$X^2
        mu1_model_2 <- lm(Y^2 ~ X + X2, data_11)
        # Y_2_11 <- predict(mu1_model_2, data)
        
        data_10$X2 <- data_10$X^2
        mu0_model_2 <- lm(Y^2 ~ X+ X2, data_10)
        # Y_2_10 <- predict(mu0_model_2, data)
        
        data_1$X2 <- data_1$X^2
        data_0$X2 <- data_0$X^2
        
        
        # check term1 -- done! 
        {
          # term1 <- N / N_R * mean(R / pi_A_X * (Y_2_11 - (mu1_X)^2) )
          # term1 <- 1 / p_R_1 * mean(R * A / pi_A_X^2 * epsi_1^2)
          # the previous two are wrong formula
          
          term1_id <-  R / p_R_1 * (A / pi_A_X * epsi_1) 
          mean(term1_id)
          var(term1_id)
          mean(term1_id^2)
          
          # Equation 3, same as mean(term1_id^2)
          term1_eq3 <- 1 / p_R_1 * mean(data_1$A /  (predict(pi_A_model, data_1, type = "response") )^2 * 
                                          (data_1$Y - predict(mu1_model, data_1))^2 )
          all.equal(mean(term1_id^2), term1_eq3)
          
          # Equation 5, similar
          term1_eq5 <- 1 / p_R_1 * mean( 1 /  (predict(pi_A_model, data_1, type = "response") ) * 
                                           (predict(mu1_model_2, data_1) - (predict(mu1_model, data_1))^2 ) )
          
          # Equation 6 using true value of V(Y|X,R=1,A=1), similar
          term1_eq6 <- 1 / p_R_1 * mean( 1 /  (predict(pi_A_model, data_1, type = "response") ))
          
          
          # Use density to check equation
          {
            # # density for X|R=1 
            # p_X <- ecdf(data_1$X) # get the empirical cdf for X|R=1
            # 
            # data_1_new <- data_1[order(data_1$X), ]
            # sort_X_1 <- data_1_new$X
            # p_X_cdf <- p_X(sort_X_1)
            # p_X_cdf1 <- c(0,  p_X_cdf[1: (N_R - 1)])
            # p_X_df <-  p_X_cdf -  p_X_cdf1
            # 
            # f_X <- data_1_new$A /(predict(pi_A_model, data_1_new, type = "response"))^2  *  (data_1_new$Y - predict(mu1_model, data_1_new)) ^2 
            # term1 <- 1 / p_R_1 * sum(f_X * p_X_df ) 
            # # this is the right formula
            
            
            
            # # density for X|R=1,A=1
            # # p_X <- ecdf(data_11$X) # get the empirical cdf for X|R=1, A=1
            # # data_11_new <- data_11[order(data_11$X), ]
            # # sort_X_11 <- sort(data_11_new$X)
            # # p_X_cdf <- p_X(sort_X_11)
            # # p_X_cdf1 <- c(0,  p_X_cdf[1: (nrow(data_11) - 1)])
            # # p_X_df <-  p_X_cdf -  p_X_cdf1
            # p_X <- ecdf(data_1$X)
            # 
            # data_1_new <- data_1[order(data_1$X), ]
            # sort_X_1 <- data_1_new$X
            # p_X_cdf <- p_X(sort_X_1)
            # p_X_cdf1 <- c(0,  p_X_cdf[1: (N_R - 1)])
            # p_X_df <-  p_X_cdf -  p_X_cdf1
            # 
            # f1_X <- 1 / (predict(pi_A_model, data_1_new, type = "response")) * (predict(mu1_model_2, data_1_new) - (predict(mu1_model, data_1_new))^2)
            # term1 <-  1 / p_R_1 * sum(f1_X * p_X_df)   
          } 
        }
        
        
        # check term2 -- done!
        {
          term2_id <-  R / p_R_1 * ((1 - A) * q_X * epsi_0 / (q_X * (1 - pi_A_X) + r_X)  )
          mean(term2_id)
          var(term2_id)
          mean(term2_id^2)
          
          q_X_1 <- predict(pi_R_model, data_1, type = "response") / (1 - predict(pi_R_model, data_1, type = "response"))
          
          # Equation 3, same as mean(term2_id^2)
          term2_eq3 <- 1 / p_R_1 * 
            mean( (1 - data_1$A) * q_X_1^2 / (q_X_1 * (1 - predict(pi_A_model, data_1, type = "response")) + r_X)^2 * 
            (data_1$Y - predict(mu0_model, data_1))^2 
                )
          all.equal(mean(term2_id^2), term2_eq3)

          # Equation 5, similar
          term2_eq5 <- 1 / p_R_1 * 
            mean( (1 - predict(pi_A_model, data_1, type = "response")) * q_X_1^2 / (q_X_1* (1 - predict(pi_A_model, data_1, type = "response")) +  r_X)^2 * 
                    (predict(mu0_model_2, data_1) - (predict(mu0_model, data_1))^2) 
                )
          
          # Equation 6 using true value of V(Y|X,R=1,A=0), similar
          term2_eq6 <- 1 / p_R_1 * 
            mean( (1 - predict(pi_A_model, data_1, type = "response")) * q_X_1^2 / (q_X_1* (1 - predict(pi_A_model, data_1, type = "response")) +  r_X)^2  
            )
          term2_eq6
          
        }
        
        
        # check term3 -- done!
        {
          # term3 <- N / N_R * ((sd(data_11$Y))^2 - mean(R / pi_A_X * (Y - epsi_1)^2))
          # the previous one is the wrong formula
          
          term3_id <-  R / p_R_1 * ((mu1_X - 3.8) - (mu0_X - 3))
          mean(term3_id) # not equals 0
          var(term3_id)
          mean(term3_id^2)
          # variance around 0
          
          # Equation 5, similar
          mean( R / p_R_1^2 *  (predict(mu1_model, data) - 3.8)^2 )
          1 / p_R_1 * mean((predict(mu1_model, data_1) - 3.8)^2)
          term3_1_eq5 <- 1 / p_R_1 * mean(predict(mu1_model, data_1)^2 - 3.8^2)
          term3_1_eq5
          
          mean( R / p_R_1^2 *  (predict(mu0_model, data) - 3)^2 )
          1 / p_R_1 * mean((predict(mu0_model, data_1) - 3)^2)
          term3_2_eq5 <- 1 / p_R_1 * mean(predict(mu0_model, data_1)^2 - 3^2)
          term3_2_eq5
          
          mean( R / p_R_1^2 *  (predict(mu1_model, data) - 3.8) * (predict(mu0_model, data) - 3))
          1 / p_R_1 * mean((predict(mu1_model, data_1) - 3.8) * (predict(mu0_model, data_1) - 3))
          term3_3_eq5 <-  1 / p_R_1 * mean(predict(mu1_model, data_1) * predict(mu0_model, data_1) - 3.8*3)
          term3_3_eq5
          
          term3_eq5 <- term3_1_eq5 + term3_2_eq5 - 2* term3_3_eq5
          
          # Equation 6 using the true model of V(mu_a(X)|R=1)
          term3_eq6 <-  1 / p_R_1 * 4 + 1 / p_R_1 * 4 - 2 *  1 / p_R_1 * 4 
        }
        
        
        # check term4
        {
          # term4 <- N / N_R * ((sd(data_10$Y))^2 - mean(R / pi_A_X * (Y - epsi_0)^2))
          # the previous wrong formula
          
          term4_id <-  (1 - R) / p_R_1 * (r_X * q_X * epsi_0 / (q_X * (1 - pi_A_X) + r_X)  )
          mean(term4_id) # around 0
          var(term4_id)
          mean(term4_id^2)
          
          q_X_0 <- predict(pi_R_model, data_0, type = "response") / (1 - predict(pi_R_model, data_0, type = "response"))
          
          # original formula, same 
          mean((1-R) /  p_R_1^2 * r_X^2 * q_X^2 / (q_X * (1 -  pi_A_X) + r_X)^2 * (Y - predict(mu0_model, data, type = "response"))^2)
          
          # Equation 3, same 
          term4_eq3 <- (1 - p_R_1) / p_R_1^2 * mean(
            r_X^2 * q_X_0^2 / (q_X_0 * (1 - predict(pi_A_model, data_0, type = "response")) + r_X)^2 * 
              (data_0$Y - predict(mu0_model, data_0, type = "response"))^2
          )
          all.equal(term4_eq3, mean(term4_id^2))
          
          # Equation 5, similar
          term4_eq5 <- (1 - p_R_1) / p_R_1^2 * mean(
            r_X^2 * q_X_0^2 / (q_X_0 * (1 - predict(pi_A_model, data_0, type = "response")) + r_X)^2 * 
              (predict(mu0_model_2, data_0, type = "response") - (predict(mu0_model, data_0, type = "response"))^2)
          )
          
          # Equation 6 using true value of V(Y|X,R=1,A=1), similar
          term4_eq6 <- (1 - p_R_1) / p_R_1^2 * mean(
            r_X^2 * q_X_0^2 / (q_X_0 * (1 - predict(pi_A_model, data_0, type = "response")) + r_X)^2 
          )
          term4_eq6
          
        }
        
        
        est_form <- term1_eq5 + term2_eq5 + term3_eq5 + term4_eq5 
        est_form_true <- term1_eq6 + term2_eq6 + term3_eq6 + term4_eq6 
        
      }
      
      
      ### Estimation process based on IF
      {
        IF_part1 <- N / N_R * (R * (mu1_X - mu0_X + A * epsi_1 / pi_A_X )) 
        IF_part2 <- N / N_R * ( 
          (R * (1- A) + (1 - R) * r_X) * (q_X * epsi_0) / (q_X * (1 - pi_A_X) + r_X)
        )  
        
        IF <- IF_part1 - IF_part2
        IF_est <- mean(IF)
          
        IF_est_var <- mean(IF^2) 
        IF_est_var
      }
      
      
      ### verify EVVE formula
      {
        value1 <- var(data_11$Y)
        value2 <- mean( (predict(mu1_model_2, data_1) - (predict(mu1_model, data_1))^2 ) )
        value3 <- var(predict(mu1_model, data_1))
        value1 - value2
        value3
        
        value1 <- var(data_10$Y)
        value2 <- mean( (predict(mu0_model_2, data_1) - (predict(mu0_model, data_1))^2 ) )
        value3 <- var(predict(mu0_model, data_1))
        value1 - value2
        value3
        
      }
      
      
      return(list(data.frame(point_est_3, point_est_4, IF_est),
                  data.frame(est_3, est_4, IF_est_var)))
    }
    
    psi(data=data_full)
    
    
    
  }
  
}








###############################################
####### Nonlinear optimization function #######
####### Default RCT Sample Size - Schuler #######
# Need to estimate nu for IF based simple case
{
  
  # Power function
  power_f <- function(n1, # treatment
                      n2, # control
                      alpha = 0.05,
                      beta = 0.8,
                      sigma1,
                      sigma2,
                      delta){
    zalpha <- qnorm(1-alpha/2)
    res <- pnorm(-zalpha + delta / sqrt(sigma1^2/n1 + sigma2^2/n2))
    return(res)
  }
  
  # sample size calculation 
  sample_size <- function(pi_A, 
                          alpha = 0.05,
                          beta = 0.8,
                          sigma1,
                          sigma2,
                          delta){
    ratio <-  (1 - pi_A) / pi_A # n2/n1 = control : treatment
    zalpha <- qnorm(1-alpha/2)
    zbeta <- qnorm(beta)
    n1 <- (ratio * sigma1^2 + sigma2^2) / ratio * (zalpha + zbeta)^2 / delta^2
    return(data.frame(n1, n2 = ratio*n1))
  }
  
  
  ####### for pi_A = 0.5 #######
  sample_size(pi_A = 0.5, # n2/n1 = control : treatment
              alpha = 0.05,
              beta = 0.8,
              sigma1 = sqrt(5),
              sigma2 = sqrt(5),
              delta = 0.8)
  
  power_f(n1 = 123,
          n2 = 123,
          alpha = 0.05,
          beta = 0.8,
          sigma1 = sqrt(5),
          sigma2 = sqrt(5),
          delta = 0.8)
  
  
  ####### for pi_A = 0.2 #######
  sample_size(pi_A = 0.2, # n2/n1 = control : treatment
              alpha = 0.05,
              beta = 0.8,
              sigma1 = sqrt(5),
              sigma2 = sqrt(5),
              delta = 0.8)
  #       n1       n2
  # 76.64922 306.5969
  
  power_f(n1 = 77,
          n2 = 307,
          alpha = 0.05,
          beta = 0.8,
          sigma1 = sqrt(5),
          sigma2 = sqrt(5),
          delta = 0.8)
  
  
  ####### for pi_A = 0.4 #######
  sample_size(pi_A = 0.4, # n2/n1 = control : treatment
              alpha = 0.05,
              beta = 0.8,
              sigma1 = sqrt(5),
              sigma2 = sqrt(5),
              delta = 0.8)
  #       n1       n2
  # 102.199 153.2984
  
  power_f(n1 = 103, # treatment group
          n2 = 154, # control group
          alpha = 0.05,
          beta = 0.8,
          sigma1 = sqrt(5),
          sigma2 = sqrt(5),
          delta = 0.8)
  
  ####### for pi_A = 0.6 #######
  sample_size(pi_A = 0.6, # n2/n1 = control : treatment
              alpha = 0.05,
              beta = 0.8,
              sigma1 = sqrt(5),
              sigma2 = sqrt(5),
              delta = 0.8)
  #       n1       n2
  # 153.2984 102.199
  
  power_f(n1 = 154, # treatment group
          n2 = 103, # control group
          alpha = 0.05,
          beta = 0.8,
          sigma1 = sqrt(5),
          sigma2 = sqrt(5),
          delta = 0.8)
  
  
  ####### for pi_A = 0.8 #######
  sample_size(pi_A = 0.8, # n2/n1 = control : treatment
              alpha = 0.05,
              beta = 0.8,
              sigma1 = sqrt(5),
              sigma2 = sqrt(5),
              delta = 0.8)
  #       n1       n2
  # 153.2984 102.199
  
  power_f(n1 = 307, # treatment group
          n2 = 77, # control group
          alpha = 0.05,
          beta = 0.8,
          sigma1 = sqrt(5),
          sigma2 = sqrt(5),
          delta = 0.8)
  
}


## previous thing
{
  
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
  
}





