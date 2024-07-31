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
    gamma1 = 1 # the ratio of two kappa
    # var_fun = var_fun # given function to estimate variance
  ){
    
    # Step 0: given propensity score function - skip, let's use 1/2
    
    # Step 1: estimate V(Y|X,R=0) using data_EC
    # var_X_EC <- var_fun(data_EC)
    var_X_EC <- sd(data_EC$X)^2
    
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
    var_X_EC <- sd(data_EC$X)^2
    
    # Step 2: estimate V(Y|X,R=1,A=0)
    var_X_RCT_0 <- r * var_X_EC
    
    # Step 3: estimate kappa 0 and kappa 1
    kappa_0 <- mean(var_X_RCT_0)
    kappa_1 <- gamma1 * kappa_0
    
    # Step 4: estimate V(Y|R=1,A=1), V(Y|R=1, A=0) 
    #  Step 4 relies on the d(X) function 
    var_EC <- var(data_EC$Y)
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
    data_RCT <- f_gen_RCT(N_RCT = 500,  pi_A = 0.5, tau = 0.8)
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
      
      res <- 1 / N_R * sum(epsi_10_r^2) / 
             (1 / N_E * sum(epsi_0_r^2) )
      
      return(res)
      
    }
      
    # Variance estimation based on EC data
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
      var_X_EC <- sd(data_EC$X)^2
      
      # Step 2: estimate V(Y|X,R=1,A=0)
      var_X_RCT_0 <- r * var_X_EC
      
      # Step 3: estimate kappa 0 and kappa 1
      kappa_0 <- mean(var_X_RCT_0)
      kappa_1 <- gamma1 * kappa_0
      
      # Step 4: estimate V(Y|R=1,A=1), V(Y|R=1, A=0) 
      #  Step 4 relies on the d(X) function 
      var_EC <- var(data_EC$Y)
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
    
    # est_var_simple <- Est_var_simple(data_EC, 
    #                                  N_RCT = 200)
    # est_var_simple
    # 
   
  
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





