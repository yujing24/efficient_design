rm(list=ls())
library(dplyr)
library(tidyr)
library(paletteer)

# Set tau = 0.4 for power and tau = 0 for type I error
pi_A_list <- c(0.5, 0.6, 0.7, 0.8, 0.9)


###############################################
####### Functions of generating EC data #######
{
  
  ### main paper simulation and sensitivity analysis 1 ###
  f_gen_EC_Y <- function(x1,
                         x2,
                         beta0_EC,
                         beta1_EC,
                         beta2_EC, 
                         sdY_EC){
    mu <- beta0_EC + beta1_EC * x1 + beta2_EC*x2  
    return(rnorm(1, mean = mu, sd = sdY_EC))
  }
  
  f_gen_EC_suff <- function(seed_EC = 123,
                            N_EC = 1000,
                            muX1_EC = 1, 
                            sdX1_EC = 1,
                            pX2_EC = 0.5,
                            beta0_EC = 1,
                            beta1_EC = 0.5,
                            beta2_EC = -1,
                            sdY_EC = 1){
    
    set.seed(seed_EC)
    
    # generate the pre-treatment covariates
    X1_EC <- rnorm(N_EC, mean = muX1_EC, sd = sdX1_EC)
    X2_EC <- rbinom(N_EC, size = 1, prob = pX2_EC)
    
    Y_EC <- mapply(f_gen_EC_Y,
                   x1 = X1_EC,
                   x2 = X2_EC, 
                   beta0_EC,
                   beta1_EC,
                   beta2_EC, 
                   sdY_EC)
    
    return(data.frame(A = rep(0, N_EC), X1 = X1_EC, X2 = X2_EC, Y = Y_EC))
    
  }
  
  f_gen_EC_insuff <- function(seed_EC = 123,
                              N_EC = 60,
                              muX1_EC = 1.2, 
                              sdX1_EC = sqrt(1.5),
                              pX2_EC = 0.7,
                              beta0_EC = 1,
                              beta1_EC = 0.5,
                              beta2_EC = -1,
                              sdY_EC = 1){
    
    set.seed(seed_EC)
    
    # generate the pre-treatment covariates
    X1_EC <- rnorm(N_EC, mean = muX1_EC, sd = sdX1_EC)
    X2_EC <- rbinom(N_EC, size = 1, prob = pX2_EC)
    
    Y_EC <- mapply(f_gen_EC_Y,
                   x1 = X1_EC,
                   x2 = X2_EC, 
                   beta0_EC,
                   beta1_EC,
                   beta2_EC, 
                   sdY_EC)
    
    return(data.frame(A = rep(0, N_EC), X1 = X1_EC, X2 = X2_EC, Y = Y_EC))
    
  }
  
  
  ### sensitivity analysis 2 ###
  f_gen_EC_Y_sensi2 <- function(x1,
                         x2,
                         beta0_EC,
                         beta1_EC,
                         beta2_EC, 
                         epi_EC,
                         sdY_EC){
    mu <- beta0_EC + beta1_EC * x1 + beta2_EC*x2  
    epi <- rnorm(1, mean = 0, sd = sdY_EC)
    
    return(mu +  epi_EC * x1^2 * epi)
  }
  
  f_gen_EC_sensi2 <- function(seed_EC = 123,
                             N_EC = 1000,
                             muX1_EC = 1, 
                             sdX1_EC = 1,
                             pX2_EC = 0.5,
                             beta0_EC = 1,
                             beta1_EC = 0.5,
                             beta2_EC = -1,
                             epi_EC = 0.4,
                             sdY_EC = 1){
    
    set.seed(seed_EC)
    
    # generate the pre-treatment covariates
    X1_EC <- rnorm(N_EC, mean = muX1_EC, sd = sdX1_EC)
    X2_EC <- rbinom(N_EC, size = 1, prob = pX2_EC)
    
    Y_EC <- mapply(f_gen_EC_Y_sensi2,
                   x1 = X1_EC,
                   x2 = X2_EC, 
                   beta0_EC,
                   beta1_EC,
                   beta2_EC,
                   epi_EC,
                   sdY_EC)
    
    return(data.frame(A = rep(0, N_EC), X1 = X1_EC, X2 = X2_EC, Y = Y_EC))
    
  }
  
}

####################################################
####### 1. Sample Size (Z-test) - true value #######
{
  
  ### Functions ###
  {
    # Power function for verifying
    power_f <- function(n1, # treatment
                        n2, # control
                        alpha = 0.05,
                        beta = 0.2,
                        sigma1,
                        sigma2,
                        delta){
      zalpha <- qnorm(1-alpha/2)
      res1 <- pnorm(-zalpha + delta / sqrt(sigma1^2/n1 + sigma2^2/n2))
      res2 <- 1 - pnorm(zalpha - delta / sqrt(sigma1^2/n1 + sigma2^2/n2)) +
        pnorm(- zalpha - delta / sqrt(sigma1^2/n1 + sigma2^2/n2))
      return(c(res1, res2))
    }
    
    # Sample size calculation 
    size_diff_in_means <- function(pi_A, 
                                   alpha = 0.05,
                                   beta = 0.2,
                                   sigma1,
                                   sigma2,
                                   delta){
      # 1: treat
      # 2: control
      # n2/n1 = control : treatment
      res <- c()
      for(i in 1:length(pi_A)){
        ratio <- pi_A[i]/(1 - pi_A[i]) # n2/n1 = control : treatment
        zalpha <- qnorm(1-alpha/2)
        zbeta <- qnorm(1 - beta)
        n1 <- (sigma1^2 + ratio * sigma2^2) * (zalpha + zbeta)^2 / delta^2
        res <- rbind(res, c(pi_A[i], ceiling(n1), n2 = ceiling(n1 / ratio), ceiling(n1) + ceiling(n1 / ratio)))
      }
      res <- as.data.frame(res)
      colnames(res) <- c("pi_A", "N_t", "N_c", "size")
      return(res)
    }
    
  }
  
  ## Simulation section in the main paper
  std_size <- size_diff_in_means(pi_A = pi_A_list, 
                                 alpha = 0.05,
                                 beta = 0.2,
                                 sigma1 = sqrt(1.3),
                                 sigma2 = sqrt(1.3), 
                                 delta = 0.4)
  std_size
  
  # # verify for pi_A = 0.5
  # power_f(n1 = 128,
  #         n2 = 128,
  #         alpha = 0.05,
  #         beta = 0.2,
  #         sigma1 = sqrt(1.3),
  #         sigma2 = sqrt(1.3),
  #         delta = 0.4)
 
  ## Sensitivity analysis 1
  std_size_sensi1 <- size_diff_in_means(pi_A = pi_A_list,
                                         alpha = 0.05,
                                         beta = 0.2,
                                         sigma1 = sqrt(1.6),
                                         sigma2 = sqrt(1.3), 
                                         delta = 0.4)
  std_size_sensi1
  
  ## Sensitivity analysis 2
  std_size_sensi2 <- size_diff_in_means(pi_A = pi_A_list,
                                         alpha = 0.05,
                                         beta = 0.2,
                                         sigma1 = sqrt(1.524),
                                         sigma2 = sqrt(1.524),
                                         delta = 0.4)
  std_size_sensi2
  
}

#######################################################################
####### 2. RCT Sample Size based on AIPW Estimator - true value #######
{
  
  ### Functions ###
  {
    ## Functions for simulation section in the main paper and sensitivity analysis 1
    eval_aipw <- function(n, 
                          alpha, 
                          beta,
                          pi_A, # given propensity score
                          tau,
                          var_X_EC = 1, 
                          var_EC = 1.5,
                          r0_M = 1.3/1.5,
                          r1_M = 1.3/1.5,
                          r = 0.8, # the ratio of variance in EC and variance in control of RCT
                          gamma1 = 1,
                          gamma = 1 
    ){
      
      # Step 1: Estimate V(Y|X, R=0, A=0) and V(Y|X, R=0, A=0)
      
      # Step 2: estimate V(Y|R=1, A=1), V(Y|R=1, A=0) 
      var_RCT_0 <- var_EC * r0_M
      var_RCT_1 <- var_EC * r1_M
      
      # Step 3: estimate V(Y|X, R=1, A=0) using r(X)
      var_X_RCT_0 <- r * var_X_EC  
      
      # Step 4: estimate kappa 0 and kappa 1
      kappa_0 <- mean(var_X_RCT_0)
      kappa_1 <- gamma1 * kappa_0
      
      term1 <-  (var_RCT_1 + ((1 - pi_A) / pi_A) * kappa_1)
      term2 <-  (var_RCT_0 + (pi_A / (1 - pi_A)) * kappa_0)  
      term3 <-  -2 * gamma * sqrt((var_RCT_0 - kappa_0) * (var_RCT_1 - kappa_1))
      
      nu2 <- (term1 + term2 + term3)
      
      sqrt_nu2 <- sqrt(nu2)
      
      res <- 1 - beta - pnorm(qnorm(alpha/2) + sqrt(n) * tau / sqrt_nu2) -
        pnorm(qnorm(alpha/2) - sqrt(n) * tau / sqrt_nu2)
      
      return(res)
      
    }
    
    # Grid search function for sample size of the RCT-only AIPW  
    size_aipw <- function(start = 1,
                          end = 500,
                          alpha = 0.05,
                          beta = 0.2,
                          pi_A = 0.5,
                          tau = 0.4,
                          var_X_EC = 1, 
                          var_EC = 1.5,
                          r0_M = 1.3/1.5,
                          r1_M = 1.3/1.5,
                          r = 0.8, 
                          gamma1 = 1,
                          gamma = 1){
      
      size_length <- end - start + 1
      AIPW_size <- matrix(NA, nrow = size_length, ncol = 2)
      
      for(i in start:end){
        power_diff <- eval_aipw(n = i,
                                alpha = alpha,
                                beta = beta,
                                pi_A = pi_A,
                                tau = tau,
                                var_X_EC = var_X_EC, 
                                var_EC = var_EC,
                                r0_M = r0_M,
                                r1_M = r1_M,
                                r = r, 
                                gamma1 = gamma1,
                                gamma = gamma)
        AIPW_size[i, ] <- c(i, power_diff)
        if(power_diff <= 0){break}
      }
      return(list(i,
                  AIPW_size))
      
    }
    
    
    ## Functions for sensitivity analysis 2
    eval_aipw_sensi2 <- function(n, 
                                 alpha, 
                                 beta,
                                 pi_A, # given propensity score
                                 tau,
                                 var_X_EC = function(x){0.16*x^4}, 
                                 var_EC = 2.1,
                                 r0_M = 1.524/2.1,
                                 r1_M = 1.524/2.1,
                                 r = function(x){3.2/x^2}, 
                                 dist_para = list(mu_R = 1,
                                                  sigma_R = 1,
                                                  p_R = 0.5,
                                                  mu_E = 1,
                                                  sigma_E = 1,
                                                  p_E = 0.5),
                                 sample_size = 5*10^5,
                                 gamma1 = 1,
                                 gamma = 1 
    ){
      
      # Step 1: Estimate V(Y|X, R=0, A=0) and V(Y|X, R=0, A=0)
      
      # Step 2: estimate V(Y|R=1, A=1), V(Y|R=1, A=0) 
      var_RCT_0 <- var_EC * r0_M
      var_RCT_1 <- var_EC * r1_M
      
      # Step 3: estimate V(Y|X, R=1, A=0) using r(X)
      var_X_RCT_0 <- function(x){
        return(r(x) * var_X_EC(x))
      }
      
      # Step 4: estimate kappa 0 and kappa 1
      X1_RCT_sample <- rnorm(sample_size, mean = dist_para$mu_R, sd = dist_para$sigma_R)
      kappa_0 <- mean(var_X_RCT_0(X1_RCT_sample))
      kappa_1 <- gamma1 * kappa_0
      
      term1 <-  (var_RCT_1 + ((1 - pi_A) / pi_A) * kappa_1)
      term2 <-  (var_RCT_0 + (pi_A / (1 - pi_A)) * kappa_0)  
      term3 <-  -2 * gamma * sqrt((var_RCT_0 - kappa_0) * (var_RCT_1 - kappa_1))
      
      nu2 <- (term1 + term2 + term3)
      
      sqrt_nu2 <- sqrt(nu2)
      
      res <- 1 - beta - pnorm(qnorm(alpha/2) + sqrt(n) * tau / sqrt_nu2) -
        pnorm(qnorm(alpha/2) - sqrt(n) * tau / sqrt_nu2)
      
      return(res)
      
    }
    
    # Grid search function for sample size of the RCT-only AIPW 
    size_aipw_sensi2 <- function(start = 1,
                                 end = 600,
                                 alpha = 0.05,
                                 beta = 0.2,
                                 pi_A = 0.5,
                                 tau = 0.4,
                                 var_X_EC = function(x){0.16*x^4}, 
                                 var_EC = 2.1,
                                 r0_M = 1.524/2.1,
                                 r1_M = 1.524/2.1,
                                 r = function(x){3.2/x^2}, 
                                 dist_para = list(mu_R = 1,
                                                  sigma_R = 1,
                                                  p_R = 0.5,
                                                  mu_E = 1,
                                                  sigma_E = 1,
                                                  p_E = 0.5),
                                 sample_size = 5*10^5,
                                 gamma1 = 1,
                                 gamma = 1){
      
      
      size_length <- end - start + 1
      AIPW_size <- matrix(NA, nrow = size_length, ncol = 2)
      
      for(i in start:end){
        power_diff <- eval_aipw_sensi2(n = i,
                                       alpha,
                                       beta,
                                       pi_A,
                                       tau,
                                       var_X_EC = var_X_EC, 
                                       var_EC = var_EC,
                                       r0_M = r0_M,
                                       r1_M = r1_M,
                                       r = r, 
                                       dist_para = dist_para,
                                       sample_size = sample_size,
                                       gamma1 = gamma1,
                                       gamma = gamma)
        AIPW_size[i, ] <- c(i, power_diff)
        if(power_diff <= 0){break}
      }
      return(list(i,
                  AIPW_size))
      
    }
    
    
  }
  
  ## Simulation section in the main paper
  aipw_size <- data.frame(pi = pi_A_list, 
                          size = rep(NA, 5),
                          ratio = rep(NA, 5))
  for(i in 1:5){
    aipw_size[i, 2] <- size_aipw(start = 1,
                                 end = 500,
                                 alpha = 0.05,
                                 beta = 0.2,
                                 pi_A = pi_A_list[i], 
                                 tau = 0.4)[[1]]
    aipw_size[i, 3] <- round((std_size$size[i] - aipw_size[i, 2]) / std_size$size[i], 4) 
  }
  aipw_size
  

  ## Sensitivity analysis 1
  aipw_size_sensi1 <- data.frame(pi = pi_A_list, 
                                size = rep(NA, 5),
                                ratio = rep(NA, 5))
  for(i in 1:5){
    aipw_size_sensi1[i, 2] <- size_aipw(start = 1,
                                        end = 500,
                                        alpha = 0.05,
                                        beta = 0.2,
                                        pi_A = pi_A_list[i],
                                        tau = 0.4,
                                        var_X_EC = 1, 
                                        var_EC = 1.5,
                                        r0_M = 1.3/1.5,
                                        r1_M = 1.6/1.5,
                                        r = 0.8,
                                        gamma1 = 1,
                                        gamma = 0.6/sqrt(0.4))[[1]]
    
    aipw_size_sensi1[i, 3] <- round((std_size_sensi1$size[i] - aipw_size_sensi1[i, 2]) / std_size_sensi1$size[i], 4) 
  }
  aipw_size_sensi1
  
  ## Sensitivity analysis 2
  aipw_size_sensi2 <- data.frame(pi = pi_A_list, 
                                 size = rep(NA, 5),
                                 ratio = rep(NA, 5))
  for(i in 1:5){
    aipw_size_sensi2[i, 2] <- size_aipw_sensi2(start = 1,
                                               end = 600,
                                               alpha = 0.05,
                                               beta = 0.2,
                                               pi_A = pi_A_list[i], 
                                               tau = 0.4)[[1]]
    
    aipw_size_sensi2[i, 3] <- round((std_size_sensi2$size[i] - aipw_size_sensi2[i, 2]) / std_size_sensi2$size[i], 4) 
  }
  aipw_size_sensi2
  
}

#########################################################################
####### 3.1 hybrid design - true tuning parameter and true models #######
{
 
  ### Functions ###
  {
    ## Functions for simulation section in the main paper and sensitivity analysis 1
    eval_hybrid_true <- function(N_RCT, 
                                 alpha, 
                                 beta,
                                 tau,
                                 pi_A, 
                                 N_EC = 1000, # 1000 for sufficient case and 60 for insufficient case
                                 var_X_EC = 1, 
                                 var_EC = 1.5, # sufficient  = 1.5, insufficient = 1.585
                                 r0_M = 1.3/1.5, # V(Y|R=1, A=0) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                 r1_M = 1.3/1.5, # V(Y|R=1, A=1) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                 r = 0.8, # the ratio of variance in EC and variance in control of RCT
                                 d_X = 1, # the ratio of f(X|R=1) with f(X|R=0), 1 for sufficient, NULL for insufficient
                                 dist_para = list(mu_R = 1,
                                                  sigma_R = 1,
                                                  p_R = 0.5,
                                                  mu_E = 1.2,
                                                  sigma_E = sqrt(1.5),
                                                  p_E = 0.7), # if var_X_EC or r is function, then we must have dist_para
                                 sample_size = 5*10^5,
                                 gamma1 = 1, # the ratio of two kappa
                                 gamma = 1 # correlation
    ){
      
      n <- N_RCT + N_EC
      r_R <- N_RCT / N_EC
      
      # Step 1: Estimate V(Y|X, R=0, A=0) and V(Y|X, R=0, A=0)
      
      # Step 2: estimate V(Y|R=1, A=1), V(Y|R=1, A=0) 
      var_RCT_0 <- var_EC * r0_M
      var_RCT_1 <- var_EC * r1_M
      
      # Step 3: estimate V(Y|X, R=1, A=0) using r(X)
      var_X_RCT_0 <- r * var_X_EC  
      
      # Step 4: estimate kappa 0 and kappa 1
      kappa_0 <- mean(var_X_RCT_0)
      kappa_1 <- gamma1 * kappa_0
      
      # Step 5 - 7: d(X), gamma, variance
      if(length(d_X) < 1){
        
        pmf_bern <- function(x2, p){
          return(x2 * p + (1 - x2) * (1 - p))
        }
        
        d_X_true <- function(x1, x2, 
                             mu_R, sigma_R,
                             p_R,
                             mu_E, sigma_E,
                             p_E
        ){
          
          res <- dnorm(x1, mean = mu_R, sd = sigma_R) * pmf_bern(x2, p_R) / dnorm(x1, mean = mu_E, sd = sigma_E) / pmf_bern(x2, p_E) 
          return(res)
          
        }
        
        dX_x1_R <- rnorm(sample_size, mean = dist_para$mu_R, sd = dist_para$sigma_R)
        dX_x2_R <- rbinom(sample_size, size = 1, prob = dist_para$p_R)
        d_X_R <-  d_X_true(x1 = dX_x1_R,
                           x2 = dX_x2_R, 
                           dist_para$mu_R, dist_para$sigma_R, dist_para$p_R,
                           dist_para$mu_E, dist_para$sigma_E, dist_para$p_E)
        
        dX_x1_E <- rnorm(sample_size, mean = dist_para$mu_E, sd = dist_para$sigma_E)
        dX_x2_E <- rbinom(sample_size, size = 1, prob = dist_para$p_E)
        d_X_E <-  d_X_true(x1 = dX_x1_E,
                           x2 = dX_x2_E, 
                           dist_para$mu_R, dist_para$sigma_R, dist_para$p_R,
                           dist_para$mu_E, dist_para$sigma_E, dist_para$p_E)
        
        # step 
        term1 <- kappa_1 / pi_A
        term2 <- mean(((1 - pi_A) * var_X_RCT_0) / ((1 - pi_A) + r / d_X_R / r_R)^2)
        term3 <- (var_RCT_1 - kappa_1) + (var_RCT_0 - kappa_0) -
          2 * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0))
        term4 <- mean(r^2  * var_X_EC / r_R / ((1 - pi_A) + r / d_X_E / r_R)^2)
        
      }else{
        
        term1 <- kappa_1 / pi_A
        term2 <- ((1 - pi_A) * var_X_RCT_0) / ((1 - pi_A) + r / d_X / r_R)^2
        term3 <- (var_RCT_1 - kappa_1) + (var_RCT_0 - kappa_0) -
          2 * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0))
        term4 <- r^2  * var_X_EC / r_R / ((1 - pi_A) + r / d_X / r_R)^2
        
      }
      
      nu2 <- (term1 + term2 + term3 + term4)
      
      sqrt_nu2 <- sqrt(nu2)
      
      res <- 1 - beta - pnorm(qnorm(alpha/2) + sqrt(N_RCT) * tau / sqrt_nu2) -
        pnorm(qnorm(alpha/2) - sqrt(N_RCT) * tau / sqrt_nu2)
      
      return(res)
    } 
    
    size_hybrid_true <- function(start = 1,
                                 end = 500,
                                 alpha = 0.05,
                                 beta = 0.2,
                                 pi_A = 0.5,
                                 tau = 0.4,
                                 N_EC = 1000, # 1000 for sufficient case and 60 for insufficient case
                                 var_X_EC = 1, 
                                 var_EC = 1.5, # sufficient  = 1.5, insufficient = 1.585
                                 r0_M = 1.3/1.5, # V(Y|R=1, A=0) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                 r1_M = 1.3/1.5, # V(Y|R=1, A=1) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                 r = 0.8, # the ratio of variance in EC and variance in control of RCT
                                 d_X = 1, # the ratio of f(X|R=1) with f(X|R=0), 1 for sufficient, NULL for insufficient
                                 dist_para = list(mu_R = 1,
                                                  sigma_R = 1,
                                                  p_R = 0.5,
                                                  mu_E = 1.2,
                                                  sigma_E = sqrt(1.5),
                                                  p_E = 0.7),
                                 sample_size = 5*10^5,
                                 gamma1 = 1, # the ratio of two kappa
                                 gamma = 1 # correlation
    ){
      
      size_length <- end - start + 1
      Hybrid_size <- c()
      
      for(i in start:end){
        
        power_diff <- eval_hybrid_true(N_RCT = i,
                                       alpha = alpha,
                                       beta = beta,
                                       tau = tau,
                                       pi_A = pi_A,
                                       N_EC = N_EC, 
                                       var_X_EC = var_X_EC, 
                                       var_EC = var_EC,
                                       r0_M = r0_M, 
                                       r1_M = r1_M,
                                       r = r, 
                                       d_X = d_X,
                                       dist_para = dist_para,
                                       sample_size = sample_size,
                                       gamma1 = gamma1,
                                       gamma = gamma)
        Hybrid_size <- rbind(Hybrid_size, c(i, power_diff))
        if(power_diff <= 0){break}
        
      }
      
      return(list(i,
                  Hybrid_size))
      
    }
    
    ## Functions for sensitivity analysis 2
    eval_hybrid_sensi2 <- function(N_RCT, 
                                   alpha, 
                                   beta,
                                   pi_A,
                                   tau, 
                                   N_EC = 1000, # 1000 for sufficient case and 60 for insufficient case
                                   var_X_EC = function(x){0.16*x^4}, 
                                   var_EC = 2.1,
                                   r0_M = 1.524/2.1, # V(Y|R=1, A=0) / V(Y|R=0)
                                   r1_M = 1.524/2.1, # V(Y|R=1, A=1) / V(Y|R=0)
                                   r = function(x){3.2/x^2}, # the ratio of variance in EC and variance in control of RCT
                                   d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
                                   dist_para = list(mu_R = 1,
                                                    sigma_R = 1,
                                                    p_R = 0.5,
                                                    mu_E = 1,
                                                    sigma_E = 1,
                                                    p_E = 0.5),
                                   sample_size = 5*10^5,
                                   gamma1 = 1, # the ratio of two kappa
                                   gamma = 1 # correlation
    ){
      
      n <- N_RCT + N_EC
      r_R <- N_RCT / N_EC
      
      # Step 1: Estimate V(Y|X, R=0, A=0) and V(Y|X, R=0, A=0)
      
      # Step 2: estimate V(Y|R=1, A=1), V(Y|R=1, A=0) 
      var_RCT_0 <- var_EC * r0_M
      var_RCT_1 <- var_EC * r1_M
      
      # Step 3: estimate V(Y|X, R=1, A=0) using r(X)
      var_X_RCT_0 <- function(x){
        return(r(x) * var_X_EC(x))
      }
      
      # Step 4: estimate kappa 0 and kappa 1
      X1_RCT_sample <- rnorm(sample_size, mean = dist_para$mu_R, sd = dist_para$sigma_R)
      kappa_0 <- mean(var_X_RCT_0(X1_RCT_sample))
      kappa_1 <- gamma1 * kappa_0
      
      # Step 5: d(X)
      
      # Step 6: gamma
      
      # Step 7: 
      X1_EC_sample <- rnorm(sample_size, mean = dist_para$mu_E, sd = dist_para$sigma_E)
      
      term1 <- kappa_1 / pi_A
      term2 <- mean(((1 - pi_A) * var_X_RCT_0(X1_RCT_sample)) / ((1 - pi_A) + r(X1_RCT_sample) / d_X / r_R)^2)
      
      term3 <- (var_RCT_1 - kappa_1) + (var_RCT_0 - kappa_0) -
        2 * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0) )
      
      term4 <- mean(r(X1_EC_sample)^2  * var_X_EC(X1_EC_sample) / r_R / ((1 - pi_A) + r(X1_EC_sample) / d_X / r_R)^2)
      
      nu2 <- (term1 + term2 + term3 + term4)
      
      sqrt_nu2 <- sqrt(nu2)
      
      res <- 1 - beta - pnorm(qnorm(alpha/2) + sqrt(N_RCT) * tau / sqrt_nu2) -
        pnorm(qnorm(alpha/2) - sqrt(N_RCT) * tau / sqrt_nu2)
      
      return(res)
    }
    
    size_hybrid_sensi2 <- function(start = 1,
                                   end = 500,
                                   alpha = 0.05,
                                   beta = 0.2,
                                   pi_A = 0.5,
                                   tau = 0.4,
                                   N_EC = 1000, # 1000 for sufficient case and 60 for insufficient case
                                   var_X_EC = function(x){0.16*x^4}, 
                                   var_EC = 2.1,
                                   r0_M = 1.524/2.1, # V(Y|R=1, A=0) / V(Y|R=0)
                                   r1_M = 1.524/2.1, # V(Y|R=1, A=1) / V(Y|R=0)
                                   r = function(x){3.2/x^2}, # the ratio of variance in EC and variance in control of RCT
                                   d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
                                   dist_para = list(mu_R = 1,
                                                    sigma_R = 1,
                                                    p_R = 0.5,
                                                    mu_E = 1,
                                                    sigma_E = 1,
                                                    p_E = 0.5),
                                   sample_size = 5*10^5,
                                   gamma1 = 1, # the ratio of two kappa
                                   gamma = 1 # correlation
    ){
      
      
      size_length <- end - start + 1
      Hybrid_size <- c()
      
      for(i in start:end){
        
        power_diff <- eval_hybrid_sensi2(N_RCT = i,
                                         alpha = alpha,
                                         beta = beta,
                                         tau = tau,
                                         pi_A = pi_A,
                                         N_EC = N_EC, # 1000 for sufficient case and 60 for insufficient case
                                         var_X_EC = var_X_EC, 
                                         var_EC = var_EC,
                                         r0_M = r0_M, # V(Y|R=1, A=0) / V(Y|R=0)
                                         r1_M = r1_M, # V(Y|R=1, A=1) / V(Y|R=0)
                                         r = r, # the ratio of variance in EC and variance in control of RCT
                                         d_X = d_X, # the ratio of f(X|R=1) with f(X|R=0)
                                         dist_para = dist_para,
                                         sample_size = sample_size,
                                         gamma1 = gamma1, # the ratio of two kappa
                                         gamma = gamma)
        Hybrid_size <- rbind(Hybrid_size, c(i, power_diff))
        if(power_diff <= 0){break}
      }
      
      return(list(i,
                  Hybrid_size))
      
    }
  }
  
  ## Simulation section - sufficient EC (sample size 1000)
  hd_size_suff <- data.frame(pi = pi_A_list,
                             size = rep(NA, 5),
                             ratio = rep(NA, 5))
  for(i in 1:5){
    hd_size_suff[i, 2] <- size_hybrid_true(start = 10,
                                             end = 500,
                                             alpha = 0.05,
                                             beta = 0.2,
                                             pi_A = pi_A_list[i],
                                             tau = 0.4)[[1]]
    hd_size_suff[i, 3] <- round((std_size$size[i] - hd_size_suff[i, 2]) / std_size$size[i], 4) 
  }
  hd_size_suff
  
  ## Simulation section - insufficient EC (sample size 60)
  hd_size_insuff_60 <- data.frame(pi = pi_A_list, 
                                  size = rep(NA, 5),
                                  ratio = rep(NA, 5))
  for(i in 1:5){
    hd_size_insuff_60[i, 2] <-  size_hybrid_true(start = 10,
                                                 end = 500,
                                                 alpha = 0.05,
                                                 beta = 0.2,
                                                 pi_A = pi_A_list[i], # pi = 0.5
                                                 tau = 0.4,
                                                 N_EC = 60, # 1000 for sufficient case and 60 for insufficient case
                                                 var_X_EC = 1, 
                                                 var_EC = 1.585, # sufficient  = 1.5, insufficient = 1.585
                                                 r0_M = 1.3/1.585, # V(Y|R=1, A=0) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                                 r1_M = 1.3/1.585, # V(Y|R=1, A=1) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                                 r = 0.8, # the ratio of variance in EC and variance in control of RCT
                                                 d_X = NULL, # the ratio of f(X|R=1) with f(X|R=0), 1 for sufficient, NULL for insufficient
                                                 dist_para = list(mu_R = 1,
                                                                  sigma_R = 1,
                                                                  p_R = 0.5,
                                                                  mu_E = 1.2,
                                                                  sigma_E = sqrt(1.5),
                                                                  p_E = 0.7),
                                                 sample_size = 5*10^5, 
                                                 gamma1 = 1, # the ratio of two kappa
                                                 gamma = 1 # correlation
    )[[1]]
    hd_size_insuff_60[i, 3] <- round((std_size$size[i] - hd_size_insuff_60[i, 2]) / std_size$size[i], 4) 
  }
  hd_size_insuff_60
  
 
  ## Simulation section - insufficient EC (sample size 30)
  hd_size_insuff_30 <- data.frame(pi = pi_A_list, 
                                  size = rep(NA, 5),
                                  ratio = rep(NA, 5))
  for(i in 1:5){
    hd_size_insuff_30[i, 2] <- size_hybrid_true(start = 10,
                                                end = 500,
                                                alpha = 0.05,
                                                beta = 0.2,
                                                pi_A = pi_A_list[i], # pi = 0.5
                                                tau = 0.4,
                                                N_EC = 30, # 1000 for sufficient case and 60 for insufficient case
                                                var_X_EC = 1, 
                                                var_EC = 1.585, # sufficient  = 1.5, insufficient = 1.585
                                                r0_M = 1.3/1.585, # V(Y|R=1, A=0) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                                r1_M = 1.3/1.585, # V(Y|R=1, A=1) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                                r = 0.8, # the ratio of variance in EC and variance in control of RCT
                                                d_X = NULL, # the ratio of f(X|R=1) with f(X|R=0), 1 for sufficient, NULL for insufficient
                                                dist_para = list(mu_R = 1,
                                                                 sigma_R = 1,
                                                                 p_R = 0.5,
                                                                 mu_E = 1.2,
                                                                 sigma_E = sqrt(1.5),
                                                                 p_E = 0.7),
                                                sample_size = 5*10^5,
                                                gamma1 = 1, # the ratio of two kappa
                                                gamma = 1 # correlation
    )[[1]]
    hd_size_insuff_30[i, 3] <- round((std_size$size[i] - hd_size_insuff_30[i, 2]) / std_size$size[i], 4) 
  }
  hd_size_insuff_30

  
  ## Sensitivity section 1 - sufficient EC (sample size 1000)
  hd_size_sensi1 <- data.frame(pi = pi_A_list,
                             size = rep(NA, 5),
                             ratio = rep(NA, 5))
  for(i in 1:5){
    hd_size_sensi1[i, 2] <- size_hybrid_true(start = 10,
                                             end = 500,
                                             alpha = 0.05,
                                             beta = 0.2,
                                             pi_A = pi_A_list[i],
                                             tau = 0.4,
                                             N_EC = 1000, 
                                             var_X_EC = 1, 
                                             var_EC = 1.5, 
                                             r0_M = 1.3/1.5, 
                                             r1_M = 1.6/1.5, 
                                             r = 0.8, 
                                             d_X = 1, 
                                             dist_para = list(mu_R = 1,
                                                              sigma_R = 1,
                                                              p_R = 0.5,
                                                              mu_E = 1.2,
                                                              sigma_E = sqrt(1.5),
                                                              p_E = 0.7),
                                             sample_size = 5*10^5,
                                             gamma1 = 1,
                                             gamma = 0.6/sqrt(0.4)
    )[[1]]
    hd_size_sensi1[i, 3] <- round((std_size_sensi1$size[i] - hd_size_sensi1[i, 2]) / std_size_sensi1$size[i], 4) 
  }
  hd_size_sensi1
  

  ## Sensitivity section 2 - sufficient EC (sample size 1000)
  hd_size_sensi2 <- data.frame(pi = pi_A_list,
                               size = rep(NA, 5),
                               ratio = rep(NA, 5))
  for(i in 1:5){
    hd_size_sensi2[i, 2] <-  size_hybrid_sensi2(start = 10,
                                                end = 500,
                                                alpha = 0.05,
                                                beta = 0.2,
                                                pi_A = pi_A_list[i],
                                                tau = 0.4)[[1]]
    hd_size_sensi2[i, 3] <- round((std_size_sensi2$size[i] - hd_size_sensi2[i, 2]) / std_size_sensi2$size[i], 4) 
  }
  hd_size_sensi2
  
  
}

#########################################################################################
####### 3.2 hybrid design - non-informative tuning parameter and estimated models #######
{

  ### Functions ###
  {
    eval_hybrid_noninfo <- function(N_RCT,
                                    alpha, 
                                    beta,
                                    tau,
                                    pi_A, 
                                    data_EC, # available EC data
                                    r0_M = 1, 
                                    r1_M = 1,
                                    r = 1,
                                    d_X = 1,
                                    dist_para = list(mu_R = 1,
                                                     sigma_R = 1,
                                                     p_R = 0.5,
                                                     mu_E = 1.2,
                                                     sigma_E = sqrt(1.5),
                                                     p_E = 0.7),
                                    sample_size = 5*10^5,
                                    gamma1 = 1, # the ratio of two kappa
                                    gamma = 1 # correlation
    ){
      
      N_EC <- nrow(data_EC)
      n <- N_RCT + N_EC
      r_R <- N_RCT / N_EC
      
      # Step 1: Estimate V(Y|X, R=0, A=0) and V(Y|X, R=0, A=0)
      mu0_model <- lm(Y ~ X1 + X2, data_EC)
      var_X_EC <- mean((data_EC$Y - predict(mu0_model, data_EC))^2)
      var_EC <- (sd(data_EC$Y))^2
      
      # Step 2: estimate V(Y|R=1, A=1), V(Y|R=1, A=0) 
      var_RCT_0 <- var_EC * r0_M
      var_RCT_1 <- var_EC * r1_M
      
      # Step 3: estimate V(Y|X, R=1, A=0) using r(X)
      var_X_RCT_0 <- r * var_X_EC
      
      # Step 4: estimate kappa 0 and kappa 1
      kappa_0 <- mean(var_X_RCT_0)
      kappa_1 <- gamma1 * kappa_0
      
      # Step 5 - 7: d(X), gamma, variance
      if(length(d_X) < 1){
        
        pmf_bern <- function(x2, p){
          return(x2 * p + (1 - x2) * (1 - p))
        }
        
        d_X_true <- function(x1, x2, 
                             mu_R, sigma_R,
                             p_R,
                             mu_E, sigma_E,
                             p_E
        ){
          
          res <- dnorm(x1, mean = mu_R, sd = sigma_R) * pmf_bern(x2, p_R) / dnorm(x1, mean = mu_E, sd = sigma_E) / pmf_bern(x2, p_E) 
          return(res)
          
        }
        
        dX_x1_R <- rnorm(sample_size, mean = dist_para$mu_R, sd = dist_para$sigma_R)
        dX_x2_R <- rbinom(sample_size, size = 1, prob = dist_para$p_R)
        d_X_R <-  d_X_true(x1 = dX_x1_R,
                           x2 = dX_x2_R, 
                           dist_para$mu_R, dist_para$sigma_R, dist_para$p_R,
                           dist_para$mu_E, dist_para$sigma_E, dist_para$p_E)
        
        dX_x1_E <- rnorm(sample_size, mean = dist_para$mu_E, sd = dist_para$sigma_E)
        dX_x2_E <- rbinom(sample_size, size = 1, prob = dist_para$p_E)
        d_X_E <-  d_X_true(x1 = dX_x1_E,
                           x2 = dX_x2_E, 
                           dist_para$mu_R, dist_para$sigma_R, dist_para$p_R,
                           dist_para$mu_E, dist_para$sigma_E, dist_para$p_E)
        
        term1 <- kappa_1 / pi_A
        term2 <- mean(((1 - pi_A) * var_X_RCT_0) / ((1 - pi_A) + r / d_X_R / r_R)^2)
        term3 <- (var_RCT_1 - kappa_1) + (var_RCT_0 - kappa_0) -
          2 * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0) )
        term4 <- mean(r^2  * var_X_EC / r_R / ((1 - pi_A) + r / d_X_E / r_R)^2)
        
      }else{
        
        term1 <- kappa_1 / pi_A
        term2 <- ((1 - pi_A) * var_X_RCT_0) / ((1 - pi_A) + r / d_X / r_R)^2
        term3 <- (var_RCT_1 - kappa_1) + (var_RCT_0 - kappa_0) -
          2 * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0))
        term4 <- r^2  * var_X_EC / r_R / ((1 - pi_A) + r / d_X / r_R)^2
        
      }
      
      nu2 <- (term1 + term2 + term3 + term4)
      
      sqrt_nu2 <- sqrt(nu2)
      
      res <- 1 - beta - pnorm(qnorm(alpha/2) + sqrt(N_RCT) * tau / sqrt_nu2) -
        pnorm(qnorm(alpha/2) - sqrt(N_RCT) * tau / sqrt_nu2)
      
      return(res)
    }
    
    
    size_hybrid_noninfo <- function(start = 1,
                                    end = 500,
                                    alpha = 0.05,
                                    beta = 0.2,
                                    tau = 0.4,
                                    pi_A = 0.5,
                                    data_EC, # available EC data
                                    r0_M = 1, 
                                    r1_M = 1,
                                    r = 1,
                                    d_X = 1,
                                    dist_para = list(mu_R = 1,
                                                     sigma_R = 1,
                                                     p_R = 0.5,
                                                     mu_E = 1.2,
                                                     sigma_E = sqrt(1.5),
                                                     p_E = 0.7),
                                    sample_size = 5*10^5,
                                    gamma1 = 1, # the ratio of two kappa
                                    gamma = 1 # correlation
    ){
      
      size_length <- end - start + 1
      Hybrid_size <- c()
      
      for(i in start:end){
        
        power_diff <- eval_hybrid_noninfo(N_RCT = i,
                                          alpha = alpha,
                                          beta = beta,
                                          tau = tau,
                                          pi_A = pi_A,
                                          data_EC = data_EC,
                                          r0_M = r0_M,
                                          r1_M = r1_M, 
                                          r = r,
                                          d_X = d_X, 
                                          dist_para = dist_para,
                                          sample_size = sample_size,
                                          gamma1 = gamma1, 
                                          gamma = gamma)
        Hybrid_size <- rbind(Hybrid_size, c(i, power_diff))
        if(power_diff <= 0){break}
        
      }
      
      return(list(i,
                  Hybrid_size))
      
    }
    
  }

  ### Use M = 2000 EC data to calculate the smallest sample size N_RCT 
  # Since for difference EC data, the calculated sample size is different
  # We generate 2000 EC data and calculate the average RCT sample size
  M <- 2000 # replicates times
  
  ###### Simulation section - sufficient EC (sample size 1000) ######
  {
    
    set.seed(123)
    Index_all <- sample(1:50000000, M, replace=F)
    N_RCT_res_non_suff <- matrix(NA, nrow = M, ncol = 5)
    for(i in 1:M){
      print(i)
      data_EC <- f_gen_EC_suff(seed_EC = Index_all[i],
                               N_EC = 1000)
      for(j in 1:5){
        N_RCT_res_non_suff[i, j] <- size_hybrid_noninfo(start = 1,
                                                end = 500,
                                                alpha = 0.05,
                                                beta = 0.2,
                                                tau = 0.4,
                                                pi_A = pi_A_list[j],
                                                data_EC = data_EC)[[1]]
      }
      
    }
    
    # average N_RCT for M = 2000 EC data
    hd_size_non_suff <- data.frame(pi = pi_A_list,
                               size = ceiling(apply(N_RCT_res_non_suff, 2, mean)))
    hd_size_non_suff
    # hd_size_suff$ratio <- (std_size$size - hd_size_suff$size) / std_size$size
    # hd_size_suff 
  }
  
  ###### Simulation section - insufficient EC (sample size 60) ######
  {
    
    set.seed(234)
    Index_all <- sample(1:50000000, M, replace=F)
    N_RCT_res_non_60 <- matrix(NA, nrow = M, ncol = 5)
    for(i in 1:M){
      print(i)
      data_EC <- f_gen_EC_insuff(seed_EC = Index_all[i],
                                 N_EC = 60)
      for(j in 1:5){
        N_RCT_res_non_60[i, j] <- size_hybrid_noninfo(start = 1,
                                                end = 500,
                                                alpha = 0.05,
                                                beta = 0.2,
                                                tau = 0.4,
                                                pi_A = pi_A_list[j],
                                                data_EC = data_EC)[[1]]
      }
      
    }
    
    # average N_RCT for M = 2000 EC data
    hd_size_non_insuff_60 <- data.frame(pi = pi_A_list,
                                        size = ceiling(apply(N_RCT_res_non_60, 2, mean)))
    hd_size_non_insuff_60
    
    # hd_size_non_insuff_60$ratio <- (std_size$size - hd_size_non_insuff_60$size) / std_size$size
    
  }
 
  ###### Simulation section - insufficient EC (sample size 30) ######
  {

    set.seed(345)
    Index_all <- sample(1:50000000, M, replace=F)
    N_RCT_res_non_30 <- matrix(NA, nrow = M, ncol = 5) 
    for(i in 1:M){
      print(i)
      data_EC <- f_gen_EC_insuff(seed_EC = Index_all[i],
                                 N_EC = 30)
      for(j in 1:5){
        N_RCT_res_non_30[i, j] <- size_hybrid_noninfo(start = 1,
                                                  end = 500,
                                                  alpha = 0.05,
                                                  beta = 0.2,
                                                  tau = 0.4,
                                                  pi_A = pi_A_list[j],
                                                  data_EC = data_EC)[[1]]
      }
    }
    
    # average N_RCT for M = 2000 EC data
    hd_size_non_insuff_30 <- data.frame(pi = pi_A_list,
                                        size = ceiling(apply(N_RCT_res_non_30, 2, mean)))
    hd_size_non_insuff_30

    # hd_size_non_insuff_30$ratio <- (std_size$size - hd_size_non_insuff_30$size) / std_size$size
    
  }
  
  ###### Sensitivity section 1 - sufficient EC (sample size 1000) ######
  # same as Simulation section - sufficient EC (sample size 1000)
  hd_size_non_suff
  
  ###### Sensitivity section 2 - sufficient EC (sample size 1000)  ######
  {
    
    set.seed(456)
    Index_all <- sample(1:50000000, M, replace=F)
    N_RCT_res_non_sensi2 <- matrix(NA, nrow = M, ncol = 5) 
    for(i in 1:M){
      print(i)
      data_EC <- f_gen_EC_sensi2(seed_EC = Index_all[i],
                                 N_EC = 1000)
      for(j in 1:5){
        N_RCT_res_non_sensi2[i,j] <- size_hybrid_noninfo(start = 1,
                                                      end = 500,
                                                      alpha = 0.05,
                                                      beta = 0.2,
                                                      tau = 0.4,
                                                      pi_A = pi_A_list[j],
                                                      data_EC = data_EC)[[1]]
      }
      
    }
    
    # average N_RCT for M = 2000 EC data
    hd_size_non_sensi2 <- data.frame(pi = pi_A_list,
                                        size = ceiling(apply(N_RCT_res_non_sensi2, 2, mean)))
    hd_size_non_sensi2
    
    # hd_size_non_sensi2$ratio <- (std_size$size - hd_size_non_sensi2$size) / std_size$size
    
  }
  
}

###########################################################################
####### 4.1 Single arm case - true tuning parameter and true models #######
####### only for sufficient EC case #######
{
 
  ### Functions ###
  {
    ## Functions for simulation section in the main paper and sensitivity analysis 1
    eval_sa_true <- function(N_RCT, 
                             alpha, 
                             beta,
                             tau,
                             N_EC = 1000, # 1000 for sufficient case and 60 for insufficient case
                             var_X_EC = 1, 
                             var_EC = 1.5, # sufficient  = 1.5, insufficient = 1.585
                             r0_M = 1.3/1.5, # V(Y|R=1, A=0) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                             r1_M = 1.3/1.5, # V(Y|R=1, A=1) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                             r = 0.8, # the ratio of variance in EC and variance in control of RCT
                             d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
                             dist_para = list(mu_R = 1,
                                              sigma_R = 1,
                                              p_R = 0.5,
                                              mu_E = 1.2,
                                              sigma_E = sqrt(1.5),
                                              p_E = 0.7),
                             sample_size =  5*10^5,
                             gamma1 = 1, # the ratio of two kappa
                             gamma = 1 # correlation
    ){
      
      n <- N_RCT + N_EC
      r_R <- N_RCT / N_EC
      
      # Step 1: Estimate V(Y|X, R=0, A=0) and V(Y|X, R=0, A=0)
      
      # Step 2: estimate V(Y|R=1, A=1), V(Y|R=1, A=0) 
      var_RCT_0 <- var_EC * r0_M
      var_RCT_1 <- var_EC * r1_M
      
      # Step 3: estimate V(Y|X, R=1, A=0) using r(X)
      var_X_RCT_0 <- r * var_X_EC  
      
      # Step 4: estimate kappa 0 and kappa 1
      kappa_0 <- mean(var_X_RCT_0)
      kappa_1 <- gamma1 * kappa_0
      
      # Step 5 - 7: d(X), gamma, variance
      if(length(d_X) < 1){
        
        pmf_bern <- function(x2, p){
          return(x2 * p + (1 - x2) * (1 - p))
        }
        
        d_X_true <- function(x1, x2, 
                             mu_R, sigma_R,
                             p_R,
                             mu_E, sigma_E,
                             p_E
        ){
          
          res <- dnorm(x1, mean = mu_R, sd = sigma_R) * pmf_bern(x2, p_R) / dnorm(x1, mean = mu_E, sd = sigma_E) / pmf_bern(x2, p_E) 
          return(res)
          
        }
        
        dX_x1_E <- rnorm(sample_size, mean = dist_para$mu_E, sd = dist_para$sigma_E)
        dX_x2_E <- rbinom(sample_size, size = 1, prob = dist_para$p_E)
        d_X_E <-  d_X_true(x1 = dX_x1_E,
                           x2 = dX_x2_E, 
                           dist_para$mu_R, dist_para$sigma_R, dist_para$p_R,
                           dist_para$mu_E, dist_para$sigma_E, dist_para$p_E)
        
        term1 <- kappa_1
        term3 <- (var_RCT_1 - kappa_1) + (var_RCT_0 - kappa_0) - 
          2 * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0))
        term4 <- mean((d_X_E * r_R)^2 / r_R  * var_X_EC)
        
      }else{
        term1 <- kappa_1
        term3 <- (var_RCT_1 - kappa_1) + (var_RCT_0 - kappa_0) - 
          2 * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0) )
        term4 <- (d_X * r_R)^2 / r_R  * var_X_EC 
        
      }
      
      nu2 <- (term1 + term3 + term4)
      
      sqrt_nu2 <- sqrt(nu2)
      
      res <- 1 - beta - pnorm(qnorm(alpha/2) + sqrt(N_RCT) * tau / sqrt_nu2) -
        pnorm(qnorm(alpha/2) - sqrt(N_RCT) * tau / sqrt_nu2)
      
      return(res)
    }
    
    size_sa_true <- function(start = 1,
                             end = 500,
                             alpha = 0.05,
                             beta = 0.2,
                             tau = 0.4,
                             N_EC = 1000, # 1000 for sufficient case and 60 for insufficient case
                             var_X_EC = 1, 
                             var_EC = 1.5, # sufficient  = 1.5, insufficient = 1.585
                             r0_M = 1.3/1.5, # V(Y|R=1, A=0) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                             r1_M = 1.3/1.5, # V(Y|R=1, A=1) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                             r = 0.8, # the ratio of variance in EC and variance in control of RCT
                             d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
                             dist_para = list(mu_R = 1,
                                              sigma_R = 1,
                                              p_R = 0.5,
                                              mu_E = 1.2,
                                              sigma_E = sqrt(1.5),
                                              p_E = 0.7),
                             sample_size =  5*10^5,
                             gamma1 = 1, # the ratio of two kappa
                             gamma = 1 # correlation
    ){
      
      size_length <- end - start + 1
      SA_size <- c()
      
      for(i in start:end){
        
        power_diff <- eval_sa_true(N_RCT = i,
                                   alpha = alpha,
                                   beta = beta,
                                   tau = tau,
                                   N_EC = N_EC,
                                   var_X_EC = var_X_EC, 
                                   var_EC = var_EC, 
                                   r0_M = r0_M, 
                                   r1_M = r1_M, 
                                   r = r,
                                   d_X = d_X,
                                   dist_para = dist_para,
                                   sample_size = sample_size,
                                   gamma1 = gamma1, 
                                   gamma = gamma)
        
        SA_size <- rbind(SA_size, c(i, power_diff))
        if(power_diff <= 0){break}
      }
      
      return(list(i,
                  SA_size))
      
    }
    
    # Function for verifying the sufficient feasibility condition (under true setting)
    eval_sa_verify_true <- function(alpha = 0.05, 
                                    beta = 0.2,
                                    tau = 0.4,
                                    N_EC = 60, # 1000 for sufficient case and 60 for insufficient case
                                    var_X_EC = 1, 
                                    d_X = NULL, # the ratio of f(X|R=1) with f(X|R=0)
                                    dist_para = list(mu_R = 1,
                                                     sigma_R = 1,
                                                     p_R = 0.5,
                                                     mu_E = 1.2,
                                                     sigma_E = sqrt(1.5),
                                                     p_E = 0.7),
                                    sample_size =  5*10^5,
                                    gamma1 = 1, # the ratio of two kappa
                                    gamma = 1 # correlation
    ){
      
      # Step 6: estimate d(X)
      if(length(d_X) < 1){
        
        pmf_bern <- function(x2, p){
          return(x2 * p + (1 - x2) * (1 - p))
        }
        
        d_X_true <- function(x1, x2, 
                             mu_R, sigma_R,
                             p_R,
                             mu_E, sigma_E,
                             p_E
        ){
          
          res <- dnorm(x1, mean = mu_R, sd = sigma_R) * pmf_bern(x2, p_R) / dnorm(x1, mean = mu_E, sd = sigma_E) / pmf_bern(x2, p_E) 
          return(res)
          
        }
        
        dX_x1_E <- rnorm(sample_size, mean = dist_para$mu_E, sd = dist_para$sigma_E)
        dX_x2_E <- rbinom(sample_size, size = 1, prob = dist_para$p_E)
        d_X_E <-  d_X_true(x1 = dX_x1_E,
                           x2 = dX_x2_E, 
                           dist_para$mu_R, dist_para$sigma_R, dist_para$p_R,
                           dist_para$mu_E, dist_para$sigma_E, dist_para$p_E)
        
        condition <- mean((d_X_E)^2 * var_X_EC) * (qnorm(1 - beta) - qnorm(alpha/2))^2 / tau^2
        
      }else{
        
        condition <- (d_X)^2 * var_X_EC * (qnorm(1 - beta) - qnorm(alpha/2))^2 / tau^2
        
      }
      
      if(N_EC > condition){
        return(list(TRUE, condition))
      }else{
        return(list(FALSE, condition))
      }
      
    }
   
    ## Functions for sensitivity analysis 2
    eval_sa_sensi2_true <- function(N_RCT, 
                                    alpha, 
                                    beta,
                                    tau,
                                    N_EC = 1000, # 1000 for sufficient case and 60 for insufficient case
                                    var_X_EC = function(x){0.16*x^4}, 
                                    var_EC = 2.1,
                                    r0_M = 1.524/2.1, # V(Y|R=1, A=0) / V(Y|R=0)
                                    r1_M = 1.524/2.1, # V(Y|R=1, A=1) / V(Y|R=0)
                                    r = function(x){3.2/x^2}, # the ratio of variance in EC and variance in control of RCT
                                    d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
                                    dist_para = list(mu_R = 1,
                                                     sigma_R = 1,
                                                     p_R = 0.5,
                                                     mu_E = 1,
                                                     sigma_E = 1,
                                                     p_E = 0.5),
                                    sample_size = 5*10^5,
                                    gamma1 = 1, # the ratio of two kappa
                                    gamma = 1 # correlation
    ){
      
      n <- N_RCT + N_EC
      r_R <- N_RCT / N_EC
      
      # Step 1: Estimate V(Y|X, R=0, A=0) and V(Y|X, R=0, A=0)
      
      # Step 2: estimate V(Y|R=1, A=1), V(Y|R=1, A=0) 
      var_RCT_0 <- var_EC * r0_M
      var_RCT_1 <- var_EC * r1_M
      
      # Step 3: estimate V(Y|X, R=1, A=0) using r(X)
      var_X_RCT_0 <- function(x){
        return(r(x) * var_X_EC(x))
      }
      
      # Step 4: estimate kappa 0 and kappa 1
      X1_RCT_sample <- rnorm(sample_size, mean = dist_para$mu_R, sd = dist_para$sigma_R)
      kappa_0 <- mean(var_X_RCT_0(X1_RCT_sample))
      kappa_1 <- gamma1 * kappa_0
      
      # Step 5: d(X)
      
      # Step 6: gamma
      
      # Step 7: 
      X1_EC_sample <- rnorm(sample_size, mean = dist_para$mu_E, sd = dist_para$sigma_E)
      
      term1 <- kappa_1
      term3 <- (var_RCT_1 - kappa_1) + (var_RCT_0 - kappa_0) -  
        2 * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0) )
      term4 <- mean((d_X * r_R)^2 / r_R  * var_X_EC(X1_EC_sample))
      
      nu2 <- (term1 + term3 + term4)
      
      sqrt_nu2 <- sqrt(nu2)
      
      res <- 1 - beta - pnorm(qnorm(alpha/2) + sqrt(N_RCT) * tau / sqrt_nu2) -
        pnorm(qnorm(alpha/2) - sqrt(N_RCT) * tau / sqrt_nu2)
      
      return(res)
    }
    
    size_sa_sensi2_true <- function(start = 1,
                                    end = 500,
                                    alpha = 0.05,
                                    beta = 0.2,
                                    tau = 0.4,
                                    N_EC = 1000, # 1000 for sufficient case and 60 for insufficient case
                                    var_X_EC = function(x){0.16*x^4}, 
                                    var_EC = 2.1,
                                    r0_M = 1.524/2.1, # V(Y|R=1, A=0) / V(Y|R=0)
                                    r1_M = 1.524/2.1, # V(Y|R=1, A=1) / V(Y|R=0)
                                    r = function(x){3.2/x^2}, # the ratio of variance in EC and variance in control of RCT
                                    d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
                                    dist_para = list(mu_R = 1,
                                                     sigma_R = 1,
                                                     p_R = 0.5,
                                                     mu_E = 1,
                                                     sigma_E = 1,
                                                     p_E = 0.5),
                                    sample_size = 5*10^5,
                                    gamma1 = 1, # the ratio of two kappa
                                    gamma = 1 # correlation
    ){
      
      size_length <- end - start + 1
      SA_size <- c()
      
      for(i in start:end){
        
        power_diff <- eval_sa_sensi2_true(N_RCT = i,
                                          alpha = alpha,
                                          beta = beta,
                                          tau = tau,
                                          N_EC = N_EC,
                                          var_X_EC = var_X_EC, 
                                          var_EC = var_EC, 
                                          r0_M = r0_M, 
                                          r1_M = r1_M, 
                                          r = r,
                                          d_X = d_X,
                                          dist_para = dist_para,
                                          sample_size = sample_size,
                                          gamma1 = gamma1, 
                                          gamma = gamma)
        
        SA_size <- rbind(SA_size, c(i, power_diff))
        if(power_diff <= 0){break}
      }
      
      return(list(i,
                  SA_size))
      
    }
    
  }
  
  ###### Simulation section - sufficient EC (sample size 1000) ######
  sa_size_suff <- size_sa_true(start = 1,
                   end = 500,
                   alpha = 0.05,
                   beta = 0.2,
                   tau = 0.4)[[1]]
  sa_ratio_suff <- round((std_size$size - sa_size_suff) / std_size$size, 4)
  sa_size_suff
  sa_ratio_suff
  
  ###### Simulation section - insufficient EC (sample size 60) ######
  eval_sa_verify_true(N_EC = 60)
  
  ###### Simulation section - insufficient EC (sample size 30) ######
  eval_sa_verify_true(N_EC = 30)
   
  ###### Sensitivity section 1 - sufficient EC (sample size 1000) ######
  sa_size_sensi1 <- size_sa_true(start = 1,
                                 end = 500,
                                 alpha = 0.05,
                                 beta = 0.2,
                                 tau = 0.4,
                                 N_EC = 1000,
                                 var_X_EC = 1, 
                                 var_EC = 1.5, 
                                 r0_M = 1.3/1.5, # V(Y|R=1, A=0) / V(Y|R=0), 
                                 r1_M = 1.6/1.5, # V(Y|R=1, A=1) / V(Y|R=0), 
                                 r = 0.8, # the ratio of variance in EC and variance in control of RCT
                                 d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
                                 dist_para = list(mu_R = 1,
                                                  sigma_R = 1,
                                                  p_R = 0.5,
                                                  mu_E = 1.2,
                                                  sigma_E = sqrt(1.5),
                                                  p_E = 0.7),
                                 sample_size =  5*10^5,
                                 gamma1 = 1, # the ratio of two kappa
                                 gamma = 0.6/sqrt(0.4) # correlation
                    )[[1]]
  sa_ratio_sensi1 <- round((std_size_sensi1$size - sa_size_sensi1) / std_size_sensi1$size, 4)
  sa_size_sensi1 
  sa_ratio_sensi1 
  
  ###### Sensitivity section 2 - sufficient EC (sample size 1000) ######
  sa_size_sensi2 <- size_sa_sensi2_true(start = 1,
                                        end = 500,
                                        alpha = 0.05,
                                        beta = 0.2,
                                        tau = 0.4)[[1]]
  sa_ratio_sensi2 <- round((std_size_sensi2$size - sa_size_sensi2) / std_size_sensi2$size, 4)
  sa_size_sensi2
  sa_ratio_sensi2
  
}

###########################################################################################
####### 4.2 Single arm case - non-informative tuning parameter and estimated models #######
####### only for sufficient EC case #######
{
  ### Functions ###
  {
    eval_sa_noninfo <- function(N_RCT, 
                                alpha, 
                                beta,
                                tau,
                                data_EC,
                                r0_M = 1, # V(Y|R=1, A=0) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                r1_M = 1, # V(Y|R=1, A=1) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                r = 1, # the ratio of variance in EC and variance in control of RCT
                                d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
                                dist_para = list(mu_R = 1,
                                                 sigma_R = 1,
                                                 p_R = 0.5,
                                                 mu_E = 1.2,
                                                 sigma_E = sqrt(1.5),
                                                 p_E = 0.7),
                                sample_size = 5*10^5,
                                gamma1 = 1, # the ratio of two kappa
                                gamma = 1 # correlation
    ){
      
      N_EC <- nrow(data_EC)
      n <- N_RCT + N_EC
      r_R <- N_RCT / N_EC
      
      # Step 1: Estimate V(Y|X, R=0, A=0) and V(Y|X, R=0, A=0)
      mu0_model <- lm(Y ~ X1 + X2, data_EC)
      var_X_EC <- mean((data_EC$Y - predict(mu0_model, data_EC))^2)
      var_EC <- (sd(data_EC$Y))^2
      
      # Step 2: estimate V(Y|R=1, A=1), V(Y|R=1, A=0) 
      var_RCT_0 <- var_EC * r0_M
      var_RCT_1 <- var_EC * r1_M
      
      # Step 3: estimate V(Y|X, R=1, A=0) using r(X)
      var_X_RCT_0 <- r * var_X_EC
      
      # Step 4: estimate kappa 0 and kappa 1
      kappa_0 <- mean(var_X_RCT_0)
      kappa_1 <- gamma1 * kappa_0
      
      # Step 5 - 7: d(X), gamma, variance
      if(length(d_X) < 1){
        
        pmf_bern <- function(x2, p){
          return(x2 * p + (1 - x2) * (1 - p))
        }
        
        d_X_true <- function(x1, x2, 
                             mu_R, sigma_R,
                             p_R,
                             mu_E, sigma_E,
                             p_E
        ){
          
          res <- dnorm(x1, mean = mu_R, sd = sigma_R) * pmf_bern(x2, p_R) / dnorm(x1, mean = mu_E, sd = sigma_E) / pmf_bern(x2, p_E) 
          return(res)
          
        }
        
        dX_x1_E <- rnorm(sample_size, mean = dist_para$mu_E, sd = dist_para$sigma_E)
        dX_x2_E <- rbinom(sample_size, size = 1, prob = dist_para$p_E)
        d_X_E <-  d_X_true(x1 = dX_x1_E,
                           x2 = dX_x2_E, 
                           dist_para$mu_R, dist_para$sigma_R, dist_para$p_R,
                           dist_para$mu_E, dist_para$sigma_E, dist_para$p_E)
        
        term1 <- kappa_1
        term3 <- (var_RCT_1 - kappa_1) + (var_RCT_0 - kappa_0) - 
          2 * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0))
        term4 <- mean((d_X_E * r_R)^2 / r_R  * var_X_EC)
        
      }else{
        term1 <- kappa_1
        term3 <- (var_RCT_1 - kappa_1) + (var_RCT_0 - kappa_0) - 
          2 * gamma * sqrt( (var_RCT_1 - kappa_1) * (var_RCT_0 - kappa_0) )
        term4 <- (d_X * r_R)^2 / r_R  * var_X_EC 
        
      }
      
      nu2 <- (term1 + term3 + term4)
      
      sqrt_nu2 <- sqrt(nu2)
      
      res <- 1 - beta - pnorm(qnorm(alpha/2) + sqrt(N_RCT) * tau / sqrt_nu2) -
        pnorm(qnorm(alpha/2) - sqrt(N_RCT) * tau / sqrt_nu2)
      
      return(res)
    }
    
    size_sa_noninfo <- function(start = 1,
                                end = 500,
                                alpha = 0.05,
                                beta = 0.2,
                                tau = 0.4,
                                data_EC,
                                r0_M = 1, # V(Y|R=1, A=0) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                r1_M = 1, # V(Y|R=1, A=1) / V(Y|R=0), 1.3/1.5 for sufficient and 1.3/1.585 for insufficient
                                r = 1, # the ratio of variance in EC and variance in control of RCT
                                d_X = 1, # the ratio of f(X|R=1) with f(X|R=0)
                                dist_para = list(mu_R = 1,
                                                 sigma_R = 1,
                                                 p_R = 0.5,
                                                 mu_E = 1.2,
                                                 sigma_E = sqrt(1.5),
                                                 p_E = 0.7),
                                sample_size = 5*10^5,
                                gamma1 = 1, # the ratio of two kappa
                                gamma = 1 # correlation
    ){
      
      size_length <- end - start + 1
      SA_size <- c()
      
      for(i in start:end){
        
        power_diff <- eval_sa_noninfo(N_RCT = i,
                                      alpha = alpha,
                                      beta = beta,
                                      tau = tau,
                                      data_EC = data_EC,
                                      r0_M = r0_M, 
                                      r1_M = r1_M, 
                                      r = r,
                                      d_X = d_X,
                                      dist_para = dist_para,
                                      sample_size = sample_size,
                                      gamma1 = gamma1, 
                                      gamma = gamma)
        
        SA_size <- rbind(SA_size, c(i, power_diff))
        if(power_diff <= 0){break}
      }
      
      return(list(i,
                  SA_size))
      
    }
    
    # Function for verifying the sufficient feasibility condition (under default setting)
    eval_sa_verify_noninfo <- function(alpha = 0.05, 
                                       beta = 0.2,
                                       tau = tau_use,
                                       N_EC = 60, 
                                       data_EC, 
                                       d_X = 1, 
                                       dist_para = NULL,
                                       sample_size =  5*10^5,
                                       gamma1 = 1, 
                                       gamma = 1 
    ){
      
      mu0_model <- lm(Y ~ X1 + X2, data_EC)
      var_X_EC <- mean((data_EC$Y - predict(mu0_model, data_EC))^2)
      
      # Step 6: estimate d(X)
      if(length(d_X) < 1){
        
        pmf_bern <- function(x2, p){
          return(x2 * p + (1 - x2) * (1 - p))
        }
        
        d_X_true <- function(x1, x2, 
                             mu_R, sigma_R,
                             p_R,
                             mu_E, sigma_E,
                             p_E
        ){
          
          res <- dnorm(x1, mean = mu_R, sd = sigma_R) * pmf_bern(x2, p_R) / dnorm(x1, mean = mu_E, sd = sigma_E) / pmf_bern(x2, p_E) 
          return(res)
          
        }
        
        dX_x1_E <- rnorm(sample_size, mean = dist_para$mu_E, sd = dist_para$sigma_E)
        dX_x2_E <- rbinom(sample_size, size = 1, prob = dist_para$p_E)
        d_X_E <-  d_X_true(x1 = dX_x1_E,
                           x2 = dX_x2_E, 
                           dist_para$mu_R, dist_para$sigma_R, dist_para$p_R,
                           dist_para$mu_E, dist_para$sigma_E, dist_para$p_E)
        
        condition <- mean((d_X_E)^2 * var_X_EC) * (qnorm(1 - beta) - qnorm(alpha/2))^2 / tau^2
        
      }else{
        
        condition <- (d_X)^2 * var_X_EC * (qnorm(1 - beta) - qnorm(alpha/2))^2 / tau^2
        
      }
      
      if(N_EC > condition){
        return(list(TRUE, condition))
      }else{
        return(list(FALSE, condition))
      }
      
    }
    
  }
  
  ### Use M = 2000 EC data to calculate the smallest sample size N_RCT 
  # Since for difference EC data, the calculated sample size is different
  # We generate 2000 EC data and calculate the average RCT sample size
  M <- 2000
  
  ###### Simulation section - sufficient EC (sample size 1000) ######
  {
   
    set.seed(987)
    Index_all <- sample(1:50000000, M, replace=F)
    sa_res_suff <- rep(NA, M)
    for(i in 1:M){
      print(i)
      data_EC <- f_gen_EC_suff(seed_EC = Index_all[i],
                               N_EC = 1000)
      sa_res_suff[i] <- size_sa_noninfo(start = 1,
                                        end = 500,
                                        alpha = 0.05,
                                        beta = 0.2,
                                        tau = 0.4,
                                        data_EC = data_EC)[[1]]
    }
    
    mean(sa_res_suff) |> ceiling()
    
  }
  
  ###### Sensitivity section 1 - sufficient EC (sample size 1000) ######
  # same as Simulation section - sufficient EC (sample size 1000)
  
  ###### Sensitivity section 2 - sufficient EC (sample size 1000) ######
  {
   
    set.seed(876)
    Index_all <- sample(1:50000000, M, replace=F)
    sa_res_sensi2 <- rep(NA, M)
    for(i in 1:M){
      print(i)
      data_EC <- f_gen_EC_sensi2(seed_EC = Index_all[i],
                                 N_EC = 1000)
      sa_res_sensi2[i] <- size_sa_noninfo(start = 1,
                                      end = 500,
                                      alpha = 0.05,
                                      beta = 0.2,
                                      tau = 0.4,
                                      data_EC = data_EC)[[1]]
    }
    
    mean(sa_res_sensi2) |> ceiling()
    
  }
  
}
