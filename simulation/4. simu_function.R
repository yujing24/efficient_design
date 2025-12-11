library(nloptr)
library(dplyr)

setwd("/code/simulation")
source("3. est_function.R")
source("4. data_gen_function.R")

## Assessment functions for simulation section in the main paper
{
  # Function for assessment in sufficient EC case (difference-in-means, RCT-only AIPW and hybrid design)
  f_assessment_suff <- function(n, # RCT sample size
                                N_EC = 1000, # EC data sample size 
                                M = 2000, # simulation times
                                alpha = 0.05, # type I error for two-sided hypothesis test
                                tau = 0, # true ATT, 0 for type I error and 0.4 for power 
                                seed = 123, # random seed for generating trial data
                                pi_A = 0.5 # treatment allocation in RCT data
  ){
    
    N <- N_EC + n # total size 
    
    set.seed(seed)
    Index_all <- sample(1:50000000, M, replace=F)
    
    Result <- c()
    for(i in 1:M){
      data_EC <- f_gen_EC(seed_EC = Index_all[i],
                          N_EC = N_EC)
      data_RCT <- f_gen_RCT(seed_RCT = Index_all[i],
                            N_RCT = n,
                            pi_A = pi_A,
                            tau = tau)
      data_full <-  f_gen_full(data_EC, data_RCT)
      
      res_EC <-  f_est_hybrid(data = data_full)
      res_AIPW <- f_est_AIPW(data = data_RCT)
      res_std <- f_est_naive(data = data_RCT)
      
      res <- cbind(res_EC, res_AIPW, res_std)
      
      Result <- rbind(Result, res)
    }
    
    Z_IF_EC <- (Result$Est_EC - 0) / sqrt(Result$Est_EC_var / N)
    res_IF_EC <- (abs(Z_IF_EC) > qnorm(1 - alpha/2)) * 1
    
    Z_IF_AIPW <- (Result$Est_AIPW - 0) / sqrt(Result$Est_AIPW_var / n)
    res_IF_AIPW <- (abs(Z_IF_AIPW) > qnorm(1 - alpha/2)) * 1
    
    res_naive <- (Result$p_value < 0.05) *1
    
    res <- data.frame(power_EC = sum(res_IF_EC) / M,
                      power_AIPW = sum(res_IF_AIPW) / M,
                      power_naive =  sum(res_naive) / M)
    
    return(res)
    
  }
  
  # Function for assessment in sufficient EC case (single-arm design)
  f_assessment_sa <- function(n, # RCT sample size
                              N_EC = 1000, # EC data sample size 
                              M = 2000, # simulation times
                              alpha = 0.05, # type I error for two-sided hypothesis test
                              tau = 0, # true ATT, 0 for type I error and 0.4 for power 
                              seed = 122, # random seed for generating trial data
                              pi_A = 1 # treatment allocation in RCT data
  ){
    
    N <- N_EC + n # total size 
    
    set.seed(seed)
    Index_all <- sample(1:50000000, M, replace=F)
    
    Result <- c()
    for(i in 1:M){
      data_EC <- f_gen_EC(seed_EC = Index_all[i],
                          N_EC = N_EC)
      data_RCT <- f_gen_RCT(seed_RCT = Index_all[i],
                            N_RCT = n,
                            pi_A = pi_A,
                            tau = tau)
      data_full <-  f_gen_full(data_EC, data_RCT)
      
      res_SA <- f_est_sa(data_full)
      
      res <- res_SA
      
      Result <- rbind(Result, res)
    }
    
    Z_IF_SA <- (Result$Est_SA - 0) / sqrt(Result$Est_SA_var / N)
    res_IF_SA <- (abs(Z_IF_SA) > qnorm(1 - alpha/2)) * 1
    
    res <- data.frame(power_SA = sum(res_IF_SA) / M)
    
    return(res)
    
  }
  
  # Function for assessment in insufficient EC case (difference-in-means, RCT-only AIPW and hybrid design)
  f_assessment_insuff <- function(n, # RCT sample size
                                  N_EC = 60, # EC data sample size 
                                  M = 2000, # simulation times
                                  alpha = 0.05, # type I error for two-sided hypothesis test
                                  tau = 0, # true ATT, 0 for type I error and 0.4 for power 
                                  seed = 123, # random seed for generating trial data
                                  pi_A = 0.5 # treatment allocation in RCT data
  ){
    
    N <- N_EC + n # total size 
    
    set.seed(seed)
    Index_all <- sample(1:50000000, M, replace=F)
    
    Result <- c()
    for(i in 1:M){
      data_EC <- f_gen_EC(seed_EC = Index_all[i],
                          N_EC = N_EC,
                          muX1_EC = 1.2, 
                          sdX1_EC = sqrt(1.5),
                          pX2_EC = 0.7,
                          beta0_EC = 1,
                          beta1_EC = 0.5,
                          beta2_EC = -1,
                          sdY_EC = 1)
      data_RCT <- f_gen_RCT(seed_RCT = Index_all[i],
                            N_RCT = n,
                            pi_A = pi_A,
                            tau = tau,
                            muX1_RCT = 1, 
                            sdX1_RCT = 1,
                            pX2_RCT = 0.5,
                            beta0_RCT = 1,
                            beta1_RCT = 0.5,
                            beta2_RCT = -1,
                            sdY_RCT = sqrt(0.8)
      )
      data_full <-  f_gen_full(data_EC, data_RCT)
      
      res_EC <-  f_est_hybrid(data = data_full)
      res_AIPW <- f_est_AIPW(data = data_RCT)
      res_std <- f_est_naive(data = data_RCT)
      
      res <- cbind(res_EC, res_AIPW, res_std)
      
      Result <- rbind(Result, res)
    }
    
    Z_IF_EC <- (Result$Est_EC - 0) / sqrt(Result$Est_EC_var / N)
    res_IF_EC <- (abs(Z_IF_EC) > qnorm(1 - alpha/2)) * 1
    
    Z_IF_AIPW <- (Result$Est_AIPW - 0) / sqrt(Result$Est_AIPW_var / n)
    res_IF_AIPW <- (abs(Z_IF_AIPW) > qnorm(1 - alpha/2)) * 1
    
    res_naive <- (Result$p_value < 0.05) *1
    
    res <- data.frame(power_EC = sum(res_IF_EC) / M,
                      power_AIPW = sum(res_IF_AIPW) / M,
                      power_naive =  sum(res_naive) / M)
    
    return(res)
    
  }
  
}


## Assessment functions for sensitivity analysis 1
{
  f_assessment_sensi1 <- function(n, # RCT sample size
                                  N_EC = 1000, # EC data sample size 
                                  M = 2000, # simulation times
                                  alpha = 0.05, # type I error for two-sided hypothesis test
                                  tau1 = 0, 
                                  tau2 = 0, # tau1 = 0.3, tau2 = 0.2 for power analysis
                                  seed = 123, # random seed for generating trial data
                                  pi_A = 0.5 # treatment allocation in RCT data
  ){
    
    N <- N_EC + n
    
    set.seed(seed)
    Index_all <- sample(1:50000000, M, replace=F)
    
    Result <- c()
    for(i in 1:M){
      data_EC <- f_gen_EC(seed_EC = Index_all[i],
                          N_EC = N_EC)
      data_RCT <- f_gen_RCT_sensi1(seed_RCT = Index_all[i],
                                   N_RCT = n,
                                   pi_A = pi_A,
                                   tau1 = tau1,
                                   tau2 = tau2)
      data_full <-  f_gen_full(data_EC, data_RCT)
      
      res_EC <- f_est_hybrid(data_full)
      res_AIPW <- f_est_AIPW(data_RCT)
      res_std <- f_est_naive(data = data_RCT)
      
      res <- cbind(res_EC, res_AIPW, res_std)
      
      Result <- rbind(Result, res)
    }
    
    Z_IF_EC <- (Result$Est_EC - 0) / sqrt(Result$Est_EC_var / N)
    res_IF_EC <- (abs(Z_IF_EC) > qnorm(1 - alpha/2)) * 1
    
    Z_IF_AIPW <- (Result$Est_AIPW - 0) / sqrt(Result$Est_AIPW_var / n)
    res_IF_AIPW <- (abs(Z_IF_AIPW) > qnorm(1 - alpha/2)) * 1
    
    res_naive <- (Result$p_value < 0.05) *1
    
    res <- data.frame(power_EC = sum(res_IF_EC) / M,
                      power_AIPW = sum(res_IF_AIPW) / M,
                      power_naive =  sum(res_naive) / M)
    
    return(res)
    
  }
  
  f_assessment_sa_sensi1 <- function(n, # RCT sample size
                                     N_EC = 1000, # EC data sample size 
                                     M = 2000, # simulation times
                                     alpha = 0.05, # type I error for two-sided hypothesis test
                                     tau1 = 0, 
                                     tau2 = 0, # tau1 = 0.3, tau2 = 0.2 for power analysis
                                     seed = 122, # random seed for generating trial data
                                     pi_A = 1 # treatment allocation in RCT data
  ){
    
    N <- N_EC + n
    
    set.seed(seed)
    Index_all <- sample(1:50000000, M, replace=F)
    
    Result <- c()
    for(i in 1:M){
      data_EC <- f_gen_EC(seed_EC = Index_all[i],
                          N_EC = N_EC)
      data_RCT <- f_gen_RCT_sensi1(seed_RCT = Index_all[i],
                                   N_RCT = n,
                                   pi_A = pi_A,
                                   tau1 = tau1,
                                   tau2 = tau2)
      data_full <- f_gen_full(data_EC, data_RCT)
      
      res_SA <- f_est_sa(data_full)
      
      res <- res_SA
      
      Result <- rbind(Result, res)
      
    }
    
    Z_IF_SA <- (Result$Est_SA - 0) / sqrt(Result$Est_SA_var / N)
    res_IF_SA <- (abs(Z_IF_SA) > qnorm(1 - alpha/2)) * 1
    
    res <- data.frame(power_SA = sum(res_IF_SA) / M)
    
    return(res)
    
    
  }
}


## Assessment functions for sensitivity analysis 2
{
  
  # Function for assessment in sufficient EC case (difference-in-means, RCT-only AIPW and hybrid design)
  f_assessment_sensi2 <- function(n, # RCT sample size
                                N_EC = 1000, # EC data sample size 
                                M = 2000, # simulation times
                                alpha = 0.05, # type I error for two-sided hypothesis test
                                tau = 0, # true ATT, 0 for type I error and 0.4 for power 
                                seed = 123, # random seed for generating trial data
                                pi_A = 0.5 # treatment allocation in RCT data
  ){
    
    N <- N_EC + n # total size 
    
    set.seed(seed)
    Index_all <- sample(1:50000000, M, replace=F)
    
    Result <- c()
    for(i in 1:M){
      data_EC <- f_gen_EC_sensi2(seed_EC = Index_all[i],
                                  N_EC = N_EC)
      data_RCT <- f_gen_RCT_sensi2(seed_RCT = Index_all[i],
                                  N_RCT = n,
                                  pi_A = pi_A,
                                  tau = tau)
      data_full <-  f_gen_full(data_EC, data_RCT)
      
      res_EC <-  f_est_hybrid(data = data_full)
      res_AIPW <- f_est_AIPW(data = data_RCT)
      res_std <- f_est_naive(data = data_RCT)
      
      res <- cbind(res_EC, res_AIPW, res_std)
      
      Result <- rbind(Result, res)
    }
    
    Z_IF_EC <- (Result$Est_EC - 0) / sqrt(Result$Est_EC_var / N)
    res_IF_EC <- (abs(Z_IF_EC) > qnorm(1 - alpha/2)) * 1
    
    Z_IF_AIPW <- (Result$Est_AIPW - 0) / sqrt(Result$Est_AIPW_var / n)
    res_IF_AIPW <- (abs(Z_IF_AIPW) > qnorm(1 - alpha/2)) * 1
    
    res_naive <- (Result$p_value < 0.05) *1
    
    res <- data.frame(power_EC = sum(res_IF_EC) / M,
                      power_AIPW = sum(res_IF_AIPW) / M,
                      power_naive =  sum(res_naive) / M)
    
    return(res)
    
  }
  
  # Function for assessment in sufficient EC case (single-arm design)
  f_assessment_sa_sensi2 <- function(n, # RCT sample size
                              N_EC = 1000, # EC data sample size 
                              M = 2000, # simulation times
                              alpha = 0.05, # type I error for two-sided hypothesis test
                              tau = 0, # true ATT, 0 for type I error and 0.4 for power 
                              seed = 122, # random seed for generating trial data
                              pi_A = 1 # treatment allocation in RCT data
  ){
    
    N <- N_EC + n # total size 
    
    set.seed(seed)
    Index_all <- sample(1:50000000, M, replace=F)
    
    Result <- c()
    for(i in 1:M){
      data_EC <- f_gen_EC_sensi2(seed_EC = Index_all[i],
                                N_EC = N_EC)
      data_RCT <- f_gen_RCT_sensi2(seed_RCT = Index_all[i],
                                    N_RCT = n,
                                    pi_A = pi_A,
                                    tau = tau)
      data_full <-  f_gen_full(data_EC, data_RCT)
      
      res_SA <- f_est_sa(data_full)
      
      res <- res_SA
      
      Result <- rbind(Result, res)
    }
    
    Z_IF_SA <- (Result$Est_SA - 0) / sqrt(Result$Est_SA_var / N)
    res_IF_SA <- (abs(Z_IF_SA) > qnorm(1 - alpha/2)) * 1
    
    res <- data.frame(power_SA = sum(res_IF_SA) / M)
    
    return(res)
    
  }
  
}
