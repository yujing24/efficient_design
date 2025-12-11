library(nloptr)
library(dplyr)

setwd("/code/simulation")
source("4. simu_function.R")

n_list <- seq(6, 500, by = 1)
n_list_sensi <- seq(6, 600, by = 1)


### Simulation section - sufficient EC ###
{
  
  ### Type I error
  {
    # Type I error for RCT-only designs (difference-in-means, RCT-only AIPW), and hybrid design.
    Res05_typeI <- parallel::mclapply(n_list,
                                      f_assessment_suff, 
                                      mc.cores = 32,
                                      N_EC = 1000, 
                                      M = 2000, 
                                      alpha = 0.05, 
                                      tau = 0, # 0.4 for power
                                      seed = 123, 
                                      pi_A = 0.5 
    )
    saveRDS(Res05_typeI, "results/res_suff/res_05_TypeI_suff.RData")
    
    
    Res06_typeI <- parallel::mclapply(n_list,
                                      f_assessment_suff, 
                                      mc.cores = 32,
                                      N_EC = 1000, 
                                      M = 2000, 
                                      alpha = 0.05, 
                                      tau = 0, # 0.4 for power
                                      seed = 234, 
                                      pi_A = 0.6 
    )
    saveRDS(Res06_typeI, "results/res_suff/res_06_TypeI_suff.RData")
    
    Res07_typeI <- parallel::mclapply(n_list,
                                      f_assessment_suff, 
                                      mc.cores = 32,
                                      N_EC = 1000, 
                                      M = 2000, 
                                      alpha = 0.05, 
                                      tau = 0, # 0.4 for power
                                      seed = 345, 
                                      pi_A = 0.7 
    )
    saveRDS(Res07_typeI, "results/res_suff/res_07_TypeI_suff.RData")
    
    
    Res08_typeI <- parallel::mclapply(n_list,
                                      f_assessment_suff, 
                                      mc.cores = 32,
                                      N_EC = 1000, 
                                      M = 2000, 
                                      alpha = 0.05, 
                                      tau = 0, # 0.4 for power
                                      seed = 456, 
                                      pi_A = 0.8
    )
    saveRDS(Res08_typeI, "results/res_suff/res_08_TypeI_suff.RData")
    
    Res09_typeI <- parallel::mclapply(n_list,
                                      f_assessment_suff, 
                                      mc.cores = 32,
                                      N_EC = 1000, 
                                      M = 2000, 
                                      alpha = 0.05, 
                                      tau = 0, # 0.4 for power
                                      seed = 567, 
                                      pi_A = 0.9 
    )
    saveRDS(Res09_typeI, "results/res_suff/res_09_TypeI_suff.RData")
    
    # Type I error for single-arm design. 
    Res_sa_typeI <- parallel::mclapply(n_list,
                                       f_assessment_sa, 
                                       mc.cores = 32,
                                       N_EC = 1000, 
                                       M = 2000, 
                                       alpha = 0.05, 
                                       tau = 0, # 0.4 for power
                                       seed = 122, 
                                       pi_A = 1 
    )
    saveRDS(Res_sa_typeI, "results/res_suff/res_sa_TypeI_suff.RData")
  }
  
  ### Power
  {
    # Power for RCT-only designs (difference-in-means, RCT-only AIPW), and hybrid design.
    Res05_power <- parallel::mclapply(n_list,
                                      f_assessment_suff, 
                                      mc.cores = 32,
                                      N_EC = 1000, 
                                      M = 2000, 
                                      alpha = 0.05, 
                                      tau = 0.4, 
                                      seed = 123, 
                                      pi_A = 0.5 
    )
    saveRDS(Res05_power, "results/res_suff/res_05_power_suff.RData")
    
    
    Res06_power <- parallel::mclapply(n_list,
                                      f_assessment_suff, 
                                      mc.cores = 32,
                                      N_EC = 1000, 
                                      M = 2000, 
                                      alpha = 0.05, 
                                      tau = 0.4, 
                                      seed = 234, 
                                      pi_A = 0.6 
    )
    saveRDS(Res06_power, "results/res_suff/res_06_power_suff.RData")
    
    
    Res07_power <- parallel::mclapply(n_list,
                                      f_assessment_suff, 
                                      mc.cores = 32,
                                      N_EC = 1000, 
                                      M = 2000, 
                                      alpha = 0.05, 
                                      tau = 0.4, 
                                      seed = 345, 
                                      pi_A = 0.7 
    )
    saveRDS(Res07_power, "results/res_suff/res_07_power_suff.RData")
    
    
    Res08_power <- parallel::mclapply(n_list,
                                      f_assessment_suff, 
                                      mc.cores = 32,
                                      N_EC = 1000, 
                                      M = 2000, 
                                      alpha = 0.05, 
                                      tau = 0.4, 
                                      seed = 456, 
                                      pi_A = 0.8
    )
    saveRDS(Res08_power, "results/res_suff/res_08_power_suff.RData")
    
    
    Res09_power <- parallel::mclapply(n_list,
                                      f_assessment_suff, 
                                      mc.cores = 32,
                                      N_EC = 1000, 
                                      M = 2000, 
                                      alpha = 0.05, 
                                      tau = 0.4, 
                                      seed = 567, 
                                      pi_A = 0.9 
    )
    saveRDS(Res09_power, "results/res_suff/res_09_power_suff.RData")
    
    
    # Power for single-arm design. 
    Res_sa_power <- parallel::mclapply(n_list,
                                       f_assessment_sa, 
                                       mc.cores = 32,
                                       N_EC = 1000, 
                                       M = 2000, 
                                       alpha = 0.05, 
                                       tau = 0.4, 
                                       seed = 122, 
                                       pi_A = 1 
    )
    saveRDS(Res_sa_power, "results/res_suff/res_sa_power_suff.RData")
  }
  
  
}

### Simulation section - insufficient EC (60) ###
{
  ### Type I error
  {
    
    # Type I error for RCT-only designs (difference-in-means, RCT-only AIPW), and hybrid design.
    Res05_typeI_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 60, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 123, 
                                             pi_A = 0.5 
    )
    saveRDS(Res05_typeI_insuff, "results/res_insuff_60/res_05_TypeI_insuff.RData")
    
    
    Res06_typeI_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 60, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 234, 
                                             pi_A = 0.6 
    )
    saveRDS(Res06_typeI_insuff, "results/res_insuff_60/res_06_TypeI_insuff.RData")
    
    
    Res07_typeI_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 60, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 345, 
                                             pi_A = 0.7 
    )
    saveRDS(Res07_typeI_insuff, "results/res_insuff_60/res_07_TypeI_insuff.RData")
    
    
    Res08_typeI_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 60, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 456, 
                                             pi_A = 0.8
    )
    saveRDS(Res08_typeI_insuff, "results/res_insuff_60/res_08_TypeI_insuff.RData")
    
    
    Res09_typeI_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 60, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 567, 
                                             pi_A = 0.9 
    )
    saveRDS(Res09_typeI_insuff, "results/res_insuff_60/res_09_TypeI_insuff.RData")
    
  }
  
  ### Power
  {
    
    # Power for RCT-only designs (difference-in-means, RCT-only AIPW), and hybrid design.
    Res05_power_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 60, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 123, 
                                             pi_A = 0.5 
    )
    saveRDS(Res05_power_insuff, "results/res_insuff_60/res_05_power_insuff.RData")
    
    
    Res06_power_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 60, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 234, 
                                             pi_A = 0.6 
    )
    saveRDS(Res06_power_insuff, "results/res_insuff_60/res_06_power_insuff.RData")
    
    
    Res07_power_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 60, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 345, 
                                             pi_A = 0.7 
    )
    saveRDS(Res07_power_insuff, "results/res_insuff_60/res_07_power_insuff.RData")
    
    
    Res08_power_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 60, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 456, 
                                             pi_A = 0.8
    )
    saveRDS(Res08_power_insuff, "results/res_insuff_60/res_08_power_insuff.RData")
    
    Res09_power_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 60, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 567, 
                                             pi_A = 0.9 
    )
    saveRDS(Res09_power_insuff, "results/res_insuff_60/res_09_power_insuff.RData")
    
  }
}

### Simulation section - insufficient EC (30) ###
{
  ### Type I error
  {
    
    # Type I error for RCT-only designs (difference-in-means, RCT-only AIPW), and hybrid design.
    Res05_typeI_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 30, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 123, 
                                             pi_A = 0.5 
    )
    saveRDS(Res05_typeI_insuff, "results/res_insuff_30/res_05_TypeI_insuff.RData")
    
    
    Res06_typeI_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 30, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 234, 
                                             pi_A = 0.6 
    )
    saveRDS(Res06_typeI_insuff, "results/res_insuff_30/res_06_TypeI_insuff.RData")
    
    
    Res07_typeI_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 30, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 345, 
                                             pi_A = 0.7 
    )
    saveRDS(Res07_typeI_insuff, "results/res_insuff_30/res_07_TypeI_insuff.RData")
    
    
    Res08_typeI_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 30, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 456, 
                                             pi_A = 0.8
    )
    saveRDS(Res08_typeI_insuff, "results/res_insuff_30/res_08_TypeI_insuff.RData")
    
    
    Res09_typeI_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 30, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 567, 
                                             pi_A = 0.9 
    )
    saveRDS(Res09_typeI_insuff, "results/res_insuff_30/res_09_TypeI_insuff.RData")
    
  }
  
  ### Power
  {
    
    # Power for RCT-only designs (difference-in-means, RCT-only AIPW), and hybrid design.
    Res05_power_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 30, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 123, 
                                             pi_A = 0.5 
    )
    saveRDS(Res05_power_insuff, "results/res_insuff_30/res_05_power_insuff.RData")
    
    
    Res06_power_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 30, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 234, 
                                             pi_A = 0.6 
    )
    saveRDS(Res06_power_insuff, "results/res_insuff_30/res_06_power_insuff.RData")
    
    
    Res07_power_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 30, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 345, 
                                             pi_A = 0.7 
    )
    saveRDS(Res07_power_insuff, "results/res_insuff_30/res_07_power_insuff.RData")
    
    
    Res08_power_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 30, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 456, 
                                             pi_A = 0.8
    )
    saveRDS(Res08_power_insuff, "results/res_insuff_30/res_08_power_insuff.RData")
    
    
    Res09_power_insuff <- parallel::mclapply(n_list,
                                             f_assessment_insuff, 
                                             mc.cores = 32,
                                             N_EC = 30, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 567, 
                                             pi_A = 0.9 
    )
    saveRDS(Res09_power_insuff, "results/res_insuff_30/res_09_power_insuff.RData")
    
  }
}

### Sensitivity analysis 1 - sufficient EC ###
{
  ### Type I error
  {
    
    # Type I error for RCT-only designs (difference-in-means, RCT-only AIPW), and hybrid design.
    Res05_typeI_sensi1 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi1, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau1 = 0, 
                                             tau2 = 0, 
                                             seed = 123, 
                                             pi_A = 0.5 
    )
    saveRDS(Res05_typeI_sensi1, "results/res_sensi1/res_05_TypeI_sensi1.RData")
    
    
    Res06_typeI_sensi1 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi1, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau1 = 0, 
                                             tau2 = 0, 
                                             seed = 234, 
                                             pi_A = 0.6 
    )
    saveRDS(Res06_typeI_sensi1, "results/res_sensi1/res_06_TypeI_sensi1.RData")
    
    
    Res07_typeI_sensi1 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi1, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau1 = 0, 
                                             tau2 = 0, 
                                             seed = 345, 
                                             pi_A = 0.7 
    )
    saveRDS(Res07_typeI_sensi1, "results/res_sensi1/res_07_TypeI_sensi1.RData")
    
    
    Res08_typeI_sensi1 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi1, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau1 = 0, 
                                             tau2 = 0, 
                                             seed = 456, 
                                             pi_A = 0.8
    )
    saveRDS(Res08_typeI_sensi1, "results/res_sensi1/res_08_TypeI_sensi1.RData")
    
    
    Res09_typeI_sensi1 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi1, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau1 = 0, 
                                             tau2 = 0, 
                                             seed = 567, 
                                             pi_A = 0.9 
    )
    saveRDS(Res09_typeI_sensi1, "results/res_sensi1/res_09_TypeI_sensi1.RData")
    
    
    # Type I error for single-arm design. 
    Res_sa_typeI_sensi1 <- parallel::mclapply(n_list_sensi,
                                              f_assessment_sa_sensi1, 
                                              mc.cores = 32,
                                              N_EC = 1000, 
                                              M = 2000, 
                                              alpha = 0.05, 
                                              tau1 = 0, 
                                              tau2 = 0, 
                                              seed = 122, 
                                              pi_A = 1 
    )
    saveRDS(Res_sa_typeI_sensi1, "results/res_sensi1/res_sa_TypeI_sensi1.RData")
    
  }
  
  ### Power
  {
    
    # Power for RCT-only designs (difference-in-means, RCT-only AIPW), and hybrid design.
    Res05_power_sensi1 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi1, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau1 = 0.3, 
                                             tau2 = 0.2, 
                                             seed = 123, 
                                             pi_A = 0.5 
    )
    saveRDS(Res05_power_sensi1, "results/res_sensi1/res_05_power_sensi1.RData")
    
    
    Res06_power_sensi1 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi1, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau1 = 0.3, 
                                             tau2 = 0.2,
                                             seed = 234, 
                                             pi_A = 0.6 
    )
    saveRDS(Res06_power_sensi1, "results/res_sensi1/res_06_power_sensi1.RData")
    
    
    Res07_power_sensi1 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi1, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau1 = 0.3, 
                                             tau2 = 0.2,
                                             seed = 345, 
                                             pi_A = 0.7 
    )
    saveRDS(Res07_power_sensi1, "results/res_sensi1/res_07_power_sensi1.RData")
    
    
    Res08_power_sensi1 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi1, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau1 = 0.3, 
                                             tau2 = 0.2,
                                             seed = 456, 
                                             pi_A = 0.8
    )
    saveRDS(Res08_power_sensi1, "results/res_sensi1/res_08_power_sensi1.RData")
    
    
    Res09_power_sensi1 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi1, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau1 = 0.3, 
                                             tau2 = 0.2,
                                             seed = 567, 
                                             pi_A = 0.9 
    )
    saveRDS(Res09_power_sensi1, "results/res_sensi1/res_09_power_sensi1.RData")
    
    
    # Power for single-arm design. 
    Res_sa_power_sensi1 <- parallel::mclapply(n_list_sensi,
                                              f_assessment_sa_sensi1, 
                                              mc.cores = 32,
                                              N_EC = 1000, 
                                              M = 2000, 
                                              alpha = 0.05, 
                                              tau1 = 0.3, 
                                              tau2 = 0.2,
                                              seed = 122, 
                                              pi_A = 1 
    )
    saveRDS(Res_sa_power_sensi1, "results/res_sensi1/res_sa_power_sensi1.RData")
    
  }
  
}

### Sensitivity analysis 2 - sufficient EC ###
{
  ### Type I error
  {
    # Type I error for RCT-only designs (difference-in-means, RCT-only AIPW), and hybrid design.
    Res05_typeI_sensi2 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi2, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 123, 
                                             pi_A = 0.5 
    )
    saveRDS(Res05_typeI_sensi2, "results/res_sensi2/res_05_TypeI_sensi2.RData")
    
    
    Res06_typeI_sensi2 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi2, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 234, 
                                             pi_A = 0.6 
    )
    saveRDS(Res06_typeI_sensi2, "results/res_sensi2/res_06_TypeI_sensi2.RData")
    
    Res07_typeI_sensi2 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi2, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 345, 
                                             pi_A = 0.7 
    )
    saveRDS(Res07_typeI_sensi2, "results/res_sensi2/res_07_TypeI_sensi2.RData")
    
    
    Res08_typeI_sensi2 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi2, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 456, 
                                             pi_A = 0.8
    )
    saveRDS(Res08_typeI_sensi2, "results/res_sensi2/res_08_TypeI_sensi2.RData")
    
    Res09_typeI_sensi2 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi2, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0, # 0.4 for power
                                             seed = 567, 
                                             pi_A = 0.9 
    )
    saveRDS(Res09_typeI_sensi2, "results/res_sensi2/res_09_TypeI_sensi2.RData")
    
    # Type I error for single-arm design. 
    Res_sa_typeI_sensi2 <- parallel::mclapply(n_list_sensi,
                                              f_assessment_sa_sensi2, 
                                              mc.cores = 32,
                                              N_EC = 1000, 
                                              M = 2000, 
                                              alpha = 0.05, 
                                              tau = 0, # 0.4 for power
                                              seed = 122, 
                                              pi_A = 1 
    )
    saveRDS(Res_sa_typeI_sensi2, "results/res_sensi2/res_sa_TypeI_sensi2.RData")
  }
  
  ### Power
  {
    # Power for RCT-only designs (difference-in-means, RCT-only AIPW), and hybrid design.
    Res05_power_sensi2 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi2, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 123, 
                                             pi_A = 0.5 
    )
    saveRDS(Res05_power_sensi2, "results/res_sensi2/res_05_power_sensi2.RData")
    
    
    Res06_power_sensi2 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi2, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 234, 
                                             pi_A = 0.6 
    )
    saveRDS(Res06_power_sensi2, "results/res_sensi2/res_06_power_sensi2.RData")
    
    Res07_power_sensi2 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi2, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 345, 
                                             pi_A = 0.7 
    )
    saveRDS(Res07_power_sensi2, "results/res_sensi2/res_07_power_sensi2.RData")
    
    
    Res08_power_sensi2 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi2, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 456, 
                                             pi_A = 0.8
    )
    saveRDS(Res08_power_sensi2, "results/res_sensi2/res_08_power_sensi2.RData")
    
    Res09_power_sensi2 <- parallel::mclapply(n_list_sensi,
                                             f_assessment_sensi2, 
                                             mc.cores = 32,
                                             N_EC = 1000, 
                                             M = 2000, 
                                             alpha = 0.05, 
                                             tau = 0.4, 
                                             seed = 567, 
                                             pi_A = 0.9 
    )
    saveRDS(Res09_power_sensi2, "results/res_sensi2/res_09_power_sensi2.RData")
    
    # Power for single-arm design. 
    Res_sa_power_sensi2 <- parallel::mclapply(n_list_sensi,
                                              f_assessment_sa_sensi2, 
                                              mc.cores = 32,
                                              N_EC = 1000, 
                                              M = 2000, 
                                              alpha = 0.05, 
                                              tau = 0.4, 
                                              seed = 122, 
                                              pi_A = 1 
    )
    saveRDS(Res_sa_power_sensi2, "results/res_sensi2/res_sa_power_sensi2.RData")
  }
}
