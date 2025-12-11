rm(list=ls())
library(dplyr)
library(MASS)
library(haven)
library(nloptr)

setwd("/code/case_study")

data <- read_sas("chapter15_example.sas7bdat")
data_trt <- filter(data, THERAPY == "DRUG")  
data_crl <- filter(data, THERAPY == "PLACEBO")

### Get the week 4's outcome
{
  data_4 <-  filter(data, VISIT == 4)
  data_4$X1 <- data_4$basval
  data_4$X2 <- ifelse(data_4$GENDER == "F", 1, 0)
  data_4$Y <- data_4$change
  
  data_4_trt <- filter(data_4, THERAPY == "DRUG") # size = 84
  data_4_crl <- filter(data_4, THERAPY == "PLACEBO") # size = 88
  
  hist(data_4_trt$change)
  hist(data_4_crl$change)
  
  mean(data_4_trt$change) - mean(data_4_crl$change) 
  
  data_4_crl$GENDER <- as.factor(data_4_crl$GENDER)
  data_4_crl$POOLINV <- as.factor(data_4_crl$POOLINV)
  
}

###  Trial design  
## Stage 1: Generate one bootstrap control data as the EC data before trial design
{
  
  N_4_crl <- nrow(data_4_crl)
  
  ## Draw data from data_4_crl as the EC data (sufficient 1000 or insufficient 60)
  # set.seed(123)
  # EC_index_suff <- sample(N_4_crl, size = 1000, replace = T)
  # EC_4_suff <- data_4_crl[EC_index_suff, ]
  # saveRDS(EC_4_suff, "sufficient_EC.RData")
  EC_4_suff <- readRDS("sufficient_EC.RData")
    
  # set.seed(234)
  # EC_index_insuff <- sample(N_4_crl, size = 60, replace = T)
  # EC_4_insuff <- data_4_crl[EC_index_insuff, ]
  # saveRDS(EC_4_insuff, "insufficient_EC.RData")
  EC_4_insuff <- readRDS("insufficient_EC.RData")
  
  # RCT data:  data_4_trt and data_4_crl
  
  # pi_A in the original data
  N_4_trt <- nrow(data_4_trt)
  pi_A_orig <- N_4_trt / (N_4_trt + N_4_crl)
  pi_A_orig
  
  # delta: treatment effect
  delta2 <- -1
  # target effect size
  tau_use <- delta2
}

## Stage 2: Calculate the sample size
{
  
  ## 1.1 RCT-only with difference-in-means - with sufficient EC data (size = 1000)
  {
    
    # sample size calculation 
    sample_size <- function(pi_A = pi_A_orig, 
                            alpha = 0.05,
                            beta = 0.2,
                            var_RCT_0,
                            var_RCT_1,
                            tau){
      ratio <- pi_A / (1 - pi_A) # n2/n1 = control : treatment
      zalpha <- qnorm(alpha/2)
      zbeta <- qnorm(1 - beta)
      N_T <- (var_RCT_1 + ratio * var_RCT_0) * (zbeta - zalpha)^2 / tau^2
      return(data.frame(N_T, N_C = N_T / ratio))
    }
    
    sample_size(pi_A = pi_A_orig, # n2/n1 = control : treatment
                alpha = 0.05,
                beta = 0.2,
                var_RCT_0 = var(EC_4_suff$Y), # treat
                var_RCT_1 = var(EC_4_suff$Y), # control
                tau = tau_use) |> ceiling() |> sum()
    # 458
    
  }
  
  ## 1.2 RCT-only with AIPW - with sufficient EC data (size = 1000)
  {
    
    eval_aipw_apl <- function(n,
                              alpha, 
                              beta,
                              pi_A, 
                              tau,
                              data_EC = EC_4_suff,
                              r0_M = 1, 
                              r1_M = 1,
                              r = 1,
                              gamma1 = 1,
                              gamma = 1
    ){
      
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
      
      term1 <- (((1 - pi_A) / pi_A) * kappa_1 + var_RCT_1)
      term2 <- ((pi_A / (1 - pi_A)) * kappa_0 + var_RCT_0)  
      term3 <- -2 * gamma * sqrt((var_RCT_0 - kappa_0) * (var_RCT_1 - kappa_1))
      
      nu2 <- (term1 + term2 + term3)
      
      sqrt_nu2 <- sqrt(nu2)
      
      res <- 1 - beta - pnorm(qnorm(alpha/2) + sqrt(n) * tau / sqrt_nu2) -
        pnorm(qnorm(alpha/2) - sqrt(n) * tau / sqrt_nu2)
      
      return(res)
    }
    
    size_aipw_apl <- function(start = 1,
                              end = 500,
                              alpha = 0.05,
                              beta = 0.2,
                              pi_A = pi_A_orig,
                              tau = tau_use,
                              data_EC = EC_4_suff,
                              r0_M = 1, 
                              r1_M = 1,
                              r = 1,
                              gamma1 = 1,
                              gamma = 1){
      
      
      size_length <- end - start + 1
      AIPW_size <- matrix(NA, nrow = size_length, ncol = 2)
      
      for(i in start:end){
        power_diff <- eval_aipw_apl(n = i,
                                    alpha = alpha,
                                    beta = beta,
                                    pi_A = pi_A,
                                    tau = tau,
                                    data_EC = data_EC,
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
    
    size_aipw_apl(data_EC = EC_4_suff)[[1]]
    
    # 421
  }
  
  ## 1.3 hybrid experimental design - with sufficient EC data (size = 1000)
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
                                    dist_para = NULL,
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
                                    tau = tau_use,
                                    pi_A = 0.5,
                                    data_EC, 
                                    r0_M = 1, 
                                    r1_M = 1,
                                    r = 1,
                                    d_X = 1,
                                    dist_para = NULL,
                                    sample_size = 5*10^5,
                                    gamma1 = 1, 
                                    gamma = 1 
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
    
    size_hybrid_noninfo(start = 1,
                        end = 500,
                        alpha = 0.05,
                        beta = 0.2,
                        tau = tau_use,
                        pi_A = pi_A_orig,
                        data_EC = EC_4_suff)[[1]]
    
    # 238
  }
  
  ## 1.4 hybrid experimental design - with insufficient EC data (size = 60)
  {
    
    size_hybrid_noninfo(start = 1,
                        end = 500,
                        alpha = 0.05,
                        beta = 0.2,
                        tau = tau_use,
                        pi_A = pi_A_orig,
                        data_EC = EC_4_insuff)[[1]]
    
    # 377
    
  }
  
  ## 1.5 single arm design - with sufficient EC data (size = 1000)
  {
    
    eval_sa_noninfo <- function(N_RCT, 
                                alpha, 
                                beta,
                                tau,
                                data_EC,
                                r0_M = 1, 
                                r1_M = 1, 
                                r = 1, 
                                d_X = 1, 
                                dist_para = NULL,
                                sample_size = 5*10^5,
                                gamma1 = 1, 
                                gamma = 1 
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
                                r0_M = 1,
                                r1_M = 1,
                                r = 1,
                                d_X = 1, 
                                dist_para = NULL,
                                sample_size = 5*10^5,
                                gamma1 = 1, 
                                gamma = 1 
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
    
    size_sa_noninfo(start = 1,
                    end = 500,
                    alpha = 0.05,
                    beta = 0.2,
                    tau = tau_use,
                    data_EC =  EC_4_suff)[[1]]
    
    # 118
    
  }
  
  ## 1.6 single arm design - with insufficient EC data (size = 60)
  {
    
    eval_sa_verify_true <- function(alpha = 0.05, 
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
    
    eval_sa_verify_true(data_EC =  EC_4_insuff)
    
  }
  
}

## Stage 3: Use the EC data and repeatedly generate RCT data to do estimation 
{
  # repeat times
  M <- 2000
  
  ## difference in means estimation
  {
    
    f_est_diffinmeans <- function(n, 
                                  pi_A = pi_A_orig,
                                  M = 200,
                                  seed = 123,
                                  data_orig_crl = data_4_crl,
                                  data_orig_trt = data_4_trt
    ){
      
      N_4_crl <- nrow(data_orig_crl)
      N_4_trt <- nrow(data_orig_trt)
      
      ## estimation functions 
      # Naive two-sample z-test estimator
      f_est_naive <- function(data = data_RCT){
        
        N_R <- nrow(data)
        
        data_1 <- filter(data, A == 1)
        N_t <- nrow(data_1)
        data_0 <- filter(data, A == 0)
        N_c <- nrow(data_0)
        
        mu_1 <- mean(data_1$change)
        mu_0 <- mean(data_0$change)
        sigma_1 <- sd(data_1$change)
        sigma_0 <- sd(data_0$change)
        
        est_naive <- mu_1 - mu_0
        est_naive_var <- (sigma_1^2 / N_t + sigma_0^2 / N_c)
        
        return(data.frame(est_naive,
                          est_naive_var))
      }
      
      set.seed(seed)
      Index_all <- sample(1:50000000, M, replace = F)
      
      Result <- c()
      for(i in 1:M){
        print(i)
        set.seed(Index_all[i])
        
        RCT_trt_index <- sample(N_4_trt, size = floor(n * pi_A), replace = T)
        RCT_crl_index <- sample(N_4_crl, size = floor(n * (1 - pi_A)) + 1, replace = T)
        
        data_RCT_trt <- data_4_trt[RCT_trt_index, ] 
        data_RCT_crl <- data_4_crl[RCT_crl_index, ] 
        
        data_RCT <- rbind(data_RCT_trt, data_RCT_crl)
        
        data_RCT$A <- ifelse(data_RCT$THERAPY == "DRUG", 1, 0)
        
        res_std <- f_est_naive(data = data_RCT)
        
        res <- cbind(res_std)
        
        Result <- rbind(Result, res)
      }
      
      return(Result)
      
    }
    
    res_diffinmeans <- f_est_diffinmeans(n = 458,
                                         M = 2000,
                                         seed = 123) # z-test
  }
  
  ## RCT-only AIPW estimation
  {
    
    f_est_AIPW <- function(n, 
                           pi_A = pi_A_orig,
                           M = 2000,
                           seed = 123,
                           data_orig_crl = data_4_crl,
                           data_orig_trt = data_4_trt
    ){
      
      N_4_crl <- nrow(data_orig_crl)
      N_4_trt <- nrow(data_orig_trt)
      
      ## estimation functions 
      f_est_AIPW_estpi <- function(data = data_RCT){
        
        N_R <- nrow(data)
        
        pi_A_model <- glm(A ~ X1 + X2, data, family = "binomial")
        
        data_1 <- filter(data, A == 1)
        mu1_model <- lm(Y ~ X1 + X2, data_1)
        
        data_0 <- filter(data, A == 0)
        mu0_model <- lm(Y ~ X1 + X2, data_0)
        
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
          IF_est_tau <- mean(IF)
          IF_est_tau
        }
        
        ## Estimator 2 - Variance estimation - based on IF
        {
          IF_est_var <- mean(IF^2)
          IF_est_var
        }
        
        return(data.frame(Est_AIPW = IF_est_tau,  
                          Est_AIPW_var = IF_est_var / N_R))
      }
      
      set.seed(seed)
      Index_all <- sample(1:50000000, M, replace = F)
      
      Result <- c()
      for(i in 1:M){
        print(i)
        set.seed(Index_all[i])
        
        RCT_trt_index <- sample(N_4_trt, size = floor(n * pi_A), replace = T)
        RCT_crl_index <- sample(N_4_crl, size = floor(n * (1 - pi_A)) + 1, replace = T)
        
        data_RCT_trt <- data_4_trt[RCT_trt_index, ] 
        data_RCT_crl <- data_4_crl[RCT_crl_index, ] 
        
        data_RCT <- rbind(data_RCT_trt, data_RCT_crl)
        
        data_RCT$A <- ifelse(data_RCT$THERAPY == "DRUG", 1, 0)
        
        res_AIPW <- f_est_AIPW_estpi(data = data_RCT)
        
        Result <- rbind(Result, res_AIPW)
      }
      
      return(Result)
      
    }
    
    res_AIPW <- f_est_AIPW(n = 421,
                           M = 2000,
                           seed = 234) 
    
  }
  
  ## Hybrid design with sufficient EC data
  {
    
    f_est_hybrid <- function(n, 
                             pi_A = pi_A_orig,
                             M = 2000,
                             seed = 123,
                             data_EC = EC_4_suff,
                             data_orig_crl = data_4_crl,
                             data_orig_trt = data_4_trt
    ){
      
      N_EC <- nrow(data_EC)
      data_EC$R <- 0
      N_4_crl <- nrow(data_orig_crl)
      N_4_trt <- nrow(data_orig_trt)
      
      ## estimation functions 
      # Use the default formula, assume the constant conditional variance
      f_r <- function(data, 
                      N_R,
                      N_E,
                      mu1_model,
                      mu0_model){
        data_10_r <- filter(data, R == 1 & A == 0)
        epsi_10_r <- data_10_r$change -  predict(mu0_model, data_10_r)
        
        data_0_r <- filter(data, R == 0)
        epsi_0_r <- data_0_r$change -  predict(mu0_model, data_0_r)
        
        # res <- 1 / N_R * sum(epsi_10_r^2) / mean(epsi_0_r^2) 
        res <- mean(epsi_10_r^2) / mean(epsi_0_r^2) 
        
        return(res)
        
      }
      
      # Est_EC
      f_est_full_estpi <- function(data = data_full){
        
        data_11 <- filter(data,  R == 1 & A == 1)
        mu1_model <- lm(Y ~ X1 + X2, data_11)
        
        # change the estimation method for mu0 model
        data_10 <- filter(data,  R == 1 & A == 0)
        data_0 <- filter(data, R == 0)
        mu0_data <- rbind(data_10, data_0)
        mu0_model <- lm(Y ~ X1 + X2, mu0_data)
        
        data_1 <- filter(data, R == 1)
        pi_A_model <- glm(A ~ X1 + X2, data_1, family = "binomial")
        pi_R_model <- glm(R ~ X1 + X2, data, family = "binomial")
        
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
          IF_part1 <- N / N_R * (data$R * (mu1_X - mu0_X + data$A * epsi_1 / pi_A_X )) 
          IF_part2 <- N / N_R * ( 
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
        
        
        return(data.frame(Est_EC = IF_est_tau,  
                          Est_EC_var = IF_est_var / N))
      }
      
      
      set.seed(seed)
      Index_all <- sample(1:50000000, M, replace = F)
      
      Result <- c()
      for(i in 1:M){
        print(i)
        set.seed(Index_all[i])
        
        RCT_trt_index <- sample(N_4_trt, size = floor(n * pi_A), replace = T)
        RCT_crl_index <- sample(N_4_crl, size = floor(n * (1 - pi_A)) + 1, replace = T)
        
        data_RCT_trt <- data_4_trt[RCT_trt_index, ]
        data_RCT_crl <- data_4_crl[RCT_crl_index, ]
        data_RCT <- rbind(data_RCT_trt, data_RCT_crl)
        data_RCT$R <- 1
          
        data_RCT$A <- ifelse(data_RCT$THERAPY == "DRUG", 1, 0)
        data_EC$A <- 0
        data_full <- rbind(data_EC, data_RCT)
        
        res_EC <- f_est_full_estpi(data_full)
    
        Result <- rbind(Result, res_EC)
      }
      
      return(Result)
      
    }
    
    res_hybrid_suff <- f_est_hybrid(n = 238,
                                    M = 2000,
                                    seed = 345,
                                    data_EC = EC_4_suff) 
    
  }
  
  ## Hybrid design with insufficient EC data
  {
    
    res_hybrid_insuff <- f_est_hybrid(n = 377,
                                      M = 2000,
                                      seed = 456,
                                      data_EC = EC_4_insuff) 
    
  }
  
  ## Single-arm design with sufficient EC data
  {
    
    f_est_sa <- function(n = 118, 
                         M = 2000,
                         seed = 123,
                         data_EC = EC_4_suff,
                         data_orig_crl = data_4_crl,
                         data_orig_trt = data_4_trt
    ){
      
      N_EC <- nrow(data_EC)
      data_EC$R <- 0
      N_4_crl <- nrow(data_orig_crl)
      N_4_trt <- nrow(data_orig_trt)
      
      ## estimation functions 
      # Use the default formula, assume the constant conditional variance
      f_r <- function(data, 
                      N_R,
                      N_E,
                      mu1_model,
                      mu0_model){
        data_10_r <- filter(data, R == 1 & A == 0)
        epsi_10_r <- data_10_r$change -  predict(mu0_model, data_10_r)
        
        data_0_r <- filter(data, R == 0)
        epsi_0_r <- data_0_r$change -  predict(mu0_model, data_0_r)
        
        # res <- 1 / N_R * sum(epsi_10_r^2) / mean(epsi_0_r^2) 
        res <- mean(epsi_10_r^2) / mean(epsi_0_r^2) 
        
        return(res)
        
      }
      
      # Est_SA
      f_est_sa_estpi <- function(data = data_full){
        
        data_1 <- filter(data, R == 1)
        mu1_model <- lm(Y ~ X1 + X2, data_1)
        
        # change the estimation method for mu0 model
        data_0 <- filter(data, R == 0)
        mu0_model <- lm(Y ~ X1 + X2, data_0)
        
        pi_R_model <- glm(R ~ X1 + X2, data, family = "binomial")
        
        N_R <- nrow(data_1)
        N_E <- nrow(data_0)
        
        mu1_X <- predict(mu1_model, data)
        mu0_X <- predict(mu0_model, data)
        pi_A_X <- 1
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
          IF_part1 <- N / N_R * (data$R * (mu1_X - mu0_X + epsi_1)) 
          IF_part2 <- N / N_R * (((1 - data$R)) * q_X * epsi_0 )  
          
          IF <- IF_part1 - IF_part2
          IF_est_tau <- mean(IF)
          IF_est_tau
        }
        
        ## Estimator 1 - Variance estimation - based on IF
        {
          IF_est_var <- mean(IF^2)
          IF_est_var
        }
        
        
        return(data.frame(Est_SA = IF_est_tau,  
                          Est_SA_var = IF_est_var / N))
      }
      
      set.seed(seed)
      Index_all <- sample(1:50000000, M, replace = F)
      
      Result <- c()
      for(i in 1:M){
        print(i)
        set.seed(Index_all[i])
        
        RCT_trt_index <- sample(N_4_trt, size = n, replace = T)
        
        data_RCT_trt <- data_4_trt[RCT_trt_index, ]
        data_RCT <- data_RCT_trt
        data_RCT$R <- 1
        
        data_RCT$A <- ifelse(data_RCT$THERAPY == "DRUG", 1, 0)
        data_EC$A <- 0
        data_full <- rbind(data_EC, data_RCT)
        
        res_EC <- f_est_sa_estpi(data_full)
        
        Result <- rbind(Result, res_EC)
      }
      
      return(Result)
      
    }
    
    res_sa_suff <- f_est_sa(n = 118,
                            M = 2000,
                            seed = 567,
                            data_EC = EC_4_suff) 
    
  }
  
  
  apply(res_diffinmeans, 2, mean)
  apply(res_AIPW, 2, mean)
  apply(res_hybrid_suff, 2, mean)
  apply(res_hybrid_insuff, 2, mean)
  apply(res_sa_suff, 2, mean)
  
  wald <- function(mean_value, var_value){
    left <- mean_value - qnorm(0.975) * sqrt(var_value)
    right <-  mean_value + qnorm(0.975) * sqrt(var_value)
    return(c(left, right))
  }
  
  temp_diffinmeans <- apply(res_diffinmeans, 2, mean)
  c(temp_diffinmeans, wald(temp_diffinmeans[1], temp_diffinmeans[2])) |> round(digits = 2)
  
  temp_AIPW <- apply(res_AIPW, 2, mean)
  c(temp_AIPW, wald(temp_AIPW[1], temp_AIPW[2])) |> round(digits = 2)
  
  temp_hybrid_insuff <- apply(res_hybrid_insuff, 2, mean)
  c(temp_hybrid_insuff, wald(temp_hybrid_insuff[1], temp_hybrid_insuff[2])) |> round(digits = 2)
  
  temp_hybrid_suff <- apply(res_hybrid_suff, 2, mean)
  c(temp_hybrid_suff, wald(temp_hybrid_suff[1], temp_hybrid_suff[2])) |> round(digits = 2)
  
  temp_sa_suff <- apply(res_sa_suff, 2, mean)
  c(temp_sa_suff, wald(temp_sa_suff[1], temp_sa_suff[2])) |> round(digits = 2)
 
   
}

