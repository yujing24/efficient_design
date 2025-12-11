## Data generation functions for simulation section in the main paper
{
  # Generate outcome Y based on linear regression in EC data
  f_gen_EC_Y <- function(x1,
                         x2,
                         beta0_EC,
                         beta1_EC,
                         beta2_EC, 
                         sdY_EC){
    mu <- beta0_EC + beta1_EC * x1 + beta2_EC*x2  
    return(rnorm(1, mean = mu, sd = sdY_EC))
  }
  
  # Generate EC data
  f_gen_EC <- function(seed_EC = 123,
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
  
  # Generate outcome Y based on linear regression in RCT data
  f_gen_RCT_Y <- function(x1,
                          x2,
                          a,
                          beta0_RCT,
                          beta1_RCT,
                          beta2_RCT, 
                          tau,
                          sdY_RCT){
    mu <- beta0_RCT + beta1_RCT * x1 + beta2_RCT * x2 + tau * a
    return(mu + rnorm(1, mean = 0, sd = sdY_RCT))
  }
  
  # Generate RCT data
  f_gen_RCT <- function(seed_RCT = 123,
                        N_RCT, # sample size of RCT data
                        pi_A, # propensity score in RCT
                        tau, # treatment effect
                        muX1_RCT = 1, 
                        sdX1_RCT = 1,
                        pX2_RCT = 0.5,
                        beta0_RCT = 1,
                        beta1_RCT = 0.5,
                        beta2_RCT = -1,
                        sdY_RCT = sqrt(0.8)
  ){
    
    set.seed(seed_RCT)
    
    N_crl <- ceiling(N_RCT * (1 - pi_A)) 
    N_trt <- N_RCT - N_crl 
    
    # generate the pre-treatment covariates
    X1_RCT <- rnorm(N_RCT, mean = muX1_RCT, sd = sdX1_RCT)
    X2_RCT <- rbinom(N_RCT, size = 1, prob = pX2_RCT)
    
    # generate A by using pi
    A_trt <- rep(1, N_trt)
    A_crl <- rep(0, N_crl)
    A <- c(A_trt, A_crl)
    
    Y_RCT <- mapply(f_gen_RCT_Y,
                    x1 = X1_RCT,
                    x2 = X2_RCT,
                    a = A,
                    beta0_RCT,
                    beta1_RCT,
                    beta2_RCT, 
                    tau,
                    sdY_RCT)
    
    return(data.frame(A = A, X1 = X1_RCT, X2 = X2_RCT,  Y = Y_RCT))
    
  }
  
  # combine data_RCT and data_EC
  f_gen_full <- function(data_EC,
                         data_RCT){
    
    data_full <- rbind(data_RCT, data_EC)
    data_full$R <- c(rep(1, nrow(data_RCT)), rep(0, nrow(data_EC)) )
    return(data_full)
  }
  
}


## Additional data generation functions for sensitivity analysis 1
{
  # Generate outcome Y based on linear regression in RCT data
  f_gen_RCT_Y_sensi1 <- function(x1,
                                 x2,
                                 a,
                                 beta0_RCT,
                                 beta1_RCT,
                                 beta2_RCT, 
                                 tau1,
                                 tau2,
                                 sdY_RCT
  ){
    mu <- beta0_RCT + beta1_RCT * x1 + beta2_RCT * x2 + tau1 * a * x1 + tau2 * a * x2
    return(mu + rnorm(1, mean = 0, sd = sdY_RCT))
  }
  
  # Generate RCT data
  f_gen_RCT_sensi1 <- function(seed_RCT = 123,
                               N_RCT, # sample size of RCT data
                               pi_A, # propensity score in RCT
                               tau, # treatment effect
                               muX1_RCT = 1, 
                               sdX1_RCT = 1,
                               pX2_RCT = 0.5,
                               beta0_RCT = 1,
                               beta1_RCT = 0.5,
                               beta2_RCT = -1,
                               tau1,
                               tau2,
                               sdY_RCT = sqrt(0.8)
  ){
    
    set.seed(seed_RCT)
    
    N_crl <- ceiling(N_RCT * (1 - pi_A)) 
    N_trt <- N_RCT - N_crl 
    
    # generate the pre-treatment covariates
    X1_RCT <- rnorm(N_RCT, mean = muX1_RCT, sd = sdX1_RCT)
    X2_RCT <- rbinom(N_RCT, size = 1, prob = pX2_RCT)
    
    # generate A by using pi
    A_trt <- rep(1, N_trt)
    A_crl <- rep(0, N_crl)
    A <- c(A_trt, A_crl)
    
    Y_RCT <- mapply(f_gen_RCT_Y_sensi1,
                    x1 = X1_RCT,
                    x2 = X2_RCT,
                    a = A,
                    beta0_RCT,
                    beta1_RCT,
                    beta2_RCT, 
                    tau1,
                    tau2,
                    sdY_RCT)
    
    return(data.frame(A = A, X1 = X1_RCT, X2 = X2_RCT,  Y = Y_RCT))
    
  }
}


## Additional data generation functions for sensitivity analysis 2
{
  # Generate outcome Y based on linear regression in EC data
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
  
  # Generate EC data
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
  
  # Generate outcome Y based on linear regression in RCT data
  f_gen_RCT_Y_sensi2 <- function(x1,
                                x2,
                                a,
                                beta0_RCT,
                                beta1_RCT,
                                beta2_RCT, 
                                tau,
                                epi_RCT,
                                sdY_RCT){
    
    mu <- beta0_RCT + beta1_RCT * x1 + beta2_RCT * x2 + tau * a
    epi <- rnorm(1, mean = 0, sd = sdY_RCT)
    
    return(mu +  epi_RCT * x1 * epi)
  }
  
  # Generate RCT data
  f_gen_RCT_sensi2 <- function(seed_RCT = 123,
                                N_RCT, # sample size of RCT data
                                pi_A, # propensity score in RCT
                                muX1_RCT = 1, 
                                sdX1_RCT = 1,
                                pX2_RCT = 0.5,
                                beta0_RCT = 1,
                                beta1_RCT = 0.5,
                                beta2_RCT = -1,
                                tau = tau,
                                epi_RCT = 0.8,
                                sdY_RCT = sqrt(0.8)
  ){
    
    set.seed(seed_RCT)
    
    N_crl <- ceiling(N_RCT * (1 - pi_A)) 
    N_trt <- N_RCT - N_crl 
    
    # generate the pre-treatment covariates
    X1_RCT <- rnorm(N_RCT, mean = muX1_RCT, sd = sdX1_RCT)
    X2_RCT <- rbinom(N_RCT, size = 1, prob = pX2_RCT)
    
    # generate A by using pi
    A_trt <- rep(1, N_trt)
    A_crl <- rep(0, N_crl)
    A <- c(A_trt, A_crl)
    
    Y_RCT <- mapply(f_gen_RCT_Y_sensi2,
                    x1 = X1_RCT,
                    x2 = X2_RCT,
                    a = A,
                    beta0_RCT,
                    beta1_RCT,
                    beta2_RCT, 
                    tau,
                    epi_RCT,
                    sdY_RCT)
    
    return(data.frame(A = A, X1 = X1_RCT, X2 = X2_RCT, Y = Y_RCT))
    
  }
}

