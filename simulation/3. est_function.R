
# Naive two-sample z-test estimator
f_est_naive <- function(data = data_RCT,
                        delta = 0){
  A <- data$A
  X1 <- data$X1
  X2 <- data$X2
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

# IF-based estimator for design based on RCT-based AIPW
f_est_AIPW <- function(data = data_RCT){
  A <- data$A
  X1 <- data$X1
  X2 <- data$X2
  Y <- data$Y
  
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
                    Est_AIPW_var = IF_est_var))
}

# Use the r(X) formula provided in Gao 2025's Remark 1
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

# IF-based estimator for hybrid design
f_est_hybrid <- function(data = data_full){
  R <- data$R
  A <- data$A
  X1 <- data$X1
  X2 <- data$X2
  Y <- data$Y
  
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
                    Est_EC_var = IF_est_var))
}

# IF-based estimator for single-arm design
f_est_sa <- function(data = data_full){
  R <- data$R
  A <- data$A
  X1 <- data$X1
  X2 <- data$X2
  Y <- data$Y
  
  # mu1_model
  data_1 <- filter(data,  R == 1)
  mu1_model <- lm(Y ~ X1 + X2, data_1)
  
  # mu0 model
  data_0 <- filter(data, R == 0)
  mu0_model <- lm(Y ~ X1 + X2, data_0)
  
  # R's model
  pi_R_model <- glm(R ~ X1 + X2, data, family = "binomial")
  
  N_R <- nrow(data_1)
  N_E <- nrow(data_0)
  
  mu1_X <- predict(mu1_model, data)
  mu0_X <- predict(mu0_model, data)
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
                    Est_SA_var = IF_est_var))
}