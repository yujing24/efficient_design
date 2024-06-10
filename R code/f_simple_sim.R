rm(list=ls())


f_gen_sim <- function(
    seed = ,
    n = 200, # sample size 
    muX = 1, 
    sdX = 1
    ){
  
  # generate the pre-treatment covariates
  X <- rnorm(n, mean = muX, sd = sdX)
  
}